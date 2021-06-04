/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gnuradio/io_signature.h>

#include <openemissions/histogram_solve.h>

// with c++17 this can be removed in favor of std::apply
#include <boost/fusion/functional/invocation/invoke.hpp>
// #include <boost/fusion/adapted/std_array.hpp> // boost ~1.63
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/algorithm/iteration/fold.hpp>

namespace gr {
namespace openemissions {

namespace detail {

// a square volume with a compile-time number of dimensions
// i think these were standardised or are on a standards track, somewhere, not sure ...
// probably in boost somewhere at least
template <typename T, size_t N>
class square_ctd : public std::vector<T>
{
public:
  square_ctd(size_t edgelength = 0, const T & value = {})
  : d_edgelength(0)
  {
    resize(edgelength, value);
  }

  void resize(size_t edgelength, const T & value = {})
  {
    size_t size = 1;
    for (size_t dim = 0; dim < N; dim ++) {
      size *= edgelength;
    }
    std::vector<T>::reserve(size);
    std::vector<T>::resize(size, value);
  }

  size_t size() const
  {
    return d_edgelength;
  }

  T & operator[](const size_t * coords)
  {
    return std::vector<T>::operator[](coords2coord(coords));
  }

  const T & operator[](const size_t * coords) const
  {
    return std::vector<T>::operator[](coords2coord(coords));
  }

private:
  size_t d_edgelength;

  constexpr size_t coords2coord(const size_t * coords, size_t skip_idx = ~0)
  {
    size_t coord = 0;
    for (size_t idx = 0; idx < N + (skip_idx == ~0 ? 0 : 1); idx ++) {
      if (skip_idx == ~0 || idx != skip_idx) {
        coord = (coord * d_edgelength) + coords[idx];
      }
    }
    return coord;
  }
};

}

template <typename freq_type, typename... Doubles>
class histogram_solve_impl : public histogram_solve<freq_type, Doubles...>
{
private:
  double d_min;
  double d_max;
  std::function<double(Doubles...)> d_expr;
  std::vector<std::function<std::array<double, sizeof...(Doubles)>(Doubles...)>> d_extrema;
  size_t d_output_idx;
  size_t d_param_output_idx;
  size_t d_nbuckets;

  double d_bucket_to_value_coeff;
  double d_bucket_to_value_offset;
  double d_value_to_bucket_coeff;
  double d_value_to_bucket_offset;

  std::vector<double> d_unk_vec;
  std::vector<double> d_pop_vec;
  std::vector<bool> d_pop_mask;

  bool d_warned;

  // sparse_dot is used to calculate each output value and could likely be optimised
  class sparse_dot {
  public:
      std::vector<size_t> index_sequences;
      std::vector<double> proportions;
      
      void clear() {
          index_sequences.clear();
      }
  
      void add_sequence(size_t *sequence, size_t skip_idx = ~0, double proportion = 1) {
          // each sequence lists inputs that combine for this result
          for (size_t w = 0; w < sizeof...(Doubles)/* + 1*/; w ++) {
            if (w != skip_idx/* && w != result_idx()*/) {
              index_sequences.push_back(sequence[w]);
            }
          }
          proportions.push_back(proportion);
      }
  
      template <typename T>
      T calc(const T ** knowns, size_t offset, size_t output_idx = ~0, bool * zero_density = 0) {
          T result = 0;
          bool found_zero_density = true;
          for (size_t sequence = 0, proportion = 0;
               sequence < index_sequences.size();
               proportion ++)
          {
              T product = 1;
              size_t knownidx = 0;
              for (size_t varnum = 0; varnum < sizeof...(Doubles) + 1; varnum ++) {
                if (varnum != output_idx) {
                    if (varnum != result_idx()) {
                        // the number of outcomes made by each combination of inputs
                        // is the combined product of the number of outcomes they hold
                        found_zero_density = false;
                        T outcomes = knowns[knownidx][offset + index_sequences[sequence]];
                        sequence ++;
                        if (outcomes == 0) {
                          product = 0;
                        } else {
                          product *= outcomes;
                        }
                    }
                    knownidx ++;
                }
              }
              result += product * proportions[proportion];
          }
          if (0 != zero_density) {
              *zero_density = found_zero_density;
          }
          return result;
      }
  };

  // if the result is the output unknown,
  // each output item is the sum of the products of sequences of input items,
  // different for each output item.
  // otherwise the components made by each possible unknown value are broken out
  std::vector<std::vector<sparse_dot>> d_unknown_result_map;

public:
  /*
   * The private constructor
   */
  histogram_solve_impl(double min, double max, const std::function<double(Doubles...)> & expr, size_t output_idx, size_t nbuckets, const std::vector<std::function<std::array<double, sizeof...(Doubles)>(Doubles...)>> & extrema)
    : gr::sync_block("histogram_solve",
            gr::io_signature::make(sizeof...(Doubles), sizeof...(Doubles), sizeof(freq_type) * nbuckets),
            gr::io_signature::make(1, output_is_forward_solving(output_idx) ? 1 : 2, sizeof(freq_type) * nbuckets)),
      d_min(min),
      d_max(max),
      d_expr(expr),
      d_output_idx(output_idx),
      d_param_output_idx(output_idx == result_idx() ? ~0 : var_idx_to_param_idx(output_idx)),
      d_nbuckets(nbuckets)
  {
    update_coeffs();
  }
  
  /*
   * Our virtual destructor.
   */
  ~histogram_solve_impl()
  {
  }

  bool start() override
  {
    d_warned = false;
    return sync_block::start();
  }

  int
  work(int noutput_items,
       gr_vector_const_void_star &input_items,
       gr_vector_void_star &output_items)
  {
    const freq_type **in = reinterpret_cast<const freq_type**>(input_items.data());
    freq_type *out = reinterpret_cast<freq_type*>(output_items[0]);
    //memset(out, 0, sizeof(freq_type) * noutput_items);

    std::vector<double> calculated_result_histogram(d_nbuckets);
    std::map<double, double> calculated_result_ranges;
    std::vector<double> calculated_result_independent_probabilities(d_nbuckets);

    for (size_t item = 0; item < noutput_items; item ++) {

      if (is_forward_solving()) {
        // forward solving
        for (size_t unknown_bucket = 0; unknown_bucket < d_nbuckets; unknown_bucket ++) {
          freq_type & unknown_proportion = out[unknown_bucket];

          unknown_proportion = d_unknown_result_map[0][unknown_bucket].calc(in, item * d_nbuckets, d_output_idx);
        }
      } else {
        // reverse solving
        
        double total_outcomes = 0;
        double max_sample = 0;
        d_unk_vec.resize(d_nbuckets);
        d_pop_vec.resize(d_nbuckets);
        d_pop_mask.resize(d_nbuckets);

        for (size_t result_bucket = 0; result_bucket < d_nbuckets; result_bucket ++) {
            const double & sample = in[result_idx()][item * d_nbuckets + result_bucket];
            if (sample > max_sample) {
                max_sample = sample;
            }
        }


        size_t unknowns_leaving_model = 0;
        for (size_t unknown_bucket = 0; unknown_bucket < d_nbuckets; unknown_bucket ++) {
            double max_population = 0;
            bool this_unknown_leaves_model = false;
            for (size_t result_bucket = 0; result_bucket < d_nbuckets; result_bucket ++) {
                bool found_zero;
                double & calculated_result_population = d_pop_vec[result_bucket];
                calculated_result_population = d_unknown_result_map[unknown_bucket][result_bucket].calc(in, item * d_nbuckets, d_output_idx, &found_zero);
                bool leaves_model;
                if (calculated_result_population > max_population) {
                    max_population = calculated_result_population;
                    leaves_model = false;
                } else {
                    leaves_model =
                            found_zero &&
                            in[result_idx()][item * d_nbuckets + result_bucket] > 0;
                    this_unknown_leaves_model |= leaves_model;
                }
                d_pop_mask[result_bucket] = leaves_model;
            }
            if (this_unknown_leaves_model) {
                unknowns_leaving_model ++;
            }
            double & unknown_result_outcomes = d_unk_vec[unknown_bucket];
            unknown_result_outcomes = max_population;
            for (size_t result_bucket = 0; result_bucket < d_nbuckets && unknown_result_outcomes != 0; result_bucket ++) {
                // We interpret the calculated histogram as a population,
                // and count how many samplings give the measured bucket.
                //
                // The equation for the count of samplings that give a single bucket:
                // Combinations(sampcount, sampbinheight) * pow(popbinheight, sampbinheight) * pow(popcount - popbinheight, sampcount - sampbinheight)
                //
                // The equation for the count of samplings that give a specific histogram:
                // sum[ Combinations(sampcountremaining, sampbinheight) * pow(popbinheight, sampbinheight) ]
                //
                // We're considering the whole histogram, for every unknown.
                //
                // Since we care about the ratio between the unknowns, the only thing that changes is popbinheight and the combinations terms can be factored out and ignored.
                // The other terms are divided by their maximums to prevent overflow
                //
                // TODO: speed up
                //  - we could consider only dense bins
                //  - we could approximate the numerical result
                //  - we could only consider the last item in the buffer
                //  - we could update based on a user-provided frequency
                //  the user can also decimate the incoming histograms.

                double sampbinheight = in[result_idx()][item * d_nbuckets + result_bucket] / max_sample;
                //if (d_pop_mask[result_bucket]) {
                    double factor = pow(d_pop_vec[result_bucket] / max_population, sampbinheight);
                    unknown_result_outcomes *= factor;
                //} else {
                //    // this condition happens if some of the result data is unmodelled
                //    // this could happen if the bounds are too small, or if sampled data
                //    // is stretched ...
                //    // a wrong answer will result.  there are various ways to try to fix it,
                //    // but the best answer is achieved if the user avoids the situation.
                //}
            }
            total_outcomes += unknown_result_outcomes;
        }
        /*
        // basically, if data from one of the inputs leaves the bounds of a histogram,
        // then this block produces quite wrong results.
        // but the check below seems to end up being roughly unrelated to that.
        if (unknowns_leaving_model >= d_nbuckets / 4 && !d_warned) {
            GR_LOG_WARN(this->d_logger, "Result has density in unmodelled areas for " + std::to_string(unknowns_leaving_model * 100 / d_nbuckets) + "% of the unknown values.");
            GR_LOG_WARN(this->d_logger, "Is background noise included in equation?  Do bounds of histogram surround signal?");
            GR_LOG_WARN(this->d_logger, "When calculating models, only data within histograms is included.  It is expected that data coming from outside histogram bounds is nonpresent or simulated.");
            GR_LOG_WARN(this->d_logger, "Further warnings from this block suppressed.");
            d_warned = true;
        }
        */
        for (size_t unknown_bucket = 0; unknown_bucket < d_nbuckets; unknown_bucket ++) {
            out[unknown_bucket] = d_unk_vec[unknown_bucket];
        }
      }

      out += d_nbuckets;
    }
  
    return noutput_items;
  }

  void set_min(double min) override
  {
    d_min = min;
    update_coeffs();
  } 

  double min() const override
  {
    return d_min;
  }

  void set_max(double max) override
  {
    d_max = max;
    update_coeffs();
  }

  double max() const override
  {
    return d_max;
  }

  void set_expr(const std::function<double(Doubles...)> &expr) override
  {
    d_expr = expr;
  }

  const std::function<double(Doubles...)> & expr() const override
  {
    return d_expr;
  }

  void set_extrema(const std::vector<std::function<std::array<double, sizeof...(Doubles)>(Doubles...)>> &extrema) override
  {
    d_extrema = extrema;
  }

  const std::vector<std::function<std::array<double, sizeof...(Doubles)>(Doubles...)>> & extrema() const override
  {
    return d_extrema;
  }

  size_t output_idx() const override
  {
    return d_output_idx;
  }

  size_t nbuckets() const override
  {
    return d_nbuckets;
  }

private:
  struct map_approximation_work
  {
    struct point_work {
      // variable names could be clearer.  these two are result values.
      double value;
      size_t bucket;
      // whereas these are input parameter values.  stored so many times to ease design, likely could be reduced.
      std::array<double, sizeof...(Doubles)> parameters;
      std::vector<std::array<double, sizeof...(Doubles)>> extrema_values;
    };
    detail::square_ctd<point_work, sizeof...(Doubles) - 1> last_ndplane;
    detail::square_ctd<point_work, sizeof...(Doubles) - 1> this_ndplane;
    size_t input_buckets[sizeof...(Doubles)];
  };

  void update_coeffs()
  {
    double range = (double)d_max - d_min;

    d_bucket_to_value_coeff = range / d_nbuckets;
    d_bucket_to_value_offset = d_min;
    d_value_to_bucket_coeff = d_nbuckets / range;
    d_value_to_bucket_offset = -d_value_to_bucket_coeff * d_min;

    d_unknown_result_map.clear();
    d_unknown_result_map.resize(d_nbuckets, std::vector<sparse_dot>{d_nbuckets});

    map_approximation_work interpolation;
    size_t * bucket_ptr = interpolation.input_buckets;
    size_t & outermost_bucket = *bucket_ptr;
    interpolation.this_ndplane.resize(d_nbuckets + 1);
    interpolation.last_ndplane.resize(d_nbuckets + 1);
    for (outermost_bucket = 0; outermost_bucket <= d_nbuckets; outermost_bucket ++) {
      update_expr_map(interpolation, bucket_ptr, outermost_bucket > 0, bucket_to_value(outermost_bucket));
      interpolation.this_ndplane.swap(interpolation.last_ndplane);
    }
  }

  // recursive template function walks each variable, trying all combinations to see what result they produce
  // this function has so many arguments that i'd like it if this algorithm were eventually bundled into a class with member variables
  // this recursion to enumerate a parameter pack can probably also be replaced with something from boost or the c++ std
  template <typename... Params>
  void update_expr_map(map_approximation_work & interpolation, size_t *bucket_ptr, bool do_interpolation, Params... params)
  {
    bucket_ptr ++;
    size_t & bucket = *bucket_ptr;
    for (bucket = 0; bucket <= d_nbuckets; bucket ++) {
      update_expr_map(interpolation, bucket_ptr, do_interpolation && bucket > 0, bucket_to_value(bucket), params...);
    }
  }

  void update_expr_map(map_approximation_work & interpolation, size_t *bucket_ptr, bool do_interpolation, Doubles... params)
  {
    double value = d_expr(params...);
    interpolation.this_ndplane[interpolation.input_buckets + 1].value = value;
    interpolation.this_ndplane[interpolation.input_buckets + 1].bucket = value_to_bucket(value);
    interpolation.this_ndplane[interpolation.input_buckets + 1].parameters = { params... };
    interpolation.this_ndplane[interpolation.input_buckets + 1].extrema_values.resize(d_extrema.size());
    for (size_t w = 0; w < d_extrema.size(); w++) {
        auto extrema_values = d_extrema[w](params...);
        interpolation.this_ndplane[interpolation.input_buckets + 1].extrema_values[w] = extrema_values;
    }

    if (! do_interpolation) {
      return;
    }

    // prev_buckets and input_buckets describe a volume of output values, with the first axis held by last_ndplane vs this_ndplane.
    // the goal is to stretch them over the result buckets with appropriate linearly interpolated densities.
    // the fraction of the volume that covers a bucket, is its weight.
    // todo: unify prev_buckets and difference_indices to organise a little

    size_t prev_buckets[sizeof...(Doubles)];
    for (size_t w = 0; w < sizeof...(Doubles); w ++) {
      // this could have been tracked next to bucket_ptr to not loop over the coords every call here?
      prev_buckets[w] = interpolation.input_buckets[w] - 1;
    }

    // the simplest solution involves finding the minimum and maximum result for the bucket.
    // this is a little inaccurate but ensures all possibilities are included.
    size_t difference_indices[sizeof...(Doubles)];
    for (size_t w = 0; w < sizeof...(Doubles); w ++) {
      difference_indices[w] = interpolation.input_buckets[w] - 1;
    }
    std::array<double, sizeof...(Doubles)> & minimum_params = interpolation.last_ndplane[difference_indices + 1].parameters;
    std::array<double, sizeof...(Doubles)> & maximum_params = interpolation.this_ndplane[interpolation.input_buckets + 1].parameters;
    double minimum, maximum;
    minimum = maximum = interpolation.last_ndplane[difference_indices + 1].value;
    while (true) {
      // walk all corner coordinates via base 2 incrementation
      size_t idx = 0;
      while (idx < sizeof...(Doubles)) {
        difference_indices[idx] ++;
        if (difference_indices[idx] <= interpolation.input_buckets[idx]) {
          break;
        }
        difference_indices[idx] = interpolation.input_buckets[idx] - 1;
        idx ++;
      }
      if (idx == sizeof...(Doubles)) {
        break;
      }

      // store min/max
      auto & ndplane = difference_indices[0] ? interpolation.last_ndplane : interpolation.this_ndplane;
      auto & point_work = ndplane[difference_indices + 1];
      double & value = point_work.value;
      if (value < minimum) {
        minimum = value;
      }
      if (value > maximum) {
        maximum = value;
      }
      for (size_t w = 0; w < point_work.extrema_values.size(); w ++) {
        bool within_range = true;
        for (size_t x = 0; x < sizeof...(Doubles) && within_range; x++) {
          double & extreme_parameter = point_work.extrema_values[w][x];
          if (extreme_parameter < minimum_params[x] || extreme_parameter > maximum_params[x]) {
            within_range = false;
          }
        }
        if (within_range) {
                     // i had boost 1.53 so converted to a fusion::vector.  use of invoke directly on std::array was probably added around boost 1.63 .
          boost::fusion::vector<Doubles...> parameters;
          boost::fusion::fold(parameters, 0, [&](int idx, double & val) { val = point_work.extrema_values[w][idx]; return idx + 1; });
                                 // can be replaced in c++17 with std::apply
          double extreme_result = boost::fusion::invoke(d_expr, parameters);

          if (extreme_result < minimum) {
            minimum = extreme_result;
          }
          if (extreme_result > maximum) {
            maximum = extreme_result;
          }
        }
      }
    }

    double result_range = maximum - minimum;

    auto & result_maps = is_reverse_solving() ? d_unknown_result_map[prev_buckets[d_param_output_idx]] : d_unknown_result_map[0];

    // for now values are given a very simple linear interpolation between min and max

    size_t first_bucket = value_to_bucket(minimum);
    size_t last_bucket = value_to_bucket(maximum);
    size_t head_bucket = first_bucket < 0 ? 0 : first_bucket;
    size_t tail_bucket = last_bucket >= d_nbuckets ? d_nbuckets - 1 : last_bucket;

    for (size_t cur_bucket = head_bucket; cur_bucket <= tail_bucket; cur_bucket ++) {
      double density;
      if (cur_bucket == last_bucket) {
        if (cur_bucket == first_bucket) {
          density = 1.0;
        } else {
          density = (maximum - bucket_to_value(last_bucket)) / result_range;
        }
      } else if (cur_bucket == first_bucket) {
        density = (bucket_to_value(first_bucket + 1) - minimum) / result_range;
      } else {
        density = d_bucket_to_value_coeff / result_range;
      }
      // prev_buckets could be the same variable as difference_indices
      if (density > 0) {
        result_maps[cur_bucket].add_sequence(prev_buckets, d_param_output_idx, density);
      }
    }
  }

  inline double bucket_to_value(size_t bucket)
  {
    return bucket * d_bucket_to_value_coeff + d_bucket_to_value_offset;
  }

  inline size_t value_to_bucket(double value)
  {
    return value * d_value_to_bucket_coeff + d_value_to_bucket_offset;
  }
  
  // index types.  variable index is what is passed to the constructor for the output.
  // with variable indices, all variables are in order.  the output is at d_output_idx.  the result is at result_idx().

  // input indices index the block inputs and have no output in them.
  // with input indices, the total is one fewer, and the output is skipped.
  size_t var_idx_to_input_idx(size_t var_idx) const
  {
    return var_idx < d_output_idx ? var_idx : var_idx - 1;
  }

  // parameter indices index the expression parameters.
  // with parameter indices, the total is one fewer, and the result is skipped.
  // when the result is the output (forward-solving), the parameter index happens to be the same as the input index.
  size_t var_idx_to_param_idx(size_t var_idx) const
  {
    return var_idx < result_idx() ? var_idx : var_idx - 1;
  }

  size_t result_input_idx() const
  {
    size_t constexpr var_idx = result_idx();
    if (var_idx == 0) {
      return var_idx;
    } else {
      return var_idx_to_input_idx(var_idx);
    }
  }

  // providing the result index with this member function
  // gives some space towards letting the user set their own variable/input order in the UI some day
  static constexpr size_t result_idx()
  {
    return 0;
  }

  static constexpr size_t first_parameter_idx()
  {
    return result_idx() == 0 ? 1 : 0;
  }

  inline bool is_forward_solving() const
  {
    return output_is_forward_solving(d_output_idx);
  }

  inline bool is_reverse_solving() const
  {
    return output_is_reverse_solving(d_output_idx);
  }

  static constexpr bool output_is_forward_solving(size_t output_idx)
  {
    return result_idx() == output_idx;
  }

  static constexpr bool output_is_reverse_solving(size_t output_idx)
  {
    return result_idx() != output_idx;
  }
};

template <typename freq_type, typename... Doubles>
typename histogram_solve<freq_type, Doubles...>::sptr
histogram_solve<freq_type, Doubles...>::make(double min, double max, const std::function<double(Doubles...)> & expr, size_t output_idx, size_t nbuckets, const std::vector<std::function<std::array<double, sizeof...(Doubles)>(Doubles...)>> & extrema)
{
  return gnuradio::make_block_sptr<histogram_solve_impl<freq_type, Doubles...>>(
    min, max, expr, output_idx, nbuckets, extrema);
}

template class histogram_solve<double,   double>;
template class histogram_solve<float,    double>;
template class histogram_solve<uint64_t, double>;
template class histogram_solve<double,   double, double>;
template class histogram_solve<float,    double, double>;
template class histogram_solve<uint64_t, double, double>;
template class histogram_solve<double,   double, double, double>;
template class histogram_solve<float,    double, double, double>;
template class histogram_solve<uint64_t, double, double, double>;

} /* namespace openemissions */
} /* namespace gr */

