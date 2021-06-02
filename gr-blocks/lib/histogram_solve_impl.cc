/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gnuradio/io_signature.h>

#include <openemissions/histogram_solve.h>

//#include <lapacke/lapacke.h>

namespace gr {
namespace openemissions {

namespace detail {

// a square volume with a compile-time number of dimensions
// i think these were standardised or are on a standards track, somewhere, not sure ...
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

  // i was getting confused mapping indices, so this is just done here for now
  T & operator()(const size_t * coords, size_t skip_idx = ~0)
  {
    return std::vector<T>::operator[](coords2coord(coords, skip_idx));
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
              proportions.push_back(proportion);
            }
          }
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
                        product *= knowns[knownidx][offset + index_sequences[sequence]];
                        found_zero_density = false;
                        sequence ++;
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
  histogram_solve_impl(double min, double max, std::function<double(Doubles...)> expr, size_t output_idx, size_t nbuckets)
    : gr::sync_block("histogram_solve",
            gr::io_signature::make(sizeof...(Doubles), sizeof...(Doubles), sizeof(freq_type) * nbuckets),
            gr::io_signature::make(1, output_is_forward_solving(output_idx) ? 1 : 2, sizeof(freq_type) * nbuckets)),
      d_min(min),
      d_max(max),
      d_expr(expr),
      d_output_idx(output_idx),
      d_param_output_idx(var_idx_to_param_idx(output_idx)),
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
            for (size_t result_bucket = 0; result_bucket < d_nbuckets; result_bucket ++) {
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
                    unknown_result_outcomes *= pow(d_pop_vec[result_bucket] / max_population, sampbinheight);
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
            out[unknown_bucket] = d_unk_vec[unknown_bucket]; // / total_outcomes;
        }
#if 0
        /*
          If a conditional probability occurs in a finite set of mutually exclusive conditions,
          then the unconditional probability is: 

          P(A) = P(A | C_1) * P(C_1) + ... + P(A | C_n) * P(C_n)

          This can be applied to our conditions here.  We calculate the probability of each
          result assuming each possible unknown value.  The unknown values form mutually exclusive
          conditions:

          P(R_x) = P(R_x | Unk_1) * P(Unk_1) + ... + P(R_x | Unk_n) * P(Unk_n)

          This same relation can be stated as a simple matrix equation:

              .--- the error is in how this vector is built.
              v    We sampled from P(R_x | Unk_real), not P(R_x) which here would be independent of Unk_real.
          [ P(R_1) ]   [ P(R_1 | Unk_1) ... P(R_1 | Unk_n) ]   [ P(Unk_1) ]
          [   ...  ] = [       ...      ...        ...     ] x [    ...   ]
          [ P(R_n) ]   [ P(R_n | Unk_1) ... P(R_n | Unk_n) ]   [ P(Unk_n) ]
          
          The matrix of conditionals in the middle, calculated from d_result_map,
          can be inverted and multiplied by the result vector on the left, which is
          the result input histogram, to calculate the probabilities of all possible
          unknown values on the right.
        */
        
        d_mat.resize(d_nbuckets * d_nbuckets);
        d_vec.resize(d_nbuckets);
        d_mat_p.resize(d_nbuckets);
        double * conditional = d_mat.data();
        double total = 0;
        // my openblas library crashes when used with row major data, so this is definitely col major
        for (size_t unknown_bucket = 0; unknown_bucket < d_nbuckets; unknown_bucket ++) {
          d_vec[unknown_bucket] = in[result_input_idx()][item * d_nbuckets + unknown_bucket];
          for (size_t result_bucket = 0; result_bucket < d_nbuckets; result_bucket ++) {
            *conditional = d_unknown_result_map[unknown_bucket][result_bucket].calc(in, item * d_nbuckets, result_idx());
            total += *conditional;
            conditional ++;
          }
        }
        for (size_t result_bucket = 0; result_bucket < d_nbuckets; result_bucket ++) {
          std::cerr << "[\t" << d_vec[result_bucket] << "\t]\t";
          if (result_bucket == d_nbuckets / 2) {
            std::cerr << "=";
          }
          for (size_t unknown_bucket = 0; unknown_bucket < d_nbuckets; unknown_bucket ++) {
            std::cerr << "\t" << d_mat[result_bucket + unknown_bucket * d_nbuckets] << "=" <<
                d_unknown_result_map[unknown_bucket][result_bucket].calc(in, item * d_nbuckets, result_idx());
          }
          if (result_bucket == d_nbuckets / 2) {
            std::cerr << "\t*\tX";
          }
          std::cerr << std::endl;
        }

        // solve
        //int status = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', d_nbuckets, d_nbuckets, 1, d_mat.data(), d_nbuckets, d_vec.data(), d_nbuckets);
        int status = LAPACKE_dgesv(LAPACK_COL_MAJOR, d_nbuckets, 1, d_mat.data(), d_nbuckets, d_mat_p.data(), d_vec.data(), d_nbuckets);
        if (status < 0) {
          throw std::runtime_error("argument " + std::to_string(-status) + " to lapack has illegal value");
        } else if (status > 0) {
          throw std::runtime_error("matrix factor is exactly singular, U(" + std::to_string(status) + ") == 0");
          //throw std::runtime_error("could not compute least squares solution, A triangular factor diagonal " + std::to_string(status) + " == 0");
        }

        for (size_t bucket = 0; bucket < d_nbuckets; bucket ++) {
          out[bucket] = d_vec[bucket] * total;
        }
        
#endif
        
#if 0
        double measured_result_total = 0;
        for (size_t result_bucket = 0; result_bucket < d_nbuckets; result_bucket ++) {
          measured_result_total += in[result_input_idx()][result_bucket];
        }

        // we form distributions for each potential value of the unknown variable
        for (size_t unknown_bucket = 0; unknown_bucket < d_nbuckets; unknown_bucket ++) {
          double unknown_proportion;

          // given the potential value, we form the histogram of the result variable
          // this histogram is placed into calculated_result_histogram

          double calculated_result_total = 0;

          for (size_t result_bucket = 0; result_bucket < d_nbuckets; result_bucket ++) {
            double frequency = d_unknown_result_map[unknown_bucket][result_bucket].calc(in, item * d_nbuckets);
            calculated_result_total += frequency;
            calculated_result_histogram[result_bucket] = frequency;
          }

          /*
          
          std::map<double, double> frequency_density; // for calculating P(bucket value)
          double last_frequency;
          for (size_t result_bucket = 0; result_bucket < d_nbuckets; result_bucket ++) {
            double frequency = d_unknown_result_map[unknown_bucket][result_bucket].calc(in, item * d_nbuckets);
            calculated_result_total += frequency;
            calculated_result_histogram[result_bucket] = frequency;
            if (result_bucket > 0) {
              double range_min = last_frequency;
              double range_max = frequency;
              if (frequency < last_frequency) {
                std::swap(range_min, range_max);
              }
              double range_density = 1 / (range_max - range_min);
              // first we store the changes, then later integrate them to look up the overlaps
              frequency_density[range_min] += range_density;
              frequency_density[range_max] -= range_density;
            }
            last_frequency = frequency;
          }
          double range_cumulation = 0;
          for (std::pair<double, double> const & frequency_diff : frequency_density) {
            // integrate
            range_cumulation += frequency_diff.second;
            frequency_diff.second = range_cumulation;
          }
          std::vector<double> result_probabilities(d_nbuckets, 0);
          {
            size_t result_bucket = 1;
            double frequency = calculated_result_histogram[result_bucket];
            if (result_bucket > 0) {
              double range_min = last_frequency;
              double range_max = frequency;
              if (frequency < last_frequency) {
                std::swap(range_min, range_max);
              }
              auto range_begin_end = frequency_density.upper_bound(range_min);
              auto range_begin_begin = range_begin_end;
              range_begin_begin --;
              auto range_end_end = frequency_density.upper_bound(range_max);
              auto range_end_begin = range_end_end;
              range_end_begin --;

              double begin_integral = (range_begin_end->second - range_begin_begin->second) * (range_min - range_begin_begin->first) / (range_begin_end->first - range_begin_begin->first) + range_begin_begin->second;
              double end_integral = (range_end_end->second - range_end_begin->second) * (range_min - range_end_begin->first) / (range_end_end->first - range_end_begin->first) + range_end_begin->second;

              double range_probability = (end_integral - begin_integral) / range_cumulation;
              double half_result_probability = range_probability / 2; // one for left bucket, one for right bucket
              result_probabilities[result_bucket - 1] += half_result_probability;
              result_probabilities[result_bucket] += half_result_probability;
            }
            last_frequency = frequency;
          }
          */

          unknown_proportion = measured_result_total;

          // in the qa_test, for output bucket 0 == output value -3,
          // the result total is zero because there is no density in the overlapping region
          // adding -3 to the addend makes it completely nonintersecting with the sum.
          // it only intersects when the scalar addend is >= 0

          // this might relate to meaning of calculated_result_total that is left out of the algorithm.
          // areas that overlap are much more valuable than those that don't.

          for (size_t result_bucket = 0; unknown_proportion != 0 && result_bucket < d_nbuckets; result_bucket ++) {
    
            double calculated_result_frequency = calculated_result_histogram[result_bucket];
            freq_type const & measured_result_frequency = in[result_input_idx()][result_bucket];
    
            double population_frequency = calculated_result_frequency;
            double population_count = calculated_result_total;
            double sample_frequency = measured_result_frequency;
            double sample_count = measured_result_total;
            
            if (population_frequency <= 0) {
                if (sample_frequency > 0) {
                    // no population here, but samples found.
                    // result is technically zero assuming accurate population.
                    // maybe population/sample could be inverted?
                }
                continue;
            } else if (population_frequency >= population_count) {
                // all population here, but not all samples here.
                // result is technically zero assuming accurate population
                // maybe population/sample could be inverted?
                unknown_proportion = sample_frequency;
                break;
            }
    
            double population_proportion = population_count > 0 ? population_frequency / population_count : 0;
            double sample_proportion = sample_count > 0 ? sample_frequency / sample_count : 0;

                // output value goes up as sample total goes down.
                //
                // sample total is proportional to sample mean, variance, and deviation^2
                // this scales the sample z scale by 1 / sqrt(sample total)
                // and the likelihood by sqrt(sample total)
                // and the final proportion by sqrt(sample total) ^ buckets
                // -> ther sample total that is contributing multiple times.
                //    the sample total is condensing the z scores, and that factor has an exponential
                //    relationship with the result, due to repeated multiplication
                //    the event count of the result drops many orders of magnitude below 1.0
    
            // binomial distribution, might be correct
            double theorised_sample_mean = sample_count * population_proportion;
            double theorised_sample_variance = theorised_sample_mean * (1 - population_proportion);
            double theorised_sample_deviation = sqrt(theorised_sample_variance);
            double sample_z = (sample_frequency - theorised_sample_mean) / theorised_sample_deviation;
    
            // likelihood of sampling a bucket like this, or a less-similar one.  hits 100% at precise similarity.
            // this is a bucket analogous to the result bucket
            // it might end up being helpful to show these histograms here, and to test them in qa.  they could be cached in a member variable.
            double sample_probability_given_unknown = erfc(-abs(sample_z) * M_SQRT1_2);

            // the formula is P(Unknown given Histogram) = P(Histogram given Unknown) * P(Unknown) / P(Histogram)
            //double unknown_probability_given_sample = sample_probability_given_unknown / d_result_independent_probability[result_bucket] / d_nbuckets;
    
            // old: the unknown likelihood is the likelihood of all the sample likelihoods happening, given the population
            unknown_proportion *= sample_probability_given_unknown; //unknown_probability_given_sample;
            // old: sample_count -= sample_frequency;
          }

          if (std::is_integral<freq_type>::value) {
            out[unknown_bucket] = round(unknown_proportion);
          } else {
            out[unknown_bucket] = unknown_proportion;
          }

        }
#endif
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

  size_t output_idx() const override
  {
    return d_output_idx;
  }

  size_t nbuckets() const override
  {
    return d_nbuckets;
  }

private:
  struct map_linear_interpolation
  {
    detail::square_ctd<double, sizeof...(Doubles) - 1> last_ndplane_values;
    detail::square_ctd<double, sizeof...(Doubles) - 1> this_ndplane_values;
    detail::square_ctd<size_t, sizeof...(Doubles) - 1> last_ndplane_buckets;
    detail::square_ctd<size_t, sizeof...(Doubles) - 1> this_ndplane_buckets;
    size_t input_buckets[sizeof...(Doubles)];
  };

  void update_coeffs()
  {
    double range = (double)d_max - d_min;

    // bucket values are slid by half a bucket to represent the middle of the bucket
    // rather than the start
    d_bucket_to_value_coeff = range / d_nbuckets;
    d_bucket_to_value_offset = d_min + /*slide up*/d_bucket_to_value_coeff / 2;
    d_value_to_bucket_coeff = d_nbuckets / range;
    d_value_to_bucket_offset = -d_value_to_bucket_coeff * d_min - /*slide down*/0.5;

    d_unknown_result_map.clear();
    d_unknown_result_map.resize(d_nbuckets, std::vector<sparse_dot>{d_nbuckets});

    map_linear_interpolation interpolation;
    size_t * bucket_ptr = interpolation.input_buckets;
    size_t & outermost_bucket = *bucket_ptr;
    interpolation.last_ndplane_values.resize(d_nbuckets);
    interpolation.last_ndplane_buckets.resize(d_nbuckets);
    interpolation.this_ndplane_values.resize(d_nbuckets);
    interpolation.this_ndplane_buckets.resize(d_nbuckets);
    for (outermost_bucket = 0; outermost_bucket <= d_nbuckets; outermost_bucket ++) {
      update_expr_map(interpolation, bucket_ptr, outermost_bucket > 0, bucket_to_value(outermost_bucket));
      interpolation.this_ndplane_buckets.swap(interpolation.last_ndplane_buckets);
      interpolation.this_ndplane_values.swap(interpolation.last_ndplane_values);
    }
  }

  // recursive template function walks each variable, trying all combinations to see what result they produce
  // this function has so many arguments that i'd like it if this algorithm were eventually bundled into a class with member variables
  template <typename... Params>
  void update_expr_map(map_linear_interpolation & interpolation, size_t *bucket_ptr, bool do_interpolation, Params... params)
  {
    bucket_ptr ++;
    size_t & bucket = *bucket_ptr;
    for (bucket = 0; bucket <= d_nbuckets; bucket ++) {
      update_expr_map(interpolation, bucket_ptr, do_interpolation && bucket > 0, bucket_to_value(bucket), params...);
    }
  }

  // WIP maybe? this is just sampling midpoints, but it could maybe do sub-bucket ranges by outputing weights for each sequence
  void update_expr_map(map_linear_interpolation & interpolation, size_t *bucket_ptr, bool do_interpolation, Doubles... params)
  {
    double value = d_expr(params...);
    interpolation.this_ndplane_values[interpolation.input_buckets + 1] = value;
    interpolation.this_ndplane_buckets[interpolation.input_buckets + 1] = value_to_bucket(value);

    if (! do_interpolation) {
      return;
    }

    size_t prev_buckets[sizeof...(Doubles)];
    for (size_t w = 0; w < sizeof...(Doubles); w ++) {
      // this could have been tracked next to bucket_ptr to not loop over the coords every call here
      prev_buckets[w] = interpolation.input_buckets[w] - 1;
    }

    // todo: update with the densities implied by change since the preceding ndplane

    // prev_buckets and buckets describe a volume of output values, with the first axis held by last_ndplane vs this_ndplane.
    // - the goal is to stretch them over the result buckets with appropriate linearly interpolated densities, like the second
    //   chunk of commented code misapplied farther up this file.  the fraction of the volume that covers a bucket, is its weight.
    //   this will likely simply mean iterating over every output bucket within the ranges, just like this outer recursive iteration.

    size_t result = interpolation.last_ndplane_buckets[prev_buckets + 1];

    if (result >= 0 && result < d_nbuckets) {
      if (is_reverse_solving()) {
        d_unknown_result_map[prev_buckets[d_param_output_idx]][result].add_sequence(prev_buckets, d_param_output_idx);
      } else {
        d_unknown_result_map[0][result].add_sequence(prev_buckets);
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
  // all variables are in order.  the output is at d_output_idx.  the result is at result_idx().

  // input indices index the block inputs and have no output in them.  they are one fewer, and the output is skipped.
  size_t var_idx_to_input_idx(size_t var_idx) const
  {
    return var_idx < d_output_idx ? var_idx : var_idx - 1;
  }

  // parameter indices index the expression parameters.  they are one fewer, and the result is skipped.
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

  // having the result index be provided by this member function
  // provides some space towards letting the user set their own variable/input order in the UI
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
histogram_solve<freq_type, Doubles...>::make(double min, double max, std::function<double(Doubles...)> expr, size_t output_idx, size_t nbuckets)
{
  return gnuradio::make_block_sptr<histogram_solve_impl<freq_type, Doubles...>>(
    min, max, expr, output_idx, nbuckets);
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

