/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gnuradio/io_signature.h>

#include <openemissions/histogram_solve.h>

namespace gr {
namespace openemissions {

/*
 * notes on phone.
 * forward and reverse solve are different cases.
 * for both, we output equation result.
 * for reverse solve, we use normal distribution to make distribution of likelihood of sampling from one result distribution to get the other
 */

template <typename freq_type, typename... Doubles>
class histogram_solve_impl : public histogram_solve<freq_type, Doubles...>
{
private:
  double d_min;
  double d_max;
  std::function<double(Doubles...)> d_expr;
  size_t d_output_idx;
  size_t d_nbuckets;

  double d_bucket_to_value_coeff;
  double d_bucket_to_value_offset;
  double d_value_to_bucket_coeff;
  double d_value_to_bucket_offset;

  // sparse_dot is used to calculate each output value and could likely be optimised
  class sparse_dot {
  public:
      std::vector<size_t> index_sequences;
      
      void clear() {
          index_sequences.clear();
      }
  
      void add_sequence(size_t *sequence, size_t output_idx) {
          for (size_t i = 0; i < sizeof...(Doubles) + 1; i ++) {
            if (i != output_idx) {
              index_sequences.push_back(sequence[i]);
            }
          }
      }
  
      template <typename T>
      T calc(const T ** knowns, size_t offset) {
          T result = 0;
          for (size_t sequence = 0; sequence < index_sequences.size();) {
              T product = 1;
              for (size_t coeff = 0; coeff < sizeof...(Doubles); coeff ++, sequence ++) {
                product *= knowns[coeff][offset + index_sequences[sequence]];
              }
              result += product;
          }
          return result;
      }
  };

  // each output item is the sum of the products of sequences of input items,
  // different for each output item
  std::vector<sparse_dot> d_result_map;

  // this breaks the result map out into the components made by each possible unknown value
  std::vector<std::vector<sparse_dot>> d_unknown_result_map;

public:
  /*
   * The private constructor
   */
  histogram_solve_impl(double min, double max, std::function<double(Doubles...)> expr, size_t output_idx, size_t nbuckets)
    : gr::sync_block("histogram_solve",
            gr::io_signature::make(sizeof...(Doubles), sizeof...(Doubles), sizeof(freq_type) * nbuckets),
            gr::io_signature::make(1, 1, sizeof(freq_type) * nbuckets)),
      d_min(min),
      d_max(max),
      d_expr(expr),
      d_output_idx(output_idx),
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

  int
  work(int noutput_items,
       gr_vector_const_void_star &input_items,
       gr_vector_void_star &output_items)
  {
    const freq_type **in = reinterpret_cast<const freq_type**>(input_items.data());
    freq_type *out = reinterpret_cast<freq_type*>(output_items[0]);
    //memset(out, 0, sizeof(freq_type) * noutput_items);

    std::vector<double> hist_cache(d_nbuckets);

    for (size_t item = 0; item < noutput_items; item ++) {

      if (d_output_idx == result_idx()) {
        // forward solving
        for (size_t unknown_bucket = 0; unknown_bucket < d_nbuckets; unknown_bucket ++) {
          freq_type & unknown_proportion = out[unknown_bucket];

          double unknown = unknown_bucket * d_bucket_to_value_coeff + d_bucket_to_value_offset;
          unknown_proportion = d_result_map[unknown_bucket].calc(in, item * d_nbuckets);
        }
      } else {
        // reverse solving
        double measured_result_total = 0;
        for (size_t result_bucket = 0; result_bucket < d_nbuckets; result_bucket ++) {
          measured_result_total += in[result_input_idx()][result_bucket];
        }

        // we form distributions for each potential value of the unknown variable
        for (size_t unknown_bucket = 0; unknown_bucket < d_nbuckets; unknown_bucket ++) {
          double unknown_proportion;

          // given the potential value, we form the histogram of the result variable
          // this histogram is placed into hist_cache

          double calculated_result_total = 0;
          for (size_t result_bucket = 0; result_bucket < d_nbuckets; result_bucket ++) {
            double frequency = d_unknown_result_map[unknown_bucket][result_bucket].calc(in, item * d_nbuckets);
            calculated_result_total += frequency;
            hist_cache[result_bucket] = frequency;
          }

          unknown_proportion = calculated_result_total;

          // in the qa_test, for output bucket 0 == output value -3,
          // the result total is zero because there is no density in the overlapping region
          // adding -3 to the addend makes it completely nonintersecting with the sum.
          // it only intersects when the scalar addend is >= 0

          // this might relate to meaning of calculated_result_total that is left out of the algorithm.
          // areas that overlap are much more valuable than those that don't.

          for (size_t result_bucket = 0; unknown_proportion != 0 && result_bucket < d_nbuckets; result_bucket ++) {
            // calculate potential result bucket and compare with reality
    
            double calculated_result_frequency = hist_cache[result_bucket];
            freq_type const & measured_result_frequency = in[result_input_idx()][result_bucket];
    
            double population_frequency = calculated_result_frequency;
            double population_count = calculated_result_total;
            double sample_frequency = measured_result_frequency;
            double sample_count = measured_result_total;
    
            // i'm thinking about the meaning of population here,
            // in the context of having 0 measured population.
            // can we be clear on what the distributions are of, here?
            //  we vor
            
            if (population_frequency <= 0) {
                if (sample_frequency > 0) {
                    unknown_proportion = 0;
                }
                continue;
            } else if (population_frequency >= population_count) {
                if (sample_frequency < sample_count) {
                    unknown_proportion = 0;
                }
                continue;
            }
    
            double population_proportion = population_count > 0 ? population_frequency / population_count : 0;
            double sample_proportion = sample_count > 0 ? sample_frequency / sample_count : 0;
    
            // binomial distribution, might be correct
            double theorised_sample_mean = (sample_count + 1) * population_proportion;
            double theorised_sample_variance = theorised_sample_mean * (1 - population_proportion);
            double theorised_sample_deviation = sqrt(theorised_sample_variance);
            double sample_upper_z = (sample_frequency + 1 - theorised_sample_mean) / theorised_sample_deviation;
            double sample_lower_z = (sample_frequency + 0 - theorised_sample_mean) / theorised_sample_deviation;
    
            // likelihood of this result bucket happening, given this unknown value
            // this is a bucket analogous to the result bucket
            // it might end up being helpful to show these histograms here, and to test them in qa.  they could be cached in a member variable.
            double sample_likelihood = 0.5 * (erfc(-sample_upper_z * M_SQRT1_2) - erfc(-sample_lower_z * M_SQRT1_2));
    
            // the unknown likelihood is the likelihood of all the sample likelihoods happening, given the population
            unknown_proportion *= sample_likelihood;
          }

          out[unknown_bucket] = unknown_proportion;
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

  size_t output_idx() const override
  {
    return d_output_idx;
  }

  size_t nbuckets() const override
  {
    return d_nbuckets;
  }

private:
  void update_coeffs()
  {
    double range = (double)d_max - d_min;

    // bucket values are slid by half a bucket to represent the middle of the bucket
    // rather than the start
    d_bucket_to_value_coeff = range / d_nbuckets;
    d_bucket_to_value_offset = d_min + /*slide up*/d_bucket_to_value_coeff / 2;
    d_value_to_bucket_coeff = d_nbuckets / range;
    d_value_to_bucket_offset = -d_value_to_bucket_coeff * d_min - /*slide down*/0.5;

    d_result_map.clear();
    d_result_map.resize(d_nbuckets);
    d_unknown_result_map.clear();
    d_unknown_result_map.resize(d_nbuckets, std::vector<sparse_dot>{d_nbuckets});
    size_t buckets[sizeof...(Doubles) + 1];
    update_expr_map(buckets, buckets + first_parameter_idx());
  }

  // LATER maybe? this is just sampling midpoints, but it could maybe do sub-bucket ranges by outputing weights for each sequence
  void update_expr_map(size_t *bucket_head, size_t *bucket_ptr, Doubles... params)
  {
    size_t & result = bucket_head[result_idx()];
    result = d_expr(params...) * d_value_to_bucket_coeff + d_value_to_bucket_offset;

    if (result >= 0 && result < d_nbuckets) {
      d_result_map[result].add_sequence(bucket_head, result_idx());
      if (d_output_idx != result_idx()) {
        d_unknown_result_map[bucket_head[d_output_idx]][result].add_sequence(bucket_head, d_output_idx);
      }
    }
  }
  // recursive template function walks each variable, trying all combinations
  template <typename... Params>
  void update_expr_map(size_t *bucket_head, size_t *bucket_ptr, Params... params)
  {
    size_t & bucket = *bucket_ptr;
    for (bucket = 0; bucket < d_nbuckets; bucket ++) {
      double value = bucket * d_bucket_to_value_coeff + d_bucket_to_value_offset;
      update_expr_map(bucket_head, bucket_ptr + 1, value, params...);
    }
  }

  size_t var_idx_to_input_idx(size_t var_idx) const
  {
    return var_idx < d_output_idx ? var_idx : var_idx - 1;
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

  static constexpr size_t result_idx()
  {
    return 0;
  }

  static constexpr size_t first_parameter_idx()
  {
    return 1;
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

