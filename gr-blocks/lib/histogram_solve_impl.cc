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
  std::vector<sparse_dot> d_output_map;

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
      d_output_idx(output_idx > 0 ? output_idx - 1 : sizeof...(Doubles)),
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

    for (size_t item = 0; item < noutput_items; item ++) {
      double total = 0;
      for (size_t unknown_bucket = 0; unknown_bucket < d_nbuckets; unknown_bucket ++) {
        double unknown = unknown_bucket * d_bucket_to_value_coeff + d_bucket_to_value_offset;
        freq_type & sample = out[unknown_bucket];
        sample = d_output_map[unknown_bucket].calc(in, item * d_nbuckets);
        total += sample;
      }
      if (std::is_floating_point<freq_type>::value) {
          total = 1.0 / total;
          for (size_t unknown_bucket = 0; unknown_bucket < d_nbuckets; unknown_bucket ++) {
            out[unknown_bucket] *= total;
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
    return d_output_idx < sizeof...(Doubles) ? d_output_idx + 1 : 0;
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

    d_output_map.clear();
    d_output_map.resize(d_nbuckets);
    size_t buckets[sizeof...(Doubles) + 1];
    update_expr_map(buckets, buckets);
  }

  // LATER? this is just sampling midpoints, but it could maybe do sub-bucket ranges by outputing weights for each sequence
  void update_expr_map(size_t *bucket_head, size_t *bucket_ptr, Doubles... params)
  {
    size_t & result = *bucket_ptr;
    result = d_expr(params...) * d_value_to_bucket_coeff + d_value_to_bucket_offset;

    if (result >= 0 && result < d_nbuckets) {
      d_output_map[bucket_head[d_output_idx]].add_sequence(bucket_head, d_output_idx);
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

