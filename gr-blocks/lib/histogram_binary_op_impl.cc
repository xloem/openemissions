/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gnuradio/io_signature.h>

#include <openemissions/histogram_binary_op.h>

namespace gr {
namespace openemissions {

namespace detail {

// sparse_dot is used to calculate each output value and could likely be optimised
class sparse_dot {
public:
    std::vector<size_t> index_pairs;
    
    void clear() {
        index_pairs.clear();
    }

    void add_pair(size_t idx1, size_t idx2) {
        index_pairs.push_back(idx1);
        index_pairs.push_back(idx2);
    }

    template <typename T>
    T calc(const T * left, const T * right) {
        T result = 0;
        for (size_t pair = 0; pair < index_pairs.size(); pair += 2) {
            result += left[index_pairs[pair]] * right[index_pairs[pair+1]];
        }
        return result;
    }
};
}

template <typename freq_type>
class histogram_binary_op_impl : public histogram_binary_op<freq_type>
{
private:
  double d_min;
  double d_max;
  std::function<double(double, double)> d_op;
  size_t d_nbuckets;

  double d_bucket_to_value_coeff;
  double d_bucket_to_value_offset;
  double d_value_to_bucket_coeff;
  double d_value_to_bucket_offset;

  // each output item is the sum of the products of pairs of input items,
  // different for each output item
  std::vector<detail::sparse_dot> d_output_map;

  // items that leave the bounds of the histograms are collected as a single quantity
  detail::sparse_dot d_overflow_map;

  bool d_warned;

public:
  /*
   * The private constructor
   */
  histogram_binary_op_impl(double min, double max, std::function<double(double,double)> op, size_t nbuckets)
  : gr::sync_block("histogram_binary_op",
                   gr::io_signature::make(2, 2, sizeof(freq_type) * nbuckets),
                   gr::io_signature::make(1, 1, sizeof(freq_type) * nbuckets)),
    d_min(min),
    d_max(max),
    d_op(op),
    d_nbuckets(nbuckets)
  {
    update_coeffs();
  }
  
  /*
   * Our virtual destructor.
   */
  ~histogram_binary_op_impl()
  {
  }

  bool start() override
  {
    d_warned = false;
    return sync_block::start();
  }

  // Where all the action really happens
  int
  work(int noutput_items,
       gr_vector_const_void_star &input_items,
       gr_vector_void_star &output_items)
  {
    const freq_type *in1 = reinterpret_cast<const freq_type*>(input_items[0]);
    const freq_type *in2 = reinterpret_cast<const freq_type*>(input_items[1]);
    freq_type *out = reinterpret_cast<freq_type*>(output_items[0]);
  
    for (size_t item = 0; item < noutput_items; item ++, in1 += d_nbuckets, in2 += d_nbuckets, out += d_nbuckets) {
      freq_type total = 0;
      freq_type inclusion = 0;

      freq_type overflow = d_overflow_map.calc(in1, in2);
      freq_type theoretical_total = 1;

      if (!std::is_floating_point<freq_type>::value) {
        freq_type left_total = 0, right_total = 0;
        for (size_t w = 0; w < d_nbuckets; w ++) {
          left_total += in1[w];
          right_total += in2[w];
        }
        theoretical_total = left_total * right_total;
      }

      for (size_t bucket = 0; bucket < d_nbuckets; bucket ++) {
        freq_type proportion = d_output_map[bucket].calc(in1, in2);
        // assuming overflow, rescale to theoretical total
        out[bucket] = proportion * theoretical_total / (theoretical_total - overflow);
      }
      if (! d_warned && overflow) {
        GR_LOG_WARN(this->d_logger, "possible operation result " + std::to_string(overflow) + " outside histogram range [" + std::to_string(d_min) + ", " + std::to_string(d_max) + ")");
        GR_LOG_WARN(this->d_logger, "Further range warnings from this block suppressed.");
        d_warned = true;
      }
    }
  
    // Tell runtime system how many output items we produced.
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

  void set_op(const std::function<double(double,double)> &op) override
  {
    d_op = op;
  }

  const std::function<double(double,double)> & op() const override
  {
    return d_op;
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

    // operations are performed in advance by making bucket mappings
    d_output_map.clear();
    d_output_map.resize(d_nbuckets);
    d_overflow_map.clear();

    for (size_t bucket1 = 0; bucket1 < d_nbuckets; bucket1 ++) {
      double value1 = bucket1 * d_bucket_to_value_coeff + d_bucket_to_value_offset;
      for (size_t bucket2 = 0; bucket2 < d_nbuckets; bucket2 ++) {
        double value2 = bucket2 * d_bucket_to_value_coeff + d_bucket_to_value_offset;

        double result = d_op(value1,  value2);

        ssize_t bucket3 = result * d_value_to_bucket_coeff + d_value_to_bucket_offset;
        if (bucket3 >= 0 && bucket3 < d_nbuckets) {
          d_output_map[bucket3].add_pair(bucket1, bucket2);
        } else {
          d_overflow_map.add_pair(bucket1, bucket2);
          d_warned = false;
        }
      }
    }
  }
};

template <typename T>
typename histogram_binary_op<T>::sptr
histogram_binary_op<T>::make(T min, T max, const std::function<double(double,double)> & op, size_t nbuckets)
{
  return gnuradio::make_block_sptr<histogram_binary_op_impl<T>>(
    min, max, op, nbuckets);
}

template class histogram_binary_op<double>;
template class histogram_binary_op<float>;
template class histogram_binary_op<uint64_t>;

} /* namespace openemissions */
} /* namespace gr */

