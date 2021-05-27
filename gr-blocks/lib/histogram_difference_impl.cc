/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gnuradio/io_signature.h>

#include <openemissions/histogram_difference.h>

namespace gr {
namespace openemissions {

template <typename freq_type>
class histogram_difference_impl : public histogram_difference<freq_type>
{
private:
  double d_min;
  double d_max;
  size_t d_nbuckets;

  double d_bucket_to_value_coeff;
  double d_bucket_to_value_offset;
  double d_value_to_bucket_coeff;
  double d_value_to_bucket_offset;

  bool d_warned;

public:
  /*
   * The private constructor
   */
  histogram_difference_impl(double min, double max, size_t nbuckets)
  : gr::sync_block("histogram_difference",
                   gr::io_signature::make(2, 2, sizeof(freq_type) * nbuckets),
                   gr::io_signature::make(1, 1, sizeof(freq_type) * nbuckets)),
    d_min(min),
    d_max(max),
    d_nbuckets(nbuckets)
  {
    update_coeffs();
  }
  
  /*
   * Our virtual destructor.
   */
  ~histogram_difference_impl()
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
      memset(out, 0, d_nbuckets * sizeof(freq_type));

      for (size_t bucket1 = 0; bucket1 < d_nbuckets; bucket1 ++) {
        double value1 = bucket1 * d_bucket_to_value_coeff + d_bucket_to_value_offset;
        for (size_t bucket2 = 0; bucket2 < d_nbuckets; bucket2 ++) {
          double value2 = bucket2 * d_bucket_to_value_coeff + d_bucket_to_value_offset;

          freq_type proportion = in1[bucket1] * in2[bucket2];
          total += proportion;

          double difference = value1 - value2;
          ssize_t bucket3 = difference * d_value_to_bucket_coeff + d_value_to_bucket_offset;
          if (bucket3 >= 0 && bucket3 < d_nbuckets) {
            out[bucket3] += proportion;
            inclusion += proportion;
          } else if (! d_warned) {
            GR_LOG_WARN(this->d_logger, "possible difference " + std::to_string(difference) + " outside histogram range [" + std::to_string(d_min) + ", " + std::to_string(d_max) + ")");
            GR_LOG_WARN(this->d_logger, "Further range warnings from this block suppressed.");
            d_warned = true;
          }
        }
      }

      if (inclusion != total) {
        inclusion = total / inclusion;
        for (size_t bucket = 0; bucket < d_nbuckets; bucket ++) {
          out[bucket] *= inclusion;
        }
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

  size_t nbuckets() const override
  {
    return d_nbuckets;
  }

private:
  void update_coeffs()
  {
        double range = (double)d_max - d_min;
        d_bucket_to_value_coeff = range / d_nbuckets;
        d_bucket_to_value_offset = d_min * d_bucket_to_value_coeff;
        d_value_to_bucket_coeff = d_nbuckets / range;
        d_value_to_bucket_offset = d_nbuckets * -d_min / range;
  }
};

template <typename T>
typename histogram_difference<T>::sptr
histogram_difference<T>::make(T min, T max, size_t nbuckets)
{
  return gnuradio::make_block_sptr<histogram_difference_impl<T>>(
    min, max, nbuckets);
}

template class histogram_difference<double>;
template class histogram_difference<float>;
template class histogram_difference<uint64_t>;

} /* namespace openemissions */
} /* namespace gr */

