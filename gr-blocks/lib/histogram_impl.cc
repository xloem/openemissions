/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gnuradio/io_signature.h>
#include <openemissions/histogram.h>

#include <type_traits>

namespace gr {
namespace openemissions {

template <typename input_type, typename freq_type>
class histogram_impl : public histogram<input_type, freq_type>
{
  private:
  input_type d_min;
  input_type d_max;
  std::vector<uint64_t> d_histogram;
  const size_t d_vinlen;
  double d_value_to_bucket_coeff;
  double d_value_to_bucket_offset;
  uint64_t d_total;

  public:
  /*
   * The private constructor
   */
  histogram_impl(input_type min, input_type max, size_t nbuckets, size_t vinlen)
  : gr::sync_block("histogram",
                   gr::io_signature::make(1, 1, sizeof(input_type) * vinlen),
                   gr::io_signature::make(1, 1, sizeof(freq_type) * nbuckets)),
    d_min(min),
    d_max(max),
    d_histogram(nbuckets, 0),
    d_vinlen(vinlen),
    d_total(0)
  {
    update_coeff();
  }

  /*
   * Our virtual destructor.
   */
  ~histogram_impl()
  { }

  // Where all the action really happens
  int work(int noutput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items)
  {
    const input_type *in = reinterpret_cast<const input_type *>(input_items[0]);
    freq_type *out = reinterpret_cast<freq_type *>(output_items[0]);
    bool warned = false;

    for (int sample = 0; sample < noutput_items * d_vinlen; sample += d_vinlen)
    {
      for (int subsample = sample; subsample < sample + d_vinlen; subsample ++) {
        input_type value = in[subsample];
        if (value >= d_min && value < d_max) {
          size_t bucket = value * d_value_to_bucket_coeff + d_value_to_bucket_offset;
          d_histogram[bucket] ++;
          d_total ++;
        } else if (!warned) {
          GR_LOG_WARN(this->d_logger, "value " + std::to_string(value) + " outside histogram range [" + std::to_string(d_min) + ", " + std::to_string(d_max) + ")");
          warned = true;
        }
      }
      
      if (std::is_floating_point<freq_type>::value) {
        for (size_t bucket = 0; bucket < d_histogram.size(); bucket ++) {
          *out = freq_type(d_histogram[bucket]) / d_total;
          *out ++;
        } 
      } else {
        memcpy(out, d_histogram.data(), sizeof(freq_type) * d_histogram.size());
        out += d_histogram.size();
      }
    }
  
    // Tell runtime system how many output items we produced.
    return noutput_items;
  }

  void set_min(input_type min) override
  {
    d_min = min;
    update_coeff();
  }

  input_type min() const override
  {
    return d_min;
  }

  void set_max(input_type max) override
  {
    d_max = max;
    update_coeff();
  }

  input_type max() const override
  {
    return d_max;
  }

  size_t nbuckets() const override
  {
    return d_histogram.size();
  }

  uint64_t total() const override
  {
    return d_total;
  }

private:
  inline void update_coeff()
  {
    input_type range = d_max - d_min;
    double nbuckets = d_histogram.size();
    d_value_to_bucket_coeff = nbuckets / range;
    d_value_to_bucket_offset = -nbuckets * d_min / range;
  }
};

template <typename input_type, typename freq_type>
typename histogram<input_type, freq_type>::sptr
histogram<input_type, freq_type>::make(input_type min, input_type max, size_t nbuckets, size_t vinlen)
{
  return gnuradio::make_block_sptr<histogram_impl<input_type, freq_type>>(min, max, nbuckets, vinlen);
}

template class histogram<double, double>;
template class histogram<float, double>;
template class histogram<int64_t, double>;
template class histogram<int32_t, double>;
template class histogram<int16_t, double>;
template class histogram<int8_t, double>;
template class histogram<uint64_t, double>;
template class histogram<uint32_t, double>;
template class histogram<uint16_t, double>;
template class histogram<uint8_t, double>;

template class histogram<double, float>;
template class histogram<float, float>;
template class histogram<int64_t, float>;
template class histogram<int32_t, float>;
template class histogram<int16_t, float>;
template class histogram<int8_t, float>;
template class histogram<uint64_t, float>;
template class histogram<uint32_t, float>;
template class histogram<uint16_t, float>;
template class histogram<uint8_t, float>;

template class histogram<double, uint64_t>;
template class histogram<float, uint64_t>;
template class histogram<int64_t, uint64_t>;
template class histogram<int32_t, uint64_t>;
template class histogram<int16_t, uint64_t>;
template class histogram<int8_t, uint64_t>;
template class histogram<uint64_t, uint64_t>;
template class histogram<uint32_t, uint64_t>;
template class histogram<uint16_t, uint64_t>;
template class histogram<uint8_t, uint64_t>;

} /* namespace openemissions */
} /* namespace gr */
