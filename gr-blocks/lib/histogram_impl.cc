/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gnuradio/io_signature.h>
#include <pmt/pmt.h>

#include <openemissions/histogram.h>

#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

static size_t pmt_hash(pmt::pmt_t const &val)
{
  if (val->is_symbol() || val->is_null()) {
    return std::hash<pmt::pmt_t>{}(val);
  } else if (val->is_bool()) {
    return std::hash<bool>{}(val == pmt::PMT_T);
  } else if (val->is_pair()) {
    return pmt_hash(pmt::car(val)) ^ pmt_hash(pmt::cdr(val));
  } else if (val->is_integer()) {
    return std::hash<long>{}(pmt::to_long(val));
  } else if (val->is_real()) {
    return std::hash<double>{}(pmt::to_double(val));
  } else if (val->is_complex()) {
    const std::complex<double> cmplx = pmt::to_complex(val);
    return std::hash<double>{}(cmplx.real()) ^ std::hash<double>{}(cmplx.imag());
  } else if (val->is_vector()) {
    size_t len = pmt::length(val);
    size_t result = 0;
    for (size_t i = 0; i < len; i ++) {
      result ^= pmt_hash(pmt::vector_ref(val, i));
    }
    return result;
  }
  throw pmt::notimplemented("hash pmt (?)", val);
}

namespace std {
template <>
class hash<std::vector<pmt::pmt_t>>
{
public:
  size_t operator()(std::vector<pmt::pmt_t> const &val) const
  {
    size_t result = val.size();
    for (size_t w = 0; w < val.size(); w ++) {
      result ^= pmt_hash(val[w]);
    }
    return result;
  }
};
}

namespace gr {
namespace openemissions {

template <typename input_type, typename freq_type>
class histogram_impl : public histogram<input_type, freq_type>
{
private:
  input_type d_min;
  input_type d_max;
  size_t d_nbuckets;
  const size_t d_vinlen;

  struct histogram_t {
    std::vector<uint64_t> d_overall;
    uint64_t d_total;
  };

  std::unordered_map<pmt::pmt_t, size_t> d_prop_tag_idxs;
  using prop_iterator = typename std::unordered_map<pmt::pmt_t, size_t>::iterator;
  std::vector<pmt::pmt_t> d_prop_tags;

  std::unordered_map<std::vector<pmt::pmt_t>, histogram_t> d_histograms;
  using hist_iterator = typename std::unordered_map<std::vector<pmt::pmt_t>, histogram_t>::iterator;
  histogram_t * d_histogram;

  double d_value_to_bucket_coeff;
  double d_value_to_bucket_offset;
  std::vector<gr::tag_t> d_work_tags;

public:
  /*
   * The private constructor
   */
  histogram_impl(input_type min, input_type max, size_t nbuckets, size_t vinlen, const std::vector<std::string> & prop_tag_keys, const std::string & filename)
  : gr::sync_block("histogram",
                   gr::io_signature::make(1, 1, sizeof(input_type) * vinlen),
                   gr::io_signature::make(1, 1, sizeof(freq_type) * nbuckets)),
    d_min(min),
    d_max(max),
    d_nbuckets(nbuckets),
    d_vinlen(vinlen),
    d_prop_tags(prop_tag_keys.size(), pmt::PMT_NIL),
    d_histogram(nullptr)
  {
    update_coeff();
    for (size_t w = 0; w < prop_tag_keys.size(); w ++) {
      d_prop_tag_idxs[pmt::string_to_symbol(prop_tag_keys[w])] = w;
    }
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
    std::vector<freq_type *> outs(
      reinterpret_cast<freq_type **>(output_items.data() + 1),
      reinterpret_cast<freq_type **>(output_items.data() + output_items.size())
    );
    bool warned = false;

    this->get_tags_in_window(d_work_tags, 0, 0, noutput_items);
    std::vector<gr::tag_t>::iterator tag_it = d_work_tags.begin();

    uint64_t abs_offset = this->nitems_read(0);
    int sample_tail = noutput_items * d_vinlen;

    for (int subsample = 0; subsample < sample_tail; abs_offset ++)
    {
      while (tag_it != d_work_tags.end() && tag_it->offset <= abs_offset) {
        gr::tag_t & tag = *tag_it;
        prop_iterator prop_it = d_prop_tag_idxs.find(tag.key);
        if (prop_it != d_prop_tag_idxs.end()) {
          // found a property tag
          size_t idx = prop_it->second;
          d_prop_tags[idx] = tag.value;
          d_histogram = nullptr;
        }
        ++ tag_it;
      }

      if (nullptr == d_histogram) {
        d_histogram = &d_histograms[d_prop_tags];
        if (d_histogram->d_overall.empty()) {
          d_histogram->d_overall.reserve(d_nbuckets);
          d_histogram->d_overall.resize(d_nbuckets);
          d_histogram->d_total = 0;
        }
      }

      for (int subtail = subsample + d_vinlen; subsample < subtail; subsample ++) {
        input_type value = in[subsample];
        if (value >= d_min && value < d_max) {
          size_t bucket = value * d_value_to_bucket_coeff + d_value_to_bucket_offset;
          d_histogram->d_overall[bucket] ++;
          d_histogram->d_total ++;
        } else if (!warned) {
          GR_LOG_WARN(this->d_logger, "value " + std::to_string(value) + " outside histogram range [" + std::to_string(d_min) + ", " + std::to_string(d_max) + ")");
          warned = true;
        }
      }
      
      if (std::is_floating_point<freq_type>::value) {
        for (size_t bucket = 0; bucket < d_nbuckets; bucket ++) {
          *out = freq_type(d_histogram->d_overall[bucket]) / d_histogram->d_total;
          out ++;
        } 
      } else {
        memcpy(out, d_histogram->d_overall.data(), sizeof(freq_type) * d_nbuckets);
        out += d_nbuckets;
      }
    }

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
    return d_nbuckets;
  }

  uint64_t total() const override
  {
    if (nullptr == d_histogram) {
        return 0;
    } else {
        return d_histogram->d_total;
    }
  }

private:
  inline void update_coeff()
  {
    input_type range = d_max - d_min;
    d_value_to_bucket_coeff = d_nbuckets / range;
    d_value_to_bucket_offset = d_nbuckets * -d_min / range;
  }
};

template <typename input_type, typename freq_type>
typename histogram<input_type, freq_type>::sptr
histogram<input_type, freq_type>::make(input_type min, input_type max, size_t nbuckets, size_t vinlen, const std::vector<std::string> & prop_tag_keys, const std::string & filename)
{
  return gnuradio::make_block_sptr<histogram_impl<input_type, freq_type>>(min, max, nbuckets, vinlen, prop_tag_keys, filename);
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
