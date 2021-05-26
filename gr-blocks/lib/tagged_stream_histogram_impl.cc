/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gnuradio/io_signature.h>
#include <pmt/pmt.h>

#include <openemissions/tagged_stream_histogram.h>

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
class tagged_stream_histogram_impl : public tagged_stream_histogram<input_type, freq_type>
{
private:
  input_type d_min;
  input_type d_max;
  size_t d_nbuckets;
  const size_t d_vinlen;

  struct histogram_t {
    std::vector<uint64_t> d_histograms;
    std::vector<uint64_t> d_totals;
    using iterator = typename std::vector<uint64_t>::iterator;
  };

  std::unordered_map<pmt::pmt_t, size_t> d_prop_tag_idxs;
  using prop_iterator = typename std::unordered_map<pmt::pmt_t, size_t>::iterator;
  std::vector<pmt::pmt_t> d_prop_tags;

  std::unordered_map<std::vector<pmt::pmt_t>, histogram_t> d_histograms;
  using hist_iterator = typename std::unordered_map<std::vector<pmt::pmt_t>, histogram_t>::iterator;
  histogram_t * d_histogram;
  size_t d_nitems;

  double d_value_to_bucket_coeff;
  double d_value_to_bucket_offset;
  std::vector<std::vector<gr::tag_t>> d_work_tags;
  std::vector<std::vector<gr::tag_t>::iterator> d_work_tag_its;
  std::vector<int64_t> d_work_input_offsets;

public:
  /*
   * The private constructor
   */
  tagged_stream_histogram_impl(input_type min, input_type max, size_t nbuckets, size_t vinlen, const std::vector<std::string> & prop_tag_keys, const std::string & len_tag_key, const std::string & filename)
  : gr::tagged_stream_block("tagged_stream_histogram",
                            gr::io_signature::make(1, 2, sizeof(input_type) * vinlen),
                            gr::io_signature::make(1, 2, sizeof(freq_type) * nbuckets),
                            len_tag_key),
    d_min(min),
    d_max(max),
    d_nbuckets(nbuckets),
    d_vinlen(vinlen),
    d_prop_tags(prop_tag_keys.size(), pmt::PMT_NIL),
    d_histogram(nullptr),
    d_nitems(0)
  {
    update_coeff();
    for (size_t w = 0; w < prop_tag_keys.size(); w ++) {
      d_prop_tag_idxs[pmt::string_to_symbol(prop_tag_keys[w])] = w;
    }
  }

  /*
   * Our virtual destructor.
   */
  ~tagged_stream_histogram_impl()
  { }

  // Where all the action really happens
  int work(int noutput_items,
           gr_vector_int& ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items)
  {
    const input_type *in = reinterpret_cast<const input_type *>(input_items[0]);
    const input_type *in2 = input_items.size() > 1 ? reinterpret_cast<const input_type *>(input_items[1]) : nullptr;
    size_t nitems_min = in2 ? std::min(ninput_items[0], ninput_items[1]) : ninput_items[0];
    size_t nitems_max = in2 ? std::max(ninput_items[0], ninput_items[1]) : ninput_items[0];
    freq_type *out = reinterpret_cast<freq_type *>(output_items[0]);
    freq_type *out2 = output_items.size() > 1 ? reinterpret_cast<freq_type *>(output_items[1]) : nullptr;
    bool warned = false;

    // grab the tags, and prepare iterators to walk them in parallel
    d_work_tags.resize(input_items.size());
    d_work_tag_its.resize(input_items.size());
    d_work_input_offsets.resize(input_items.size());
    int64_t output_offset = this->nitems_written(0);
    for (size_t input = 0; input < input_items.size(); input ++) {
        this->get_tags_in_window(d_work_tags[input], input, 0, ninput_items[input]);
        d_work_tag_its[input] = d_work_tags[input].begin();
        d_work_input_offsets[input] = this->nitems_read(input);
    }

    // size up histograms if needed; keeps output consistent
    if (nitems_max > d_nitems) {
      d_nitems = nitems_max;
      for (hist_iterator it = d_histograms.begin(); it != d_histograms.end(); it ++) {
        it->second.d_histograms.reserve(d_nitems * d_nbuckets);
        it->second.d_histograms.resize(d_nitems * d_nbuckets, 0);
        it->second.d_totals.reserve(d_nitems);
        it->second.d_totals.resize(d_nitems, 0);
      }
    }
    
    // loop through each item of the histogram, each packet sample
    for (int item = 0, sample = 0, hist_offset = 0; item < d_nitems; item ++, hist_offset += d_nbuckets)
    {
      // process the tags to match were we are in iterating the histogram
      for (size_t input = 0; input < input_items.size(); input ++) {
        uint64_t item_offset = d_work_input_offsets[input] + item;
        for (
          std::vector<gr::tag_t>::iterator & tag_it = d_work_tag_its[input];
          tag_it != d_work_tags[input].end() && tag_it->offset <= item_offset;
          tag_it ++
        ) {
          gr::tag_t & tag = *tag_it;
          prop_iterator prop_it = d_prop_tag_idxs.find(tag.key);
          // if a property tag is found, it changes the histogram
          if (prop_it != d_prop_tag_idxs.end()) {
            size_t idx = prop_it->second;
            d_prop_tags[idx] = tag.value;
            d_histogram = nullptr;
          }
        }
      }

      // set histogram from property tags if changed
      if (nullptr == d_histogram) {
        d_histogram = &d_histograms[d_prop_tags];
        if (0 == d_histogram->d_totals.size()) {
          d_histogram->d_histograms.reserve(d_nitems * d_nbuckets);
          d_histogram->d_histograms.resize(d_nitems * d_nbuckets, 0);
          d_histogram->d_totals.reserve(d_nitems);
          d_histogram->d_totals.resize(d_nitems, 0);
        }
      }

      // if this packet has this sample, add it to the histogram
      if (item < nitems_min) {
        for (int subtail = sample + d_vinlen; sample < subtail; sample ++) {
          input_type value = in[sample];
          size_t count = in2 ? in2[sample] : 1;
          if (value >= d_min && value < d_max) {
            size_t bucket = value * d_value_to_bucket_coeff + d_value_to_bucket_offset;
            d_histogram->d_histograms[hist_offset + bucket] += count;
            d_histogram->d_totals[item] += count;
          } else if (!warned) {
            GR_LOG_WARN(this->d_logger, "value " + std::to_string(value) + " outside histogram range [" + std::to_string(d_min) + ", " + std::to_string(d_max) + ")");
            warned = true;
          }
        }
      }
      
      // output the histogram
      if (std::is_floating_point<freq_type>::value) {
        size_t bucket = hist_offset;
        size_t tail = bucket + d_nbuckets;
        for (; bucket < tail; bucket ++) {
          *out = freq_type(d_histogram->d_histograms[bucket]) / d_histogram->d_totals[item];
          out ++;
        } 
      } else {
        // WARNING: the d_histograms type must be the same as freq_type
        static_assert(std::is_floating_point<freq_type>::value || sizeof(freq_type) == sizeof(d_histogram->d_histograms[0]));
        memcpy(out, &d_histogram->d_histograms[hist_offset], sizeof(freq_type) * d_nbuckets);
        out += d_nbuckets;
      }

      // output the counter
      if (out2) {
        *out2 = d_histogram->d_totals[item];
        out2 ++;
      }
    }

    return d_nitems;
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

protected:
  int calculate_output_stream_length(const gr_vector_int &ninput_items) {
      size_t length = std::max((size_t)ninput_items[0], d_nitems);
      if (ninput_items.size() > 1) {
        length = std::max(length, (size_t)ninput_items[1]);
      }
      return length;
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
typename tagged_stream_histogram<input_type, freq_type>::sptr
tagged_stream_histogram<input_type, freq_type>::make(input_type min, input_type max, size_t nbuckets, size_t vinlen, const std::vector<std::string> & prop_tag_keys, const std::string & len_tag_keys, const std::string & filename)
{
  return gnuradio::make_block_sptr<tagged_stream_histogram_impl<input_type, freq_type>>(min, max, nbuckets, vinlen, prop_tag_keys, len_tag_keys, filename);
}

template class tagged_stream_histogram<double, double>;
template class tagged_stream_histogram<float, double>;
template class tagged_stream_histogram<int64_t, double>;
template class tagged_stream_histogram<int32_t, double>;
template class tagged_stream_histogram<int16_t, double>;
template class tagged_stream_histogram<int8_t, double>;
template class tagged_stream_histogram<uint64_t, double>;
template class tagged_stream_histogram<uint32_t, double>;
template class tagged_stream_histogram<uint16_t, double>;
template class tagged_stream_histogram<uint8_t, double>;

template class tagged_stream_histogram<double, float>;
template class tagged_stream_histogram<float, float>;
template class tagged_stream_histogram<int64_t, float>;
template class tagged_stream_histogram<int32_t, float>;
template class tagged_stream_histogram<int16_t, float>;
template class tagged_stream_histogram<int8_t, float>;
template class tagged_stream_histogram<uint64_t, float>;
template class tagged_stream_histogram<uint32_t, float>;
template class tagged_stream_histogram<uint16_t, float>;
template class tagged_stream_histogram<uint8_t, float>;

template class tagged_stream_histogram<double, uint64_t>;
template class tagged_stream_histogram<float, uint64_t>;
template class tagged_stream_histogram<int64_t, uint64_t>;
template class tagged_stream_histogram<int32_t, uint64_t>;
template class tagged_stream_histogram<int16_t, uint64_t>;
template class tagged_stream_histogram<int8_t, uint64_t>;
template class tagged_stream_histogram<uint64_t, uint64_t>;
template class tagged_stream_histogram<uint32_t, uint64_t>;
template class tagged_stream_histogram<uint16_t, uint64_t>;
template class tagged_stream_histogram<uint8_t, uint64_t>;

} /* namespace openemissions */
} /* namespace gr */
