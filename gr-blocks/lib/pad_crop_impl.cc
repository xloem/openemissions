/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gnuradio/io_signature.h>

#include <openemissions/pad_crop.h>

namespace gr {
namespace openemissions {

class pad_crop_impl : public pad_crop
{
private:
  int d_itemsize;
  int d_length;
  bool d_pad;
  bool d_crop;

public:
  /*
   * The private constructor
   */
  pad_crop_impl(size_t itemsize,
                size_t length,
                const std::string& len_tag_key="packet_len",
                bool pad = true,
                bool crop = true)
  : gr::tagged_stream_block("pad_crop",
                   gr::io_signature::make(1, 1, itemsize),
                   gr::io_signature::make(1, 1, itemsize),
                   len_tag_key),
    d_itemsize(itemsize),
    d_length(length),
    d_pad(pad),
    d_crop(crop)
  {
  }

  int calculate_output_stream_length(const gr_vector_int &ninput_items) override
  {
    if (d_pad && d_length > ninput_items[0]) {
      return d_length;
    }
    if (d_crop && d_length < ninput_items[0]) {
      return d_length;
    }
    return ninput_items[0];
  }

  ~pad_crop_impl()
  {
  }

  int work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items) override
  {
    noutput_items = ninput_items[0];

    if (d_crop && noutput_items > d_length) {
      noutput_items = d_length;
    }

    std::vector<gr::tag_t> tags;
    get_tags_in_range(tags, 0, nitems_read(0), nitems_read(0) + noutput_items);
    for (size_t i = 0; i < tags.size(); i++) {
      gr::tag_t & t = tags[i];
      if (t.offset >= d_length) {
        break;
      }
      t.offset += nitems_written(0) - nitems_read(0);
      add_item_tag(0, t);
    }

    int input_size = noutput_items * d_itemsize;
    memcpy(output_items[0], input_items[0], input_size);

    if (d_pad && noutput_items < d_length) {
      int extra = d_length - noutput_items;
      memcpy((char*)output_items[0] + input_size, 0, extra * d_itemsize);
      noutput_items += extra;
    }

    return noutput_items;
  }

  void set_length(size_t length)
  {
    d_length = length;
  }

  size_t length() const
  {
    return d_length;
  }

  void set_pad(bool pad)
  {
    d_pad = pad;
  }
  
  bool pad() const
  {
    return d_pad;
  }

  void set_crop(bool crop)
  {
    d_crop = crop;
  }

  bool crop() const
  {
    return d_crop;
  }
};

pad_crop::sptr
pad_crop::make(size_t itemsize,
               size_t length,
               const std::string& len_tag_key,
               bool pad,
               bool crop)
{
  return gnuradio::make_block_sptr<pad_crop_impl>(
    itemsize, length, len_tag_key, pad, crop);
}

} /* namespace openemissions */
} /* namespace gr */

