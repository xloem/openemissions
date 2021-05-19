/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */


#ifndef INCLUDED_OPENEMISSIONS_PAD_CROP_H
#define INCLUDED_OPENEMISSIONS_PAD_CROP_H

#include <openemissions/api.h>
#include <gnuradio/tagged_stream_block.h> 

namespace gr {
namespace openemissions { 

class OPENEMISSIONS_API pad_crop : virtual public tagged_stream_block
{
public:
  typedef std::shared_ptr<pad_crop> sptr;

  static sptr make(size_t itemsize,
                   size_t length,
                   const std::string& len_tag_key="packet_len",
                   bool pad = true,
                   bool crop = true);

  virtual void set_length(size_t length) = 0;
  virtual size_t length() const = 0;

  virtual void set_pad(bool pad) = 0;
  virtual bool pad() const = 0;

  virtual void set_crop(bool crop) = 0;
  virtual bool crop() const = 0;
}; 

} // namespace digital
} // namespace gr

#endif /* INCLUDED_OPENEMISSIONS_PAD_CROP_H */
