/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_OPENEMISSIONS_HISTOGRAM_PDF_H
#define INCLUDED_OPENEMISSIONS_HISTOGRAM_PDF_H

#include <openemissions/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
namespace openemissions {

/*!
 * \brief Per-sample histograms
 * \ingroup openemissions
 *
 */
template <typename input_type, typename freq_type>
class OPENEMISSIONS_API histogram : virtual public gr::sync_block {
  public:
  typedef std::shared_ptr<histogram> sptr;

  /*!
   * \brief Return a shared_ptr to a new instance of
   *openemissions::histogram.
   * \param[in] min Smallest value at the start of the histogram
   * \param[in] max Tail of the histogram; values at this level will be excluded
   * \param[in] nbuckets Vector length of the output; number of histogram buckets
   * \param[in] vinlen Vector length of the input; samples will be consolidated together
   */
  static sptr make(input_type min, input_type max, size_t nbuckets = 1024, size_t vinlen = 1, const std::vector<std::string> & prop_tag_keys = {}, const std::vector<std::string> & len_tag_keys = {}, const std::string & filename = "");

  virtual void set_min(input_type min) = 0;
  virtual input_type min() const = 0;

  virtual void set_max(input_type max) = 0;
  virtual input_type max() const = 0;

  virtual size_t nbuckets() const = 0;

  virtual uint64_t total() const = 0;
};

} // namespace openemissions
} // namespace gr

#endif /* INCLUDED_OPENEMISSIONS_HISTOGRAM_PDF_H */
