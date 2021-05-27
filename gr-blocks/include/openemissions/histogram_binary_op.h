/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_OPENEMISSIONS_HISTOGRAM_DESUMMATION_H
#define INCLUDED_OPENEMISSIONS_HISTOGRAM_DESUMMATION_H

#include <openemissions/api.h>
#include <gnuradio/sync_block.h>

#include <functional>

namespace gr {
namespace openemissions {

/*!
 * \brief Histogram binary operation.
 * \ingroup openemissions
 *
 */
template <typename freq_type>
class OPENEMISSIONS_API histogram_binary_op : virtual public gr::sync_block
{
public:
  typedef std::shared_ptr<histogram_binary_op> sptr;

  /*!
   * \brief Return a shared_ptr to a new instance of openemissions::histogram_binary_op.
   * \param[in] min Smallest value at the start of the histograms
   * \param[in] max Tail of the histograms; values at this level are excluded
   * \param[in] op Callback that performs the binary operation
   * \param[in] nbuckets Vector length of the histograms; number of buckets
   */
  static sptr make(freq_type min, freq_type max, const std::function<double(double,double)> & op, size_t nbuckets = 1024);

  virtual void set_min(double min) = 0;
  virtual double min() const = 0;
  
  virtual void set_max(double max) = 0;
  virtual double max() const = 0;

  virtual void set_op(const std::function<double(double,double)> & op) = 0;
  virtual const std::function<double(double,double)> & op() const = 0;

  virtual size_t nbuckets() const = 0;
};

} // namespace openemissions
} // namespace gr

#endif /* INCLUDED_OPENEMISSIONS_HISTOGRAM_DESUMMATION_H */

