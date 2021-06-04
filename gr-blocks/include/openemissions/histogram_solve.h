/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_OPENEMISSIONS_HISTOGRAM_BINARY_SOLVE_H
#define INCLUDED_OPENEMISSIONS_HISTOGRAM_BINARY_SOLVE_H

#include <openemissions/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
namespace openemissions {

/*!
 * \brief Solves an equation of histograms
 * \ingroup openemissions
 *
 */
template <typename freq_type, typename... Doubles>
class OPENEMISSIONS_API histogram_solve : virtual public gr::sync_block
{
public:
  typedef std::shared_ptr<histogram_solve> sptr;

  /*!
   * \brief Return a shared_ptr to a new instance of openemissions::histogram_solve.
   * \param[in] min Smallest value at the start of the histograms
   * \param[in] max Tail of the histograms; values at this level are excluded
   * \param[in] expr Callback that evaluates an expression to solve.
   * \param[in] output_idx Parameter index to solve for.  Index 0 is the expression result.
   * \param[in] nbuckets Vector length of the histograms; number of buckets
   * \param[in] extrema Callbacks that evaluate to minima and maxima parameters in buckets
   */
  static sptr make(double min, double max, const std::function<double(Doubles...)> & expr, size_t output_idx = 0, size_t nbuckets = 1024, const std::vector<std::function<std::array<double, sizeof...(Doubles)>(Doubles...)>> & extrema = {});

  virtual void set_min(double min) = 0;
  virtual double min() const = 0;
  virtual void set_max(double max) = 0;
  virtual double max() const = 0;
  virtual void set_expr(const std::function<double(Doubles...)> &expr) = 0;
  virtual const std::function<double(Doubles...)> & expr() const = 0;
  virtual void set_extrema(const std::vector<std::function<std::array<double, sizeof...(Doubles)>(Doubles...)>> &extrema) = 0;
  virtual const std::vector<std::function<std::array<double, sizeof...(Doubles)>(Doubles...)>> & extrema() const = 0;
  virtual size_t output_idx() const = 0;
  virtual size_t nbuckets() const = 0;
};

} // namespace openemissions
} // namespace gr

#endif /* INCLUDED_OPENEMISSIONS_HISTOGRAM_BINARY_SOLVE_H */

