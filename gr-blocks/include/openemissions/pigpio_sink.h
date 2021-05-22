/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */


#ifndef INCLUDED_OPENEMISSIONS_PIGPIO_SINK_H
#define INCLUDED_OPENEMISSIONS_PIGPIO_SINK_H

#include <openemissions/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
namespace openemissions {

/*!
 * \brief Raspberry Pi GPIO pin sink
 * \ingroup openemissions
 *
 */
template <class T>
class OPENEMISSIONS_API pigpio_sink : virtual public gr::sync_block
{
public:
  typedef std::shared_ptr<pigpio_sink<T>> sptr;

  /*!
   * \param[in] samp_rate Samples per second of input
   * \param[in] pin GPIO to send output on
   * \param[in] level Pin is changed when input crosses this value
   * \param[in] address Host of pigpiod to connect to
   * \param[in] wave_buffer_percent Portion of the device dma resources to allocate for each waveform
   * \param[in] hardware_clock_frequency For pins with a hardware clock, frequency to set the clock to
   */
  static sptr make(double samples_per_second,
                   unsigned pin = 4,
                   T level = 0,
                   const std::string &address = "127.0.0.1:8888",
                   int wave_buffer_percent = 25,
                   unsigned hardware_clock_frequency = 30000000,
                   unsigned pad_milliamps = 0);

  virtual void set_pin(unsigned pin) = 0;
  virtual unsigned pin() const = 0;

  virtual void set_level(T level) = 0;
  virtual T level() const = 0;

  virtual void set_sample_rate(double samples_per_second) = 0;
  virtual double sample_rate() const = 0;

  virtual void set_address(const std::string &address) = 0;
  virtual const std::string &address() const = 0;

  virtual void set_hardware_clock_frequency(unsigned hardware_clock_frequency) = 0;
  virtual unsigned hardware_clock_frequency() const = 0;

  virtual void set_wave_buffer_percent(int wave_buffer_percent) = 0;
  virtual int wave_buffer_percent() const = 0;

  virtual void set_pad_milliamps(unsigned pad_milliamps) = 0;
  virtual unsigned pad_milliamps() const = 0;
};

} // namespace openemissions
} // namespace gr

#endif /* INCLUDED_OPENEMISSIONS_PIGPIO_SINK_H */

