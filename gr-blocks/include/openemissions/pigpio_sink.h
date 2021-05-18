/* -*- c++ -*- */
/* 
 * Copyright 2021 <+YOU OR YOUR COMPANY+>.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */


#ifndef INCLUDED_OPENEMISSIONS_PIGPIO_SINK_H
#define INCLUDED_OPENEMISSIONS_PIGPIO_SINK_H

#include <openemissions/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace openemissions {

    /*!
     * \brief <+description of block+>
     * \ingroup openemissions
     *
     */
    template <class T>
    class OPENEMISSIONS_API pigpio_sink : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<pigpio_sink> sptr;

      /*!
       * \param[in] pin GPIO to send output on
       * \param[in] level Pin is changed when input crosses this value
       * \param[in] samples_per_second Sample rate of input
       * \param[in] address Host of pigpiod to connect to
       * \param[in] hardware_clock_frequency For pins with a hardware clock, frequency to set the clock to
       */
      static sptr make(double samp_rate, unsigned pin = 4, T level = 0, const std::string &address = "127.0.0.1:8888", int wave_buffer_percent = 50, unsigned hardware_clock_frequency = 30000000);

      virtual void set_pin(unsigned pin) = 0;
      virtual unsigned pin() const = 0;

      virtual void set_level(T level) = 0;
      virtual T level() const = 0;

      virtual void set_samples_per_second(double samples_per_second) = 0;
      virtual double samples_per_second() const = 0;

      virtual void set_address(const std::string &address) = 0;
      virtual const std::string &address() const = 0;

      virtual void set_hardware_clock_frequency(unsigned hardware_clock_frequency) = 0;
      virtual unsigned hardware_clock_frequency() const = 0;

      virtual void set_wave_buffer_percent(int wave_buffer_percent) = 0;
      virtual int wave_buffer_percent() const = 0;
    };

  } // namespace openemissions
} // namespace gr

#endif /* INCLUDED_OPENEMISSIONS_PIGPIO_SINK_H */

