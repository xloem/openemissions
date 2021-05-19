/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gnuradio/io_signature.h>

#include <openemissions/pigpio_sink.h>
#include "pigpiod.h"

namespace gr {
namespace openemissions {

using namespace detail;

template <class T>
class pigpio_sink_impl : public pigpio_sink<T>, public pigpiod::wave_sender
{
private:
  T d_level;
  double d_sample_rate;
  unsigned d_hardware_clock_frequency;
  int d_wave_buffer_percent;

  pigpiod::sptr d_server;

public:
  /*
   * The private constructor
   */
  pigpio_sink_impl(double samp_rate,
                   unsigned pin,
                   T level,
                   const std::string &address,
                   int wave_buffer_percent,
                   int hardware_clock_frequency)
  : gr::sync_block("pigpio_sink",
                   gr::io_signature::make(
                     1 /* min inputs */,
                     1 /*max inputs */,
                     sizeof(T)),
                   gr::io_signature::make(0, 0, 0)),
    pigpiod::wave_sender{
        .d_pin = 0,
        .d_accumulated_us = 0
    },
    d_level(level),
    d_sample_rate(samp_rate),
    d_hardware_clock_frequency(0),
    d_wave_buffer_percent(wave_buffer_percent)
  {
    set_hw(pin, address, hardware_clock_frequency);
  }

  /*
   * Our virtual destructor.
   */
  ~pigpio_sink_impl()
  {
    stop();
  }

  bool start() override
  {
    gr::thread::scoped_lock lk(d_server->d_mtx);
    pigpiothrow(set_mode(d_server->d_handle, d_pin, PI_OUTPUT));
    d_server->d_wave_senders[d_pin] = this;
    
    return sync_block::start();
  }

  bool stop() override
  {
    {
      gr::thread::scoped_lock lk(d_server->d_mtx);
      d_server->d_wave_senders.erase(d_pin);
    }
    pulses_filled();

    return sync_block::stop();
  }

  // Where all the action really happens
  int work(int noutput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items) override
  {
    const T *in = reinterpret_cast<const T*>(input_items[0]);
    bool last_state = false;
    bool might_send = (0 == d_accumulated_us);
    int last_us = 0;

    for (int sample = 0; sample <= noutput_items; ++ sample) {
      bool send_last_pulse;
      bool state;
      if (sample == noutput_items) {
        send_last_pulse = true;
        state = last_state;
      } else {
        state = in[sample] > d_level;
        if (sample == 0) {
          send_last_pulse = false;
        } else {
          send_last_pulse = (state != last_state);
        }
      }

      if (send_last_pulse) {
        gr::thread::scoped_lock lk(d_server->d_mtx);
        d_pulses.resize(d_pulses.size() + 1);
        gpioPulse_t & pulse = d_pulses.back();

        if (last_state) {
          pulse.gpioOn = d_pin;
          pulse.gpioOff = 0;
        } else {
          pulse.gpioOn = 0;
          pulse.gpioOff = d_pin;
        }

        double target_us = (sample - (double)0.5) * 1000000 / d_sample_rate;
        pulse.usDelay = target_us - last_us;
        d_accumulated_us += pulse.usDelay;
        last_us = target_us;
      }
    }

    if (might_send) {
      pulses_filled();
    }

    // Tell runtime system how many output items we produced.
    return noutput_items;
  }

  void set_pin(unsigned pin) override
  {
    set_hw(pin, address(), d_hardware_clock_frequency);
  }

  unsigned pin() const override
  {
    return d_pin;
  }

  void set_level(T level) override
  {
    d_level = level;
  }

  T level() const override
  {
    return d_level;
  }

  void set_sample_rate(double sample_rate) override
  {
    d_sample_rate = sample_rate;
  }

  double sample_rate() const override
  {
    return d_sample_rate;
  }

  void set_address(const std::string &address) override
  {
    set_hw(d_pin, address, d_hardware_clock_frequency);
  }

  const std::string &address() const override
  {
    return d_server->d_address;
  }

  void set_hardware_clock_frequency(unsigned hardware_clock_frequency) override
  {
    set_hw(d_pin, address(), hardware_clock_frequency);
  }

  unsigned hardware_clock_frequency() const override
  {
    return d_hardware_clock_frequency;
  }

  void set_wave_buffer_percent(int wave_buffer_percent) override
  {
    d_wave_buffer_percent = wave_buffer_percent;
    if (d_wave_buffer_percent > 50) {
      d_wave_buffer_percent = 50;
    }
    if (d_wave_buffer_percent < 1) {
      d_wave_buffer_percent = 1;
    }
  }

  int wave_buffer_percent() const override
  {
    return d_wave_buffer_percent;
  }

private:
  void pulses_filled()
  {
    gr::thread::scoped_lock lk(d_server->d_mtx);
    uint32_t min_us = ~0;
  
    for (auto & item : d_server->d_wave_senders) {
      pigpiod::wave_sender &sink = *item.second;

      if (0 == sink.d_accumulated_us) {
        return;
      }

      if (sink.d_accumulated_us < min_us) {
        min_us = sink.d_accumulated_us;
      }
    }

    pigpiothrow(wave_add_new(d_server->d_handle));

    for (auto & item : d_server->d_wave_senders) {
      pigpiod::wave_sender &sink = *item.second;

      uint32_t us_accumulation = 0;
      size_t pulse_idx = 0;
      for (
        us_accumulation = 0, pulse_idx = 0;
        us_accumulation < min_us;
        us_accumulation += sink.d_pulses[pulse_idx].usDelay, ++pulse_idx
      ) {
        if (pulse_idx >= sink.d_pulses.size()) {
          throw std::logic_error("iterated off end of pulse list, should have small enough min_us to prevent");
        }
      }

      // this is now the last pulse in the wave
      gpioPulse_t & pulse = sink.d_pulses[pulse_idx];

      // shorten the pulse to waveform length
      uint32_t extra_us = us_accumulation - min_us;
      pulse.usDelay -= extra_us;

      // add the waveform
      pigpiothrow(wave_add_generic(d_server->d_handle, pulse_idx, sink.d_pulses.data()));
      
      if (us_accumulation > min_us) {
        // adjust the rest of the pulse to be part of the next waveform
        pulse = {
          .gpioOn = 0,
          .gpioOff = 0,
          .usDelay = extra_us
        };
        -- pulse_idx;
      }

      // remove the added waveform from the queue
      sink.d_pulses.erase(sink.d_pulses.begin(), sink.d_pulses.begin() + pulse_idx);

      sink.d_accumulated_us -= min_us;
    }

    // send waveform
    int waveform = pigpiothrow(wave_create_and_pad(d_server->d_handle, d_wave_buffer_percent));
    pigpiothrow(wave_send_using_mode(d_server->d_handle, waveform, PI_WAVE_MODE_ONE_SHOT_SYNC));
    d_server->d_waveforms_in_flight.push_back(waveform);

    // delete waves that have completed
    // this loop assumes that the wave just sent prevents wave_tx_at from returning a no-wave-sent error
    while (d_server->d_waveforms_in_flight.front() != pigpiothrow(wave_tx_at(d_server->d_handle))) {
      pigpiothrow(wave_delete(d_server->d_handle, d_server->d_waveforms_in_flight.front()));
      d_server->d_waveforms_in_flight.pop_front();
    }
  }

  void set_hw(int pin, const std::string &address, unsigned hardware_clock_frequency)
  {
    if (d_server) {
      if (pin != d_pin || address != this->address()) {
        {
          gr::thread::scoped_lock lk(d_server->d_mtx);
          if (pigpiod::get(address)->d_wave_senders.count(pin)) {
            throw std::runtime_error("pin " + std::to_string(pin) + " already has output bound on " + address);
          }
          d_server->d_wave_senders.erase(d_pin);
        }
        pulses_filled();
      }
    }

    pigpiod::sptr new_server = pigpiod::get(address);
    d_pin = pin;

    {
      gr::thread::scoped_lock lk(new_server->d_mtx);

      pigpiothrow(set_mode(new_server->d_handle, d_pin, PI_OUTPUT));

      d_server = new_server;
      d_server->d_wave_senders[d_pin] = this;

      if (d_accumulated_us) {
        d_pulses.push_back({
          .gpioOn = 0,
          .gpioOff = 0,
          .usDelay = d_accumulated_us
        });
      }

      int r = hardware_clock(d_server->d_handle, pin, hardware_clock_frequency);
      if (r == PI_BAD_HCLK_FREQ && hardware_clock_frequency > 375000000) {
        r = hardware_clock(d_server->d_handle, pin, 375000000);
      }
      if (r == PI_BAD_HCLK_FREQ && hardware_clock_frequency > 250000000) {
        r = hardware_clock(d_server->d_handle, pin, 250000000);
      }
      d_hardware_clock_frequency = hardware_clock_frequency;
    }
  }
};

template <typename T>
typename pigpio_sink<T>::sptr
pigpio_sink<T>::make(double samp_rate, unsigned pin, T cutoff, const std::string& address, int wave_buffer_percent, unsigned hardware_clock_frequency)
{
  return gnuradio::make_block_sptr<pigpio_sink_impl<T>>(
    samp_rate, pin, cutoff, address, wave_buffer_percent, hardware_clock_frequency);
}

template class pigpio_sink<float>;
template class pigpio_sink<int>;
template class pigpio_sink<short>;
template class pigpio_sink<char>;

} /* namespace openemissions */
} /* namespace gr */

