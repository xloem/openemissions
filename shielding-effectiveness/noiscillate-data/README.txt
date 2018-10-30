2018-08-31

Although the files showed appropriate increases in noise level, my attempts to measure the difference
between off-state and on-state inside a given file failed.  In order to debug it, I ran a live test, and
found
- my radio's LNA doesn't function above 43 dB anymore
- my noise generator no longer produces any noise at all.

I lowered the LNA gain and switched to a different noise generator that I have verified work.  I've
verified the 10 hz oscillation shows on a waterfall chart.  I made 4 new recordings, although I did not
watch the chart while specifically making the recordings.

The 2018-08-31 recordings are made with no antenna attached to the noise generator.  In the 'distant' ones,
the receiver is at least a couple feet from the noise generator and I expect the noise signal to be
nonobvious.  In the other two, the antenna is right against the noise generator's output, and the noise
signal is quite apparent in osmocom_fft .

Karl

ADDENDUM 2018-08-31:
I made a mistake recording the output of the commands, so I reran the tests with the output saved properly.
I took the time to make a third recording at a much greater distance (other side of the room, maybe at least
10 feet away).  I verified gqrx picked up the noise for the close recording, but again I did not verify
simultaneously with the recording.


2018-08-21

I used rtl_sdr to make 3 raw iq recordings from an rtl-sdr while running noiscillate at 40 hz on
a noise generator.

One without the noise generator active.  One with no antenna attached.  And one with an antenna.
I had to move the setup a ways to connect it to the big antenna.

I was running on a raspberry pi, so there could have been samples dropped due to write delay.

Karl

2018-10-15

I'm attempting to measure how measured wavelength of an emitter changes from recording to recording.
2018-10-15_noiscillate-40hz-500mhz-active-close_broken-rtlsdr-3-hot.iq
    tuned to 500 mhz, emitter is 40hz, equipment is broken and antenna must be placed right against emitter to see it; recording with dongle hot
    sample rate set to 2400000 Hz
      Periods: 23998
      Wavelength Avg: 0.02500093659276854
      Frequency Avg: 39.998501507709413
      Wavelength STD: 4.958333333333333e-05
2018-10-15_noiscillate-40hz-500mhz-active-close_broken-rtlsdr-3-cold.iq
  I placed the dongle in the fridge for a bit, then tried recording with it, the same setup:
      Periods: 23998
      Wavelength Avg: 0.025000955952827437
      Frequency Avg: 39.998470533959996
      Wavelength STD: 3.1250000000000001e-05
    More accuracy, but mean wavelength has changed slightly, perhaps due to PPM.
