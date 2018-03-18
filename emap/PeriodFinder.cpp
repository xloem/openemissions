#include "PeriodFinder.hpp"

#include <iostream>

// TODO:
// Okay ... it looks like the video signal I thought I was finding was just noise in memory.
// Due to a bug in countPrecision, it was only calculating for the highest frequency.
//
// I thought I fixed the bug, but then it ran so slowly it didn't really do anything ...
// I have yet to demonstrate that convolving the complex stream with itself will show the refresh rate.
// At this low sample-rate, the absolute magnitude or the frequency domain may be needed, or peaks from distant parts of the spectrum that would reveal underlying shapes.
//
// I could perhaps check for this information quickly in gnuradio by performing an fft or a waterfall that actually covered 60 Hz
//      This would be done by decimating the signal down enough that a reasonably-sized fft would include 60 Hz.
// Gnuradio I believe has a peak detector that could be used afterwards, but I'd probably need to write my own blocks to find the precise period
// and then amplify it.

// TWO RUNTIME ERRORS:
// -- occasional segfault, somehow strengths[] is getting cleared! data = 0 at crash
// -- when totalBest is stored, strengths seem to become corrupted !

// TODO:
// Right now PeriodFinder enumerates all periods.
//
// Ideally it would find a highly outlying period, and only collect enough samples or use enough precisionFrequency to identify this period.
//
// Such functionality should be broken into pluggable chunks so this can be a useful tool.  The user is smarter than the heuristic.
//
// Let's factor out the period correlation function, and the buffer accumulation.  This class can be turned into a single function that
// finds a precisely high period.  Generally-high periods can be found with a separate function that perhaps uses downsample + fft.

// "LEVELS" of period finding:
// 1. general: perhaps use downsampled data to find wide swaths of periods
// 2. precise & memory intensive: find arrangement of periods over a long IQ recording that maximizes an evaluation function
// 3. unlimited precision: do #2 over and over, but don't store the data because we already know the frequency within sample resolution
//                         observe drift over time to increase precision of frequency

PeriodFinder::PeriodFinder(double minFrequency, double maxFrequency, double precisionFrequency, Source & source, double tunedHertz)
: window(GUIWindow::create()), rangeHz(minFrequency, maxFrequency), precisionHz(precisionFrequency), source(source), tunedHz(tunedHertz)
{
	init(source.sampleHertz());
	ready();
}

PeriodFinder::~PeriodFinder()
{
	done();
}

void PeriodFinder::init(double rate)
{
	expectedSamplingHz = rate;
	longestPeriod = expectedSamplingHz / rangeHz.first + 1;
	shortestPeriod = expectedSamplingHz / rangeHz.second;

	/*values.resize(longestPeriod - shortestPeriod);
	weights.resize(values.size());
	offsets.resize(values.size());
	for (unsigned i = 0; i < values.size(); ++ i) {
		values[i].set_size(shortestPeriod + i);
		weights[i].set_size(values[i].size());
		weights[i].zeros();
		offsets[i] = 0;
	}*/
}

/*
Memory could be conserved by having a concept of signal peaks, and narrowing down those peaks.
Alternatively we could just record the information for a wide variety of periods, and display those.

It seems we should start with the latter, to make a nice display that would aid in deciding what determines a peak.
*/

#include <itpp/stat/misc_stat.h>
double evaluateWaveform(itpp::cvec & waveform, double scale, typename itpp::cvec::value_type mean)
{
  // TODO: make this pluggable; waveform plug-in should probably handle recording the peak too
  return (itpp::sum(itpp::abs(waveform - mean))) / waveform.size() / scale;
}

void PeriodFinder::receiveQuadrature(itpp::cvec const & newData, double samplingHertz, double tunedHertz, double dBGain, double unixSecondsCompleted, class Source & source)
{
	if (&source != &this->source || samplingHertz != expectedSamplingHz || tunedHertz != this->tunedHz)
		return;

  // TODO: we probably don't need to store ALL of this old data ...
  // for sure!  it's crashing the prog due to memory exhaustion
	data = itpp::concat(data, newData);

	itpp::vec strengths(window->size().first);
	double totalBestValue  = -0.0/1.0;
  double totalBestPeriod = 0;
	itpp::vec totalBest;

  double shortLength = floor(shortestPeriod);
  double totalCount = floor(data.size() / shortLength);
  double longCount = floor((shortestPeriod - shortLength) * totalCount);

  itpp::cvec waveform(shortLength);
  typename itpp::cvec::value_type sampleSum;
  std::vector<int> periodStarts(totalCount);
  bool needToInitializeWaveform = true;


  for (;;) {
    // each iteration of this outer loop is a different fractional period to check, defined by this formula:
    double precisePeriod = shortLength + longCount / totalCount;

    if (precisePeriod > longestPeriod)
      break;

    double totalLength = shortLength * totalCount + longCount;

    // TODO: maybe not use countPrecision, increment longCount by one, and use a different heuristic for finding imprecise periods?
    double countPrecision = 1 / (1 / precisePeriod - precisionHz / expectedSamplingHz) - precisePeriod;

    if (totalLength <= data.size()) {

      // TODO: average waveform over entirety of periodPrecision, if countPrecision > 1
  
      // determine the average waveform for this period
      if (needToInitializeWaveform) {

        for (double i = 0; i < totalCount; ++ i)
          periodStarts[i] = i * totalLength / totalCount;

        sampleSum = 0;
        for (double j = 0; j < shortLength; ++ j) {
          waveform[j] = 0;
          for (double i = 0; i < totalCount; ++ i)
          {
            auto sample = data[periodStarts[i] + j];
            sampleSum += sample;
            waveform[j] += sample;
          }
        }

        needToInitializeWaveform = false;

      } else {
        // already have data from the last run of this shortLength; can just adjust it

        for (int i = 0; i < totalCount; ++ i) {
          int periodStart = i * totalLength / totalCount;
          if (periodStart != periodStarts[i]) {
            for (int j = 0; j < shortLength; ++ j) {
              auto delta = data[periodStart + j] - data[periodStarts[i] + j];
              sampleSum += delta;
              waveform[j] += delta;
            }
            periodStarts[i] = periodStart;
          }
        }
      }
  
      // measure and store the strength of this waveform
      double denom = totalCount * shortLength;
      double value = evaluateWaveform(waveform, denom, sampleSum / denom);
			strengths[strengths.size() * (precisePeriod - shortestPeriod) / (longestPeriod - shortestPeriod)] = value;

      if (value > totalBestValue) {
        totalBestValue = value;
        totalBestPeriod = precisePeriod;
        // TODO: uncommenting the following seems to corrupt things somehow ...
        //totalBest = itpp::abs(waveform);
      }

      // advance to the next period
      longCount += countPrecision;

    } else { // !(totalLength <= data.size())

      // longCount stretched us longer than we have data
      -- totalCount;

      // TODO: we don't actually have to reinitialize here; we can just remove the values being dropped
      needToInitializeWaveform = true;

    }

    if (longCount >= totalCount) {
      do {
        longCount -= totalCount;
        ++ shortLength;
      } while (longCount >= totalCount);

      waveform.set_size(shortLength);
      needToInitializeWaveform = true;
    }
  }


  /*

	for (double shortLength = floor(shortestPeriod); shortLength < longestPeriod; ++ shortLength) {
    double longCount = 0;
    int totalCount = data.size() / shortLength;

		double bestPeriod;
		double bestValue = -0.0/1.0;

		itpp::cvec waveform(shortLength);
		for (double increment = 0; increment < increments; ++ increment) {
			waveform.zeros();
			for (double i = 0; i < increments; ++ i) {
				int offset = i * ((shortLength - 1) * increments + increment) / (increments - 1);
				waveform += data.get(offset, offset + waveform.size());
			}
			double value = evaluateWaveform(waveform);
			double period = shortLength + increment / increments;
			strengths[strengths.size() * (period - shortestPeriod) / (longestPeriod - shortestPeriod)] = value;
			if (value > bestValue) {
				bestIncrement = increment;
				bestValue = value;
				bestPeriod = period;
				if (bestValue > totalBestValue) {
					totalBestValue = bestValue;
					totalBest = itpp::abs(waveform);
				}
			}
		}
		if (bestIncrement > 0 && bestIncrement < increments - 1) {
			strengths[strengths.size() * (bestPeriod - shortestPeriod) / (longestPeriod - shortestPeriod)] = bestValue;
		}
	}
  */

	//window->setLines({strengths, totalBest});
  window->setLines({strengths});
  std::cout << "Peak = " << (1 / totalBestPeriod * samplingHertz) << " Hz   " << (totalBestPeriod / samplingHertz) << " s" << std::endl;
  std::cout << "Value = " << totalBestValue << std::endl;
		/*
but really we should choose X - X/period to begin
and stretch it X/period increments to the right.
This increases correctly.
Just confused because my math was wrong.

Okay, so ...
		*/
	//}
	
/*	if (&source == &this->source && samplingHertz != expectedSamplingHz)
		init(samplingHertz);

	for (int row = 0, len = shortestPeriod; len < longestPeriod; ++ row, ++ len) {
		auto valueRow = values[row];
		auto weightRow = weights[row];
		auto & offset = offsets[row];

		int dataOffset = 0;
		weights[*/
/*
Well, my brain isn't quite right here.
The plan is to sum the complex fourier transform over a range of different sizes ... it may
not be well thouight out.
What I was implementing was to sum the actual data over a range of different periods ... it
will probably work, but I think the fourier approach is better.
The goal is basically to _perform_ a fourier transform, but over very low frequencies.
I'm interested in the power at very precise frequencies of whatever arbitrary repeating
signal there may be at those frequencies.  Step #1 is to recreate a video signal from
display emanations.  These signals repeat at roughly 60 hz, but in reality it is some precise
value unique to the particular display.
So if we just want to _find_ such repeating data, we can look for peaks in autocorrelations
that repeat at a variety of frequencies.  A range of frequencies as finely sampled as we
can manage.
We have data coming in at 2 MHz, radio data.  So I can just fill repeating buckets at varying
lengths (1/2MHz increments).
A problem I considered had something to do with phase shift ... my preferred approach is to
sum many subsequent instances of these periods.  The paper used an "inferior" approach which
correlated an instance with one occuring far in the future.  This makes echos look like
signals and it is not robust in the face of noise.  It is better to sum many signals.
But if my sampling windows are off from the actualy frequency by even a little, my sum
will go to zero as we phase by the actual signal.
One approach in that scenario is to sum the spectral content instead of the raw signal.
A spectrogram created with a much higher sampling rate than the ~60 Hz rate of interest will
"blur" the interesting features and they should show through the sum despite small phase
differences.
I was trying to think of how to do this without making arbitrary spectrogram buckets, by
perhaps taking one big huge fourier transform.  Somehow I thought I should take a fourier
transform of subsequent sections of data equal to the period length of interest -- a range
of subsequent fourier transforms -- and sum them.  The hope is then that the summed fourier
transform would represent one taken of the underlying signal for that precise period.
Which it would!
But if I use the complex fourier transform, I imagine we'd have a similar problem with
phase offsets.
If I instead took the absolute value of the fourier transform, then the positioning information of the signal components would be lost (the phase information) and the data would not
be reusable to extract the actual underlying signal.  It seems better to create reusable data.

Options
	- Sum raw data
	- Sum buckets of fourier transforms (spectrograms, contain original signal)
	- Sum complex fourier transform of raw data (contains same data as raw data)
	- Sum abs fourier transforms of raw data (seems more direct than the spectrogram)
	- Work with fourier transforms of huge chunks of data (promising but uncertain)

Okay, so I take the absolute value of the fourier transform of this data.
I want information on the data that repeats at precisely this period

Basically I have bursts of high frequency data.  These bursts happen in a pattern that repeats at a low frequency.

If I take a fourier of a big long sample of this data, these high frequency repeating bursts ... I don'[ quite have the intuition for this, but they should make relatively low frequency components.

So my data has low frequency components, even though its fourier transform is about high frequency information ...  But these components may be very weak.  We can recover them by finding the precise frequency at which they are at, and summing the data.  I don't really have the sinusoidal intuition here.

So I'll use a repeated summing approach.

So, basically what's gotten scrambled in my head, is this midway conclusion:
- Spectrogram breaks the data into buckets.  We could make these buckets be the period length
of interest.  Then we need only one transform per period length, and we don't have to worry
about issues where the bucket sizes don't quite line up.

- On the other hand, using buckets could be just as good ... but we have to be willing to
come up with a good bucket size.  The advantage of a spectrogram is we only need to break
the signal into buckets once, and make one spectrogram of it, and then we just resample that
spectrogram to handle different period lengths.
If we use windows that are the same size as the periods, then we have to take geometrically
more fourier transforms.  n^2/2 as many.  If I knew more about fourier transforms, or
spent more time to figure it out, there's probably a way to re-use a lot of the calculation
that gets done.  But I think using spectrogram buckets could be the way to go fo rnow.
The other disadvantage of the spectrogram is that it is much less precise.
That's more of a problem.  Precision could help a lot ...

I guess I'll try the N^2 fourier transforms for now?
NOTE THIS THOUGH!  We can't use these transforms to reconstruct the actual signal, although
they do give us interesting data on it such as perhaps a likely clock frequency, so we
shouldn't be using the fourier transform here unless we need to.
If we can identify maxima without using the fourier transform, it could be valuable.

Perhaps the spectrogram approach could be valuable here, because we could switch techniques
once we want to become more precise.  We'd switch to a technique that does not lose the
phase information, and precisely align the periods such that they align in phase.  Then the
signal may be reconstructed.

So, two-phased approach:
PeriodFinder: finds signals using some kind of rudimentary spectral decomposition.  wants to
              look for correlation peaks in a spectrogram.  Probably a simpler approach exists.
PeriodSignalRecorder? :  First narrows down a period found by PeriodFinder
                         Then spews out a detailed waveform of it.
I guess PeriodFinder could do this too, once it has a period identified.

Okay, so I guess I should compare "correlation peaks in a spectrogram" to "simpler approach exists".  QUITE POSSIBLE that these low-frequency peaks would show up in a simple fourier transform, for example.

So, let's compare these two approaches.  Well, I feel pretty confident the spectrogram approach will work.  I'm quite curious how non-sinusoidal repeating waveforms show up in a fourier transform.  Hrm, didn't I learn about this already, once, long ago?
It shows up as harmonics.  A repeating waveform at the frequency of a sinusoid, not shaped like the sinusoid, has a component for the sinusoid, and successive components at integral fraction multiples of the sinusoid.  The differing strengths of these harmonic peaks represent the shape of the signal.

I could imagine processing a fourier transform, knowing I'm looking for a particular peak, to actually sum the harmonics together.  This could be a valuable way to process a fourier transform, and I imagine it representing better the peak of interest.

I guess a meaningful approach might be to play around in e.g. matlab for a bit.
Make a repeating signal of random data.
Look at the fourier of the random data.
See how it goes.

But obviously the spectrogram approach will work.  The N^2 fourier approach will work.
And this approach of summing harmonic peaks should work too.  I guess you'd want to look at the weighted average of the resulting summed fourier window and verify that it was way over to the edge.

---

After thinking about this for a bit, I'm settling into the spectrogram approach.
Buckets of fourier transforms (spectrograms) for blurriness.
Then a precise compare for honing right in.
I was in pretty weird states of mind throughout this.  Made this decision suddenly.  But it does seem good and will move the pro0ject forward.

----

Okay, I sent an e-mail to myself with at least somewhat more clarity.  Use basic summation to find the peaks; we'll skip the blurry approach.  Peaks are one sample wide, so we check all possible samples.
Accumulate the data in the peak finder.
Use a separate thread for radio from processing.
*/
/*		
		int idx, weight;
		for (idx = 0, weight = 0; idx < data.size() - len; idx += len, ++ weight) {
			
		}
	}
*/
}
