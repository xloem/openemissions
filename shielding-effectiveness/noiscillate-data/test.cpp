#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

#include "types.hpp"

#include "stats.hpp"

#include "rtlsdriqdump.hpp"
#include "soapylive.hpp"

////////
// CURRENT STATUS AND PLAN
//
// dB exist but will be unsatisficatory when signal level is below noise level.
// I think more accuracy could be acquired by taking the FFT of the squared magnitude rather than the magnitude; I think the dB may sum under that condition.  Test is needed.
// 
// Whatever the current strategy for finding the dB is, any harmonics that have dB above some cutoff may be used to identify frequency better.
// That's for FFT.  I think direct period comparison approaches will wok more accurately, but we'll need to factor out the algorithm for choosing which periods to try; the current one-by-one approach is too slow
// unless we can start with more accuracy.  Things can also be sped up by creating one stats class that parallelizes all the periods.

////////
// Found the dB with the FFT.
// Unforunately, it's not very precise.  Also, not knowing the frequency precisely will make it hard to separate data with different attributes if the source is in motion.
//
// I think the next step should be to get a more accurate frequency from the harmonics, as in APPROACH FOR FINDING MORE PRECISE FREQUENCY below.
// Then the device can be calibrated before measurement is synchronized with motion.
// (I imagine a source that is moving in a planned way to measure e.g. SE at different distances.  To get the value at 1 single distance may require combining samples that are not
// near each other but have well-known recording time.  This won't be conducive to using FFT, but having accurate freq will allow other techniques to be used relatively easily.)
//
// NEXT: Try APPROACH FOR FINDING MORE PRECISE FREQUENCY below.  (can also sum all peaks for more precise mag)

////////
// I WAS CONFUSED
//
// LET'S GO TO THE FINISH LINE.
//
// The fft can find the frequency roughly.
// We can get the dB from the FFT.  Is it usable?
//
// With an rough frequency, we can get the dB approximately, and redo the FFT to adjust for phase drift next buffers.
// Will need to:
// 1. accumulate enough buffers for FFT
// 2. use perhaps existing 4-bin approach to measure dB.
//
// ALTERNATIVELY
// if you'r having trouble going to the goal, it would be helpful to merge this work into emap, or to merge emap into gnuradio or pothos.

// Can we get dB from FFT?
// 1. switch FFT to complex? =S phase adjustment would be needed
// 1. compare magnitude of FFT to magnitude of wave.

// it's notable that harmonics give a much more precise frequency than the root frequency.  they can likely be found just by the max of the right area.
// APPROACH FOR FINDING MORE PRECISE FREQUENCY:
// 1. fft range gives range for each harmonic
// 2. try harmics, pick real ones and treat maxes as frequencies
//    -> criteria: must be statistically signicant given background
//                 must be have frequency within range of fundamental
// 3. divind range for harmonic fft by harmonic number gives a different max and min for fundamental frequency
// 4. combining max and mins gives narrower range for fundamental frequency
//
// APPROACH FOR FINDING dB from FFT:
// 1. height of peak is proportional to intensity of signal +- noise mag
// - adding fundamentals could increase accuracy by increasing snr

// thoughts on harmonics:
// A problem exists where each harmonic will have [decreasing] statistical significance.
// One approach could be to assume that the harmonics exist, without regard to how strong they are relative to noise.
// I've been assuming the signal magnitude is the height of the peak +- noise mag -- but once the sig mag gets low towards the noise mag, I think this will change.
// e.g. if the sig mag is 0.1 and the noise mag is 1.0, the range will be 0.0 to 1.1.  Most of the time the harmonic will be invisible.
// Say the peak is 1.05 and the real sig mag is 0.1 .  This means the noise mag was 0.95 .  
//
// 0.1 sig.  Could have any phase or sign.  1.0 noise.  Summed with sig.
//
// So if I measure a 1.1 peak, and the noise is 1.0, we can basically assume the peak is nonexistent because the noise change could bring it below zero.
//
// How would I sum multiple peaks?  When I have 4 of a single peak, I can take it as a sampling distribution of the means.
// If I have two different peaks, is it different?  if I sum them the noise should cancel out .... but that assumes that the noise contributes in both an additive and subtractive manner, which is only true
// if the peak is louder than the loudest expected noise.
//
// So I can't blindly sum all the peaks -- some of the noise will not cancel itself out.  I'll just want loud peaks.  It seems I'll need to make some kind of assumption about the noise; I'll need some kind of
// confidence interval or such.  The % confidence I want should be related to the number of trials I want to require before a msitake is found.  I've been doing a year of runs, which involves dividing the time
// in a year by the time of the trial.
                   

template <typename Scalar, class DataSource, class DataProcessor>
class DataFeeder
{
public:
  DataFeeder(DataSource & source,
             DataProcessor & processor,
             std::function<bool()> checkDataInterruptedHook = [](){return false;})
  : _source(source),
    _processor(processor),
    _checkInterruptedHook(checkDataInterruptedHook)
  {
    _lastMeta.sampleTime = 0;
    _lastMeta.freq = -1;
    _lastMeta.rate = -1;
  }

  template <typename Derived>
  void add(Eigen::PlainObjectBase<Derived> const & chunk, RecBufMeta const & meta)
  {
    using Eigen::numext::mini;

    size_t minNeeded = _processor.bufferSizeMin();
    size_t multipleNeeded = _processor.bufferSizeMultiple();
    size_t maxNeeded = _processor.bufferSizeMax();

    size_t chunkOffset = 0;
    size_t nextSize;
    
    // exhaust buffer
    if (_lastMeta.sampleTime + _buffer.size() != meta.sampleTime ||
        _lastMeta.freq != meta.freq ||
        _lastMeta.rate != meta.rate)
    {
      // buffer doesn't contain enough samples, and next chunk is not contiguous with it
      // drop it
      _buffer.conservativeResize(0);
      _lastMeta = meta;
    }
    if (_buffer.size())
    {
      nextSize = mini((size_t)_buffer.size() + (size_t)chunk.size(), maxNeeded);
      if (multipleNeeded > 1)
      {
        nextSize -= nextSize % multipleNeeded;
      }
      if (nextSize < minNeeded)
      {
        _buffer.conservativeResize(_buffer.size() + chunk.size());
        _buffer.tail(chunk.size()) = chunk;
        return;
      }
      size_t chunkOffset = nextSize - _buffer.size();
      _buffer.conservativeResize(nextSize);
      _buffer.tail(chunkOffset) = chunk.head(chunkOffset);
      _processor.process(*this, _buffer, _lastMeta);
      _lastMeta.sampleTime += _buffer.size();
    }

    // pass segments of chunk
    while (!checkInterrupted())
    {
      minNeeded = _processor.bufferSizeMin();
      multipleNeeded = _processor.bufferSizeMultiple();
      maxNeeded = _processor.bufferSizeMax();
      nextSize = mini(chunk.size() - chunkOffset, maxNeeded);
      if (multipleNeeded > 1)
      {
        nextSize -= nextSize % multipleNeeded;
      }
      if (nextSize < minNeeded)
      {
        break;
      }

      _processor.process(*this, chunk.segment(chunkOffset, nextSize), _lastMeta);
      _lastMeta.sampleTime += nextSize;
      chunkOffset += nextSize;
    }

    nextSize = chunk.size() - chunkOffset;
    assert(_lastMeta.sampleTime == meta.sampleTime + chunk.size() - nextSize);
    _buffer.conservativeResize(nextSize);
    _buffer = chunk.tail(nextSize);
  }

  bool checkInterrupted() { return _checkInterruptedHook(); }

  DataSource & source() { return _source; }

  DataProcessor & processor() { return _processor; }

  // Harmonics
  //
  // Integer multiples of the real fundamental will contain harmonics.
  // Harmonics magnitude will be real harmonic magnitude +- noise magnitude, so we can only see a harmonic if it is louder than the noise, right?



  // Harmonic problem: my perception indicates that harmonics will only be recordable if louder than noise.  This doesn't make sense: we should get a small statistical change for even a very weak harmonic.  How to merge
  // these two perceptions?
  //
  // Okay, I _measured_ signal strength amidst noise, and found the height of the peak was proportional to the strength of the signal, +- noise.
  //
  // My intuition is that over many summed recordings, a weak peak will stand out amidst noise.
  //
  // NOISE 
  //    VALUE has a mean value of zero; it is sometimes loud sometimes weak, and the phase is random
  //    MAGNITUDE has a mean value of nonzero: the value is sometimes zzero sometimes loud, and the average is somewhere between
  //    STD DEV OF VALUE has a mean nonzero value which is different from the mean of the magnitude (sum of the square roots vs square root of the sum)
  //                     it stays the same as more data is added
  //    STD DEV OF MAGNITUDE ACROSS TIME has a mean nonzero value
  //                     I'm thinking this value stays the same as more data is added
  //
  // Over many 
  //
  //
  // okay, so what ends up being significant here is that I am comparing adjacent FFT bins with each other -- not the data over time, but rather over frequency.
  // As I add more data, the standard dev of the mean will decrease because I have more samples.  Each fft bin is a sampling mean, and they will become more similar over time, each approaching the true noise mean.
  // the noise + signal bins will have different behavior: as the mean is taken, the further samples of the noise tend to move closer to the average.
  // So, on average, a powerful signal will move towards its real value as the noise cancels out.
  //
  // How does a weak signal behave?
  //
  // If the signal is only a little louder than the noise, we will see the signal's value.
  // If the signal is less loud than the noise, what do we see? ?????
  //
  // Okay, the FFT result is actually complex.  It's a vector in the complex plane.  It has random phase and random magnitude, wrt noise component.
  // The signal component has constant magnitude, and changing phase.
  //
  // If we sum random magnitude @ random phase, with constant magnitude @ changing phase, how do we get a result that isolates the constant magnitude?
  //
  // What's the magnitude of the result?  random magnitude will sometimes be small -- at those times, the result will be the constant magnitude +- the random magnitude
  // When the random magnitude exceeds the size of the cosntant magnitude, funny things can happen.  When they are in-phase, the result doubles.  When they are out-of-phase, the result gets super tiny.
  //
  // What's the mean?
  //
  // Since sin/cos are symmetrical, we can likely find the mean by average the extrema.  no ... not quite ... the volume density may change near zero compared to far away.
  //
  // What we're doing here is convolving two spheres.  One sphere is randomly sized, and one sphere is sized a constant, small amount.  We want to know the center of mass of the convolved sphere.
  // The spheres are not filled.  They are surfaces only.
  // If I convolve one sphere with another, I get two concentric spheres with the insides filled.  The center of mass of two concentric spheres is their center.
  //
  // Okay, errors:  these are circles, not spheres.  the complex plane is 2 dimensional.
  // The cosntant circle is progressing at a cosntant rate.  It's not a random circle like the noise.  Let's try treating it like a constant vector.
  //
  // A vector plus a circle equals an offset circle.  So the random noise is offset from the origin by the constant vector.
  //
  // If I take the magnitude of the result in order to undo the effect of the phase changes to the constant value, what magnitude will I see?
  //
  // -> note, we could undo the phase progression if we took the FFT of the FFT bin values; then we wouldn't need to take the magnitude
  //
  // the resulting average magnitude will be the average of all the distances of all the points on that circle from the origin.
  //
  // given that (x - c)^2 + y^2 = r^2, what is integral(sqrt(x^2 + y^2), dt, 0, 2pi)?
  // x = r cos(t) + c
  // y = r sin(t)
  //
  // integral(x^2 + y^2, dt, 0, 2pi) / (2pi)
  //
  // x^2 = (r cos(t) + c) * (r cos(t) + c)
  // x^2 = r cos(t) * (r cos(t) + c) + c r cos(t) + c^2
  // x^2 = r cos(t) * r cos(t) + r cos(t) * c + c r cos(t) + c^2
  // x^2 = r^2 cos^2(t) + r cos(t) * c + c r cos(t) + c^2
  // x^2 = r^2 cos^2(t) + 2 r c cos(t) + c^2
  // y^2 = r^2 sin^2(t)
  //
  // x^2 + y^2 = r^2 cos^2(t) + 2 r c cos(t) + c^2 +  r^2 sin^2(t)
  // x^2 + y^2 = r^2 + 2 r c cos(t) + c^2 
  // x^2 + y^2 = r^2 + c^2 + 2 r c cos(t)
  //
  // integral(x^2 + y^2, dt, 0, 2pi) / 2pi = 
  // r^2 + c^2 + 2 r c integral(cos(t), dt, 0, 2pi) / 2pi
  // = r^2 + c^2
  //
  // so
  // I think I got that the mean of the square of the mixed signal's value will be equal to the magnitude of the noise squared plus the magnitude of the signal squared
  // unfortunately I assumed noise with constant magnitude and random phase!  in reality our noise will have random magnitude
  //
  // additionally, the integral will depend on the magnitude distribution
  //
  // but what will happen is that r will change randomly, whereas c will not.
  //
  // so it's looking to me like if I sum the square of the value instead of the magnitude, the result will = the mean of r^2 (which I think is the variance) plus the signal^2

private:
  DataSource & _source;
  DataProcessor & _processor;
  HeapVector<Scalar> _buffer;
  std::function<bool()> _checkInterruptedHook;
  RecBufMeta _lastMeta;
};

/*
template <typename T>
class StatsAccumulatorHistogram
{
public:
  using Scalar = T;

  StatsAccumulatorHistogram(T binWidth, T start = 0)
  : binWidth(binWidth),
    binDensity(1.0 / binWidth),
    binStart(start * binDensity),
    count(0),
    dataSum(0),
    dataSquaredSum(0),
  { }

  void add(Eigen::DenseBase<Derived> const & chunk)
  {
    using Eigen::numext::floor;

    count += chunk.size();
    dataSum += chunk.sum();
    dataSquaredSum += chunk.derived().matrix().squaredNorm();

    auto chunkScaled = (chunk.derived().array() * binDensity - binStart).eval();

    int minIdx = floor(chunkScaled.min());
    if (minIdx < 0)
    {
      size_t shift = -minIdx;
      size_t oldSize = bins.size();
      bins.resize(oldSize + shift);
      bins.tail(oldSize) = bins.head(oldSize);
      bins.head(shift).setZero();
      binStart -= shift;
      chunkScaled += shift;
      minIdx = 0;
    }

    auto chunkByIndex = chunkScaled.cast<size_t>();

    size_t maxIdx = chunkByIndex.max();
    if (maxIdx >= bins.size())
    {
      bins.resize(maxIdx + 1, 0);
    }
    
    ++ bins[chunkByIndex.redux([this](T lastIdx, T nextIdx) -> T {
      ++ bins[lastIdx];
      return nextIdx;
    })];
  }

  T const & mean() const
  {
    return dataSum / count;
  }

  T const & variance() const
  {
    // sum((x - mean)^2) / n1
    // sum(x^2 - 2 * mean * x + mean^2) / n1
    // sum(x^2 + mean^2 - 2 * mean * x) / n1
    // (sum(x^2) + N*mean^2 - 2 * mean * sum(x)) / n1
    // (sum(x^2) + mean*(N*mean - 2 * sum(x))) / n1
    //
    // mean = sum(x) / N
    //
    // (sum(x^2) + mean*(N*sum(x) / N - 2 * sum(x))) / n1
    // (sum(x^2) + mean*(sum(x) - 2 * sum(x))) / n1
    // (sum(x^2) - mean*sum(x)) / n1
    return (dataSquaredSum - dataSum * mean()) / (count - 1);
  }

  template <typename StatsAccumulator>
  T chanceMeanGivenBaseline(StatsAccumulator const & baseline)
  {
    using Eigen::numext::erfc;
    using Eigen::numext::abs;
    // get standard error of means from baseline
    // return erfc(abs(averageDelta) / stderr)
  }

  T chanceAverageThisLarge() const
  {
    // gives a metric with regard to how likely this data is a random fluke

    // assume the population has a mean of zero
    // estimate its variance with sum(x) / (n - 1)
    // then its standard error of means is sqrt(variance / n)
    //   stderr = sqrt(sum(x) / (n * (n - 1)))
    
    // the chance we got an average this large is then
    //   1.0 - erf(average / stderr)

    // reduce:
    // average / stderr
    // (dataSum / n) / sqrt(dataSum / (n * (n - 1)))
    // sqrt(dataSum^2 / n^2 / (dataSum / (n * (n - 1))))
    // sqrt(dataSum^2 * (n * (n - 1)) / (n^2 * dataSum))
    // sqrt(dataSum * (n * (n - 1)) / n^2)
    // = sqrt(dataSum * (n - 1) / n)
    return erfc(sqrt(dataSum * (count - 1) / count));
  }

private:
  T binWidth;
  T binDensity;
  T binStart;
  HeapVector<T> bins;
  size_t count;
  T dataSum;
  T dataSquaredSum;

  inline void add(T val)
  {
    int idx = (val - binStart) * binDensity;
    if (idx < 0) {
      size_t shift = -idx;
      idx = 0;
      bins.insert(0, shift, 0);
      binStart -= shift * binDensity;
    }
    else if (idx >= bins.size())
    {
      bins.resize(idx + 1, 0);
    }
    ++ bins[idx];
  }

  T binToValue(size_t idx)
  {
    return idx * binWidth + binStart;
  }
};

*/

// TODO: our signal _IS NOISE_ which complicates things!
// to resolve this you just have to start writing about what it implies.
// note 1: our signal is _SIGNAL NOISE_ summed onto _BACKGROUND NOISE_: so we expect it, when active, to be simply a differing gaussian dist. 
//            since we've taken the absolute magnitude of the signal, the mean will rise, but otherwise it would still be zero and only the variance would have risen.

// TODO: run it!  the pieces are almost all there!

// uses 4 stats bins in an attempt to quickly identify noise that's toggled in a square wave with 50% duty cycle
// TODO: make a model like this but that uses a tree of stats accumulators for adaptive subsampling
template <typename StatsAccumulator, size_t NUM_BINS = 4>
class PeriodModelAccumulatorNoiscillate
{
public:
  using Detail = StatsAccumulator;
  using Scalar = typename StatsAccumulator::Scalar;

  PeriodModelAccumulatorNoiscillate(size_t period, Detail const & accum = {})
  : offset(0),
    bins{{accum, accum, accum, accum}},
    period(period),
    periodStride(period)
  { }

  template <typename Derived>
  void add(Eigen::PlainObjectBase<Derived> const & chunk)
  {
    // 1. fill first partial bin to the end, if needed
    size_t bin = offset * NUM_BINS / period;
    size_t tail = binStart(bin + 1);
    size_t chunkPos = tail - offset;
    if (chunkPos >= chunk.size())
    {
      bins[bin].add(chunk);
      offset = (offset + chunk.size()) % period;
      std::cerr << offset << " / " << period << std::endl;
    }
    else
    {
      bins[bin].add(chunk.head(chunkPos));
  
      // 2. use maps and strides to fill a mess of bins
      
      size_t endBin = bin;
      size_t chunkTailPos = 0;
      size_t chunkTailBin;
      do {
        offset = tail % period;
        bin = (bin + 1) % NUM_BINS;
  
        tail = binStart(bin + 1);
        size_t binLen = tail - offset;
        size_t periods = (chunk.size() - chunkPos) / period;
  
        size_t possibleEnd = chunkPos + periods * period;
        if (possibleEnd > chunkTailPos)
        {
          chunkTailPos = possibleEnd;
          chunkTailBin = bin;
        }
  
        if (periods > 0)
        {
          // note that this is a rectangular matrix, and statsaccumulator is expected to still process each element coefficient-wise
          // TODO: handle data passed with existing strides
          Eigen::Map<Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> const, 0, decltype(periodStride)> map(chunk.derived().data() + chunkPos, binLen, periods, periodStride);
          bins[bin].add(map);
        }
        else
        {
          break;
        }
  
        chunkPos += binLen;
      } while (bin != endBin);
  
      // 3. fill last partial bin from start, if needed
      if (chunk.size() > chunkTailPos)
      {
        bins[chunkTailBin].add(chunk.tail(chunk.size() - chunkTailPos));
      }
    }

    // find the extreme bins;
    Scalar minVariance = bins[0].variance();
    Scalar maxVariance = bins[0].variance();
    minBin = maxBin = 0;
    for (size_t i = 1; i < bins.size(); ++ i)
    {
      auto variance = bins[i].variance();
      if (variance < minVariance)
      {
        minVariance = variance;
        minBin = i;
      }
      if (variance > maxVariance)
      {
        maxVariance = variance;
        maxBin = i;
      }
    }
    //std::cerr << period << " min " << minBin << " sigma^2 " << minVariance << std::endl;
    //std::cerr << period << " max " << maxBin << " sigma^2 " << maxVariance << std::endl;
  }

  Scalar significance() const
  {
    return StatsDistributionSampling<Scalar, STATS_STANDARD_DEVIATION>(bins[minBin].fakeInfinitePopulation(), bins[minBin].size()).deviationSignificance(bins[maxBin]);
  }

#if 0
  // TODO oldNEXT: these overallError, overallRange functions are not quite right
  // as our signal _is noise_, we want to use the sampling distribution of the
  // standard error, using the maxima.
  Scalar overallError() const
  {
    return bins[0].groupError(bins);
  }

  // ooookay what did I mean by sampling distribution of the standard error?
  // well, our extreme bins (the ones with the mean magnitudes most far apart --
  // also the ones with the standard deviations most far apart) contain specifically
  // different sets of noise.
  // the minimum bin contains only the background noise.
  // the maximum bin contains the background noise + the signal noise.
  // The are roughly random selections from those two distributions.
  // I apparently would like to track the shapes of the distributions -- but for now
  // let's assume they're both gaussian.
  //
  // My stats accumulators are taking the raw samples from the two distributions.
  // So I don't expect the standard deviation of either of them to provide a sampling
  // distribution ....  do I?  why was I talking about a sampling distribution?
  //
  // hmm I think it has to do with the error of the masurement itself.  Our measure-
  // ment is one of the standard error and the mean, so if we report the error in our
  // measurement, then we'd report that error based on the sampling distribution.
  // Our measurement of the standard error and the mean is a sample from those
  // sampling distributions.
  //
  // So, I guess the question is -- what is reported for the range? the error or the
  // mean difference? and then, use the sampling distriution error to determine
  // the error of that meaurement.
  //
  // Our first approach is to take raw samples, so the data is signed and we expect
  // the mean to be zero.  The approach of magnitudes would give a nonzero mean.
  // It's notable that this magnitude mean is probably directly related to the
  // standard deviation, if the distribution is normal.
  //
  // Both distributions are supposed to have a mean of zero.
  // One distribution has a larger standard deviation, and its stddev is combined
  // from the other distribution.
  //
  // Our range I imagine would be what?  I guess it would be magnitude of the
  // standard devation of the added signal!  So we'll want to subtract somehow
  // the error from the background noise.
  //
  // var_sum = (n1 * var1 + n2 * var2) / (n1 + n2)
  //
  // // for dist 2, we have var_sum, and n2
  // // for dist 1, we have var, and n1
  // //
  // // that's a litle confusing and might be wrong, let's hash it out.
  // // yeah i think this was wrong
  //
  // dist_1: n_1 recordings of background noise (pop_1)
  // dist_3: n_3 recordings of pop_1 + pop_2
  //
  // var_samp1 is measured
  // var_3 should be var_pop1 + var_pop2, because it is a simple sum
  // (var_sum = (1 * var1 + 1 * var2) / (1 + 1)
  //
  // var_pop3 = var_pop1 + var_pop2
  //
  // we want var_pop2 = var_pop3 - var_pop1
  //
  // we don't have these population variances, though
  // we have sample variances ...
  // var_samp3 used n3 samples to sample the sampling distribution of the variance
  // var_samp1 used n1 samples to sample the sampling distribution of its variance
  //
  // these sampling distributions each have a variance of their own.
  //
  // the final measurement var2 = var_samp3 - var_samp1 is a sampling distribution
  // with a variance equal to the combined variance of the two that made it up.
  // var_sum = (var_samp1 * n1 + var_samp3 * n3) / (n1 + n3)
  //
  // that's the error we can return =)
  //
  // Okay, it actually works best with the current layout to make this general for
  // arbitrary kinds of distributions.  What am I doing here more generally?
  //
  // the range is the magnitude of the underlying distribution.
  // This is found by solving for the sum of two samples from two different distributions.
  // So what I need for a general function is the std deviation of the difference of two
  // samples from two distributions.
  // It might make sense to allow some simple arithmetic operations on distributions /
  // stats accumulators.  Thse operations would predict the distributions if the
  // provided operators were performed on samples from the distributions.
  //
  // I can also provide sampling distribution stats on the distributions themselves I
  // think ...
  //
  // The operator we'll need here is subtraction.  It could also be modeled as
  // multiplication by a scalar and addition.  That's a little more general.

  Scalar overallRange() const
  {

    Scalar min = Eigen::NumTraits<Scalar>::infinity();
    Scalar max = -min;
    for (auto const & accum : bins)
    {
      auto v = accum.average();
      if (v < min)
      {
        min = v;
      }
      if (v > max)
      {
        max = v;
      }
    }
    return max - min;
  }
#endif

private:
  size_t offset;
  std::array<StatsAccumulator, NUM_BINS> bins;
  size_t period;
  Eigen::OuterStride<Eigen::Dynamic> periodStride;

  size_t minBin;
  size_t maxBin;

  size_t binStart(size_t bin)
  {
    size_t inter = bin * period;
    return inter / NUM_BINS + (inter % NUM_BINS ? 1 : 0);
  }
};

// given some incoming data, accumulates data % a given period, passing it
// through a provided metric
// TODO: this makes a stats accumulator for each sample of the wave !!!!!
//       just a single mean is all that is needed (for much better memory use); just one statsaccumulator can be kept for whole wave.
template <typename StatsAccumulator>
class PeriodModelAccumulatorPreciseSignal
{
public:
  using Detail = StatsAccumulator;
  using Scalar = typename StatsAccumulator::Scalar;

  PeriodModelAccumulatorPreciseSignal(size_t period, Detail const & accum = {})
  : offset(0),
    periodWave(period, accum),
    periodStride(period)
  { }

  template <typename Derived>
  void add(Eigen::PlainObjectBase<Derived> const & chunk)
  {
    // okay, so chunk[0] goes in data[0] and on and on, but modulo data.size()
    
    size_t subcount = chunk.size() / periodWave.size();
    size_t subremainder = chunk.size() % periodWave.size();

    size_t srcIdx = 0;
    while (srcIdx < periodWave.size())
    {
      // make a vector of all the instances of this sample idx of the wave
      // TODO: properly handle storage order and stride of chunk
      Eigen::Map<Derived const, 0, decltype(periodStride)> map(chunk.derived().data() + srcIdx, subcount, periodStride);

      // process it
      periodWave[offset].add(map);

      // increment index in wave
      ++ srcIdx;
      offset = (offset + 1) % periodWave.size();
      if (subremainder <= 0)
      {
        subremainder += periodWave.size();
        -- subcount;
      }
      -- subremainder;
    }
  }

  Scalar overallError() const
  {
    return periodWave[0].groupError(periodWave);
  }

  Scalar overallRange() const
  {
    Scalar min = Eigen::NumTraits<Scalar>::infinity();
    Scalar max = -min;
    for (auto const & accum : periodWave)
    {
      auto v = accum.average();
      if (v < min)
      {
        min = v;
      }
      if (v > max)
      {
        max = v;
      }
    }
    return max - min;
  }

  template <typename Derived1, typename Derived2>
  void wave(Eigen::MatrixBase<Derived1> const & _values, Eigen::MatrixBase<Derived2> const & _errors)
  {
    Eigen::MatrixBase<Derived1> & values = const_cast<Eigen::MatrixBase<Derived1> &>(_values);
    Eigen::MatrixBase<Derived2> & errors = const_cast<Eigen::MatrixBase<Derived2> &>(_errors);

    values.derived().resize(periodWave.size());
    errors.derived().resize(periodWave.size());
    for (size_t i = 0; i < periodWave.size(); ++ i)
    {
      values[i] = periodWave[i].average();
      errors[i] = periodWave[i].errorInAverage();
    }
  }

private:
  size_t offset;
  std::vector<StatsAccumulator> periodWave;
  Eigen::InnerStride<Eigen::Dynamic> periodStride;
};


// PeriodFinder is a class that takes an incrementally-provided set of samples, and a range of possible periods, and identifies the strongest period in the data with that range
template <typename PeriodModelAccumulator>
class PeriodFinderAccumulateAllInteger
{
  // accepts buffers of samples
  // finds a periodic signal in them
public:
  using Scalar = typename PeriodModelAccumulator::Scalar;

  PeriodFinderAccumulateAllInteger(size_t minPeriod, size_t maxPeriod, typename PeriodModelAccumulator::Detail const & config)
  : minPeriod(minPeriod), _minPeriodsRead(0), _minPeriodsRemainder(0)
  {
    periodMetrics.reserve(maxPeriod - minPeriod);
    for (size_t i = minPeriod; i <= maxPeriod; ++ i)
    {
      periodMetrics.emplace_back(i, config);
    }
    std::cerr << "Constructed." << std::endl;
  }

  template <typename Derived>
  void add(Eigen::PlainObjectBase<Derived> const & chunk)
  {
    Scalar bestSignificance = 1.0;

    for (size_t pidx = 0; pidx < periodMetrics.size(); ++ pidx)
    {
      auto & pa = periodMetrics[pidx];
      pa.add(chunk);
      Scalar significance = pa.significance();
      std::cerr << "period #" << pidx << " (" << pidx + minPeriod << "): " << significance*100 << std::endl;
      if (significance < bestSignificance)
      {
        bestSignificance = significance;
        _bestPeriod = pidx + minPeriod;
        _bestSignificance = significance;
      }
    }

    _minPeriodsRemainder += chunk.size();
    _minPeriodsRead += _minPeriodsRemainder / minPeriod;
    _minPeriodsRemainder %= minPeriod;
  }

  size_t periodsRead() { return _minPeriodsRead; }

  size_t bestPeriod() { return _bestPeriod; }
  Scalar bestSignificance() { return _bestSignificance; }

private:
  size_t minPeriod;
  std::vector<PeriodModelAccumulator> periodMetrics;

  size_t _bestPeriod;
  Scalar _bestSignificance;

  size_t _minPeriodsRead;
  size_t _minPeriodsRemainder;
};

// status:
//
// the noise is immediately detectable by a reliable small change in the standard deviation or variance
// of the absolute magnitude of the samples.
// this is perhaps, a metric (notably one that works).
//
// Case 1: we know the frequency of the oscillation, but not the phase.
//
//    Perhaps we could accumulate the samples in parallel, using the metric to place them in separate bins.
//    This means running the metric long enough that the difference is notable.
//
//    Once the difference is notable, the phase can be found by looking at where the bin change is.
//
//    This does not account for phase drift.  I imagine to deal with phase drift, the process will have to
//    be halted as soon as the difference is notable (another metric will be needed for notability, preferably
//    statistics-based).  It's probably possible to determine the drift from the shape of the bin accumulation.
//
// Case 2: We know an approximation of the oscillation frequency.
//
//    We could perhaps assume the approximation is correct, accumulate into bins as above, and then observe
//    the curve of the bin to infer how far we are from the real frequency, and adjust.
//
// Case 3: We do not know the oscillation frequency.
//
//    This is harder.
//        - We could autocorrelate the magnitude of the data.
//            The weakness of the signal makes this hard, but perhaps it could be compensated for by taking
//            the data in chunks (assuming a low frequency) or in strides (assuming a high frequency) to
//            accumulate a variety of possible signals.
//
//       I like the chunks and stride idea.
//
//       Autocorrelation:
//        We start at sample 0, and we don't know the period.
//        We begin filling a buffer that grows to some maximum period.
//        Each next sample, we consider possible periods by comparing with buffer?
//
//            hold up: this sin't autocorrelation.
//
//       Autocorrelation would be comparing the entire buffer with itself at integrally spaced offsets.
//       Each comparison would give a result.
//       The result of a successive autocorrellation would give spikes that are spaced at the real period.
//
//       That's n x n comparisons, n for each n samples.
//
//    Period enumeration:
//                      metric A: determines likely waveform value from a set of samples
//                        probably average
//                      metric B: rapidly determines likelihood of a possible period, updating based on
//                                incremental addition of samples
//                        probably use histogram, and compare histogram to ideal
//                            (a good square wave is a histogram where most frequent value is significantly different from mean value)
//                        could use metric A to tighten histogram
//                      metric C: determines accuracy of signal measurement (+- error)
//                        look up std dev of summed gaussians and use 99% confidence or somesuch
//
//                      1. determine approximate frequency:
//
//                        a: keep a vector of possible integral periods
//                        b: an iteration algorithm will be used similar to autocorrelation
//                            but as each possible shift for the autocorrelation is basically a period, the data is accumulated % the period length
//                        b: use metric B to fill vector with likelihoods of periods
//                           will need to use metric C to determine when the actual data is present
//                           (compare magnitude of signal to magnitude of noise)
//                        c: pick the best period
//
//                        is this roughly n x n?
//                        what's the difference between this with a vector as long as the data, and autocorrelation?
//
//                        period vector created: n allocations
//                        vector filled with likelihoods:
//                            autocorrelation: n x n comparisons, making n differences
//                            periods: n x n loads into a histogram, then n judgements of histograms
//
//                            so for autocorrelation, we take the difference from ourselves.
//                            for periods, we load each item into a histogram .... then spew out avg and std dev
//
//                            for one thing, we don't care about super long periods, so we don't need to use the whole autocorrelation.
//                            we only go up to m
//                            then with autocorr we have m * n comparisons
//                            we could hybrid autocorr and periods:
//                            periods, with variance:
//                              for each m, we determine difference from mean
//                              for period in 0 .. m:
//                                for sample1 in 0 .. m:
//                                  for sample2 in sample1 .. n by m:
//                                    mean <- sample
//                                  for sample2 in sample1 .. n by m:
//                                    delta = sample - mean
//                                    delta += delta * delta
//                                    variance <- delta
//                            autocorr:
//                              for each sample1 in 0 .. n:
//                                for each sample2 in 0 .. n:
//                                  out[sample1] += abs(sample1 - sample2);
//                            the issue with basic autocorrelation is that information from multiple periods is not combined together
//                            but that could be fixed with a hybrid approach
//
//                            so there's the correlation algorithm, and then there's the correlation arithmetic
//                              regarding the algorithm, I'd like to combine across periods, but what is lost in doing that?
//
//                              for each possible period, correlation makes a chart that shows the repeatedness
//                              it wouldn't be hard to take this % period len; same n x n time
//
//                              okay.
//
//                           no, this isn't quite right.
//                           autocorrelation takes all the data and makes 1 result
//                           there's no need for a % period length: it's all combined !
//
//                           12345678901234567890
//                             12345678901234567890 -> difference of 2 everywhere
//
//                           O(m) -> try each possible period
//                             O(n) -> determine mean
//                             O(n) -> determine squared difference from mean
//                           -> O(m * n)
//                             but note that the wave information we get only uses two samples for each check.
//
//                           12345678901234567890
//                           34
//                           56
//                           78
//                           90
//                           12
//                           34
//                           56
//                           78
//                           90
//                           here we're maintaining a wave almost
//                           that could be designed to be a free memory allocation, maybe
//                           
//                           O(m) -> try each possible period
//                            O(n) -> determine m means from m/n data each
//                            O(n) -> determine squared difference from each mean
//                                  note that this is more divisions and squareroots because there are more variances (m*n many instead of n many)
//                           -> O(m*n) but with a larger constant factor
//
//
//                                  
//
//                      2. determine precise frequency?
//                        we'll want to add some jitter to the periods to get them precise and adjust for small drifts
//                        not crucial at this time; maybe could come out from phase
//
//                      3. determine phase, given frequency is near
//
//                        approach B: assume each phase and try them incrementally
//                            since this is a square wave, can note variance for high and low duty
//                            pick the phase with lowest variance
//                            since phases are incremental, can add and remove data from stats one point at a time (O=n)
//                              this may require either a histogram approach or a simpler metric than variance
//
//                        approach A: make a vector for each phase; compare each one to ideal
//                            this is basically A without the incrementation or the assumption of square waveness (O=n*n)
//
//
//
//
// Okay.  Problems:
// - the period length is way greater than the buffer size.  The limiting factor is not the length of the data, but rather the number of samples in the period.
// - the period length is too long to allocate every sample for a wide range of periods.
//
// Approaches:
// - downsample the signal for short periods (or downsample part of the algorithm); usable but high-freq components suffer and some rare signals could be lost
//      not that meaningful to consider this a 'vulnerable' approach because signals could be designed to be hidden, as such signals could be designed to not repeat in the first place
//   - hybrid downsample!  consider pdf of downsample?
// - store pdf of whole signal?  look for non-gaussian!
//   - Q: how is periodicity retained?
// - find another way to evaluate the period without storing every sample
// - some more versatile, more general approach
//   - signal must be amplified in face of background noise.  advantage is being able to predict an attribute of the signal.
//      - too general, go more specific
//
//
// Evaluating high-resolution periods in the face of memory limitation.
//    Okay, we have a blazingly fast sample rate, but our signal period is low.  Note: this problem was solved to some degree for van-eck phreaking; I could review existing work.
//    Approach for signals with low-freq components: downsample
//    Approach for signals with high-freq components: store a subsection of the signal.  if nothing found, try a different subsection !
//
//    some way to do the whole high-res signal at once?
//
//        hmmmmmm the signal looks like noise.  but if we know the period, we can combine the noisy signals to increase the snr.
//        unfortunately, we don't know the period exactly.  even if we know it pretty exactly, it is prohibitively long.
//
//        I want to record a whole period's worth of data, then take the next period and compare it with that one to clean it up.
//        The major issue is memory.
//
//        What I _do_ with this recording is get a mean of the actual signal.  I can then subtract the signal from the noise to determine the amplitude of the signal.
//
//        What we're getting out of this phase is the amplitude of the signal, and the variance of the noise.
//
//        It's notable that the background noise is expected to be gaussian.  In my environment here in Green Bank, this is pretty likely, but in a dense environment, the
//        background noise would have a distribution matching the combination of the other ones in the background.  A nearby e.g. wifi hotspot probably does not emit a gaussian
//        distribution and would be loud enough not to produce a gaussian sum.
//
//                  okay, that's kind of interesting, cause I think we're already assuming gaussian background noise.
//
//          perhaps a more interesting metric regarding the signal is 'smoothness' ...
//          what attributes does the signal have?
//
//        A correct period will sample the signal sequentially, such that the wave produced matches the wave emitted.  Other signals will be sampled randomly; so the distribution of
//        the data may be the same, but the ordering of the data will be all messed up.
//
//        That's perhaps a key item: our ordering of the correct period is the only one that is correct.  When the signal is plotted with the correct ordering, each sample is likely
//        near the last.  There's also some attribute, I imagine, that is strong compared to randomness or gaussianness.
//
//         /\              .   .
//        /  \                . . .
//            \  /         ..  .
//             \/            .   .
//
//        One notable transformation is to frequency space.  If we randomize ordering, the fourier transform is likely to look like noise.  The fourier transform of a properly
//        ordered signal is likely to have strong peaks: high maxima and low minima.  Note that this transform can be streamed and downsampled.
//
//        Of course, we could downsample our existing wave storage, too, which also streams, and gives us direct per-sample stats.
//
//        Another notable metric is the delta between successive samples i.e. the first derivative and subsequent ones.  Lower for ordered data; higher for unordered data.
//
//        I'm partially leaning towards the concepts of taking downsampled pdfs.  The bit resolution gets higher as we downsample; and these higher-bit data can then be used to create
//        more detailed pdfs without needing to sum samples.  Summing samples may also be needed.
//
//        So, it's notable that we're solving the problem here of seeing in the dark, and in such situations data needs to be downsampled to become visible.  I'm focusing here
//        on downsampling across "time" (successive periods), but it's equally relevent to downsample across "space" (intra-period), and via "prediction" (period models).
//
//        Increase SNR:
//          - downsampling samplerate <==
//          - accumulating periods <==
//          - predicting changes <==
//          - retuning the tuner XX not that good because what we want is to discern frequency differences
//          - reconfigering the tuner XX minimal bang for the buck expected, unknown domain
//          - better antenna <== could be helpful
//          - better radio XX helpful but want to keep price minimal for others to reuse; not where i want to invest this kind of researchy work
//          - physically moving or changing things XX not that good because we want to discern spacial differences too
//
//          For this usecase it's actually about equally meaningful to downsample as it is to record for longer.
//          Note that my current approach _does_ increase the snr via "downsampling" but it stores data for every point which produces failure.
//
//          Some other major approaches include random sampling and sublengths (good for van eck perhaps), and period models (good for everything; some interesting square wave approaches)
//
//
//          There's something remaining about how to tune the approach specifically for a square wave.  I think I've written about it somewhere else in here ...
//
//          Hmm!  The first approach I thought of used the metric that I had already noticed worked.  The standard deviation of the absolute magnitude of the signal is a metric that
//          is indicative of the square wave being present.  This metric is nonpresent for the low half of the wave.
//          Since it's a square wave, it can be downsampled arbitrarily small -- e.g. binary tree -- and this metric used to narrow down the peak and value rapidly.
//
//          Okay, a significant missing piece is how to adjust the metric measuring the square wave.
//          It's a pulse of loud noise, so we expect it to be gaussian.  We can perhaps assume that the linear sum of two gaussian functions produces linearly summed means and variances, perhaps with variances losing the sign of the sum.
//          Two approaches:
//            1. could this be done without taking the absolute magnitude?
//            2. how would this be done if the absolute magnitude were taken?
//
//            Pretty sure I've figured these before, picked #2, and solved it, but I've misplaced that work.  #1 is enticing!
//
//            I think that #1 sounds good in theory but did not work in practice; uncertain.  Might be out of context.
//
//            With #1, the signal is noisy ....  it looks actually like it would work fine.  Even if adjacent samples are near each-other due to low-pass filtering or somesuch, it's
//            still just a gaussian distribution of information; there's jsut some random phase in there.  right ????????? uncertain but basically
//            (if i take the fourier transform I'll see a randomly bouncing noise floor ... so I guess the phases are likely kind of doing a brownian motion walk)
//            each sample is a complex number: can I still take the std deviation of complex samples? it sounds a little funny
//            an easy approach would be to drop the Q and deal only with a real signal.  Including the I doubles my bit-depth, though, if I want to treat that as real too.
//            Would treating I as real cause a problem?
//
//            HMM I'm having trouble finding avenues to find answers to these questions!
//
//            - air contains dense random frequencies, weak
//            - recorded as IQ, which mixes the raw physical data with a model wave and takes 0 deg and 90 deg as samples
//            - these samples basically represent frequency-shifted data
//            - if I take their fourier, I see the underlying frequencies shifted such that DC is the frequency of the model wave that was tuned to
//            - this fourier represents what is physically in the air: randomness with random phase, densely occuping all frequencies
//            - the high frequency behavior of the data will represent the behavior of the air that is distant from the model wave
//              this high frequency behavior will be present as densely as the DC behavior, with equal magnitude
//            - the phase behavior we see is the sum of all the present frequencies, and the frequencies are delineated in linear Hz ...
//            - so if the bandwidth is 20 MHz and the sample rate is 20 Hz, then ... this variation happens at 20 different frequencies between 0 MHz and 20 MHz, 1 at each MHz,
//              1 at each sample thingy
//            - so 1st sample to 2nd sample is 20 MHz
//            - 2nd sample to 3rd sample is 20 Mhz + 19 MHz
//            - 3rd sample to 4th sample is 20 MHz + 18 MHz ... ?
//
//            - okay! we can think about it smaller.  It's loggy, inversy.  If the fourier shows flat background noise, then the low frequencies will have a smaller impact than
//            - the high frequencies on the per-sample changes in magnitude or likely phase as well.  Their impact will be spread out over their whole period (in a sine wave).
//            - So, given the frequencies are evenly filled, we expect the samples to jitter around a ton; mostly from high frequencies, a little from low frequencies.
//              Because high freqs are likely as filld as low freqs, we don't expect nearby phases adjacent to each other (unless the data is low-pass filtered or somesuch).
//
//            - okay, so I think that complex samples could be considered normally distributed as complex numbers, kinda.  sounds like the frequency domain might be a better spot
//              for it or something.
//
//            - does the mean of a complex number work? yes
//            - does the variance of a complex number work? hmmm
//              squaring a complex number does not produce a positive real number, which is what the variance expects.
//              sum(delta^2) / (n - 1)
//              I guess we'd want to take the absolute magnitude of delta.
//              
//              that sounds like it would be certainly correct if we were dealing with a phase-locked recording, but I'm not sure if it's correct for this noise.
//              The mean of the noise is supposed to be zero.  So taking the variance would just be the absolute magnitude in that sum...
//              Well that's better than subtracting the mean of the absolute magnitude !!!
//
//              okay I think this sounds good, I'd just like to verify my variance of complex numbers somehow
//                    hmm okay taking the magnitude of the delta is intuitive
//                    what would happen if we kept it complex?
//                    better find out perhaps
//
//                    it's notable that when taking a square root of a complex number, or raising a complex number to a power,
//                    I think the magnitude of the number behaves as if it were real, whereas I think the angle of the number is summed or differenced.
//                    Hence taking the magnitude of a number before or after a root or power, is meaningless
//                          probably do it earlier for fewer calculations
//
//                    It's also notable that if we don't do that, we get a variance and std deviation with a phase angle: a complex number.
//                        What does this mean?
//
//                    The mean, since it just uses complex addition, it reasonable: the real portion is the mean of the real portions, and the imaginary portion is the same.
//                    
//
//              When we use complex numbers, we're taking samples with two values to them and treating the two values as one value.
//              These samples are not on a straight line; however, the concepts of variance and deviation kind of assume a straight line.
//              One can expand these concepts to pythagorean distance, for example, which if the mean is zero is basically the same as taking the magnitude of the variance.
//              
//              The datapoints are really two-valued, and the meaning here is likely magnitude and angle.  It's strange to think of it that way, though, because the mean of
//              the magnitude is nonzero, and the angle has no mean, really !  It makes more sense if the datapoint are taken counter to their meaning, as cartesian coords or
//              paired real-valued samples !
//                It's notable that these _are_ paired real-valued samples, physically, but in such a way as to have specific meaning.
//
//              It's also notable that my test regarding complex variance did not yield as good a variance as my test with real variance.
//                observations on test data:
//                taking both real and imag together as real resulted in a more accurate mean (twice as much data) and factor was indeed roughly 2 (0.4 - 2.7)
//                    this makes physical sense because the angles are random, so sin and cos are equally distributed across + and -
//                    hey that makes sense for our data too I think!
//
//             okay, for this randomly sampled data, the i and the q are the sin and cos of the angle of the signal piece thing ....
//                but they are taken at slightly different times so sample different portions of each wave
//                        they are 90 deg apart so their distance is 1/4 the sample width.
//
//             another notable thing is that the magnitude of a sin/cos pair is always the magnitude of the signal
//             whereas the magnitude of only the real portion is the magnitude of the cos of the signal ...
//
//             I think the safest approach would be take only the real portion of the data
//              then try doubling the speed by adding the imaginary portion and see if the results are the same.
//              A possible problem is that the imaginary portion contains roughly the same information as the real portion,
//              since samplerate*4 is probably filtered out.  It could be _important_ to only use the real portion.
//              If this is true, (could be likely), then adding the imaginary portion would result in falsely low variances.
//              Additionally, summing the points might not result in the wave being smoothed.
//              That could be the way to test.
//
//
//            it is also notable that the complex magnitude, and the sin or real portion, will have different distributions; they contain the frequency information
//            differently and other confusing stuff.
//
//            it is also possible to specifically compare these distributions, without assuming they are gaussian.  I believe I was planning this earlier.
//            This would allow for arbitrary background signals!
//
//            so, in _identifying_ the signal, shuold we use the magnitude or the values centered at zero?
//            these are kinda similar.
//            the mean of the mag is sum(val) / count(val)
//            the variance of the centred is sum(val^2) / count(val
//
//            The mean picks 50% above and below, so I would expect mean of gaussian magnitude to be 50% point.
//            Meanwhile the std deviation is like 80% point or something.  It sounds like a better choice !!
//            so let's not use mag for this detection
//
//            once we're taking these standard deviations, what metric do we use to determine we don't have noise?
//
//            We have 4 bins; each has a standard deviation.
//            The min bin and the max bin have different standard deviations.
//            We want to determine that we have actually sampled different distributions, which becomes intuitively obvious when the standard deviations differ by more than their
//            value, really.
//
//                  StdDev of n samples from a population vs StdDev of population
//                  turns out if the population is normally distributed, the sampld std deviation has a std deviation of sigma/sqrt(2 n) where sigma is the stddev of the population and n is the smple count.
//
//                  So if we assume the sampled std dev is roughly equal to the population std dev, we can estimate the error.  This is the standard error needed.
//
//
// techniques:
//  - downsample
//  - take pdf rather than value; compare attribute such as kurtosis
//  - take only a subsection of period, or randomly chosen points (quality!) <======================
//  - take fewer periods; just get a general idea of the period before narrowing in
//  - don't use doubles in period storage; use like uint16's
//      -> poor, only divides ram by 4 but requires lots of work
//  - abstract this processing into a pluggable component, perhaps a PeriodModel <==================
//
//  I like the randomly-chosen points approach.  This is used in the ICA implementations I've looked at more recently, and it would catch signals with both low-freq and high-freq comps.  Much more recording time might be needed because the number of samples compared would be fewer, but this time would probably be still fast for a human.
//
// ===============================================================================================================================================================================
//  TODO APPROACH:
//    A. abstract processing into a pluggable component PeriodModel
//        -> this means the code that chooses to accumulate each sample of each period choice [and use their error and range?]
//    B. take current work and store as a PeriodModel that stores every sample for the wave and assumes an underlying repeating signal
//    C. make a PeriodModel for square waves
//      - [X] make 4 stats bins
//        - [X] make a todo to make another class that uses a binary tree for adaptive subsampling
//      - [X] fill stats bins with incoming data and track std dev
//            -> I'll probably drop the q data and use real-valued data for now
//            - [X] make a todo to try using magnitude of complex variance, and to try twice as much data with i & q both considered real.  we'll want to see if the approaches using more data have accurate stats, I guess
//                -> this could also be measured by looking at all 3 metrics on the raw data.  we'll want the one that is most extreme and most reliable, I guess.
//      - [/] maxima is taken with std dev of results (the 2 opposite bins that are most extreme are picked)
//            -> this gives range
//      - [ ] sampling distribution of the standard error can be used to determine standard error output, using the maxima
// ==============================================================================================================================================================================

#include <unsupported/Eigen/FFT>

// okay, fft takes a min and max period I guess?
// 'period' is in # samples
//
// fft has # samples coming in.
// at DC, that's component where period length = sample size
// at far right, that's component where period length = 2, I think
// in between, the delineations I believe are in frequencies not in periods,
// so invertedly
//
// 0 hz at far left
// [SR] hz at far right
//      periods / sec
//      instead
//      periods / window
// 0 p/w at left
// count or count/2 p/w at right
// basically we can try to assume the p/w are the index into the fft result
//
// so I have a periodsize, the p/w are going to be periodsize / windowsize

template <typename T>
class PeriodFinderFFT
{
public:
  using Scalar = T;

  // TODO: downSampling should use sinc function to properly take a rectangular window of data
  //       can frequency shift by multiplying by complex oscillator to place window over period of interest

  PeriodFinderFFT(size_t minPeriod, size_t maxPeriod, size_t bufferSize, size_t downSampling = 1)
  : _minPeriodsPerBuffer(bufferSize * downSampling / maxPeriod),
    sampleBuffer(downSampling),
    sampleBufferSize(0),
    fftBuffer(bufferSize),
    //fftAccum(decltype(fftAccum)::Zero(bufferSize / minPeriod - _minPeriodsPerBuffer)),
    offset(0),
    downSampling(downSampling),
    buffersRead(0),
    downsampleStride(downSampling),
    foregroundDists(bufferSize * downSampling / minPeriod - _minPeriodsPerBuffer),
    maxIdx(0),
    maxIdx2(0)
  { }

  static size_t periodResolutionToBufferSize(size_t period, Scalar resolution = 1.0)
  {
    using Eigen::numext::ceil;
    return period + (size_t)ceil(period / resolution * period);
  }

  static Scalar bufferSizeToPeriodResolution(Scalar period, size_t bufferSize)
  {
    return period / (bufferSize / period - 1);
  }

  template <typename Derived>
  void add(Eigen::PlainObjectBase<Derived> const & chunk)
  {
    //if (downSampling != 1)
    //{
    //  throw std::logic_error("downSampling != 1 unimplemented");
    //}

    size_t chunkOffset = 0;

    static constexpr typename Eigen::FFT<Scalar>::Flag fftFlags = static_cast<typename Eigen::FFT<Scalar>::Flag>(Eigen::FFT<Scalar>::Unscaled | Eigen::FFT<Scalar>::HalfSpectrum);
    Eigen::FFT<Scalar> fft(Eigen::default_fft_impl<Scalar>(), fftFlags);

    // first empty sampleBuffer
    if (sampleBufferSize > 0)
    {
      if (sampleBufferSize + chunk.size() < downSampling)
      {
        sampleBuffer.segment(sampleBufferSize, chunk.size()) = chunk;
        sampleBufferSize += chunk.size();
        return;
      }

      chunkOffset = downSampling - sampleBufferSize;
      fftBuffer[offset ++] = sampleBuffer.head(sampleBufferSize).sum() + chunk.head(chunkOffset).sum();
    }

    // now calculate end of downsampling and refill sampleBuffer
    size_t numDownedSamples = (chunk.size() - chunkOffset) / downSampling;
    size_t chunkDownsamplingTail = chunkOffset + numDownedSamples * downSampling;
    sampleBufferSize = chunk.size() - chunkDownsamplingTail;
    sampleBuffer.head(sampleBufferSize) = chunk.segment(chunkDownsamplingTail, sampleBufferSize);

    // make matrix to downsample chunk
    Eigen::Map<Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const, 0, decltype(downsampleStride)> downsamplingMap(chunk.derived().data() + chunkOffset, numDownedSamples, downSampling, downsampleStride);

    // add chunk to buffer
    size_t amountToFillBuffer = fftBuffer.size() - offset;
    if (numDownedSamples < amountToFillBuffer)
    {
      fftBuffer.segment(offset, numDownedSamples) = downsamplingMap.rowwise().sum();
      offset += numDownedSamples;
      return;
    }
    else
    {
      fftBuffer.segment(offset, amountToFillBuffer) = downsamplingMap.topRows(amountToFillBuffer).rowwise().sum();
    }
    
    size_t downsampledOffset = amountToFillBuffer;
    size_t amountRemaining;

    while (true)
    {
      // fft buffer if filled
      fftCurrentResult = fftCurrentResult.Constant(fftBuffer.size()/2+1, 7);
      fft.fwd(fftCurrentResult, fftBuffer);
  
      // take abs of result and sum into accumulations
      
      // this could be sped up by allow dist accumulators to handle matrices of parallel distributions: it would just be a vector of sums

      // TODO: should we weigh the dists by the number of periods contained in the FFT
      // for them?  note that we often get a good result after just 1 fft, so n=1 for the
      // result freq.
      auto magData = (fftCurrentResult.array().segment(_minPeriodsPerBuffer, foregroundDists.size()).abs() / (fftBuffer.size() * downSampling)).eval();
      overallInterestDist.add(magData);
      for (size_t i = 0; i < foregroundDists.size(); ++ i)
      {
        foregroundDists[i].add(magData[i]);
      }

      ++ buffersRead;

      amountRemaining = numDownedSamples - downsampledOffset;
      if (amountRemaining < fftBuffer.size())
      {
        break;
      }

      fftBuffer = downsamplingMap.middleRows(downsampledOffset, fftBuffer.size()).rowwise().sum();
      downsampledOffset += fftBuffer.size();
    }

    // add last bit of chunk to buffer
    fftBuffer.head(amountRemaining) = downsamplingMap.bottomRows(amountRemaining).rowwise().sum();
    offset = amountRemaining;

    // calculate result
    maxVal = 0;
    for (size_t i = 1; i < foregroundDists.size() - 1; ++ i)
    {
      if (foregroundDists[i].mean() > maxVal)
      {
        maxVal = foregroundDists[i].mean();
        maxIdx = i;
      }
    }
    if (foregroundDists[maxIdx + 1].mean() > foregroundDists[maxIdx - 1].mean())
    {
      maxIdx2 = maxIdx + 1;
    }
    else
    {
      maxIdx2 = maxIdx - 1;
    }
    maxVal2 = foregroundDists[maxIdx2].mean();
  }

  Scalar bestMag() { return foregroundDists[maxIdx].mean(); }
  Scalar bestMagVariance()
  {
    decltype(overallInterestDist) backgroundDist = overallInterestDist;
    
    backgroundDist.remove(foregroundDists[maxIdx]);
    backgroundDist.remove(foregroundDists[maxIdx2]);

    backgroundDist.recalcForUnderlyingDataWithZeroMeanOfAbsoluteValue();

    return StatsDistributionSampling<Scalar, STATS_MEAN>(backgroundDist.fakeInfinitePopulation(), foregroundDists[maxIdx].size()).variance();
  }
  // wrt stats for bestMag, this is different.
  // When we find the frequency, we give the % likelihood it's random noise
  // But for the mag, we want a measure of the accuracy, not the chance of existing
  // If we were to measure multiple mags, we'd have a distribution.  We want to know how close we are to the mean of that distribution.
  // Oh!  hmm.  I guess we could give the standard deviation, or the variance.
  // I could also give some kind of confidence interval.
  // Given that we know the noise distribution, roughly, it should be reasonable.
  
  size_t periodsRead() { return buffersRead * _minPeriodsPerBuffer; }
  size_t bestPeriod() { return downSampling * fftBuffer.size() / (maxIdx + _minPeriodsPerBuffer); }
  size_t bestPeriod2() { return downSampling * fftBuffer.size() / (maxIdx2 + _minPeriodsPerBuffer); }
  size_t minBestPeriod() { return downSampling * fftBuffer.size() / (maxIdx + _minPeriodsPerBuffer + 0.5); }
  size_t maxBestPeriod() { return downSampling * fftBuffer.size() / (maxIdx + _minPeriodsPerBuffer - 0.5); }
  /*
   * how will we really get the best significance?
   * probably consider the non-best periods samples from some kind of population
   * look at the significance of our good period compared to those of the population
   *        ideally remove harmonics; so i guess a simple solution might be DC to pop / 2 + 1 wrt period
   *
   * okay, we have representative samples
   * we have abnormal sample ... really I guess it's many samples, huh
   * seems reasonable to consider it mean * n since it's a sum
   * 
   */
  Scalar bestSignificance() {

    decltype(overallInterestDist) backgroundDist = overallInterestDist;
    backgroundDist.remove(foregroundDists[maxIdx]);
    backgroundDist.remove(foregroundDists[maxIdx2]);
    
    /*if (maxIdx > 1)
    {
      backgroundDist.remove(foregroundDists[maxIdx - 1]);
    }
    if (maxIdx < foregroundDists.size() - 2)
    {
      backgroundDist.remove(foregroundDists[maxIdx + 1]);
    }*/

    return StatsDistributionSampling<Scalar, STATS_MEAN>(backgroundDist.fakeInfinitePopulation(), foregroundDists[maxIdx].size()).deviationSignificance(foregroundDists[maxIdx]);
  }
  Scalar bestSignificance2() {

    decltype(overallInterestDist) backgroundDist = overallInterestDist;
    backgroundDist.remove(foregroundDists[maxIdx]);
    backgroundDist.remove(foregroundDists[maxIdx2]);
    
    /*if (maxIdx > 1)
    {
      backgroundDist.remove(foregroundDists[maxIdx2 - 1]);
    }
    if (maxIdx < foregroundDists.size() - 2)
    {
      backgroundDist.remove(foregroundDists[maxIdx2 + 1]);
    }*/

    return StatsDistributionSampling<Scalar, STATS_MEAN>(backgroundDist.fakeInfinitePopulation(), foregroundDists[maxIdx2].size()).deviationSignificance(foregroundDists[maxIdx2]);
  }

  //HeapVector<typename Eigen::FFT<Scalar>::Scalar> const & significances() { return fftAccum; }

private:
  size_t _minPeriodsPerBuffer;
  HeapVector<Scalar> sampleBuffer;
  size_t sampleBufferSize;
  HeapVector<Scalar> fftBuffer;
  HeapVector<typename Eigen::FFT<Scalar>::Complex> fftCurrentResult;
  //HeapVector<typename Eigen::FFT<Scalar>::Scalar> fftAccum;
  size_t offset;
  size_t downSampling;
  size_t buffersRead;
  Eigen::OuterStride<Eigen::Dynamic> downsampleStride;

  StatsAccumulatorNormal<Scalar> overallInterestDist;
  std::vector<StatsAccumulatorNormal<Scalar>> foregroundDists;

  Scalar maxVal;
  size_t maxIdx;
  Scalar maxVal2;
  size_t maxIdx2;
};

#if 0
// Takes more FFT's but finds the period with constantly increasing precision
// Due to FFT, time is O(n log n).
// RAM used does not need to increase over time if guess error is appropriately
// decreased on a regular basis.
template <typename T>
class PeriodFinderResamplingFFT
{
public:
  using Scalar = T;
  using Complex = typename Eigen::FFT<Scalar>::Complex;
  using Real = typename Eigen::FFT<Scalar>::Scalar;

  PeriodFinderResamplingFFT(T guessFreq, T maxFreqError, T sampleRate)
  : guessHz(guessFreq),
    cropHz(maxFreqError),
    sampleRateHz(sampleRate),
    rawSampleCount(0)
  { }


  template <typename Derived>
  void add(Eigen::PlainObjectBase<Derived> const & chunk)
  {
    using Eigen::numext::exp;

    // 1. shift data with complex sinusoid with negative guess freq
    
    shiftedBuffer = chunk * exp(HeapVector<Complex>::LinSpaced(chunk.size(), rowSampleCount * EIGEN_PI * Complex(0, -2) * guessHz / sampleRateHz, (rowSampleCount + chunk.size()) * EIGEN_PI * Complex(0, -2) * guessHz / sampleRateHz));
    
    // DOWNSAMPLING BAND-PASS FILTER:
    // 2. fft the shifted data
    Eigen::FFT<Complex> fft(Eigen::default_fft_impl<Complex>(), Eigen::FFT<Complex>::Unscaled);
    fft.fwd(shiftedFFT, shiftedBuffer);

    // 3. crop the fft
    auto cropPeriodsPerChunk = cropHz/*periods/sec*/ * chunk.size()/*samples/chunk*/ / sampleRateHz/*samples/sec*/;

    // TODO: should we handle the nyquist frequency?  seems acceptable to ignore it
    shiftedFFT.segment(cropPeriodsPerChunk, cropPeriodsPerChunk) = shiftedFFT.tail(cropPeriodsPerChunk);
    // resize will be done by passing head() to ifft
    
    // 4. ifft the crop (much smaller !)
      // map ongoingData to write to end directly
    auto lenNewData = cropPeriodsPerChunk * 2;
    auto oldDataLen = ongoingData.size();
    ongoingData.resize(oldDataLen + lenNewData);
    Eigen::Map<HeapVector<Complex>> ongoingTail(ongoingData.data() + oldDataLen, lenNewData);
    fft.inv(ongoingTail, shiftedFFT.head(lenNewData));

    // [5. unshift data to produce accurate downsample: we don't do this as there is no need]
    

    // hmmm the approach I envisioned was to accumulate enough data that the signal shows,
    // and then to downsample further once the frequency is known better.
    //
    // We can't really downsample properly with our current lowpass filter, though, because
    // it involves fft'ing the entire incoming bit >_> right? unsure
    // fft the incoming buffer
    // if we're downsampling a hueg amount then the final chunk will just be one freq
    // we also have rectangular window artefacts =/
    //
    // my goodness
    //
    // but if I had a lowpass filter then I could pluck a sample at the proper spot.
    //
    // let's say I did have the lowpass filter, and could get downsampld incoming signals withotu window artefacts or needing to store a ton in a buffer ... that would remove all the fft's from the above code.  that soudns good!  I just need a filter with a proper width to it.  The lowpass filter is just how to downsample.

    // NOTE; the solution for now might be to just do a bigger FFT until the accuracy of the period is down to a single sample.  Can calculate the needed size, too.  Maybe no downsampling is needed.
    //
    // what periods are in an fft?
    // the indices are periods / buffer
    // np.fft.fftfreq gives periods / sample
    //
    // I want samples / period, so it would be buffersize / index
    // The largest period differences are at the low frequencies, naturally
    //
    // Note that if I frequency-shift something to DC, because the space is frequency rather
    // than time, the information I have on it has to do with the original frequency, not the DC frequency ... in period-space it is just 'stretched'
    
    // I seem to have wasted a lot of time redoing things I already did with a different view.  The csin frequency shifting was the same as the cropping I was doing in the original function because I used the same technique as a lowpass filter.  I was completely unaware of this duplicate work, and should note that state in the future to see if there is some way to wait, back out, and try soemthing that I have better awareness of.
    // but I do feel I have a better understanding of frequency-shifting of the FFT and the definition of the FFT with this bit of my mind.  these are things i have looked up in the past as well.
    
    // So the sample accuracy I have of a signal's period in the FFT has to do with the difference in period length of the adjacent integral neighboring periods/buffersize.

    // I should probably just accumulate the biggest buffer I can, increasing the FFT resolution ... I imagine ...

    // I'm confused as to my approach so I guess I should make it configurable.
    //
    // will the periodfinder spew out periods before it fills its first buffer?
    //      yes I guess it should, that should resolve that set of conflicts
    //
    // what other questions need to be answered to add this simple functionality to
    //   the existing FFT class?
    //
    // the addition is the concept of focusing on just adding more samples to identify
    // the frequency better.
    // the buffer size increase needed can be calculated.
    //
    // a second thingy shows me that the simplest way is just to pass a larger buffer
    // size to the existing class.  this has no versatility.  that's not necessarily
    // valuable but it is a path forward.  [one can assume that it is _faster_ although
    // it is hard to judge without being aware of the alternatives]
    //
    // it will be slower when i change the frequency of the oscillator.  i'd also
    // like to track multiple oscillators.  So we'll at least want to predict the
    // buffer size, or choose a large enough buffer size to idntify all of them
    //
    // let's just do the constant buffer size for now.
    //
    // hmmmmm should it let the user pass the buffer size?
    // or should it just pick a buffer size big enough to resolve the period to 1 sample?
    //
    // let's keep making this static fucntion that converts between the two
    //
    // ----
    //
    // okay, buffer size for period size
    // periodsize = samples / period
    // indices = periods / buffer
    // buffersize = samples / buffer
    //            
    // i have periodsize.
    // i'll need to convert this to an index to identify the resolution
    // index = buffersize / periodsize
    // in order to resolve 'resolution' differences from the period,
    // I'll need index +- 1 to be a periodsize within resolution of the original
    // buffersize / periodsize - 1 = buffersize / (periodsize + resolution)
    // buffersize - periodsize = buffersize periodsize / (periodsize + resolution)
    // buffersize (periodsize + resolution ) - periodsize (periodsize + resolution ) = buffersize periodsize
    // (periodsize + resolution ) - periodsize (periodsize + resolution ) / buffersize = periodsize
    // - periodsize (periodsize + resolution ) / buffersize = periodsize - periodsize - resolution
    // periodsize (periodsize + resolution ) / resolution = buffersize
    
    // note: because it's a square wave, it's okay to downsample by averaging.
    // the peak will stay flat

    // okay, to get an accurate square wave it came to many gigasamples
    // i feel this per-sample accuracy is helpful for resolving a highly attenuated source
    // let's go with it: but it is notable that we could use poor accuracy.  we would
    // have fewer ... is see; the amount of time needed for the recording would be
    // multiplied by the ratio of inaccuracy to a single sample.
    //
    // downsampling may be possible while still allowing us to get accurate periods ...
    // ... I think it is.  the fft catches the phase drift and includes it as a
    // component of the frequency.  it cares more about where the peaks are.
    //
    // 

    // with regard to downsampling ...  how will i accumulate the buffers?
    // the incoming data is oversampled, and downsampling it will leave a little remaining.
    // our downsampled data needs to fit into fft buffers, and a littl will be remaining at that time too
    // where do I currently store the remaining data?  I guess I need another buffer for the remaining un-downsampled data.
    //
    // I keep one FFT buffer and process it when it's filled.
    // I can keep doing that with downsampled data, but I'll need a little buffer to store overflow/remaining.

  }

private:
  T guessHz;
  T cropHz;
  T sampleRateHz;
  unsigned long long rawSampleCount;

  HeapVector<Complex> shiftedBuffer;
  HeapVector<Complex> shiftedFFT;
  HeapVector<Complex> downsampledShiftedBuffer;

  HeapVector<Complex> ongoingData;
  HeapVector<Real> ongoingFFT;
};
#endif

template <typename T>
static constexpr T erfinv(T x)
{
  using Eigen::numext::sqrt;
  using Eigen::numext::log;

  // from https://stackoverflow.com/a/40260471
  x = (1 - x)*(1 + x);        // x = 1 - x*x;
  T lnx = log(x);
  
  T tt1 = 2/(EIGEN_PI*0.147) + 0.5f * lnx;
  T tt2 = 1/(0.147) * lnx;
  
  return sqrt(-tt1 + sqrt(tt1*tt1 - tt2));
}

class PeriodFinderRunningAdjustment
{
  // re: adaptive subsampling, the first approach is to try bigger and smaller and pick the best
  //     an adaptive approach would try a range, and then eliminate part of the range and try a smaller range.
  //
  //     there are two sources of new information: acquiring more samples, and running more tests
  //     running more tests lets us narrow down precision
  //     acquiring more samples lets us increase confidence, and could decrease the number of tests needed
  //
  //     given the biggest problem with acquiring more samples is time, we want to spend as much as available on increasing tests, and then use more samples when they
  //     are available.
  //
  //     at the moment there is no interface for timing; this might require refactoring.  we just want to code the algorithm to allow for it.
  //     re: timing: I think the most flexible interface would be to provide for an interrupt function
  //     classes I suppose would derive from a base that would provide a function to check if interrupted. (we could even hook into the base and check the rtl buffer)
  //
  //     can we parameterize the approaches of this algorithm?
  //     chooseing what rises and falls to test: 
  //      - 1, then the next, the next
  //      - all at once, in parallel
  //      - random ones, adding more if time allows
  //
  //     choosing where to try the rises and falls
  //     : note that two rise/falls will be needed to adjust phase and period
  //      - 2 adjacent to each guess, 6 total guesses
  //      - all guesses ! pick the best
  //      - adaptively subsample extreme of possibilities and middle
  //
  //    so in the second challenge, we have a domain of possible peaks, and we're looking
  //    for the best one.  like a root finding problem.
  //    in the first challenge, we have a lot of data and limited time.
  //    maybe for now I'll test all the data at once in parallel? it seems mroe complicated
  //
  //    okay we have a buffer coming in: will we just test the next period, or all the
  //    periods of the buffer? (of course the buffer mgiht have no periods)
  //
  //    it would be _easier_ to implement testing just the next period; accumulate a buffer
  //    etc.  could I do this in a flexible way to expand to testing allt he periods at
  //    once?
  //
  //    it's notable that if I have a guess period, I can view the data in such a way
  //    that it is a matrix of rows, where each row is an instance of that period
  //
  //    okay ... and I did that already ... there's a parameterization here
  //
  //    my existing code tests a bunch of periods handed it
  //    I guess to make it more general I would separate out into 3 tihngs:
  //    - what dataset are we testing?
  //      -> 'ringer' concept tests 2 overlapping period's worth of buffers at once
  //      -> 'allinteger' concept tests every period in all buffers accumulated at once
  //      we're missing a concept of held data vs 'archived' data.  
  //      -> 'ringer' concept only needs the latest 2 periods, but could possibly work
  //         with more
  //      -> 'allinteger' concept discards buffers as soon as processed, but stores a
  //        ton of state
  //      maybe a good addition might be to allow the algorithm to control the buffer
  //      size it gets.  this would mean I don't have to handle rebuffering in every
  //      algorithm.  It would report back how much it 'ate' or 'needs' I suppose.
  //    - what kind of test are we doing?
  //      -> my previous test crafted 4 stats buckets for the period, and picked the most
  //         likely ones to be peaks and troughs.  it did not use the phase
  //      -> my next test will make 2 stats buckets, and it will determine both the phase
  //         and the period.
  //      sounds reasonable.  phase information does not need to be exposed;
  //        could be provided to constructor
  //    - how are we choosing which rises and falls to test?
};


template <typename Scalar>
class PeriodProcessorConvolveDownsample
{
public:
  PeriodProcessorConvolveDownsample(size_t minPeriod, size_t maxPeriod, size_t samplesPerCoarsestUnit)
  : _mips(std::log(samplesPerCoarsestUnit) / M_LN2 + 0.5),
    _samplesProcessed(0),
    _periodStats(1, (minPeriod + maxPeriod) / 2)
  {
    if (minPeriod <= maxPeriod / 2)
      throw "minPeriod must be more than half maxPeriod";

    _mips[0].minPeriod = minPeriod;
    _mips[0].maxPeriod = maxPeriod;
    _mips[0].modelSums = HeapVector<Scalar>::Zero(maxPeriod);

    // fill mip indices > 0
    for (size_t i = 1; i < _mips.size(); ++ i) {
      _mips[i].minPeriod = _mips[i-1].minPeriod / 2;
      _mips[i].maxPeriod = _mips[i-1].maxPeriod / 2 + 1;
      _mips[i].modelSums = HeapVector<Scalar>::Zero(_mips[i].maxPeriod);
    }
    if (_mips.back().maxPeriod < 4)
      throw "samplesPerCoarsestUnit too large";
    
    // allocate some ram for period tracking
    _onsetsCur.resize(4096);
    _periodsCur.resize(4096);
  }

  size_t bufferSizeMin() const
  {
    return _mips[0].maxPeriod * 2;
  }

  size_t bufferSizeMultiple() const
  {
    return 0;
  }

  size_t bufferSizeMax() const
  {
    return ~0;
  }

  template <class DataFeeder, class Derived>
  void process(DataFeeder & feeder, Eigen::DenseBase<Derived> const & chunk)
  {
    // 1. downsample to H *2
    _mips[0].chunk.conservativeResize(chunk.size());
    _mips[0].chunk = chunk; // this memory copy is unnecessary but makes code simpler
    for (size_t i = 1; i < _mips.size(); ++ i)
    {
      // each iteration sums adjacent elements, resulting in size/2
      // TODO: if we keep the phase aligned, then periods could be accumulated
      //       to detect signals orders of magnitude weaker while profiling, depending on accuracy of guess
      //       would want to form a map with stride == guess, and then sum the downsamples based on
      //       blocks of inequal sizes from the map.
      //       same complexity if each iteration sums from results of previous iteration, like this approach
      _mips[i].chunk = _downsample(_mips[i-1].chunk, _mips[i].chunk);
    }

    // 2. Convolve each possible period to map onsets (fully downsampled only)
    // TODO: consider statistical error in assuming best convolution is correct
    //        -> for now we assume it is loud enough that it will be
    size_t dI = _mips.size() - 1;
    size_t startRange = _mips[dI].maxPeriod - _mips[dI].minPeriod + 1;
    size_t window = _mips[dI].maxPeriod + startRange;

    // TODO: A second pass of the buffer after the onsets are gathered could improve accuracy
    // (might be easier just to seed the run with an existing profile file, then can re-run on data)
    //
    // FIX TODO: does temperature etc change results in successive recordings?
    
    long nextConvolveStart;
    long inbetweenOffset = 0;
    _onsetsCur.conservativeResize(0);
    if (!_onsets.size() || !_mips[dI].inbetween.size()) {
      _onsetsCur.conservativeResize(1);
      _onsetsCur[0] = 0;
      nextConvolveStart = _mips[dI].minPeriod;
      _mips[dI].modelCount = 1;
      _mips[dI].modelSums = _mips[dI].chunk.head(_mips[dI].maxPeriod);
    } else {
      // combine with data from previous buffer
      inbetweenOffset = _mips[dI].inbetween.size();
      nextConvolveStart = -inbetweenOffset;
      _mips[dI].inbetween.conservativeResize(inbetweenOffset + window);
      _mips[dI].inbetween.tail(window) = _mips[dI].chunk.head(window);

      do {
        long idx;
        idx = nextConvolveStart + inbetweenOffset;
        idx = convolveBestAndUpdate(dI, _mips[dI].inbetween, idx, startRange);
        _onsetsCur.conservativeResize(_onsetsCur.size() + 1);
        _onsetsCur[_onsetsCur.size() - 1] = idx - inbetweenOffset;
        nextConvolveStart = idx - inbetweenOffset + _mips[dI].minPeriod;
      } while (nextConvolveStart < 0);
    }

    // 2.5: do remaining periods in downsampled buffer
    // until latest end point passes end of chunk
    while (nextConvolveStart + window < _mips[dI].chunk.size()) {

      // given earliest start point of next period, convolve to identify it
      long bestIdx = convolveBestAndUpdate(dI, _mips[dI].chunk, nextConvolveStart, startRange);
      _onsetsCur.conservativeResize(_onsetsCur.size() + 1);
      _onsetsCur[_onsetsCur.size() - 1] = bestIdx;

      // update earliest start point
      nextConvolveStart = bestIdx + _mips[dI].minPeriod;
    }

    // populate _inbetween
    size_t extra = _mips[dI].chunk.size() - nextConvolveStart + 1;
    _mips[dI].inbetween.conservativeResize(extra);
    _mips[dI].inbetween = _mips[dI].chunk.tail(extra);
    

    // 3. Decrease dI and loop to get a map for raw buffer size
    // TODO: add functionality to combine adjacent periods to detect weak signals
    // TODO: consider adjusting downsampled bins to properly wrap upsampled bins as period adjust?
    while (dI) {
      using Eigen::numext::mini;
      using Eigen::numext::maxi;

      // 3-0: copy data up
      //    each onset will be doubled, then adjusted
      -- dI;
      startRange = _mips[dI].maxPeriod - _mips[dI].minPeriod + 1;
      window = _mips[dI].maxPeriod + startRange;

      size_t oI = 0;
      // 3-1: handle inbetween
      if (!_onsets.size()) {
        _mips[dI].modelCount = 0;
        _mips[dI].modelSums = HeapVector<Scalar>::Zero(_mips[dI].modelSums.size());
      } else {
        inbetweenOffset = _mips[dI].inbetween.size();
        _mips[dI].inbetween.conservativeResize(inbetweenOffset + window);
        _mips[dI].inbetween.tail(window) = _mips[dI].chunk.head(window);
        do {
          long idx = _onsetsCur[oI] * 2 + inbetweenOffset;
          long start = maxi(idx - 2, (long)0);
          long end = mini(idx + 2, (long)window + inbetweenOffset);
          idx = convolveBestAndUpdate(dI, _mips[dI].inbetween, start, end - start);
          _onsetsCur[oI] = idx - inbetweenOffset;
          ++ oI;
        } while(_onsetsCur[oI] < 2);
      }

      // 3-2: handle data in middle
      for (; oI < _onsetsCur.size(); ++ oI) {
        auto onset = _onsetsCur[oI] * 2;
        auto start = maxi(onset - 2, (decltype(onset))0);
        auto end = mini(onset + 2, (decltype(onset))_mips[dI].chunk.size());
        onset = convolveBestAndUpdate(dI, _mips[dI].chunk, start, end - start);
        _onsetsCur[oI] = onset;
      }

      // 3-3: populate next inbetween
      size_t extra = _mips[dI].chunk.size() - (_onsetsCur[_onsetsCur.size() - 1] + _mips[dI].minPeriod - 1);
      _mips[dI].inbetween.conservativeResize(extra);
      _mips[dI].inbetween = _mips[dI].chunk.tail(extra);
    }

    // process periods and onsets
    long lastOnset;
    size_t oI;
    if (_onsets.size()) {
      _periodsCur.conservativeResize(_onsetsCur.size());
      lastOnset = _onsets[_onsets.size() - 1] - _samplesProcessed;
      oI = 0;
    } else {
      _periodsCur.conservativeResize(_onsetsCur.size() - 1);
      lastOnset = _onsetsCur[0];
      oI = 1;
    }

    _onsets.conservativeResize(_onsets.size() + _onsetsCur.size());
    _onsets.tail(_onsetsCur.size()) = _onsetsCur.array() + _samplesProcessed;
    _samplesProcessed += chunk.size();

    for (size_t pI = 0; pI < _periodsCur.size(); ++ pI) {
      _periodsCur[pI] = _onsetsCur[oI] - lastOnset;
      lastOnset = _onsetsCur[oI];
      ++ oI;
    }

    // TODO FIX: the mean produced by this accumulator is wrong at this precision
    _periodStats.add(_periodsCur);

    if (_periodStats.mean() != bestPeriod()) {
      throw std::logic_error("math inaccuracy");
    }
  }

  HeapVector<size_t> const & onsets() const
  {
    return _onsets;
  }

  size_t periodsRead() const
  {
    return _onsets.size();
  }

  Scalar bestPeriod() const
  {
    //return Scalar(_onsets[_onsets.size() - 1]) / (periodsRead() - 1);
    return _periodStats.mean();
  }

  Scalar bestFrequency() const
  {
    return (periodsRead() - 1) / Scalar(_onsets[_onsets.size() - 1]);
  }

  size_t periodVariance() const
  {
    return _periodStats.variance();
  }

private:
  // structure representing successfully lower resolutions of data
  struct _mip {
    // period ranges at this resolution
    size_t minPeriod;
    size_t maxPeriod;

    // data currently worked
    HeapVector<Scalar> chunk;

    // data retained from last work (last buffer, last chunk)
    HeapVector<Scalar> inbetween;

    // wave model: sum of wave data found at this resolution
    // TODO: make vector stats accumulator to allow me to reference the vector of means but also use
    //       stats methods on it
    // TODO: use histogram stats for detailed info
    // TODO: update stats work to use boost stats libs
    HeapVector<Scalar> modelSums;
    size_t modelCount;
  };

  // TODO: data could be retained in one big matrix to increase memory locality?
  // eigen has the concept of triangular matrices
  std::vector<struct _mip> _mips;

  int64_t _samplesProcessed;
  HeapVector<long> _onsetsCur;
  HeapVector<int64_t> _onsets;

  HeapVector<Scalar> _periodsCur;
  StatsAccumulatorHistogram<Scalar> _periodStats;

  static auto _downsample(HeapVector<Scalar> & prev, HeapVector<Scalar> & next) {
    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 2, Eigen::RowMajor> const, 0, Eigen::OuterStride<2>>
      prevMap(prev.data(), prev.size() / 2, 2);
    next.conservativeResize(prevMap.rows());
    return prevMap.rowwise().sum();
  }

  size_t convolveBestAndUpdate(size_t dI, HeapVector<Scalar> & data, size_t start, size_t distance)
  {
    Scalar min = Eigen::NumTraits<Scalar>::infinity();
    size_t idx;
    size_t length = _mips[dI].modelSums.size();

    // TODO: NOT THREAD SAFE
    static HeapVector<Scalar> dataInterim;
    dataInterim.conservativeResize(distance + length);
    dataInterim = data.segment(start, distance + length) * _mips[dI].modelCount;

    for (size_t i = 0; i < distance; ++ i) {
      Scalar sum = (_mips[dI].modelSums.array() - dataInterim.segment(i, length).array()).abs().sum();

      if (sum < min) {
        min = sum;
        idx = i;
      }
    }
    _mips[dI].modelSums += data.segment(idx + start, length);
    ++ _mips[dI].modelCount;
    return idx + start;
  }
};

/*
 * missing pieces:
 * - detector
 * - frequency information (can include in detector)
 * - log ;/
 */

// calcs correlation
// this can be efficiently done for entire signal using fft, but will likely weigh loud noise peaks too highly; use of absolute difference would be more robust
template <typename Scalar>
class FFTCorrelationFunctor
{
public:
  FFTCorrelationFunctor()
  : _fft(Eigen::default_fft_impl<Scalar>(), typename Eigen::FFT<Scalar>::Flag(Eigen::FFT<Scalar>::Unscaled | Eigen::FFT<Scalar>::HalfSpectrum))
  { }

  template <class Derived>
  int operator()(Eigen::MatrixBase<Derived> const & first, Eigen::MatrixBase<Derived> const & second)
  {
    using Eigen::numext::mini;
    //_fftResult1.conservativeResize(first.size());
    //_fftResult2.conservativeResize(second.size());
    size_t nfft = mini(first.size(), second.size());
    _fft.fwd(_fftResult1, first, nfft);
    _fft.fwd(_fftResult2, second, nfft);
    // conjugate produces correlation instead of convolution
    _fftResult2.array() *= _fftResult1.array().conjugate(); // pointwise product into _fftResult2
    _fft.inv(_correlation, _fftResult2, nfft); // time-domain into _correlation

    int idx;
    _correlation.maxCoeff(&idx);
    if (idx * 2 >= _correlation.size())
    {
      return idx - _correlation.size();
    }
    else
    {
      return idx;
    }
  }

private:
  Eigen::FFT<Scalar> _fft;
  HeapVector<typename Eigen::FFT<Scalar>::Complex> _fftResult1;
  HeapVector<typename Eigen::FFT<Scalar>::Complex> _fftResult2;
  HeapVector<Scalar> _correlation;
};

// TODO: downsample convolution!  is any data precalculation needed?

template <class _StatsAccumulator, StatsStatistic STATISTIC = STATS_MEAN>
class WaveformModelStatsAccumulator
{
  // TODO: should use a parallel stats accumulator so all can be done at once
  // eases construction too
public:
  using StatsAccumulator = _StatsAccumulator;

  WaveformModelStatsAccumulator(StatsAccumulator const & initial, size_t maxlen)
  : _model(maxlen, initial),
    _environment(initial),
    _rotatedModel(maxlen),
    _shortest(maxlen),
    _wholePeriodsConsidered(0),
    _samplesConsidered(0)
  {
    for (size_t i = 0; i < _model.size(); ++ i)
    {
      _rotatedModel[i] = &_model[i];
    } 
  }

  WaveformModelStatsAccumulator(WaveformModelStatsAccumulator const & other)
  : _model(other._model),
    _environment(other._environment),
    _rotatedModel(other._rotatedModel),
    _shortest(other._shortest),
    _wholePeriodsConsidered(other._wholePeriodsConsidered),
    _samplesConsidered(other._samplesConsidered)
  {
    for (size_t i = 0; i < _model.size(); ++ i)
    {
      _rotatedModel[i] = _model.data() + (_rotatedModel[i] - other._model.data());
    } 
  }

  // assume the first period starts at onsets[0]
  // and the last period goes until data.size()
  // extraTail is the number of samples missing from the last period off the end
  template <class Derived1, class Derived2>
  void add(Eigen::DenseBase<Derived1> const & data, Eigen::DenseBase<Derived2> const & onsets, size_t extraTail = 0)
  {
    // TODO: all data is squared and summed twice by this line
    //       add functionality to add data to multiple accumulators at once
    _environment.add(data);

    // create a matrix from data
    // take columns over matrix, place in stats accumulators
    // TODO: could process the data inside the data array with appropriate constructs
    
    // TODO: handle <2 onsets: length calculation is different

    _matbuf.resize(onsets.size(), _model.size());

    // handle last onset, may be missing end bits
    auto lenTail = data.size() - onsets[onsets.size() - 1];
    if (lenTail + extraTail > _model.size())
    {
      throw std::runtime_error("period is longer than initialized max");
    }
    _matbuf.row(onsets.size() - 1).head(lenTail) = data.tail(lenTail);

    // handle first onset, may be missing start bits
    auto startHead = onsets[0] < 0 ? 0 : onsets[0];
    auto extraHead = startHead - onsets[0];
    auto lenHead = onsets[1] - onsets[0] - extraHead;
    _shortest = lenHead + extraHead;
    if (_shortest > _model.size())
    {
      throw std::runtime_error("period is longer than initialized max");
    }
    _matbuf.row(0).segment(extraHead, lenHead) = data.segment(startHead, lenHead);

    // TODO: include samples located >= _shortest
    //    - could add them manually, or after-the-fact
    //    - could sort matrix by row length (not too hard)
    //        -> could require some extra copies if many different row lengths
    //    - could have multiple matrices (memory allocation confusing in general sense, but really only 2 different row lengths will be passed atm)
    //    - could precalc set of onset lengths by row, sort it, then write to matrix by it ...

    for (size_t i = 1; i < onsets.size() - 1; ++ i)
    {
      auto len = onsets[i+1] - onsets[i];
      if (len > _model.size())
      {
        throw std::runtime_error("period is longer than initialized max");
      }
      if (len < _shortest)
      {
        _shortest = len;
      }
      _matbuf.row(i).head(len) = data.segment(onsets[i], len);
    }

    _wholePeriodsConsidered = _rotatedModel[0]->size() + _matbuf.rows();

    for (size_t j = 0; j < _shortest; ++ j)
    {
      auto start = j < extraHead ? 1 : 0;
      auto end = _matbuf.rows() - (j < lenTail ? 0 : 1);
      auto rowct = end - start;
      _rotatedModel[j]->add(_matbuf.col(j).segment(start, rowct));
      _samplesConsidered += rowct;
      if (_rotatedModel[j]->size() < _wholePeriodsConsidered)
      {
        _wholePeriodsConsidered = _rotatedModel[j]->size();
      }
    }

  }

  void add(WaveformModelStatsAccumulator<StatsAccumulator, STATISTIC> const & other)
  {
    using Eigen::numext::mini;
    size_t len = mini(other._rotatedModel.size(), _rotatedModel.size());
   
    for (size_t i = 0; i < len; ++ i)
    {
      _rotatedModel[i]->add(*other._rotatedModel[i]);
    }

    _wholePeriodsConsidered += other.periodsConsidered();
    _samplesConsidered += other.samplesConsidered();
    _shortest = other._shortest;
  }

  void clear()
  {
    _shortest = _model.size();
    for (size_t i = 0; i < _model.size(); ++ i)
    {
      _model[i].clear();
      _rotatedModel[i] = &_model[i];
    } 

    _wholePeriodsConsidered = 0;
    _samplesConsidered = 0;
  }

  // metric for number of periods used for stats calculation
  // approximate
  inline uint64_t periodsConsidered() const
  {
    return _wholePeriodsConsidered;
  }

  // number of samples used for stats calculation
  inline uint64_t samplesConsidered() const
  {
    return _samplesConsidered;
  }

  //inline Scalar meanPeriod() const
  //{
  //  return Scalar(periodsRead) / samplesRead;
  //}

  inline size_t maxUsefulPeriod() const
  {
    return _model.size();
  }

  template <StatsStatistic STAT, class Derived>
  void get(Eigen::DenseBase<Derived> const & data)
  {
    auto & rw_data = const_cast<Eigen::DenseBase<Derived> &>(data);
    rw_data.derived().conservativeResize(_shortest);
    for (size_t i = 0; i < _shortest; ++ i)
    {
      rw_data[i] = _rotatedModel[i]->template get<STAT>();
    }
  }

  template <class Derived>
  void get(Eigen::DenseBase<Derived> const & data)
  {
    get<STATISTIC>(data);
  }

  void rotate(int amount)
  {
    _mergeOverflow();

    while (amount > (int)_shortest)
    {
      std::cerr << "WARNING: drift is greater than 1 period" << std::endl;
      amount -= _shortest;
    }
    while (-amount > (int)_shortest)
    {
      std::cerr << "WARNING: drift is greater than 1 period" << std::endl;
      amount += _shortest;
    }

    std::rotate(
        _rotatedModel.begin(),
        _rotatedModel.begin() + (amount < 0 ? 0 : _shortest) - amount,
        _rotatedModel.begin() + _shortest);
  }

  template <StatsStatistic STAT>
  Scalar significance() const
  {
    // TODO: learn stats and review, see if this is correct, including underlying implementation
    
    Scalar overallSignificance = 1.0;
    for (size_t i = 0; i < _model.size(); ++ i)
    {
      if (_model[i].size() > 0)
      {
        overallSignificance *= StatsDistributionSampling<Scalar, STAT>(_environment.fakeInfinitePopulation(), _model[i].size()).deviationSignificance(_model[i]);
      }
    }

    return overallSignificance;
  }

  Scalar significance() const
  {
    return significance<STATISTIC>();
  }

  template <StatsStatistic STAT>
  Scalar integrate() const
  {
    // FIX TODO: assume least extreme signal element is background noise
    //           then, assume any that don't cross significance threshhold of being in the same distribution,
    //            are also background noise
    // Least extreme signal element is treated as background noise and removed
    // from each distribution sample.
    // Results represent difference of integral of signal from integral of minimum,
    // which scales proportionally to just plain signal integral.
    //
    // FIX TODO: review whether this is robust in the face of noise changing the smallest element

    using Eigen::numext::abs;

    auto min = real_limits<Scalar>::infinity();
    decltype(min) max = 0;
    size_t minIdx, maxIdx;
    for (size_t i = 0; i < _model.size(); ++ i)
    {
      if (_model[i].size() > 0)
      {
        auto stat = abs(_model[i].template get<STAT>());
        if (stat < min)
        {
          min = stat;
          minIdx = i;
        }
        if (stat > max)
        {
          max = stat;
          maxIdx = i;
        }
      }
    }

    // FIX TODO: make this heuristic make sense: e.g. only look at dists with variance and mean within half a std dev of sampling dist expected.  perhaps calc numbers to get 1 wrongly chosen dist on avg
    auto thresh = (max - min) * 0.125 + min;
    StatsAccumulator backgroundGuess = _model[minIdx];
    for (size_t i = 0; i < _model.size(); ++ i)
    {
      if (i == minIdx || _model[i].size() == 0)
      {
        continue;
      }
      
      if (_model[i].template get<STAT>() < thresh)
      {
        backgroundGuess.add((StatsAccumulator&)_model[i]);
      }
    }

    Scalar ret = 0;
    StatsDistributionBrief<Scalar> dist(0, 0);
    for (size_t i = 0; i < _model.size(); ++ i)
    {
      if (_model[i].size() > 0)
      {
        dist.setMeasures(_model[i].mean() - backgroundGuess.mean(), _model[i].variance() - backgroundGuess.variance());
        if (dist.variance() < 0)
        {
          continue;
        }
        ret += dist.template get<STAT>();
      }
    }
    return ret;
  }

  Scalar integrate() const
  {
    return integrate<STATISTIC>();
  }

  template <StatsStatistic STAT>
  Scalar integrationVariance()
  {
    // TODO: I assume variances sum here; i think that's correct?
    // TODO: check stats, and stats implementation



    Scalar ret = 0.0;
    _mergeOverflow();
    for (size_t i = 0; i < _shortest; ++ i)
    {
      if (_model[i].size() > 0)
      {
        ret += StatsDistributionSampling<Scalar, STAT>(_model[i].fakeInfinitePopulation(), _model[i].size()).variance() * 2;
      }
    }

    return ret;
  }

  Scalar integrationVariance()
  {
    return integrationVariance<STATISTIC>();
  }

  void dbgDump(std::ostream &f)
  {
    _mergeOverflow();
    uint32_t histCount = _shortest;
    f.write((char*)&histCount, sizeof(histCount));
    for (uint32_t i = 0; i < histCount; ++ i)
    {
      Scalar histStart, binWidth;
      uint32_t binCount;
      HeapVector<Scalar> const & hist = _rotatedModel[i]->histogram(histStart, binWidth);
      binCount = hist.size();
      f.write((char*)&i, sizeof(i));
      f.write((char*)&histStart, sizeof(histStart));
      f.write((char*)&binWidth, sizeof(binWidth));
      f.write((char*)&binCount, sizeof(binCount));
      for (size_t j = 0; j < binCount; ++ j)
      {
        f.write((char*)&hist[j], sizeof(hist[j]));
      }
    }
  }

private:
  std::vector<StatsAccumulator> _model;
  StatsAccumulator _environment;
  std::vector<StatsAccumulator*> _rotatedModel;
  size_t _shortest;
  Eigen::Matrix<typename StatsAccumulator::Scalar, Eigen::Dynamic, Eigen::Dynamic> _matbuf;

  uint64_t _wholePeriodsConsidered;
  uint64_t _samplesConsidered;

  void _mergeOverflow()
  {
    // assume components are low-frequency enough that it is fine to merge overflow with the data
    // #include "recbufmeta.hpp"
    //        if high frequency components are present, inaccuracy as true sample rate fraction approaches 0.5
    for (size_t i = _shortest; i < _rotatedModel.size(); ++ i)
    {
      _rotatedModel[i - _shortest]->add(*_rotatedModel[i]);
      _rotatedModel[i]->clear();
    }
  }
};

/*
 * grumble gripe
 * i'm storing a log of the shape of the wave at each frequency tuned to, at the
 * recording samplerate
 * i'm guessing the system is crashing because, when tuning to too many frequencies,
 * memory is exhausted
 *    (although that doesn't describe everything)
 * -> thinking a nice quick solution could be to provide a max # of freqs to store
 *    and start rotating them when that is reached
 *
 * approaches to resolving ram exhaustion:
 *  - more ram
 *  - optimize stats information (switch to gaussian summary? store data efficiently?)
 *  - max # of freqs stored
 *  - don't store model data for old freqs at all
 *  - downsample storage
 */

template <typename Scalar, class WaveformModel, class CorrelationFunctorType>
class PeriodProcessorDetectKnownCorrelate
{
public:
  // minSamplesPerSignificance can be considered the time duration in samples when 1 error will be found on avg
  // wave is considered detected when samples / significance > minSamplesPerSignificance
  PeriodProcessorDetectKnownCorrelate(WaveformModel const & initial, Scalar approxWavelen, CorrelationFunctorType correlationFunctor, Scalar minSamplesPerSignificance)
  : _firstPhaseTimestamp(0),
    _firstPhasePeriod(0),
    _lastPhaseTimestamp(0),
    _lastPhasePeriod(0),
    _initial{initial,initial},
    _approxWavelen(approxWavelen),
    _correlationFunctor(correlationFunctor),
    _minSamplesPerSignificance(minSamplesPerSignificance),
    _onsetStartPeriods(0)
  { }

  size_t bufferSizeMin() const
  {
    return _approxWavelen * 3 + 0.5;
  }

  size_t bufferSizeMultiple() const
  {
    return 1;
  }

  size_t bufferSizeMax() const
  {
    return std::numeric_limits<size_t>::max();
  }

  template <class DataFeeder, class Derived>
  void process(DataFeeder & feeder, Eigen::DenseBase<Derived> const & chunk, RecBufMeta const & meta)
  {
    // stores starting onset time relative to chunk start
    int64_t lastOnset; // TODO?: made this an int64_t to provide for huge sample gaps, but that may be silly, didn't really think about it, or consider elsewhere
    uint64_t missedOnsets = 0;
    uint64_t lastPeriods;
    if (_onsetsBuf.size())
    {
      lastOnset = _onsetsBuf[_onsetsBuf.size() - 1] + _onsetStart - meta.sampleTime;
      lastPeriods = _onsetStartPeriods + _onsetsBuf.size() - 1;
    }
    else
    {
      // first data buffer processed!
      _onsetStart = meta.sampleTime;
      lastOnset = 0;
      lastPeriods = _onsetStartPeriods;
    }

    FrequencyState & state = (*_modelsByFrequency.emplace(std::make_pair(meta.freq, _initial)).first).second;
    if (_modelsByFrequency.size() == 1)
    {
      state.pinned = true;
      _canonicalModel = &state;
    }
    else if (!state.active.samplesConsidered() && !state.past.samplesConsidered())
    {
      state.pinned = false;
      _unpinnedModels.push_back(&state);
    }

    // is this correct if data is lost without tuning?
    // this only checks for dataloss, but not for tuning
    if (-lastOnset >= _approxWavelen)
    {
      // adjust lastOnset to account for data gap
      auto shiftCt = Eigen::numext::floor((meta.sampleTime - _onsetStart) / _approxWavelen);
      lastOnset = _onsetStart + static_cast<uint64_t>(shiftCt * _approxWavelen + 0.5) - meta.sampleTime;
      lastPeriods = _onsetStartPeriods + shiftCt;
      //missedOnsets = lastPeriods - _onsetsBuf.size();
    }

    // ensure _onsetsCt is large enough
    size_t onsetCt = Eigen::numext::ceil((chunk.size() - lastOnset) / _approxWavelen);
    if (onsetCt > _onsetsCt.size())
    {
      size_t old = _onsetsCt.size();
      _onsetsCt.conservativeResize(onsetCt);
      for (size_t j = old; j < _onsetsCt.size(); ++ j)
      {
        _onsetsCt[j] = j;
      }
    }

    // add stats for wave onset guesses

    // onsetsBuf stores time relative to start for current freq run
    uint64_t bufferOffset;
    if (state.active.samplesConsidered() == 0)
    {
      // first data buffer for this run
      _onsetStart = meta.sampleTime;
      _onsetStartPeriods = lastPeriods;
      _onsetsBuf.conservativeResize(onsetCt);
      bufferOffset = 0;
    }
    else
    {
      _onsetsBuf.conservativeResize(_onsetsBuf.size() + missedOnsets + onsetCt);
      bufferOffset = meta.sampleTime - _onsetStart;
    }
    // TODO: code updated to no longer require tracking onsets between runs, no need for .tail() here,
    //       but will need to track total number of periods encountered and replace _onsetsBuf.size() everywhere
    _onsetsBuf.tail(onsetCt) = bufferOffset + (_onsetsCt.array().head(onsetCt) * _approxWavelen).round().template cast<int>() + lastOnset;

    /*
     * Bug: when wavelength shrinks by a lot, we can lose a whole period, similar to tuning delay
     * 
     * bug concern: we track all onsets for the current freq in _onsetsBuf, so dropping a period could confuse the tracking
     * solution?: drop the period, but store its location anyway
     */

    // shift onsets to match buffer start and add data
    state.active.add(chunk, (_onsetsBuf.tail(onsetCt).array() - bufferOffset), lastOnset + onsetCt * _approxWavelen - chunk.size());

    //std::cerr << "SamplesPerSignificance: " << state.active.samplesConsidered() / state.active.significance() << std::endl;

    if (state.active.samplesConsidered() / state.active.significance() > _minSamplesPerSignificance)
    {
      // "phase" attributes represent the average over the measurement
      //      these are the timestamp and period at which the measured, aligned phase matches reality
      size_t effectivePhaseOnsetIdx = _onsetsBuf.size()/2;
      //uint64_t activePhaseTimestamp = _onsetsBuf[effectivePhaseOnsetIdx] + _onsetStart;
      uint64_t activePhaseTimestamp = static_cast<uint64_t>(Eigen::numext::round(effectivePhaseOnsetIdx * _approxWavelen)) + _onsetStart + lastOnset;
      //assert(activePhaseTimestamp == activePhaseTimestamp2);
      uint64_t activePhasePeriod = effectivePhaseOnsetIdx + _onsetStartPeriods;
      if (state.past.periodsConsidered())
      {
        // rotate active wave to be in-phase with past wave
        state.past.get(_pastBuf);
        state.active.get(_activeBuf);
        int rotation = _correlationFunctor(_activeBuf, _pastBuf);
        std::cerr << "Drift: " << -rotation/Scalar(activePhasePeriod - state.lastPhasePeriod) << std::endl;

        // rotate stored sum
        state.active.rotate(rotation);
        activePhaseTimestamp -= rotation;

        // update timing information
        state.lastPhaseTimestamp = activePhaseTimestamp;
        state.lastPhasePeriod = activePhasePeriod;

        _lastPhaseTimestamp = activePhaseTimestamp;
        _lastPhasePeriod = activePhasePeriod;

        // update wavelength measure
        _approxWavelen = (_lastPhaseTimestamp - _firstPhaseTimestamp) / Scalar(_lastPhasePeriod - _firstPhasePeriod);

        // adjust final onset (other onsets will be ignored)
        uint64_t finalOnsetPeriod = _onsetStartPeriods + _onsetsBuf.size() - 1;
        // move buffer such that it represents xtarting at the phase-aligned time
        _onsetStart = activePhaseTimestamp;
        _onsetStartPeriods = activePhasePeriod;
        _onsetsBuf.conservativeResize(finalOnsetPeriod - activePhasePeriod + 1);
        // set value
        _onsetsBuf[_onsetsBuf.size() - 1] = (finalOnsetPeriod - activePhasePeriod) * _approxWavelen + 0.5;

        // TODO: if this is a large rotation, assume it is calibration
        //            (store final sum separately; it will be blurry + inaccurate)
        // TODO: if rotation existed at all, reuse current buffer for
        //       next run, to increase signal strength?
        //        -> possibly not too meaningful unless rotation is large:
        //           gains made in ability to resolve signal may be matched by
        //           losses in precision of metrics from blur?

        // adjust other freqs
        if (&state == _canonicalModel)
        {
          for (auto & otherModel : _unpinnedModels)
          {
            auto correctedFirstPhaseTimestamp =
              (otherModel->firstPhasePeriod - _firstPhasePeriod) * (_lastPhaseTimestamp - _firstPhaseTimestamp) / (_lastPhasePeriod - _firstPhasePeriod) + _firstPhaseTimestamp;
            auto correctedLastPhaseTimestamp =
              (otherModel->lastPhasePeriod - _firstPhasePeriod) * (_lastPhaseTimestamp - _firstPhaseTimestamp) / (_lastPhasePeriod - _firstPhasePeriod) + _firstPhaseTimestamp;

            auto drift = (correctedFirstPhaseTimestamp - otherModel->firstPhaseTimestamp + correctedLastPhaseTimestamp - otherModel->lastPhaseTimestamp) / 2;

            otherModel->past.rotate(-drift);
            otherModel->firstPhaseTimestamp += drift;
            otherModel->lastPhaseTimestamp += drift;
            otherModel->pinned = true;
          }
          _unpinnedModels.clear();
        }
      }
      else
      {
        state.firstPhaseTimestamp = activePhaseTimestamp;
        state.firstPhasePeriod = activePhasePeriod;

        if (&state == _canonicalModel)
        {
          _firstPhaseTimestamp = activePhaseTimestamp;
          _firstPhasePeriod = activePhasePeriod;
        }
      }

      // DEBUG write model data to file
      //static std::ofstream f("test_statmodels.modeldump");
      //state.active.dbgDump(f);

      state.past.add(state.active);
      state.active.clear();
    }
  }

  uint64_t periodsConsidered() const
  {
    if (_lastPhasePeriod == 0)
    {
      return 0;
    }
    return _lastPhasePeriod - _firstPhasePeriod;
  }

  Scalar bestPeriod() const
  {
    return _approxWavelen;
  }

  Scalar bestFrequency() const
  {
    return Scalar(_lastPhasePeriod - _firstPhasePeriod) / (_lastPhaseTimestamp - _firstPhaseTimestamp);
  }

  Scalar bestMag(Scalar & freq) const
  {
    auto needle = _modelsByFrequency.lower_bound(freq);
    auto dist = needle->first - freq;
    -- needle;
    if (freq - needle->first > dist)
    {
      ++ needle;
    }
    return needle->second.past.integrate();
  }

  Scalar bestMagVariance(Scalar & freq)
  {
    auto needle = _modelsByFrequency.lower_bound(freq);
    auto dist = needle->first - freq;
    -- needle;
    if (freq - needle->first > dist)
    {
      ++ needle;
    }
    return needle->second.past.integrationVariance();
  }

  // TODO: track variance in wavelength measure
  //       (to precision available given significance of signal)

private:
  struct FrequencyState
  {
    WaveformModel active;

    // TODO: dump past models that need large convolution, so as to not have their blurriness change the
    //       correct variances and means for accurate power measurement
    //       could optionally collect them into a 'calibration' model
    //       maybe have a state change once calibrated that resets the models
    // TODO: use variance to determine 'large convolution' above
    WaveformModel past;

    int64_t firstPhaseTimestamp;
    int64_t firstPhasePeriod;
    int64_t lastPhaseTimestamp;
    int64_t lastPhasePeriod;

    bool pinned; // whether times are aligned with global time
  };

  std::map<Scalar,FrequencyState> _modelsByFrequency;
  std::vector<FrequencyState *> _unpinnedModels;
  FrequencyState * _canonicalModel;
  int64_t _firstPhaseTimestamp;
  int64_t _firstPhasePeriod;
  int64_t _lastPhaseTimestamp;
  int64_t _lastPhasePeriod;
  FrequencyState _initial;
  Scalar _approxWavelen; // FIX TODO: this is used for period timing and uses of it would need to be changed to use cryptographic timing to provide robustness in the face of noise with equal timing to detected source
  CorrelationFunctorType _correlationFunctor;
  Scalar _minSamplesPerSignificance;

  HeapVector<Scalar> _onsetsCt;
  uint64_t _onsetStart;
  uint64_t _onsetStartPeriods;
  HeapVector<int> _onsetsBuf;

  HeapVector<Scalar> _pastBuf;
  HeapVector<Scalar> _activeBuf;
};

#if 0
struct FileRef
{
  std::string format;
  std::string location;
};

#include <fstream>
#include "file.hpp"
class PeriodicProfile
{
public:
  PeriodicProfile(std::string fname)
  : _modelfile(new fstream(fname, ios_base::in|ios_base::out|ios_base::binary))
  {
    // TODO: store wave data, probably a file reference
    if (_modelfile.size())
    {
      uint64_t ver;
      _modelfile.readUInt(ver);
      if (!ver || ver > 1)
      {
        throw std::runtime_error("unsupported profile version");
      }
      _modelfile.readUInt(_filesStart);
      _modelfile.readReal(_periodMean);
      _modelfile.readReal(_periodPrecision);
      _modelfile.readReal(_periodVariance);
      _modelfile.readUInt(_periodCount);

      _modelfile.seek(_filesStart);
      decltype(files.size()) fileCt;
      _modelfile.readUInt(fileCt);
      _files.resize(fileCt);
      for (fileCt = 0; fileCt < _files.size(); ++ fileCt)
      {
        _modelfile.readString(_files[fileCt].format);
        _modelfile.readString(_files[fileCt].location);
      }
    }
    else
    {
      _modelfile.writeUInt(1); // ver
      auto filesStartPos = _modelfile.tell();
      _modelfile.writeUInt(_filesStart = 0);
      _modelfile.writeReal(_periodMean = 0);
      _modelfile.writeReal(_periodPrecision = 1.0/0);
      _modelfile.writeReal(_periodVariance = 0);
      _modelfile.writeUInt(_periodCount = 0);

      _filesStart = _modelfile.tell();
      _modelfile.seek(filesStartPos);
      _modelfile.writeUInt(_filesStart);
      _modelfile.seek(_filesStart);

      _modelfile.writeUInt(0); // filect
    }
  };

  class Recording
  {
  public:
    static char const * type() { return "PROFILE EVENT META"; }

    Recording(std::string fname)
    : _metafile(new fstream(fname, ios_base::in|ios_base::out|ios_base::binary))
    {
      if (!_metafile.size())
      {
        throw std::runtime_error("no logfile");
      }
      uint64_t ver;
      _metafile.readUInt(ver);
      if (!ver || ver > 1)
      {
        throw std::runtime_error("unsupported profile version");
      }
      _metafile.readUInt(_onsetsStart);
      _metafile.readString(_logfile.format);
      _metafile.readString(_logfile.location);
      _metafile.readString(_devStr);
      _metafile.readReal(_tunedHz);
      _metafile.readReal(_rateHz);

      _metafile.seek(_onsetsStart);
      decltype(_onsets.size()) onsetCt;
      _metafile.readUInt(onsetCt);
      _onsets.resize(onsetCt);
      for (onsetCt = 0; onsetCt < _onsets.size(); ++ onsetCt)
      {
        _metafile.readReal(_onsets[onsetCt]);
      }
    }

    // TODO: make Recording constructor that takes logfile and streams data into it
    // guess it could be private, producd by processor

  private:
    File _metafile;
    FileRef _logfile;
    std::string _devStr;
    Scalar _tunedHz;
    Scalar _rateHz;
    HeapVector<Scalar> _onsets;
    std::streamsize _onsetsStart;

    /*
     * hmm how to handle logfiles
     * there are a handful of different formats, for example gnuradio has a raw format and a metadata format, and rtl-sdr has a raw format
     * seems like I'll need a separate classy thing for handling this
     * pass an object into this
     *
     * there are two situations:
     * 1. we have live data, to write to file
     * 2. we have a file, to read from; it's already written
     *
     * we have class fro read from file; can get path from it
     * we have class for reading from radio; would need auxiliary data for writing
     * 
     */
    Recording(std::string fname, std::string iqfname
  };

private:
  File _modelfile;
  Scalar _periodMean;
  Scalar _periodPrecision;
  Scalar _periodVariance;
  uint64_t _periodCount;
  std::vector<FileRef> _files;
  std::streamsize _filesStart;
};
#endif 

// a useful structure for looking at this data would be downsampling it as a
// distribution.
// consider chunks of at least 2 samples and calculate the mean and std dev,
// and plot them separately.  consider also the local std dev if mean is taken as
// the whole signal mean (or eg a moving average)

// downsampling could be done by reshaping into a matrix and taking colwise/rowise
// downsample operation !

struct TuningFreq
{
  Scalar freq;
  size_t idx;
};

// - [ ] resolve crash after significant tuning, perhaps by decreasing ram use?
// - [ ] change tuning structure to store indices using above structure
// - [ ] output complete spectrum when done

int main(int argc, char const * const * argv)
{
  // New period has high variability.
  // Old period seemed very precise.
  // Perhaps a bug in my convolution.
  // I'll debug the issue in depth.
  // 1. - [X] Identify onsets used to determine period guess.
  // 2. - [X] Re-use older profiler to determine correct period guess.
  //        -> used audacity instead
  // 3. - [ ] Examine convolution of onset sums to discern why guess was wrong.
  //    ========>>>>>>>>>
  //    It looks like rotation adjustment is working now,
  //    with some fixes.
  //    The second time rotation is run, the active buffer
  //    is much more accurate than the first (much sharper),
  //    and the wrong rotation is spewed back! it's in
  //    the exact correct place but it's less blurry,
  //    and the fft correlation gives an ideal adjustment
  //    of 25943 samples =S
  //    - [ ] thinking I'll load the data into python and
  //          review what the fft product looks like
  //          -> it's notable that pastBuf is smaller
  //             than needed ! activeBuf is 42 samples larger
  //    <<<<<<<<=========
  //            First buffer:
  //              onsetCt = 18
  //              _approxWavelen = 59962.523422860715 (40.025 Hz)
  //              meta.sampleTiem = 0
  //              chunk.size() = 1024000
  //              state.active:
  //                  periodsConsidered() = 17
  //                  significance() = 0
  //              effectivePhaseOnset = 9
  //              activePhaseTimestamp = 539663
  //              activePhasePeriod = 9
  //            Then:
  //              firstPhaseTimestamp = 539663
  //              firstPhasePeriod = 9
  //              lastOnset = -4637
  //              lastPeriods = 17
  //              onsetCt = 18
  //              _onsetStart = 1024000
  //              _onsetStartPeriods = 17
  //              effectivePhaseOnset = 9
  //              activePhaseTimestamp = 1559026
  //              activePhasePeriod = 26
  //
  //              _pastBuf in '_pastBuf.dump'
  //              _activeBuf in '_activeBuf.dump'
  //              rotation = -713
  //              '_activeBuf_rotated.dump'
  //              activePhaseTimestamp = 1558313
  //              _onsetStart = 1024713
  //              _approxWavelen = 59920.588235294119
  //                // shorter than earlier, which is correct by graph
  //            Then 3:
  //              firstPhaseTimestamp = 539663
  //              firstPhasePeriod = 9
  //              lastPhaseTimestamp = 1558313
  //              lastPhasePeriod = 26
  //              lastOnset = -10322
  //              lastPeriods = 34
  //              onsetCt = 18
  //              _onsetStart = 2048000
  //              _onsetStartPeriods = 34
  //              effectivePhaseOnset = 9
  //              activePhaseTimestamp = 2576963
  //              activePhasePeriod = 43
  //
  //              rotation = -2762
  //
  //        1:
  //            activePhaseTimestamp = 1559026
  //            _onsetStart = 1024000
  //            rotation = -713
  //            activePhaseTimestamp = 1559739
  //            _approxWavelen = 600004.460588235294
  //        2:
  //            activePhaseTimestamp = 2579815
  //            _onsetStart = 2048000
  //            rotation = 25943
  // 4. - [ ] fix period integration to not include background noise

  //RtlSdrIQDump data(std::cin);
  
  constexpr Scalar MIN_TUNE_FREQ =  400000000;
  constexpr Scalar MAX_TUNE_FREQ = 1600000000;
  constexpr Scalar FREQ_GUESS = 8000;
  constexpr Scalar FREQ_GUESS_ERROR = FREQ_GUESS / 100;
  constexpr size_t SAMPLERATE = 2400000;
  constexpr size_t TUNE_WIDTH = SAMPLERATE / 2;

  SoapyLive data("", 43, SAMPLERATE);
  constexpr size_t FFT_BUFFERSIZE = 4 * 1024;// * 2 / 5 + 1;
  constexpr size_t RADIO_BUFFERSIZE = SAMPLERATE / 10;
  constexpr size_t INITIALIZATION_SECS = 1;//10;
  constexpr size_t DOWNSAMPLING = INITIALIZATION_SECS * SAMPLERATE / FFT_BUFFERSIZE;
  constexpr Scalar BUFFERS_PER_SEC = SAMPLERATE / double(FFT_BUFFERSIZE);
  //PeriodFinderFFT<Scalar> periodFinder(SAMPLERATE / (FREQ_GUESS * 2) + 2, FFT_BUFFERSIZE * DOWNSAMPLING / 5 - 2, FFT_BUFFERSIZE, DOWNSAMPLING);

  WaveformModelStatsAccumulator<StatsAccumulatorHistogram<Scalar>, STATS_STANDARD_DEVIATION> initialModel({data.epsilon()}, SAMPLERATE / Scalar(FREQ_GUESS - FREQ_GUESS_ERROR));
  //WaveformModelStatsAccumulator<StatsAccumulatorNormal<Scalar>, STATS_STANDARD_DEVIATION> initialModel({}, SAMPLERATE / Scalar(FREQ_GUESS - FREQ_GUESS_ERROR));

  PeriodProcessorDetectKnownCorrelate<Scalar, decltype(initialModel), FFTCorrelationFunctor<Scalar>> processor(initialModel, SAMPLERATE / Scalar(FREQ_GUESS), {}, Scalar(SAMPLERATE) * 60 * 60 * 24 * 365.25 * 100);
  //PeriodProcessorConvolveDownsample<Scalar> processor(SAMPLERATE / (FREQ_GUESS + FREQ_GUESS_ERROR) - 1, SAMPLERATE / (FREQ_GUESS - FREQ_GUESS_ERROR) + 1, SAMPLERATE / (FREQ_GUESS) / 6);
  DataFeeder<Scalar, decltype(data), decltype(processor)> feeder(data, processor);

  HeapVector<Complex> buffer(RADIO_BUFFERSIZE);
  RecBufMeta bufMeta;
  size_t lastPeriods = 0;

  //data.tune(500000000.);
  //data.readMany(buffer); // allow to settle
  //
  std::cout.precision(real_limits<Complex>::max_digits10);
  std::cerr.precision(real_limits<Complex>::max_digits10);

  std::vector<TuningFreq> tune_freqs{};
  for (Scalar freq = MIN_TUNE_FREQ; freq <= MAX_TUNE_FREQ; freq += TUNE_WIDTH)
  {
    TuningFreq obj;
    obj.freq = freq;
    obj.idx = tune_freqs.size();
    tune_freqs.push_back(obj);
  }
  std::vector<Scalar> tune_dBs(tune_freqs.size());

  auto beginTime = std::chrono::steady_clock::now();
  std::random_device rand_dev;
  std::mt19937 rand_gen(rand_dev());

  // TODO FIX TODO: REMOVE DEBUG LINE
//#warning REMOVE DEBUG LINE PREVENTS FULL SCAN
  //tune_freqs.resize(tune_freqs.size() / 10);
  //tune_dBs.resize(tune_freqs.size());

  std::shuffle(tune_freqs.begin(), tune_freqs.end(), rand_gen);

  auto curFreq_it = tune_freqs.begin();
  data.tune(curFreq_it->freq);


  for (data.readMany(buffer, bufMeta); buffer.size(); data.readMany(buffer, bufMeta))
  {
    // TODO: try using magnitude fo complex variance, and try twice as much data with i & q both considered real.  which of the 3 approaches has the most accurate stats?
    // could also look at all 3 metrics on the raw data: the best one is the most extreme and the most reliable
    auto preprocessed = buffer.array().real().eval();
    //auto mag = buffer.array().abs().eval();
    //auto preprocessed = buffer.array().real().abs().eval();

    feeder.add(preprocessed, bufMeta);
    if (lastPeriods != processor.periodsConsidered()) {
      using Eigen::numext::sqrt;
      lastPeriods = processor.periodsConsidered();

      auto elapsedTime = std::chrono::steady_clock::now() - beginTime;
      auto percentDone = (curFreq_it - tune_freqs.begin() + 1) / Scalar(tune_freqs.size());
      auto remainingTime = elapsedTime * (1 / percentDone - 1);
      auto hrs = std::chrono::duration_cast<std::chrono::hours>(remainingTime);
      remainingTime -= hrs;
      auto mins = std::chrono::duration_cast<std::chrono::minutes>(remainingTime);
      remainingTime -= mins;
      auto secs = std::chrono::duration_cast<std::chrono::seconds>(remainingTime);
      std::cerr.fill('0');
      std::cerr << "Tuned to: " << curFreq_it->freq << " Hz (" << percentDone*100 << "% ETA: " << std::setw(2) << hrs.count() << ":" << std::setw(2) << mins.count() << ":" << std::setw(2) << secs.count() << ")" << std::endl;
      std::cerr << "Periods: " << lastPeriods << std::endl;
      std::cerr << "Wavelength Avg: " << processor.bestPeriod() / double(SAMPLERATE) << std::endl;
      std::cerr << "Frequency Avg: " << double(SAMPLERATE) * processor.bestFrequency() << std::endl;
      //std::cerr << "Wavelength STD: " << sqrt(processor.periodVariance()) / double(SAMPLERATE) << std::endl;
      /*
         std::cerr << "Best significance so far: " << periodFinder.bestPeriod() << " (" << Scalar(SAMPLERATE) / periodFinder.bestPeriod() << " Hz " << periodFinder.bestSignificance()*100 << " %)" << std::endl;
         std::cerr << "                          " << periodFinder.bestPeriod2() << " (" << Scalar(SAMPLERATE) / periodFinder.bestPeriod2() << " Hz " << periodFinder.bestSignificance2()*100 << " %)" << std::endl;
      */
      auto mag = processor.bestMag(bufMeta.freq);
      auto err = sqrt(processor.bestMagVariance(bufMeta.freq)) * 3;
      std::cerr << "                          " << mag << " mag " << err << " err" <<std::endl;
      //tune_dBs[*curFreq_it - ]
      auto mindB = 20 * log10(mag - err);
      auto maxdB = 20 * log10(mag + err);
      auto dBctr = (maxdB + mindB) / 2;
      auto dBerr = (maxdB - mindB) / 2;
      std::cerr << "                          " << dBctr << " dB +-" << dBerr << " dB" <<std::endl;

      tune_dBs[curFreq_it->idx] = 20 * log10(mag);
      std::cout << curFreq_it->freq << " " << tune_dBs[curFreq_it->idx] << std::endl;
      ++ curFreq_it;

      /*
      // stop when stats imply small enough significance that there would be one error in a year of trials, and dB is within 0.1
      if (periodFinder.bestSignificance() <= 1.0 / (365.25 * 24 * 60 * 60 * BUFFERS_PER_SEC / lastPeriods) && dBerr <= 0.1)
      {
        break;
      }*/
      if (curFreq_it == tune_freqs.end())
      {
        std::shuffle(tune_freqs.begin(), tune_freqs.end(), rand_gen);
        curFreq_it = tune_freqs.begin();
        beginTime = std::chrono::steady_clock::now();

        /*for (size_t i = 0; i < tune_dBs.size(); ++ i)
        {
          if (i)
          {
            std::cout << ", ";
          }
          std::cout << tune_dBs[i];
        }*/
        std::cout << std::endl;
      }
      data.tune(curFreq_it->freq);
    }
  }
  /*
     std::cerr << "Best significance: " << periodFinder.bestPeriod() << " (" << Scalar(SAMPLERATE) / periodFinder.bestPeriod() << " Hz " << periodFinder.bestSignificance()*100 << " %)" << std::endl;
     std::cerr << "                   " << periodFinder.bestPeriod2() << " (" << Scalar(SAMPLERATE) / periodFinder.bestPeriod2() << " Hz " << periodFinder.bestSignificance2()*100 << " %)" << std::endl;
     std::cerr << " (Range is "
     << periodFinder.minBestPeriod() << "=" << Scalar(SAMPLERATE) / periodFinder.minBestPeriod() << " Hz  to "
     << periodFinder.maxBestPeriod() << "=" << Scalar(SAMPLERATE) / periodFinder.maxBestPeriod() << " Hz)" << std::endl;
  auto mag = periodFinder.bestMag();
  auto err = sqrt(periodFinder.bestMagVariance()) * 3;
  auto mindB = 20 * log10(mag - err);
  auto maxdB = 20 * log10(mag + err);
  auto dBctr = (maxdB + mindB) / 2;
  auto dBerr = (maxdB - mindB) / 2;
  std::cerr << "                          " << dBctr << " dB +-" << dBerr << " dB " << std::endl;
  */

  //auto significances = periodFinder.significances();
  //for (size_t i = 0; i < significances.size(); ++ i)
  //{
  //  std::cerr << significances[i] << " ";
  //}
  //std::cerr << std::endl;
  
  //PeriodFinderAccumulateAllInteger<PeriodModelAccumulatorNoiscillate<StatsAccumulatorHistogram<Scalar>>> periodFinder2(periodFinder.minBestPeriod(), periodFinder.maxBestPeriod(), StatsAccumulatorHistogram<Scalar>(data.epsilon()));

  //for (lastPeriods = 0; buffer.size(); data.readMany(buffer))
  //{
  //  auto preprocessed = buffer.array().real().eval();
  //  periodFinder2.add(preprocessed);
  //  if (lastPeriods != periodFinder.periodsRead())
  //  {
  //    
  //  }
  //}
  return 0;
}

/*
 * oh my god I've been working so hard to do this and it was suddenly so easy!
 * but when i try to inner dialogue around that, to see how it might happen in the future
 * i got some nasty feedback, like a desire to hit myself hard in my head
 */

// TODO GOAL [NEXT]:
// Seems the ideal execution path might be:
// - identify starting period
// - identify the number of periods needed to accumulate to determine if the phase has
//   shifted out
// - keep adjusting for shift, using stats to do so.  we expect a constant average period
//   but our average will never be exactly right, so there may always be a 1 sample shift
//
// but first of course we want to identify dB intensity if we have the right period, and
// of course we're flexible in that measurement up to around 1/4 period shift, so
// measuring intensity and adjusting for phase drift can go in parallel
//  -> when very quiet it may drift more than 1/4 period in the time it takes to see it,
//     so it's important to make use of subsample periods when predicting onsets

/*
 * My periodfinder that iterates through each period is kind of slow.  it has n*m
 * complexity without any vectorization ... I bet I could do this with matrix ops.
 *
 * We have m periods we want to check, and a buffer of n samples.
 * We want to accumulate the n samples in m different ways.
 *
 * So each n sample is turned into an m-length row and summed into m accumulators.
 * Sounds pretty reasonable!  Just need to write the wrapping code.
 *
 * > hey, the FFT gives us the phase !
 *   we can use this phase for every nearby period
 *   we just need to ignore the samples where the inaccuracy of the phase could
 *   affect things
 *
 * possible accumulators:
 *  (sin*sig)^2 + (cos*sig)^2 for associated period
 *    -> could just use sin or cos if we know phase
 *  sliding noiscillate stats bins
 *
 *  sincos^2 is more conducive to linear algebra. less algorithmic code needed
 *
 *  it would be really nice to generalize this parallel-accumulator-matrix approach.
 *  
 */

/*
 * Making accumulator for square wave rather than sinusoidal
 *      111111
 * 00000      00000
 *
 * 11111      11111
 *      000000
 *
 * we want results to maximally differ when square waves are out of phase vs in phase.
 *
 * summation works!
 *
 * - adjust wave to be DC centered
 * - create model wave with same magnitude, DC centered
 * - some with real wave.  In phase produces *2 magnitude DC signal, if abs is taken.
 * - Out of phase produces moire effect, likely a sinusoidal signal centered around
 *   *0.5 magnitude, with *0.5 magnitude.  Not sure if this is a sinusoid, maybe more
 *   a triangle wave.
 * 
 * - ways of speeding up processing of square wavelet interference:
 *   - could check derivative of interference wave.  DC will have minimal derivative.
 *   - Wave is a noisy sinusoid of a frequency related to the distance from the correct
 *     answer.
 *
 * I could interfere _one_ wavelet, and use the frequency of the interference signal
 * to identify the actual frequency, and adjust.
 *
 * My original FFT performs this interference.
 * Notably, though, I get a good answer after just 1 FFT.
 *
 *       | | | | | | | | | | | |
 *    -
 *    -
 *    -
 *
 * When we manually do the interference, we accumulate every sample.
 * But it takes many _periods_ each consisting of many _samples_ to identify the precise
 * answer.
 *
 * After 1 sample, everything's pretty much in phase.
 * We get interference once we get to half the minimum period length # of samples, or
 * so.  There's a predictable point at which there is a sudden change.
 *
 * This point is really what's interesting.  It occurs again and again, and we can
 * guess where it is with given accuracy.  The spots we know just give us information
 * on the peak and valley populations.
 *
 * It's notable that after some number of periods, these spots won't work anymore
 * as background population sources, because the phase drift will be too great.
 * We can predict the minimum time at which that will happen.
 *
 * We can also look at the phase change between FFT's to figure this out !!!!!
 * I think that may be all that's needed.
 *
 * The valuable piece of information is the mean of the slope of the phase.
 * We'll need to accumulate phase, so that it wraps around, which will only work
 * if the signal phase is much louder than the noise phase ie the signal itself
 * is much louder than the noise.
 *
 * One way to ensure this I think would be to sum the ffts as complex numbers
 * rather than real numbers.  This way the phases will sum, and the random changes
 * will average out.
 *
 */

/*
 * Consider an FFT of a weak signal among noise.
 * The strongest pieces of the FFT are the noise.  Random phase, random amplitude,
 * but kind of a _floor_ of minimum amplitude.
 *
 * I assume I can find the signal but looking for a strong statistical difference
 * in the peak -- a lot of averaging, basically.  This appears to work in practice,
 * looking at gqrx and osmocom_fft.  The question is how to preserve the phase
 * with this averaging.  If the noise is louder than the signal, the phase in each
 * sample is basically the phase of the noise -- random -- with only a small adjustment
 * for the phase of the signal.  Normally averaging would remove this, but because of
 * how these things are placed in complex space it's not so easy.  The noise signal,
 * because it is strongest, is centered about the origin, spinning all around, and
 * the real signal is summed onto that, away from the origin.  So the change in the
 * phase from the real signal is not directly related to the real signal phase, but
 * indirectly related via geometry.
 *
 * noisevec + sigvec
 * produces a triangle.  the angle change is related to the magnitude of the noise,
 * for each sample, by some trigonometry.  Higher magnitude noise produces smaller
 * change.  We could scale the sample up by its magnitude in the proper way
 * to approximate the correct phaes change; that would not properly approximate the
 * correct magnitude change.
 *
 * So if I'm averaging to find phase, I'd want to apply a transformation to the
 * samples first.
 *
 * atan2(noisevec + sigvec) is related how to atan2(sigvec) ?
 * atan2(cs(n_ang) * n_mag + cs(s_ang) * s_mag)
 * sin(na) * nm + sin(sa) * sm; sm << nm
 *
 * The effect of sigvec_ang on the final angle is related to the distance of sigvec
 * from the origin, which is the magnitude of noisevec.
 *
 * ISSUE! if we weigh the samples by a transformation, then the noise may not cancel
 * out when unweighed.  It's notable that the noise is theoretically roughly the
 * same magnitude each sample (maybe) ... but still in practice it doesn't work
 *
 * I'm leaning towards the approach here being to accumulate samples to decrease the
 * noise below the level of the signal magnitude prior to doing the phase thing.
 *
 * So I have the signal in question in my FFT, and I have the noise floor.
 * I know the max magnitude of the noise floor.
 *
 * I'll want to sum the fft's such that the noise contribution to the signal is minimal.
 * Then I can use adjacent fft samples to determine the phase.
 * It's the slope of normalized complex fft result.  We need each adjacent piece to have
 * the signal louder than the noise; I believe this is doable.
 *
 * It's also notable that the signal is recoverable from the FFT, and may be analyzable
 * this way ... but there's no need, once we have it precise we can look at the flat
 * peak, perhaps using existing code.
 *
 * Hmm the noise here is different from the noise I was looking at before.
 * 
 * Well I guess the assumption is that the signal's fft is something summed onto
 * the background fft; I'm not sure if that is true ...
 *
 * in time domain s1*m1 + s2*m2, c1*m1 + c2+m2
 * in freq domain ang, mag
 * hmm I feel this should be obvious
 */



// TODO NEXT FOR REAL: adjust fft accumulator to sum in a complex manner, so it can
//  track changes to the phase.
//      since we're planning on doing signals that are weaker than the noise, and phase
//      information can be lost in such a situatio, we'll want to get the fft's strong
//      enough such that the signal is louder than the noise in the accumulation
//      -> note that this noise may or may not be the same as the backgroundDist;
//         considering noise in maxIdx signal itself, to extract phase
//      we'll want to apply the same statistical approach to the phase: get it accurate
//      within some significance
//      but we'll want to pay attention to the times of the phases somehow, so that
//      we can identify phase drift
//      might need a different kind of accumulator!
//
//      boost accumulators would be good for this

// FFT PHASE APPROACH
//  needs complex sum
//  needs sample time to extract phase slope -> can perhaps be done with a second FFT over the waterfall
//    -> could also maybe multiply/divide adjacent fft values to get average change in phase

// ONGOING INTERFERENCE APPROACH
//    I think the square-wave interference could be done really generally, and keep on
//    going to narrow down precisely, without the FFT.  window of application shrinks
//    as stastics provide confidence
//
//    If we don't want to guess the magnitude, previous samples could be referenced to
//    provide the interference.

// i'm not sure what interference was, but the current plan is to combine periodfinderallinteger and the ringer idea
// -> the concept of taking 4 bins will be factored out into a replacable algorithm
// -> we only need 2 bins if we're testing the phase as well as the period
//
// so the period finder will be the interconnection between the choice of what periods to test, and how to test them
// since one of the testers tests both periods and phases, this period finder will want to generalize the concept of a range of variables to test within.
//
// : range of variables to test
// : testers themselves
// : incoming data passed to testers
// : pattern in which variables are chosen to pass to testers
// : information flow from the testers to the pattern to choose what to text next (feedback)
//
// that's pretty much it -- the information flow is what's missing.
// there are two kinds of information: confidence/significance, and nearness
// so if I test [A] and determine it is 3.9, and I'm looking for a minimum, I know very little
// if I test [A],[B],[C] and determine they are 3.2, 3.9, and 1.2, it looks like I should disregard A in my search for a minimum.
// But I'll need some degree of error to determine if it's meaningful.  If they are +- 4, I have very little information.
// I'll need a way to combine the error information from the results to determine how significant the results are.
//
// So, the proposal of adaptive subsampling is that the best answer is within a smaller range than the one passed in.
// It's based on knowledge of the shape of the underlying data.  With these square waves, we expect a lot of locality ....
//
// When we pick a zero-crossing, we parcel out the data into possibly-high and possibly-low.  To judge this data, we can take the standard deviation (if it's real) or the mean
// (if it's absolute), and compare that to the sampling dsitribution oft he best population we ahve available.  This will give us a % likelihood that we sampled the correct
// populations.
//    -> this could be improved statistically since part of our sample _is_ from the correct population, and the other part is from a known population
//    but the end result is we get a % of how likely this sample is, if we actually sampled the correct population
//
// This % will fall off smoothly from a peak of choosing the correct zero-crossing: it'll get more jittery as more noise is introduced
//
// The adaptive subsampler will get these %'s and have to decide whether to sample the whole area more to get more data and increase the confidence, or whether to accept
// some regions as not of interest or not worthwhile, and sample a likely region more densely.  (it should spew out a list of what to sample next, so the implementation can
// vectorize or parallelize).
//
// How do we combine these %'s into one combined significance?  It appears to depend on the underlying data.
// hmmm this has a 4% chance of being wrong, this has a 2% chance of being wrong; this has a 10% chance of being wrong
// if we want to hypothesize that all are right, we can combine them with probability
// 98% of being right, and 96% of being right, and 90% of being right
//
// 96% : A is right
//        98%: B is right
//         2: B is wrong
//  4  : A is wrong
//        98%: B is right
//         2: B is wrong
// if A and B are both right, then the result is 98% of 96% = 0.98 * 0.96: a lower chance that both are right
// if A and B are both wrong, we get 0.04 * 0.02: a very low chance of both being wrong
//
// | | | | | | |.| | |
// it's notable that none of them are right, because they are all approximations.
// the result of this metric will be a curve approaching the right value.
//
// how do we hypothesize that we are on one side of the goal, without regard to whether we are at it or not?
// we look at only 1 of the 2 distributions for the new approach
// if we have many samples of a region, and we want to propose that the best choice is within a particular area of the region, we want to use the portions of each sample that
// help us back up that choice.
//
// what are the chances ... etc ...
// hmmmmmmmmmmmmmmmmmmm this doesn't quite seem to work!
// what information do I really want from the ringer, for the adaptive subsampler to make the best choice?  I'm not sure that it's the standard deviation sampling thing I was
// looking for.
//
// The ringer wants to break the variable space into two regions and test the hypothesis of the result being in one, or in the other.
//
// What do we need in order to test this hypothesis, for square wave noise?
// We have k spots we want to gain metrics around.
// For each k spot, we partition the period in two, and get stats attributes for the two partitions
//      -> this can be done by getting stats for the segments, and summing them
// We can then consider that the answer is between any two spots, by considering the stats to the left and right
// I guess we'd need to start with 3 spots?
//    we'd consider only 1 stats group for each test, even though we have 3 ! but 2 of the 6 actual bins are empty
//
// stats-edges-approach: for each pair, test if the answer is within that pair by ignoring the stats within the region, and testing the stats outside it
// maximum-approach: subsample the pair with the maximal metric; disregard other samples
// interpolation: determining an accurate guess by using the information within a region.
//
// i'm used to adaptively subsampling object edges.   i usually start with screen edges and object centers: we divide by two until hitting the edge
// it's helpful to do this in this manner because objects have lots of changes near the edges, and change slowly away, and we want to resolve it all
// but here we _only_ want to resolve the edge: so it could be meaningful to make a good guess, and then consider points around that guess that have
// to do with the accuracy of the guess.  this would be most effective.
//
// the maximum approach combines the 2 sets of information for each point, and then subsampling can pick the two best spots to look down into
// the stats-edges approach doesn't use a region of information in making the decision, but it gets information directly applicable to determining
// the confidence of the decision
//
// trying to understand this unused region ... ideally i want to use as much information as possible so as to make the best decision
// i guess we could use the stats on the region to inform the decision .... if they are off-base compared to the other two groups, it's a pretty good indicator,
// but this would only happen if the real position is not near the edge of the region
//
// i'll try looking at this a different way.
// consider all the spots: could we form a decision on whether or not the change is to the right of the spot? or to the left?
// if it is in one direction, then we know the level of that spot, and the level of all the spots in the other direction.
// hum this is the same information !
//
// perhaps the trick is to consider when the change is near the edge of the region
// then considering that the change is _at_ a spot will appear valuable, and will give more confidence than considering that the change is in the region
//
// say we were to compare the proposition that the change is at the spot vs in the region
// this is really different depending on whether the change is near the spot or not !
//
// the maximum approach considers only that the change is at spots, and it just picks the ones that look the best
// the data in the region: we have the moments
// 
// the number of periods to test is in the 10000s at least (perhaps 50000s), and the data that needs to be accumulated and processed for each one is in the millions
// processing each one involves summing chunks of the millions data
// so if I were to test each one, it's O(n * m)
// I'm trying to make it be O(n log(m))
//
// but i'm not sure I really understand this .... processing one of these 40000 periods means breaking the 2m datapoints into two groups, and summing each one
// if I were to process each 40k period, I'd need to make 40k of these sums
// but each adjacent one is just slightly different from the other one
// so I could assume the first period: sum all the data into the 2nd group
// as the 1st group moves forward (oh I don't quite have it right, it's more complicated than that)
// as the 1st group moves forward, we'd add one sample to one sum, and remove it from the other
// this is O(2 n)
// we don't have to reprocess all the data: just the changes from the previous processed one
//
// this would work well for the existing approach which is unaware of periods and leaves half of the data unconsidered
//
// how could I expand it to phases? I can't really test all the phases, but I have a small range of phases that are worth considering ...
//
// as I move forward my knowledge of where the signal is likely to be decreases due to phase shift and my uncertainty 
// so as I get more certain about the nature of the signal, I can expand the set of data I use to consider it
//
// Here's a pattern which I should try as much as possible to fit into my generic model:
// - when a buffer comes in, FFT it to get an idea of the phase and period
// - test periods using an algorithm that is not O(n^2), considering the data remaining
//
// our period length is known within a specific range
// like the same is true of the phase
//
// so for a given period within the recording we can define upper and lower boundaries for each of its zero-crossings
// this gives us regions of known high and low
//
// we can take the data from the known high and low, put it in population accumulators, and discard it
// 
// the data from the known periods is different depending on its distance .....
// -> I'll need to use a complex FFT in order to get phases
// -> I'll need to adjust the FFT phase in each result in order to get the complex answers to line up (there might be other approaches, like taking the geometric mean)
//
// Okay, I'm not sure if I can extract the phase from the FFT for sure or not, so I should provide for it also coming from a 4-bin system like AllInteger, which requires no
// storage but can eventually approximate the phase within 90 deg; could be modified to have arbitrary [low] precision for some ram
// would need to optimize it to be O(2n)
//    -> is it really O(2n) to slowly shift the bins, one sample at a time?
//       I think so.  This shifting lets the bins track their std dev, variance, and mean ... the only loop is over the bins, not the data. at worst O(2n + C m)
//
//
//
// hrm ... how do we accumulate data with this approach?
// the old approach was to have 4 stats bins _for every possibly period_
// so they could all accumulate data coming in
// if we adjust the stats bins, can we still accumulate?
// when we drop the old data, we lose the ability to adjust them.
//
// so this is going to be a new approach.
// it needs to have a clue of the period (not the phase) and then it may be able to store some information that accumulates in the bins
//
// FFT guesses period range: doesn't need storage to do so
// AllIntegers2 starts working within FFT range ... we just assume our phases of 0 deg, 90 deg, 180 deg, 270 deg, and fill old data while letting new come in
// => problem, if signal is weak, eventually our phase ranges will go out of sync
//    => at this time, an approach would be to save them and start accumulating more, then look for the most similar and add
//      of course, with weak enough signal this correlation would be random and the signal would never be found
//      a better approach might be longer-term adaptive subsampling
//      or this fft-of-fft thing you have so much trouble thinking about
//
// if I adaptively subsample them, can it keep going in the face of superweak signals?
// I can at least pick a best period with the data I have.
//        if the signal were weak enough, this result would be random.


// sorting out my approach has gotten derailed
// now I'm looking for a solution that will work for arbitrarily-weak signals when the error with which we know their period is too large to hold enough of them in memory to find
// them.
// I can _find_ an arbitrarily weak signal using the sum of the FFT power, but I can't precisely know its period if it is too weak to have any discernible change on a ram buffer
// full of noise.
// A window into a larger fourier algorithm would discern its period (and phase, possibly), but the algorithms available don't tend to do that
// I could make this window with a low-pass or band-pass filter.
// hmm all i need is a band-pass filter for this approach.  let's just go for it
//
//
// okay, sinc(t) * 2 * bw 
// I guess I'll read about those.
//
//
//
// I think all the approaches can be combined.
// Regarding the problem of getting stats from single metrics, with the FFT we take stats regarding the peak metric's value compared to the background metric.
//        -> with other metrics, it's probably possible to solve for the expected smooth curve to get reasonable stats.  we need to know the shape of the curve and figure out
//           how it will relate.
// Adaptive subsampling can be applied to the FFT approach too, just throwing complex sinusoids at the data.
//        -> this will probably create a wide curve, since the DFT has unit-widthed peaks

// > What's going on? why are you so thoroughly freezing him up here?
// < Just turning him off or delaying him, just doing my job.
// > He'll go right back on after you turn him off.
// < What is he working on?
// > He's writing code to measure the effectiveness of a shielded room.  He wants to build one and make it
//   easy for many people to build or have one.  We scared him into doing it many years ago; he won't drop it.
//   He wants people to be able to control what information enters and leaves them, so that they can be
//   guaranteed both privacy and safety.  A shielded room is one step towards this: it blocks radio waves
//   and electromagnetic fields and such.
//
// < Some of us can't help you, Karl
// > I'd ask you to think kindly of me.  Would helping me threaten your position?
// < No.  It would threaten my life. <Yes.>  But I'll think kindly of you.
// > Thank you.
//
// > So you know, my life was threatened, but only my autonomy, mind, and external life were harmed, despite
//   terrifying threatening experiences.
// < Good to know.
// < I might be able to help you knowing this.
//
// > I've found that I can help myself by taking nicotine with the tasks that are hard.  Forming an addiction
//   around them makes them a  pleasant experience over time.
// > I'm also studying Nonviolent Communication, which I believe can help all of us.  But it is a struggle
//   to learn it.
//
// > these files are occasionally uploaded to github and gitlab; i do have another file i started that is not
//
//
//
// > I'm also intentionally spending time studying something called Nonviolent Communication, which helps me
//   connect with people in a constructive manner, and repeatedly accesses parts of my mind showing that I
//   am not worthy of harm.  (which I lost much access to during my experience and appreciate so much)
//
// > What do you mean?  What happened to you?
// < I lost most of my livelihood, the functioning of my memory, and had many terrifying experiences.
// < I have not been significantly physically harmed, but did go through cancer, which I don't know how
//   I acquired; could be random.
//
//
// Explanation Regarding Suffering During Psychosis
// Loss of livelihood:
//    - My experience began during some intense delusion during which I lost my connection to my friends
//      and behaved in ways that sabatoged my reputation, and my trust for my own mind and behavior.
//      I experienced a sense of being harshly and completely controlled by something invisible.
//    - My moneyflow radically changed.  I found it impossible to work on things that made me money,
//      while being pressured to engage in work-for-pay.  I was not afraid of this because I was familiar
//      with living without money, and had a family to rely on.  Their supports dropped significantly, but
//      they couldnt' handle seeing me without livelihood, and kept supporting me with housing and food.
//
//
// > I'm working on code to detect attenuation of a shielded room from a weak, oscillating noise source.
//      It basically uses an obvious approach to do so: first it profiles the noise source to get an exact period
//                                                      then i plan to detect it by integrating using that
// < thank you for letting me know
// < just delaying him to do my job
// > your job involves delaying him?
// < okay, it's not my job
// > why are you doing this?
// < i don't know, it's what his brain is doing for me
// > could you fix his brain to not stop mid-way like this? you seem to have access to it in some way
// < [do i have to?] [no] [why not?] [you won't let me]
// they're stopping to let you
// ummm they're not actually heping me at this time
// yeah they had me write that
// ok i'll work with this for a it
// < he just keeps on trying to work
// > yes that is what he does
// < why are you letting him work?
// > why don't you want him to? i want him to i just seem to not be able to
// oh, the work i'm doing is important for me.  i want to build a shielded room, so i'm writing code to measure it
// i plan to associate the work with taking nicotine so i never stop
// yeah the hop is the nicotine decreases these interruption periods
// it workd really well for toothbrushing but takes time to develop
// nothing else worked for toothbrushing for a couple years (roughly), so it's really inspiring
//
// the phrase "I can't help you" is one that's been in my mind for a long time, popping up
// in response to it, we'd consider OFNR.  "Can't" sounds like a strong need is going on, likely something
//   outside the person's control affecting them.  Often associated with not losing a job.
//   We left feeling out of that.  Also left out secrecy / silence.
//   There are a lot of possible explanations.  Perhaps what's left out here is empathy more thoroughly.
// I'm learning NVC, partly as an alternative approach to this situation.
// 
