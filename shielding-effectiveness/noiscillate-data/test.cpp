#include <iostream>

#include "types.hpp"

#include "stats.hpp"

#include "rtlsdriqdump.hpp"

// options:
// - keep using eigen buffer as a stretchy buffer and accept memory reallocations twice per process
// - use one size of eigen buffer, and track the size of it in the code
// - use a map to adjust the perceived size of eigen buffer
// - use a vector instead of eigen buffer, and map the vector when passing it in
// - eigen has a wrapper for std::vector !
                   

template <typename Scalar, class DataProcessor>
class DataFeeder
{
public:
  DataProcessor(DataProcessor & processor, std::function<bool()> checkDataInterruptedHook = [](){return false;})
  : processor(processor),
    checkInterruptedHook(hook)
  { }

  template <typename Derived>
  void add(Eigen::PlainObjectBase<Derived> const & chunk)
  {
    using min = Eigen::numext::mini;

    size_t minNeeded = processor.bufferSizeMin();
    size_t multipleNeeded = processor.bufferSizeMultiple();
    size_t maxNeeded = processor.bufferSizeMax();
    
    // exhaust buffer
    if (buffer.size() + chunk.size() < minNeeded)
    {
      buffer.conservativeResize(buffer.size() + chunk.size());
      buffer.segment(buffer.size(), chunk.size()) = chunk;
      return
    }
    auto chunkOffset = min(buffer.size() + chunk.size(), maxNeeded) % multipleNeeded - buffer.size();
    buffer.conservativeResize(buffer.size() + chunkOffset);
    buffer.segment(buffer.size(), chunkOffset) = chunk.head(chunkOffset);
    processor.process(*this, buffer);

    // pass segments of chunk
    size_t nextSize;
    while (!checkInterruped())
    {
      minNeeded = processor.bufferSizeMin();
      multipleNeeded = processor.bufferSizeMultiple();
      maxNeeded = processor.bufferSizeMax();
      nextSize = min(chunk.size() - chunkOffset, maxNeeded) % multipleNeeded;
      if (nextSize < minNeeded)
      {
        break;
      }

      process.process(*this, chunk.segment(chunkOffset, nextSize));
      chunkOffset += nextSize;
    }

    buffer = chunk.tail(nextSize);
  }

  bool checkInterrupted() { return checkInterruptedHook(); }

private:
  DataProcessor & processor;
  std::vector<Scalar> buffer;
  HeapVector<Scalar> buffer;
  std::function<bool()> checkInterruptedHook;
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
      std::cout << offset << " / " << period << std::endl;
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
    //std::cout << period << " min " << minBin << " sigma^2 " << minVariance << std::endl;
    //std::cout << period << " max " << maxBin << " sigma^2 " << maxVariance << std::endl;
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
    std::cout << "Constructed." << std::endl;
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
      std::cout << "period #" << pidx << " (" << pidx + minPeriod << "): " << significance*100 << std::endl;
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

  // TODO: this needlessly performs a wide FFT and only looks at a fraction of it.
  // instead one could apply a properly-widthed band-pass filter to the data and sparsely
  // sample it, to tighten up the FFT (perhaps freq shift it and apply low-pass filter).
  // i guess that would speed things up by roughly bufferSize / (maxPeriod-minPeriod)
  // (if correct that would be a lot !!)

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
      auto magData = fftCurrentResult.array().segment(_minPeriodsPerBuffer, foregroundDists.size()).abs().eval();
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

class PeriodFinderRunningAdjustment : public PeriodFinderBase
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

// a useful structure for looking at this data would be downsampling it as a
// distribution.
// consider chunks of at least 2 samples and calculate the mean and std dev,
// and plot them separately.  consider also the local std dev if mean is taken as
// the whole signal mean (or eg a moving average)

// downsampling could be done by reshaping into a matrix and taking colwise/rowise
// downsample operation !

int main(int argc, char const * const * argv)
{
  RtlSdrIQDump data(std::cin);

  constexpr size_t SAMPLERATE = 2048000;
  constexpr Scalar FREQ_GUESS = 10;
  constexpr Scalar FREQ_GUESS_ERROR = 5;
  constexpr size_t FFT_BUFFERSIZE = 4 * 1024;// * 2 / 5 + 1;
  constexpr size_t RADIO_BUFFERSIZE = 2048000 / 2;
  constexpr size_t INITIALIZATION_SECS = 10;
  constexpr size_t DOWNSAMPLING = INITIALIZATION_SECS * SAMPLERATE / FFT_BUFFERSIZE;
  constexpr double BUFFERS_PER_SEC = SAMPLERATE / double(FFT_BUFFERSIZE);
  PeriodFinderFFT<Scalar> periodFinder(SAMPLERATE / (FREQ_GUESS * 2) + 2, FFT_BUFFERSIZE * DOWNSAMPLING / 5 - 2, FFT_BUFFERSIZE, DOWNSAMPLING);

  HeapVector<Complex> buffer(RADIO_BUFFERSIZE);
  size_t lastPeriods = 0;
  for (data.readMany(buffer); buffer.size(); data.readMany(buffer))
  {
    // TODO: try using magnitude fo complex variance, and try twice as much data with i & q both considered real.  which of the 3 approaches has the most accurate stats?
    // could also look at all 3 metrics on the raw data: the best one is the most extreme and the most reliable
    //auto preprocessed = buffer.array().real().eval();
    //auto mag = buffer.array().abs().eval();
    auto preprocessed = buffer.array().real().abs().eval();

    periodFinder.add(preprocessed);
    if (lastPeriods != periodFinder.periodsRead()) {
      lastPeriods = periodFinder.periodsRead();

      // stop when stats imply small enough significance that there would be one error in a year of trials
      if (periodFinder.bestSignificance2() <= 1.0 / (365.25 * 24 * 60 * 60 * BUFFERS_PER_SEC / lastPeriods))
      {
        break;
      }

      std::cout << "Best significance so far: " << periodFinder.bestPeriod() << " (" << Scalar(SAMPLERATE) / periodFinder.bestPeriod() << " Hz " << periodFinder.bestSignificance()*100 << " %)" << std::endl;
      std::cout << "                          " << periodFinder.bestPeriod2() << " (" << Scalar(SAMPLERATE) / periodFinder.bestPeriod2() << " Hz " << periodFinder.bestSignificance2()*100 << " %)" << std::endl;
    }
  }
  std::cout << "Best significance: " << periodFinder.bestPeriod() << " (" << Scalar(SAMPLERATE) / periodFinder.bestPeriod() << " Hz " << periodFinder.bestSignificance()*100 << " %)" << std::endl;
  std::cout << "                   " << periodFinder.bestPeriod2() << " (" << Scalar(SAMPLERATE) / periodFinder.bestPeriod2() << " Hz " << periodFinder.bestSignificance2()*100 << " %)" << std::endl;
  std::cout << " (Range is "
     << periodFinder.minBestPeriod() << "=" << Scalar(SAMPLERATE) / periodFinder.minBestPeriod() << " Hz  to "
     << periodFinder.maxBestPeriod() << "=" << Scalar(SAMPLERATE) / periodFinder.maxBestPeriod() << " Hz)" << std::endl;

  //auto significances = periodFinder.significances();
  //for (size_t i = 0; i < significances.size(); ++ i)
  //{
  //  std::cout << significances[i] << " ";
  //}
  //std::cout << std::endl;
  
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





// When converting PeriodFinderBase to DataProcessor with a passed template algorithm,
// I found a desire to make the API nicer.  I wanted to make an interrupt() function
// instead of registerCheckInterrupted().  But since we're still single-threaded,
// there is nobody to call the interrupt() function.  Instead the DataProcessor will
// need to call out.
// A lightweight signals/slots implementation might fix this: we could register
// a signal and the user could provide a handler to determine things.
//
// I'd like it to work for multithreaded things too.  The nice thing about the hook
// is that the hook could call a signal, or it could check a condition variable
// set concurrently ...
//
// it's also notable that setting a function pointer has comparable overhead to
// using a virtual function.  It's actually easier to use a std::function from
// other code, though.
//
// Really I have streaming going on here, and I'll be plugging the class into
// a stream.  The process may be interrupted either with new data, or termination
// due to the result being satistfactory or a time limit exhausted, I suppose.
//
// why did I do this conversion?  it allows the algorithm functions to be
// inlined and called with fastcall/etc .. it allows the 
// it separates the idea of how much data is coming in. it lets the central
// class define an api function that calls the algorithm's function, without using
// virtual functions
//
// i figured this was the best choice, but got stuck
