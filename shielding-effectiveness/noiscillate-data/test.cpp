#include <Eigen/Core>

#include <iostream>
#include <vector>

using Scalar = double;
using Complex = std::complex<Scalar>;

template <typename Type>
using HeapVector = Eigen::Matrix<Type, Eigen::Dynamic, 1>;
template <typename Type, size_t size>
using StackVector = Eigen::Matrix<Type, Eigen::Dynamic, 1, 0, size>;

class RtlSdrIQDump
{
public:
  RtlSdrIQDump(std::istream & stream)
  : stream(stream)
  { }
  
  Complex readOne()
  {
    std::complex<uint8_t> val;
    stream.read(reinterpret_cast<char*>(&val), 2);
    return (Complex(val.real(), val.imag()) - Complex(127.5, 127.5)) / 127.5;
  }

  template <typename Derived>
  void readMany(Eigen::MatrixBase<Derived> const & _vec)
  {
    std::vector<uint8_t> val(2 * _vec.size());
    stream.read(reinterpret_cast<char*>(&val[0]), val.size());
    Eigen::Map<Eigen::Array<uint8_t, 2, Eigen::Dynamic>> eigenRawData(&val[0], 2, stream.gcount() / 2);
    auto eigenScalarData = eigenRawData.cast<Scalar>() / 127.5 - 1.0;
    Eigen::MatrixBase<Derived> & vec = const_cast<Eigen::MatrixBase<Derived> &>(_vec);
    vec.derived().resize(eigenScalarData.cols());
    vec.real() = eigenScalarData.row(0);
    vec.imag() = eigenScalarData.row(1);
  }

private:
  std::istream & stream;
};

// maintains a running expected value and error given incoming samples
// assumes gaussian input
template <typename T>
class StatsAccumulatorGaussianPct
{
public:
  using Scalar = T;

  StatsAccumulatorGaussianPct(T pct)
  : cutoff(erfinv(pct)),
    mean(0),
    dataSum(0),
    dataSquaredSum(0),
    count(0)
  { }

  StatsAccumulatorGaussianPct(StatsAccumulatorGaussianPct const & other) = default;

  template <typename Derived>
  void add(Eigen::DenseBase<Derived> const & chunk)
  {
    count += chunk.size();
    // note: sum() and squaredNorm() do the right component-wise thing for vectors and for
    //       matrices: i.e. treating every entry as simple an item to add to the stats
    dataSum += chunk.sum();
    dataSquaredSum += chunk.derived().matrix().squaredNorm();
    /*
     * variance = sum(delta * delta) / (n - 1)
     * delta[i] = data[i] - mean
     *
     * sum((data[i] - mean) * (data[i] - mean)) / (n - 1)
     * sum(data[i] * data[i] - 2 * data[i] * mean + mean * mean) / (n - 1)
     * (sum(data[i] * data[i]) - 2 * mean * sum(data[i]) + mean * mean * n) / (n - 1)
     *
     * mean = sum(data[i]) / n
     *
     * variance = (dataSquaredSum - dataSum * dataSum / n) / (n - 1)
     */

    mean = dataSum / count;

    // actual variance calculation moved to error() call
  }

  T const & expectedValue() const
  {
    return mean;
  }

  T error() const
  {
    using Eigen::numext::sqrt;

    T variance = (dataSquaredSum - dataSum * mean) / (count - 1);

    // the variance of an averaged set of samples is the variance of one
    // sample divided by the number of samples taken
    return sqrt(variance / count) * cutoff;
  }

  // takes a bunch of statsaccumulators, subtracts their differing means,
  // and determines their error as if treated as one dataset
  T groupError(std::vector<StatsAccumulatorGaussianPct<T>> const & group) const
  {
    using Eigen::numext::sqrt;
    /*
     * delta[i][j] = data[i][j] - mean[i]
     * variance = sum(delta * delta) / sum(n - 1)
     *
     * sum(data[i][j] * data[i][j] - 2 * data[i][j] * mean[i] + mean[i] * mean[i]) / (sum(n[i]) - 1)
     * (sum(data[i][j] * data[i][j]) - 2 * sum(data[i][j] * mean[i]) + sum(mean[i] * mean[i])) / (sum(n[i]) - 1)
     * (sum(dataSquaredSum) - 2 * sum(data[i][j] * mean[i]) + sum(mean[i] * mean[i])) / (sum(n[i]) - 1)
     *
     * mean[i] = sum_j(data[i][j]) / n[i]
     *
     * (sum_i(dataSquaredSum) - 2 * sum_i(sum_j(data[i][j]) * mean[i]) + sum_i(mean[i] * mean[i])) / (sum_i(n[i]) - 1)
     * (sum_i(dataSquaredSum) - 2 * sum_i(sum_j(data[i][j]) * sum_j(data[i][j])/n[i]) + sum_i(sum_j(data[i][j]) / n[i] * sum_j(data[i][j]) / n[i])) / (sum_i(n[i]) - 1)
     *
     * (sum_i(dataSquaredSum[i]) - 2 * sum_i(dataSum[i] * dataSum[i] / n[i]) + sum_i(dataSum[i] * dataSum[i] / n[i] / n[i])) / (sum_i(n[i]) - 1)
     *
     * (sum_i(dataSquaredSum[i]) - 2 * sum_i(dataSum[i] * mean[i]) + sum_i(mean[i] * mean[i])) / (sum_i(n[i]) - 1)
     * (sum_i(dataSquaredSum[i] - 2 * dataSum[i] * mean[i] + mean[i] * mean[i]) / (sum_i(n[i]) - 1)
     * (sum_i(dataSquaredSum[i] + mean[i] * (mean[i] - 2 * dataSum[i])) / (sum_i(n[i]) - 1)
     */

    T numerator = 0;
    T denominator = 0;
    for (auto const & accum : group)
    {
      numerator += accum.dataSquaredSum + accum.mean * (accum.mean - 2 * accum.dataSum);
      denominator += accum.count;
    }
    T variance = numerator / (denominator - 1);
    return sqrt(variance / denominator) * cutoff;
  }

private:
  T cutoff;
  T mean;
  T dataSum;
  T dataSquaredSum;
  T count;

  // from https://stackoverflow.com/a/40260471
  static constexpr T erfinv(T x)
  {
    using Eigen::numext::sqrt;
    using Eigen::numext::log;

    x = (1 - x)*(1 + x);        // x = 1 - x*x;
    T lnx = log(x);
    
    T tt1 = 2/(EIGEN_PI*0.147) + 0.5f * lnx;
    T tt2 = 1/(0.147) * lnx;
    
    return sqrt(-tt1 + sqrt(tt1*tt1 - tt2));
  }
};

// TODO: our signal _IS NOISE_ which complicates things!
// to resolve this you just have to start writing about what it implies.
// note 1: our signal is _SIGNAL NOISE_ summed onto _BACKGROUND NOISE_: so we expect it, when active, to be simply a differing gaussian dist. 
//            since we've taken the absolute magnitude of the signal, the mean will rise, but otherwise it would still be zero and only the variance would have risen.

// TODO: run it!  the pieces are almost all there!

// uses 4 stats bins in an attempt to quickly identify noise that's toggled in a square wave with 50% duty cycle
template <typename StatsAccumulator, size_t NUM_BINS = 4>
class PeriodModelAccumulatorNoiscillate
{
public:
  using Detail = StatsAccumulator;

  PeriodModelAccumulatorNoiscillate(size_t period, Detail const & accum = {})
  : offset(0),
    bins{accum, accum, accum, accum},
    period(period),
    periodStride(period)
  { }

  template <typename Derived>
  void add(Eigen::PlainObjectBase<Derived> const & chunk)
  {
    // okay .... we want to store lengths of chunk into the bins
    //
    // note that in my first example, the period length is longer than the chunk
    // length.  but should code for the opposite too.

    // which bin does each sample go in?
    // period is length of 4 bins
    // length of 1 bin is period / 4, max remainder is 3
    //
    // (offset % period) / (period / 4) is bin number
    // (offset % period) * 4 / period is bin number
    //    
    // (offset % period) % (period / 4) would be bin offset except it doesn't deal
    // with rounding error
    //
    // TODO NEXT
    //
    // Okay, so, given an offset into a period, the bin is pretty precisely
    // binnum = offset * numbins / periodlength
    //
    // given a bin, the start is
    // offset = ceil(binnum * periodlength / numbins)
    //
    // an integer ceil can be taken via
    // inter = binnum * periodlength
    // offset = inter / numbins + (inter % numbins ? 1 : 0)

    // our offset will be likely partway through a bin (but not always)
    // 1. fill first partial bin to the end, if needed
    size_t bin = offset * NUM_BINS / period;
    size_t tail = binStart(bin + 1);
    size_t chunkPos = tail - offset;
    if (chunkPos >= chunk.size())
    {
      bins[bin].add(chunk);
      return;
    }
    bins[bin].add(chunk.head(chunkPos));
    offset = tail % period;
    bin = (bin + 1) % NUM_BINS;

    // 2. use maps and strides to fill a mess of bins
    
    size_t chunkTail = chunkPos + period;
    if (chunkTail > chunk.size())
    {
      // TODO: make sure we have the right periodCount for every map made
      //       depending on how that is calculated, this line could be changed, as atm
      //       it doesn't pick a tail that's aligned with period boundaries
      chunkTail = chunk.size();
    }
    // also need to find total length

    // okay, we start at chunkPos in chunk, and I guess the appropriate way is
    // to use the outer stride to look at the whole chunk 4 different times.
    // need to keep the total size within proper limits

    // NOTE it is fine to pass a rectangular matrix to statsAccumulator.add(); it will Do The Right Thing

    // TODO
    //  - [ ] set periodCount
    //  - [ ] loop over all bins
    //  - [ ] add map content to bin
    Eigen::Map<Derived const, 0, decltype(periodStride>) map(chunk.derived().data() + chunkPos, period, periodCount, periodStride);

    // 3. fill last partial bin from start, if needed
  }

private:
  size_t offset;
  StatsAccumulator bins[NUM_BINS];
  size_t period;
  Eigen::OuterStride<Eigen::Dynamic> periodStride;

  size_t binStart(size_t bin)
  {
    size_t inter = bin * period;
    return inter / NUM_BINS + (inter % numbins ? 1 : 0);
  }
};

// given some incoming data, accumulates data % a given period, passing it
// through a provided metric
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
      auto v = accum.expectedValue();
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
      values[i] = periodWave[i].expectedValue();
      errors[i] = periodWave[i].error();
    }
  }

private:
  size_t offset;
  std::vector<StatsAccumulator> periodWave;
  Eigen::InnerStride<Eigen::Dynamic> periodStride;
};



// PeriodFinder is a class that takes an incrementally-provided set of samples, and a range of possible periods, and identifies the strongest period in the data with that range
template <typename PeriodModelAccumulator>
class PeriodFinder
{
  // accepts buffers of samples
  // finds a periodic signal in them
public:
  using Scalar = typename PeriodModelAccumulator::Scalar;

  PeriodFinder(size_t minPeriod, size_t maxPeriod, typename PeriodModelAccumulator::Detail const & config)
  : minPeriod(minPeriod)
  {
    periodMetrics.reserve(maxPeriod - minPeriod);
    for (size_t i = minPeriod; i <= maxPeriod; ++ i)
    {
      std::cout << "Emplacing test period " << i << " of " << maxPeriod << " ..." << std::endl;
      periodMetrics.emplace_back(i, config);
    }
    std::cout << "Constructed."
  }

  // TODO: <==========================
  // this finds the best error
  // but with some changes, it could output a statistical metric regarding
  // the likelihood that the measured signal is not noise.
  //
  // resolving that would also resolve the need to pass a weird config object to the constructor
  template <typename Derived>
  void add(Eigen::PlainObjectBase<Derived> const & chunk)
  {
    Scalar bestRatio = 0;

    for (size_t pidx = 0; pidx < periodMetrics.size(); ++ pidx)
    {
      auto & pa = periodMetrics[pidx];
      pa.add(chunk);
      Scalar error = pa.overallError();
      Scalar range = pa.overallRange();
      Scalar ratio = range / error;
      if (ratio > bestRatio)
      {
        bestRatio = ratio;
        _bestPeriod = pidx + minPeriod;
        _bestError = error;
        _bestRange = range;
      }
    }
  }

  Scalar bestPeriod() { return _bestPeriod; }
  Scalar bestError() { return _bestError; }
  Scalar bestRange() { return _bestRange; }

private:
  size_t minPeriod;
  std::vector<PeriodModelAccumulator> periodMetrics;

  size_t _bestPeriod;
  Scalar _bestError;
  Scalar _bestRange;
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
//      - [ ] make 4 stats bins
//        - [ ] make a todo to make another class that uses a binary tree for adaptive subsampling
//      - [ ] fill stats bins with incoming data and track std dev
//            -> I'll probably drop the q data and use real-valued data for now
//            - [ ] make a todo to try using magnitude of complex variance, and to try twice as much data with i & q both considered real.  we'll want to see if the approaches using more data have accurate stats, I guess
//                -> this could also be measured by looking at all 3 metrics on the raw data.  we'll want the one that is most extreme and most reliable, I guess.
//      - [ ] maxima is taken with std dev of results (the 2 opposite bins that are most extreme are picked)
//            -> this gives range
//      - [ ] sampling distribution of the standard error can be used to determine standard error output, using the maxima
// ==============================================================================================================================================================================

int main()
{
  RtlSdrIQDump data(std::cin);

  PeriodFinder<PeriodModelAccumulatorPreciseSignal<StatsAccumulatorGaussianPct<Scalar>>> periodFinder(2048000 * 40 - 1, 2048000 * 40 + 1, StatsAccumulatorGaussianPct<Scalar>(95.0));

  HeapVector<Complex> buffer(2048000 * 2);
  while (buffer.size())
  {
    data.readMany(buffer);
    auto mag = buffer.array().abs().eval();

    periodFinder.add(mag);
    std::cout << "Best period so far: " << periodFinder.bestPeriod() << " Hz (" << periodFinder.bestRange() << " +-" << periodFinder.bestError() << ")" << std::endl;
  }
  return 0;
}
