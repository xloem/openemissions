#include "types.hpp"

#include <unsupported/Eigen/SpecialFunctions>

template <typename T> static constexpr T erfinv(T x);

using StatsOptions = int;
enum StatsOptionsBit
{
  STATS_SAMPLE =        1 << 0,
  STATS_INFINITE =      1 << 1,
  STATS_ASSUME_NORMAL = 1 << 2
};

enum StatsStatistic
{
  STATS_MEAN,
  STATS_VARIANCE,
  STATS_STANDARD_DEVIATION
};

using Eigen::numext::sqrt;

// I know the CRTP is wordy for the derived class, but it lets us have variadic
// inline functions, and these functions will be called in inner loops.  It's
// worth it.
template <typename T, typename Derived, StatsOptions OPTIONS>
class StatsDistribution
{
public:
  using Scalar = T;

  Derived & derived() { return *static_cast<Derived *>(this); }
  Derived const & derived() const { return *static_cast<Derived const *>(this); }

  Scalar mean() { return derived().mean(); }
  Scalar mean() const { return derived().mean(); }

  Scalar variance() { return derived().variance(); }
  Scalar variance() const { return derived().variance(); }

  Scalar standardDeviation() { return sqrt(variance()); }
  Scalar standardDeviation() const { return sqrt(variance()); }

  template <StatsStatistic STATISTIC>
  Scalar get() const
  {
    switch(STATISTIC)
    {
    case STATS_MEAN: return mean();
    case STATS_VARIANCE: return variance();
    case STATS_STANDARD_DEVIATION: return standardDeviation();
    default: throw std::logic_error("invalid statistic");
    }
  }

  static constexpr StatsOptions options() { return OPTIONS; }
  static constexpr bool isSample() { return OPTIONS & STATS_SAMPLE; }
  static constexpr bool isPopulation() { return !isSample(); }
  static constexpr bool isInfinite() { return OPTIONS & STATS_INFINITE; }
  static constexpr bool isFinite() { return !isInfinite(); }
  static constexpr bool assumeNormal() { return OPTIONS & STATS_ASSUME_NORMAL; }
};

template <typename T, typename Derived, StatsOptions OPTIONS = STATS_INFINITE>
class StatsDistributionByMeasures : public StatsDistribution<T, Derived, OPTIONS>
{
public:
  using Scalar = T;

  Scalar const & mean() const { return _mean; }
  Scalar const & variance() const { return _variance; }

protected:
  StatsDistributionByMeasures(T mean, T variance)
  : _mean(mean),
    _variance(variance)
  { }

  StatsDistributionByMeasures(StatsDistributionByMeasures const &) = default;

  void setMeasures(Scalar mean, Scalar variance)
  {
    setMean(mean);
    setVariance(variance);
  }

  void setMean(Scalar mean)
  {
    _mean = mean;
  }

  void setVariance(Scalar variance)
  {
    _variance = variance;
  }

  Scalar deviationSignificance(Scalar value)
  {
    using Eigen::numext::abs;
    using Eigen::numext::erfc;

    if (this->assumeNormal())
    {
      return erfc(abs(value - mean()) / this->standardDeviation());
    }
    else
    {
      throw std::logic_error("unimplemented");
    }
  }

private:
  Scalar _mean;
  Scalar _variance;
};

template <typename T, typename Derived, StatsOptions OPTIONS>
class StatsDistributionBySums : public StatsDistributionByMeasures<T, Derived, OPTIONS>
{
public:
  using Scalar = T;

  size_t size() const { return _count; }

  Scalar sum1() const { return _sum1; }
  Scalar sum2() const { return _sum2; }
  Scalar sum4() const
  {
    if (this->assumeNormal()) throw std::logic_error("not tracking 4th moment");
    return _sum4;
  }

  Scalar moment1() const { return sum1() / size(); }
  Scalar moment2() const { return sum2() / size(); }
  Scalar moment4() const { return sum4() / size(); }

  // For data which has had its absoute value taken before the statistics were taken,
  // if that data had a zero mean,
  // this alters the object such that the mean and variance will represent those of the data prior
  // to the absolute value having been taken.
  void recalcForUnderlyingDataWithZeroMeanOfAbsoluteValue()
  {
    this->setMean(0);
    this->setVariance(this->sum2() / (size() - (this->isSample() ? 1 : 0)));
  }

protected:
  StatsDistributionBySums()
  : StatsDistributionByMeasures<T, Derived, OPTIONS>(0, 0),
    _count(0),
    _sum1(0),
    _sum2(0),
    _sum4(0)
  { }

  StatsDistributionBySums(size_t count, Scalar sum1/*sum[x]*/, Scalar sum2/*sum[x^2]*/, Scalar sum4 = 0)
  : _count(count),
    _sum1(sum1),
    _sum2(sum2),
    _sum4(sum4)
  { }

  StatsDistributionBySums(StatsDistributionBySums const &) = default;

  void addSums(size_t count, Scalar sum1/*sum[x]*/, Scalar sum2/*sum[x^2]*/, Scalar sum4/*sum[x^4]*/ = 0)
  {
    setSums(_count + count, _sum1 + sum1, _sum2 + sum2, _sum4 + sum4);
  }

  void setSums(size_t count, Scalar sum1/*sum[x]*/, Scalar sum2/*sum[x^2]*/, Scalar sum4/*sum[x^4]*/ = 0)
  {
    _count = count;
    _sum1 = sum1;
    _sum2 = sum2;
    if (!this->assumeNormal()) _sum4 = sum4;

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

    this->setMean(this->sum1() / size());
    this->setVariance((this->sum2() - this->sum1() * this->mean()) / (size() - (this->isSample() ? 1 : 0)));
  }

private:
  size_t _count;
  Scalar _sum1; // sum[x]
  Scalar _sum2; // sum[x^2]
  Scalar _sum4; // sum[x^4], only used if not assuming normal?
};

//template <typename T, StatsStatistic STATISTIC, StatsOptions OPTIONS = STATS_INFINITE | STATS_ASSUME_NORMAL>
template <typename T, StatsStatistic STATISTIC, StatsOptions OPTIONS = STATS_INFINITE | STATS_ASSUME_NORMAL>
class StatsDistributionSampling : public StatsDistributionByMeasures<T, StatsDistributionSampling<T, STATISTIC, OPTIONS | STATS_INFINITE>, OPTIONS>
{
public:
  using Scalar = T;

  static_assert(!(OPTIONS & STATS_SAMPLE), "A sampling distribution is not a sample.");

  template <typename Derived, StatsOptions OPTIONS2>
  StatsDistributionSampling(StatsDistribution<T,Derived,OPTIONS2> const & sampled, size_t sampleSize)
  : StatsDistributionByMeasures<T, StatsDistributionSampling<T, STATISTIC, OPTIONS | STATS_INFINITE>, OPTIONS | STATS_INFINITE>(
      getMean(sampled, sampleSize),
      getVariance(sampled, sampleSize)
    )
  {
    if (!this->assumeNormal())
    {
      throw std::logic_error("unimplemented");
    }
    if (sampled.isSample())
    {
      throw std::logic_error("sampled distribution is not a population");
    }
  }

  Scalar standardError() { return this->standardDeviation(); }
  Scalar standardError() const { return this->standardDeviation(); }

  template <typename Derived, StatsOptions OPTIONS2>
  Scalar deviationSignificance(StatsDistribution<Scalar, Derived, OPTIONS2> const & outliers)
  {
    return StatsDistributionByMeasures<T, StatsDistributionSampling<T, STATISTIC, OPTIONS | STATS_INFINITE>, OPTIONS>::deviationSignificance(outliers.template get<STATISTIC>());
  }

private:
  template <typename Derived, StatsOptions OPTIONS2>
  static Scalar getMean(StatsDistribution<Scalar,Derived,OPTIONS2> const & sampled, size_t sampleSize)
  {
    switch (STATISTIC)
    {
    case STATS_STANDARD_DEVIATION:
    case STATS_VARIANCE:
      // TODO: inaccurate?
    case STATS_MEAN:
      return sampled.template get<STATISTIC>();
    default:
      throw std::logic_error("unimplemented");
    }
  }

  template <typename Distribution>
  [[deprecated("still needlessly assuming non-normal distributions are normal")]]
  static Scalar getVariance(Distribution const & sampled, size_t sampleSize)
  {
    switch (STATISTIC)
    {
    case STATS_MEAN:
      return sampled.variance() / sampleSize;
    case STATS_VARIANCE:
      if (sampled.assumeNormal())
      {
        // distribution is chi-square with n-1 degrees of freedom
        // chi^2 = (n-1)s^2 / sigma^2 where sigma^2 is pop variance, and s^2 is sample variance
        // variance is 2 * degrees of freeedom = 2 * (n - 1)
        //
        // s^2 = chi^2 * sigma^2 / (n-1)
        // variance of s^2 = 2 * (n - 1) * sigma^2 / (n - 1)
        //                 = 2 * sigma^2
        //
        // TODO: use an actual chi-squared distribution
        return 2 * sampled.variance();
      }
      else
      {
        throw std::logic_error("unimplemented");
      }
    case STATS_STANDARD_DEVIATION:
      if (sampled.assumeNormal())
      {
        return sampled.variance() / (2 * sampleSize);
      }
      else
      {
        // reduce:
        // sqrt of:
        //  (mu4 - mu2 ^2) / (4 N mu2)
        //
        //  (sum4/N - sum2*sum2/N/N) / (4 N sum2/N)
        //  (sum4/N - sum2*sum2/N/N) / (4 sum2)
        //  (sum4 - sum2*sum2/N) / (4 N sum2)
        //  (sum4/sum2 - sum2/N) / (4 N)
        //  (sum4/(sum2 * N) - sum2) / 4
        //auto derived = sampled.derived();
        //return ((derived.sum4() / derived.sum2() * sampleSize) - derived.sum2()) / 4;
        //throw std::logic_error("must assume the population is normal for now; formula from schaum's outline doesn't work");
        // TODO: REMOVE THIS LINE AND THE DEPRECATION ATTRIBUTE
        return sampled.variance() / (2 * sampleSize);
      }
    default:
      throw std::logic_error("unimplemented");
    }
  }
};

#include <memory>
#define EIGEN_FFTW_DEFAULT
#include <unsupported/Eigen/FFT>

template <typename T, StatsOptions OPTIONS = STATS_SAMPLE>
class StatsAccumulatorHistogram : public StatsDistributionBySums<T, StatsAccumulatorHistogram<T, OPTIONS>, OPTIONS>
{
public:
  using Scalar = T;

  static_assert(!(OPTIONS & STATS_ASSUME_NORMAL), "Use StatsAccumulatorNormal to assume normal.");
  //static_assert(!(OPTIONS & STATS_INFINITE), "StatsAccumulators cannot be infinite.");

  StatsAccumulatorHistogram(Scalar binWidth, Scalar initial = 0)
  : binWidth(binWidth),
    binDensity(1.0 / binWidth),
    binStart(Eigen::numext::floor(initial * binDensity))
  { }

  StatsAccumulatorHistogram(StatsAccumulatorHistogram const &) = default;

  [[deprecated("learn statistics and improve library")]]
  StatsAccumulatorHistogram<T, (OPTIONS & ~(STATS_SAMPLE)) | STATS_INFINITE> const &
  fakeInfinitePopulation() const
  {
    return reinterpret_cast<decltype(fakeInfinitePopulation())>(*this);
  }

  void clear()
  {
    this->setSums(0, 0, 0);
    for (size_t i = 0; i < bins.size(); ++ i)
      bins[i] = 0;
  }

  HeapVector<T> const & histogram(T & start, T & width)
  {
    start = binStart * binWidth;
    width = binWidth;
    return bins;
  }

  template <typename Derived>
  Eigen::DenseBase<Derived> const & add(Eigen::DenseBase<Derived> const & chunk, size_t repeat = 1)
  {
    using Eigen::numext::floor;

    auto chunkScaled = (chunk.derived().array() * binDensity - binStart + 0.5).floor().eval();

    int minIdx = chunkScaled.minCoeff();
    int maxIdx = chunkScaled.maxCoeff();
    if (minIdx < 0)
    {
      size_t shift = -minIdx;
      growByStart(shift);
      chunkScaled += shift;
      minIdx = 0;
      maxIdx += shift;
    }
    growToEnd(maxIdx + 1);
    stretchNonzero(minIdx, maxIdx);

    this->addSums(chunk.size() * repeat, chunk.sum() * repeat, chunk.derived().matrix().squaredNorm() * repeat, (chunk.derived().array().abs2() * chunk.derived().array().abs2()).sum() * repeat);

    bins[chunkScaled.redux([this, repeat](T lastIdx, T nextIdx) -> T {
      bins[lastIdx] += repeat;
      return nextIdx;
    })] += repeat;

    return chunk;
  }

  void add(StatsAccumulatorHistogram<T,OPTIONS> & group)
  {
    if (group.size() == 0)
    {
      return;
    }
    growBoth<OPTIONS>(group);
    stretchNonzero(group.minNonzero, group.maxNonzero);
    this->addSums(group.size(), group.sum1(), group.sum2(), group.sum4());
    bins += group.bins;
    assert(this->size() == bins.sum());
  }

  void remove(StatsAccumulatorHistogram<T,OPTIONS> & subgroup)
  {
    this->setSums(this->size() - subgroup.size(), this->sum1() - subgroup.sum1(), this->sum2() - subgroup.sum2(), this->sum4() - subgroup.sum4());
    growBoth(subgroup);
    bins = (bins - subgroup.bins).round();
    assert(this->size() == bins.sum());
    if (this->size() == 0)
    {
      return;
    }
    while (bins[minNonzero] == 0)
    {
      ++ minNonzero;
    }
    while (bins[maxNonzero] == 0)
    {
      -- maxNonzero;
    }
  }

  // Assuming this distribution is made by summing data with another distribution,
  // estimate a distribution of the unsummed data (the data summed with the other to
  // produce this).
  template <StatsOptions OPTIONS2>
  StatsAccumulatorHistogram<T, OPTIONS2> withRemovedFromSum(StatsAccumulatorHistogram<T, OPTIONS2> & componentDistribution)
  {
    // TODO: this function does not succeed when tested.  it's notable that it appears that performing this calculation in the forward direction is much more robust.
    //        -> mean of result looks reasonable; variance does not.   might have been testing with too small a dataset and with artefacts in recording not accounted for
    // TODO: negative values in resulting histogram should be removed

    // Convolution:
    //  - f conv g (n) = sum_m f(m) g(n - m)
    //                 = sum_m f(n - m) g(m)
    //  - equal to pointwise product of ffts [Convolution Theorem]
    //  - fft of convolution is called Cross-Energy Density Spectrum
    //  - done of histograms, produces distribution of sum after sampling each distribution once
    //  - in time domain, discrete index 0 of input and output is value 0.0
    //  - can be reversed using pointwise quotient
    // Correlation:
    //  - f corr g (n) = sum_m f*(m) g(m + n)
    //  - equal to pointwise product of fft 1 and complex conjugate of fft 2 [Correlation Theorem]
    //  - done of histograms, produces distribution of difference after sampling each distribution once
    //  - in time domain, discrete index of output is shift required to make two input buffers match
    //        -> taking the other conjugate results in negative output indices
    // Buffers are circular, so "negative" values are at tail, with -1 at end.

    if (binWidth != componentDistribution.binWidth)
      throw std::logic_error("distributions have mismatching bin widths");

    // uses convolution theorem in reverse to deconvolve one from another
    std::cerr << "Component bins sum: " << componentDistribution.bins.sum() << std::endl;

    // 1. size these bins to match and hold result
    // this - comp = ret
    // the smallest item in ret could be the smallest item in this minus the largest item in comp
    auto smallestIdx = minNonzero - (componentDistribution.maxNonzero + componentDistribution.binStart);
    auto smallest = smallestIdx + binStart;
    if (smallest > 0)
    {
      smallest = 0;
      smallestIdx = -binStart;
    }
    if (smallestIdx < 0)
    {
      growByStart(-smallestIdx);
      smallestIdx = 0;
    }
    // the largest item in ret could be the largest item in this minus the smallest item in comp
    auto largestIdx = maxNonzero - (componentDistribution.minNonzero + componentDistribution.binStart);
    auto largest = largestIdx + binStart;
    if (largest < 0)
    {
      largest = 0;
      largestIdx = -binStart;
    }
    if (largestIdx >= bins.size())
    {
      growToEnd(largestIdx + 1);
    }

    growBoth(componentDistribution);

    if (largestIdx != bins.size() - 1)
    {
      largestIdx = bins.size() - 1;
      largest = largestIdx + binStart;
    }
    if (binStart != smallest)
    {
      smallest = binStart;
      smallestIdx = 0;
    }

    size_t bufsize = largest - smallest + 1;
    assert(bufsize == bins.size());
    assert(smallest == componentDistribution.binStart);
    HeapVector<Scalar> buf(bufsize), otherbuf(bufsize);

    // stats histogram buffers for convolution fft align value zero to index zero
    buf.head(largest + 1) = bins.segment(-binStart, largest + 1);
    buf.tail(-smallest) = bins.segment(smallest - binStart, -smallest);
    otherbuf.head(largest + 1) = componentDistribution.bins.segment(-componentDistribution.binStart, largest + 1);
    otherbuf.tail(-smallest) = componentDistribution.bins.segment(smallest - componentDistribution.binStart, -smallest);
    
    // 2. fft both
    Eigen::FFT<Scalar> fft (Eigen::default_fft_impl<Scalar>(),
       typename Eigen::FFT<Scalar>::Flag(Eigen::FFT<Scalar>::Unscaled | Eigen::FFT<Scalar>::HalfSpectrum));

    HeapVector<typename Eigen::FFT<Scalar>::Complex> fftResult, fftOther;
    StatsAccumulatorHistogram<T, OPTIONS2> ret (binWidth, binDensity, smallest);
    std::cerr << "Our bins sum: " << buf.sum() << std::endl;
    std::cerr << "Component bins sum: " << otherbuf.sum() << std::endl;
    std::cerr << "Product sum: " << (buf.array() * otherbuf.array()).sum() << std::endl;

    fft.fwd(fftResult, buf);
    fft.fwd(fftOther, otherbuf);

    // 3. take the quotient/product
    // TODO: is it worthwhile to use sparse division trick to map data from a compressed sparse array
    //          https://stackoverflow.com/questions/36226818/eigen-c-sparse-matrix-elementwise-product-and-divide
    //   but note that some inf arise from small nonzero values in denominator
    fftResult.array() /= fftOther.array();
    fftResult.array() = (fftResult.array() - fftResult.array()).isNaN().select(0, fftResult.array());

    std::cerr << "fftResult sum: " << fftResult.sum() << std::endl;

    // 4. ifft quotient/product
    fft.inv(buf, fftResult, bufsize);
    ret.bins.conservativeResize(bufsize);
    ret.bins.head(-smallest) = buf.tail(-smallest);
    ret.bins.tail(largest) = buf.head(largest);

    std::cerr << "ret.bins sum: " << ret.bins.sum() << std::endl;

    // 5. get effective sample size to be the smaller one (could also sqrt I guess but then we'd have a non-integral samplesize)

    size_t smallCount, largeCount;
    if (this->size() < componentDistribution.size())
    {
      smallCount = this->size();
      largeCount = componentDistribution.size();
    }
    else
    {
      largeCount = this->size();
      smallCount = componentDistribution.size();
    }
    ret.bins *= smallCount / ret.bins.sum();

    // 6. return a histogram using result as bins
    // rebuild sums & count
    if (!sumCoeffs)
    {
      if (componentDistribution.sumCoeffs)
      {
        sumCoeffs = componentDistribution.sumCoeffs;
      }
      else
      {
        generateSumCoeffs();
        componentDistribution.sumCoeffs = sumCoeffs;
      }
    }
    ret.sumCoeffs = sumCoeffs;
    ret.updateSums(smallCount);
    
    return ret;
  }

private:
  T binWidth;
  T binDensity;
  T binStart;
  HeapVector<T> bins;

  size_t minNonzero, maxNonzero;

  std::shared_ptr<HeapVectors<T,3>> sumCoeffs;

  StatsAccumulatorHistogram(T binWidth, T binDensity, T binStart)
  : binWidth(binWidth),
    binDensity(binDensity),
    binStart(binStart)
  { } 

  template <StatsOptions OPTIONS2>
  void growBoth(StatsAccumulatorHistogram<T,OPTIONS2> & other)
  {
    if (binStart < other.binStart)
    {
      other.growByStart(other.binStart - binStart);
    }
    else if (other.binStart < binStart)
    {
      growByStart(binStart - other.binStart);
    }

    if (bins.size() > other.bins.size())
    {
      other.growToEnd(bins.size());
    }
    else if (other.bins.size() > bins.size())
    {
      growToEnd(other.bins.size());
    }
  }

  void stretchNonzero(size_t possibleMinIdx, size_t possibleMaxIdx)
  {
    assert(possibleMinIdx < bins.size());
    assert(possibleMaxIdx < bins.size());
    if (this->size() == 0)
    {
      minNonzero = possibleMinIdx;
      maxNonzero = possibleMaxIdx;
      return;
    }
    assert(minNonzero < bins.size());
    assert(maxNonzero < bins.size());
    if (possibleMinIdx < minNonzero)
    {
      minNonzero = possibleMinIdx;
    }
    if (possibleMaxIdx > maxNonzero)
    {
      maxNonzero = possibleMaxIdx;
    }
  }
  
  void growByStart(size_t shift)
  {
    sumCoeffs.reset();
    size_t oldSize = bins.size();
    assert(oldSize + shift <= 384); // TODO: remove for non-rtl-sdr
    bins.conservativeResize(oldSize + shift);
    bins.tail(oldSize) = bins.head(oldSize).eval();
    bins.head(shift).setZero();
    assert(binStart - shift >= -256); // TODO: remove for non-rtl-sdr
    binStart -= shift;
    minNonzero += shift;
    maxNonzero += shift;
  }

  void growToEnd(size_t size)
  {
    if (size > bins.size())
    {
      size_t oldSize = bins.size();
      assert(size <= 384); // TODO: remove for non-rtl-sdr
      sumCoeffs.reset();
      bins.conservativeResize(size);
      bins.tail(size - oldSize).setZero();
    }
  }

  void generateSumCoeffs()
  {
    sumCoeffs.reset(new HeapVectors<T,3>());
    sumCoeffs->resize(bins.size(), 3);
    for (size_t i = 0; i < sumCoeffs->rows(); ++ i)
    {
      T value = (binStart + i) * binWidth;
      for (size_t j = 0; j < 3; ++ j)
      {
        (*sumCoeffs)(i, j) = value;
        value *= value;
      }
    }
  }

  void updateSums(size_t size)
  {
    // bins has n rows
    // sumCoeffs has n rows and 3 columns
    auto sums = (sumCoeffs->array().colwise() * bins.array()).colwise().sum().eval();
    this->setSums(size, sums[0], sums[1], sums[2]);
  }
};

template <typename T, StatsOptions OPTIONS = STATS_SAMPLE | STATS_ASSUME_NORMAL>
class StatsDistributionBrief : public StatsDistributionByMeasures<T, StatsDistributionBrief<T, OPTIONS>, OPTIONS>
{
  using Measures = StatsDistributionByMeasures<T, StatsDistributionBrief<T, OPTIONS>, OPTIONS>;
public:
  using Scalar = T;

  StatsDistributionBrief(Scalar mean, Scalar variance)
  : Measures(mean, variance)
  { }

  using Measures::setMeasures;
  using Measures::setMean;
  using Measures::setVariance;
};

template <typename T, StatsOptions OPTIONS = STATS_SAMPLE | STATS_ASSUME_NORMAL>
class StatsAccumulatorNormal : public StatsDistributionBySums<T, StatsAccumulatorNormal<T, OPTIONS>, OPTIONS>
{
public:
  using Scalar = T;

  static_assert(OPTIONS & STATS_ASSUME_NORMAL, "StatsAccumulatorNormal must assume normal.");
  //static_assert(!(OPTIONS & STATS_INFINITE), "StatsAccumulators cannot be infinite.");

  StatsAccumulatorNormal()
  { }

  void clear()
  {
    this->setSums(0, 0, 0);
  }

  template <typename Derived>
  void add(Eigen::DenseBase<Derived> const & chunk, size_t repeat = 1)
  {
    // note: sum() and squaredNorm() do the right component-wise thing for vectors and for
    //       matrices: i.e. treating every entry as simply an item to add to the stats
    this->addSums(chunk.size() * repeat, chunk.sum() * repeat, chunk.derived().matrix().squaredNorm() * repeat);
  }

  void add(Scalar const & value, size_t repeat = 1)
  {
    auto square = value * value;
    this->addSums(repeat, value * repeat, square * repeat);
  }

  void add(StatsAccumulatorNormal<T, OPTIONS> const & group)
  {
    this->addSums(group.size(), group.sum1(), group.sum2());
  }

  void remove(StatsAccumulatorNormal<T, OPTIONS> const & subgroup)
  {
    this->setSums(this->size() - subgroup.size(), this->sum1() - subgroup.sum1(), this->sum2() - subgroup.sum2());
  }

  [[deprecated("learn statistics and improve library")]]
  StatsAccumulatorNormal<T, (OPTIONS & ~(STATS_SAMPLE)) | STATS_INFINITE> const &
  fakeInfinitePopulation() const
  {
    return reinterpret_cast<decltype(fakeInfinitePopulation())>(*this);
  }

#if 0
  T errorInAverage() const
  {
    using Eigen::numext::sqrt;

    // the variance of an averaged set of samples is the variance of one
    // sample divided by the number of samples taken
    //
    // this is the same thing as the variance of the sampling distribution of means, and hence represents the error in the expected value
    return sqrt(variance() / count) * cutoff;
  }

  // TODO AFTER NEXT: rather than this confusing function, provide a method
  // to get the whole sampling distribution of the error, and it can be queried
  // for its error.
  T errorInErrorInData() const
  {
    using Eigen::numext::sqrt;

    // 1. first we get the variance of the sampling distribution of standard deviation
    // std dev2 = std dev1 / sqrt(2 N)
    // variance2 = variance1 / (2 N)
    // 2. multiply it by cutoff, since errorInData is multiplied by cutoff
    // 3. turn that into another standard error by taking the square root and multiplying that by cutoff
    //
    T samplingVariance = sqrt(variance() / (2 * count) * cutoff);

    // it's notable that these metrics are roughly multiples of each other
    // errorInData = sqrt(v) * cutoff
    // errorInAverage = sqrt(v / count) * cutoff
    //                = errorInData / sqrt(count)
    // samplingVariance = sqrt(v / (2 * count)) * cutoff
    //                  = errorInData / sqrt(2*count)
    //                  = errorInAverage / sqrt(2)

    return sqrt(samplingVariance) * cutoff;
  }

  //T const & average() const
  //{
    // sampling std dev = std dev / sqrt(count)
    //                  = sqrt(variance / count)
  //}

  // takes a bunch of statsaccumulators, subtracts their differing means,
  // and determines their error as if treated as one dataset
  template <template Container>
  T groupErrorInAverage(Container<StatsAccumulatorNormal<T>> const & group) const
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
    for (StatsAccumulatorNormal<T> const & accum : group)
    {
      numerator += accum.dataSquaredSum + accum.mean * (accum.mean - 2 * accum.dataSum);
      denominator += accum.count;
    }
    T variance = numerator / (denominator - 1);
    return sqrt(variance / denominator) * cutoff;
  }

  StatsAccumulatorNormal<T> operator*(T scalar) const
  {
    auto absScalar = abs(scalar);
    return {mean * scalar, dataSum * absScalar, dataSquaredSum * absScalar, count};
  }

  StatsAccumulatorNormal<T> operator+(StatsAccumulatorNormal<T> const & other)
  {
    // what's the distribution of the sum of two distributions, given that we only have
    // a sample for both populations?
    //
    // two distributions summed have a sum mean and a sum variance ...
    // it's clear that if both have the same count, then the count of the final is this
    // count.
    // it also seems valuable to see if we could model the datapoints being summed under
    // such conditions.
    //
    // so, it's notable that the population variance sums, so the correct resulting variance is best approximated by the sum of the component variances, withou regard to the count.
    // this is likely true of the mean too.
    // 
    // it's also notable that each of these variances is a sample from the sampling distribution of variance
    // the count would be used to determine this.  we're looking here for the error in the metrics themselves, not in the data.  so i'll need to add methods for that.
    //
    // each distribution has a standard error with regard to its mean and its variance
    // what is the standard error of the summed mean and variance?
    //  the variance of a summed statistic is the summed variance
    //  because the variance of each statistic includes the count, the combined count is included
    //  in the sum
    // TODO: this should be pretty obvious
    //          add methods for standard error of mean and variance
    //          update data for them properly in this sum
    //          we might be able to share the count to make it work, or we might need
    //          a separate count to track this.
    //          -> I got really confused working towards resolving the above line,
    //             and ended up not resolving it.
    //
    // TODO: provide the body of this function
    //    The distribution spewed out by this function doesn't use the normal
    //    statistic bodies unless the proper count is solved for to make them
    //    work out.
    //    I got very confused, but I'm guessing a way to solve this would be to
    //    create a new class that wraps the two parent classes or stores special
    //    data.
    //    One might also be able to jump straight to the statistics in question,
    //    which would be the error of the sampling distribution of the error of
    //    the sum, and just plain the error of the sum (which is the average of the sampling distribution of the error, naturally).
    //
    //    okay ... I want the missing distribution.
    //    I have a distribution resulting from its sum, to repeat again to help
    //    think:
    //    A + B = C
    //    I have A and C, what is B?
    //
    //    B is a distribution produced by an appropriate transform on the particular
    //    distribution class.  For histogram, pointwise division to make a pdf with
    //    minimum count.  For normal, the variances will actually subtract,
    //    because shat_A + shat_B = shat_C.
    //    I wonder what approach there might be to generalize the transform
    //    of solving for a sum from a distribution ...
    //
    //    There are two approaches to subtracting a distribution:
    //      - find me a distribution that would be given if we sampled
    //        these two and subtracted the result
    //      - find me a distribution that would make this distribution if
    //        we sampled it and one other and added
    //
    //      why are these different for normal but not for histogram?
    //
    //      I convolved two distributions and found that pointwise product
    //      produced approximate variance sum, wherease pointwise division
    //      produced approximate variance difference.
    //
    //      I created three A + B = C normal distributions and found that
    //      the variance summed with the data operators.
    //
    //      I'm thinking the operators here will need to two things:
    //      - whether one distribution is being combined with, or removed from.
    //        the other distribution
    //      - the sign of the underlying data operation
    //
    //      So we could calculate A - B to figure out a distribution resulting
    //         from the subtraction of samples from A and B
    //         -> data operation is a - b
    //         -> distribution operation is combination
    //      Or we could calculate A - B to figure out a distribution used 
    //         to sum with B to produce A
    //         -> data operation is a + c
    //            i.e. a - b
    //         -> distribution operation is removal of B from A
    //      So there are two different signs when distributions are added to
    //        or subtracted from each oher: the sign of the mean, and the sign
    //        of the variance.  The mean takes the sign of the data operation,
    //        and the variance takes the sign of the combination or removal.
    //
    // It would be slick to just put operator+ and operator- on dists but
    // won't work.
    //
    // We could do it in two steps.  First negate the distribution to get the
    // data sign right, then add or subtract to get the variance sign.
    // This would work and be nice to do, but it would be confusing later
    // because the meaning is not fully held by the symbology.
    //
    // I guess variance addition/subtraction should be done with
    // combine/removal functions.  The only operators available to distributions
    // should be unary.

    // calculate metrics for the sum of two distributions
    // 
    // because the signal is noise, in the final sum we're going to be looking
    // at errorInData as the value in question.
    //
    // That might need some changes to the code, but we'll then want to look at
    // the error in that metric, the errorInErrorInData, for which we'll need
    // the variance of the sampling distribution of standard deviations.
    //
    // I think the solution to the 'changes in code' question is to make a method
    // that returns the distribution of the error.
    //
    // First we need the distribution of the sum.
    //
    // ====
    // OKAY yeesh
    // The means will add, and so will the variances.
    // The mean is actually a member variable.  If streaming happens, this mean
    // is calculated by dataSum / count.  Have to figure out if I want to
    // provide for streaming.
    //
    // variances should add.  I could make a new class for this ... that's an
    // option.
    // (A1 - B1 * mean1) / C1 + (A2 - B2 * mean1) / C2
    // C2 (A1 - B1 * mean) / C1C2 + C1 (A2 - B2 * mean) / C1C2
    //
    // To preserve the variance, an equivalently sampled imaginary distribution
    // can be produced, that has
    // dataSquaredSum = (count2 - 1) * dataSquaredSum1 + (count1 - 1) * dataSquaredSum
    // count = (count2 - 1) * (count1 - 1) + 1
    //
    // and then the mean and dataSum
    // dataSum * mean = (count2 - 1) * dataSum1 * mean1 + (count1 - 1) * dataSum2 * mean2
    
    // This production of an imaginary distribution may be doable, but we're
    // going to need the variance of the standard deviation of this distribution.
    // This is supposed to have a correct count ...
    //
    // what's the sampling distribution of the standard deviation of a sum of two distributions?
    // this sounds really familiar I may have just looked this up
    // is this the sampling distribution of differences and sum? it is !!
    // or is it?  
    // THE VARIANCE OF THE STATISTIC OF THE DIFFERENCE IS EQUAL TO THE SUM OF
    // THE INDIVIDUAL VARIANCES OF THE STATISTIC
    
    // I seem to be having an issue holding this task in my working memory.
  }

  StatsAccumulatorNormal<T> operator-() const
  {
    return *this * -1;
  }
#endif

private:

  StatsAccumulatorNormal(size_t count, T sum1, T sum2)
  : StatsDistributionBySums<T, StatsAccumulatorNormal<T, OPTIONS>, OPTIONS>(count, sum1, sum2)
  { }
};


/*
 * Hmm okay ... I'm a little confused by variance.
 *
 * I have a population of samples.  There are two models for it: a constant plus
 * noise, and noise plus noise.
 *
 * There are two error metrics:
 *    the variance / standard error of the data
 *    the variance / standard error of the mean of the data
 *
 * The first one represents how widely the data spreads around the mean.
 * For the constant case, this would represent the loudness of the noise.
 * The the combined case, this would be used to find the loudness of the noise.
 *
 * The second one represents how accurate the mean is.
 * For the constant case, this would represent the error in the measurement,
 * which is naturally related to the loudness of the noise but naturally
 * decreases as the sample size increases; which the variance does not.
 * For the noise-noise case, this doesn't seem as useful. ... there
 * we're concerned with the variances themselves, rather than the mean.
 *
 * Perhaps it would be safest to use statistical words in my methods.
 *
 */


// Task:
//    Goal: Given samples from normal distributions A and A + B, identify
//        both the standard error and the standard error of the standard error,
//        of distribution B.
//
//   If we call A + B distribution C,
//   then B = C - A
//
//   mu_B = mu_C - mu_A
//   sigma^2_B = sigma^2_C + sigma^2_D
//
//   (sigma^2)_sigma_B = (sigma^2)_sigma_C + (sigma^2)_sigma_A
//   
//  now I'll need to find (sigma^2)_sigma_C ... these are the variances of the
//  error.  I've calculated them already, but I should step back and decide how
//  to organize them.
//
//    Goal 2: Generalize the StatsAccumulator class to handle distributions
//    of sampled data, sampling distributions of their statistics, and sampling
//    distributions of their sums.
//    Ideally make it general enough to continue finding the sampling distributions of further statistics.
//
//    So a major issue here is what to store when we do a sum.
//    A sum of two distributions produces a distribution with all its attributes
//    defined:
//      - mean is the sum of means
//      - variance is the sum of abs(variance)s
//      - sampling dist can also be made:
//        - variance is sum of variances of individual sampling dists
//        - mean is sum of means of individual sampling dists
//
//    So a separate class might be appropriate, although we could solve between
//      the sampling distributions of differences and the sampling distributions
//      of the mean and variance etc to find the right underlying values, maybe
//
//    Shouldn't be hard because the definitions are simple.
//    The variances add in absoluet value, and the means add.
//    So we want to track variance and mean, and hold those as primary members.
//
//    There's an issue here ... to get the sampling distribution of the standard
//    error, we divide the variance by the count ... but a sampling distribution
//    doesn't really have a count.
//    
//    So the question is, what is the sampling distribution of the standard error
//    of a sampling distribution of the standard error?
//    Maybe kind of.  I seem kind of co
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//    Provide code general enough such that a statistic on a statistic
//    can be treated via the same paths as a statistic on the data itself.
//
//    One approach is to parameterize which statistic is pulled from the class
//    that tracks the data.  Another approach is to provide a method which returns
//    another object, that has the same methods but gives the statistics of the
//    statistic, rather than the statistics of the data.
//
//    The second approach is more general and likely better.  However, I'm having
//    trouble holding the information in my mind well enough to understand if I
//    can do both tasks with one class, which could be preferably, or if I will
//    need a second class with different code to handle the statistic of the
//    statistic.  It is certainly possible to do it in a completely general way,
//    but I don't want to do anything that requires writing a compile-time
//    expression evaluator.  This would take weeks at least.
//
//    It's notable that the statistics of a normal distribution I believe have
//    a normal distribution, so this should be pretty doable.  I just need to
//    figure out which member variables to store the crucial information in, to
//    hold this.  Another question is how to handle data streamed in ...
//
//    To manage memory, the current sample distribution class doesn't store every
//    data point, it just accumulates information as data comes in, in smaller
//    buffers.  It keeps running totals and updates the stats each buffer, without
//    having to hold the previous buffers. (obviously this was done by algebraic
//    manipulation of the expressions for the statistics of mean and variance)
//
//    Okay, I'm off in a 'diatribe', now.  I'm having trouble continuing this goal.  Let's update it.
//
//    
//
//    sampling distribution stats:
//      stderr of sigma^2 = sigma^2 sqrt(2/N)
//      stderr of sigma = sigma / sqrt(2 N)
//      stderr of mean = sigma / sqrt(N)
//      variance of sigma^2 = sigma^2^2 2 / N
//      variance of sigma = sigma^2 / (2 N)
//      variance of mean = sigma^2 / N
//
//
//    The variance of a sampling distribution of a distribution of differences
//    Is just the sum of the difference of the sampling distributions.
//
//    It's notable that a distribution of differences, has a variance itself
//    equal to the sum of the variances of the distributions of the different
//    populations.
//
//    So it seems
//    var(dist(var(a + b))) = var(dist(var(a)) + dist(var(b)))
//
//    of course we can't just add the underlying distributions because this
//    statistic works only for the variance, and a different addition would
//    be needed for the mean, for example.
//
//    If we are a sum statistic, to calculate our sampling distributions, it looks
//    like we will need to access the distributions that we are a sum of.
//
//    why can't we just use the sampling distribution stats on our variance and
//    mean ???
//
//    hmmm maybe it all works out =)
//
//    variance of sigma^2 = sigma^2^2 2 / N
//    variance of sum of sigma^2 = sigma1^4 2 / N1 + sigma2^4 2 / N2
//
//    sigma^2 of the sum = (sigma1^2 * N1 + sigma2^2 * N2) / (N1 + N2)
//    variance would be = big^2 2 / N
//
//    I see.  We'd have to calculate count to make this equation be correct.
//    And if we added new statistics, they might be wrong if they use a count.
//    Of course there are very few statistics to add to a normal distribution.
//    But even like say the variance of sigma instead of sigma^2 ... this may
//    work out or it may not!  I don't know!
//
//    I'm thinking since things are hard right now, I won't solve it out.
//    Let's see what other options there are.
//
//    I could make the fucntion return an objet that references its parents.
//
//    An _easier_ approach would be to make a simple class that just stores
//    the metrics, but doesn't update when the source classes are updated.kjj
//
//    Say I make a new class, should it wrap the parent classes, or just
//    store data?
//
//    References:
//        - object is created at the start.  memory allocation is delicate.
//        - execution path involves calculation every request, from the references
//        - only stats requested are calculated
//        - can generically apply to nested references
//
//    Flat data:
//        - lightweight, simple to make, memory allocation straightforward
//        - all stats prepared when object made, must be made every request (light)
//        - i'm too confused to determine immediately if this can do nested references.  do we need them?
//      
//   A way to have flat data more normal could be to have it be the only way
//   to access the things in question.
//
//   There are a few kinds of normal distributions.  Sampled data.  Population data, infinite or finite.
//
//   Hmm, so everything we have is based on one of those distributions.
//   I'm basically sampling infinite data here.
//

// Okay ... my algebraic idea wasn't that great
//
// There are two things we want atm, function + noise and noise + noise
// for function + noise, we want the mean and std error, to extract function,
// as we can multiply sample the same spot on function repeatedly.
//
// for noise + noise, we want the std error, and std error error, because
// we care about the magnitude of the noise, and it can be found statistically
// pretty accurately.
//
// To find noise + noise, we need to create the missing noise distribution.
// Its variance is the difference of the sum and the component.
// Its mean is zero.
// We'll also want the variance of its standard error, which can be calcualted.
//
// So, what structure do I provide to find these?
// Maybe a different stats accumulator.
// 
// 1. Stats Accum for mean and error
//
// 2. Stats Accum for noise difference?  It would accumulate, and reference
//    another accum to get the difference?
// 
// Since we dont' know which bin is the maximum, mgiht be better to have
// a class that references two others and finds the particular metrics we need.
// ^================== ANSWER PERHAPS
//
// There's another way of looking at this that is more accurate/powerful but
// harder and slower.  The distributions may not be gaussian and the constant
// signals may not be constant.  I can estimate each distribution by making
// a histogram, and then convolve the distributions to perform arithmetic on them.
//
// Would this solution solve both problems at once?
// Would it at least be usable for noise+noise ?
//
// I think this could work.
//
// I'd make a stats accumulator that accepts 8-bit samples, or uses 256 histogram buckets.
// It fills them as it goes along, and it integrates to identify the mean and std
// dev.  It can now use more advanced metrics like kurtosis.
//
// I can identify a missing distribution by de-convolving it.
// I FFT the sum distribution and the individual distribution, and divide point-
// wise the larger by the smaller.  I can then apply the results to this.
//
// What about the sampling distribution?  How do I make use of this?
// I think I can use the same functions from the book, but I'm going to need to
// understand the more advanced moments of the data.
//
// Specifically, I'll need the variance of the sampling distribution of the standard error ... or maybe some other metric.


// [ ] check the definition of sampling distribution of standard error for non-gaussian distribution, read about moments
//
// what do we want from these distributions?
// (y'know it's kind of interesting that if I did this pointwise division, I could have a sample count in the result, which is weird because I'm not sure what that means, but I could totally integrate it and see the # of 'samples')
// So once we have the addend distribution, we basically want to know how wide it is, and the standard error will work fine for this.  This tells us how loud the noise is.
// But then we want to know how accurate this measure of loudness is.  So I'll want the sampling error ... it could probably be composed based on the source distributions, which kind of brings us back to square 1.

// 2 ways to get standard error of standard error: from direct formula including number of samples, or from difference or sum of others.  B does not have a sample count; we'll want to use a difference or sum.  I'm not sure if this error here should be larger or smaller than the two its made from, but I guess I'd expect it to be larger here.
// I'm thinking the standard error of the standard error of B would be similar to A + C
// But the plain old error of B would be more like C - A

// => The sampling distributions of sums and differences give the variance of
//    the calculated sums and differences of _the_statistics_themselves_: so
//    for example you can determine the variance of the distribution of the
//    sum of two variances.  It does _not_ mean that distribution is one of
//    sums or differences of other distributions, but could be used in the
//    calculation of such.
//
// Basically I'm concerned with the accuracy of my estimated distribution.
// I can probably use the standard deviation, and then use the standard error
// of the standard deviation.  It just involves some moments.
//
// A moment is just the power the data is raised to before being summed.
// So my dataSum is the first moment (times count)
// and dataSquaredSum is the second moment (times count)
// The variance is the second moment about the mean, because the mean
// is subtracted prior to the power.
//
// So, of course, the standard error of the standard deviation requires the
// number of samples to calculate.  (I COULD JUST USE THE MINIMUM)
// I'm curious how it would go with the histograms.
//
// Histogram A has Na samples.
// Histogram C has Nc samples.
// These are found via integration or tracking.
//
// We find Histogram B by deconvolving A and from C.
// Ha, this totally removes the count information. (also the integral of
// the spectrum is different from that of the data, but on the same order). 
// I'm thinking this convolution likely applies to the pdf, whose integral is
// 1.0 .
//
// We can probably just use the minimum count.
//
// Pointwise product and division does appear to work.




// i'm thinking about calculating the error.
// The approach I came up with was to assume the distribution really was normal
// about zero, and to calculate the likelihood that, if sampling such a
// distribution, we would get a mean closer to zero than the one we calculated.
//
// This sounds like a pretty meaningful approach; it just leaves out all
// the information we have that indicates the ways in which the distribution
// is not normal.
//
// Really, perhaps what I'd want ideally, is to compare this distribution
// with a distribution that represents the background noise.  What's the
// likelihood of sampling that distribution, and instead getting this one?
//
// Since the mean is our value metric, it's probably meaningful to use it
// as the statistic used to determine this.
//
// But i guess if we're going to assume that our statistics from our sampled
// data are not quite right, we're going to need to assume something about that
// not-rightness.  Assuming normal is a reasonable approach to this.`
//
// What normal distribution do we assume in figuring our likely mean?
// We could also look at the distribution of the mean inferred from this
// particular distribution.
//
// The standard error for means is just sigma / sqrt(N).
//
// If we assume the mean is zero, than what is the standard deviation?
// We could infer it by expansion.

// variance = sum(x - mean) / (n - 1)
//          = sum(x) / (n - 1)
//          = 

// I'm just cnosidering a concern here ... what if the real data has a DC offset?
// Notably, our radio cannot represent zero, at least not with thei nterface I've
// made for it, so there will always be a small DC offset in the actual recordings.
// If I applied tihs metric to straight samples, it wouldn't work.  However,
// it should work for the checking that underlying noise.
//
// Say there were a constant signal, added to noise, added to a DC offset.
// 
// We'd then have some accumulated averages for the wave; we'd want to find the
// one that's most extreme from the noise average, and see its likelihood.
// So it's more meaningful to look at the average of the background noise
// and compare.
//
// Constant signal case:
// - stats accum for whole recording: get standard error of mean, use it to
//   compare extreme mean for sampled offset into signal.
// 
// Combined noise case:
// - stats accum for minimum: get standard error of standard deviation, use it to
//   compare standard deviation of inferred signal (or of maximum, really)
//
// so really we want to reference a different distribution and get the standard error from its mean.

// okayyyyyyy
// I have two distributions, A and B.  
// Y'know, if B doesn't exist then it should be zero!  All zero!
// But I created B by remove samples of one distribution from another.
//
// hmmmm my function is called approximateSignificanceOfMeanGivenBaseline
// it compares mean to baseline mean, and divides by the standard error
// of the means sampling distribution of the baseline.
// So, the mean difference divided by the standard error is the z-score.
// That does show how our mean differs from a baseline.
//
// I don't tihnk I even needed to convolve to use this function.
// I can apply it to the standard deviation sampling distributions of the
// combined signal and the baseline signal.
//
// The sampling distribution of the standard deviation of the baseline signal
// show what standard deviations we are likely to get.
// Then the standard deviation of the combined signal is the standard deviation
// we sampled specially.  Sounds good, I think ...
// Well, the thing is that we are calling two sampling distributions.
// First I get the sampling distribution of the standard deviation.
// Then I get the sampling distribution of its mean.
// So the standard error I'm dividing by is the standard error of the means of
// the standard deviations ... as opposed to the standard error of just the
// standard deviations, which would be directly appropriate since that
// is what I'm comparing.  This is only correct if the means
// sampling distribution of th standard deviations sampling distribution is
// the same .... and it won't be.
//
// So
//

// Okay, I'm considering it good to treat A as a population, look at its
// sampling distribution, and consider how likely B's statistic is given
// that sampling distribution.
//
// Of course, A and B are both just samples, and I think chapter 16 "analysis
// of variance" may cover this.
//
// deviation: difference from mean
// [total] variation: sum of square of deviation?
// standard deviation: sqrt variation

// okay, got it compiled
// I'm using a fakeInfinitePopulation function in c++-style to remind me
// that I'm fudging the statistics
// it casts a template param to create a view of a class as infinite pop
// but this breaks static_cast in the derived() function because he
// derived class still believes it is a finite sample
// I guess I should make a lightweight view
// hrm that means the view can't call derived functions =S which might not
// be needed but is unneccessarily limiting

// okay there could be a delay bug after the first buffer unless I fix this
// I'm pushing to the back of these heapvectors but they prform a memory
// allocation every time I do that, unlike std::vectors which only do it
// occasionally but aren't usable for quick math unless wrapped in eigen
// to be used this way
// maybe I'll use std::vectors but wrap them in eigen vectors

// Accumulator problem:
// I'm using a histogram stat, I think, to calculate a mean period
// the mean that spews out of the stat accumulator is incorrect, off by a little
//
// - is this from error accumulating many integers in a double?
// - if not, what is it from?
// - will this affect the calculated model waveform, which uses accumulation like this too?
//    -> errors could be considered noise ... just decreases strength of signal
// - [ ] add todos noting this error, keep going

// thoughts regarding fixing accumulation precision:
// note: accumulating many floats will lose precision, almost always
// - could use integers to sum, uint64_t, and throw on overflow
//   samples from hardware have finite precision, which would make this reasonable
// - could use rational numbers, only issue is there are unlikely to be fast libraries for doing calculations
//   with them.  the code needs to be able to complete the analysis of a buffer before the next buffer arrives,
//   basically
// - could use integers combined with rational numbers ....
//
// - i could make a class for a vector of rational numbers.
//   there are only a few operations i need:
//   - summation
//   - division
//   - small integer powers
//   basic arithmetic.
//    -> note: it could make sense to make a vector stats accumulator with this; would also speed things up a lot
//   probably not needed yet.
// let's try throwing rational numbers at eigen and see what happens
// could we be backendy in a ... the approach karl desires is unreasonable for his goal i think
// yeah makes sense, that should be done if required, the work involved to test it will be needed anyway
//
// since i'll likely be coding some new math backend, i should switch to boost and to gnuradio
// and merge my work with more general approaches.
//
// what kinds of generality will i want to code in a way that combines with other things I could write,
// or to make it usable by other people?
// - algorithmic plugging?
//    i'd like to separate out the common pieces of my algorithms
//    so that when i am researching data analysis approaches using heuristics, I can replace different
//    pieces of the heuristics without having to start over
//    i'd also like to parameterize these things so that intelligent code can try different approaches
//    i was doing this previously for finding weak signals without profiling.
//    i came up with many common algorithms, broke them into 3 conceptual pieces, to make the handling
//    pluggable.  i'm pretty sure my framework would have provided for my current approach to plug right in.
// X- root integration?
// - -> use gnuradio, soapy, pothos, and/or osmosdr
// - -> use BOOST stats
//
// - -> use BOOST accumulators
//          see plan.txt for jumping-off points
// - -> make a vector math wrapper to plug in vector math libs or hw acceleration
//
// pluggable concepts for analysis frameworks:
//   algorithms for:
//    - finding optima by minimizing or maximizing a parameter
//    - given data desired, and data available, what data to provide
//          -> math takes microseconds, which are in limited supply
//          -> we want to perform calculations on data that is judged to be valuable
//   norms:
//    - providing for pluggable access to algorithm parameters
//    - specifying what data is needed next
//    - specifying what data is needed to start
//
// so, current convolution approach fits into some of these, and may indicate more
// we have a huge stream of incoming data offered in dense, large buffers
// -> buffers need to be a certain size to be useful
// we downsample the data (pluggable algorithm)
//   and then convolve the downsampled data (pluggable algorithm)
//
// downsampling reduces the data available by summarizing it
//
// we then upsample and repeat the downsampling, using the results of the earlier downsampling to inform
// what data to look at
//
// karl's tree approach likely achieved a similar thing, it just had a bottom-up view rather than a top-down
// one.  plusses and minuses to both
// or something like that
//
// so this periodfinder that convolves produces PROPOSALS of periods to try
// it has a TECHNIQUE for checking them, which is one-by-one
// and it has a TECHNIQUE for PRODUCING NEW PROPOSALS
//
// it also has a TECHNIQUE for MASSAGING THE DATA to be have DENSER MEANING, so that the proposal space is
// smaller.
//
// in general, profile and detection have some of these things ...
// profile tries to DETERMINE the WAVELENGTH and the ONSET TIMES
// detector NEEDS and DETERMINES the WAVELENGTH
// detector DETERMINES the ONSET TIMES and the POWER
//
// i'd like a more general class to let me replug some of these concepts or algorithms
// under the assumption i will run into more problems
// and to make the code reusable for other purposes
//
// tree approach was adaptive subsamping
// i think it kept statistics for chunks of the data, it GUESSED possible ONSET TIMES and used EVALUATIONS
// of those guesses to GUESS NEW onset times
// it likely used DENSE ACTIVE DATA BUFFERS and METRICS of old data
//    ... long vectors of metrics sometimes are used, summarizing old data
//
// maybe i'll just make a general class
// and use it to separate things out more and more as i encounter new approaches
// maybe try to put 2 different approaches in it
// if everything has a base class, then i can move concepts up towards the base as i use them
