#include <iostream>
#include <fstream>

#include "types.hpp"

#include "stats.hpp"

#include "rtlsdriqdump.hpp"

int main(int argc, char const * const * argv)
{
  StatsAccumulatorHistogram<Scalar> overall(RtlSdrIQDump::epsilon());

  {
    HeapVector<Complex> buffer(2048000 / 2);
    RecBufMeta meta;
    std::ifstream stream(argv[1]);
    RtlSdrIQDump data(stream);
    while (true)
    {
      data.readMany(buffer, meta);
      if (!buffer.size()) break;
      auto processed = buffer.real().eval();
      overall.add(processed);
    }
  }

  HeapVector<Complex> buffer(2048000 / 40 / 32);
  RecBufMeta meta;
  std::ifstream stream(argv[1]);
  RtlSdrIQDump data(stream);
  auto overallSTD = overall.standardDeviation();
  while (true)
  {
    data.readMany(buffer, meta);

    if (! buffer.size()) break;

    StatsAccumulatorHistogram<Scalar> stats(data.epsilon());

    auto processed = buffer.real().eval();

    stats.add(processed);

    //Scalar sample[2] = { stats.standardDeviation(), stats.mean() };
    Scalar sample = stats.standardDeviation() - overallSTD;

    std::cout.write(reinterpret_cast<char*>(&sample),sizeof(sample));
  }

  return 0;
}
