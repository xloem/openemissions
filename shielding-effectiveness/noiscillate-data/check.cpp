#include <iostream>

#include "types.hpp"
#include "stats.hpp"
#include "rtlsdriqdump.hpp"

int main()
{
  HeapVector<Complex> buffer(2048000 / 2);
  RtlSdrIQDump data(std::cin);

  while (true)
  {
    data.readMany(buffer);
    if (!buffer.size()) break;
    auto processed = buffer.real().array().abs();
    StatsAccumulatorHistogram<Scalar> stats(RtlSdrIQDump::epsilon());
    stats.add(processed);
    std::cout << stats.standardDeviation() << std::endl;
  }

  return 0;
}
