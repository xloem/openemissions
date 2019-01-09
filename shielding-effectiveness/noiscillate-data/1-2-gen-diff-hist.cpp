// 1-2-gen-diff-hist.cpp
//
// Given two histogram noise profiles of the same environment, with one emitter different betweeen them,
// produces an estimated histogram noise profile of the emitter itself.

#include <algorithm>

#include "processornoiseprofile.hpp"
#include "rootapplication.hpp"

using Scalar = double;

class T12GenDiff : public RootApplication
{
public:
  T12GenDiff(Int_t * argc, char ** argv)
  : RootApplication(
      "1-2-gen-diff-hist",
      "Loads two histogram profiles of the same environment and produces\n"
      "a profile of one emitter that changed between the recordings.",
      argc, argv)
  {
    regArg({{"-o"}, {"emitter.noisep"}}, {"output file to write difference in"},
      [this](std::string, std::string outfn)
      {
        _outfn = outfn;
      }, 1 // required
    );
    regArg({{}, {"1.noisep"},{"2.noisep"}}, {"two noise histograms to compare"},
      [this](std::string, std::string fname1, std::string fname2)
      {
        _infn1 = fname1;
        _infn2 = fname2;
      }, 1 // required
    );
  }

  virtual void Run(Bool_t retrn) override
  {
    ProcessorNoiseProfile<Scalar, StatsAccumulatorHistogram<Scalar>, STATS_VARIANCE> inp1(_infn1, {0});
    ProcessorNoiseProfile<Scalar, StatsAccumulatorHistogram<Scalar>, STATS_VARIANCE> inp2(_infn2, {0}, inp1.environment(), {}, inp1.srcDescs());

    // find emitter from symmetric difference
    std::set<int64_t> emits1{inp1.emitters().begin(), inp1.emitters().end()};
    std::set<int64_t> emits2{inp2.emitters().begin(), inp2.emitters().end()};

    std::cerr << "Emitters in " << _infn1 << ":";
    for (auto emitter : emits1)
    {
      std::cerr << " " << emitter;
    }
    std::cerr << std::endl;
    std::cerr << "Emitters in " << _infn2 << ":";
    for (auto emitter : emits2)
    {
      std::cerr << " " << emitter;
    }
    std::cerr << std::endl;

    std::vector<int64_t> emitterVec(emits1.size()+emits2.size());
    auto tail = std::set_symmetric_difference(emits1.begin(), emits1.end(), emits2.begin(), emits2.end(), emitterVec.begin());
    emitterVec.resize(tail - emitterVec.begin());

    if (emitterVec.size() != 1)
    {
      std::cerr << "Profiles must differ by exactly 1 active emitter." << std::endl;
      Terminate(-1);
    }

    auto emitter = emitterVec[0];

    decltype(inp1) * procEmitting, * procSilent;
    std::string fnEmitting, fnSilent;
    if (emits1.count(emitter))
    {
      procEmitting = &inp1;
      procSilent = &inp2;
      fnEmitting = _infn1;
      fnSilent = _infn2;
    }
    else
    {
      procEmitting = &inp2;
      procSilent = &inp1;
      fnEmitting = _infn2;
      fnSilent = _infn1;
    }
    std::cerr << "Producing a profile of emitter " << emitter << " by calculating " << fnEmitting << " - " << fnSilent << std::endl;

    ProcessorNoiseProfile<Scalar, StatsAccumulatorHistogram<Scalar>, STATS_VARIANCE> outp(_outfn, {0}, inp1.environment(), {emitter}, {"1-2-gen-diff-hist",fnSilent,fnEmitting,inp1.srcDescs()[0],inp1.srcDescs()[1]});

    // - [ ] check this to make sure value is actually propagated ..
    //   (for some reason default constructor of accumulator is called)
    outp.fromStatsPredicate(decltype(outp)::StatsAccumulator::fromRemovedFromSum<decltype(outp)::StatsAccumulator::options()>, *procEmitting, *procSilent);

    outp.write();

    // - [X] check for difference functions present in processornoiseprofile
    //    - [X] likely add a difference func to processornoiseprofile
    // - [ ] verify correctness of withRemovedFromSum
    // - [ ] chart result
    // - [X] write result to file
    // - [ ] fix compile-time errors
    // - [ ] fix run-time errors
    // - [ ] test behavior
  }

private:
  std::string _infn1;
  std::string _infn2;
  std::string _outfn;
};

int main(int argc, char ** argv)
{
  T12GenDiff app(&argc, argv);
  app.GetOptions(&argc, argv);
  app.Run(false);
}
