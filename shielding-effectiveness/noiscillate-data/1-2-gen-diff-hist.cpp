// 1-2-gen-diff-hist.cpp
//
// Given two histogram noise profiles of the same environment, with one emitter different betweeen them,
// produces an estimated histogram noise profile of the emitter itself.

#include <algorithm>

#include <TCanvas.h>

#include "chartprofile.hpp"
#include "processornoiseprofile.hpp"
#include "rootapplication.hpp"

using Scalar = double;

class T12GenDiff : public RootApplication
{
public:
  using Profile = ProcessorNoiseProfile<Scalar, StatsAccumulatorHistogram<Scalar>, STATS_VARIANCE>;

  T12GenDiff(Int_t * argc, char ** argv)
  : RootApplication(
      "1-2-gen-diff-hist",
      "Loads two histogram profiles of the same environment that differ by\n"
      "only one emitter, and produces a profile of that emitter.",
      argc, argv),
    _inchart1("No Emitter"),
    _inchart2("With Emitter"),
    _outchart("Emitter")
  {
    regArg({{"-o"}, {"emitter.noisep"}}, {"output file to write difference in"},
      [this](std::string, std::string outfn)
      {
        _outfn = outfn;
      }, 1 // required
    );
    regArg({{}, {"1.noisep"}}, {"first noise histogram to compare"},
      [this](std::string, std::string fname1)
      {
        _infn1 = fname1;
      }, 1 // required
    );
    regArg({{}, {"2.noisep"}}, {"second noise histogram to compare"},
      [this](std::string, std::string fname2)
      {
        _infn2 = fname2;
      }, 1 // required
    );
  }

  virtual void GetOptions(Int_t * argc, char **argv) override
  {
    RootApplication::GetOptions(argc, argv);
    if (quiet < 1)
    {
      _canvas.reset(new TCanvas());
      _canvas->ToggleToolTips();
      _canvas->cd();
      _canvas->Divide(1,2);
      _inchart1.pad() = _canvas->GetPad(1);
      _inchart2.pad() = _canvas->GetPad(1);
      _outchart.pad() = _canvas->GetPad(2);
    }
  }

  virtual void Run(Bool_t retrn) override
  {
    ProcessorNoiseProfile<Scalar, StatsAccumulatorHistogram<Scalar>, STATS_VARIANCE> inp1(_infn1, {0});
    std::cout << "== Source identification strings ==" << std::endl;
    for (auto str : inp1.srcDescs())
    {
      std::cout << str << std::endl;
    }
    std::cout << "===================================" << std::endl;
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

    _inchart1.prepPoints(*procSilent, quiet > 0);
    _inchart2.prepPoints(*procEmitting, quiet > 0);
    _inchart1.paint(*procSilent);
    _inchart2.paint(*procEmitting);
    _inchart1.finalize();
    _inchart2.finalize();
    _inchart1.paint(*procSilent);
    _inchart2.paint(*procEmitting);
    tick();

    ProcessorNoiseProfile<Scalar, StatsAccumulatorHistogram<Scalar>, STATS_VARIANCE> outp(_outfn, {0}, inp1.environment(), {emitter}, {"1-2-gen-diff-hist",fnSilent,fnEmitting,inp1.srcDescs()[0],inp1.srcDescs()[1]});

    // - [ ] check this to make sure value is actually propagated ..
    //   (for some reason default constructor of accumulator is called)
    outp.fromStatsPredicate(decltype(outp)::StatsAccumulator::fromRemovedFromSum<decltype(outp)::StatsAccumulator::options()>, *procEmitting, *procSilent);

    _outchart.prepPoints(outp, quiet > 0);
    _outchart.paint(outp);
    _outchart.finalize();
    _outchart.paint(outp);

    outp.write();

    // - [X] check for difference functions present in processornoiseprofile
    //    - [X] likely add a difference func to processornoiseprofile
    // - [ ] chart result
    //      plan:
    //        - [ ] split canvas vertically into two pads
    //        - [ ] draw the original profiles in pad 1
    //        - [ ] draw the difference profile in pad 2
    // - [ ] verify correctness of withRemovedFromSum
    // - [X] write result to file
    // - [ ] fix compile-time errors
    // - [ ] fix run-time errors
    // - [ ] test behavior
  }

private:
  std::string _infn1;
  std::string _infn2;
  std::string _outfn;

  std::unique_ptr<TCanvas> _canvas;
  RootChartProfile<Profile> _inchart1;
  RootChartProfile<Profile> _inchart2;
  RootChartProfile<Profile> _outchart;
};

int main(int argc, char ** argv)
{
  T12GenDiff app(&argc, argv);
  app.GetOptions(&argc, argv);
  app.Run(false);
}
