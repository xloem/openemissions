#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <freesdr/Block.hpp>

#include <freesdr/LinAlgImplementationArmadillo.hpp>
#include <freesdr/StreamGraphImplementationPothos.hpp>

// A first test block: to produce a ROOT waterfall from a Pothos graph input in a modular manner

class Waterfall : public freesdr::Block<freesdr::LinAlgImplementationArmadillo, freesdr::StreamGraphImplementationPothos, freesdr::InputStreamType<double>>
{
public:
  static std::string blockCategory() { return "test"; }
  static std::sting blockName() { return "waterfall"; }

  virtual void work(freesdr::Block::InputBuffer<double> & inputBuf)
  {
  };
};

BlockRegistration<Waterfall> waterfallReg;

TEST_CASE( "Armadillo/Pothos input block", "[blockarmpothinput]" )
{
  // TODO: make these calls modular
  Pothos::ScopedInit pothos;
  Pothos::Topology topology;
  auto waterfall = Pothos::BlockRegistry::make("/test/waterfall");
  auto source = Pothos::BlockRegistry::make("/blocks/vector_source", Pothos::DType(typeid(double)));
  source.call("setElements", std::vector<double>({-1,0,1}));
  source.call("setMode", "REPEAT");
  topology.connect(source, 0, fft, 0);
  topology.commit();
  topology.waitInactive();
  REQUIRE( 1 == 1 );
}
