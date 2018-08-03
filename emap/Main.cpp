#include "Main.hpp"

#include "Source.hpp"
#include "Oscilloscope.hpp"
#include "Waterfall.hpp"
#include "PeriodFinder.hpp"
#include "PeriodViewer.hpp"

//#include "FLAC_XMP.hpp"

#include <iostream>
#include <stdexcept>

Main * Main::_instance = 0;

Main::Main()
{
	if (_instance)
		throw std::invalid_argument("Duplicate Main instance");
	_instance = this;

	auto stoppedFuture = stopped.get_future();

	std::vector<std::unique_ptr<Source>> sources;
	for (SourceType const * source : Source::sources())
	{
		sources.emplace_back(source->construct());
	}
	std::cout << sources.size() << " source(s) created." << std::endl;

	//FLAC_XMP flacXmp("test.flac", "test.xmp");

	if (sources.size() > 0) {
		//sources[0]->tuneHertz(100000000);
		//sources[0]->tuneHertz(76750000); // my thinkpad's pixel clock
		sources[0]->tuneHertz(147800000); // my purism's advertised pixel clock
		
		Oscilloscope scope0(*sources[0].get());
		Waterfall waterfall0(*sources[0].get());
		//PeriodFinder periods0(39.09, 40.01, 0.01, *sources[0].get(), sources[0]->hertz());
		//PeriodFinder periods0(55, 65, 0.01, *sources[0].get(), sources[0]->hertz());
		//PeriodFinder periods0(300000, 1200000, *sources[0].get(), sources[0]->hertz());
		//PeriodViewer period0(60.1, *sources[0].get());

		stoppedFuture.wait();
	}
}

Main::~Main()
{
	_instance = 0;
}

void Main::stop()
{
	std::move(stopped).set_value();
}

Main & Main::instance()
{
	if (!_instance)
		throw std::invalid_argument("Not instantiated");
	return *_instance;
}

int main(int argc, char const * const * argv)
{
	Main();
}
