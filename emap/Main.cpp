#include "Main.hpp"

#include "RtlSdr.hpp"
#include "SDL.hpp"

#include <iostream>
#include <stdexcept>

Main * Main::_instance = 0;

Main::Main()
{
	if (_instance)
		throw std::invalid_argument("Duplicate Main instance");
	_instance = this;

	auto stoppedFuture = stopped.get_future();

	SDLWindow window;

	window.draw(itpp::mat("0 0.5; 0.25 1"));
	window.draw(itpp::vec("0 1 0"));

	std::vector<RtlSdr> radios;
	radios = RtlSdr::createForAllDevices();
	std::cout << radios.size() << " radio(s) created." << std::endl;

	stoppedFuture.wait();
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
