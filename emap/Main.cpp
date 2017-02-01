#include "Main.hpp"

#include "RtlSdr.hpp"
#include "GUI.hpp"

#include <iostream>
#include <stdexcept>

Main * Main::_instance = 0;

Main::Main()
{
	if (_instance)
		throw std::invalid_argument("Duplicate Main instance");
	_instance = this;

	auto stoppedFuture = stopped.get_future();

	auto window = GUIWindow::create();

	window->draw(itpp::mat("0 0.5; 0.25 1"));
	window->draw(itpp::vec("0 1 0"));

	std::vector<std::unique_ptr<Source>> sources;
	for (SourceType const * source : Source::sources())
	{
		sources.emplace_back(source->construct());
	}
	std::cout << sources.size() << " source(s) created." << std::endl;

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
