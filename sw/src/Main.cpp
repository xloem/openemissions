#include "RtlSdr.hpp"
#include "SDL.hpp"
#include <iostream>

class Main
{
public:
	Main() {
		SDLWindow window;

		window.draw({0,1,0});

		radios = RtlSdr::createForAllDevices();
		std::cout << radios.size() << " radio(s) created." << std::endl;

		SDLEvents::loop();
	}
	~Main() {
	}
private:
	std::vector<RtlSdr> radios;
};

int main(int argc, char const * const * argv)
{
	Main();
}
