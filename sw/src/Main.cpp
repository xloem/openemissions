#include "RtlSdr.hpp"
#include <iostream>

class Main
{
public:
	Main() {
		radios = RtlSdr::createForAllDevices();
		std::cout << radios.size() << " radio(s) created." << std::endl;
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
