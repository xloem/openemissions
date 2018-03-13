#pragma once

#include "Source.hpp"

#include <itpp/base/vec.h>

class Spectrum : public Destination
{
public:
	Spectrum(Source & source);
	~Spectrum();
	void receiveQuadrature(itpp::cvec const & data, double samplingHertz, double tunedHertz, double dBGain, double unixSecondsCompleted, class Source & source);

private:
	Source & source;
	std::unique_ptr<class GUIWindow> window;
};
