#pragma once

#include "Source.hpp"
#include "Math.hpp"

class Spectrum : public Destination
{
public:
	Spectrum(Source & source);
	~Spectrum();
	void receiveQuadrature(cvec const & data, double samplingHertz, double tunedHertz, double dBGain, double unixSecondsCompleted, class Source & source);

private:
	Source & source;
	std::unique_ptr<class GUIWindow> window;
  FFT fft;
};
