#pragma once

#include "Oscilloscope.hpp"

#include "Math.hpp"

class Waterfall : public Oscilloscope
{
public:
	Waterfall(Source & source);
	void receiveQuadrature(cvec const & data, double samplingHz, double tunedHertz, double gainDB, double unixSecondsCompleted, Source & source);

private:
  FFT fft;
};
