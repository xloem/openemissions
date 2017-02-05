#pragma once

#include "Oscilloscope.hpp"

class Waterfall : public Oscilloscope
{
public:
	Waterfall(Source & source);
	void receiveQuadrature(itpp::cvec const & data, double samplingHz, double tunedHertz, double gainDB, double unixSecondsCompleted, Source & source);
};
