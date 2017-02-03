#pragma once

#include "Oscilloscope.hpp"

class Waterfall : public Oscilloscope
{
public:
	Waterfall(Source & source);
	void receive(itpp::cvec const & data, double secondsDuration, double tunedHertz, double gainDB, double unixSecondsCompleted, Source & source);

private:
	itpp::mat data;
};
