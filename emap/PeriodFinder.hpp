#pragma once

#include "GUI.hpp"
#include "Source.hpp"

#include <itpp/base/mat.h>
#include <itpp/base/vec.h>

class PeriodFinder : public Destination
{
public:
	PeriodFinder(double minFrequency, double maxFrequency, double precisionFrequency, Source & source, double tunedHertz);
	~PeriodFinder();
	void receiveQuadrature(itpp::cvec const & data, double samplingHertz, double tunedHertz, double dBGain, double unixSecondsCompleted, class Source & source);

private:
	void init(double rate);

	std::unique_ptr<GUIWindow> window;

	std::pair<double,double> rangeHz;
  double precisionHz;
	Source & source;
	double tunedHz;

	double expectedSamplingHz;
	double longestPeriod, shortestPeriod;

	itpp::cvec data;
  typename itpp::cvec::value_type dataSum;

	std::vector<itpp::cvec> values;
	std::vector<itpp::vec> weights;
	std::vector<int> offsets;
};
