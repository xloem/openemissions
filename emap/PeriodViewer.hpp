#pragma once

#include "Source.hpp"

#include <mutex>

class PeriodViewer : public Destination
{
public:
	PeriodViewer(double frequency, Source & source);
	~PeriodViewer();
	void receiveQuadrature(cvec const & data, double samplingHertz, double tunedHertz, double dBGain, double unixSecondsCompeted, class Source & source);

	double tuneHertz(double);
	double hertz();

private:
	std::unique_ptr<class GUIWindow> window;
	std::mutex mtx;
	Source & source;
	double frequency;
	
	double expectedSamplingHz;
	double expectedTuningHz;

	double periodLen;
	double currentOffset;
	
	cvec waveform;
	vec weight;
};
