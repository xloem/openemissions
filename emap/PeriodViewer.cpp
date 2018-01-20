#include "PeriodViewer.hpp"

#include <cmath>
#include <itpp/base/operators.h>
#include <itpp/signal/transforms.h>

#include "GUI.hpp"

PeriodViewer::PeriodViewer(double frequency, Source & source)
: window(GUIWindow::create()), source(source), frequency(frequency), expectedSamplingHz(-1)
{
	ready();
}

PeriodViewer::~PeriodViewer()
{
	done();
}

double PeriodViewer::tuneHertz(double frequency)
{
	std::lock_guard<std::mutex> lk(mtx);
	this->frequency = frequency;
	this->expectedSamplingHz = -1;
	return frequency;
}

double PeriodViewer::hertz()
{
	std::lock_guard<std::mutex> lk(mtx);
	return frequency;
}

void PeriodViewer::receiveQuadrature(itpp::cvec const & data, double samplingHertz, double tunedHertz, double dBGain, double unixSecondsCompleted, class Source & source)
{
	if (&source != &this->source)
		return;

	{
		std::lock_guard<std::mutex> lk(mtx);
		if (samplingHertz != expectedSamplingHz || tunedHertz != expectedTuningHz) {
			expectedSamplingHz = samplingHertz;
			expectedTuningHz = tunedHertz;
			periodLen = samplingHertz / frequency;
			currentOffset = 0;
			waveform.set_size(periodLen);
			waveform.zeros();
			weight.set_size(periodLen);
			weight.zeros();
		}
	}

	for (double i = 0, wavePos = currentOffset; i < data.size(); ++ i) {
		int intWavePos = wavePos;
		waveform[intWavePos] += data[i];
		weight[intWavePos] = weight[intWavePos] + 1.0;
		++ wavePos;
		if (wavePos >= waveform.size())
			wavePos -= waveform.size();
	}

	window->setLines({elem_div(itpp::abs(waveform), weight)});

	/*
	add vector to other vector with fractional offset
	I'm trying to sum the periods given that the actual period length is not an integral
	multiple of samples.
	This shouldn't be hard.
	B ut it is, right now.
	currentOffset is a double that represents somehow how where we are in the SAMPLES
	is offset in the ACTUAL WAVE.
	So it loops.
	Right now currentOffset is 0, representing that sample 0 and wave-sample 0 precisely
	coincide.  We are beginning recording.
	If periodLen is say pi, then after we hit the 3rd sample, we have a fractional sample
	that could be skipped, or could be placed as the next 1st sample.
	So, after we hit the 3rd sample, we could imagine ourselves as recording the 1st
	sample of the next wave, which is at pi.
	This pi sample maps to sample #4.
	the actual offset is (pi - 3).
	so when we hit the 4th sample, we could store (pi - 3) as the currentOffset.
	This means that the value in the data is equal to the position in the wave
	well
	the situation is that the value in the data belongs at wave position ...
	confused
	our things are shifted.
	the wave start is now 0.14159 more than the sample start
	in theory the sample is precisely sampled at a given moment
	so its value belongs, roughly, at some spot in the wave
	and if we snap those spots to the wave spots, the phase problems this causes will go away on average, because some will be high and some will be low
	where does the value at 4 belong?
	well, the "wave address" of 4 is 4 - pi, or roughly 0.85
	so if I want to find the "wave address" from the sample number, I can add (4 - pi)
	then I can round the wave address to yeeeeeeeeesh obviously ?? how confused my brain is
	not so obvious now!
	all I had to do was start implementing it
*/
}
