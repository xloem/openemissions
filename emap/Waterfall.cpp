#include "Waterfall.hpp"

#include <itpp/signal/transforms.h>

Waterfall::Waterfall(Source & source)
: Oscilloscope(source)
{ }

#include <iostream>

void Waterfall::receive(itpp::cvec const & datavec, double secondsDuration, double tunedHertz, double gainDB, double unixSecondsCompleted, Source & source)
{
	if (&source != &this->source)
		return;

	itpp::cvec fourier = itpp::fft(datavec);
	itpp::vec mag = itpp::abs(fourier) / 64.0;

	window->addRow(mag);
}
