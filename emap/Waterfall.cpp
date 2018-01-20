#include "Waterfall.hpp"

#include <itpp/signal/transforms.h>

Waterfall::Waterfall(Source & source)
: Oscilloscope(source)
{ }

#include <iostream>

void Waterfall::receiveQuadrature(itpp::cvec const & datavec, double samplingHz, double tunedHertz, double gainDB, double unixSecondsCompleted, Source & source)
{
	if (&source != &this->source)
		return;

	itpp::vec row = itpp::abs(itpp::fft(datavec)) * (128.0 / datavec.size());

	row = itpp::concat(row,row.split(row.size()/2));

	window->addRow(row);

	window->setText(std::to_string(tunedHertz) + " Hz");
}
