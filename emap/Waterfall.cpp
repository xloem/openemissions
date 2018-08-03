#include "Waterfall.hpp"

Waterfall::Waterfall(Source & source)
: Oscilloscope(source, "Waterfall")
{ }

#include <iostream>

void Waterfall::receiveQuadrature(cvec const & datavec, double samplingHz, double tunedHertz, double gainDB, double unixSecondsCompleted, Source & source)
{
	if (&source != &this->source)
		return;

	vec row = abs(fft.execute(datavec)) * (256.0 / datavec.size());

  size_t half = row.size() / 2;
  row.head(half).swap(row.tail(half));

	window->addRow(row);

	window->setText(std::to_string(tunedHertz) + " Hz");
}
