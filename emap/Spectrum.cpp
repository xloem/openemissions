#include "Spectrum.hpp"

#include "GUI.hpp"

Spectrum::Spectrum(Source & source)
: source(source), window(window)
{
	ready();
}

Spectrum::~Spectrum()
{
	done();
}

void Spectrum::receiveQuadrature(itpp::cvec const & data, double samplingHertz, double tunedHz, double dBGain, double unixSecondsCompleted, class Source & source)
{
	if (&source != &this->source)
		return;

	itpp::vec fft = itpp::abs(itpp::fft(data)) * (128.0 / data.size());

	fft = itpp::concat(fft,fft.split(fft.size()/2));

	window->setLines(fft);
}
