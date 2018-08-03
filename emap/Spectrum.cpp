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

void Spectrum::receiveQuadrature(cvec const & data, double samplingHertz, double tunedHz, double dBGain, double unixSecondsCompleted, class Source & source)
{
	if (&source != &this->source)
		return;

	vec fft = abs(fft.execute(data)) * (128.0 / data.size());

	concat(fft, tail(fft, fft.size()/2));

	window->setLines(fft);
}
