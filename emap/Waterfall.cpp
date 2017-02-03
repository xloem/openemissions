#include "Waterfall.hpp"

#include <itpp/signal/transforms.h>

Waterfall::Waterfall(Source & source)
: Oscilloscope(source)
{ }

void Waterfall::receive(itpp::cvec const & datavec, double secondsDuration, double tunedHertz, double gainDB, double unixSecondsCompleted, Source & source)
{
	if (&source != &this->source)
		return;

	itpp::cvec fourier = itpp::fft(datavec);
	itpp::vec mag = itpp::abs(fourier) / 64.0;

	int size = 8192;
	if (mag.size() < size)
		size = mag.size();

	if (data.cols() == 0)
		data.set_size(1, size);
	else
		data.set_size(data.rows() + 1, size, true);

	data.set_row(data.rows() - 1, mag.get(mag.size() / 2 - size / 2, mag.size() / 2 - size / 2 + size - 1));
	//data.set_row(data.rows() - 1, mag.get(mag.size() - size, mag.size() - 1));
	//data.set_row(data.rows() - 1, mag.get(0, size - 1));

	window->draw(data);
}
