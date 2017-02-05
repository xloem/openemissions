#include "Oscilloscope.hpp"

Oscilloscope::Oscilloscope(Source & source)
: source(source), window(GUIWindow::create())
{
	ready();
}

Oscilloscope::~Oscilloscope()
{
	done();
}

void Oscilloscope::receiveQuadrature(itpp::cvec const & data, double samplingHz, double tunedHertz, double gainDB, double unixSecondsCompleted, Source & source)
{
	if (&source != &this->source)
		return;

	static std::vector<itpp::vec> drawVecs(2);
	drawVecs[0].set_size(std::min(1024,data.size()));
	drawVecs[1].set_size(std::min(1024,data.size()));
	for (int i = 0; i < drawVecs[0].size(); ++ i) {
		drawVecs[0][i] = data[i].real();
		drawVecs[1][i] = data[i].imag();
	}
	window->setLines(drawVecs);
	window->setText(std::to_string(tunedHertz) + " Hz");
}

