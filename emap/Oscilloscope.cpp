#include "Oscilloscope.hpp"

Oscilloscope::Oscilloscope(Source & source, std::string title)
: source(source), window(GUIWindow::create(title))
{
	ready();
}

Oscilloscope::~Oscilloscope()
{
	done();
}

void Oscilloscope::receiveQuadrature(cvec const & data, double samplingHz, double tunedHertz, double gainDB, double unixSecondsCompleted, Source & source)
{
	if (&source != &this->source)
		return;

	static std::vector<vec> drawVecs(2);
	drawVecs[0].resize(std::min(1024L,data.size()));
	drawVecs[1].resize(std::min(1024L,data.size()));
	for (int i = 0; i < drawVecs[0].size(); ++ i) {
		drawVecs[0][i] = data[i].real();
		drawVecs[1][i] = data[i].imag();
	}
	window->setLines(drawVecs);
	window->setText(std::to_string(tunedHertz) + " Hz");
}

