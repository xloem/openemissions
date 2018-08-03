#pragma once

#include "GUI.hpp"
#include "Source.hpp"

class Oscilloscope : public Destination
{
public:
	Oscilloscope(Source & source, std::string title = "Oscilloscope");
	~Oscilloscope();
	void receiveQuadrature(cvec const & data, double samplingHz, double tunedHertz, double gainDB, double unixSecondsCompleted, Source & source);

protected:
	Source & source;
	std::unique_ptr<GUIWindow> window;
};
