#pragma once

#include "GUI.hpp"
#include "Source.hpp"

class Oscilloscope : public Destination
{
public:
	Oscilloscope(Source & source);
	~Oscilloscope();
	void receive(itpp::cvec const & data, double secondsDuration, double tunedHertz, double gainDB, double unixSecondsCompleted, Source & source);

protected:
	Source & source;
	std::unique_ptr<GUIWindow> window;
};
