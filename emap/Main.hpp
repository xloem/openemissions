#pragma once

#include <future>

class Main
{
public:
	Main();
	~Main();
	
	void stop();
	static Main & instance();

private:
	static Main * _instance;

	std::promise<void> stopped;
};
