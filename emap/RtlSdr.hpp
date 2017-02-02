#pragma once

#include "Source.hpp"

#include <mutex>
#include <thread>
#include <vector>

#include <itpp/base/vec.h>

class RtlSdr : public Source
{
friend class RtlSdrType;
public:
	~RtlSdr();

	double tuneHertz(double);
	double hertz();
	std::pair<double,double> hertzRange();

	double setDB(double);
	double dB();
	std::pair<double,double> dBRange();

	double setSampleHertz(double);
	double sampleHertz();
	std::pair<double,double> sampleHertzRange();

private:
	RtlSdr(uint32_t index);

	bool running;
	std::mutex mtx;
	struct rtlsdr_dev *dev;
	std::vector<int> gains;

	std::thread th;
	void run();
};
