#pragma once

#include "Source.hpp"

#include "DriverThread.hpp"

#include <mutex>
#include <vector>


class RtlSdr : public Source, public DriverThread<std::vector<uint8_t>>
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

	std::vector<uint8_t> && construct();
	bool fill(std::vector<uint8_t> &);
	bool process(std::vector<uint8_t> &);

	struct rtlsdr_dev *dev;
	std::vector<int> gains;

	std::vector<uint8_t> commands;
	std::mutex commandsMtx;
};
