#pragma once

#include "Source.hpp"

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

	struct rtlsdr_dev *dev;
	std::vector<int> gains;

	itpp::cvec data;

	static void readAsyncCb(unsigned char * buf, unsigned len, void * ctx);
};
