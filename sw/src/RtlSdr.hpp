#pragma once

#include "Source.hpp"

#include <vector>
#include <itpp/base/vec.h>

class RtlSdr : public Source
{
public:
	static std::vector<RtlSdr> createForAllDevices();
	
	RtlSdr(uint32_t index);

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
	struct rtlsdr_dev *dev;
	std::vector<int> gains;

	itpp::cvec data;

	static void readAsyncCb(unsigned char * buf, uint32_t len, RtlSdr * ctx);
}


extern "C" {
#include <rtl-sdr.h>
}
#include <system_error>

std::vector<RtlSdr> RtlSdr::createForAllDevices()
{
	uint32_t count = rtlsdr_get_device_count();
	std::vector<RtlSdr> ret;
	for (uint32_t index = 0; index < count; ++ index)
	{
		ret.emplace_back(index);
	}
	return ret;
}

static void rtlErrAny(int r)
{
	throw std::system_error(r, std::generic_category());
}

static int rtlErrNonzero(int r)
{
	if (r)
		rtlErrAny(r);
	return r;
}

static int rtlErrZero(int r)
{
	if (!r)
		rtlErrAny(r);
	return r;
}

static int rtlErrNonpositive(int r)
{
	if (r <= 0)
		rtlErrAny(r);
	return r;
}

static int rtlErrNegative(int r)
{
	if (r < 0)
		rtlErrAny(r);
	return r;
}

RtlSdr::RtlSdr(uint32_t index)
{
	rtlErrNonzero( rtlsdr_open(&dev, index) );

	// sample rate
	setSampleHertz(2400000);

	// manual gain mode
	rtlErrNonzero( rtlsdr_set_tuner_gain_mode(dev, 1) );
	rtlErrNonzero( rtlsdr_set_agc_mode(dev, 0) );
	
	gains.resize( rtlErr( rtlsdr_get_tuner_gains(dev, 0) );
	rtlErrNonpositive( rtlsdr_get_tuner_gains(dev, &gains[0]) );

	// start
	rtlErrNonzero( rtlsdr_read_async(dev, rtlReadSyncCb, this, 0, 0) );
}

RtlSdr::~RtlSdr()
{
	rtlErrNonzero( rtlsdr_cancel_async(dev) );
	rtlErrNonzero( rtlsdr_close(dev) );
	dev = 0;
}

void RtlSdr::readAsyncCb(unsigned char *buf, uint32_t len, RtlSdr * ctx)
{
	ctx->data.set_size(len / 2);
	for (size_t i = 0; i < len / 2; ++ i)
		ctx->data[i].real = (buf[i*2] - 127.5d) / 127.5d;
		ctx->data[i].imag = (buf[i*2+1] - 127.5d) / 127.5d;
	ctx->dispatch(ctx->data);
}

double RtlSdr::tuneHertz(double freq)
{
	int directSampling = rtlErrNegative( rtlsdr_get_direct_sampling(dev);
	if (freq > 28800000) {
		if (directSampling)
			rtlErrNonzero( rtlsdr_set_direct_sampling(dev, 0) );
		rtlErrNonzero( rtlsdr_set_center_freq(dev, freq) );
	} else {
		if (!directSampling)
			rtlErrNonzero( rtlsdr_set_direct_sampling(dev, 1) );
		rtlErrNonzero( rtlsdr_set_center_freq(dev, freq) );
	}
	return hertz();
}

double RtlSdr::hertz()
{
	return rtlErrZero(rtlsdr_get_center_freq(dev));
}

std::pair<double,double> RtlSdr::hertzRange()
{
	//return std::make_pair<double,double>(24000000, 1850000000);
	return std::make_pair<double,double>(0, 1850000000);
}

double RtlSdr::setDb(double dB)
{
	int closest = gains[0];
	double dist = 1.0/0.0d;
	for (int & i : gains) {
		double gain = i / 10.0d;
		if (abs(gain - dB) < dist) {
			closest = i;
			dist = abs(gain - dB);
		}
	}
	rtlErrNonzero(rtlsdr_set_tuner_gain(dev, closest));
	return dB();
}

double RtlSdr::dB()
{
	return rtlErrZero(rtlsdr_get_tuner_gain(dev)) / 10.0d;
}

std::pair<double,double> RtlSdr::dBRange()
{
	std::pair<double,double> range = std::make_pair<double,double>(1.0/0.0d, -1.0/0.0d);
	for (int & i : gains) {
		double gain = i / 10.0d;
		if (gain < range.left)
			range.left = gain;
		if (gain > range.right)
			range.right = gain;
	}
}

double RtlSdr::setSampleHertz(double rate)
{
	rtlErrNonzero( rtlsdr_set_sample_rate(dev, rate) );
	return sampleHertz();
}

double RtlSdr::sampleHertz()
{
	return rtlErrZero( rtlsdr_get_sample_rate(dev) );
}

std::pair<double, double> RtlSdr::sampleHertzRange()
{
	return std::make_pair<double, double>(225001, 3200000);
}
