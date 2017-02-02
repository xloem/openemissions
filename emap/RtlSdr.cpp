#include "RtlSdr.hpp"

#include <rtl-sdr.h>
#include <system_error>

class RtlSdrType : public SourceType
{
public:
	RtlSdrType(std::string name, uint32_t index)
	: SourceType(name), index(index)
	{ }
	std::unique_ptr<Source> construct() const
	{
		return std::unique_ptr<Source>(new RtlSdr(index));
	}

private:
	uint32_t index;
};

static class RtlSdrEnumerator
{
public:
	RtlSdrEnumerator()
	{
		uint32_t numDevices = rtlsdr_get_device_count();
		std::string name;
		char serial[256];
		for (uint32_t dev = 0; dev < numDevices; ++ dev) {
			name = rtlsdr_get_device_name(dev);
			rtlsdr_get_device_usb_strings(dev, 0, 0, serial);
			types.emplace_back(name + " " + serial + "(" + std::to_string(dev) + ")", dev);
		}
		for (auto & type : types) {
			Source::_register(&type);
		}
	}
private:
	std::vector<RtlSdrType> types;
} rtlSdrEnumerator;

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

RtlSdr::RtlSdr(uint32_t index)
: running(true)
{
	rtlErrNonzero( rtlsdr_open(&dev, index) );

	// sample rate
	setSampleHertz(2400000);

	// manual gain mode
	rtlErrNonzero( rtlsdr_set_tuner_gain_mode(dev, 1) );
	rtlErrNonzero( rtlsdr_set_agc_mode(dev, 0) );
	
	gains.resize( rtlErrNonpositive( rtlsdr_get_tuner_gains(dev, 0) ) );
	rtlErrNonpositive( rtlsdr_get_tuner_gains(dev, &gains[0]) );

	// Default tuning
	auto hzRange = hertzRange();
	tuneHertz((hzRange.first + hzRange.second) / 2);

	// Default gain
	setDB(dBRange().second);

	// Seems to need this
	rtlsdr_reset_buffer(dev);

	// start
	th = std::thread(&RtlSdr::run, this);
}

RtlSdr::~RtlSdr()
{
	{
		std::lock_guard<std::mutex> devLock(mtx);
		running = false;
	}
	th.join();
	
	rtlErrNonzero( rtlsdr_close(dev) );
	dev = 0;
}

void RtlSdr::run()
{
	int n_read;
	std::vector<uint8_t> buf(16 * 32 * 512);
	itpp::cvec data;

	for (;;) {
		{
			std::lock_guard<std::mutex> devLock(mtx);
			if (!running) break;
			rtlErrNonzero( rtlsdr_read_sync(dev, &buf[0], buf.size(), &n_read) );
		}

		data.set_size(n_read / 2);
		for (int i = 0; i < data.size(); ++ i) {
			data[i].real( (buf[i*2] - 127.5) / 127.5 );
			data[i].imag( (buf[i*2+1] - 127.5) / 127.5 );
		}
		dispatch(data, data.size() / sampleHertz(), hertz(), dB(), 0);
	}
}

double RtlSdr::tuneHertz(double freq)
{
	{
		std::lock_guard<std::mutex> devLock(mtx);
		rtlErrNonzero( rtlsdr_set_center_freq(dev, freq) );
	}
	return hertz();
}

double RtlSdr::hertz()
{
	std::lock_guard<std::mutex> devLock(mtx);
	return rtlsdr_get_center_freq(dev);
}

std::pair<double,double> RtlSdr::hertzRange()
{
	return std::make_pair<double,double>(24000000, 1850000000);
}

double RtlSdr::setDB(double dB)
{
	int closest = gains[0];
	double dist = 1.0/0.0;
	for (int & i : gains) {
		double gain = i / 10.0;
		if (std::abs(gain - dB) < dist) {
			closest = i;
			dist = std::abs(gain - dB);
		}
	}
	{
		std::lock_guard<std::mutex> devLock(mtx);
		rtlErrNonzero(rtlsdr_set_tuner_gain(dev, closest));
	}
	return this->dB();
}

double RtlSdr::dB()
{
	std::lock_guard<std::mutex> devLock(mtx);
	return rtlErrZero(rtlsdr_get_tuner_gain(dev)) / 10.0;
}

std::pair<double,double> RtlSdr::dBRange()
{
	std::pair<double,double> range = std::make_pair<double,double>(1.0/0.0, -1.0/0.0);
	for (int & i : gains) {
		double gain = i / 10.0;
		if (gain < range.first)
			range.first = gain;
		if (gain > range.second)
			range.second = gain;
	}
	return range;
}

double RtlSdr::setSampleHertz(double rate)
{
	{
		std::lock_guard<std::mutex> devLock(mtx);
		rtlErrNonzero( rtlsdr_set_sample_rate(dev, rate) );
	}
	return sampleHertz();
}

double RtlSdr::sampleHertz()
{
	std::lock_guard<std::mutex> devLock(mtx);
	return rtlErrZero( rtlsdr_get_sample_rate(dev) );
}

std::pair<double, double> RtlSdr::sampleHertzRange()
{
	return std::make_pair<double, double>(225001, 3200000);
}
