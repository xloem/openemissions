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

#include <future>

enum SerializedCommandIDs
{
	RTLSDR_TUNE,
	RTLSDR_GAIN,
	RTLSDR_RATE,
	RTLSDR_QUIT,
	RTLSDR_NO_COMMAND = 0xff
};

class SerializedCommand
{
public:
	SerializedCommand(int ID, double arg)
	: _ID(ID), _arg(arg)
	{ }

	static void initialize(std::vector<uint8_t> & data)
	{
		data.assign({RTLSDR_NO_COMMAND});
	}

	double send(std::vector<uint8_t> & data, std::mutex & mtx)
	{
		std::promise<double> promise;
		std::future<double> future = promise.get_future();

		_promise = &promise;

		{
			std::lock_guard<std::mutex> lk(mtx);
			int idx = data.size() - 1;
			data.resize(idx + sizeof(SerializedCommand) + 1);
			*reinterpret_cast<SerializedCommand *>(&data[idx]) = *this;
			data[idx + sizeof(SerializedCommand)] = RTLSDR_NO_COMMAND;
		}
		
		return future.get();
	}

	static int get(std::vector<uint8_t> & data, SerializedCommand * & command, size_t last = 0)
	{
		if (data[last] == RTLSDR_NO_COMMAND)
			return 0;
		command = reinterpret_cast<SerializedCommand*>(&data[last]);
		return last + sizeof(SerializedCommand);
	}

	int ID() { return _ID; }
	double arg() { return _arg; }
	std::promise<double> && promise() { return std::move(*_promise); }

private:
	int _ID;
	double _arg;
	std::promise<double> * _promise;

};
RtlSdr::RtlSdr(uint32_t index)
: DriverThread([]{
	constexpr int BUFFER_MS = 64;

	// ~ 2 MHz default sampling rate, 2 bytes per sample
	std::vector<uint8_t> buf(BUFFER_MS * (4 << 10));

	SerializedCommand::initialize(buf);
	buf.resize(buf.capacity());

	return std::move(buf);
})
{
	SerializedCommand::initialize(commands);

	rtlErrNonzero( rtlsdr_open(&dev, index) );

	// manual gain mode
	rtlErrNonzero( rtlsdr_set_tuner_gain_mode(dev, 1) );
	rtlErrNonzero( rtlsdr_set_agc_mode(dev, 0) );
	
	gains.resize( rtlErrNonpositive( rtlsdr_get_tuner_gains(dev, 0) ) );
	rtlErrNonpositive( rtlsdr_get_tuner_gains(dev, &gains[0]) );

	// sample rate
	rtlErrNonzero( rtlsdr_set_sample_rate(dev, 2400000) );

	// Default tuning is mid-range
	auto hzRange = hertzRange();
	rtlErrNonzero( rtlsdr_set_center_freq(dev, (hzRange.first + hzRange.second) / 2) );

	// Default gain is maximum
	rtlErrNonzero( rtlsdr_set_tuner_gain(dev, dBRange().second * 10) );

	// Seems to need this
	rtlsdr_reset_buffer(dev);

	start();

}

RtlSdr::~RtlSdr()
{
	SerializedCommand(RTLSDR_QUIT, 0).send(commands, commandsMtx);

	join();
	
	rtlErrNonzero( rtlsdr_close(dev) );
	dev = 0;
}

bool RtlSdr::fill(std::vector<uint8_t> & buf)
{
	SerializedCommand * command;
	for (size_t idx = 0; (idx = SerializedCommand::get(buf, command, idx));) {
		try {
			double result = 0;
			switch (command->ID()) {
			case RTLSDR_TUNE:
				rtlErrNonzero( rtlsdr_set_center_freq(dev, command->arg()) );
				break;
			case RTLSDR_GAIN: 
				rtlErrNonzero( rtlsdr_set_tuner_gain(dev, command->arg()) );
				break;
			case RTLSDR_RATE:
				rtlErrNonzero( rtlsdr_set_sample_rate(dev, command->arg()) );
				break;
			case RTLSDR_QUIT:
				command->promise().set_value(result);
				buf.clear();
				return false;
			default:
				throw std::logic_error("corrupt commands buffer");
			}
			command->promise().set_value(result);
		} catch(std::logic_error) {
			throw;
		} catch(...) {
			command->promise().set_exception(std::current_exception());
		}
	}

	int n_read;
	rtlErrNonzero( rtlsdr_read_sync(dev, &buf[0], buf.size(), &n_read) );
	buf.resize(n_read);

	return true;
}

bool RtlSdr::process(std::vector<uint8_t> & buf)
{
	if (buf.empty())
		return false;

	itpp::cvec data;

	data.set_size(buf.size() / 2);
	for (int i = 0; i < data.size(); ++ i) {
		data[i].real( (buf[i*2] - 127.5) / 127.5 );
		data[i].imag( (buf[i*2+1] - 127.5) / 127.5 );
	}
	dispatchQuadrature(data, sampleHertz(), hertz(), dB(), 0);

	{
		std::lock_guard<std::mutex> commandsLk(commandsMtx);
		if (commands.size() > buf.size())
			buf = commands;
		else
			memcpy(&buf[0], &commands[0], commands.size());
		SerializedCommand::initialize(commands);
	}
	return true;
}

double RtlSdr::tuneHertz(double freq)
{
	SerializedCommand(RTLSDR_TUNE, freq).send(commands, commandsMtx);
	return hertz();
}

double RtlSdr::hertz()
{
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
	SerializedCommand(RTLSDR_TUNE, closest).send(commands, commandsMtx);
	return this->dB();
}

double RtlSdr::dB()
{
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
	SerializedCommand(RTLSDR_RATE, rate).send(commands, commandsMtx);
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
