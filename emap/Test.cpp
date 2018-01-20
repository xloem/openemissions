#include "Main.hpp"

#include "GUI.hpp"

#include <itpp/base/random.h>
#include <itpp/signal/transforms.h>
#include <itpp/stat/misc_stat.h>

#include <iostream>
#include <map>


Main * Main::_instance = 0;

struct SignalGraders
{
	static std::vector<double> ideal_karls(int period, itpp::cvec data)
	{
		itpp::vec partial(period);
		partial *= 0.0;
		int i;
		for (i = 0; i < data.size() / period; ++ i)
			partial += itpp::abs(data(i * period, (i + 1) * period));
		partial /= i;
		return { itpp::sum(itpp::abs(partial - itpp::mean(partial))) };
	}

	static std::vector<double> fourier_harmonics(int period, itpp::cvec data)
	{
		itpp::cvec total = itpp::fft(data) / data.size();
		itpp::vec partial(period);
		partial *= 0.0;
		int i;
		for (i = 0; i < data.size() / period; ++ i)
			partial += itpp::abs(total(i * period, (i + 1) * period));
		partial /= i;
		std::vector<double> ret;
		ret.push_back(partial(0) / itpp::mean(partial(1,partial.size())));
		double weightedIndex = 0;
		double weightedTotal = 0;
		for (i = 0; i < partial.size(); ++ i) {
			weightedIndex += partial(i) * i;
			weightedTotal += partial(i);
		}
		ret.push_back(weightedIndex / weightedTotal / partial.size());
		return ret;
	}
};

Main::Main()
{
	_instance = this;
	auto stoppedFuture = stopped.get_future();

	itpp::Complex_Normal_RNG rand(0, 1);

	auto win1 = GUIWindow::create();

	constexpr int signal_length = 256;//2048;
	constexpr int signal_count = 256;

	itpp::cvec signal = rand(signal_length) * 0.1;
	itpp::cvec data(signal_length * signal_count);
	for (int i = 0; i < signal_count; ++ i)
		data.set_subvector(i * signal_length, signal + rand(signal_length));

	std::map<std::string,std::vector<double>(*)(int, itpp::cvec)> graders = {{"ideal_karl",SignalGraders::ideal_karls}, {"fourier_harmonic",SignalGraders::fourier_harmonics}};
	
	for (auto const & kv : graders) {
		auto good = kv.second(signal_length, data);
		auto bad = kv.second(signal_length * 3.141592653589793238/3, data);
		std::cout << kv.first << ": ";
		for (auto d : ret)
			std::cout << d << ", ";
		std::cout << std::endl;
	}

	stoppedFuture.wait();
}

Main::~Main()
{
	_instance = 0;
}

void Main::stop()
{
	std::move(stopped).set_value();
}

Main & Main::instance()
{
	return *_instance;
}

int main( int argc, char const * const * argv)
{
	Main();
}
