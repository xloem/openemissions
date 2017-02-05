#include "Source.hpp"

#include <map>
#include <mutex>
#include <set>

std::mutex receiversMtx;
std::set<Destination *> receivers;
std::map<std::string, SourceType *> sources;

void Destination::ready()
{
	std::lock_guard<std::mutex> receiversLk(receiversMtx);
	receivers.insert(this);
}

void Destination::done()
{
	std::lock_guard<std::mutex> receiversLk(receiversMtx);
	receivers.erase(this);
}

Destination::~Destination()
{
	std::lock_guard<std::mutex> receiversLk(receiversMtx);
	if (receivers.count(this))
		throw std::logic_error("done() never called");
}

SourceType::SourceType(std::string name)
: name(name)
{ }

void Source::dispatchQuadrature(itpp::cvec const & data, double samplingHertz, double tunedHertz, double dBGain, double unixSecondsCompleted)
{
	std::lock_guard<std::mutex> receiversLk(receiversMtx);
	for (Destination * receiver : receivers)
		try {
			receiver->receiveQuadrature(data, samplingHertz, tunedHertz, dBGain, unixSecondsCompleted, *this);
		} catch (std::invalid_argument) { }
}

std::vector<SourceType const *> Source::sources()
{
	std::vector<SourceType const *> ret;
	for (auto const & it : ::sources) {
		ret.push_back(it.second);
	}
	return ret;
}

Source::~Source()
{ }

void Source::_register(SourceType * source)
{
	::sources[source->name] = source;
}

void Source::_deregister(std::string name)
{
	::sources.erase(name);
}
