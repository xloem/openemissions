#include "Source.hpp"

#include <map>

std::vector<Destination *> receivers;
std::map<std::string, SourceType *> sources;

SourceType::SourceType(std::string name)
: name(name)
{ }

void Source::dispatch(itpp::cvec const & data, double secondsDuration)
{
	for (Destination * receiver : receivers)
		try {
			receiver->receive(data, secondsDuration, *this);
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

void Source::_register(SourceType * source)
{
	::sources[source->name] = source;
}

void Source::_deregister(std::string name)
{
	::sources.erase(name);
}
