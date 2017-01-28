#include "Source.hpp"

std::vector<Destination *> Source::receivers;

void Source::dispatch(itpp::cvec const & data)
{
	for (Destination * receiver : receivers)
		try {
			receiver->receive(data, *this);
		} catch (std::invalid_argument) { }
}
