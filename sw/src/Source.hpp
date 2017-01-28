#pragma once

#include <itpp/base/vec.h>
#include <stdexcept>

class Destination
{
friend class Source;
private:
	virtual void receive(itpp::cvec const & data, class Source & source) { throw std::invalid_argument("unimplemented"); }
};

class Source
{
public:
	static void subscribe(Destination & destination);

	virtual double tuneHertz(double) { throw std::invalid_argument("unimplemented"); }
	virtual double hertz() { throw std::invalid_argument("unimplemented"); }
	virtual std::pair<double,double> hertzRange() { throw std::invalid_argument("unimplemented"); }

	virtual double setDb(double) { throw std::invalid_argument("unimplemented"); }
	virtual double dB() { throw std::invalid_argument("unimplemented"); }
	virtual std::pair<double,double> dBRange() { throw std::invalid_argument("unimplemented"); }

protected:
	void dispatch(itpp::cvec const & data);

	static std::vector<Destination *> receivers;
};
