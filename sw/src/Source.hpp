#pragma once

#include <itpp/base/vec.h>
#include <stdexcept>

class Destination
{
friend class Source;
private:
	virtual void receive(cvec const & data, Source & source) = 0;
}

class Source
{
public:
	static void subscribe(Destination destination);

	virtual double tuneHertz(double) { throw std::invalid_argument("unimplemented"); }
	virtual double hertz() { throw std::invalid_argument("unimplemented"); }
	virtual std::pair<double,double> hertzRange() { throw std::invalid_argument("unimplemented"); }

	virtual double setDb(double) { throw std::invalid_argument("unimplemented"); }
	virtual double dB() { throw std::invalid_argument("unimplemented"); }
	virtual std::pair<double,double> dBRange() { throw std::invalid_argument("unimplemented"); }

protected:
	void dispatch(cvec const & data);

	// TODO static list of receivers
}

void Source::dispatch(cvec const & data)
{
	// TODO dispatch to receiver list
}
