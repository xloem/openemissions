#pragma once

#include <itpp/base/vec.h>
#include <memory>
#include <stdexcept>
#include <string>

class Destination
{
friend class Source;
private:
	virtual void receive(itpp::cvec const & data, double secondsDuration, class Source & source) { throw std::invalid_argument("unimplemented"); }
};

class SourceType
{
public:
	std::string const name;
	virtual std::unique_ptr<class Source> construct() const = 0;

protected:
	SourceType(std::string name);
};

class Source
{
public:
	static std::vector<SourceType const *> sources();

	static void subscribe(Destination & destination);

	virtual double tuneHertz(double) { throw std::invalid_argument("unimplemented"); }
	virtual double hertz() { throw std::invalid_argument("unimplemented"); }
	virtual std::pair<double,double> hertzRange() { throw std::invalid_argument("unimplemented"); }

	virtual double setDb(double) { throw std::invalid_argument("unimplemented"); }
	virtual double dB() { throw std::invalid_argument("unimplemented"); }
	virtual std::pair<double,double> dBRange() { throw std::invalid_argument("unimplemented"); }

	static void _register(SourceType *  source);
	static void _deregister(std::string name);

protected:
	void dispatch(itpp::cvec const & data, double secondsDuration);
};
