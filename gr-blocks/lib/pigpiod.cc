#include "pigpiod.h"

namespace gr {
namespace openemissions {
namespace detail {

pigpiod::sptr pigpiod::get(std::string const & address)
{
  gr::thread::scoped_lock lk(s_mtx);
  if (s_servers.count(address)) {
    return s_servers[address].lock();
  } else {
    sptr result(new pigpiod(address));
    s_servers[address] = result;
    return result;
  }
}

pigpiod::~pigpiod()
{
  gr::thread::scoped_lock lk(s_mtx);

  s_servers.erase(d_address);

  pigpio_stop(d_handle);
}

pigpiod::pigpiod(std::string const & address)
: d_address(address),
  d_handle(start(address))
{
  pigpiothrow(wave_clear(d_handle));
}

int pigpiod::start(std::string const & address)
{
  int portIdx = address.find(':');
  std::string addr, port;
  if (portIdx >= 0) {
    port = address.substr(portIdx + 1);
  } else {
    port = "8888";
  }
  if (portIdx > 0) {
    addr = address.substr(0, portIdx);
  } else {
    addr = "127.0.0.1";
  }
  return pigpiothrow(pigpio_start(addr.c_str(), port.c_str()));
}

int pigpiothrow(int result_code)
{
  if (result_code < 0) {
    throw std::runtime_error(pigpio_error(result_code)); 
  }
  switch (result_code) {
  case PI_WAVE_NOT_FOUND:
    throw std::runtime_error("Transmited wave not found");
  case PI_NO_TX_WAVE:
    throw std::runtime_error("No wave being transmitted");
  default:
    return result_code;
  }
}

boost::unordered_map<std::string, boost::weak_ptr<pigpiod>> pigpiod::s_servers;
gr::thread::mutex pigpiod::s_mtx;

}
}
}
