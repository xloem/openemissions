#include <gnuradio/thread/thread.h>

#include <pigpiod_if2.h>

#include <boost/unordered_map.hpp>

namespace gr {
namespace openemissions {
namespace detail {

int pigpiothrow(int result_code);

class pigpiod
{
public:
  typedef boost::shared_ptr<pigpiod> sptr;

  static sptr get(std::string const & address);

  ~pigpiod();

  std::string const d_address;

  const int d_handle;

  struct wave_sender {
    unsigned d_pin;
    unsigned d_accumulated_us;
    std::vector<gpioPulse_t> d_pulses;
  };

  boost::unordered_map<unsigned, wave_sender *> d_wave_senders;
  std::list<int> d_waveforms_in_flight;

  gr::thread::mutex d_mtx;
  gr::thread::condition_variable d_cond;

private:
  /*needs_lock*/ pigpiod(std::string const & address);

  static int start(std::string const & address);

  static boost::unordered_map<std::string, boost::weak_ptr<pigpiod>> s_servers;
  static gr::thread::mutex s_mtx;
};

}
}
}
