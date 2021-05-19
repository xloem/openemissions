#include <gnuradio/thread/thread.h>

#include <pigpiod_if2.h>

#include <boost/unordered_map.hpp>

namespace gr {
  namespace openemissions {
    namespace detail {

      static int pigpiothrow(int result_code)
      {
        switch (result_code) {
        case PI_INIT_FAILED     : throw std::runtime_error("gpioInitialise failed");
        case PI_BAD_USER_GPIO   : throw std::runtime_error("GPIO not 0-31");
        case PI_BAD_GPIO        : throw std::runtime_error("GPIO not 0-53");
        case PI_BAD_MODE        : throw std::runtime_error("mode not 0-7");
        case PI_BAD_LEVEL       : throw std::runtime_error("level not 0-1");
        case PI_BAD_PUD         : throw std::runtime_error("pud not 0-2");
        case PI_BAD_PULSEWIDTH  : throw std::runtime_error("pulsewidth not 0 or 500-2500");
        case PI_BAD_DUTYCYCLE   : throw std::runtime_error("dutycycle outside set range");
        case PI_BAD_TIMER       : throw std::runtime_error("timer not 0-9");
        case PI_BAD_MS          : throw std::runtime_error("ms not 10-60000");
        case PI_BAD_TIMETYPE    : throw std::runtime_error("timetype not 0-1");
        case PI_BAD_SECONDS     : throw std::runtime_error("seconds < 0");
        case PI_BAD_MICROS      : throw std::runtime_error("micros not 0-999999");
        case PI_TIMER_FAILED    : throw std::runtime_error("gpioSetTimerFunc failed");
        case PI_BAD_WDOG_TIMEOUT: throw std::runtime_error("timeout not 0-60000");
        case PI_NO_ALERT_FUNC   : throw std::runtime_error("PI_NO_ALERT_FUNC");
        case PI_BAD_CLK_PERIPH  : throw std::runtime_error("clock peripheral not 0-1");
        case PI_BAD_CLK_SOURCE  : throw std::runtime_error("PI_BAD_CLK_SOURCE");
        case PI_BAD_CLK_MICROS  : throw std::runtime_error("clock micros not 1, 2, 4, 5, 8, or 10");
        case PI_BAD_BUF_MILLIS  : throw std::runtime_error("buf millis not 100-10000");
        case PI_BAD_DUTYRANGE   : throw std::runtime_error("dutycycle range not 25-40000");
        case PI_BAD_SIGNUM      : throw std::runtime_error("signum not 0-63");
        case PI_BAD_PATHNAME    : throw std::runtime_error("can't open pathname");
        case PI_NO_HANDLE       : throw std::runtime_error("no handle available");
        case PI_BAD_HANDLE      : throw std::runtime_error("unknown handle");
        case PI_BAD_IF_FLAGS    : throw std::runtime_error("ifFlags > 4");
        case PI_BAD_CHANNEL     : throw std::runtime_error("DMA channel not 0-15");
        case PI_BAD_SOCKET_PORT : throw std::runtime_error("socket port not 1024-32000");
        case PI_BAD_FIFO_COMMAND: throw std::runtime_error("unrecognized fifo command");
        case PI_BAD_SECO_CHANNEL: throw std::runtime_error("DMA secondary channel not 0-15");
        case PI_NOT_INITIALISED : throw std::runtime_error("function called before gpioInitialise");
        case PI_INITIALISED     : throw std::runtime_error("function called after gpioInitialise");
        case PI_BAD_WAVE_MODE   : throw std::runtime_error("waveform mode not 0-3");
        case PI_BAD_CFG_INTERNAL: throw std::runtime_error("bad parameter in gpioCfgInternals call");
        case PI_BAD_WAVE_BAUD   : throw std::runtime_error("baud rate not 50-250K(RX)/50-1M(TX)");
        case PI_TOO_MANY_PULSES : throw std::runtime_error("waveform has too many pulses");
        case PI_TOO_MANY_CHARS  : throw std::runtime_error("waveform has too many chars");
        case PI_NOT_SERIAL_GPIO : throw std::runtime_error("no bit bang serial read on GPIO");
        case PI_BAD_SERIAL_STRUC: throw std::runtime_error("bad (null) serial structure parameter");
        case PI_BAD_SERIAL_BUF  : throw std::runtime_error("bad (null) serial buf parameter");
        case PI_NOT_PERMITTED   : throw std::runtime_error("GPIO operation not permitted");
        case PI_SOME_PERMITTED  : throw std::runtime_error("one or more GPIO not permitted");
        case PI_BAD_WVSC_COMMND : throw std::runtime_error("bad WVSC subcommand");
        case PI_BAD_WVSM_COMMND : throw std::runtime_error("bad WVSM subcommand");
        case PI_BAD_WVSP_COMMND : throw std::runtime_error("bad WVSP subcommand");
        case PI_BAD_PULSELEN    : throw std::runtime_error("trigger pulse length not 1-100");
        case PI_BAD_SCRIPT      : throw std::runtime_error("invalid script");
        case PI_BAD_SCRIPT_ID   : throw std::runtime_error("unknown script id");
        case PI_BAD_SER_OFFSET  : throw std::runtime_error("add serial data offset > 30 minutes");
        case PI_GPIO_IN_USE     : throw std::runtime_error("GPIO already in use");
        case PI_BAD_SERIAL_COUNT: throw std::runtime_error("must read at least a byte at a time");
        case PI_BAD_PARAM_NUM   : throw std::runtime_error("script parameter id not 0-9");
        case PI_DUP_TAG         : throw std::runtime_error("script has duplicate tag");
        case PI_TOO_MANY_TAGS   : throw std::runtime_error("script has too many tags");
        case PI_BAD_SCRIPT_CMD  : throw std::runtime_error("illegal script command");
        case PI_BAD_VAR_NUM     : throw std::runtime_error("script variable id not 0-149");
        case PI_NO_SCRIPT_ROOM  : throw std::runtime_error("no more room for scripts");
        case PI_NO_MEMORY       : throw std::runtime_error("can't allocate temporary memory");
        case PI_SOCK_READ_FAILED: throw std::runtime_error("socket read failed");
        case PI_SOCK_WRIT_FAILED: throw std::runtime_error("socket write failed");
        case PI_TOO_MANY_PARAM  : throw std::runtime_error("too many script parameters (> 10)");
        case PI_SCRIPT_NOT_READY: throw std::runtime_error("script initialising");
        case PI_BAD_TAG         : throw std::runtime_error("script has unresolved tag");
        case PI_BAD_MICS_DELAY  : throw std::runtime_error("bad MICS delay (too large)");
        case PI_BAD_MILS_DELAY  : throw std::runtime_error("bad MILS delay (too large)");
        case PI_BAD_WAVE_ID     : throw std::runtime_error("non existent wave id");
        case PI_TOO_MANY_CBS    : throw std::runtime_error("No more CBs for waveform");
        case PI_TOO_MANY_OOL    : throw std::runtime_error("No more OOL for waveform");
        case PI_EMPTY_WAVEFORM  : throw std::runtime_error("attempt to create an empty waveform");
        case PI_NO_WAVEFORM_ID  : throw std::runtime_error("no more waveforms");
        case PI_I2C_OPEN_FAILED : throw std::runtime_error("can't open I2C device");
        case PI_SER_OPEN_FAILED : throw std::runtime_error("can't open serial device");
        case PI_SPI_OPEN_FAILED : throw std::runtime_error("can't open SPI device");
        case PI_BAD_I2C_BUS     : throw std::runtime_error("bad I2C bus");
        case PI_BAD_I2C_ADDR    : throw std::runtime_error("bad I2C address");
        case PI_BAD_SPI_CHANNEL : throw std::runtime_error("bad SPI channel");
        case PI_BAD_FLAGS       : throw std::runtime_error("bad i2c/spi/ser open flags");
        case PI_BAD_SPI_SPEED   : throw std::runtime_error("bad SPI speed");
        case PI_BAD_SER_DEVICE  : throw std::runtime_error("bad serial device name");
        case PI_BAD_SER_SPEED   : throw std::runtime_error("bad serial baud rate");
        case PI_BAD_PARAM       : throw std::runtime_error("bad i2c/spi/ser parameter");
        case PI_I2C_WRITE_FAILED: throw std::runtime_error("i2c write failed");
        case PI_I2C_READ_FAILED : throw std::runtime_error("i2c read failed");
        case PI_BAD_SPI_COUNT   : throw std::runtime_error("bad SPI count");
        case PI_SER_WRITE_FAILED: throw std::runtime_error("ser write failed");
        case PI_SER_READ_FAILED : throw std::runtime_error("ser read failed");
        case PI_SER_READ_NO_DATA: throw std::runtime_error("ser read no data available");
        case PI_UNKNOWN_COMMAND : throw std::runtime_error("unknown command");
        case PI_SPI_XFER_FAILED : throw std::runtime_error("spi xfer/read/write failed");
        case PI_BAD_POINTER     : throw std::runtime_error("bad (NULL) pointer");
        case PI_NO_AUX_SPI      : throw std::runtime_error("no auxiliary SPI on Pi A or B");
        case PI_NOT_PWM_GPIO    : throw std::runtime_error("GPIO is not in use for PWM");
        case PI_NOT_SERVO_GPIO  : throw std::runtime_error("GPIO is not in use for servo pulses");
        case PI_NOT_HCLK_GPIO   : throw std::runtime_error("GPIO has no hardware clock");
        case PI_NOT_HPWM_GPIO   : throw std::runtime_error("GPIO has no hardware PWM");
        case PI_BAD_HPWM_FREQ   : throw std::runtime_error("invalid hardware PWM frequency");
        case PI_BAD_HPWM_DUTY   : throw std::runtime_error("hardware PWM dutycycle not 0-1M");
        case PI_BAD_HCLK_FREQ   : throw std::runtime_error("invalid hardware clock frequency");
        case PI_BAD_HCLK_PASS   : throw std::runtime_error("need password to use hardware clock 1");
        case PI_HPWM_ILLEGAL    : throw std::runtime_error("illegal, PWM in use for main clock");
        case PI_BAD_DATABITS    : throw std::runtime_error("serial data bits not 1-32");
        case PI_BAD_STOPBITS    : throw std::runtime_error("serial (half) stop bits not 2-8");
        case PI_MSG_TOOBIG      : throw std::runtime_error("socket/pipe message too big");
        case PI_BAD_MALLOC_MODE : throw std::runtime_error("bad memory allocation mode");
        case PI_TOO_MANY_SEGS   : throw std::runtime_error("too many I2C transaction segments");
        case PI_BAD_I2C_SEG     : throw std::runtime_error("an I2C transaction segment failed");
        case PI_BAD_SMBUS_CMD   : throw std::runtime_error("SMBus command not supported by driver");
        case PI_NOT_I2C_GPIO    : throw std::runtime_error("no bit bang I2C in progress on GPIO");
        case PI_BAD_I2C_WLEN    : throw std::runtime_error("bad I2C write length");
        case PI_BAD_I2C_RLEN    : throw std::runtime_error("bad I2C read length");
        case PI_BAD_I2C_CMD     : throw std::runtime_error("bad I2C command");
        case PI_BAD_I2C_BAUD    : throw std::runtime_error("bad I2C baud rate, not 50-500k");
        case PI_CHAIN_LOOP_CNT  : throw std::runtime_error("bad chain loop count");
        case PI_BAD_CHAIN_LOOP  : throw std::runtime_error("empty chain loop");
        case PI_CHAIN_COUNTER   : throw std::runtime_error("too many chain counters");
        case PI_BAD_CHAIN_CMD   : throw std::runtime_error("bad chain command");
        case PI_BAD_CHAIN_DELAY : throw std::runtime_error("bad chain delay micros");
        case PI_CHAIN_NESTING   : throw std::runtime_error("chain counters nested too deeply");
        case PI_CHAIN_TOO_BIG   : throw std::runtime_error("chain is too long");
        case PI_DEPRECATED      : throw std::runtime_error("deprecated function removed");
        case PI_BAD_SER_INVERT  : throw std::runtime_error("bit bang serial invert not 0 or 1");
        case PI_BAD_EDGE        : throw std::runtime_error("bad ISR edge value, not 0-2");
        case PI_BAD_ISR_INIT    : throw std::runtime_error("bad ISR initialisation");
        case PI_BAD_FOREVER     : throw std::runtime_error("loop forever must be last command");
        case PI_BAD_FILTER      : throw std::runtime_error("bad filter parameter");
        case PI_BAD_PAD         : throw std::runtime_error("bad pad number");
        case PI_BAD_STRENGTH    : throw std::runtime_error("bad pad drive strength");
        case PI_FIL_OPEN_FAILED : throw std::runtime_error("file open failed");
        case PI_BAD_FILE_MODE   : throw std::runtime_error("bad file mode");
        case PI_BAD_FILE_FLAG   : throw std::runtime_error("bad file flag");
        case PI_BAD_FILE_READ   : throw std::runtime_error("bad file read");
        case PI_BAD_FILE_WRITE  : throw std::runtime_error("bad file write");
        case PI_FILE_NOT_ROPEN  : throw std::runtime_error("file not open for read");
        case PI_FILE_NOT_WOPEN  : throw std::runtime_error("file not open for write");
        case PI_BAD_FILE_SEEK   : throw std::runtime_error("bad file seek");
        case PI_NO_FILE_MATCH   : throw std::runtime_error("no files match pattern");
        case PI_NO_FILE_ACCESS  : throw std::runtime_error("no permission to access file");
        case PI_FILE_IS_A_DIR   : throw std::runtime_error("file is a directory");
        case PI_BAD_SHELL_STATUS: throw std::runtime_error(" bad shell return status");
        case PI_BAD_SCRIPT_NAME : throw std::runtime_error("bad script name");
        case PI_BAD_SPI_BAUD    : throw std::runtime_error("bad SPI baud rate, not 50-500k");
        case PI_NOT_SPI_GPIO    : throw std::runtime_error("no bit bang SPI in progress on GPIO");
        case PI_BAD_EVENT_ID    : throw std::runtime_error("bad event id");
        case PI_CMD_INTERRUPTED : throw std::runtime_error("Used by Python");
        case PI_NOT_ON_BCM2711  : throw std::runtime_error("not available on BCM2711");
        case PI_ONLY_ON_BCM2711 : throw std::runtime_error("only available on BCM2711");
        default                 : return result_code;
        }
      }

      class pigpiod {
      public:
        typedef boost::shared_ptr<pigpiod> sptr;

        static sptr get(std::string const & address)
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
  
        ~pigpiod()
        {
          gr::thread::scoped_lock lk(s_mtx);
  
          s_servers.erase(d_address);
  
          pigpio_stop(d_handle);
        }
  
        std::string const d_address;
  
        const int d_handle;
  
        boost::unordered_map<unsigned, gr::block *> d_sinks;
        std::list<int> d_waveforms_in_flight;
  
        gr::thread::mutex d_mtx;
  
      private:
        /*needs_lock*/ pigpiod(std::string const & address)
        : d_address(address),
          d_handle(start(address))
        { }

        static int start(std::string const & address)
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
  
        static boost::unordered_map<std::string, boost::weak_ptr<pigpiod>> s_servers;
        static gr::thread::mutex s_mtx;
      };
    }
  }
}
