#include "tincanterm.h"
#include "impl.h"
#include "codec.h"

#include "Arduino.h"

// State for encoding or decoding a character
struct CodingState
{
  // whether or not data is being coded
  bool active;

  // storage
  unsigned wire;

  // time to code next character
  unsigned long start;

  // which bit is being read
  unsigned bit;
};

#ifdef NDEBUG
#define debug(arg) 
#else
void debug(char message) {
  implLocalSend(message);
}
void debug(char const * message) {
  while (message[0]) {
    debug(message[0]);
    message = message + 1;
  }
}
#endif

int main()
{
  struct CodingState recvState;
  struct CodingState sendState;

  recvState.active = false;
  sendState.active = false;

  implInit();

  // Begin in resting state
  debug("-> LOW init"); debug(LINEBREAK);
  implRemoteSend(false);

  while (true)
  {
    // check time once per loop
    unsigned long msNow = implMillis();

    // Sending
    if (! sendState.active) {
      // not in the process of sending

      // Check for new data to send
      if (implLocalAvail()) {

        // Load character
        unsigned char character = implLocalRecv();
        debug("-> sending character ");
          if (character >= 100) debug('0' + ((character / 100) % 10));
          if (character >= 10) debug('0' + ((character / 10) % 10));
          debug('0' + (character % 10));
          debug(LINEBREAK);
         
        // Start sending character
        debug("-> HIGH start"); debug(LINEBREAK);
        implRemoteSend(true);

        sendState.active = true;
        sendState.wire = asciiToWire(character);
        sendState.start = msNow + BIT_MS;
        sendState.bit = 0;

      }

    } else if (msNow >= sendState.start) {
      // we're sending and it's time to send the next bit

      // an extra bit is sent to return to the resting state
      if (sendState.bit >= WIRE_BITS + 1) {

        // done sending
        sendState.active = false;

      } else {
  
        // Send the bit
        bool bitValue = (sendState.wire >> sendState.bit) & 1;
        debug("-> Bit ");
          debug('0' + sendState.bit);
          debug(": ");
          debug(bitValue ? "HIGH" : "LOW");
          debug(LINEBREAK);
        implRemoteSend(bitValue);
  
        // Update state
        sendState.start += BIT_MS;
        sendState.bit += 1;
      }
    }

    // Receiving
    if (! recvState.active) {
      // not in the process of receiving

      // Check for starting HIGH
      if (implRemoteRecv()) {
        debug("<- HIGH start"); debug(LINEBREAK);

        // Start receiving character
        recvState.active = true;
        recvState.wire = 0;
        recvState.start = msNow + SETTLE_MS + BIT_MS;
        recvState.bit = 0;
      }

    } else if (msNow >= recvState.start) {
      // we're receiving and it's time to read the next bit

      if (recvState.bit >= WIRE_BITS) {

        // done receiving
        recvState.active = false;
        
        // decode and forward to user
        implLocalSend(wireToAscii(recvState.wire));

      } else {
      
        // Read the bit
        bool bitValue = implRemoteRecv();
  
        debug("<- Bit ");
          debug('0' + recvState.bit);
          debug(": ");
          debug(bitValue ? "HIGH" : "LOW");
          debug(LINEBREAK);
        // Since wire starts all 0, we only need to change for HIGH
        if (bitValue)
          recvState.wire = recvState.wire | (1 << recvState.bit);
  
        // Update state
        recvState.start += BIT_MS;
        recvState.bit += 1;
      }
    }
  }

  implDestroy();
  return 0;
}
