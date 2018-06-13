#include "tincanterm.h"

#include "codec.h"
#include "debug.h"
#include "impl.h"

// State for encoding or decoding a character
struct CodingState
{
  // whether or not data is being coded
  bool active;

  // storage
  unsigned wire;

  // time to process next bit
  unsigned long start;

  // which bit is being read
  unsigned bit;
};

struct CodingState recvState;
struct CodingState sendState;

// called when the input port changes or a deadline is hit
unsigned long handleInputState(bool priorValue, bool newValue, unsigned long time)
{
  // Receiving

  if (!recvState.active) {
    // not in the process of receiving
    
    // Check if this is the starting HIGH
    if (newValue) {
      debug("<- ");
        debug((long)INPUT_PORT);
        debug(" HI start");
        debug(" @ ");
        debug((long)time);
        debug(LINEBREAK);

      // Start receiving character
      recvState.active = true;
      recvState.wire = 0;
      recvState.start = time + SETTLE_US + BIT_US;
      recvState.bit = 0;
    } else {
      debug("(<-");
        debug((long)INPUT_PORT);
        debug(" ignored inactive ");
        debug(newValue ? "high" : "low");
        debug(" @ ");
        debug((long)time);
        debug(")");

      recvState.start = time + 60000000;
    }
  } else if (time >= recvState.start) {

    // we're receiving and it's time to read the next bit
    if (recvState.bit >= WIRE_BITS) {

      // done receiving
      recvState.active = false;
      recvState.start = time + 60000000;

      // decode and forward to user
      debug("<- received ");
        debug((long)recvState.wire);
        debug(LINEBREAK);
      char character = wireToAscii(recvState.wire);
      debug("<- decoded to ");
        debug((long)character);
        debug(LINEBREAK);
      implLocalSend(character);

    } else {

      // Read the bit
      debug("<- ");
        debug((long)INPUT_PORT);
        debug(" Bit ");
        debug((long)recvState.bit);
        debug(": ");
        debug(priorValue ? "HI" : "LO");
        debug(" @ ");
        debug((long)recvState.start);
        debug(LINEBREAK);
      // Since wire starts all 0, we only need to change for HIGH
      if (priorValue)
        recvState.wire = recvState.wire | (1 << recvState.bit);
  
      // Update state
      recvState.start += BIT_US;
      recvState.bit += 1;
    }

  } else {
    debug("(<- ");
      debug((long)INPUT_PORT);
      debug(" ignored early ");
      debug(newValue ? "high" : "low");
      debug(" @ ");
      debug((long)time);
      debug(")");
  }

  debug("(deadline = ");
    debug((long)recvState.start);
    debug(")");
  return recvState.start;
}

int main()
{
  recvState.active = false;
  sendState.active = false;

  implLocalInit();
  implRemoteInit();

  // Begin in resting state
  unsigned long startMicros = implRemoteSend(false);

  // attach handler
  implRemoteRecvHandle(handleInputState, startMicros + 60000000);

  debug("-> ");
    debug((long)OUTPUT_PORT);
    debug(" LO init @ ");
    debug((long)startMicros);
    debug(LINEBREAK);

  while (true)
  {
    // Sending
    
    // Load character
    char character = implLocalRecv();
    debug("-> sending character ");
      debug((long)character);
      debug(LINEBREAK);
     
    // Start sending character

    sendState.active = true;
    sendState.wire = asciiToWire(character);
    sendState.start = implRemoteSend(true);
    debug("-> ");
      debug((long)OUTPUT_PORT);
      debug(" HI start @ ");
      debug((long)sendState.start);
      debug(LINEBREAK);
    sendState.start += BIT_US;
    sendState.bit = 0;

    debug("-> encoding as ");
      debug((long)sendState.wire);
      debug(LINEBREAK);

    // one extra bit will be sent to return to the resting state
    while (sendState.bit <= WIRE_BITS) {

      // Send the current bit
      bool value = (sendState.wire >> sendState.bit) & 1;

      debug("-> ");
        debug((long)OUTPUT_PORT);
        debug(" Bit ");
        debug((long)sendState.bit);
        debug(": ");
        debug(value ? "HI" : "LO");
        debug(" @ ");
        debug((long)sendState.start);
        debug(LINEBREAK);

      implRemoteSendSchedule(value, sendState.start);

      // Update state
      sendState.start += BIT_US;
      sendState.bit += 1;
    }
  }

  implRemoteDestroy();
  implLocalDestroy();
  return 0;
}
