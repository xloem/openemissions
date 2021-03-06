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

  // time to code next character
  unsigned long start;

  // which bit is being read
  unsigned bit;
};

int main()
{
  struct CodingState recvState;
  struct CodingState sendState;

  recvState.active = false;
  sendState.active = false;

  implRemoteInit();
  implLocalInit();

  // Begin in resting state
  debug("-> ");
    debug((long)OUTPUT_PORT);
    debug(" LOW init");
    debug(LINEBREAK);
  implRemoteSend(false);

  while (true)
  {
    // Sending
    if (! sendState.active) {
      // not in the process of sending

      // Check for new data to send
      if (implLocalAvail()) {

        bool ignored = false;

        // Load character
        char character = implLocalRecv();

        // Check if character is squelched
        for (unsigned i = 0; i < sizeof(WIRE_IGN_CHARS); ++ i)
          if (character == WIRE_IGN_CHARS[i])
            ignored = true;

        if (!ignored) {
          debug("-> sending character ");
            debug((long)character);
            debug(LINEBREAK);
           
          // Start sending character
          debug("-> ");
            debug((long)OUTPUT_PORT);
            debug(" HIGH start");
            debug(LINEBREAK);
  
          sendState.active = true;
          sendState.wire = asciiToWire(character);
          sendState.start = implMicros() + BIT_US;
          implRemoteSend(true);
          sendState.bit = 0;
          debug("-> encoding as ");
            debug((long)sendState.wire);
            debug(LINEBREAK);
        }
      }

    } else if (implMicros() >= sendState.start) {
      // we're sending and it's time to send the next bit

      // an extra bit is sent to return to the resting state
      if (sendState.bit >= WIRE_BITS + 1) {

        // done sending
        sendState.active = false;

      } else {
  
        // Send the bit
        bool bitValue = (sendState.wire >> sendState.bit) & 1;
        implRemoteSend(bitValue);
        debug("-> ");
          debug((long)OUTPUT_PORT);
          debug(" Bit ");
          debug((long)sendState.bit);
          debug(": ");
          debug(bitValue ? "HIGH" : "LOW");
          debug(LINEBREAK);
  
        // Update state
        sendState.start += BIT_US;
        sendState.bit += 1;
      }
    }

    // Receiving
    if (! recvState.active) {
      // not in the process of receiving

      // Check for starting HIGH
      if (implRemoteRecv()) {
        debug("<- ");
          debug((long)INPUT_PORT);
          debug(" HIGH start");
          debug(LINEBREAK);

        // Start receiving character
        recvState.active = true;
        recvState.wire = 0;
        recvState.start = implMicros() + SETTLE_US + BIT_US;
        recvState.bit = 0;

      } else {
        debug(implRemoteRecv() ? "(1)" : "(0)");
      }

    } else if (implMicros() >= recvState.start) {
      // we're receiving and it's time to read the next bit

      if (recvState.bit >= WIRE_BITS) {

        // done receiving
        recvState.active = false;
        
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
        bool bitValue = implRemoteRecv();
  
        debug("<- ");
          debug((long)INPUT_PORT);
          debug(" Bit ");
          debug((long)recvState.bit);
          debug(": ");
          debug(bitValue ? "HIGH" : "LOW");
          debug(LINEBREAK);
        // Since wire starts all 0, we only need to change for HIGH
        if (bitValue)
          recvState.wire = recvState.wire | (1 << recvState.bit);
  
        // Update state
        recvState.start += BIT_US;
        recvState.bit += 1;
      }
    } else {
      debug(implRemoteRecv() ? "(1)" : "(0)");
    }
  }

  implRemoteDestroy();
  implLocalDestroy();
  return 0;
}
