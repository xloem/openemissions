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

// I want to keep handleInputState simple and easy to review.  (that probably means understanding it and leaving accurate comments)
// There is a bug at the moment such that when a high is received to start a new message, and the last message ended in low, the high
// is not processed.
// This is a logic error, but there is a deeper bug that should be resolved first, because it will change the implementation.
// Messages are not finished until the next one begins; which means last bytes will be significantly delayed.
// I'll need some way of timing things out.
// I can see a few options, there are probably more:
// - handleInputState could watchdog; trigger after a timeout automatically
// - a function could be provided to sample a value at a time
// - a function could be provided to trigger a callback at a time
// 
// Note that arduino does not provide access to the timer interrupts.  What options are there on the arduino platform?
// The current arduino implementation polls regularly so that the callback does not happen inside an ISR.  Checking for the timeout
// could happen inside the poll routine. Either solution is viable (kind of, too bad we're polling).
//
// What would the implementations look like?
// The callback is called whenever the value is changed.
//
// Timeout function:
//    We could add bits in the timeout function when needed.  The implementation could be simple
//    if we trimmed some of the code out of the handler and put it in a third function, and called the third function from both.
//    Assuming the handler is correct, we could look into the values it has recorded in order to have highly accurate sample timing.
//    How does the problem of queing the handler compare to the problem of queueing scheduled samples?
//
//
// Scheduled sample:
//    Now the handler only responds to the starting high.
//    The scheduled sampling can store a bit properly.
//    we might run into an issue if the scheduled function is called late, if the time for the next bit has already passed by the time it
//    returns the schedule may not be set.  However, this could be an issue with the existing implementation too.
//
// Watchdog handler:
//    Is this approach simpler for the application code?
//    The only concern is that the handler would need to be able to tell the difference between a change and a repeat.
//    This would mean storing the last value, or passing another value.
//
//    This is most similar to the current implementation.
//    The driver doesn't need to understand the clock: we could just set a timeout, even change it.
//    This has the plus of not adding an additional API function, although it adds a couple parameters.  
//    I think the implementation for the application might be easier because it doesn't need to examine the handler results to determine
//    the accurate value at a given time.
//
//    There's a problem here: currently the arduino implementation calls the handler outside the ISR, which means another ISR can trigger
//    before the handler is called.
//    There's a way to handle that ... but what if the timeout is adjusted in a difficult way?
//    Maybe it would be easiest just to have a one-shot timeout, and give the driver an understanding of the clock.
//
//    I'm kind of thinking what I want is to be able to use the handler to read the starting HIGH.
//    The starting HIGH determines when we want to sample all the bits (in the starting naive implementation); each bit is sampled at
//    the high time + settle time, at bit time intervals.
//    So we have a number of times we want to know the value.  We just want some function that will be called at the bit time if there
//    is no change.  I guess I just want the API to be simple.

// called when the input port changes
void handleInputState(bool value, unsigned long time)
{
  // Reciving

  if (!recvState.active) {
    // not in the process of receiving
    
    // Check if this is the starting HIGH
    if (value) {
      debug("<- ");
        debug((long)INPUT_PORT);
        debug(" HIGH start");
        debug(LINEBREAK);

      // Start receiving character
      recvState.active = true;
      recvState.wire = 0;
      recvState.start = time + SETTLE_US + BIT_US;
      recvState.bit = 0;
    } else {
      debug("(<-");
        debug((long)INPUT_PORT);
        debug(" ignored unactive ");
        debug(value ? "high" : "low");
        debug(" @ ");
        debug((long)time);
        debug(")");
    }
  } else if (time > recvState.start) {

    // we're receiving and it's time to read the next bit
    do {

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

        // bail
        break;
  
      } else {
  
        // Read the bit
        debug("<- ");
          debug((long)INPUT_PORT);
          debug(" Bit ");
          debug((long)recvState.bit);
          debug(": ");
          debug(!value ? "HIGH" : "LOW");
          debug(" @ ");
          debug((long)recvState.start);

          debug(LINEBREAK);
        // Since wire starts all 0, we only need to change for HIGH
        if (!value)
          recvState.wire = recvState.wire | (1 << recvState.bit);
    
        // Update state
        recvState.start += BIT_US;
        recvState.bit += 1;
      }

    } while (time > recvState.start);

  } else {
    debug("(<- ");
      debug((long)INPUT_PORT);
      debug(" ignored early ");
      debug(value ? "high" : "low");
      debug(" @ ");
      debug((long)time);
      debug(")");
  }
}

int main()
{
  recvState.active = false;
  sendState.active = false;

  implLocalInit();
  implRemoteInit();

  // Begin in resting state
  debug("-> ");
    debug((long)OUTPUT_PORT);
    debug(" LOW init");
    debug(LINEBREAK);
  implRemoteSend(false);

  // attach handler
  implRemoteHandle(handleInputState);

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
      debug(" HIGH start @ ");
      debug((long)sendState.start);
      debug(LINEBREAK);
    sendState.start += BIT_US;
    sendState.bit = 0;

    debug("-> encoding as ");
      debug((long)sendState.wire);
      debug(LINEBREAK);

    // an extra bit is sent to return to the resting state
    while (sendState.bit <= WIRE_BITS) {

      // Send the bit
      bool value = (sendState.wire >> sendState.bit) & 1;

      implRemoteSchedule(value, sendState.start);

      debug("-> ");
        debug((long)OUTPUT_PORT);
        debug(" Bit ");
        debug((long)sendState.bit);
        debug(": ");
        debug(value ? "HIGH" : "LOW");
        debug(" @ ");
        debug((long)sendState.start);
        debug(LINEBREAK);

      // Update state
      sendState.start += BIT_US;
      sendState.bit += 1;
    }
  }

  implRemoteDestroy();
  implLocalDestroy();
  return 0;
}
