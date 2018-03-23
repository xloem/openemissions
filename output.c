#include "impl.h"
#include "codec.h"

int main()
{
  implInitOutput();

  // Begin in resting state
  implSend(false);

  while (true)
  {
    unsigned wire;
    unsigned long start;

    // Wait for a byte and encode
    wire = asciiToWire(implGetch());

    // Mark start
    start = implMillis();
    implSend(true);

    // Send each bit
    // An extra bit is sent to return to the resting state
    for (unsigned i = 0; i < WIRE_BITS + 1; ++ i) {

      // Wait for bit time
      start += BIT_MS;
      while (implMillis() < start);

      // Send bit
      implSend((wire >> i) & 1);
    }
  }

  implDestroy();

  return 0;
}
