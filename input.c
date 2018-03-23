#include "impl.h"
#include "codec.h"

int main()
{
	implInitInput();

  while (true)
  {
    unsigned wire;
    unsigned long start;

    // Wait for starting HIGH
    while (!implRead());

    // Mark start
    start = implMillis() + SETTLE_MS;
    wire = 0;

    // Read each bit
    for (unsigned i = 0; i < WIRE_BITS; ++ i) {

      // Wait for bit time
      start += BIT_MS;
      while (implMillis() < start);
      
      // Read bit
      wire = wire | (implRead() ? (1 << i) : 0);
    }

    // Decode and display
    implDisplay(wireToAscii(wire));
  }

	implDestroy();

	return 0;
}
