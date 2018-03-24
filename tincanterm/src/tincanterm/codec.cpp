#include "codec.h"

// The intent here is that only printable characters, space, and carriage return are used.
// Everything else is turned into a '?'.

#define WIRE_LINEBREAK (WIRE_MAX - 1)
#define ASCII_MIN ' '

char wireToAscii(unsigned wire)
{
  if (wire == WIRE_LINEBREAK)
    return LINEBREAK;
  else if (wire < WIRE_MAX)
    return ' ' + wire;
  else
    return '?';
}

unsigned asciiToWire(char ascii)
{
  if (ascii == LINEBREAK)
    return WIRE_LINEBREAK;
  unsigned wire = ascii - ASCII_MIN;
  if (wire < WIRE_MAX)
    return wire;
  else
    return asciiToWire('?');
}
