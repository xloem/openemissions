#include "codec.h"

// The intent here is that only printable characters, space, and carriage return are used.
// Everything else is turned into a '?'.

#define WIRE_CR (WIRE_MAX - 1)
#define ASCII_MIN ' '

char wireToAscii(unsigned wire)
{
  if (wire == WIRE_CR)
    return '\n';
  else if (wire < WIRE_MAX)
    return ' ' + wire;
  else
    return '?';
}

unsigned asciiToWire(char ascii)
{
  if (ascii == '\n')
    return WIRE_CR;
  unsigned wire = ascii - ASCII_MIN;
  if (wire < WIRE_MAX)
    return wire;
  else
    return asciiToWire('?');
}
