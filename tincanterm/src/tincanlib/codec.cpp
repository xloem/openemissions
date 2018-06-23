#include "codec.h"

char wireToAscii(unsigned wire)
{
  if (wire < WIRE_SEQ_ASCII)
    return ' ' + wire;
  else if (wire < WIRE_MAX)
    return WIRE_CTL_CHARS[wire - WIRE_SEQ_ASCII];
  else
    return '?';
}

unsigned asciiToWire(char ascii)
{
  // non-sequential chars
  for (unsigned i = 0; i < WIRE_MAX - WIRE_SEQ_ASCII; ++ i)
    if (ascii == WIRE_CTL_CHARS[i])
      return i + WIRE_SEQ_ASCII;

  // sequential chars
  unsigned wire = ascii - ASCII_MIN;
  if (wire < WIRE_SEQ_ASCII)
    return wire;
  else
    return asciiToWire('?');
}
