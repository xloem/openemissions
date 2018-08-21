// Defines the functions that translate between wire and ascii

#pragma once

#include "config.h"

// Only printable characters, space, carriage return, and ^C are used.
// Everything else is turned into a '?'.
// Printable + space are the 95 ascii chars starting at space.

#define ASCII_MIN ' '
#define ASCII_MAX (' ' + 95)

// The remaining non-sequential chars to support
const char WIRE_CTL_CHARS[] =
{
	// may be expanded if needed
	'\n', // linebreak
	'\x03' // ctrl-C
};

// Chars that might be silently ignored (not used in codec.cpp)
const char WIRE_IGN_CHARS[] =
{
  '\r' // excess carriage return
};

#define WIRE_SEQ_ASCII (ASCII_MAX - ASCII_MIN)
#define WIRE_MAX (WIRE_SEQ_ASCII + sizeof(WIRE_CTL_CHARS)/sizeof(*WIRE_CTL_CHARS))
#define WIRE_BITS 7

char wireToAscii(unsigned wire);
unsigned asciiToWire(char ascii);
