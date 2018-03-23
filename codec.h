// Defines the functions that translate between wire and ascii

#pragma once

#define WIRE_MAX 96
#define WIRE_BITS 7

char wireToAscii(unsigned wire);
unsigned asciiToWire(char ascii);
