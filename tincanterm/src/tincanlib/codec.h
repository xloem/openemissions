// Defines the functions that translate between wire and ascii

#pragma once

#include "config.h"

#define WIRE_MAX 96
#define WIRE_BITS 7

char wireToAscii(unsigned wire);
unsigned asciiToWire(char ascii);
