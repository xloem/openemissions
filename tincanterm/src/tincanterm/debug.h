#pragma once

#include "config.h"

#if NDEBUG

// debugging not compiled in
#define debug(arg)

#else

// output debug strings
void debug(char const * message);
void debug(char character);
void debug(long number);

#endif
