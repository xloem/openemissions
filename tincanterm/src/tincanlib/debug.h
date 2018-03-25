#pragma once

#include "config.h"
#include "print.h"

#if NDEBUG
// debugging not used
#define debug(arg)
#else
// print debugging strings
#define debug(arg) print(arg)
#endif
