#pragma once

// which implementation to use
#define IMPL        PIGPIOD

// identifier for input port
#define INPUT_PORT  4

// identifier for output port
#define OUTPUT_PORT 5

// pull resistor direction for input port
#define INPUT_PULL  DOWN

// set to nonzero if input should be set output-low before reading
// this drains stored charge
#define INPUT_DRAIN 1

