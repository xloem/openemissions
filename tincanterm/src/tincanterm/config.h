#pragma once


// which implementation to use
#define IMPL        ARDUINO

// identifier for input port
#define INPUT_PORT  4

// identifier for output port
#define OUTPUT_PORT 5

// Character to use for linebreak
#define LINEBREAK '\n'

// milliseconds to spend sending each bit (inverse baud in ms)
#define BIT_MS 5

// delay to allow bit value to settle when reading
#define SETTLE_MS 0

// baud rate for local communications, if applicable
#define LOCAL_BAUD 9600

// pull resistor direction for input port
#define INPUT_PULL  DOWN

// set to nonzero if input should be set output-low before reading
// this drains stored charge
#define INPUT_DRAIN 0

// unset to enable debugging output
#define NDEBUG
