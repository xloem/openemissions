#pragma once

/*
 * Arduino sending/receiving via handmade fiber optic cables using a PHOTORESISTOR.
 * Photoresistors are very slow, so communication is like molasses.
 *
 * A length of cable is surrounded in heat shrink tubing or any other opaque material.
 * An LED is affixed to one end.  The LED is tested that it creates significant light coming out of the opposite end of the cable.
 * A photoresistor is then affixed to this opposite end.
 *
 * Dark heat shrink tubing or some other opaque material must be placed around connections to prevent external light from entering.
 *
 * Photoresistor wiring:
 *    First leg connects to +5V
 *    Second leg connects to ARDUINO PIN 4
 *    Second leg additionally connects to two 33K Ohm resistors wired in series (end-to-end), which connect to GROUND
 *      This high-resistance connection to ground acts as pull-down on the arduino pin, making a voltage divider with the photoresistor.
 *      The resistance may need to be increased if the system reads LOW/FALSE/0 when it should not
 *      The resistance may need to be decreased if the system reads HIGH/TRUE/1 when it should not
 *      This can be checked by commenting the NDEBUG definition at the bottom of this file.
 *
 * LED wiring:
 *    Negative leg connects to GROUND
 *    Positive leg connects to a 330 Ohm resistor, which connects to ARDUINO PIN 5
 */


// which implementation to use
#define IMPL        ARDUINO

// identifier for input port
#define INPUT_PORT  4

// whether to invert the input logic
#define INPUT_INVERT false

// identifier for output port
#define OUTPUT_PORT 5

// whether to invert the output logic
#define OUTPUT_INVERT false

// Character to use for linebreak
#define LINEBREAK '\n'

// microseconds to spend sending each bit (inverse baud in us)
#define BIT_US 29000

// delay to make sure bit value has settled when reading
#define SETTLE_US 13000

// pull resistor direction for input port
#define INPUT_PULL  DOWN

// baud rate for local communications, if applicable
// note that if the device has no flow control the buffer will quickly fill, causing data loss unless rate is limited elsewhere
#define LOCAL_BAUD 9600

// unset to enable debugging output
#define NDEBUG 1
