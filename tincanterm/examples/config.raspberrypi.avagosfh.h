#pragma once

/*
 * Raspberry Pi sending/receiving via handmade fiber optic cables using a AVAGO TRANSMITTERS AND RECEIVERS.
 * I used the SFH757V and SFH250V listed in hardware.txt, which made assembly very easy.
 * They can operate at speeds far in excess of what this system can handle, but the values here add a few milliseconds of settle time
 * to allow for mistakes or simplifications (my system ran fine even without any covering on the fiber).
 *
 * A length of 1mm - 2.2mm fiber optic cable is optionally surrounded in heat shrink tubing or any other opaque material.
 * One end is connected to the SFH757V and the other to the SFH250V.
 *
 * The logic is inverted, because I have heard that microcontroller pins are better at sinking current than sourcing it.
 *
 * Receiver SFH250V wiring:
 *    Back left leg connects to GROUND
 *    Back right leg connects to RASPBERRY PI BCM_GPIO PIN 4
 *    Back right leg additionally connects to a 1M Ohm resistor, which connects to +5.0V
 *
 * Transmitter SFH757V wiring:
 *    Back left leg connects to +3.3V
 *    Back right leg connects to a 33 Ohm resistor, which connects to RASPBERRY PI BCM_GPIO PIN 5
 */


// which implementation to use
#define IMPL        PIGPIO

// identifier for input port
#define INPUT_PORT  4

// whether to invert the input logic
#define INPUT_INVERT true

// identifier for output port
#define OUTPUT_PORT 5

// whether to invert the output logic
#define OUTPUT_INVERT true

// Character to use for linebreak
#define LINEBREAK '\n'

// microseconds to spend sending each bit (inverse baud in us)
#define BIT_US 425

// delay to make sure bit value has settled when reading
#define SETTLE_US 320

// pull resistor direction for input port
#define INPUT_PULL  OFF

// set to nonzero if input should be set output-low before reading
// this drains stored charge
#define INPUT_DRAIN 0

// baud rate for local communications, if applicable
// note that if the device has no flow control the buffer will quickly fill, causing data loss unless rate is limited elsewhere
//#define LOCAL_BAUD 9600

// unset to enable debugging output
#define NDEBUG 1
