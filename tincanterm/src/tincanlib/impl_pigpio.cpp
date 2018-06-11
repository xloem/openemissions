// Communication implementation for direct GPIO on the Raspberry PI

// TODO: NONFUNCTION UPDATE FOR NEW INTERFACE

#include "impl.h"

#define PIGPIO 11
#if IMPL == PIGPIO

#include <stdio.h>
#include <stdlib.h>

#define _pigpio_input_pull(pull) PI_PUD_##pull
#define pigpio_input_pull(pull) _pigpio_input_pull(pull)

#include <pigpio.h>
// for pigpio_error
#include <pigpiod_if2.h>

static int check(int e)
{
	if (e < 0) {
		fprintf(stderr, "pigpioerror: %s\n", pigpio_error(e));
		exit(e);
	}
	return e;
}

void implRemoteInit()
{
  int r;

  // pigpio
	r = gpioInitialise();
  check(r);

  // input port
  #if INPUT_DRAIN
	r = gpioSetMode(INPUT_PORT, PI_OUTPUT);
  check(r);
  #endif

	r = gpioSetMode(INPUT_PORT, PI_INPUT);
  check(r);

	r = gpioSetPullUpDown(INPUT_PORT, pigpio_input_pull(INPUT_PULL));
	check(r);

  r = gpioRead(INPUT_PORT);
  check(r);

  // output port
	r = gpioSetMode(OUTPUT_PORT, PI_OUTPUT);
  check(r);

	r = gpioWrite(OUTPUT_PORT, 1);
  check(r);
}

void implRemoteSend(bool trueOrFalse)
{
  #if OUTPUT_INVERT
  trueOrFalse = !trueOrFalse;
  #endif

	gpioWrite(OUTPUT_PORT, trueOrFalse ? 1 : 0);
}

bool implRemoteRecv()
{
  #if INPUT_DRAIN
	gpioSetMode(INPUT_PORT, PI_OUTPUT);
  gpioWrite(INPUT_PORT, 0);
  gpioSetMode(INPUT_PORT, PI_INPUT);
  #endif

  #if INPUT_INVERT
  return ! gpioRead(INPUT_PORT);
  #else
	return gpioRead(INPUT_PORT);
  #endif
}

void implRemoteDestroy()
{
	gpioTerminate();
}

#endif
