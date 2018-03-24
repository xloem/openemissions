// Communication implementation for direct GPIO on the Raspberry PI


// UNWORKING
// TODO: update to match impl_pigpiod.c

#include "impl.h"

#include <stdlib.h>

#define PIGPIO 11
#if IMPL == PIGPIO

#define _pigpio_input_pull(pull) PI_PUD_##pull
#define pigpio_input_pull(pull) _pigpio_input_pull(pull)

#include <pigpio.h>

void implInit()
{
	if (gpioInitialise() == PI_INIT_FAILED)
		exit(1);

	if (gpioSetMode(INPUT_PORT, PI_OUTPUT) != 0)
		exit(2);

	if (gpioSetMode(INPUT_PORT, PI_INPUT) != 0)
		exit(2);

	if (gpioSetPullUpDown(INPUT_PORT, pigpio_input_pull(INPUT_PULL)) != 0)
		exit(3);

	if (gpioSetMode(OUTPUT_PORT, PI_OUTPUT) != 0)
		exit(4);

	if (gpioWrite(OUTPUT_PORT, 1) != 0)
		exit(5);

	if (gpioWrite(OUTPUT_PORT, 0) != 0)
		exit(5);
	
}

void implSend(bool trueOrFalse)
{
	gpioWrite(OUTPUT_PORT, trueOrFalse ? 1 : 0);
}

bool implRead()
{
	gpioSetMode(INPUT_PORT, PI_OUTPUT);
	gpioSetMode(INPUT_PORT, PI_INPUT);
	return gpioRead(INPUT_PORT);
}

void implDestroy()
{
	gpioTerminate();
}

#endif
