// Communication implementation for daemon-managed GPIO on a Raspberry PI

#include "impl.h"

#define PIGPIOD 12
#if IMPL == PIGPIOD

#include <stdio.h>
#include <stdlib.h>

#include <pigpiod_if2.h>

#define _pigpio_input_pull(pull) PI_PUD_##pull
#define pigpio_input_pull(pull) _pigpio_input_pull(pull)

int pi = 0;

static int check(int e)
{
	if (e < 0) {
		fprintf(stderr, "pigpio daemon error: %s\n", pigpio_error(e));
		exit(e);
	}
	return e;
}

static void init()
{
	pi = pigpio_start(0, 0);
	check(pi);
}


void implInitInput()
{
  int r;

  init();

	#if INPUT_DRAIN
	r = set_mode(pi, INPUT_PORT, PI_OUTPUT);
	check(r);
	#endif

	r = set_mode(pi, INPUT_PORT, PI_INPUT);
	check(r);

	r = set_pull_up_down(pi, INPUT_PORT, pigpio_input_pull(INPUT_PULL));
	check(r);

	r = gpio_read(pi, INPUT_PORT);
	check(r);
}

void implInitOutput()
{
	int r;

  init();

	r = set_mode(pi, OUTPUT_PORT, PI_OUTPUT);
	check(r);

	r = gpio_write(pi, OUTPUT_PORT, 0);
	check(r);
}

void implSend(bool trueOrFalse)
{
	gpio_write(pi, OUTPUT_PORT, trueOrFalse ? 1 : 0);
}

bool implRead()
{
	#if INPUT_DRAIN
	set_mode(pi, INPUT_PORT, PI_OUTPUT);
  gpio_write(pi, INPUT_PORT, 0);
	set_mode(pi, INPUT_PORT, PI_INPUT);
  #endif
	return gpio_read(pi, INPUT_PORT);
}

void implDestroy()
{
	pigpio_stop(pi);
	pi = 0;
}

#endif
