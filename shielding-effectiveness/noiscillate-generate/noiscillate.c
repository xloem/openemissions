#include <pigpio.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main(int argc, const char ** argv)
{
	if (argc != 3) {
		fprintf(stderr, "Usage: %s port hertz\n", argv[0]);
		return -1;
	}
	
	unsigned port = strtoul(argv[1], 0, 10);
	unsigned hertz = strtoul(argv[2], 0, 10);

	if (gpioInitialise() == PI_INIT_FAILED) {
		fprintf(stderr, "Failed to initialise gpio library\n");
		return -2;
	}

	if (gpioSetMode(port, PI_OUTPUT)) {
		fprintf(stderr, "Failed to set pin mode\n");
		return -3;
	}

	int realRange = gpioGetPWMrealRange(port);
	gpioSetPWMrange(port, realRange);

	if (gpioSetPWMfrequency(port, hertz) == PI_BAD_USER_GPIO) {
		fprintf(stderr, "Failed to set gpio frequency\n");
		return -5;
	}

	if (gpioPWM(port, gpioGetPWMrange(port) / 2)) {
		fprintf(stderr, "Failed to set gpio PWM\n");
		return -4;
	}

	printf("pin %d\n", port);
	printf("pwm %d\n", gpioGetPWMdutycycle(port));
	printf("range %d\n", gpioGetPWMrange(port));
	printf("realrange %d\n", gpioGetPWMrealRange(port));
	printf("freq %d\n", gpioGetPWMfrequency(port));

	pause();

	printf("%s Terminated\n", argv[0]);

	gpioTerminate();
}
