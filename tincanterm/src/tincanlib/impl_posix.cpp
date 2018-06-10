// Local communication implementations for POSIX platforms

// turn off stdin line buffering

#include "impl.h"

#if unix

#include <poll.h>
#include <stdio.h>
#include <stdlib.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>

static int check(int e)
{
  if (e < 0) {
    perror("posixerror");
    exit(e);
  }
  return e;
}

static struct termios termios_orig;
static int tcgetret;

void implLocalInit()
{
  // adjust terminal for convenient use
  struct termios attrs;

  tcgetret = tcgetattr(0, &attrs);
  if (tcgetret == 0) { // there will be a nonzero error if not connected to a terminal

    termios_orig = attrs;
  
    // disable canonical mode and hence line buffering
    attrs.c_lflag &= ~ICANON;
  
    // disable local echo
    attrs.c_lflag &= ~ECHO;
  
    check(tcsetattr(0, 0, &attrs));
  };
}

unsigned long implMillis()
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec * 1000 + ts.tv_nsec / 1000000;
}

void implLocalSend(char character)
{
  write(1, &character, 1);
}

bool implLocalAvail()
{
  static struct pollfd readwait = { 0, POLLIN, 0 };
  int r = poll(&readwait, 1, 0);
  return r > 0;
}

char implLocalRecv()
{
  char character;
  int r = read(0, &character, 1);
  if (r == 1) return character;

  // some error or EOF
  implLocalDestroy();
  implRemoteDestroy();
  check(r);
  exit(0);
}

void implLocalDestroy()
{
  if (tcgetret == 0) {
    // restore original terminal settings
    tcsetattr(0, 0, &termios_orig);
  }
}

#endif
