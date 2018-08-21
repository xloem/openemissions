// Local communication implementations for POSIX platforms

#include "impl.h"

#if unix

#include <fcntl.h>
#include <poll.h>
#include <stdio.h>
#include <stdlib.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>

static int check(int e, const char *caller)
{
  if (e < 0) {
    perror(caller);
    exit(e);
  }
  return e;
}

static int pty_master;

void implLocalInit()
{
  // create pseudoterminal
  pty_master = posix_openpt(O_RDWR | O_NOCTTY);
  check(pty_master, "impl_posix posix_openpt");

  check(grantpt(pty_master), "impl_posix grantpt");
  check(unlockpt(pty_master), "impl_posix unlockpt");

  // provide to user
  printf("%s\n", ptsname(pty_master));
  fflush(stdout);
}

unsigned long implNanos()
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  unsigned long secs = ts.tv_sec;
  return secs * 1000000000 + ts.tv_nsec;
}

void implLocalSend(char character)
{
  write(pty_master, &character, 1);
}

bool implLocalAvail()
{
  static struct pollfd readwait = { pty_master, POLLIN, 0 };
  int r = poll(&readwait, 1, 0);
  return r > 0;
}

char implLocalRecv()
{
  char character;
  int r = read(pty_master, &character, 1);
  if (r == 1) return character;

  // some error or EOF
  implLocalDestroy();
  implRemoteDestroy();
  check(r, "impl_posix read");
  exit(0);
}

void implLocalDestroy()
{
  close(pty_master);
}

#endif
