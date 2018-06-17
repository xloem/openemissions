#include "debug.h"

#include "impl.h"

void print(char character)
{
  implLocalSend(character);
}

void print(char const * message)
{
  while (message[0]) {
    print(message[0]);
    message = message + 1;
  }
}

void print(long number)
{
  if (number < 0) {
    print('-');
    number = -number;
  } else if (number == 0) {
    print('0');
    return;
  }

  char buf[32];
  unsigned idx = sizeof(buf) - 1;
  
  buf[idx] = 0;

  while (number != 0 && idx != 0) {
    idx -= 1;
    buf[idx] = (number % 10) + '0';
    number /= 10;
  }

  print(&buf[idx]);
}
