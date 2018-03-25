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
  }

  if (number >= 10000) print(char('0' + ((number / 10000) % 10)));
  if (number >= 1000) print(char('0' + ((number / 1000) % 10)));
  if (number >= 100) print(char('0' + ((number / 100) % 10)));
  if (number >= 10) print(char('0' + ((number / 10) % 10)));
  print(char('0' + (number % 10)));
}
