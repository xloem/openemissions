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

  long place = 1000000000L;
  while (number < place) place /= 10;
  while (place) {
    int digit = number / place;
    print(char('0' + digit));
    number -= digit * place;
    place /= 10;
  }
}
