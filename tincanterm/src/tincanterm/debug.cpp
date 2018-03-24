#include "debug.h"

#include "impl.h"

#if ! NDEBUG

void debug(char character)
{
  implLocalSend(character);
}

void debug(char const * message)
{
  while (message[0]) {
    debug(message[0]);
    message = message + 1;
  }
}

void debug(long number)
{
  if (number < 0) {
    debug('-');
    number = -number;
  }

  if (number > 10000) debug(char('0' + ((number / 10000) % 10)));
  if (number > 1000) debug(char('0' + ((number / 1000) % 10)));
  if (number > 100) debug(char('0' + ((number / 100) % 10)));
  if (number > 10) debug(char('0' + ((number / 10) % 10)));
  debug(char('0' + (number % 10)));
}

#endif
