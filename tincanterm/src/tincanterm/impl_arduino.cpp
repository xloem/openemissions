// Communication implementation for bit banging on Arduinos

#include "impl.h"

#if defined(ARDUINO) && IMPL == ARDUINO

#undef true
#undef false

#include <Arduino.h>

#define DOWN 1
#define UP 2

void implInit()
{
  Serial.begin(LOCAL_BAUD);

  pinMode(INPUT_PORT, INPUT_PULL == UP ? INPUT_PULLUP : INPUT);
  #if INPUT_DRAIN
    #error INPUT_DRAIN unimplemented on arduino
  #endif

  pinMode(OUTPUT_PORT, OUTPUT);
}

unsigned long implMillis()
{
  return millis();
}

void implRemoteSend(bool trueOrFalse)
{
  digitalWrite(OUTPUT_PORT, trueOrFalse ? HIGH : LOW);
}

bool implRemoteRecv()
{
  return digitalRead(INPUT_PORT) == HIGH;
}

void implLocalSend(char character)
{
  Serial.write(character);
}

bool implLocalAvail()
{
  return Serial.available() > 0;
}

char implLocalRecv()
{
  return Serial.read();
}

void implDestroy()
{
}

#endif
