// Communication implementation for bit banging on Arduinos

#include "impl.h"

#if defined(ARDUINO) && IMPL == ARDUINO

#undef true
#undef false

#include <Arduino.h>

// There's also a trick inside tincanlib.h for renaming main.
#include "tincanlib.h"

#define DOWN INPUT
#define UP INPUT_PULLUP

void implInit()
{
  Serial.begin(LOCAL_BAUD);

  pinMode(INPUT_PORT, INPUT_PULL);

  pinMode(OUTPUT_PORT, OUTPUT);
}

unsigned long implMillis()
{
  return millis();
}

void implRemoteSend(bool trueOrFalse)
{
  #if OUTPUT_INVERT
  trueOrFalse = !trueOrFalse;
  #endif

  digitalWrite(OUTPUT_PORT, trueOrFalse ? HIGH : LOW);
}

bool implRemoteRecv()
{
  #if INPUT_DRAIN
  pinMode(INPUT_PORT, OUTPUT);
  digitalWrite(INPUT_PORT, LOW);
  pinMode(INPUT_PORT, INPUT_PULL);
  #endif

  #if INPUT_INVERT
  return digitalRead(INPUT_PORT) == LOW;
  #else
  return digitalRead(INPUT_PORT) == HIGH;
  #endif
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

void setup()
{
  impl_arduino_main();
}

void loop()
{
}

#endif
