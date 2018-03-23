// Communication implementation for bit banging on Arduinos

extern "C" {

#include "impl.h"

}

#if defined(ARDUINO) && IMPL == ARDUINO

#undef true
#undef false

#include <Arduino.h>

extern "C" {

static void implInit()
{
  Serial.begin(LOCAL_BAUD);
}

#define DOWN 1
#define UP 2

void implInitInput()
{
  implInit();

  pinMode(INPUT_PORT, INPUT_PULL == UP ? INPUT_PULLUP : INPUT);
  #if INPUT_DRAIN
    #error INPUT_DRAIN unimplemented on arduino
  #endif
} 

void implInitOutput()
{
  implInit();

  pinMode(OUTPUT_PORT, OUTPUT);
}

void implSend(bool trueOrFalse)
{
  digitalWrite(OUTPUT_PORT, trueOrFalse ? HIGH : LOW);
}

bool implRead()
{
  return digitalRead(INPUT_PORT) == HIGH;
}

unsigned long implMillis()
{
  return millis();
}

void implDisplay(char character)
{
  Serial.write(character);
}

char implGetch()
{
  return Serial.read();
}

void implDestroy()
{
}

}

#endif
