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

#define NOTHING_RECVD 0
static volatile unsigned long recvMicros = NOTHING_RECVD;
static volatile bool recvValue;
static implRemoteHandler recvHandler;

void implRemoteInit()
{
  pinMode(INPUT_PORT, INPUT_PULL);

  pinMode(OUTPUT_PORT, OUTPUT);

  recvHandler = 0;
}

void implLocalInit()
{
  Serial.begin(LOCAL_BAUD);
}

static void queueRecvHandler()
{
  recvMicros = micros();
  #if INPUT_INVERT
    recvValue = ! digitalRead(INPUT_PORT);
  #else
    recvValue = digitalRead(INPUT_PORT);
  #endif
}

void implRemoteHandle(implRemoteHandler handler)
{
  noInterrupts();

  recvHandler = handler;

  attachInterrupt(digitalPinToInterrupt(INPUT_PORT), queueRecvHandler, CHANGE);

  interrupts();
}

static void poll()
{
  if (recvMicros != NOTHING_RECVD) {
    noInterrupts();
    unsigned long micros = recvMicros;
    recvMicros = NOTHING_RECVD;
    interrupts();

    recvHandler(recvValue, micros);

  }
}

unsigned long implRemoteSend(bool trueOrFalse)
{
  unsigned long ret;

  poll();

  noInterrupts();

  #if OUTPUT_INVERT
  trueOrFalse = !trueOrFalse;
  #endif

  ret = micros();

  digitalWrite(OUTPUT_PORT, trueOrFalse ? HIGH : LOW);

  interrupts();

  poll();

  return ret;
}

void implRemoteSchedule(bool trueOrFalse, unsigned long time)
{
  poll();

  noInterrupts();

  #if OUTPUT_INVERT
  trueOrFalse = !trueOrFalse;
  #endif

  unsigned long us = time - micros();
  #if AVR
    unsigned long cycles = (F_CPU / 1000000) * us;
    _delay_loop_2(cycles / 4);
  #else
    delayMicroseconds(us);
  #endif
  
  digitalWrite(OUTPUT_PORT, trueOrFalse ? HIGH : LOW);

  interrupts();

  poll();
}

void implLocalSend(char character)
{
  //poll();

  Serial.write(character);

  //poll();
}

char implLocalRecv()
{
  char character;

  do {
    poll();

    character = Serial.read();
  } while (character == -1);

  poll();

  return character;
}

void implRemoteDestroy()
{
}

void implLocalDestroy()
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
