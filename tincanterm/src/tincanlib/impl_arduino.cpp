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
static volatile unsigned long recvChangeMicros = NOTHING_RECVD;
static volatile bool recvNewValue = false;
static bool recvPriorValue = false;
static implRemoteRecvHandler recvHandler;
static long recvDeadlineMicros;

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
  recvChangeMicros = micros();
  recvNewValue = digitalRead(INPUT_PORT);
  #if INPUT_INVERT
    recvNewValue = ! recvNewValue;
  #endif
}

void implRemoteRecvHandle(implRemoteRecvHandler handler, unsigned long deadline)
{
  noInterrupts();

  recvHandler = handler;

  attachInterrupt(digitalPinToInterrupt(INPUT_PORT), queueRecvHandler, CHANGE);

  interrupts();

  recvDeadlineMicros = deadline;
}

static void poll()
{
  noInterrupts();
  unsigned long changeMicros = recvChangeMicros;
  unsigned long newValue = recvNewValue;
  recvChangeMicros = NOTHING_RECVD;
  long nowMicros = micros();
  interrupts();

  if (recvDeadlineMicros - changeMicros > 0) {
    if (changeMicros != NOTHING_RECVD) {
      recvDeadlineMicros = recvHandler(recvPriorValue, newValue, changeMicros);
      recvPriorValue = newValue;
      changeMicros = NOTHING_RECVD;
    }
  }
  // TODO: seems only to work with this check for recvDeadlineMicros != 0
  // probably needed because poll() is getting called before these variables are set somehow
  if (nowMicros - recvDeadlineMicros > 0 && recvDeadlineMicros != 0) {
    recvDeadlineMicros = recvHandler(recvPriorValue, recvPriorValue, recvDeadlineMicros);
  }
  return;
  if (changeMicros != NOTHING_RECVD) {
    recvDeadlineMicros = recvHandler(recvPriorValue, newValue, changeMicros);
    recvPriorValue = newValue;
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

void implRemoteSendSchedule(bool trueOrFalse, unsigned long time)
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
