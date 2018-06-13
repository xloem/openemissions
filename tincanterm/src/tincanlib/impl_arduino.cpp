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
static implRemoteRecvHandler recvHandler = 0;
static long recvDeadlineMicros;

void implRemoteInit()
{
  pinMode(INPUT_PORT, INPUT_PULL);

  #if OUTPUT_INVERT
  digitalWrite(OUTPUT_PORT, HIGH);
  #else
  digitalWrite(OUTPUT_PORT, LOW);
  #endif

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
  // can't do anything if no handler set
  if (! recvHandler) return;

  // if this function is called recursively it could cause a stack
  // overflow.  this variable prevents this.
  static bool recursing = false;
  if (recursing) return;
  recursing = true;


  // atomically read input event
  noInterrupts();
  unsigned long changeMicros = recvChangeMicros;
  unsigned long newValue = recvNewValue;
  recvChangeMicros = NOTHING_RECVD;
  long nowMicros = micros();
  interrupts();

  // dispatch handler calls in correct order
  if (recvDeadlineMicros - changeMicros > 0) {
    if (changeMicros != NOTHING_RECVD) {
      recvDeadlineMicros = recvHandler(recvPriorValue, newValue, changeMicros);
      recvPriorValue = newValue;
      changeMicros = NOTHING_RECVD;
    }
  }
  if (nowMicros - recvDeadlineMicros > 0) {
    recvDeadlineMicros = recvHandler(recvPriorValue, recvPriorValue, recvDeadlineMicros);
  }
  if (changeMicros != NOTHING_RECVD) {
    while (recvDeadlineMicros < changeMicros) {
      recvDeadlineMicros = recvHandler(recvPriorValue, recvPriorValue, recvDeadlineMicros);
    }
    recvDeadlineMicros = recvHandler(recvPriorValue, newValue, changeMicros);
    recvPriorValue = newValue;
  }

  // reset recursion prevention
  recursing = false;
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

  // coarse delay to get precise delay cycles to fit within 16 bits 
  while (micros() < time - 800);

  // now, precise delay and send
  noInterrupts();

  #if OUTPUT_INVERT
  trueOrFalse = !trueOrFalse;
  #endif
  unsigned long start = micros();
  unsigned long us = time - start;

  delayMicroseconds(us);
  
  digitalWrite(OUTPUT_PORT, trueOrFalse ? HIGH : LOW);

  unsigned long fin = micros();

  interrupts();

  poll();
}

void implLocalSend(char character)
{
  poll();

  Serial.write(character);

  poll();
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
