#pragma once

#include "config.h"

#include <stdbool.h>

// Initialize the remote communication implementation
void implRemoteInit();

// Initialize the local communication implementation
void implLocalInit();

// Set the output port high or low and get the time it was set at
unsigned long implRemoteSend(bool trueOrFalse);

// Set the output port high or low once a time is hit
void implRemoteSendSchedule(bool trueOrFalse, unsigned long scheduleMicros);

// Set a handler to be called whenever either the input port changes or a deadline passes
// Handler should return the next deadline time in microseconds.
//
// Driver implementation note:
//    this handler could call any other function.  make sure to prevent
//    stack overflow by disallowing it to be called recursively.
//
// TODO: a more robust api might be to return an event structure from one single blocking function that waits for the next input or output event, and not expose callbacks.  This might increase driver complexity but would prevent these dangerous recursion bugs.
typedef unsigned long (*implRemoteRecvHandler)(bool priorValue, bool newValue, unsigned long eventMicros);
void implRemoteRecvHandle(implRemoteRecvHandler handler, unsigned long deadlineMicros);

// Forward a character locally
void implLocalSend(char character);

// Read a character locally, blocking
char implLocalRecv();

// Shut down the remote communication implementation
void implRemoteDestroy();

// Shut down the local communication implementation
void implLocalDestroy();
