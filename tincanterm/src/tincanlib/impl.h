#pragma once

#include "config.h"

#include <stdbool.h>

// Initialize the remote communication implementation
void implRemoteInit();

// Initialize the local communication implementation
void implLocalInit();

// Get a timestamp
unsigned long implMicros();

// Set the output port high or low
void implRemoteSend(bool trueOrFalse);

// Read the input port
bool implRemoteRecv();

// Forward a character locally
void implLocalSend(char character);

// Check if any characters are available locally
bool implLocalAvail();

// Read a character locally
char implLocalRecv();

// Shut down the remote communication implementation
void implRemoteDestroy();

// Shut down the local communication implementation
void implLocalDestroy();
