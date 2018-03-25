#pragma once

#include "config.h"

#include <stdbool.h>

// Initialize the implementation
void implInit();

// Get a timestamp
unsigned long implMillis();

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

// Shut down the implementation
void implDestroy();
