#pragma once

#include "config.h"

#include <stdbool.h>

// Initialize the implementation
void implInitInput();
void implInitOutput();

// Set the output high or low
void implSend(bool trueOrFalse);

// Read the input
bool implRead();

// Get a timestamp
unsigned long implMillis();

// Display a character locally
void implDisplay(char character);

// Wait for a character locally and return it
char implGetch();

// Shut down the implementation
void implDestroy();
