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

// Shut down the implementation
void implDestroy();
