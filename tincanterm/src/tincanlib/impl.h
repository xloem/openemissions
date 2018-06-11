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

// Set a handler to be called when the input port changes
typedef void (*implRemoteHandler)(bool trueOrFalse, unsigned long eventMicros);
void implRemoteRecvChange(implRemoteHandler handler);

// it's becoming increasingly difficult to continue due to pieces of my experience disappearing as I experience them
// i'm very slowed due to my inability to judge this small decision
// now I can't access the pieces I'm trying to use to judge it
//
// maybe I need an anonymous employee

// How do we get our handler called at the correct bit time?
// we only know the correct bit times once the handler is called
// and we know all of them.
// - api function to set an interval
//    CLARITY:
//      - pretty clear
//    SIMPLICITY:
//      - is extra func to clear interval needed?
//      - deadline doesn't need to be reset manually, which could prevent bugs if driver implementation is good
//      - driver implementation more complex due to storing interval data
//    DRIVER:
//      - need interval data
//      - need to handle 
//    APP:
//      - don't need to worry about resetting deadline
//    CLEAR but requires additional api function
// - api function to set a handler deadline
//    LESS SIMPLE because handler has dual roles ?
// - reference arg to handler deadline
//    SIMPLER because there is no additional api function
// - api function to set a general deadline
//    CLEAR due to generality; additional api function
// - api function to set a sample deadline
//    NOT AS GENERAL; HARDER TO IMPLEMENT IN DRIVER
//
// CLARITY, SIMPLICITY, EASY TO IMPLEMENT IN DRIVER, EASY TO CODE FOR IN APP

// Forward a character locally
void implLocalSend(char character);

// Read a character locally, blocking
char implLocalRecv();

// Shut down the remote communication implementation
void implRemoteDestroy();

// Shut down the local communication implementation
void implLocalDestroy();
