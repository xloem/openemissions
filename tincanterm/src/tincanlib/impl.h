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
typedef unsigned long (*implRemoteRecvHandler)(bool priorValue, bool newValue, unsigned long eventMicros);
void implRemoteRecvHandle(implRemoteRecvHandler handler, unsigned long deadlineMicros);

// it's becoming increasingly difficult to continue due to pieces of my experience disappearing as I experience them
// i'm very slowed due to my inability to judge this small decision
// now I can't access the pieces I'm trying to use to judge it
//
// maybe I need an anonymous employee

// How do we get our handler called at the correct bit time?
// we only know the correct bit times once the handler is called
// and we know all of them.
// + api function to set an interval
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
// - api function to set a handler deadline
//    CLARITY:
//      / handler has dual roles, but this could be made clear by calling to separate funcs
//    SIMPLICITY:
//      - handler has dual roles, more complex code
//      + definitely no additional callback signature needed
//    DRIVER:
//      / needs to call on a timeout, but this can be done in poll
//      + seems easy to implement: check time vs deadline, call handler early if struck
//    APP:
//      - has to reset deadline
//      -? has to interpret multiple meanings of handler function
//          -> could explore implementation in more detail
//      / can store last bit time from 
//    - I think this approach is basically the same as setting a timeout function, but it makes the code more complex by requiring
//    the handler to tell the difference between a change and a timeout, which could require extra parameters and extra comparison.
// - reference arg to handler deadline
//    SIMPLER because there is no additional api function
//    I don't think handler deadlines are the way to go UNLESS we make it clear in the handler that there is no change.
//       This would mean passing the last value, passing the reason for call, or using an unclear comparison such as checking if the
//       deadline is the bit time (could be a fluke too).
// + api function to set a general deadline
//    CLARITY:
//      + does one thing and does it well
//    SIMPLICITY
//      / does one thing and does it well
//    DRIVER:
//      / use poll function
//    APP:
//      - needs to remember samples
// + api function to set a sample deadline
//    NOT AS GENERAL; HARDER TO IMPLEMENT IN DRIVER
//    CLARITY:
//      /
//    SIMPLICITY:
//      /
//    DRIVER:
//      - needs to remember samples
//    APP:
//      /
//
// CLARITY, SIMPLICITY, EASY TO IMPLEMENT IN DRIVER, EASY TO CODE FOR IN APP
//
// I'm basically picking between a function to set a deadline, and a function to set a deadline which is called with the bit value at
// that deadline.
// At first glance the first seemss impler, api-function-wise, but the implementation of the two will be different.
//
// General deadline, app code determines sample value.
//    This will be simplest if the driver code can gaurantee that the timeout will be called prior to an event for a change that occurred
//    after the deadline.  This could be gauranteed by comparing the microseconds in the callback, and ordering the calls based on this.
//
// Sample deadline, driver code determines sample value.
//    This has the downside that sample remembering needs to be done in every implementaiton, which is a little more code.
//    The implementation would end up being almost the same as the general deadline, though.  We'd want to compare timestamps of the deadline
//    with the event.  
//    In the end they are almost the same.
// 
// The general deadline has a simpler api function.  It means bit storage 
// General deadline.
//
// Okay !!! The sample deadline doesn't actually have a more complex api function.  It can reuse the existing callback interface of the change handler.  It has a bit time and a value.  The general deadline only has an optional time.
//
// There's another issue here ... these deadlines need to be reset.  If they aren't reset fast enough, one could be missed ... ... although we
// could just call it immediately if it's late.
//
// Alternatively, if we set an interval, there is at least less handling code.  No need to reset the deadline every bit.  
// But all the interval implementation would be doing is resetting a deadline.
//
// So maybe passing the deadline as a reference param would be the way to go.  This is as fast as an interval implementation would be.
//
// How would the reference param look?
//
// Maybe I'll just try it.
//
// Hmm reference param is clearer than return value, but return value is faster because it won't require a pointer dereference.
// The speed difference there is generally negligible, but it will be at least 1 or 2 clock cycles which could affect things at whatever
// the max bitrate ends up being.
//
// I think the difference is likely negligible (1Mhz = 1 microsecond, not the order of magnitude we're concerned about now)
//
// I'm thinking I could simplify the API by _only_ using deadlines.  This would limit the utility .. for example if we wanted to average
// the bit values, it would make that impossible.  But it would simplify things a lot.
// What about noise?  Say the pin has a tendency to do sudden, brief oscillations randomly.  We could filter these out in the app code,
// in the driver code, or hope the api used by the driver allows to filter them.  If we want to filter them in app code, we'd need
// an api function that respondsn to changes.

// Forward a character locally
void implLocalSend(char character);

// Read a character locally, blocking
char implLocalRecv();

// Shut down the remote communication implementation
void implRemoteDestroy();

// Shut down the local communication implementation
void implLocalDestroy();
