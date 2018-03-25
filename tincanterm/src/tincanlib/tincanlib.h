#pragma once

// This dummy header file is just so the arduino sketch builder can find the library

#ifdef ARDUINO
  #define main impl_arduino_main
  int main();
#endif
