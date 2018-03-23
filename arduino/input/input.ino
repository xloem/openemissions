#include "../../impl_arduino.cpp"

extern "C" {

#define main _main

#include "../../codec.c"
#include "../../input.c"

}

void setup() {
  _main();
}
