#include "tincandiag.h"

#include "print.h"
#include "impl.h"

bool rand() {
  static unsigned char data[] = { 0x3a, 0x10, 0x32, 0x80, 0xf6, 0xd2, 0x93, 0xb4, 0x28, 0xd9, 0x2a, 0x47, 0x27, 0xf0, 0xb6, 0xe2 };
  static unsigned char idx = 0;
  static unsigned char bit = 0;

  ++ bit;
  if (bit >= 8) {
    bit = 0;
    ++ idx;
    if (idx >= sizeof(data))
      idx = 0;
  }

  return (data[idx] >> bit) & 1;
}

int main()
{
  implRemoteInit();
  implLocalInit();
  implRemoteSend(false);

  bool oldstate = false;

  double waitMin = 1;
  double settleMin = 1;

  unsigned long usStart, usNow;
  double usDelay;

  print("Beginning diagnostics.");

  while (true) {
  
    double toLowMin = 9e24;
    double toLowSum = 0;
    double toLowMax = 0;
  
    double toHighMin = 9e24;
    double toHighSum = 0;
    double toHighMax = 0;

    double count = 0;

    long highCount, lowCount;

    print("\r\n");

    unsigned long usTotalComparator = implMicros();
    
    while (implMicros() - usTotalComparator < 5000000)
    {
      // output summary
      print("\r___ ");
      print(long(count));
      print(" runs, low us range = ");
      print(long(toLowMin));
      print(" - ");
      print(long(toLowSum / count + 0.5));
      print(" - ");
      print(long(toLowMax));
      print(" high us range = ");
      print(long(toHighMin));
      print(" - ");
      print(long(toHighSum / count + 0.5));
      print(" - ");
      print(long(toHighMax));
      print(", apparent bias towards ");
      print(toHighSum > toLowSum ? "LOW" : "HIGH");
      print(", min bit us = ");
      if (toHighMax > waitMin) waitMin = toHighMax;
      if (toLowMax > waitMin) waitMin = toLowMax;
      long lowSettle = toLowMax - toLowMin;
      long highSettle = toHighMax - toHighMin;
      long currentSettleMin = lowSettle > highSettle ? lowSettle : highSettle;
      if (lowSettle > settleMin) settleMin = lowSettle;
      if (highSettle > settleMin) settleMin = highSettle;
      print(long(toHighMax > toLowMax ? toHighMax : toLowMax) + currentSettleMin + 1);
      print(", min settle us = ");
      print(currentSettleMin + 1);
      print("    ");

  
      // old state ends
      usStart = implMicros();
  
      bool newstate = rand();
      implRemoteSend(newstate);

      // wait for new state to settle
      usNow = usStart;
      usDelay = implMicros() - usStart;
      highCount = 0;
      lowCount = 0;
  
      while (usStart + waitMin > usNow) {
        usNow = implMicros();
        if (implRemoteRecv() != newstate) {
          print("\r_0_");
          if (newstate)
            ++ lowCount;
          else
            ++ highCount;
          usDelay = usNow - usStart;
          if (usDelay >= waitMin)
            waitMin = usDelay + 1;
          if (usDelay > 100000) {
            print(" ==== State not settled after ");
            print(long(usDelay));
            print(" us (");
            print(lowCount);
            print(" low, ");
            print(highCount);
            print(" high); is input port ");
            print((long)INPUT_PORT);
            print(" correctly connected to output port ");
            print((long)OUTPUT_PORT);
            print("? ====\n");
            break;
          }
        } else {
          if (newstate)
            ++ highCount;
          else
            ++ lowCount;
          print("\r_1_");
        }
      }
  
      if (newstate) {
        toHighSum += usDelay;
        if (usDelay > toHighMax)
          toHighMax = usDelay;
        if (oldstate != newstate && usDelay < toHighMin)
          toHighMin = usDelay;
      } else {
        toLowSum += usDelay;
        if (usDelay > toLowMax)
          toLowMax = usDelay;
        if (oldstate != newstate && usDelay < toLowMin)
          toLowMin = usDelay;
      }

      // update count
      ++ count;

      // update old state
      oldstate = newstate;
    }
  }

  implRemoteDestroy();
  implLocalDestroy();
  return 0;
}
