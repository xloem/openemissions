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
  implInit();
  implRemoteSend(false);

  bool oldstate = false;

  double waitMin = 1;
  double settleMin = 1;

  unsigned long msStart, msNow;
  double msDelay;

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

    long msTotalComparator = implMillis();
    
    while (implMillis() - msTotalComparator < 5000)
    {
      // output summary
      print("\r___ ");
      print(long(count));
      print(" runs, low ms range = ");
      print(long(toLowMin));
      print(" - ");
      print(long(toLowSum / count + 0.5));
      print(" - ");
      print(long(toLowMax));
      print(" high ms range = ");
      print(long(toHighMin));
      print(" - ");
      print(long(toHighSum / count + 0.5));
      print(" - ");
      print(long(toHighMax));
      print(", apparent bias towards ");
      print(toHighSum > toLowSum ? "LOW" : "HIGH");
      print(", min bit ms = ");
      if (toHighMax > waitMin) waitMin = toHighMax;
      if (toLowMax > waitMin) waitMin = toLowMax;
      long lowSettle = toLowMax - toLowMin;
      long highSettle = toHighMax - toHighMin;
      long currentSettleMin = lowSettle > highSettle ? lowSettle : highSettle;
      if (lowSettle > settleMin) settleMin = lowSettle;
      if (highSettle > settleMin) settleMin = highSettle;
      print(long(toHighMax > toLowMax ? toHighMax : toLowMax) + currentSettleMin + 1);
      print(", min settle ms = ");
      print(currentSettleMin + 1);
      print("    ");

  
      // old state ends
      msStart = implMillis();
  
      bool newstate = rand();
      implRemoteSend(newstate);

      // wait for new state to settle
      msNow = msStart;
      msDelay = implMillis() - msStart;
      highCount = 0;
      lowCount = 0;
  
      while (msStart + waitMin > msNow) {
        msNow = implMillis();
        if (implRemoteRecv() != newstate) {
          print("\r_0_");
          if (newstate)
            ++ lowCount;
          else
            ++ highCount;
          msDelay = msNow - msStart;
          if (msDelay >= waitMin)
            waitMin = msDelay + 1;
          if (msDelay > 100) {
            print(" ==== State not settled after ");
            print(long(msDelay));
            print(" ms (");
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
        toHighSum += msDelay;
        if (msDelay > toHighMax)
          toHighMax = msDelay;
        if (oldstate != newstate && msDelay < toHighMin)
          toHighMin = msDelay;
      } else {
        toLowSum += msDelay;
        if (msDelay > toLowMax)
          toLowMax = msDelay;
        if (oldstate != newstate && msDelay < toLowMin)
          toLowMin = msDelay;
      }

      // update count
      ++ count;

      // update old state
      oldstate = newstate;
    }
  }

  implDestroy();
  return 0;
}
