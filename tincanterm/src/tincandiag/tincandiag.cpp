#include "tincandiag.h"

#include "print.h"
#include "impl.h"

bool rand()
{
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

// the system behaves as both a receiver and a sender
// it may be self-wired or wired to another instance elsewhere

#if defined(BIT_NS)
#undef BIT_NS
#endif
const long BIT_NS = 200000;
const bool BIT_CODE[] = { true, false };

const unsigned OUTPUT_HEIGHT = 8;
const unsigned OUTPUT_WIDTH = 80;

const long DECIMATION = BIT_NS * sizeof(BIT_CODE) / OUTPUT_WIDTH;

const long WAVE_ACCUMULATION = OUTPUT_WIDTH * 8192;

char outputMatrix[OUTPUT_WIDTH][OUTPUT_HEIGHT];

void outWave(double * accumulation, long * counts)
{
  double min = 9e24;
  double max = -9e24;
  for (unsigned i = 0; i < OUTPUT_WIDTH; i += 1)
  {
    if (counts[i] == 0)
    {
      continue;
    }
    accumulation[i] = accumulation[i] / counts[i];
    if (accumulation[i] < min)
      min = accumulation[i];
    if (accumulation[i] > max)
      max = accumulation[i];
  }
  for (unsigned i = 0; i < OUTPUT_WIDTH; i += 1)
  {
    double fVal = accumulation[i];
    // extrema are treated specially to make them visually identifiable
    if (fVal == min) {
      fVal = 0;
    } else if (fVal == max) {
      fVal = OUTPUT_HEIGHT - 1;
    } else if (fVal == 0) {
      
    } else {
      fVal = ((accumulation[i] - min) * (OUTPUT_HEIGHT - 3)) / (max - min) + 1;
    }
    unsigned uVal = fVal;
    fVal -= uVal;
    char chr;
    if (fVal < 1/3.)
      chr = '_';
    else if (fVal < 2/3.)
      chr = '-';
    else
      chr = '`';
    outputMatrix[i][uVal] = chr;
  }
  for (int y = OUTPUT_HEIGHT - 1; y >= 0; y -= 1)
  {
    for (unsigned x = 0; x < OUTPUT_WIDTH; x += 1)
    {
      print(outputMatrix[x][y]);
      outputMatrix[x][y] = ' ';
    }
    print('\n');
  }
}

int main()
{
  implRemoteInit();
  implLocalInit();
  implRemoteSend(false);

  print("Beginning diagnostics.");

  for (unsigned x = 0; x < OUTPUT_WIDTH; x += 1)
  {
    for (unsigned y = 0; y < OUTPUT_HEIGHT; y += 1)
    {
      outputMatrix[x][y] = ' ';
    }
  }

  double recvWaveValues[OUTPUT_WIDTH];
  double recvWaveChanges[OUTPUT_WIDTH];
  long recvWaveCounts[OUTPUT_WIDTH];
  long recvWaveTotal = WAVE_ACCUMULATION;

  //bool const * recvCode = BIT_CODE;
  unsigned recvCodeLen = sizeof(BIT_CODE);
  //long recvBitTime = BIT_NS;
  long recvStartTime = implNanos();
  long recvNextTime = recvStartTime;
  bool recvLast = implRemoteRecv();
  long recvLastTime = recvStartTime;
 

  bool const * sendCode = BIT_CODE;
  unsigned sendCodeLen = sizeof(BIT_CODE);
  long sendBitTime = BIT_NS;
  long sendNextTime = recvStartTime + BIT_NS;
  unsigned sendIdx = 0;


  while (true)
  {

    long time = implNanos();
    if (time - sendNextTime >= 0)
    {
      implRemoteSend(sendCode[sendIdx]);

      sendNextTime += sendBitTime;
      sendIdx += 1;

      if (sendIdx >= sendCodeLen)
      {
        sendIdx = 0;
      }
    }
    else
    {
      bool recvCur = implRemoteRecv();
      if (recvWaveTotal >= WAVE_ACCUMULATION)
      {
        print("wave\r\n");
        if (recvWaveCounts[0] > 0 || recvWaveCounts[1] > 0)
        {
          outWave(recvWaveValues, recvWaveCounts);
          outWave(recvWaveChanges, recvWaveCounts);
        }

        recvWaveTotal = 0;
        for (unsigned i = 0; i < OUTPUT_WIDTH; i += 1)
        {
          recvWaveValues[i] = 0;
          recvWaveChanges[i] = 0;
          recvWaveCounts[i] = 0;
        }
      }
      while (time - recvNextTime >= 0)
      {
        recvStartTime = recvNextTime;
        recvNextTime += recvCodeLen * BIT_NS;
      }
      unsigned long recvWaveIdx = (time - recvStartTime) * OUTPUT_WIDTH / (BIT_NS * recvCodeLen);
      if (recvWaveIdx < 0 || recvWaveIdx >= OUTPUT_WIDTH) {
        print("oops");
      } else {
        // we care about the times that we change back
        // so, the time we change may oscillate, which we hope to make an average curve
        // but aeraging could also come from flutterig --
        // when we for example change from HIGH to LOW, then back to HIGH, then back to LOW
        // we mark so many HIGHs in the accumulation
        // but it's meaningful whether they wre HIGHs that have last a long time
        // or HIGHS that came after a brief LOW
        // how can we record the difference?
        // we'd like to know perhaps how many lows have preceded this high
        // i.e. how many changes there have been
        // we mark 38 highs ... where were the previous lows?
        // for 34 of them the lows were in the previous wave
        // for 2 of them the lows were in this wave, recent
        // perhaps a time spent in the bit could be meaningful to record
        // another thing that could be recorded is whether the change occured here
        // for example, rather than recording the value, I could record the value delta
        // 1, 0, or -1
        if (recvCur)
          recvWaveValues[recvWaveIdx] += 1;
        if (recvCur != recvLast)
          recvWaveChanges[recvWaveIdx] += 1.0 / (time - recvLastTime);
        recvWaveCounts[recvWaveIdx] += 1;
        recvWaveTotal += 1;
      }
      recvLast = recvCur;
      recvLastTime = time;


/*
      if ((time - recvLoopUpdateTime >= 0 || time - recvLoopDoneTime >= 0) && timesCount[0] + timesCount[1] > 0)
      {
        double avgTimes[2] = {sumTimes[0] / timesCount[0], sumTimes[1] / timesCount[1]};
        print("\r0: ");
          print(minTimes[0]);
          print("-");
          print(long(avgTimes[0] + 0.5));
          print("-");
          print(maxTimes[0]);
          print(" x");
          print(long(failsCount[0]));
          print("/");
          print(long(timesCount[0]));
        print(" 1: ");
          print(minTimes[1]);
          print("-");
          print(long(avgTimes[1] + 0.5));
          print("-");
          print(maxTimes[1]);
          print(" x");
          print(long(failsCount[1]));
          print("/");
          print(long(timesCount[1]));
        print(" bias towards ");
          if (sumTimes[0] < sumTimes[1])
          {
            print("LOW");
          }
          else if (sumTimes[1] < sumTimes[0])
          {
            print("HIGH");
          }
          else
          {
            print("NONE");
          }

        double avgTime = (sumTimes[0] + sumTimes[1]) / (timesCount[0] + timesCount[1]);
        long minTime = minTimes[0] < minTimes[1] ? minTimes[0] : minTimes[1];
        long maxTime = maxTimes[0] > maxTimes[1] ? maxTimes[0] : maxTimes[1];
        
        print(" x: ");
          print(minTime);
          print("-");
          print(long(avgTime + 0.5));
          print("-");
          print(maxTime);

        print(" jitter: ");
          print(long(maxTime - minTime));

        recvLoopUpdateTime = time + RECV_UPDATE_NS;
      }
      if (time - recvLoopDoneTime >= 0)
      {
        print("\r\n");

        minTimes[0] = 0x8fff;
        minTimes[1] = 0x8fff;
        maxTimes[0] = 0;
        maxTimes[1] = 0;
        sumTimes[0] = 0;
        sumTimes[1] = 0;
        timesCount[0] = 0;
        timesCount[1] = 0;
        failsCount[0] = 0;
        failsCount[1] = 0;

        recvLoopDoneTime = time + RECV_LOOP_NS;
        recvLastDelta = 0;
      }
      if (recvCur != recvLast)
      {
        long recvCurTime = time;
        long recvCurDelta = recvCurTime - recvLastTime;
        int idx = recvCur ? 1 : 0;

        if (recvLastDelta != 0)
        {
          if (recvCurDelta < minTimes[idx])
          {
            minTimes[idx] = recvCurDelta;
          }
          if (recvCurDelta > maxTimes[idx])
          {
            maxTimes[idx] = recvCurDelta;
          }
  
          sumTimes[idx] += recvCurDelta;
          timesCount[idx] += 1;
  
          if (recvCurDelta + recvLastDelta <= BIT_NS && recvCurDelta < recvLastDelta)
          {
            failsCount[idx] += 1;
          }
        }

        recvLast = recvCur;
        recvLastTime = recvCurTime;;
        recvLastDelta = recvCurDelta;
      }
*/
    }
  }

  implRemoteDestroy();
  implLocalDestroy();
  return 0;
}
