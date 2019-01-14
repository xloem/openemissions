#!/bin/bash

if [ "$1" == "" ]
then
  echo No deadline specified.
  guess=$(date +%s)
  guess=$((guess / 600 * 600 + 900))
  echo Try $((guess)) !
  exit -1
fi

start=$(($1))

if [ "$2" != "" ]
then
  echo == EMITTER ==
  start_emit() { raspi-gpio set 4 op dl; }
  stop_emit() { raspi-gpio set 4 op dh; }
  start_rec_emit() { true; }
  start_rec_noise() { true; }
  stop_emit
else
  start_emit() { true; }
  stop_emit() { true; }
  start_rec_emit() { ./1-1-prof-env -o ${start}-emit.noisep -n 1 0 1; }
  start_rec_noise() { ./1-1-prof-env -o ${start}-noise.noisep -n 1 0; }
  . /usr/local/bin/thisroot.sh
fi

# stages:
# 1. 00 min: noise recording starts and runs for a little more than 600
# 2. 10 min: emitter turns on
# 3. 15 min: emitter recording starts and runs for a little more than 600
# 4. 25 min: emitter turns off
# 5. 30 min: goto 1.

waitforuntil()
{
  name="$1"
  deadline=$(($2))
  while (( $(date +%s) < deadline))
  do
    left=$((deadline - $(date +%s)))
    mins=$((left / 60));
    secs=$((left % 60));
    printf "$name %0.2d:%0.2d\r" "$mins" "$secs"
  done
  printf "\n"
}

next=$((start))
stage=0
while ((next < $(date +%s) + 60*5))
do
  next=$((next + 60*30))
done
while true
do
  waitforuntil "Noise recording" $((next + 0))
  start_rec_noise # until a little more than 60*5
  waitforuntil "Emitter ON" $((next + 60*10))
  start_emit
  waitforuntil "Emitter recording" $((next + 60*15))
  start_rec_emit # until a little more than 60*20
  waitforuntil "Emitter OFF" $((next + 60*25))
  stop_emit
  next=$((next + 60*30))
done

