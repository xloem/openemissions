#!/bin/bash

# run a login terminal on tincanterm
while
  ../tincanterm | {
    read pts
    sudo setsid /sbin/getty --noclear "${pts#/dev/}" dumb
  }
do
  true
done
