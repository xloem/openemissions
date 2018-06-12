#!/bin/bash
dev=/dev/ttyACM*

stty -icanon -echo

stty -F $dev -icanon -echo min 1 9600
cat $dev &
catpid=$!

onexit() {
  stty icanon echo
  kill "$catpid"
}

trap onexit EXIT
cat > $dev
