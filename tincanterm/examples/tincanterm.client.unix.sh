#!/bin/bash

# open tty on tincanterm
tty=$(tty)
../tincanterm | {
  read pts

  # ensure first time pts is opened, it is not closed until subprocesses complete
  # if it is closeed tincanterm will close
  exec > "$pts"

  # set terminal parameters to disable buffering and echoing, and pass sigint over the line
  stty -F "$pts" -echo -icanon
  stty -F "$tty" -echo -icanon -isig

  # bind stdin/stdout
  cat "$pts" > "$tty" &
  cat "$tty" > "$pts"
  wait

  # restore terminal parameters
  stty -F "$tty" echo icanon isig
}
