#!/bin/bash

if test -z "$*"
then
  echo Usage: "$0" cmds ...
  echo
  echo Example:
  echo \$ while true\; do "$0" bash -i\; done
  exit 1
fi

export TERM=printer

tmpdir=$(mktemp --tmpdir -d tincan.XXX)
cmdin="$tmpdir"/in.fifo
cmdout="$tmpdir"/out.fifo

mkfifo "$cmdin"
mkfifo "$cmdout"

echo "../tincanterm > $cmdin < $cmdout"
../tincanterm > "$cmdin" < "$cmdout" &
tctpid=$!
echo

onexit() {
  wait $tctpid
  rm -rf "$tmpdir"
}

trap onexit EXIT

"$@" < "$cmdin" 2>&1 | tee "$cmdout"
