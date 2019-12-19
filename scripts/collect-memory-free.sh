#!/bin/bash

export LC_NUMERIC="C"
TIME_BETWEEN_COMMANDS=1
echo TIMESTAMP, DESCARTE, TOTAL, USED, FREE, AVAILABLE
while true; do
  info=`free -b | grep Mem`
  timestamp=`date +"%s"`
  echo $timestamp $info
  sleep $TIME_BETWEEN_COMMANDS
done
