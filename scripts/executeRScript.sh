#!/bin/bash

R_ALGORITHM_PATH=$1
R_EXEC_DIR=$2
TMP_DIR_PATH=$3
METADATA_DIR_PATH=$4

NUMBER_OF_TIMEOUTS=0

Rscript $R_ALGORITHM_PATH $R_EXEC_DIR $TMP_DIR_PATH > $METADATA_DIR_PATH/out.log 2> $METADATA_DIR_PATH/error.log
PROCESS_OUTPUT=$?

echo "RScript_process_output=$PROCESS_OUTPUT"
if [ $PROCESS_OUTPUT -eq 124 ]
then
  NUMBER_OF_TIMEOUTS=$(($NUMBER_OF_TIMEOUTS+1))
  echo "NUMBER OF TIMEOUTS $NUMBER_OF_TIMEOUTS"
  exit 124
elif [ $PROCESS_OUTPUT -ne 0 ]
then
  exit 1
else
  exit 0
fi
