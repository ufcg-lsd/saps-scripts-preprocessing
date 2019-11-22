#!/bin/bash

## This script pre-processes the data downloaded by the InputDownload phase

## args
if [ $# -ne 4 ]
then
  echo "Usage: $0 /tmp/teste landsat_X PPPRRR YYYY-MM-DD"
  exit 1
fi

## args
ROOT_DIR=$1
IMAGE_DATASET=$2
IMAGE_PATHROW=$3
IMAGE_DATE=$4

# global variables
SANDBOX=$(pwd)
R_EXEC_DIR=$SANDBOX/src
SCRIPTS_DIR=$SANDBOX/scripts
ELEVATION_DIR=$R_EXEC_DIR/elevation
TMP_DIR_PATH=/tmp
MAX_TRIES=1

R_ALGORITHM_PREPROCESSING_VERSION=preprocessing.R

# folders
INPUTDOWNLOAD_DIR_PATH=$ROOT_DIR/inputdownloading
PREPROCESSING_DIR_PATH=$ROOT_DIR/preprocessing

# files
TIME_FILE=$PREPROCESSING_DIR_PATH/time

echo "Step 0. Time (start)"
touch $TIME_FILE
echo "START `date +"%s"`" >> $TIME_FILE

echo "Step 1. Capture MTL and station path"

files=($INPUTDOWNLOAD_DIR_PATH/*_MTL.txt)
for file in "${files[@]}"
do
  filename="${file##*/}"
  filenameWithoutExtension="${filename%_MTL*}"
  IMAGE_NAME="$filenameWithoutExtension"
done

IMAGE_MTL_PATH=$INPUTDOWNLOAD_DIR_PATH/$IMAGE_NAME"_MTL.txt"
IMAGE_STATION_FILE_PATH=$INPUTDOWNLOAD_DIR_PATH/$IMAGE_NAME"_station.csv"

echo "MTL path: $IMAGE_MTL_PATH"
echo "Station path: $IMAGE_STATION_FILE_PATH"

echo "Step 2. Creating dados.csv for image $IMAGE_NAME"

echo "File images;MTL;File Station Weather;Path Output" > $R_EXEC_DIR/dados.csv
echo "$INPUTDOWNLOAD_DIR_PATH;$IMAGE_MTL_PATH;$IMAGE_STATION_FILE_PATH;$PREPROCESSING_DIR_PATH" >> $R_EXEC_DIR/dados.csv

echo "Step 3. Starting CPU, disk and Memory collect..."

bash $SCRIPTS_DIR/collect-cpu-usage.sh $(pidof R) | tee $PREPROCESSING_DIR_PATH/$IMAGE_NAME"_cpu_usage.txt" > /dev/null &
bash $SCRIPTS_DIR/collect-memory-usage.sh $(pidof R) | tee $PREPROCESSING_DIR_PATH/$IMAGE_NAME"_mem_usage.txt" > /dev/null &
bash $SCRIPTS_DIR/collect-disk-usage.sh $(pidof R) | tee $PREPROCESSING_DIR_PATH/$IMAGE_NAME"_disk_usage.txt" > /dev/null &

echo "Step 4. Executing R script"

for i in `seq $MAX_TRIES`
  do
  
  rm -rf $PREPROCESSING_DIR_PATH/out.log $PREPROCESSING_DIR_PATH/error.log
  bash $SCRIPTS_DIR/executeRScript.sh $R_EXEC_DIR/$R_ALGORITHM_PREPROCESSING_VERSION $R_EXEC_DIR $TMP_DIR_PATH $PREPROCESSING_DIR_PATH
  PROCESS_OUTPUT=$?
  
  echo "executeRScript_process_output=$PROCESS_OUTPUT"
  if [ $PROCESS_OUTPUT -eq 0 ]
  then
    echo "Number of tries: $i"
    break
  elif [ $PROCESS_OUTPUT -eq 124 ] && [ $i -ge $MAX_TRIES ]
  then
    exit 124
  else
    if [ $i -ge $MAX_TRIES ]
    then
      echo "Number of tries: $i"
      exit 1
    fi
  fi
done

echo "Step 5. Killing collect CPU, disk and Memory scripts"

ps -ef | grep collect-cpu-usage.sh | grep -v grep | awk '{print $2}' | xargs kill
ps -ef | grep collect-memory-usage.sh | grep -v grep | awk '{print $2}' | xargs kill
ps -ef | grep collect-disk-usage.sh | grep -v grep | awk '{print $2}' | xargs kill

echo "Step 6. Moving dados.csv"
mv $R_EXEC_DIR/dados.csv $PREPROCESSING_DIR_PATH

echo "Step 7. Generate metadata"
bash $SANDBOX/generate_metadata.sh $IMAGE_DATASET $IMAGE_PATHROW $IMAGE_DATE $INPUTDOWNLOAD_DIR_PATH $PREPROCESSING_DIR_PATH

echo "Step 8. Time (finish)"
echo "START `date +"%s"`" >> $TIME_FILE

## Exit code
# exit code `0` indicates a successful execution. Any other number indicates failure.
