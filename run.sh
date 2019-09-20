#!/bin/bash

## This script downloads the input data for later processing.

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
MAX_TRIES=2

R_ALGORITHM_PREPROCESSING_VERSION=preprocessing.R

# folders
INPUTDOWNLOAD_DATA_DIR_PATH=$ROOT_DIR/inputdownload/data
PREPROCESSING_DIR_PATH=$ROOT_DIR/preprocessing
PREPROCESSING_DATA_DIR_PATH=$PREPROCESSING_DIR_PATH/data
PREPROCESSING_METADATA_DIR_PATH=$PREPROCESSING_DIR_PATH/metadata

echo "Step 1. Capture MTL and station path"

files=($INPUTDOWNLOAD_DATA_DIR_PATH/*_MTL.txt)
for file in "${files[@]}"
do
  filename="${file##*/}"
  filenameWithoutExtension="${filename%_MTL*}"
  IMAGE_NAME="$filenameWithoutExtension"
done

IMAGE_MTL_PATH=$INPUTDOWNLOAD_DATA_DIR_PATH/$IMAGE_NAME"_MTL.txt"
IMAGE_STATION_FILE_PATH=$INPUTDOWNLOAD_DATA_DIR_PATH/$IMAGE_NAME"_station.csv"

echo "MTL path: $IMAGE_MTL_PATH"
echo "Station path: $IMAGE_STATION_FILE_PATH"

echo "Step 2. Creating dados.csv for image $IMAGE_NAME"

echo "File images;MTL;File Station Weather;Path Output" > $R_EXEC_DIR/dados.csv
echo "$INPUTDOWNLOAD_DATA_DIR_PATH;$IMAGE_MTL_PATH;$IMAGE_STATION_FILE_PATH;$PREPROCESSING_DATA_DIR_PATH" >> $R_EXEC_DIR/dados.csv

echo "Step 3. Starting CPU, disk and Memory collect..."

bash $SCRIPTS_DIR/collect-cpu-usage.sh $(pidof R) | tee $PREPROCESSING_METADATA_DIR_PATH/$IMAGE_NAME"_cpu_usage.txt" > /dev/null &
bash $SCRIPTS_DIR/collect-memory-usage.sh $(pidof R) | tee $PREPROCESSING_METADATA_DIR_PATH/$IMAGE_NAME"_mem_usage.txt" > /dev/null &
bash $SCRIPTS_DIR/collect-disk-usage.sh $(pidof R) | tee $PREPROCESSING_METADATA_DIR_PATH/$IMAGE_NAME"_disk_usage.txt" > /dev/null &

echo "Step 4. Download raster elevation"
 
IFS='_'
read -ra SPLIT_IMAGE_NAME <<< "$IMAGE_NAME"
wget -P $ELEVATION_DIR https://github.com/ThiagoWhispher/elevation-raster/raw/master/${SPLIT_IMAGE_NAME[2]}.tif

echo "Step 5. Executing R script"

for i in `seq $MAX_TRIES`
  do
  
  rm -rf $TMP_DIR_PATH/* $PREPROCESSING_METADATA_DIR_PATH/out.log $PREPROCESSING_METADATA_DIR_PATH/error.log
  bash $SCRIPTS_DIR/executeRScript.sh $R_EXEC_DIR/$R_ALGORITHM_PREPROCESSING_VERSION $R_EXEC_DIR $TMP_DIR_PATH $PREPROCESSING_METADATA_DIR_PATH
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

echo "Step 6. Killing collect CPU, disk and Memory scripts"

ps -ef | grep collect-cpu-usage.sh | grep -v grep | awk '{print $2}' | xargs kill
ps -ef | grep collect-memory-usage.sh | grep -v grep | awk '{print $2}' | xargs kill
ps -ef | grep collect-disk-usage.sh | grep -v grep | awk '{print $2}' | xargs kill

echo "Step 7. Moving dados.csv"
mv $R_EXEC_DIR/dados.csv $PREPROCESSING_METADATA_DIR_PATH

echo "Step 8. Generate metadata"
bash $SANDBOX/generate_metadata.sh $INPUTDOWNLOAD_DATA_DIR_PATH $PREPROCESSING_DATA_DIR_PATH $PREPROCESSING_METADATA_DIR_PATH

## Exit code
# exit code `0` indicates a successful execution. Any other number indicates failure.
