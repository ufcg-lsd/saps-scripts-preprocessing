#!/bin/bash

## This script pre-processes the data downloaded by the InputDownload phase

## Checking args
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

# folders
PREPROCESSING_DIR_PATH=$ROOT_DIR/preprocessing

# repo with results
REPO=http://www2.lsd.ufcg.edu.br/~thiagoyes/saps/nop-download-files/preprocessing
 
cd $PREPROCESSING_DIR_PATH

# download files
wget $REPO/dados.csv
wget $REPO/error.log
wget $REPO/out.log
wget $REPO/stage.metadata
wget $REPO/elevation.tif
wget $REPO/LC82150652015174LGN00_alb.nc
wget $REPO/LC82150652015174LGN00_EVI.nc
wget $REPO/LC82150652015174LGN00_G.nc
wget $REPO/LC82150652015174LGN00_LAI.nc
wget $REPO/LC82150652015174LGN00_NDVI.nc
wget $REPO/LC82150652015174LGN00_Rn.nc
wget $REPO/LC82150652015174LGN00_SAVI.nc
wget $REPO/LC82150652015174LGN00_TS.nc
wget $REPO/LC82150652015174LGN00_cpu_usage.txt
wget $REPO/LC82150652015174LGN00_disk_usage.txt
wget $REPO/LC82150652015174LGN00_mem_free.txt
wget $REPO/LC82150652015174LGN00_mem_usage.txt
wget $REPO/time

exit 0

## Exit code
# exit code `0` indicates a successful execution. Any other number indicates failure.
# In particular, `3` indicates that a Landsat image was not found for the given paramenters.
