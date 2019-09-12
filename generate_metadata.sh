#!/bin/bash

## Checking args
if [ $# -ne 3 ]
then
  echo "Usage: $0 input_path output_path metadata_path"
  exit 1
fi

## Capture args
INPUT_DIR_PATH=$1
OUTPUT_DIR_PATH=$2
METADATA_DIR_PATH=$3

METADATA_FILE_PATH=$METADATA_DIR_PATH/metadata.txt
rm -rf $METADATA_FILE_PATH
touch $METADATA_FILE_PATH

MTL_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_MTL.txt")
LANDSAT_SCENE_ID=$(cat $MTL_INPUT_FILE_PATH | grep LANDSAT_SCENE_ID)
TYPE_L8_INPUT_FILE=$(echo $LANDSAT_SCENE_ID | grep LC8)
TYPE_L7_INPUT_FILE=$(echo $LANDSAT_SCENE_ID | grep LE7)
TYPE_L5_INPUT_FILE=$(echo $LANDSAT_SCENE_ID | grep LT5)

B1_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B1.TIF")
B2_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B2.TIF")
B3_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B3.TIF")
B4_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B4.TIF")
B5_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B5.TIF")
B6_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B6.TIF")
B6_VCID_1_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B6_VCID_1.TIF")
B6_VCID_2_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B6_VCID_2.TIF")
B7_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B7.TIF")
B8_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B8.TIF")
B9_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B9.TIF")
B10_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B10.TIF")
B11_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_B11.TIF")
BQA_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_BQA.TIF")
STATION_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_station.csv")

RASTER_ELEVATION_OUTPUT_FILE_PATH=$(find $OUTPUT_DIR_PATH -iname "elevation.tif")
EVI_OUTPUT_FILE_PATH=$(find $OUTPUT_DIR_PATH -iname "*_EVI.nc")
LAI_OUTPUT_FILE_PATH=$(find $OUTPUT_DIR_PATH -iname "*_LAI.nc")
SAVI_OUTPUT_FILE_PATH=$(find $OUTPUT_DIR_PATH -iname "*_SAVI.nc")
NDVI_OUTPUT_FILE_PATH=$(find $OUTPUT_DIR_PATH -iname "*_NDVI.nc")
LSA_OUTPUT_FILE_PATH=$(find $OUTPUT_DIR_PATH -iname "*_alb.nc")
LST_OUTPUT_FILE_PATH=$(find $OUTPUT_DIR_PATH -iname "*_TS.nc")
RN_OUTPUT_FILE_PATH=$(find $OUTPUT_DIR_PATH -iname "*_Rn.nc")
G_OUTPUT_FILE_PATH=$(find $OUTPUT_DIR_PATH -iname "*_G.nc")

CURRENT_DATE=$(date)

echo "#Preprocessing (default) Implementation Metadata" >> $METADATA_FILE_PATH
echo "$CURRENT_DATE # Date" >> $METADATA_FILE_PATH

echo "INPUT FILES" >> $METADATA_FILE_PATH >> $METADATA_FILE_PATH
echo "$B1_INPUT_FILE_PATH # Band 1 from image" >> $METADATA_FILE_PATH
echo "$B2_INPUT_FILE_PATH # Band 2 from image" >> $METADATA_FILE_PATH
echo "$B3_INPUT_FILE_PATH # Band 3 from image" >> $METADATA_FILE_PATH
echo "$B4_INPUT_FILE_PATH # Band 4 from image" >> $METADATA_FILE_PATH
echo "$B5_INPUT_FILE_PATH # Band 5 from image" >> $METADATA_FILE_PATH
if [ "$TYPE_L5_INPUT_FILE" ]
then
  echo "$B6_INPUT_FILE_PATH # Band 6 from image" >> $METADATA_FILE_PATH
fi
if [ "$TYPE_L7_INPUT_FILE" ]
then
  echo "$B6_VCID_1_INPUT_FILE_PATH # Band 6 VCID 1 from image" >> $METADATA_FILE_PATH
  echo "$B6_VCID_2_INPUT_FILE_PATH # Band 6 VCID 2 from image" >> $METADATA_FILE_PATH
fi
if [ "$TYPE_L8_INPUT_FILE" ]
then
  echo "$B6_INPUT_FILE_PATH # Band 6 from image" >> $METADATA_FILE_PATH
fi

echo "$B7_INPUT_FILE_PATH # Band 7 from image" >> $METADATA_FILE_PATH

if [ "$TYPE_L7_INPUT_FILE" ]
then
  echo "$B8_INPUT_FILE_PATH # Band 8 from image" >> $METADATA_FILE_PATH
fi
if [ "$TYPE_L8_INPUT_FILE" ]
then
  echo "$B8_INPUT_FILE_PATH # Band 8 from image" >> $METADATA_FILE_PATH
  echo "$B9_INPUT_FILE_PATH # Band 9 from image" >> $METADATA_FILE_PATH
  echo "$B10_INPUT_FILE_PATH # Band 10 from image" >> $METADATA_FILE_PATH
  echo "$B11_INPUT_FILE_PATH # Band 11 from image" >> $METADATA_FILE_PATH
fi
echo "$BQA_INPUT_FILE_PATH # Band QA from image" >> $METADATA_FILE_PATH
echo "$MTL_INPUT_FILE_PATH # MTL from image" >> $METADATA_FILE_PATH
echo "$STATION_INPUT_FILE_PATH # Station data from image" >> $METADATA_FILE_PATH

echo "OUTPUT FILES" >> $METADATA_FILE_PATH
echo "$RASTER_ELEVATION_OUTPUT_FILE_PATH # Preprocessed Raster elevation data file path" >> $METADATA_FILE_PATH
echo "$EVI_OUTPUT_FILE_PATH # Enhanced Vegetation Index data file path" >> $METADATA_FILE_PATH
echo "$LAI_OUTPUT_FILE_PATH # Leaf Area Index data file path" >> $METADATA_FILE_PATH
echo "$SAVI_OUTPUT_FILE_PATH # Soil Adjusted Vegetation Index data file path" >> $METADATA_FILE_PATH
echo "$NDVI_OUTPUT_FILE_PATH # Normalized Different Vegetation Index data file path" >> $METADATA_FILE_PATH
echo "$LSA_OUTPUT_FILE_PATH # Land Surface Albedo data file path" >> $METADATA_FILE_PATH
echo "$LST_OUTPUT_FILE_PATH # Land Surface Temperature data file path" >> $METADATA_FILE_PATH
echo "$RN_OUTPUT_FILE_PATH # Net Radiation Balance data file path" >> $METADATA_FILE_PATH
echo "$G_OUTPUT_FILE_PATH # Ground Heat Flux data file path" >> $METADATA_FILE_PATH
