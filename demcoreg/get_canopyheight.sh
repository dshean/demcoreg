#! /bin/bash

# This script will download the global forest canopy height GEDI-Landsat data product detailed in https://doi.org/10.1016/j.rse.2020.112165
# May be better than NLCD and global bareground data for identifying canopy cover and exposed rock surfaces

# 7 regions available, see https://glad.umd.edu/dataset/gedi/ for more info and access to GEE app
# AUS:        Australasia region         [1.5 GB]
# NAFR:       North Africa region        [3.1 GB] (includes Europe and parts of Middle East)
# NAM:        North America region       [5.7 GB]
# NASIA:      North Asia region          [4.1 GB]
# SAFR:       South Africa region        [5.7 GB]
# SAM:        South America region       [6.9 GB] (includes part of Central America) 
# SASIA:      South Asia region          [4.7 GB]

# Usage
# ./get_canopyheight.sh [ region ] [ output data directory] 

###################################################################################################
# Help
Help()
{
   # Display Help
   echo "Script downloader for 2019 Global Forest Canopy Height GEDI-Landsat 30 m products from"
   echo "https://glad.umd.edu/dataset/gedi/"
   echo
   echo "Usage:"
   echo "get_canopyheight.sh [ REGIONCODE ] [ output data directory ]"
   echo
   echo "Defaults"
   echo "       region      NAM (North America)"
   echo "       datadir     $HOME/data/canopyheight"
   echo
   echo "Available (continental) regions"
   echo "       AUS:        Australasia region         [1.5 GB]"
   echo "       NAFR:       North Africa region        [3.1 GB]"
   echo "       NAM:        North America region       [5.7 GB]"
   echo "       NASIA:      North Asia region          [4.1 GB]"
   echo "       SAFR:       South Africa region        [5.7 GB]"
   echo "       SAM:        South America region       [6.9 GB]"
   echo "       SASIA:      South Asia region          [4.7 GB]"
   echo 
   echo "options:"
   echo "h     Print this Help."
   echo
}
###################################################################################################

# Get the options
while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done

###################################################################################################
set -e

# Input region
region=$1
DATADIR=$2

# Default is North America 'NAM'
if [ -z "$region" ] ; then
    echo "Using default region: North America ['NAM']"
    region='NAM'
fi

gdal_opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES"

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data/canopyheight
fi

if [ ! -e $DATADIR ] ; then
    mkdir $DATADIR
fi

cd $DATADIR

# Global 2019 grids, 30m from https://glad.geog.umd.edu/Potapov/Forest_height_2019/
url="https://glad.geog.umd.edu/Potapov/Forest_height_2019/Forest_height_2019_${region}.tif"
canopyheight_fn="Forest_height_2019_${region}.tif"

if [ ! -e $canopyheight_fn ] ; then
    echo "Downloading $canopyheight_fn"
    wget -O $canopyheight_fn $url
else
    echo "Found existing ${canopyheight_fn}"
fi

# Haven't tested this yet, checking to see if the bits are already as compressed as can be
if ! gdalinfo $canopyheight_fn | grep -q LZW ; then
    echo "Compressing $canopyheight_fn"
    gdal_translate -co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER $canopyheight_fn temp.tif
    if [ $? -eq 0 ] ; then
        mv temp.tif $canopyheight_fn
    fi
fi
###################################################################################################