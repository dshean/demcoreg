#! /bin/bash

#This script will download national land cover data
#Needed for various masking applications
#Should be better than global bareground data for identifying exposed rock surfaces

# Usage:
# ./get_nlcd.sh year

set -e

# Input year
yr=$1

# Default is 2016
if [ -z "$yr" ] ; then
    echo "Using default NLCD year 2016"
    yr=2016
fi

gdal_opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES"

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

if [ ! -e $DATADIR ] ; then
    mkdir $DATADIR
fi

cd $DATADIR

#CONUS Land Cover (NLCD) grids, 30 m from https://www.mrlc.gov/data
url="https://s3-us-west-2.amazonaws.com/mrlc/NLCD_${yr}_Land_Cover_L48_20190424.zip"
nlcd_zip_fn="NLCD_${yr}_Land_Cover_L48_20190424.zip"
nlcd_fn="NLCD_${yr}_Land_Cover_L48_20190424.tif"
nlcd_dir=${nlcd_zip_fn%.*}

if [ ! -e $nlcd_zip_fn ] ; then
    echo "Downloading $nlcd_zip_fn"
    wget -O $nlcd_zip_fn $url
fi

if [ ! -e $nlcd_fn ] ; then
    echo "Unzipping $nlcd_zip_fn"
    unzip -d $nlcd_dir $nlcd_zip_fn
    #NOTE: Default OS X unzip fails to extract: skips *.ige due to error (need PK compat. v4.5 (can do v2.1))
    #Can use 7za: brew install p7zip
    #7za x $nlcd_zip_fn
    #But tar -z should be a safe cross-platform option.
    #tar -xzvf $nlcd_zip_fn
    #Can save a lot of disk space compressing the original ArcGrid img and deleting overviews
    #Input is ~17 GB, output is 1.1 GB
    echo "Creating compressed GTiff copy"
    echo gdal_translate $gdal_opt $nlcd_dir/${nlcd_fn%.*}.img $DATADIR/temp.tif
    gdal_translate $gdal_opt $nlcd_dir/${nlcd_fn%.*}.img $DATADIR/temp.tif
    echo "Removing metadata directory and contents"
    rm -rfv $nlcd_dir
    mv -v $DATADIR/temp.tif $nlcd_fn
else
    echo "Found existing ${nlcd_fn}!"
fi
