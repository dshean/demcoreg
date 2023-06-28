#! /bin/bash

#This script will download national land cover data
#Needed for various masking applications
#Should be better than global bareground data for identifying exposed rock surfaces

# Usage:
# ./get_nlcd.sh [ region ] [ year ] [ ouput dat directory]

###################################################################################################
# Help
Help()
{
   # Display Help
   echo "Script downloader for NLCD 30 m products from"
   echo "https://www.mrlc.gov/data"
   echo
   echo "Usage:"
   echo "./get_nlcd.sh [ REGIONCODE ] [ year ] [ output data directory ]"
   echo
   echo "Defaults"
   echo "       region      CONUS (Conterminous US)"
   echo "       yr          2019"
   echo "       datadir     $HOME/data/"
   echo
   echo "Available regions"
   echo "       CONUS:        Lower 48 US states"
   echo "       AK:           Alaska"
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

echo ; echo Running: $0 $@ ; echo
# Input region
region=$1
yr=$2
DATADIR=$3

# Default is CONUS
if [ -z "$region" ] ; then
    echo "Using default region: 'CONUS'"
    region=CONUS
else
    echo "Using input region: $region"
fi

# Default is 2019
if [ -z "$yr" ] ; then
    echo "Using default NLCD year 2016" ; echo
    yr=2019
else
    echo "Using input year: $yr" ; echo
fi

gdal_opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES"

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

if [ ! -e $DATADIR ] ; then
    mkdir $DATADIR
fi

cd $DATADIR

if [ $region == "CONUS" ] ; then
    #CONUS Land Cover (NLCD) grids, 30 m from https://www.mrlc.gov/data
    url="https://s3-us-west-2.amazonaws.com/mrlc/nlcd_${yr}_land_cover_l48_20210604.zip"
elif [ $region == "AK" ] ; then
    #AK Land Cover (NLCD) grids, 30 m from https://www.mrlc.gov/data
    url="https://s3-us-west-2.amazonaws.com/mrlc/NLCD_${yr}_Land_Cover_AK_20200724.zip"
else
    echo "Region code not recognized: ${region}. Please enter 'CONUS' or 'AK'"
fi

nlcd_zip_fn=$(basename ${url})
nlcd_fn="${nlcd_zip_fn%.*}.tif"
nlcd_dir=${nlcd_zip_fn%.*}

if [[ ! -f $nlcd_zip_fn && ! -f $nlcd_fn ]] ; then
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
