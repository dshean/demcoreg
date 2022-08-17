#! /bin/bash

#David Shean
#dshean@gmail.com
# Updated by JMH, 8.16.22

#This script will download global 30-m bare earth data 
#Needed for various masking applicaitons

#Uses GNU parallel for processing

# For more information, visit https://glad.umd.edu/dataset/global-2010-bare-ground-30-m
# Percentage of bare ground cover per output grid cell is encoded as integer values (1-100). 
# Files are unsigned 8-bit data with spatial resolution of 0.00025° per pixel (approximately 30 meters per pixel at the equator). 
# Data coverage is from 80° north to 60° south. It is divided into tiles of 10 x 10 degrees. 
# The naming convention per tile refers to the latitude and longitude value of the upper left corner of the tile. 
# Tiles over ocean area are included for completeness.

# Updated links
# https://glad.umd.edu/Potapov/Bare_2010/[YY][N..S]_[XXX][E..W].tif
# Y ranges between 00–80 N
# Y ranges between 10–50 S
# X ranges between 010–180 W
# X ranges between 000–170 E
# increments by 10˚ in both lat and lon

gdal_opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES"

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data/bare
fi

# Make datadir if it doesn't exist
if [ ! -e $DATADIR ] ; then
    mkdir $DATADIR
fi

cd $DATADIR

# Final processed vrt of all tiles
be_vrt_fn=$DATADIR/bare2010.vrt

# If vrt DNE, download tiles with/without overwriting (default is overwrite)
overwrite=false
if [ ! -e $be_vrt_fn ] ; then
    # This puts things into Potapov/Bare_2010 dir
    # 51GB for 504 .tifs
    echo "Downloading tiles"
    time parallel -v wget -r -np -nH -nc --no-check-certificate -A {} https://glad.umd.edu/Potapov/Bare_2010/ ::: "[0-9][0-9][NS]_[0-9][0-9][0-9][EW].tif" &

    # move files to datadir
    echo
    echo "Moving tifs ../Potapov/Bare_2010/*tif $DATADIR"
    mv ../Potapov/Bare_2010/*tif $DATADIR
    echo
    echo "Removing empty dirs"
    rm -r ../Potapov
    #Now create a vrt for all input data
    #Set nodata as 255, as 0 is valid bare ground percentage
    echo "Building vrt of cleaned tiles: $be_vrt_fn"
    gdalbuildvrt -vrtnodata 255 $be_vrt_fn $DATADIR/*.tif
else
    echo
    echo $be_vrt_fn exists 
fi

# #~2010 global bare ground, 30 m
# #When unzipped, this is 64 GB!
# #http://landcover.usgs.gov/glc/BareGroundDescriptionAndDownloads.php

# be_zip_fn='bare2010.zip'
# be_fn='bare2010'
# be_vrt_fn='bare2010/bare2010.vrt'

# #This zipfile is 36 GB
# if [ ! -e $be_fn ] ; then
#     url='http://edcintl.cr.usgs.gov/downloads/sciweb1/shared/gtc/downloads/bare2010.zip'
#     echo "Downloading $be_zip_fn"
#     wget -O $be_zip_fn $url 
# fi

# if [ ! -d $be_fn ] ; then
#     echo "Unzipping $be_zip_fn"
#     unzip $be_zip_fn
# fi

# #These are 10x10-degree tiles at 30 m
# #Original data are uncompressed, untiled
# #Clean up
# fn_list=$(ls $DATADIR/*.tif)
# c=$(wc -c $fn_list)
# t=false
# echo "Checking $c bare2010 files"
# for fn in $fn_list
# do
#     if ! gdalinfo $fn | grep -q LZW ; then
#         t=true
#         break
#     fi
# done

# if $t ; then
#     echo "Cleaning up tiles"
#     parallel "gdal_translate $gdal_opt {} {.}_lzw.tif; mv {.}_lzw.tif {}" ::: $DATADIR/*.tif
# fi

# #Now create a vrt for all input data
# #Set nodata as 255, as 0 is valid bare ground percentage
# if [ ! -e $be_vrt_fn ] ; then
#     echo "Building vrt of cleaned tiles"
#     gdalbuildvrt -vrtnodata 255 $be_vrt_fn $DATADIR/*.tif
# fi
