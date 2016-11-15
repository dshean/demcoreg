#! /bin/bash

#David Shean
#dshean@gmail.com

#This script will download national land cover data
#Needed for various masking applicaitons

gdal_opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES"

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

cd $DATADIR

#2011 Land Use Land Cover (nlcd) grids, 30 m
#http://www.mrlc.gov/nlcd11_leg.php

nlcd_zip_fn='nlcd_2011_landcover_2011_edition_2014_10_10.zip'
nlcd_fn='nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img'

if [ ! -e $nlcd_zip_fn ] ; then
    echo "Downloading $nlcd_zip_fn"
    url='http://www.landfire.gov/bulk/downloadfile.php?TYPE=nlcd2011&FNAME=nlcd_2011_landcover_2011_edition_2014_10_10.zip'
    wget $url $nlcd_zip_fn
fi

if [ ! -e $nlcd_fn ] ; then
    echo "Unzipping $nlcd_zip_fn"
    unzip $nlcd_zip_fn
fi
