#! /bin/bash

#David Shean
#dshean@gmail.com

#This script will download national land cover data
#Needed for various masking applicaitons
#Should be better than global bareground data for identifying exposed rock surfaces

gdal_opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES"

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

cd $DATADIR

#2011 Land Use Land Cover (nlcd) grids, 30 m
#http://www.mrlc.gov/nlcd11_leg.php

nlcd_zip_fn='nlcd_2011_landcover_2011_edition_2014_10_10.zip'
nlcd_fn='nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.tif'

if [ ! -e $nlcd_zip_fn ] ; then
    echo "Downloading $nlcd_zip_fn"
    url='http://www.landfire.gov/bulk/downloadfile.php?TYPE=nlcd2011&FNAME=nlcd_2011_landcover_2011_edition_2014_10_10.zip'
    wget -O $nlcd_zip_fn $url 
fi

if [ ! -e $nlcd_fn ] ; then
    echo "Unzipping $nlcd_zip_fn"
    #unzip $nlcd_zip_fn
    #NOTE: Default OS X unzip fails to extract: skips *.ige due to error (need PK compat. v4.5 (can do v2.1))
    #Can use 7za: brew install p7zip
    #7za x $nlcd_zip_fn
    #But tar -z should be a safe cross-platform option
    tar -xzvf $nlcd_zip_fn
    #Can save a lot of disk space compressing the original ArcGrid img and deleting overviews
    #Input is ~17 GB, output is 1.1 GB
    echo "Creating compressed GTiff copy"
    echo gdal_translate $gdal_opt ${nlcd_fn%.*}.img $DATADIR/temp.tif
    gdal_translate $gdal_opt ${nlcd_fn%.*}.img $DATADIR/temp.tif
    echo "Removing $(dirname $nlcd_fn) contents"
    rm -rv $(dirname $nlcd_fn)/*
    mv -v $DATADIR/temp.tif $nlcd_fn
else 
    echo "Found existing ${nlcd_fn}!"
fi
