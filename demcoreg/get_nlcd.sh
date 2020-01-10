#! /bin/bash

#David Shean
#dshean@gmail.com

#This script will download national land cover data
#Needed for various masking applications
#Should be better than global bareground data for identifying exposed rock surfaces

gdal_opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES"

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

if [ ! -e $DATADIR ] ; then
    mkdir $DATADIR
fi

cd $DATADIR

#2016 Land Cover (NLCD) grids, 30 m
url='https://s3-us-west-2.amazonaws.com/mrlc/NLCD_2016_Land_Cover_L48_20190424.zip'
nlcd_zip_fn='NLCD_2016_Land_Cover_L48_20190424.zip'
nlcd_fn='NLCD_2016_Land_Cover_L48_20190424.tif'

if [ ! -e $nlcd_zip_fn ] ; then
    echo "Downloading $nlcd_zip_fn"
    wget -O $nlcd_zip_fn $url 
fi

if [ ! -e $nlcd_fn ] ; then
    echo "Unzipping $nlcd_zip_fn"
    unzip $nlcd_zip_fn
    #NOTE: Default OS X unzip fails to extract: skips *.ige due to error (need PK compat. v4.5 (can do v2.1))
    #Didn't run into any problems with UnZip 6.0 on OS X or Linux (-JMH)
    #Can use 7za: brew install p7zip
    #7za x $nlcd_zip_fn
    #But tar -z should be a safe cross-platform option.
    #Did run into errors with tar on Linux with stdin and multiple inputs
    #tar -xzvf $nlcd_zip_fn
    #Can save a lot of disk space compressing the original ArcGrid img and deleting overviews
    #Input is ~17 GB, output is 1.1 GB
    echo "Creating compressed GTiff copy"
    echo gdal_translate $gdal_opt ${nlcd_fn%.*}.img $DATADIR/temp.tif
    gdal_translate $gdal_opt ${nlcd_fn%.*}.img $DATADIR/temp.tif
    echo "Removing NLCD2016_spatial_metadata/ and contents"
    rm -rfv NLCD2016_spatial_metadata/*
    mv -v $DATADIR/temp.tif $nlcd_fn
else 
    echo "Found existing ${nlcd_fn}!"
fi
