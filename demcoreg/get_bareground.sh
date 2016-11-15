#! /bin/bash

#David Shean
#dshean@gmail.com

#This script will download global 30-m bare earth data 
#Needed for various masking applicaitons

#Uses GNU parallel for processing

gdal_opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES"

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

cd $DATADIR

#~2010 global bare ground, 30 m
#When unzipped, this is 64 GB!
#http://landcover.usgs.gov/glc/BareGroundDescriptionAndDownloads.php

be_zip_fn='bare2010.zip'
be_fn='bare2010'
be_vrt_fn='bare2010/bare2010.vrt'

#This zipfile is 36 GB
if [ ! -e $be_fn ] ; then
    url='http://edcintl.cr.usgs.gov/downloads/sciweb1/shared/gtc/downloads/bare2010.zip'
    echo "Downloading $be_zip_fn"
    wget $url $be_zip_fn
fi

if [ ! -d $be_fn ] ; then
    echo "Unzipping $be_zip_fn"
    unzip $be_zip_fn
fi

#These are 10x10-degree tiles at 30 m
#Original data are uncompressed, untiled
#Clean up
fn_list=$(ls bare2010/bare2010_v3/*bare2010_v3.tif)
c=$(wc -c $fn_list)
t=false
echo "Checking $c bare2010 files"
for fn in $fn_list
do
    if ! gdalinfo $fn | grep -q LZW ; then
        t=true
        break
    fi
done

if $t ; then
    echo "Cleaning up tiles"
    parallel "gdal_translate $gdal_opt {} {.}_lzw.tif; mv {.}_lzw.tif {}" ::: bare2010/bare2010_v3/*bare2010_v3.tif
fi

#Now create a vrt for all input data
#Set nodata as 255, as 0 is valid bare ground percentage
if [ ! -e $be_vrt_fn ] ; then
    echo "Building vrt of cleaned tiles"
    gdalbuildvrt -vrtnodata 255 $be_vrt_fn bare2010/bare2010_v3/*bare2010_v3.tif
fi
