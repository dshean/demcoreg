#! /bin/bash

#David Shean
#dshean@gmail.com

#This script will download global 1-km landcover daata: AVRR GLCF
#http://www.landcover.org/data/landcover/data.shtml

gdal_opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES"

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

cd $DATADIR
if [ ! -d glcf ] ; then 
    mkdir glcf
fi
cd glcf

glcf_zip_fn='AVHRR_1km_LANDCOVER_1981_1994.GLOBAL.tif.gz'
glcf_tif_fn=${glcf_zip_fn%.*}

if [ ! -e $glcf_tif_fn ] ; then
    if [ ! -e $glcf_zip_fn ] ; then
        url='ftp://ftp.glcf.umd.edu/glcf/Global_Land_Cover/Global/1km/AVHRR_1km_LANDCOVER_1981_1994.GLOBAL.tif.gz'
        echo "Downloading $glcf_zip_fn"
        wget -O $glcf_zip_fn $url 
    fi
    if [ ! -e $glcf_tif_fn ] ; then
        echo "Expanding $glcf_zip_fn"
        gunzip $glcf_zip_fn
    fi
fi

#Compress with lossless LZW
#Reduces file size from ~620 MB to ~25 MB
if ! gdalinfo $glcf_tif_fn | grep -q LZW ; then
    echo "Compressing $glcf_tif_fn"
    gdal_translate -co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER $glcf_tif_fn temp.tif
    if [ $? -eq 0 ] ; then
        mv temp.tif $glcf_tif_fn
    fi
fi
