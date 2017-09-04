#! /bin/bash

#David Shean
#dshean@gmail.com

#This script will download and prepare global RGI glacier outlines
#212248 features in 15 regions

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

echo "Downlaiding RGI glacier inventory shp"
cd $DATADIR
pwd

rgi_zip_fn='rgi60.zip'
rgi_fn='rgi60'

#This zipfile is 36 GB
if [ ! -e $rgi_zip_fn ] ; then
    url='http://www.glims.org/RGI/rgi60_files/00_rgi60.zip'
    echo "Downloading $rgi_zip_fn"
    wget -O $rgi_zip_fn $url
fi

#Check MD5 sum
#95bd0486301073bbbe26988303fdaa1d  00_rgi60.zip
orig_zip_md5=95bd0486301073bbbe26988303fdaa1d

#Should check OS with uname
if [ -x "$(command -v md5)" ]; then
    #OS X
    md5=md5
elif [ -x "$(command -v md5sum)" ]; then
    #Linux/Unix
    md5=mdfsum
else
    md5=none
fi

if [ $md5 != "none" ] ; then
    zip_md5=$($md5 $rgi_zip_fn | awk '{print $NF}')
    if [ "$zip_md5" != "$orig_zip_md5" ] ; then
        echo "MD5 checksum failed for $rgi_zip_fn:"
        echo "Expected: $orig_zip_md5"
        echo "Downloaded: $zip_md5"
        exit 1
    fi
fi

if [ ! -d $rgi_fn ] ; then
    echo "Unzipping $rgi_zip_fn"
    #Note: need subdir here
    unzip $rgi_zip_fn -d $rgi_fn
fi

if [ ! -d $rgi_fn/regions ] ; then
    mkdir $rgi_fn/regions
    #Should check to see if GNU Parallel is installed, otherwise use for loop
    if [ -x "$(command -v parallel)" ]; then
        parallel "unzip {} -d $rgi_fn/regions" ::: $(ls $rgi_fn/*.zip | grep -v ^00_)
    else
        for i in $(ls $rgi_fn/*.zip | grep -v ^00_)
        do
            unzip $i -d $rgi_fn/regions
        done
    fi
fi

cd $rgi_fn/regions
mv 00_rgi*Regions* ../

#Merge regions for global shp
#if [ ! -e rgi60_merge.shp ] ; then
echo "Merging region shp for global shp"
ogr_merge.sh rgi60_merge.shp *_rgi60_*shp
#fi

echo $rgi_fn/regions/rgi60_merge.shp
