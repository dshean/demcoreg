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

rgi_zip_fn='rgi50.zip'
rgi_fn='rgi50'

#This zipfile is 36 GB
if [ ! -e $rgi_zip_fn ] ; then
    url='http://www.glims.org/RGI/rgi50_files/rgi50.zip'
    echo "Downloading $rgi_zip_fn"
    wget -O $rgi_zip_fn $url
fi

if [ ! -d $rgi_fn ] ; then
    echo "Unzipping $rgi_zip_fn"
    #Note: need subdir here
    unzip $rgi_zip_fn -d $rgi_fn
fi

if [ ! -d $rgi_fn/regions ] ; then
    mkdir $rgi_fn/regions
    #Should check to see if GNU Parallel is installed, otherwise use for loop
    parallel "unzip {} -d $rgi_fn/regions" ::: $(ls $rgi_fn/*.zip | grep -v ^00_)
fi

cd $rgi_fn/regions
mv 00_rgi*Regions* ../

#Merge regions for global shp
#if [ ! -e rgi50_merge.shp ] ; then
echo "Merging region shp for global shp"
ogr_merge.sh rgi50_merge.shp *_rgi50_*shp
#fi

echo $rgi_fn/regions/rgi50_merge.shp
