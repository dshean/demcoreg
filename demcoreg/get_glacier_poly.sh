#! /bin/bash

#David Shean
#dshean@gmail.com

#This script will download global RGI glacier outlines 
#Needed for various masking applicaitons

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

cd $DATADIR

rgi_zip_fn='rgi50.zip'
rgi_fn='rgi50'

#This zipfile is 36 GB
if [ ! -e $rgi_zip_fn ] ; then
    url='http://www.glims.org/RGI/rgi50_files/rgi50.zip'
    echo "Downloading $rgi_zip_fn"
    wget $url $rgi_zip_fn
fi

if [ ! -d $rgi_fn ] ; then
    echo "Unzipping $rgi_zip_fn"
    #Note: need subdir here
    unzip $rgi_zip_fn -d $rgi_fn
fi

parallel "unzip {} -d $rgi_fn/regions" ::: $(ls $rgi_fn/*.zip | grep -v ^00_)

cd $rgi_fn/regions
ogr_merge.sh rgi50_merge.shp *_rgi50_*shp
