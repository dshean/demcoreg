#! /bin/bash
set -e

#David Shean
#dshean@gmail.com

#This script will download and prepare global RGI glacier outlines
#212248 features in 15 regions

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

if [ ! -d $DATADIR ] ; then
    mkdir -p $DATADIR
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
    zip_md5=$(md5 $rgi_zip_fn | awk '{print $NF}')
elif [ -x "$(command -v md5sum)" ]; then
    #Linux/Unix
    zip_md5=$(md5sum $rgi_zip_fn | awk '{print $1}')
else
    zip_md5=none
fi

if [ $zip_md5 != "none" ] ; then
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
ogrmerge.py -single -o rgi60_merge.shp *_rgi60_*shp
#fi

function clipnproj() {
    site=$1
    wkt="$2"
    proj="$3"
    if [ ! -e rgi60_merge_${site}_aea.geojson ] ; then 
        echo
        echo "Clipping to site: $site"
        echo "$wkt"
        echo "$proj"
        ogr2ogr -progress -overwrite -clipsrc "$wkt" rgi60_merge_${site}.shp rgi60_merge.shp
        ogr2ogr -progress -overwrite -t_srs "$proj" rgi60_merge_${site}_aea.shp rgi60_merge_${site}.shp
        #ogr2ogr -progress -f 'GeoJSON' rgi60_merge_${site}_aea.geojson rgi60_merge_${site}_aea.shp
    fi
}

site=CONUS
wkt='POLYGON ((-125 49, -104 49, -104 32, -125 32, -125 49))'
proj='+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
clipnproj $site "$wkt" "$proj"

site=HMA
wkt='POLYGON ((66 47, 106 47, 106 25, 66 25, 66 47))'
proj='+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
clipnproj $site "$wkt" "$proj"

echo $rgi_fn/regions/rgi60_merge.shp
