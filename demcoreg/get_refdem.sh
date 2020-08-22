#! /bin/bash

#Uses OpenTopography API to fetch reference DEM
https://portal.opentopography.org/apidocs/#/Public/getGlobalDem

#["SRTMGL3", "SRTMGL1", "SRTMGL1_E", "AW3D30", "AW3D30_E"]
demtype="SRTMGL1_E"

minlon=$1
minlat=$2
maxlon=$3
maxlat=$4

out_fn=$5

url="https://portal.opentopography.org/API/globaldem?demtype=${demtype}&south=${minlat}&north=${maxlat}&west=${minlon}&east=${maxlon}&outputFormat=GTiff" 
curl -X GET "$url" -H "accept: */*" --output $out_fn 

#Tile and compress
opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER -co NUM_THREADS=ALL_CPUS"
gdal_translate $opt $out_fn /tmp/$$.tif

if [ $? -eq 0 ] ; then
    mv -f /tmp/$$.tif $out_fn
fi
