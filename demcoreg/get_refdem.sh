#! /bin/bash

#Uses OpenTopography API to fetch reference DEM
#https://portal.opentopography.org/apidocs/#/Public/getGlobalDem

#["SRTMGL3", "SRTMGL1", "SRTMGL1_E", "AW3D30", "AW3D30_E", "NASADEM", "COP30", "COP90"]
demtype="COP30"

minlon=$1
minlat=$2
maxlon=$3
maxlat=$4

out_fn=$5

api=demoapikeyot2022

url="https://portal.opentopography.org/API/globaldem?demtype=${demtype}&south=${minlat}&north=${maxlat}&west=${minlon}&east=${maxlon}&outputFormat=GTiff&API_Key=${api}" 
echo $url
curl -X GET "$url" -H "accept: */*" --output $out_fn 

#Tile and compress
opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER -co NUM_THREADS=ALL_CPUS"
#opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER"
gdal_translate $opt $out_fn ${out_fn%.*}_lzw.tif 

#if [ $? -eq 0 ] ; then
#    mv ${out_fn%.*}_lzw.tif $out_fn
#fi
