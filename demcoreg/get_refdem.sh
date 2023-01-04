#! /bin/bash

#Uses OpenTopography API to fetch reference DEM
#https://portal.opentopography.org/apidocs/#/Public/getGlobalDem

if [ "$#" -ne 7 ]; then
    echo "Usage: $(basename $0) minlon minlat maxlon maxlat out_fn demtype epsg_code"
    echo "Example for Copernicus 30 m DEM over Utqiagvik Alaska and UTM Zone 4N"
    echo "$(basename $0) -157.5 71 -155.5 71.75 COP30_utqiagvik.tif COP30 32604"
    exit 1
fi

#OT demo key (replace with user API key)
API_key="demoapikeyot2022"

minlon=$1
minlat=$2
maxlon=$3
maxlat=$4

#Output filename
out_fn=$5

#['SRTMGL3', 'SRTMGL1', 'SRTMGL1_E', 'AW3D30', 'AW3D30_E', 'SRTM15Plus', 'NASADEM', 'COP30', 'COP90']
#demtype="COP30"
#Check against list
demtype=$6

#Example EPSG code for UTM Zone 12N (should compute automatically)
#epsg=32612
epsg=$7

url="https://portal.opentopography.org/API/globaldem?API_Key=${API_key}&demtype=${demtype}&south=${minlat}&north=${maxlat}&west=${minlon}&east=${maxlon}&outputFormat=GTiff" 
echo $url
curl -X GET "$url" -H "accept: */*" --output $out_fn 

#Tile and compress
opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER -co NUM_THREADS=ALL_CPUS"
#opt="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER"
gdal_translate $opt $out_fn ${out_fn%.*}_lzw.tif 

#if [ $? -eq 0 ] ; then
#    mv ${out_fn%.*}_lzw.tif $out_fn
#fi

#Convert to ellipsoid height (uses ASP dem_geoid utility, update to corresponding GDAL/PROJ command)
dem_geoid --reverse-adjustment --geoid EGM2008 ${out_fn%.*}_lzw.tif

#Reproject to UTM
gdalwarp -overwrite $opt -r cubic -t_srs EPSG:$epsg ${out_fn%.*}_lzw-adj.tif ${out_fn%.*}_lzw-adj_proj.tif
