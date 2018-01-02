#! /bin/bash

#http://www.hydrosheds.org/

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

echo "Downlaiding HydroBASINS and HydroSHEDS layers"
cd $DATADIR
mkdir HydroBASINS
cd HydroBASINS

pwd

#The HydroBASINS products require login, staged level 1-12 lake products on Google Drive for now
#North America
if [ ! -d hybas_lake_na_lev01-12_v1c ] ; then
    wget https://drive.google.com/open?id=1MgblI_0RepPOtFbGmSpVISIJ_6oas7Db
    unzip hybas_lake_na_lev01-12_v1c.zip
fi

#Asia
if [ ! -d hybas_lake_as_lev01-12_v1c ] ; then
    wget https://drive.google.com/open?id=1bhtiWwrYNhn58ClPpjNj6yFaxshOkb8f
    unzip hybas_lake_as_lev01-12_v1c.zip
fi

#HydroSHEDS lake polygons
if [ ! -d HydroLAKES_polys_v10_shp ] ; then
    wget https://97dc600d3ccc765f840c-d5a4231de41cd7a15e06ac00b0bcc552.ssl.cf5.rackcdn.com/HydroLAKES_polys_v10_shp.zip
    unzip HydroLAKES_polys_v10_shp.zip
fi

#HydroSHEDS lake pour points
if [ ! -d HydroLAKES_points_v10_shp ] ; then
    wget https://97dc600d3ccc765f840c-d5a4231de41cd7a15e06ac00b0bcc552.ssl.cf5.rackcdn.com/HydroLAKES_points_v10_shp.zip
    unzip HydroLAKES_points_v10_shp.zip
fi
