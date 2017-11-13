#! /bin/bash

if [ -z "$DATADIR" ] ; then
    export DATADIR=$HOME/data
fi

echo "Downlaiding Natural Earth layers"
cd $DATADIR
mkdir NaturalEarth
cd NaturalEarth

pwd

echo "50M Cultural"
wget http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/cultural/50m_cultural.zip
unzip -d 50m_cultural 50m_cultural.zip

echo "10M Cultural"
wget http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/10m_cultural.zip
unzip -d 10m_cultural 10m_cultural.zip

echo "10M Physical"
wget http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/10m_physical.zip
unzip -d 10m_physical 10m_physical.zip
