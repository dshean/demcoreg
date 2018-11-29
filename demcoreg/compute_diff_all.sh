#! /bin/bash

#Compute difference maps for all unique combinations of input files
#Input should be sorted chronologically

#Run as:
#compute_diff_all.sh *DEM_32m.tif
#compute_diff_all.sh $(ls -tr *00/dem*/*DEM_32m.tif | sort -n -t'/' -k 3)

narg=$#

set -- $@ 

echo -n > pairlist

for i; do
    shift
    for j; do
        echo $i $j >> pairlist
    done
done

npair=$(wc -l pairlist)

echo
echo "$narg inputs"
echo "$npair unique pairs"
echo

logdir=all_dh_log
if [ ! -d $logdir ] ; then
    mkdir $logdir
fi

#Compute raster differences
cat pairlist | parallel --group --colsep ' ' compute_diff.py {1} {2}

#Compute velocity pairs
#compute_diff_all.sh *hs.tif
#cat pairlist | parallel -j 12 --group --colsep ' ' 'vmap.py {1} {2}'

#Compute lagrangian elevation differences
#cat pairlist | parallel -j 12 --verbose --results $logdir --group --colsep ' ' 'vmap.py {1.}_hs.tif {2.}_hs.tif; lag_dz.py {1} {2} {1.}_hs__{2.}_hs_vmap/{1.}_hs__{2.}_hs_vmap-F.tif'
