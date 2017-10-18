#! /bin/bash

#Apply existing translation to all standard DEM products
#Input should be top-level directory (e.g., WV03_20160619_104001001D193C00_104001001E107400)

dir=$1
align=ned_align
dem=$(ls $dir/dem*/*-DEM_2m.tif)
dembase=${dem%%-*}

if ls -t $dir/dem*/*$align/*DEM.tif 1> /dev/null 2>&1 ; then 
    log=$(ls -t $dir/dem*/*$align/*.log | head -1)
    if grep -q 'Translation vector' $log ; then
        echo
        if [ -e ${dembase}-DEM_32m.tif ] ; then 
            apply_dem_translation.py ${dembase}-DEM_32m.tif $log
        fi
        if [ -e ${dembase}-DEM_8m.tif ] ; then 
            apply_dem_translation.py ${dembase}-DEM_8m.tif $log
        fi
        if [ -e ${dembase}-DEM_2m.tif ] ; then 
            ln -vsf $dir/dem*/*$align/*DEM.tif ${dembase}-DEM_2m_trans.tif
        fi
    fi
fi
