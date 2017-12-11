#! /bin/bash

#Apply existing translation to all standard DEM products
#Input should be top-level directory (e.g., WV03_20160619_104001001D193C00_104001001E107400)

dir=$1
dem=$(find $dir -name '*-DEM_2m.tif')
dembase=${dem%%-*}
echo $dembase

if ls -t ${dembase}*align/*DEM.tif 1> /dev/null 2>&1 ; then 
    log=$(ls -t ${dembase}*align/*.log | head -1)
    if grep -q 'Translation vector' $log ; then
        echo
        if [ -e ${dembase}-DEM_32m.tif ] ; then 
            apply_dem_translation.py ${dembase}-DEM_32m.tif $log
        fi
        if [ -e ${dembase}-DEM_8m.tif ] ; then 
            apply_dem_translation.py ${dembase}-DEM_8m.tif $log
        fi
        if [ -e ${dembase}-DEM_2m.tif ] ; then 
            ln -vsf $(ls $(basename $dembase)*align/*DEM.tif) ${dembase}-DEM_2m_trans.tif
        fi
    fi
else 
    echo "Unable to find output from pc_align with expected directory structure"
fi
