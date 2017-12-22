#! /bin/bash

#David Shean
#dshean@gmail.com

#This will co-register a large batch of DEMs

#Run tpfe1 or bro node
#qsub -I -q devel -lselect=1:model=bro,walltime=2:00:00
#cd dir

#Clean up all previous runs
clean=false

#Some earlier r100 xml files lack the necessary tags for toa.  Need to copy original xml files from lfe
#ssh lfe
#cd ~/hma
#nbdir=/nobackupp8/deshean/hma/sites/khumbu/rerun
#nbdir=/nobackupp8/deshean/conus_combined/sites/rainier/rerun
#pushd $nbdir
#list_dir=$(ls -Sr *00/dem*/*-DEM_32m.tif | awk -F'/' '{print $1}')
#pushd
#NOTE: these need to be tightened, can unnecessarily transfer full directory again
#for i in $list_dir; do d=$(find . -name $i); shiftc -R $d/*xml $nbdir/$i/ ; done
#for i in $list_dir; do d=$(find conus*/ -name $i); shiftc -R $d/*xml $nbdir/$i/ ; done

#Compute top-of-atmosphere reflectance 
#Uses toa.sh, toa.py, and dglib from https://github.com/dshean/dgtools 
#parallel --jobs 16 --delay 2 'toa.sh {}' ::: $list_dir

#Copy DEMs and toa to subdir
#mkdir dem_coreg
#shiftc -L -r -d --include='-DEM_2m.tif' --include='-DEM_32m.tif' --include='-DEM_8m.tif' --include 'toa.tif' *00 dem_coreg/
##rsync -aLmv --include='*-DEM_*.tif' --include '*toa.tif' --include='*/' --exclude='*' *00 dem_coreg/
##rsync -aLmv --include='*-DEM_2m.tif' --include='*-DEM_32m.tif' --include='*-DEM_8m.tif' --include '*toa.tif' --include='*/' --exclude='*' *00 dem_coreg/
#cd dem_coreg

#NOTE: Careful with QB DEMs here, need to increase max offset

#Default is for all WV/GE/QB subdir
list_dir=$(ls -Sr *00/dem*/*-DEM_32m.tif | awk -F'/' '{print $1}')
list_32m=$(ls -Sr *00/dem*/*-DEM_32m.tif)
list_32m_done=$(ls -Sr *00/dem*/*-DEM_32m_trans.tif)

njobs=10

echo "Creating done and todo lists"
#Use existing lists to determine which round we are on
if [ ! -n "$n" ] ; then
    n=$(ls coreg_round*_done.txt | cut -c 12)
    if [ -z "$n" ] ; then
        n=1
    else
        n=$((n+1))
    fi
fi

if $clean ; then
    rm coreg_round*_*.txt
    n=1
fi

done_fn=coreg_round${n}_done.txt
todo_fn=coreg_round${n}_todo.txt
echo -n > $done_fn
echo -n > $todo_fn
for i in $list_32m
do
    if echo $list_32m_done | grep -q ${i%.*} ; then
        echo $i >> $done_fn
    else 
        echo $i >> $todo_fn
    fi
done

ndone=$(cat $done_fn | wc -l)
ntodo=$(cat $todo_fn | wc -l)
echo
echo "$ndone done"
echo "$ntodo todo"
echo

list_32m=$(cat ${todo_fn})
#list_2m=$(ls -Sr *00/dem*/*-DEM_2m.tif)
list_2m=$(echo $list_32m | sed 's/-DEM_32m.tif/-DEM_2m.tif/g')

#Clean up existing masks
if $clean ; then 
    rm -v *00/*/*-DEM_32m_ref.tif *00/*/*-DEM_32m_*mask.tif *00/*/*-DEM_32m_*bareground.tif
    rm -v *00/*/hma_srtm*.tif *00/*/conus_lidar*tif
fi

#Now create masks for each 32m DEM using dem_mask.py
#Check settings for dem_mask - MODSCAG, SNODAS, TOA, etc.
#Note that parallel won't work with input for modscag username
#parallel --delay 1 'dem_mask.py {}' ::: $list_32m
#parallel --delay 1 'dem_mask.py --bareground_thresh 5 {}' ::: $list_32m
#Use TOA mask
#parallel --delay 1 'dem_mask.py --toa {}' ::: $list_32m
#SnowEx
#parallel --jobs 32 --delay 1 'dem_mask.py --filter not_forest+not_water --toa {}' ::: $list_32m

#Clean up existing pc_align runs
if $clean ; then
    rm -rv *00/*/*align *00/*/*trans.tif 
fi

#Do the co-registration
#set pc_align_wrapper threads to 2 or 4
#set pc_align-wrapper max displacement

#Reference is a csv of xyz points
if True ; then
    #For output of filter_glas.py: ICESat-1 GLAS points over reference surfaces
    list_csv=$(echo $list_2m | sed 's/DEM_2m.tif/DEM_32m_GLAH14_conus_refdemfilt_ref.csv/g')
    parallel --link --progress -j $njobs --delay 3 "dem_coreg.sh {1} {2}" ::: $list_2m ::: $list_csv
    #This should work for direct parallel command substitution, but there is still some syntax issue
    ##parallel --plus --progress -j $njobs --delay 3 "dem_coreg.sh {} {/DEM_2m.tif/DEM_32m_GLAH14_conus_refdemfilt_ref.csv}" ::: $list_2m
else
    #Reference is a gridded DEM
    #Should be existing 2-m mosaic for entire region
    #ref=/nobackup/deshean/rpcdem/lidar/conus_lidar_1m.vrt
    ref=$1
    #If we don't find it, create a reference
    if [ ! -n "$ref" ] ; then 
        ref=ref_DEM_2m_mos-tile-0.tif
        if [ ! -e $ref ] ; then 
            #For larger extents, need to run dem_mosaic_validtiles.py
            #Create mosaic, removing any low-quality products from QB02 or IK01
            #dem_mosaic -o ref_DEM_2m_mos $(echo $list_2m | tr ' ' '\n' | grep -v QB02 | grep -v IK01)
            dem_mosaic -o ref_DEM_2m_mos $(echo $list_2m | tr ' ' '\n' | grep -v '_101*_101')
            ref=ref_DEM_2m_mos-tile-0.tif
            #Coreg round number
            n=1
        fi
    fi

    #Using gridded DEM as reference
    parallel --progress -j $njobs --delay 3 "dem_coreg.sh {} $ref" ::: $list_2m
fi

#Now create new weighted average mosaics, burn into reference DEM, and rerun the co-registration
#Should use make_mos.sh - need to specify that we want to use -DEM_2m_trans.tif
