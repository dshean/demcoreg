#! /bin/bash

#David Shean
#dshean@gmail.com

#This will co-register a large batch of DEMs

#Run tpfe1 or bro node
#To check out devel bro node:
#qsub ~/bin/devel.pbs
#cd topdir

#site=scg
#site=conus
#site=hma
site=ngozumpa

#Coreg round number
n=1

njobs=10
#njobs=4
#Default is for all WV/GE/QB subdir
list_32m=$(ls -Sr *00/dem*/*-DEM_32m.tif)
list_dir=$(ls -Sr *00/dem*/*-DEM_32m.tif | awk -F'/' '{print $1}')
list_32m_done=$(ls -Sr *00/dem*/*-DEM_32m_trans.tif)

#parallel_log=dem_coreg_round2_log

echo "Creating done and todo lists"
done_list=${site}_coreg_round${n}_done.txt
todo_list=${site}_coreg_round${n}_todo.txt
echo -n > $done_list
echo -n > $todo_list
for i in $list_32m
do
    if echo $list_32m_done | grep -q ${i%.*} ; then
        echo $i >> $done_list 
    else 
        echo $i >> $todo_list
    fi
done

ndone=$(cat $done_list | wc -l)
ntodo=$(cat $todo_list | wc -l)
echo
echo "$ndone done"
echo "$ntodo todo"
echo

list_32m=$(cat ${todo_list})

#This was recovery
#list_32m=$(ls -Sr $(cat incomplete_killed.txt | sed 's#$#/dem*/*-DEM_32m.tif#' | grep -v WV02_20150702_103001004453A700_10300100456F4A00))

#If we have existing orthoimages, compute top-of-atmosphere reflectance 
#Uses toa.sh, toa.py, and dglib from https://github.com/dshean/dgtools 
#parallel --jobs 16 --delay 2 'toa.sh {}' ::: $list_dir

#Clean up existing masks
#rm */*/*-DEM_32m_ref.tif */*/*-DEM_32m_*mask.tif */*/*-DEM_32m_*perc.tif
#rm */*/hma_srtm*
#rm */*/conus_lidar*

#Now create masks for each 32m DEM
#Check settings for dem_mask - MODSCAG, SNODAS, TOA, etc.
#parallel --jobs 32 --delay 1 'dem_mask.py --toa {}' ::: $list_32m
#SnowEx
#parallel --jobs 32 --delay 1 'dem_mask.py --filter not_forest+not_water --toa {}' ::: $list_32m
#Note that parallel won't work with input for modscag username

#Clean up existing pc_align runs
#rm -r */*/*align */*/*trans.tif 

#Do the co-registration
#set pc_align_wrapper threads to 2 or 4
#set pc_align-wrapper max displacement
list_2m=$(echo $list_32m | sed 's/-DEM_32m.tif/-DEM_2m.tif/g')
#list_8m=$(ls -Sr *00/dem*/*-DEM_8m.tif)
#parallel --progress --results $parallel_log -j $njobs --delay 3 'dem_coreg.sh {}' ::: $list_2m
#echo $list_2m
parallel --progress -j $njobs --delay 3 'dem_coreg.sh {}' ::: $list_2m
#parallel --progress --results coreg_8m_log -j 24 --delay 1 'dem_coreg.sh {}' ::: $list_8m

#Now create new weighted average mosaics, burn into reference DEM, and rerun the co-registration
