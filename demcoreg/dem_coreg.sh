#! /bin/bash 

#David Shean
#dshean@gmail.com

#This will co-register a single DEM
#See dem_coreg_all.sh for batch processing

#Input should be highest res version of DEM (i.e., DEM_2m.tif)
dem=$1

if [ ! -e $dem ] ; then
    echo "Unable to find source DEM: $dem"
    exit
fi

#Define the reference DEM
#Need to create vrt with 1 arcsec over areas where 1/3 is not avail
#1-arcsec NED (30 m) for CONUS
#ref=/nobackup/deshean/rpcdem/ned1/ned1_tiles_glac24k_115kmbuff.vrt
#1/3-arcsec NED (10 m) for CONUS
#ref=/nobackup/deshean/rpcdem/ned13/ned13_tiles_glac24k_115kmbuff.vrt
#1-arcsec SRTM (30 m) for HMA
ref=/nobackup/deshean/rpcdem/hma/srtm1/hma_srtm_gl1.vrt

if [ ! -e $ref ] ; then
    echo "Unable to find ref DEM: $ref"
    exit
fi

demdir=$(dirname $dem)
dembase=$(basename $dem)
#This will be pc_align output directory
outdir=${dembase%.*}_grid_align
dembase=$(echo $dembase | awk -F'-' '{print $1}')

#This is DEM_32m reference mask output by dem_mask.py
dem_mask=$demdir/${dembase}-DEM_32m_ref.tif

if [ ! -e $dem_mask ] ; then
    echo "Unable to find reference DEM mask, need to run dem_mask.py"
    exit
fi

#Clip the reference DEM to the DEM_32m extent
warptool.py -te $dem -tr $ref -t_srs $dem -outdir $demdir $ref

refdem=$demdir/$(basename $ref)
refdem=${refdem%.*}_warp.tif
#Mask the ref using valid pixels in DEM_32m_ref.tif product
apply_mask.py -extent intersection $refdem $dem_mask
refdem_masked=${refdem%.*}_masked.tif

#point-to-point
pc_align_wrapper.sh $refdem_masked $dem

cd $demdir
log=$(ls -t $outdir/*.log | head -1)
if [ -e $log ] ; then 
    if grep -q 'Translation vector' $log ; then
        apply_dem_translation.py ${dembase}-DEM_32m.tif $log
        apply_dem_translation.py ${dembase}-DEM_8m.tif $log
        ln -sf $outdir/*DEM.tif ${dembase}-DEM_2m_trans.tif
        compute_dh.py $(basename $refdem) ${dembase}-DEM_8m_trans.tif
    fi
fi
