#! /bin/bash 

#This will co-register a single DEM
#See dem_coreg_all.sh for batch processing

#Exit if error
set -e

if [ "$#" -ne 2 ]; then
    echo "Usage: $(basename $0) source_dem ref_dem"
    exit 1
fi

#Input should be highest res version of DEM (i.e., DEM_2m.tif)
dem=$1
if [ ! -e $dem ] ; then
    echo "Unable to find source DEM: $dem"
    exit
fi

ref=$2
if [ ! -e $ref ] ; then
    echo "Unable to find ref DEM: $ref"
    exit
fi

#Set this to apply translation to original 2, 8 and 32 m DEMs
update=true

echo
echo $dem
echo $ref
echo ${ref##*.}

demdir=$(dirname $dem)
dembase=$(basename $dem)
dembase=$(echo ${dembase%.*} | awk -F'-' '{print $1}')

if [ "${ref##*.}" == "csv" ] ; then 
    refdem_masked=$ref
    #This will be pc_align output directory
    outdir=${dembase%.*}_pt_align
else
    #This will be pc_align output directory
    outdir=${dembase%.*}_grid_align

    #This is DEM_32m reference mask output by dem_mask.py
    #dem_mask=$demdir/${dembase}-DEM_32m_ref.tif
    dem_mask=$(ls $demdir/${dembase}*_ref.tif | tail -1)

    if [ ! -e $dem_mask ] ; then
        echo "Unable to find reference DEM mask, need to run dem_mask.py"
        echo "Ideally, run on low-res version of DEM:"
        echo "dem_mask.py $demdir/${dembase}-DEM_32m.tif"
        exit
    fi

    refdem=$demdir/$(basename $ref)
    refdem=${refdem%.*}_warp.tif
    refdem_masked=${refdem%.*}_masked.tif

    if [ ! -e $refdem ] ; then
        #Clip the reference DEM to the DEM_32m extent
        echo
        echo "Warping reference DEM to low-res mask extent"
        warptool.py -r near -te $dem_mask -tr $dem_mask -t_srs $dem_mask -outdir $demdir $ref
        #This does a bunch of unnecessary calculations, could do quick python call to get valid px count
        count=$(robust_stats.py $refdem | grep 'Count:' | awk '{print $2}')
        mincount=100
        if [ "$count" -lt "$mincount" ] ; then
            echo "Not enough valid pixels in clipped reference DEM: (${count} < ${mincount})" 
            exit
        fi
        echo
        echo "Applying low-res mask to high-res reference DEM"
        apply_mask.py -extent intersection $refdem $dem_mask
        count=$(robust_stats.py $refdem_masked | grep 'Count:' | awk '{print $2}')
        mincount=100
        if [ "$count" -lt "$mincount" ] ; then
            echo "Not enough valid pixels in masked reference DEM: (${count} < ${mincount})" 
            exit
        fi
        #Remove low-res check products
        rm $refdem_masked $refdem
        echo
        echo "Warping high-res reference DEM to appropriate extent"
        warptool.py -r cubic -te $dem -tr $ref -t_srs $dem -outdir $demdir $ref
    fi

    #Check if refdem has valid pixels

    #This avoids writing another copy of ref, but is slower
    #NOTE: assumes projection of $dem and $ref are identical.  Need to implement better get_extent with -t_srs option
    #dem_extent=$(~/src/demtools/get_extent.py $dem)
    #echo "Creating vrt of high-res reference DEM clipped to appropriate extent"
    #gdalbuildvrt -tr 1 1 -te $dem_extent -tap -r nearest ${refdem%.*}_warp.vrt $ref
    #refdem=${refdem%.*}_warp.vrt

    #Mask the ref using valid pixels in DEM_32m_ref.tif product
    refdem_masked=$demdir/${dembase}-$(basename $ref)
    refdem_masked=${refdem_masked%.*}_masked.tif
    if [ ! -e $refdem_masked ] ; then
        echo
        echo "Applying low-res mask to high-res reference DEM"
        apply_mask.py -out_fn $refdem_masked -extent intersection $refdem $dem_mask
    fi
fi

#Should check if refdem_masked has valid pixels

if [ -e $refdem_masked ] ; then
    #Actually run pc_align
    #Default is point-to-point ICP
    pc_align_wrapper.sh $refdem_masked $dem

    #Now apply the translation to DEMs with different resolution
    if $update ; then
        apply_trans_all.sh $demdir
    fi
fi
