#! /bin/bash 

#This will co-register a single DEM
#See dem_coreg_all.sh for batch processing

#Exit if error
set -e

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
    #Hardcode reference DEM
    #1-arcsec SRTM (30 m) for HMA
    #ref=/nobackup/deshean/rpcdem/hma/srtm1/hma_srtm_gl1.vrt
    #Round 1 after ICESat
    #ref=/nobackupp8/deshean/hma/hma1_2016dec22/hma_2m_tile_20170220/hma_2m.vrt
    #CONUS
    #Need to create vrt with 1 arcsec over areas where 1/3 is not avail
    #1-arcsec NED (30 m) for CONUS
    #ref=/nobackup/deshean/rpcdem/ned1/ned1_tiles_glac24k_115kmbuff.vrt
    #1/3-arcsec NED (10 m) for CONUS
    #ref=/nobackup/deshean/rpcdem/ned13/ned13_tiles_glac24k_115kmbuff.vrt
    #1-m lidar vrt
    #ref=/nobackup/deshean/rpcdem/lidar/conus_lidar_1m.vrt
    #2-m WV DEM mosaic, first round
    #ref=/nobackup/deshean/conus/dem2/conus_coreg1_mos_2m_tile/conus_coreg1_mos_2m.vrt
    #2-m WV DEM mosaic, second round
    #ref=/nobackup/deshean/conus/dem2/conus_coreg2_mos_2m_tile/conus_coreg2_mos_2m.vrt
    #2-m WV DEM mosaic, second round
    #ref=/nobackup/deshean/conus/dem2/conus_coreg3_mos_2m_tile/conus_coreg3_mos_2m.vrt
    #ref=/nobackup/deshean/conus/dem2/conus_coreg3_mos_2m_summer_tile/conus_coreg3_mos_2m_summer.vrt
    #Ngozumpa
    #2-m WV DEM mosaic
    #ref=/nobackupp8/deshean/hma/ngozumpa2/hma1_2016dec22/stereo/ngozumpa2_2m_ref-tile-0.tif
    #Mashel
    #ref=/nobackupp8/deshean/mashel/Rasters/dem_3ft_unitmeters-adj.tif
    #SnowEx
    #SBB
    #ref=/nobackup/deshean/snowex/aso/USCOSB20160926f1a1_dsm_vf_bf_masked.tif
    #GM
    #ref=/nobackup/deshean/snowex/aso/USCOGM20160926f1a1_dsm_1p5m_vf_bf_masked_trim.tif
    #ref=/nobackup/deshean/snowex/stereo/gm/20160925_gm_2m_trans-tile-0.tif

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
        stats=$(~/src/demtools/stats.py $refdem)
        echo
        echo $stats
        count=$(echo $stats | awk '{print $2}')
        mincount=100
        if [ "$count" -lt "$mincount" ] ; then
            echo "Not enough valid pixels in clipped reference DEM: (${count} < ${mincount})" 
            exit
        fi
        echo
        echo "Applying low-res mask to high-res reference DEM"
        apply_mask.py -extent intersection $refdem $dem_mask
        stats=$(~/src/demtools/stats.py $refdem_masked)
        echo
        echo $stats
        count=$(echo $stats | awk '{print $2}')
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
