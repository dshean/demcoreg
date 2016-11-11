#! /bin/bash

#David Shean
#dshean@gmail.com

#Wrapper for ASP utility pc_align to align two input point clouds
#Written for specific case of reference point data (ATM) and source DEM (WV), hence variable names

#To do
#Improved mask handling - make sure working with new pygeotools apply_mask.py
#If two grids are provided, clip the ref grid w/ shp, buffer src grid by max_disp, clip src grid

echo
if [ "$#" -ne 2 ] && [ "$#" -ne 3 ] ; then
    echo "Usage is $0 ref_pts.[csv|tif] src_dem.tif [trans.txt]"
    echo
    exit 1
fi

#This returns appropriate status when cmd is piped to tee
set -o pipefail

#atm=lat_lon_z_180_rock.csv
atm=$1
#dem=20100709_1534_1030010005B3D100_103001000612DC00-DEM_4x.tif
dem=$2

itrans=false
if [ ! -z "$3" ] ; then
    itrans=$3
fi

if [ ! -e $atm ] ; then
    echo "Unable to locate $atm"
    exit 1
fi

if [ ! -e $dem ] ; then
    echo "Unable to locate $dem"
    exit 1
fi

#Set this to compute trans+rot, instead of default trans
#rot=true
rot=false

#Set number of threads to use
ncpu=8

#Number of iterations (default 1000)
n_iter=2000

#Max displacement (0 is default)
#NOTE: can extract this from last column of input csv (v*dt)
#max_disp=1000
#max_disp=200
max_disp=40
#max_disp=20

pc_align_opt=''
point2dem_opt=''

#Extract info about reference point cloud to use for alignment
if [ ${atm##*.} = 'csv' ] ; then
    #ATM "resolution" should be ~10 m - finer for repeat tracks
    atm_res=10.0
    sample_pts=true
    fmt="--csv-format '1:lat 2:lon 3:height_above_datum'"
    #This x y is ECEF
    #fmt="--csv-format '1:x 2:y 3:height_above_datum'"
    pc_align_opt+=$fmt
    ref_type='point'
elif [ ${atm##*.} = 'tif' ] ; then
    #ASP PC
    if gdalinfo $atm | grep -q POINT_OFFSET ; then
        atm_res=0.5
        sample_pts=false
        ref_type='asp_pc'
    #Gridded geotif
    else
        atm_res=$(gdalinfo $atm | grep 'Pixel Size' | awk -F'[(,)]' '{print $2}')
        sample_pts=false
        ref_type='grid'
    fi
else
    echo "Unrecognized input extension:"
    echo $atm
    exit 1
fi

#Extract info about source DEM to be aligned
if [ ${dem##*.} = 'csv' ] ; then
    dem_type='point'
    #ATM "resolution" should be ~10 m - finer for repeat tracks
    dem_res=10.0
    fmt="--csv-format '1:lat 2:lon 3:height_above_datum'"
    sample_pts=false
    use_point2dem=true
    usemask=false
    #This x y is ECEF
    #fmt="--csv-format '1:x 2:y 3:height_above_datum'"
    pc_align_opt+=$fmt
elif [ ${dem##*.} = 'tif' ] ; then
    #ASP PC
    if gdalinfo $dem | grep -q POINT_OFFSET ; then
        dem_type='asp_pc'
        dem_res=0.5
        sample_pts=false
        use_point2dem=true
        #Could theoretically use mask here
        usemask=false
    #Gridded geotif
    else
        dem_type='grid'
        sample_pts=true
        use_point2dem=false
        usemask=true
        dem_res=$(gdalinfo $dem | grep 'Pixel Size' | awk -F'[(,)]' '{print $2}')
        point2dem_opt+=" --threads $ncpu --tr $dem_res"
        dem_ndv=$(gdalinfo $dem | grep 'NoData Value' | awk -F'=' '{print $2}')
        point2dem_opt+=" --nodata-value $dem_ndv"
        proj="$(proj_select.py $dem)"
        point2dem_opt+=" --t_srs '$proj'"
    fi
else
    echo "Unrecognized input extension:"
    echo $dem
    exit 1
fi

#If rotation desired, must use point2dem
if $rot ; then
    use_point2dem=true
fi

#Set this to force output of new PC and run point2dem
#Otherwise, run apply_dem_translation
if ! $use_point2dem ; then
    pc_align_opt+=" --compute-translation-only"
    #only want to mask if we're going to apply trans to orig, unmasked DEM
    usemask=true
fi

#If DEM is more dense, use as "reference" and trans-reference
#if false ; then
if [ $(echo "a=($dem_res < $atm_res); a" | bc -l) -eq 1 ] ; then
    trans_source=false
    trans_reference=true
else
    trans_source=true
    trans_reference=false
fi

#When we have an initial transform - useful to transform 
#Could also just do this with geotransform and scale z values of final grid
#Must also comment out --num-iterations below
if $itrans ; then
    itrans=$3
    pc_align_opt+=" --initial-transform $itrans"
    #Set number of iterations to 0
    n_iter=0
    #Set displacement to -1, otherwise does one iteration
    max_disp=-1
else
    pc_align_opt+=" --max-displacement $max_disp" 
    pc_align_opt+=" --num-iterations $n_iter" 
fi

#Should be pc_align ref src
pc_align_opt+=" --threads $ncpu --datum WGS_1984"

#Point-to-point works MUCH better for limited ref points in a plane
#Looks like trans-reference point-to-point wins for Thwaites test cases

align_method="point-to-point"
#align_method="point-to-plane"
pc_align_opt+=" --alignment-method $align_method"

outlier_ratio=0.75
#outlier_ratio=0.95
pc_align_opt+=" --outlier-ratio $outlier_ratio"

#Do we want to use fsaa here?  Input is already gridded at $dem_res
#point2dem_opt+='--fsaa'

#Script for sampling points
#This is Pleiades location
sample_script='/u/deshean/src/Tools/build/point_to_dem_dist'
if [ ! -x "$sample_script" ] ; then
    #Assume we're on dido
    sample_script="/Users/dshean/src/asp_tools/point_to_dem_dist"
fi

function runcmd () {
    cmd=$1
    logfile=$2
    echo | tee -a $logfile
    echo $cmd | tee -a $logfile
    echo | tee -a $logfile
    eval time $cmd | tee -a $logfile
    #This is the return status we want
    status=$?
    echo | tee -a $logfile
    return $status
}

function sample () {
    ref=$1
    dem=$2
    outdir=$3
    logfile=$4
    myout=$outdir/$(basename ${dem%.*})
    if [ "${ref##*.}" == "csv" ] ; then 
        #Run initial sampling
        echo "Sampling $dem" | tee -a $logfile
        cmd="$sample_script -o $myout $ref $dem"
        runcmd "$cmd" $logfile
        if [ ! -e ${myout}-sample.csv ] ; then
            echo "No valid overlapping points!" | tee -a $logfile
            exit 1
        else
            cmd="robust_stats.py ${myout}-sample.csv"    
            runcmd "$cmd" $logfile
        fi
    elif [ "${ref##*.}" == "tif" ] ; then
        echo "Computing elevation difference" | tee -a $logfile
        #This will print eul stats, which is what we want
        cmd="compute_dz.py $ref $dem"
        runcmd "$cmd" $logfile
        eul=${ref%.*}_$(basename ${dem%.*})_dz_eul.tif
        if [ ! -e $eul ] ; then
            echo "No valid samples!" | tee -a $logfile
            exit 1
        else
            cmd="robust_stats.py $eul"
            runcmd "$cmd" $logfile
        fi
        if [ -e ${eul%.*}_rate.tif ] ; then
            rm -v ${eul%.*}_rate.tif 
        fi
        mv -v $eul $outdir 
    fi
}

function logfile_init () {
    logfile=$1
    #Careful, this will overwrite the logfiles
    echo -n > $logfile
    date | tee -a $logfile
    echo | tee -a $logfile
    echo "Input resolution:" | tee -a $logfile
    echo "ref: " $(basename $atm) $atm_res | tee -a $logfile 
    echo "src: " $(basename $dem) $dem_res | tee -a $logfile
    echo "ref type:" $ref_type | tee -a $logfile
    echo "align_method: $align_method" | tee -a $logfile
    echo "trans_source: $trans_source" | tee -a $logfile
    echo "trans_reference: $trans_reference" | tee -a $logfile
    echo "max_disp: $max_disp" | tee -a $logfile
    echo | tee -a $logfile
}

dem_orig=$dem

#Note pc_align should append a '-' to the specified outprefix
#out=pc_align_ref_$(basename ${atm%.*})_src_$(basename ${dem%.*})
outdir=${dem_orig%.*}_align
if [ "$ref_type" == "grid" ] ; then
    outdir=${dem_orig%.*}_grid_align
fi
mkdir -pv $outdir
out=$outdir/$(basename ${dem_orig%.*})
pc_align_opt+=" -o $out"

#This masks the input DEM to the region around the points
if $usemask ; then
    #Identify mask filename
    mask=${atm%.*}_mask.tif
    #mask=${atm%.*}_masked.tif
    #mask=rainier_surfaces_for_alignment.shp
    #mask=rainier_surfaces_for_alignment_all_rivervalleys.shp
    if [ -e $mask ] ; then
        if [ "${mask##*.}" == "tif" ] ; then
            if [ ! -e ${dem%.*}_masked.tif ] ; then
                echo; echo "Applying mask to input DEM"
                cmd="apply_mask.py $dem $mask"
                #cmd="apply_mask_new.py $dem $mask"
                runcmd "$cmd" $logfile
                echo
            fi
            dem=${dem%.*}_masked.tif
        elif [ "${mask##*.}" == "shp" ] ; then
            if [ ! -e ${dem%.*}_shpclip.tif ] ; then
                echo; echo "Masking input DEM for alignment"
                cmd="clip_raster_by_shp.sh $dem $mask"
                runcmd "$cmd" $logfile
                echo
            fi
            dem=${dem%.*}_shpclip.tif
        fi
    else
        usemask=false
    fi
fi

if $trans_source ; then
    out_dem=${out}-trans_source-DEM.tif
    
    if $use_point2dem && ! $usemask ; then
        pc_align_opt+=' --save-transformed-source-points'
    fi

    logfile=${out}-trans_source.log
    logfile_init $logfile

    if $sample_pts ; then
        sample $atm $dem $outdir $logfile
    fi

    #Want to load as many ref points as possible
    #pc_align_opt+=' --max-num-reference-points 10000000 --max-num-source-points 10000000'

    #if [ ! -e ${out}-trans_source.tif ] ; then
    if [ ! -e ${out}-trans_source-end_errors.csv ] ; then
        cmd="pc_align $pc_align_opt $atm $dem"
        runcmd "$cmd" $logfile
        mv -v ${out}-beg_errors.csv ${out}-trans_source-beg_errors.csv
        mv -v ${out}-end_errors.csv ${out}-trans_source-end_errors.csv
    fi

    if [ ! -e $out_dem ] ; then
        if $use_point2dem ; then 
            #This is a hack, because we need to rerun to generate output PC using entire orig DEM, not masked area
            if $usemask ; then
                echo | tee -a $logfile
                echo "Applying pc_align transform to original input DEM" | tee -a $logfile
                echo | tee -a $logfile
                #Set number of iterations to 0
                pc_align_opt=$(echo $pc_align_opt | awk '{ gsub("--num-iterations '$n_iter'","--num-iterations 0"); print }')
                #Set displacement to -1, otherwise does one iteration
                pc_align_opt=$(echo $pc_align_opt | awk '{ gsub("--max-displacement '$max_disp'","--max-displacement -1"); print }')
                #Add final transform
                itrans=${out}-transform.txt
                pc_align_opt+=" --initial-transform $itrans"
                pc_align_opt+=' --save-transformed-source-points'
                #If point, use head -10 to generate dummy ref
                cmd="pc_align $pc_align_opt $atm $dem_orig"
                #runcmd "$cmd" $logfile
                echo $cmd | tee -a $logfile
                eval time $cmd
                mv -v ${out}-beg_errors.csv ${out}-trans_source-beg_errors.csv
                mv -v ${out}-end_errors.csv ${out}-trans_source-end_errors.csv
            fi
            cmd="point2dem $point2dem_opt ${out}-trans_source.tif"
            runcmd "$cmd" $logfile
            if [ $? -eq 0 ] ; then 
               rm ${out}-trans_source.tif
            fi
        else
            #Could be an issue to read log while writing
            cmd="apply_dem_translation.py $dem_orig $logfile"
            runcmd "$cmd" $logfile
            #Keep everything in output directory
            mv -v ${dem_orig%.*}_trans.tif $out_dem 
        fi
    fi

    if $sample_pts ; then
        sample $atm $out_dem $outdir $logfile
    fi

    date | tee -a $logfile
fi

if $trans_reference ; then
    out_dem=${out}-trans_reference-DEM.tif
    
    if $use_point2dem && ! $usemask ; then
        pc_align_opt+=' --save-inv-transformed-reference-points'
    fi

    logfile=${out}-trans_reference.log
    logfile_init $logfile

    if $sample_pts ; then
        sample $atm $dem $outdir $logfile
    fi

    #Ref should always be as dense as possible
    #Ref is now DEM, want to load as many src points as possible
    #pc_align_opt+=' --max-num-reference-points 100000 --num-source-points 100000000'
    #pc_align_opt+=' --max-num-reference-points 10000000 --max-num-source-points 10000000'

    if [ ! -e ${out}-trans_reference-end_errors.csv ] ; then
        cmd="pc_align $pc_align_opt $dem $atm"
        runcmd "$cmd" $logfile
        mv -v ${out}-beg_errors.csv ${out}-trans_reference-beg_errors.csv
        mv -v ${out}-end_errors.csv ${out}-trans_reference-end_errors.csv
    fi

    if [ ! -e $out_dem ] ; then
        if $use_point2dem ; then
            #This is a hack, because we need to rerun to generate output PC using entire orig DEM, not masked area
            if $usemask ; then
                echo | tee -a $logfile
                echo "Applying pc_align transform to original input DEM" | tee -a $logfile
                echo | tee -a $logfile
                #Set number of iterations to 0
                pc_align_opt=$(echo $pc_align_opt | awk '{ gsub("--num-iterations '$n_iter'","--num-iterations 0"); print }')
                #Set displacement to -1, otherwise does one iteration
                pc_align_opt=$(echo $pc_align_opt | awk '{ gsub("--max-displacement '$max_disp'","--max-displacement -1"); print }')
                #Add final transform
                #itrans=${out}-inverse-transform.txt
                itrans=${out}-transform.txt
                pc_align_opt+=" --initial-transform $itrans"
                pc_align_opt+=' --save-inv-transformed-reference-points' 
                #If point, use head -10 to generate dummy ref
                cmd="pc_align $pc_align_opt $dem_orig $atm"
                #Don't want to write to pollute logfile
                #runcmd "$cmd" $logfile
                echo $cmd | tee -a $logfile
                eval time $cmd
                mv -v ${out}-beg_errors.csv ${out}-trans_source-beg_errors.csv
                mv -v ${out}-end_errors.csv ${out}-trans_source-end_errors.csv
            fi
            cmd="point2dem $point2dem_opt ${out}-trans_reference.tif"  
            runcmd "$cmd" $logfile
            if [ $? -eq 0 ] ; then 
                rm ${out}-trans_reference.tif
            fi
        else
            cmd="apply_dem_translation.py $dem_orig $logfile"
            runcmd "$cmd" $logfile
            mv -v ${dem_orig%.*}_trans.tif $out_dem 
        fi
    fi

    if $sample_pts ; then
        sample $atm $out_dem $outdir $logfile
    fi
    
    date | tee -a $logfile
fi
