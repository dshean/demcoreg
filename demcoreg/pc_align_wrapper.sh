#! /bin/bash

#David Shean
#dshean@gmail.com

#Wrapper for ASP utility pc_align to align two input point clouds or DEMs

#To do
#Improved mask handling - make sure working with new pygeotools apply_mask.py

set -e

usage()
{
    echo; echo "Usage is $0 ref_pts.{csv|tif} src_dem.tif [pc_align options]"
    echo "Note that options must come AFTER the two filenames are specified"
    echo; "Also see pc_align usage"; echo; exit 1
}

if [ "$#" -lt 2 ] ; then
    usage
fi

#This returns appropriate status when cmd is piped to tee
set -o pipefail

#Reference DEM/PC
ref=$1
shift
if [ ! -e $ref ] ; then
    echo "Unable to locate $ref"
	usage
fi

#Source DEM/PC to shift
dem=$1
shift
if [ ! -e $dem ] ; then
    echo "Unable to locate $dem"
    usage
fi

## Define default parameters

#Alignment method
align_method=point-to-point
#align_method=point-to-plane

#This is maximum number of points to use
max_points=1000000
#max_points=10000000

#Initial transformation
itrans=false

#Set this to enable rotation in addition to translation 
rot=false

#Set this to enable scaling in addition to translation and rotation
scale=false

#Set number of threads to use
ncpu=4

#Number of iterations (default 1000)
n_iter=2000

#Max displacement (0 is default)
#NOTE: can extract this from last column of input csv (v*dt)
max_disp=10

#Set outlier ratio
#outlier_ratio=0.95
outlier_ratio=0.75

#Sample points and compute stats before and after coreg
#Uses raster differencing for gridded DEMs and point sampling for PC
#Can slow down processing for large inputs
sample_pts=true

## Parse any user-specfiied parameters
#See pc_align usage
#Note that there is no checking here, user must specify properly
while [ "$1" != "" ]; do
	case $1 in
		--alignment-method)	    shift
								align_method=$1
								;;
		--initial-transform)	shift
								itrans=$1
								;;
		--num-iterations )		shift
								n_iter=$1
								;;
		--max-displacement)		shift
								max_disp=$1
								;;
		--max-points)			shift
								max_points=$1
								;;
		-r | --rot )        	rot=true 
								;;
		-s | --scale )        	scale=true 
								;;
		-h | --help )           usage
								exit
								;;
		* )                     usage
								exit 1
	esac
	shift
done

#ATM "resolution" should be ~10 m - finer for repeat tracks
#ref_res=10.0
#fmt="--csv-format '1:lat 2:lon 3:height_above_datum'"
#This x y is ECEF
#fmt="--csv-format '1:x 2:y 3:height_above_datum'"

#ICESat-1 along-track spacing
ref_res=140.0
#Format from demcoreg ICESat-1 csv
fmt="--csv-format '3:lat 4:lon 5:height_above_datum'"

pc_align_opt=''
point2dem_opt=''

#Extract info about reference point cloud to use for alignment
if [ ${ref##*.} = 'csv' ] ; then
    pc_align_opt+=$fmt
    ref_type='point'
elif [ ${ref##*.} = 'tif' ] ; then
    #ASP PC
    if gdalinfo $ref | grep -q POINT_OFFSET ; then
        #This info should now be available in PC.tif header
        ref_res=0.5
        ref_type='asp_pc'
    #Gridded geotif
    else
        ref_res=$(gdalinfo $ref | grep 'Pixel Size' | awk -F'[(,)]' '{print $2}')
        ref_type='grid'
    fi
else
    echo "Unrecognized input extension:"
    echo $ref
    exit 1
fi

#Extract info about source DEM to be aligned
if [ ${dem##*.} = 'csv' ] ; then
    dem_type='point'
    dem_res=$ref_res
    use_point2dem=true
    usemask=false
    pc_align_opt+=$fmt
elif [ ${dem##*.} = 'tif' ] ; then
    #ASP PC
    if gdalinfo $dem | grep -q POINT_OFFSET ; then
        dem_type='asp_pc'
        dem_res=0.5
        use_point2dem=true
        #Could theoretically use mask here
        usemask=false
    #Gridded geotif
    else
        dem_type='grid'
        #Can used compute_diff.py here
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

if $scale ; then
    #Use similarity transform to solve for scaling
    align_method='similarity-point-to-point'
    #align_method='similarity-least-squares'
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
if [ $(echo "a=($dem_res < $ref_res); a" | bc -l) -eq 1 ] ; then
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
    pc_align_opt+=" --initial-transform $itrans"
    #Set number of iterations to 0
    n_iter=0
    #Set displacement to -1, otherwise does one iteration
    max_disp=-1
else
    pc_align_opt+=" --max-displacement $max_disp" 
    pc_align_opt+=" --num-iterations $n_iter" 
fi

#This was test for Hexagon mapping camera
#--initial-ned-translation "0 0 -6000"

#Should be pc_align ref src
pc_align_opt+=" --threads $ncpu --datum WGS_1984"

#Use specified alignment method
pc_align_opt+=" --alignment-method $align_method"

pc_align_opt+=" --outlier-ratio $outlier_ratio"

#Script for sampling points
#This is now bundled with demcoreg
sample_script=sample_raster_at_pts.py

#This was older sampling utility, fast, but required compiling
#This is Pleiades location
#sample_script='/u/deshean/src/Tools/build/point_to_dem_dist'
#if [ ! -x "$sample_script" ] ; then
#    #Assume we're on dido
#    sample_script="/Users/dshean/src/asp_tools/point_to_dem_dist"
#fi

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
        #Older sampling formats
        #cmd="$sample_script -o $myout $ref $dem"
        #sample_fn=${myout}-sample.csv
        cmd="$sample_script $dem $ref"
        ref_basename=$(basename $ref)
        sample_fn=${dem%.*}_${ref_basename%.*}_sample.csv
        runcmd "$cmd" $logfile
        if [ ! -e $sample_fn ] ; then
            echo "No valid overlapping points!" | tee -a $logfile
            exit 1
        else
            cmd="robust_stats.py -col -1 $sample_fn"    
            runcmd "$cmd" $logfile
        fi
    elif [ "${ref##*.}" == "tif" ] ; then
        echo "Computing elevation difference" | tee -a $logfile
        #This will print diff stats, which is what we want
        cmd="compute_diff.py $ref $dem"
        runcmd "$cmd" $logfile
        diff=${ref%.*}_$(basename ${dem%.*})_diff.tif
        if [ ! -e $diff ] ; then
            echo "No valid samples!" | tee -a $logfile
            exit 1
        else
            cmd="robust_stats.py $diff"
            runcmd "$cmd" $logfile
        fi
        if [ -e ${diff%.*}_rate.tif ] ; then
            rm -v ${diff%.*}_rate.tif 
        fi
        mv -v $diff $outdir 
    fi
}

function logfile_init () {
    logfile=$1
    #Careful, this will overwrite the logfiles
    echo -n > $logfile
    date | tee -a $logfile
    echo | tee -a $logfile
    echo "Input resolution:" | tee -a $logfile
    echo "ref: " $(basename $ref) $ref_res | tee -a $logfile 
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
#out=pc_align_ref_$(basename ${ref%.*})_src_$(basename ${dem%.*})
outdir=${dem_orig%.*}_pt_align
if [ "$ref_type" == "grid" ] ; then
    outdir=${dem_orig%.*}_grid_align
fi

if [ -d $outdir ] ; then
    echo "Output directory already exists: $outdir"
    exit
    #mv $outdir ${outdir}_old
fi

mkdir -pv $outdir
out=$outdir/$(basename ${dem_orig%.*})
pc_align_opt+=" -o $out"

#This masks the input DEM to the region around the points
if $usemask ; then
    #Identify mask filename
    mask=${ref%.*}_mask.tif
    #mask=${ref%.*}_masked.tif
    #mask=rainier_surfaces_for_alignment.shp
    #mask=rainier_surfaces_for_alignment_all_rivervalleys.shp
    if [ -e $mask ] ; then
        if [ "${mask##*.}" == "tif" ] ; then
            if [ ! -e ${dem%.*}_masked.tif ] ; then
                echo; echo "Applying mask to input DEM"
                cmd="apply_mask.py -extent raster $dem $mask"
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
        sample $ref $dem $outdir $logfile
    fi

    #Want to load as many ref points as possible
    pc_align_opt+=" --max-num-reference-points $max_points --max-num-source-points $max_points"

    #if [ ! -e ${out}-trans_source.tif ] ; then
    if [ ! -e ${out}-trans_source-end_errors.csv ] ; then
        cmd="pc_align $pc_align_opt $ref $dem"
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
                cmd="pc_align $pc_align_opt $ref $dem_orig"
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
        sample $ref $out_dem $outdir $logfile
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
        sample $ref $dem $outdir $logfile
    fi

    #Ref should always be as dense as possible
    #Ref is now DEM, want to load as many src points as possible
    pc_align_opt+=" --max-num-reference-points $max_points --max-num-source-points $max_points"

    if [ ! -e ${out}-trans_reference-end_errors.csv ] ; then
        cmd="pc_align $pc_align_opt $dem $ref"
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
                cmd="pc_align $pc_align_opt $dem_orig $ref"
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
        sample $ref $out_dem $outdir $logfile
    fi
    
    date | tee -a $logfile
fi
