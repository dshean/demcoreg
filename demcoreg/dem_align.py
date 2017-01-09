#! /usr/bin/env python

#David Shean
#dshean@gmail.com

#Utility to co-register two rasters
#This was written in 2012, and needs major cleanup

#Input two uncorrected, unclipped DEMs and a mask for static regions
#Clip both DEMs using mask
#Align the clipped DEMs
#Update georef info and remove elevation offset to DEM2

#Next:
#Need better way to clip input datasets to common mask of static surfaces
#Using gdalwarp is a bit clunky
#This uses PIL to draw polygon
#http://geospatialpython.com/2011/02/clip-raster-using-shapefile.html
#At least remove the temporary clipped files
#mapnik to render shapefile for given extent?
#Define shapefile
#Convert to pbm (or binary ndarray) for specified extent

#To do:
#Clean up path - implement os.path.join(outdir, fn)

import sys
import os
import argparse

from osgeo import gdal
import numpy as np

from pygeotools.lib import iolib
from pygeotools.lib import malib
from pygeotools.lib import geolib
from pygeotools.lib import filtlib
from pygeotools.lib import warplib

from demcoreg import coreglib

def getparser():
    parser = argparse.ArgumentParser(description="Perform DEM co-registration using old algorithms")
    parser.add_argument('ref_fn', type=str, help='Reference DEM filename')
    parser.add_argument('src_fn', type=str, help='Source DEM filename to be moved')
    parser.add_argument('-mode', type=str, default='ncc', choices=['ncc', 'sad', 'nuth', 'none'], help='Type of co-registration to use')
    parser.add_argument('-mask_fn', type=str, default=None, help='Mask filename')
    parser.add_argument('-outdir', default=None, help='Output directory')
    return parser

def main(argv=None):
    parser = getparser()
    args = parser.parse_args()

    #Should check filenames exist

    #Input DEMs
    #This is the reference
    dem1_fn = args.ref_fn
    #This is the source to be moved
    dem2_fn = args.src_fn

    print("\nReference: %s" % dem1_fn)
    print("Source: %s" % dem2_fn)

    #Input mask for static areas
    #Note mask must be in same projected coordinate system as input rasters
    #For areas not covering rock, can use the intersection polygon
    mask_fn = args.mask_fn 
    if mask_fn is not None:
        print("Mask: %s\n" % mask_fn)

    mode = args.mode
    print("Mode: %s\n" % mode)

    dem1_ds = gdal.Open(dem1_fn, gdal.GA_ReadOnly)
    dem2_ds = gdal.Open(dem2_fn, gdal.GA_ReadOnly)

    dem1_gt_orig = np.array(dem1_ds.GetGeoTransform())
    dem2_gt_orig = np.array(dem2_ds.GetGeoTransform())

    #Check to see if input datasets have the same extent (already clipped) and gt
    clip = np.any(dem1_gt_orig != dem2_gt_orig) or geolib.extent_compare(geolib.ds_extent(dem1_ds), geolib.ds_extent(dem2_ds))

    dem1_clip_ds = dem1_ds
    dem2_clip_ds = dem2_ds
    
    #These are needed by older filtering code, should be removed (can get from ds)
    dem1_clip_fn = dem1_fn
    dem2_clip_fn = dem2_fn

    if clip:
        if mask_fn is not None:
            #NOTE: this does not resample to ensure both are same res
            #Should be better to clip first, then resample
            print("Clipping input DEMs to input polygon\n")
            dem1_clip_ds = geolib.clip_raster_by_shp(dem1_fn, mask_fn)
            dem2_clip_ds = geolib.clip_raster_by_shp(dem2_fn, mask_fn)
            dem1_clip_fn = dem1_clip_ds.GetFileList()[0]
            dem2_clip_fn = dem2_clip_ds.GetFileList()[0]
            #Need to be careful here, as these are not necessarily set 
            #Also, I think GDAL has trouble when source file for an open ds disappears
            #os.remove(dem1_clip_fn)
            #os.remove(dem2_clip_fn)

        #Make sure the input datasets have the same resolution/extent
        dem1_clip_ds, dem2_clip_ds = warplib.memwarp_multi([dem1_clip_ds, dem2_clip_ds], res='max', extent='intersection')

    #This will be updated geotransform for dem2
    dem2_gt = np.array(dem2_clip_ds.GetGeoTransform())

    dem1 = iolib.ds_getma(dem1_clip_ds, 1)
    dem2 = iolib.ds_getma(dem2_clip_ds, 1)

    #Compute common mask and apply to all products
    #common_mask = dem1.mask + dem2.mask
    common_mask = malib.common_mask([dem1, dem2]) 

    #Compute difference for unaligned inputs
    diff_euler = np.ma.array(dem1-dem2, mask=common_mask)
    
    print("Elevation difference stats for uncorrected input DEMs")
    diff_stats = malib.print_stats(diff_euler)
    med_bias = diff_stats[5]

    #Now load the full dem2, not just the clipped 
    #dem2_orig = iolib.ds_getma(dem2_ds, 1)
    dem2_orig = np.copy(dem2)

    #Filter the input DEMs to remove blunders and smooth before alignment
    filt = False
    
    if filt:
        print("Filtering clipped DEMs before alignment")
   
        #DEM elevation cutoff filter
        """
        perc = (1.0, 99.0)
        dem1lim = malib.calcperc(dem1,perc)
        dem2lim = malib.calcperc(dem2,perc)
        rangelim = (min(dem1lim[0], dem2lim[0]), max(dem1lim[1], dem2lim[1]))

        dem1 = np.ma.masked_outside(dem1, rangelim[0], rangelim[1])
        dem2 = np.ma.masked_outside(dem2, rangelim[0], rangelim[1])
        """
        
        #Apply range filter
        #dem1 = filtlib.range_fltr_lowresDEM(dem1, dem1_fn)
        #dem2 = filtlib.range_fltr_lowresDEM(dem2, dem2_fn)
        dem1 = filtlib.range_fltr_lowresDEM(dem1, dem1_clip_fn)
        dem2 = filtlib.range_fltr_lowresDEM(dem2, dem2_clip_fn)

        #Gaussian Smoothing filter
        #Note - there are issues with nan handling for the scipy.ndimage.gaussian_filter!
        sigma = 1
        size = 3
        dem1 = filtlib.gauss_fltr_astropy(dem1, size=size, sigma=sigma)
        dem2 = filtlib.gauss_fltr_astropy(dem2, size=size, sigma=sigma)
        
        common_mask = malib.common_mask([dem1, dem2])
    
        #Note: use original mask here, as filter will fill some holes
        diff_euler = np.ma.array(dem1-dem2, mask=common_mask)
    
        print("Elevation difference stats for filtered input DEMs")
        diff_stats = malib.print_stats(diff_euler)
        med_bias = diff_stats[5]

        #Should use a slope filter here

    dfilt = False 

    #This is implemented in filtlib now

    if dfilt:
        #Elevation differences shouldn't be more than 20-30 m/yr
        #Or just do perc clip
        #Want to limit this to static surfaces, as advecting ice could have large differences
        #perc = (2, 98)
        #difflim = malib.calcperc(diff_euler, perc)
        #Throw out everything beyond 2*mad
        #mad_sigma = 2
        #difflim = (diff_stats[5] - mad_sigma*diff_stats[6], diff_stats[5] + mad_sigma*diff_stats[6])
        #Compute appropriate range limits based on dt
        #Over static surfaces, this should be really small <1 m
        #However, could be intersection, over several years
        #rangelim = (-60, 60)
        
        #Use most conservative of two    
        #difflim = (max(difflim[0], rangelim[0]), min(difflim[1], rangelim[1]))

        #print "Excluding locations with differences beyond", difflim

        #diff_euler = np.ma.masked_outside(diff_euler, *difflim)
        #diff_euler = filtlib.mad_fltr(diff_euler, mad_sigma=2)
        diff_euler = filtlib.perc_fltr(diff_euler, (2,98))

        print("Elevation difference stats after difference filter")
        diff_stats = malib.print_stats(diff_euler)
        med_bias = diff_stats[5]

        common_mask = common_mask + diff_euler.mask
   
    #Update masks
    dem1.__setmask__(common_mask)
    dem2.__setmask__(common_mask)

    print("Computing sub-pixel offset between DEMs using mode: %s" % mode)

    #Sum of absolute differences
    if mode == "sad":
        m, int_offset, sp_offset = coreglib.compute_offset_sad(dem1, dem2)
        
        #Geotransform has negative y resolution, so don't need negative sign
        #np array is positive down
        #Projection is positive up
        xshift_m = sp_offset[1]*dem2_gt[1]
        yshift_m = sp_offset[0]*dem2_gt[5]

    #Normalized cross-correlation of clipped, overlapping areas
    elif mode == "ncc":
        prefilter = False 

        if prefilter:
            #Do LoG filter here 
            print("Applying edge-detection filter to DEMs")
            #Note, this propagates Nans and greatly reduces valid data area
            sigma = 1
            import scipy.ndimage
            dem1_LoG = malib.nanfill(dem1, scipy.ndimage.filters.gaussian_laplace, sigma) 
            dem2_LoG = malib.nanfill(dem2, scipy.ndimage.filters.gaussian_laplace, sigma) 
            m, int_offset, sp_offset, fig = coreglib.compute_offset_ncc(dem1_LoG, dem2_LoG, plot=False)
        else:
            m, int_offset, sp_offset, fig = coreglib.compute_offset_ncc(dem1, dem2, plot=False)

        write_nccfig = True
        if write_nccfig:
            dst_fn = '%s_%s_nccfig.png' % (os.path.splitext(dem2_fn)[0], os.path.splitext(os.path.split(dem1_fn)[1])[0]) 
            import matplotlib.pyplot as plt
            plt.savefig(dst_fn)

        xshift_m = sp_offset[1]*dem2_gt[1]
        yshift_m = sp_offset[0]*dem2_gt[5]
        #xshift_m = 93.6283764258
        #yshift_m = -24

    #Nuth and Kaab (2011)
    elif mode == "nuth":
        #NOTE: want to generate slope/aspect maps from the result of the above filter
        dem1_fltr_fn = os.path.splitext(dem1_clip_fn)[0]+'_fltr.tif'
        iolib.writeGTiff(dem1, dem1_fltr_fn, dem1_clip_ds)

        #These just call subprocess - update for new GDAL utility API
        dem1_slope = geolib.gdaldem_slope(dem1_fltr_fn)
        dem1_aspect = geolib.gdaldem_aspect(dem1_fltr_fn)

        common_mask = common_mask + dem1_slope.mask + dem1_aspect.mask
        dem1_slope.__setmask__(common_mask)
        dem1_aspect.__setmask__(common_mask)

        #Compute relationship between elevation difference, slope and aspect
        fit_param = coreglib.compute_offset_nuth(diff_euler, dem1_slope, dem1_aspect)
        print(fit_param)

        #fit_param[0] is magnitude of shift vector
        #fit_param[1] is direction of shift vector
        #fit_param[2] is mean bias divided by tangent of mean slope 

        xshift_m = fit_param[0]*np.sin(np.deg2rad(fit_param[1]))
        yshift_m = fit_param[0]*np.cos(np.deg2rad(fit_param[1]))

        #mean_slope = dem1_slope.mean()
        #mean_bias = fit_param[2]*np.tan(np.deg2rad(mean_slope))
        med_slope = np.ma.median(dem1_slope)
        med_bias = fit_param[2]*np.tan(np.deg2rad(med_slope))
    
        print("med bias: ", med_bias)

    #Want to compare all methods, take average?
    elif mode == "all":
        m, int_offset, sp_offset = coreglib.compute_offset_sad(dem1, dem2)
        m, int_offset, sp_offset = coreglib.compute_offset_ncc(dem1, dem2)

    #This is a hack to apply the computed median bias correction for shpclip area only
    elif mode == "none":
        print("Skipping alignment, writing out DEM with median bias over static surfaces removed")
        dst_fn = os.path.splitext(dem2_fn)[0]+'_med%0.2f.tif' % med_bias
        iolib.writeGTiff(dem2_orig + med_bias, dst_fn, dem2_ds)
        sys.exit()

    print("X shift: ", xshift_m, "Y shift: ", yshift_m)
    
    dem2_gt_shift = np.copy(dem2_gt_orig)
    #Don't know why these were subtracted before - this may be what Nuth and Kaab output needs?
    #dem2_gt_shift[0] -= xshift_m
    #dem2_gt_shift[3] -= yshift_m
    dem2_gt_shift[0] += xshift_m
    dem2_gt_shift[3] += yshift_m

    print("Original geotransform:", dem2_gt_orig)
    print("Updated geotransform:", dem2_gt_shift)

    #Write out copy of dem2 with updated geotransform
    dem2_ds_align = iolib.mem_drv.CreateCopy('', dem2_ds, 1)
    dem2_ds_align.SetGeoTransform(dem2_gt_shift)

    #Write out aligned dataset, but without offset applied
    write_align = True 
    if write_align:
        align_fn = '%s_align_x%+0.2f_y%+0.2f.tif' % (os.path.splitext(dem2_fn)[0], xshift_m, yshift_m)
        #align_fn = os.path.splitext(dem2_fn)[0]+'_align.tif'
        print("Writing out shifted dem2 (no vertical offset): %s" % align_fn)
        warplib.writeout(dem2_ds_align, align_fn)

    #Write out original eulerian difference map
    write_origdiff = True 
    if write_origdiff:
        print("Writing out original euler difference map for common intersection before alignment")
        dst_fn = '%s_%s_orig_dz_eul.tif' % (os.path.splitext(dem2_fn)[0], os.path.splitext(os.path.split(dem1_fn)[1])[0]) 
        iolib.writeGTiff(diff_euler, dst_fn, dem1_clip_ds)

    #Write out aligned eulerian difference map
    write_aligndiff = True
    if write_aligndiff:
        dem1_ds, dem2_ds = warplib.memwarp_multi([dem1_ds, dem2_ds_align], res='max', extent='intersection') 

        dem1 = iolib.ds_getma(dem1_ds, 1)
        dem2 = iolib.ds_getma(dem2_ds, 1)
        common_mask = malib.common_mask([dem1, dem2]) 

        #Compute difference for aligned inputs
        diff_euler_align = np.ma.array(dem1-dem2, mask=common_mask)
        
        print("Elevation difference stats for aligned inputs")
        diff_stats_align = malib.print_stats(diff_euler_align)
        #Use median bias
        offset = diff_stats_align[5]

        dem2_clip_ds = dem2_ds_align 
        dem2_clip_fn = align_fn

        if clip:
            if mask_fn is not None:
                print("Clipping input DEMs to input polygon\n")
                #dem1_clip_ds should remain unchanged, no need to reproduce
                dem2_clip_ds = geolib.clip_raster_by_shp(align_fn, mask_fn)
                dem2_clip_fn = dem2_clip_ds.GetFileList()[0]
                #Need to be careful here, as these are not necessarily set 
                #Also, I think GDAL has trouble when source file for an open ds disappears
                #os.remove(dem1_clip_fn)
                #os.remove(dem2_clip_fn)

            #Make sure the input datasets have the same resolution/extent
            dem1_clip_ds, dem2_clip_ds = warplib.memwarp_multi([dem1_clip_ds, dem2_clip_ds], res='max', extent='intersection')
        
            dem1_clip = iolib.ds_getma(dem1_clip_ds, 1)
            dem2_clip = iolib.ds_getma(dem2_clip_ds, 1)

            common_mask = malib.common_mask([dem1_clip, dem2_clip]) 

            #Compute difference for aligned inputs
            diff_euler_clip = np.ma.array(dem1_clip-dem2_clip, mask=common_mask)
            
            print("Elevation difference stats for aligned inputs")
            diff_stats_clip = malib.print_stats(diff_euler_clip)
            #Use median bias
            offset = diff_stats_clip[5]

            if dfilt:
                diff_euler_clip = filtlib.mad_fltr(diff_euler_clip, mad_sigma=2)        
                print("Elevation difference stats for aligned inputs after difference filter")
                diff_stats_clip = malib.print_stats(diff_euler_clip)
                offset = diff_stats_clip[5]
                common_mask = common_mask + diff_euler_clip.mask

            #Write out aligned eulerian difference map for clipped extent with offset removed
            dst_fn = '%s_%s_align_x%+0.2f_y%+0.2f_z%+0.2f_dz_eul_clip.tif' % (os.path.splitext(dem2_fn)[0], os.path.splitext(os.path.split(dem1_fn)[1])[0], xshift_m, yshift_m, offset) 
            print("Writing out aligned, clipped difference map with median vertical offset removed") 
            iolib.writeGTiff(diff_euler_clip - offset, dst_fn, dem1_clip_ds) 
        
        #Write out aligned, unclipped dem2 dataset with offset removed
        #Note: offset will come from clipped static area if possible
        #dst_fn = '%s_align_%0.2fm.tif' % (os.path.splitext(dem2_fn)[0], offset)
        dst_fn = '%s_align_x%+0.2f_y%+0.2f_z%+0.2f.tif' % (os.path.splitext(dem2_fn)[0], xshift_m, yshift_m, offset)
        print("Writing out shifted dem2 with median vertical offset removed: %s" % dst_fn)
        dem2_align = iolib.ds_getma(dem2_ds_align)
        iolib.writeGTiff(dem2_align + offset, dst_fn, dem2_ds_align) 

if __name__ == "__main__":
    main()
