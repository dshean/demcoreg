#! /usr/bin/env python

import sys
import os
import argparse

from osgeo import gdal
import numpy as np

from pygeotools.lib import iolib, malib, geolib, warplib

from demcoreg import coreglib, dem_mask

def getparser():
    parser = argparse.ArgumentParser(description="Perform DEM co-registration using old algorithms")
    parser.add_argument('ref_fn', type=str, help='Reference DEM filename')
    parser.add_argument('src_fn', type=str, help='Source DEM filename to be shifted')
    parser.add_argument('-mode', type=str, default='nuth', choices=['ncc', 'sad', 'nuth', 'none'], \
            help='Type of co-registration to use')
    parser.add_argument('-nomask', action='store_true', help='Do not apply mask to isolate static surfaces')
    parser.add_argument('-max_offset', type=float, default=100, \
            help='Maximum expected horizontal offset in meters')
    parser.add_argument('-outdir', default=None, help='Output directory')
    return parser

def main(argv=None):
    parser = getparser()
    args = parser.parse_args()

    #Should check that files exist
    dem1_fn = args.ref_fn
    dem2_fn = args.src_fn
    mode = args.mode
    max_offset_m = args.max_offset

    print("\nReference: %s" % dem1_fn)
    print("Source: %s" % dem2_fn)
    print("Mode: %s\n" % mode)

    dem1_ds = gdal.Open(dem1_fn, gdal.GA_ReadOnly)
    dem2_ds = gdal.Open(dem2_fn, gdal.GA_ReadOnly)
    #Preserve copy of original DEM 2 geotransform
    dem2_gt_orig = np.array(dem2_ds.GetGeoTransform())

    #Make sure the input datasets have the same resolution/extent
    dem1_clip_ds, dem2_clip_ds = warplib.memwarp_multi([dem1_ds, dem2_ds], res='max', extent='intersection')

    #Compute size of NCC and SAD search window in pixels 
    res = float(geolib.get_res(dem1_clip_ds, square=True)[0])
    max_offset_px = (max_offset_m/res) + 1
    #print(max_offset_px)
    pad = (int(max_offset_px), int(max_offset_px))

    #This will be updated geotransform for dem2
    dem2_gt = np.array(dem2_clip_ds.GetGeoTransform())

    #Load the arrays
    dem1 = iolib.ds_getma(dem1_clip_ds, 1)
    dem2 = iolib.ds_getma(dem2_clip_ds, 1)

    if not args.nomask:
        #Mask glaciers, vegetated slopes
        #static_mask = ~(dem_mask.get_lulc_mask(dem1_clip_ds, mask_glaciers=True, \
        #        filter='not_forest+not_water', bareground_thresh=60))
        #Mask glaciers
        static_mask = ~(dem_mask.get_icemask(dem1_clip_ds))
        dem1 = np.ma.array(dem1, mask=static_mask)
        dem2 = np.ma.array(dem2, mask=static_mask)

    #Compute difference for unaligned inputs
    print("Elevation difference stats for uncorrected input DEMs")
    #Shouldn't need to worry about common mask here, as both inputs are ma

    diff_euler = np.ma.array(dem1 - dem2)
    if diff_euler.count() == 0:
        sys.exit("No valid overlap between input DEMs")

    diff_stats = malib.print_stats(diff_euler)
    med_bias = diff_stats[5]

    print("Computing sub-pixel offset between DEMs using mode: %s" % mode)

    #Sum of absolute differences
    if mode == "sad":
        m, int_offset, sp_offset = coreglib.compute_offset_sad(dem1, dem2, pad=pad)
        
        #Geotransform has negative y resolution, so don't need negative sign
        #np array is positive down
        #GDAL coordinates are positive up
        xshift_m = sp_offset[1]*dem2_gt[1]
        yshift_m = sp_offset[0]*dem2_gt[5]

    #Normalized cross-correlation of clipped, overlapping areas
    elif mode == "ncc":
        prefilter = False 
        write_nccfig = True
        m, int_offset, sp_offset, fig = coreglib.compute_offset_ncc(dem1, dem2, \
                pad=pad, prefilter=prefilter, plot=write_nccfig)

        if write_nccfig:
            dst_fn = '%s_%s_nccfig.png' % (os.path.splitext(dem2_fn)[0], \
                    os.path.splitext(os.path.split(dem1_fn)[1])[0]) 
            fig.savefig(dst_fn)

        xshift_m = sp_offset[1]*dem2_gt[1]
        yshift_m = sp_offset[0]*dem2_gt[5]

    #Nuth and Kaab (2011)
    elif mode == "nuth":
        #Generate slope and aspect maps
        print("Computing slope and aspect")
        dem1_slope = geolib.gdaldem_mem_ds(dem1_clip_ds, processing='slope', returnma=True)
        dem1_aspect = geolib.gdaldem_mem_ds(dem1_clip_ds, processing='aspect', returnma=True)

        #Apply common mask

        #Compute relationship between elevation difference, slope and aspect
        fit_param = coreglib.compute_offset_nuth(diff_euler, dem1_slope, dem1_aspect)

        #fit_param[0] is magnitude of shift vector
        #fit_param[1] is direction of shift vector
        #fit_param[2] is mean bias divided by tangent of mean slope 
        print(fit_param)

        xshift_m = fit_param[0]*np.sin(np.deg2rad(fit_param[1]))
        yshift_m = fit_param[0]*np.cos(np.deg2rad(fit_param[1]))

        """
        #This is recomputed and applied below
        med_slope = malib.fast_median(dem1_slope)
        med_bias = fit_param[2]*np.tan(np.deg2rad(med_slope))
        print("med vertial bias: ", med_bias)
        """

    elif mode == "all":
        print("Not yet implemented")
        #Want to compare all methods, average offsets
        #m, int_offset, sp_offset = coreglib.compute_offset_sad(dem1, dem2)
        #m, int_offset, sp_offset = coreglib.compute_offset_ncc(dem1, dem2)

    #This is a hack to apply the computed median bias correction for shpclip area only
    elif mode == "none":
        print("Skipping alignment, writing out DEM with median bias over static surfaces removed")
        dst_fn = os.path.splitext(dem2_fn)[0]+'_med%0.2f.tif' % med_bias
        iolib.writeGTiff(dem2_orig + med_bias, dst_fn, dem2_ds)
        sys.exit()

    #Apply the horizontal shift to the original dataset
    dem2_ds_align = coreglib.apply_xy_shift(dem2_ds, xshift_m, yshift_m)

    #Write out aligned dataset, but without vertical offset applied
    write_align = False 
    if write_align:
        align_fn = '%s_align_x%+0.2f_y%+0.2f.tif' % (os.path.splitext(dem2_fn)[0], xshift_m, yshift_m)
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

        #Compute elevation difference
        diff_euler_align = dem1 - dem2
        
        print("Elevation difference stats for aligned inputs")
        diff_stats_align = malib.print_stats(diff_euler_align)

        #Use median for vertical bias removal
        zshift_m = diff_stats_align[5]

        #Write out aligned eulerian difference map for clipped extent with vertial offset removed
        dst_fn = '%s_%s_align_x%+0.2f_y%+0.2f_z%+0.2f_dz_eul.tif' % (os.path.splitext(dem2_fn)[0], \
                os.path.splitext(os.path.split(dem1_fn)[1])[0], xshift_m, yshift_m, zshift_m) 
        print("Writing out aligned difference map with median vertical offset removed") 
        iolib.writeGTiff(diff_euler_align - zshift_m, dst_fn, dem1_clip_ds) 
        
        #Write out aligned dem_2 with vertial offset removed
        dst_fn = '%s_align_x%+0.2f_y%+0.2f_z%+0.2f.tif' % (os.path.splitext(dem2_fn)[0], \
                xshift_m, yshift_m, zshift_m)
        print("Writing out shifted dem2 with median vertical offset removed: %s" % dst_fn)
        dem2_align = iolib.ds_getma(dem2_ds_align)
        iolib.writeGTiff(dem2_align + zshift_m, dst_fn, dem2_ds_align) 

if __name__ == "__main__":
    main()
