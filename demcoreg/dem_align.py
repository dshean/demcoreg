#! /usr/bin/env python

import sys
import os
import argparse

from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt

from pygeotools.lib import iolib, malib, geolib, warplib

from demcoreg import coreglib, dem_mask

from imview.lib import pltlib, gmtColormap
cpt_rainbow = gmtColormap.get_rainbow()

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

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.splitext(dem2_fn)[0] + '_dem_align'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outprefix = '%s_%s' % (os.path.splitext(os.path.split(dem2_fn)[-1])[0], \
            os.path.splitext(os.path.split(dem1_fn)[-1])[0]) 
    outprefix = os.path.join(outdir, outprefix)

    print("\nReference: %s" % dem1_fn)
    print("Source: %s" % dem2_fn)
    print("Mode: %s" % mode)
    print("Output: %s\n" % outprefix)

    dem1_ds = gdal.Open(dem1_fn, gdal.GA_ReadOnly)
    dem2_ds = gdal.Open(dem2_fn, gdal.GA_ReadOnly)
    #Preserve copy of original DEM 2 geotransform
    dem2_gt_orig = np.array(dem2_ds.GetGeoTransform())

    #Make sure the input datasets have the same resolution/extent
    #Use projection of source DEM
    dem1_clip_ds, dem2_clip_ds = warplib.memwarp_multi([dem1_ds, dem2_ds], \
            res='max', extent='intersection', t_srs=dem2_ds)

    #Compute size of NCC and SAD search window in pixels 
    res = float(geolib.get_res(dem1_clip_ds, square=True)[0])
    max_offset_px = (max_offset_m/res) + 1
    #print(max_offset_px)
    pad = (int(max_offset_px), int(max_offset_px))

    #This will be updated geotransform for dem2
    dem2_gt = np.array(dem2_clip_ds.GetGeoTransform())

    #Load the arrays
    dem1_orig = iolib.ds_getma(dem1_clip_ds, 1)
    dem2_orig = iolib.ds_getma(dem2_clip_ds, 1)

    if not args.nomask:
        #Mask glaciers, vegetated slopes
        #static_mask = ~(dem_mask.get_lulc_mask(dem1_clip_ds, mask_glaciers=True, \
        #        filter='not_forest+not_water', bareground_thresh=60))
        #Mask glaciers
        static_mask = ~(dem_mask.get_icemask(dem1_clip_ds))
        dem1 = np.ma.array(dem1_orig, mask=static_mask)
        dem2 = np.ma.array(dem2_orig, mask=static_mask)
    else:
        dem1 = dem1_orig
        dem2 = dem2_orig

    #Compute difference for unaligned inputs
    print("Elevation difference stats for uncorrected input DEMs")
    #Shouldn't need to worry about common mask here, as both inputs are ma
    diff_euler = dem1 - dem2

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
            dst_fn = '%s_nccfig.png' % outprefix
            fig.savefig(dst_fn)

        xshift_m = sp_offset[1]*dem2_gt[1]
        yshift_m = sp_offset[0]*dem2_gt[5]

    #Nuth and Kaab (2011)
    elif mode == "nuth":
        #Generate slope and aspect maps
        print("Computing slope and aspect")
        dem1_slope = geolib.gdaldem_mem_ds(dem1_clip_ds, processing='slope', returnma=True)
        dem1_aspect = geolib.gdaldem_mem_ds(dem1_clip_ds, processing='aspect', returnma=True)

        #Compute relationship between elevation difference, slope and aspect
        fit_param, f_nuth = coreglib.compute_offset_nuth(diff_euler, dem1_slope, dem1_aspect)

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
        dst_fn = outprefix+'_med%0.1f.tif' % med_bias
        iolib.writeGTiff(dem2_orig + med_bias, dst_fn, dem2_ds)
        sys.exit()

    #Apply the horizontal shift to the original dataset
    dem2_ds_align = coreglib.apply_xy_shift(dem2_ds, xshift_m, yshift_m)

    #String to append to output filenames
    xy_shift_str = '_%s_x%+0.1f_y%+0.1f' % (mode, xshift_m, yshift_m)

    #Write out aligned dataset, but without vertical offset applied
    write_align = False 
    if write_align:
        dst_fn = outprefix + '%s_align.tif' % (xy_shift_str)
        print("Writing out shifted dem2 (no vertical offset): %s" % dst_fn)
        warplib.writeout(dem2_ds_align, dst_fn)

    #Write out original eulerian difference map
    write_origdiff = True 
    if write_origdiff:
        print("Writing out original euler difference map for common intersection before alignment")
        dst_fn = outprefix + '_orig_dz_eul.tif' 
        iolib.writeGTiff(diff_euler, dst_fn, dem1_clip_ds)

    #Write out aligned eulerian difference map
    write_aligndiff = True
    if write_aligndiff:
        dem1_clip_ds, dem2_clip_ds = warplib.memwarp_multi([dem1_ds, dem2_ds_align], \
                res='max', extent='intersection', t_srs=dem2_ds_align) 

        dem1 = iolib.ds_getma(dem1_clip_ds, 1)
        dem2 = iolib.ds_getma(dem2_clip_ds, 1)

        #Compute elevation difference
        diff_euler_align = dem1 - dem2

        #Recompute the mask
        if not args.nomask:
            #Mask glaciers, vegetated slopes
            #static_mask = ~(dem_mask.get_lulc_mask(dem1_clip_ds, mask_glaciers=True, \
            #        filter='not_forest+not_water', bareground_thresh=60))
            #Mask glaciers
            static_mask = ~(dem_mask.get_icemask(dem1_clip_ds))
            diff_euler_align_masked = np.ma.array(diff_euler_align, mask=static_mask) 
        else:
            diff_euler_align_masked = diff_euler_align

        print("Elevation difference stats for aligned inputs")
        diff_stats_align = malib.print_stats(diff_euler_align_masked)

        #Use median for vertical bias removal
        zshift_m = diff_stats_align[5]
        xyz_shift_str = xy_shift_str+'_z%+0.1f' % zshift_m

        #Write out aligned eulerian difference map for clipped extent with vertial offset removed
        dst_fn = outprefix + '%s_align_dz_eul.tif' % xyz_shift_str
        print("Writing out aligned difference map with median vertical offset removed") 
        iolib.writeGTiff(diff_euler_align - zshift_m, dst_fn, dem1_clip_ds) 
        
        #Write out aligned dem_2 with vertial offset removed
        dst_fn = outprefix + '%s_align.tif' % xyz_shift_str
        print("Writing out shifted dem2 with median vertical offset removed: %s" % dst_fn)
        dem2_align = iolib.ds_getma(dem2_ds_align)
        iolib.writeGTiff(dem2_align + zshift_m, dst_fn, dem2_ds_align) 

        makeplot = True
        if makeplot:
            dst_fn = outprefix + '%s_nuth_plot.png' % xyz_shift_str
            print("Writing Nuth and Kaab plot: %s" % dst_fn)
            f_nuth.savefig(dst_fn, dpi=300, bbox_inches='tight', pad_inches=0)

            print("Creating final plot")
            dem1_hs = geolib.gdaldem_mem_ma(dem1_orig, dem1_clip_ds, returnma=True)
            dem2_hs = geolib.gdaldem_mem_ma(dem2_orig, dem2_clip_ds, returnma=True)
            f,axa = plt.subplots(2, 3, figsize=(11, 8.5))          
            for ax in axa.ravel():
                ax.set_facecolor('k')
                pltlib.hide_ticks(ax)
            dem_clim = malib.calcperc(dem1_orig, (2,98))
            axa[0,0].imshow(dem1_hs, cmap='gray')
            axa[0,0].imshow(dem1_orig, cmap='cpt_rainbow', clim=dem_clim, alpha=0.6)
            pltlib.add_scalebar(axa[0,0], res=res)
            axa[0,0].set_title('Reference DEM')
            axa[0,1].imshow(dem2_hs, cmap='gray')
            axa[0,1].imshow(dem2_orig, cmap='cpt_rainbow', clim=dem_clim, alpha=0.6)
            axa[0,1].set_title('Source DEM')
            axa[0,2].imshow(~(np.ma.getmaskarray(diff_euler)), clim=(0,1), cmap='gray')
            axa[0,2].set_title('Surfaces for co-registration')
            dz_clim = malib.calcperc_sym(diff_euler, (2, 98))
            im = axa[1,0].imshow(-diff_euler, cmap='RdBu', clim=dz_clim)
            axa[1,0].set_title('Elev. Diff. Before')
            im = axa[1,1].imshow(-diff_euler_align_masked, cmap='RdBu', clim=dz_clim)
            axa[1,1].set_title('Elev. Diff. After')
            im = axa[1,2].imshow(-diff_euler_align, cmap='RdBu', clim=dz_clim)
            axa[1,2].set_title('Elev. Diff. Final')
            #Tried to insert Nuth fig here
            #ax_nuth.change_geometry(1,2,1)
            #f.axes.append(ax_nuth)
            pltlib.add_cbar(axa[1,2], im, label='Elev Diff (m)')
            suptitle = '%s\nx: %+0.2fm, y: %+0.2fm, z: %+0.2fm' % \
                    (os.path.split(outprefix)[-1], xshift_m, yshift_m , zshift_m)
            f.suptitle(suptitle)
            f.tight_layout()
            plt.subplots_adjust(top=0.90)
            dst_fn = outprefix + '%s_align.png' % xyz_shift_str
            print("Writing out shifted dem2 with median vertical offset removed: %s" % dst_fn)
            f.savefig(dst_fn, dpi=300, bbox_inches='tight', pad_inches=0)
            #plt.show()

if __name__ == "__main__":
    main()
