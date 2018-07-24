#! /usr/bin/env python

#Todo
#Better outlier removal
#Check Nuth and Kaab bin median
#Implement check for empty diff

import sys
import os
import argparse
import subprocess

from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt

from pygeotools.lib import iolib, malib, geolib, warplib

from demcoreg import coreglib, dem_mask

from imview.lib import pltlib, gmtColormap
cpt_rainbow = gmtColormap.get_rainbow()

def getparser():
    parser = argparse.ArgumentParser(description="Perform DEM co-registration using old algorithms", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('ref_fn', type=str, help='Reference DEM filename')
    parser.add_argument('src_fn', type=str, help='Source DEM filename to be shifted')
    parser.add_argument('-mode', type=str, default='nuth', choices=['ncc', 'sad', 'nuth', 'none'], \
            help='Type of co-registration to use')
    parser.add_argument('-nomask', action='store_true', \
            help='By default, input DEMs are masked to limit co-registration for static surfaces. \
            Set this to use all surfaces')
    filter_choices = ['rock', 'rock+ice', 'rock+ice+water', 'not_forest', 'not_forest+not_water', 'none']
    parser.add_argument('-filter', type=str, default='not_forest', choices=filter_choices, \
            help='Define areas to use as reference surfaces for co-registration')
    parser.add_argument('-tiltcorr', action='store_true', \
            help='After preliminary translation, fit plane to residual elevation offsets and remove')
    parser.add_argument('-tol', type=float, default=0.02, \
            help='When iterative translation magnitude is below this tolerance (meters), break and write out corrected DEM') 
    parser.add_argument('-max_offset', type=float, default=100, \
            help='Maximum expected horizontal offset in meters')
    parser.add_argument('-outdir', default=None, help='Output directory')
    return parser

def get_mask(ds, dem_fn, filter='not_forest', mask_glaciers=True, bareground_thresh=60):
    #This logic needs to be cleaned up
    if filter != 'none':
        #Mask glaciers, vegetated slopes
        static_mask = dem_mask.get_lulc_mask(ds, mask_glaciers=mask_glaciers, filter=filter, bareground_thresh=bareground_thresh)
    else:
        #Mask glaciers only
        static_mask = dem_mask.get_icemask(ds)
    #Top-of-atmosphere reflectance threshold (requires orthoimage and output from toa.sh)
    toa_fn = dem_mask.get_toa_fn(dem_fn)
    #toa_fn = None
    if toa_fn is not None:
        toa_ds = warplib.memwarp_multi_fn([toa_fn,], res=ds, extent=ds, t_srs=ds, r='cubicspline')[0]
        toa_mask = dem_mask.get_toa_mask(toa_ds)
        static_mask = np.logical_and(static_mask, toa_mask)
    #Return final mask, ready to be applied
    return ~(static_mask)

def compute_offset(dem1_ds, dem2_ds, dem2_fn, mode='nuth', max_offset_m=100, remove_outliers=True, \
        apply_mask=True, filter=filter):
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
    dem1 = iolib.ds_getma(dem1_clip_ds, 1)
    dem2 = iolib.ds_getma(dem2_clip_ds, 1)

    #Compute difference for unaligned inputs
    print("Elevation difference stats for uncorrected input DEMs")
    #Shouldn't need to worry about common mask here, as both inputs are ma
    diff_euler = dem2 - dem1

    static_mask = None
    if apply_mask:
        #Need dem2_fn here to find TOA fn
        static_mask = get_mask(dem2_clip_ds, dem2_fn, filter=filter)
        dem1 = np.ma.array(dem1, mask=static_mask)
        dem2 = np.ma.array(dem2, mask=static_mask)
        diff_euler = np.ma.array(diff_euler, mask=static_mask)
        static_mask = np.ma.getmaskarray(diff_euler)

    if diff_euler.count() == 0:
        sys.exit("No overlapping, unmasked pixels shared between input DEMs")

    #Compute stats for new masked difference map
    diff_stats = malib.print_stats(diff_euler)
    dz = diff_stats[5]

    #This needs further testing
    if remove_outliers:
        med = diff_stats[5]
        nmad = diff_stats[6]
        f = 3
        rmin = med - f*nmad
        rmax = med + f*nmad
        #Use IQR
        #rmin = diff_stats[7]
        #rmax = diff_stats[8]
        diff_euler = np.ma.masked_outside(diff_euler, rmin, rmax)
        #Should also apply to original dem1 and dem2 for sad and ncc

    print("Computing sub-pixel offset between DEMs using mode: %s" % mode)

    #By default, don't create output figure
    fig = None

    #Sum of absolute differences
    if mode == "sad":
        m, int_offset, sp_offset = coreglib.compute_offset_sad(dem1, dem2, pad=pad)
        #Geotransform has negative y resolution, so don't need negative sign
        #np array is positive down
        #GDAL coordinates are positive up
        dx = sp_offset[1]*dem2_gt[1]
        dy = sp_offset[0]*dem2_gt[5]
    #Normalized cross-correlation of clipped, overlapping areas
    elif mode == "ncc":
        m, int_offset, sp_offset, fig = coreglib.compute_offset_ncc(dem1, dem2, \
                pad=pad, prefilter=False, plot=True)
        dx = sp_offset[1]*dem2_gt[1]
        dy = sp_offset[0]*dem2_gt[5]
    #Nuth and Kaab (2011)
    elif mode == "nuth":
        print("Computing slope and aspect")
        dem1_slope = geolib.gdaldem_mem_ds(dem1_clip_ds, processing='slope', returnma=True)
        dem1_aspect = geolib.gdaldem_mem_ds(dem1_clip_ds, processing='aspect', returnma=True)
        #Compute relationship between elevation difference, slope and aspect
        fit_param, fig = coreglib.compute_offset_nuth(diff_euler, dem1_slope, dem1_aspect)
        #fit_param[0] is magnitude of shift vector
        #fit_param[1] is direction of shift vector
        #fit_param[2] is mean bias divided by tangent of mean slope 
        #print(fit_param)
        dx = fit_param[0]*np.sin(np.deg2rad(fit_param[1]))
        dy = fit_param[0]*np.cos(np.deg2rad(fit_param[1]))
        #med_slope = malib.fast_median(dem1_slope)
        #dz = fit_param[2]*np.tan(np.deg2rad(med_slope))
    elif mode == "all":
        print("Not yet implemented")
        #Want to compare all methods, average offsets
        #m, int_offset, sp_offset = coreglib.compute_offset_sad(dem1, dem2)
        #m, int_offset, sp_offset = coreglib.compute_offset_ncc(dem1, dem2)
    #This is a hack to apply the computed median bias correction for shpclip area only
    elif mode == "none":
        print("Skipping alignment, writing out DEM with median bias over static surfaces removed")
        dst_fn = outprefix+'_med%0.1f.tif' % dz
        iolib.writeGTiff(dem2_orig + dz, dst_fn, dem2_ds)
        sys.exit()
    #Note: minus signs here since we are computing dz=(src-ref), but adjusting src 
    return -dx, -dy, -dz, static_mask, fig

#Defined a second main to allow recursion with new arguments for second run
def main2(args):
    #Should check that files exist
    dem1_fn = args.ref_fn
    dem2_fn = args.src_fn
    mode = args.mode
    apply_mask = not args.nomask
    filter = args.filter
    max_offset_m = args.max_offset
    tiltcorr = args.tiltcorr

    #These are tolerances (in meters) to stop iteration
    tol = args.tol
    min_dx = tol
    min_dy = tol
    min_dz = tol

    #Maximum number of iterations
    max_n = 10 
    
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
    
    dem2_ds = gdal.Open(dem2_fn, gdal.GA_ReadOnly)
    #Often the "ref" DEM is high-res lidar or similar
    #This is a shortcut to resample to match "source" DEM
    dem1_ds = warplib.memwarp_multi_fn([dem1_fn,], res=dem2_ds, extent=dem2_ds, t_srs=dem2_ds)[0]
    #dem1_ds = gdal.Open(dem1_fn, gdal.GA_ReadOnly)

    #Create a copy to be updated in place
    dem2_ds_align = iolib.mem_drv.CreateCopy('', dem2_ds, 0)
    #dem2_ds_align = dem2_ds

    #Iteration number
    n = 1
    #Cumulative offsets
    dx_total = 0
    dy_total = 0
    dz_total = 0

    #Now iteratively update geotransform and vertical shift
    while True:
        print("*** Iteration %i ***" % n)
        dx, dy, dz, static_mask, fig = compute_offset(dem1_ds, dem2_ds_align, dem2_fn, mode, max_offset_m, \
                apply_mask=apply_mask, filter=filter)
        if n == 1:
            static_mask_orig = static_mask
        xyz_shift_str_iter = "dx=%+0.2fm, dy=%+0.2fm, dz=%+0.2fm" % (dx, dy, dz)
        print("Incremental offset: %s" % xyz_shift_str_iter)

        #Should make an animation of this converging
        if fig is not None:
            dst_fn = outprefix + '_%s_iter%i_plot.png' % (mode, n)
            print("Writing offset plot: %s" % dst_fn)
            fig.gca().set_title(xyz_shift_str_iter)
            fig.savefig(dst_fn, dpi=300, bbox_inches='tight', pad_inches=0.1)

        #Apply the horizontal shift to the original dataset
        dem2_ds_align = coreglib.apply_xy_shift(dem2_ds_align, dx, dy, createcopy=False)
        dem2_ds_align = coreglib.apply_z_shift(dem2_ds_align, dz, createcopy=False)

        dx_total += dx
        dy_total += dy
        dz_total += dz
        print("Cumulative offset: dx=%+0.2fm, dy=%+0.2fm, dz=%+0.2fm" % (dx_total, dy_total, dz_total))

        #Fit plane to residuals and remove
        #Might be better to do this after converging
        """
        if tiltcorr:
            print("Applying planar tilt correction")
            gt = dem2_ds_align.GetGeoTransform()
            #Need to compute diff_euler here
            #Copy portions of compute_offset, create new function 
            vals, resid, coeff = geolib.ma_fitplane(diff_euler_align, gt, perc=(4, 96))
            dem2_ds_align = coreglib.apply_z_shift(dem2_ds_align, -vals, createcopy=False)
        """

        n += 1
        print("\n")
        #If magnitude of shift in all directions is less than tol
        #if n > max_n or (abs(dx) <= min_dx and abs(dy) <= min_dy and abs(dz) <= min_dz):
        #If magnitude of shift is less than tol
        dm = np.sqrt(dx**2 + dy**2 + dz**2)
        if n > max_n or dm < tol:
            break

    #String to append to output filenames
    xyz_shift_str_cum = '_%s_x%+0.2f_y%+0.2f_z%+0.2f' % (mode, dx_total, dy_total, dz_total)
    if tiltcorr: 
        xyz_shift_str_cum += "_tiltcorr"

    #Compute original elevation difference
    if True:
        dem1_clip_ds, dem2_clip_ds = warplib.memwarp_multi([dem1_ds, dem2_ds], \
                res='max', extent='intersection', t_srs=dem2_ds) 
        dem1_orig = iolib.ds_getma(dem1_clip_ds, 1)
        dem2_orig = iolib.ds_getma(dem2_clip_ds, 1)
        diff_euler_orig = dem2_orig - dem1_orig
        if not apply_mask:
            static_mask_orig = np.ma.getmaskarray(diff_euler_orig)
        diff_euler_orig_compressed = diff_euler_orig[~static_mask_orig]
        diff_euler_orig_stats = np.array(malib.print_stats(diff_euler_orig_compressed))

        #Write out original eulerian difference map
        print("Writing out original euler difference map for common intersection before alignment")
        dst_fn = outprefix + '_orig_dz_eul.tif' 
        iolib.writeGTiff(diff_euler_orig, dst_fn, dem1_clip_ds)
        
    #Compute final elevation difference
    if True:
        dem1_clip_ds_align, dem2_clip_ds_align = warplib.memwarp_multi([dem1_ds, dem2_ds_align], \
                res='max', extent='intersection', t_srs=dem2_ds_align) 
        dem1_align = iolib.ds_getma(dem1_clip_ds_align, 1)
        dem2_align = iolib.ds_getma(dem2_clip_ds_align, 1)
        diff_euler_align = dem2_align - dem1_align
        if not apply_mask:
            static_mask = np.ma.getmaskarray(diff_euler_align)
        diff_euler_align_compressed = diff_euler_align[~static_mask]
        diff_euler_align_stats = np.array(malib.print_stats(diff_euler_align_compressed))

        #Fit plane to residuals and remove
        if tiltcorr:
            print("Applying planar tilt correction")
            gt = dem1_clip_ds_align.GetGeoTransform()
            #Need to apply the mask here, so we're only fitting over static surfaces
            #Note that the origmask=False will compute vals for all x and y indices, which is what we want 
            #Should offer option for polynomial of arbitrary order
            #Also, want a better robust fit - maybe throw out more outliers
            vals, resid, coeff = geolib.ma_fitplane(np.ma.array(diff_euler_align, mask=static_mask), \
                    gt, perc=(12.5, 87.5), origmask=False)
            #Remove planar offset from difference map
            diff_euler_align -= vals
            #Remove planar offset from aligned dem2
            #Note: dimensions of ds and vals will be different as vals are computed for clipped intersection
            #Recompute planar offset for dem2_ds_align extent
            xgrid, ygrid = geolib.get_xy_grids(dem2_ds_align)
            vals = coeff[0]*xgrid + coeff[1]*ygrid + coeff[2] 
            dem2_ds_align = coreglib.apply_z_shift(dem2_ds_align, -vals, createcopy=False)
            if not apply_mask:
                static_mask = np.ma.getmaskarray(diff_euler_align)
            diff_euler_align_compressed = diff_euler_align[~static_mask]
            diff_euler_align_stats = np.array(malib.print_stats(diff_euler_align_compressed))
            print("Creating fitplane plot")
            fig, ax = plt.subplots(figsize=(6, 6))
            fitplane_clim = malib.calcperc(vals, (2,98))
            im = ax.imshow(vals, cmap='cpt_rainbow', clim=fitplane_clim)
            res = float(geolib.get_res(dem2_clip_ds, square=True)[0])
            pltlib.add_scalebar(ax, res=res)
            pltlib.hide_ticks(ax)
            pltlib.add_cbar(ax, im, label='Fit plane residuals (m)')
            fig.tight_layout()
            dst_fn1 = outprefix + '%s_align_dz_eul_fitplane.png' % xyz_shift_str_cum
            print("Writing out figure: %s" % dst_fn1)
            fig.savefig(dst_fn1, dpi=300, bbox_inches='tight', pad_inches=0.1)
   
        #Compute higher-order fits?
        #Could also attempt to model along-track and cross-track artifacts

        #Write out aligned eulerian difference map for clipped extent with vertial offset removed
        dst_fn = outprefix + '%s_align_dz_eul.tif' % xyz_shift_str_cum
        print("Writing out aligned difference map with median vertical offset removed") 
        iolib.writeGTiff(diff_euler_align, dst_fn, dem1_clip_ds) 
          
    #Write out aligned dem_2 with vertial offset removed
    if True:
        dst_fn2 = outprefix + '%s_align.tif' % xyz_shift_str_cum
        print("Writing out shifted dem2 with median vertical offset removed: %s" % dst_fn2)
        #Might be cleaner way to write out MEM ds directly to disk
        dem2_align = iolib.ds_getma(dem2_ds_align)
        iolib.writeGTiff(dem2_align, dst_fn2, dem2_ds_align) 
        dem2_ds_align = None

    #Create output plot
    if True:
        print("Creating final plot")
        dem1_hs = geolib.gdaldem_mem_ma(dem1_orig, dem1_clip_ds, returnma=True)
        dem2_hs = geolib.gdaldem_mem_ma(dem2_orig, dem2_clip_ds, returnma=True)
        f, axa = plt.subplots(2, 3, figsize=(11, 8.5))          
        for ax in axa.ravel()[:-1]:
            ax.set_facecolor('k')
            pltlib.hide_ticks(ax)
        dem_clim = malib.calcperc(dem1_orig, (2,98))
        axa[0,0].imshow(dem1_hs, cmap='gray')
        axa[0,0].imshow(dem1_orig, cmap='cpt_rainbow', clim=dem_clim, alpha=0.6)
        res = float(geolib.get_res(dem1_clip_ds, square=True)[0])
        pltlib.add_scalebar(axa[0,0], res=res)
        axa[0,0].set_title('Reference DEM')
        axa[0,1].imshow(dem2_hs, cmap='gray')
        axa[0,1].imshow(dem2_orig, cmap='cpt_rainbow', clim=dem_clim, alpha=0.6)
        axa[0,1].set_title('Source DEM')
        axa[0,2].imshow(~static_mask_orig, clim=(0,1), cmap='gray')
        axa[0,2].set_title('Surfaces for co-registration')
        dz_clim = malib.calcperc_sym(diff_euler_orig_compressed, (5, 95))
        im = axa[1,0].imshow(diff_euler_orig, cmap='RdBu', clim=dz_clim)
        pltlib.add_cbar(axa[1,0], im, label=None)
        axa[1,0].set_title('Elev. Diff. Before (m)')
        im = axa[1,1].imshow(diff_euler_align, cmap='RdBu', clim=dz_clim)
        pltlib.add_cbar(axa[1,1], im, label=None)
        axa[1,1].set_title('Elev. Diff. After (m)')

        #Tried to insert Nuth fig here
        #ax_nuth.change_geometry(1,2,1)
        #f.axes.append(ax_nuth)

        bins = np.linspace(dz_clim[0], dz_clim[1], 128)
        axa[1,2].hist(diff_euler_orig_compressed, bins, color='g', label='Before', alpha=0.5)
        axa[1,2].hist(diff_euler_align_compressed, bins, color='b', label='After', alpha=0.5)
        axa[1,2].axvline(0, color='k', linewidth=0.5, linestyle=':')
        axa[1,2].set_xlabel('Elev. Diff. (m)')
        axa[1,2].set_ylabel('Count (px)')
        axa[1,2].set_title("Source - Reference")
        #axa[1,2].legend(loc='upper right')
        #before_str = 'Before\nmean: %0.2f\nstd: %0.2f\nmed: %0.2f\nnmad: %0.2f' % tuple(diff_euler_orig_stats[np.array((3,4,5,6))])
        #after_str = 'After\nmean: %0.2f\nstd: %0.2f\nmed: %0.2f\nnmad: %0.2f' % tuple(diff_euler_align_stats[np.array((3,4,5,6))])
        before_str = 'Before\nmed: %0.2f\nnmad: %0.2f' % tuple(diff_euler_orig_stats[np.array((5,6))])
        axa[1,2].text(0.05, 0.95, before_str, va='top', color='g', transform=axa[1,2].transAxes) 
        after_str = 'After\nmed: %0.2f\nnmad: %0.2f' % tuple(diff_euler_align_stats[np.array((5,6))])
        axa[1,2].text(0.65, 0.95, after_str, va='top', color='b', transform=axa[1,2].transAxes) 

        suptitle = '%s\nx: %+0.2fm, y: %+0.2fm, z: %+0.2fm' % (os.path.split(outprefix)[-1], dx_total, dy_total, dz_total)
        f.suptitle(suptitle)
        f.tight_layout()
        plt.subplots_adjust(top=0.90)

        dst_fn = outprefix + '%s_align.png' % xyz_shift_str_cum
        print("Writing out figure: %s" % dst_fn)
        f.savefig(dst_fn, dpi=300, bbox_inches='tight', pad_inches=0.1)
       
        #Removing residual planar tilt can introduce additional slope/aspect dependent offset 
        #Want to run another round of main dem_align after removing planar tilt
        if tiltcorr:
            print("\n Rerunning after applying tilt correction \n")
            #Create copy of original arguments
            import copy
            args2 = copy.copy(args)
            #Use aligned, tilt-corrected DEM as input src_fn for second round
            args2.src_fn = dst_fn2
            #Assume we've already corrected most of the tilt during first round (also prevents endless loop)
            args2.tiltcorr = False
            main2(args2)

def main(argv=None):
    parser = getparser()
    args = parser.parse_args()
    main2(args)

if __name__ == "__main__":
    main()
