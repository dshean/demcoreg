#! /usr/bin/env python

"""
Compute difference between two rasters
"""

import sys
import os
import re
import argparse
from datetime import timedelta

import numpy as np

from pygeotools.lib import timelib, iolib, malib, warplib

def getparser():
    parser = argparse.ArgumentParser(description="Compute difference between two rasters")
    parser.add_argument('-outdir', default=None, help='Output directory')
    parser.add_argument('-tr', default='max', help='Output resolution (default: %(default)s)')
    parser.add_argument('-te', default='intersection', help='Output extent (default: %(default)s)')
    parser.add_argument('-t_srs', default='first', help='Output projection (default: %(default)s)')
    parser.add_argument('-rate', action='store_true', help='Attempt to generate elevation change rate products (dz/dt) in m/yr (default: %(default)s)')
    parser.add_argument('fn1', type=str, help='Raster filename 1')
    parser.add_argument('fn2', type=str, help='Raster filename 2')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    #This is output ndv, avoid using 0 for differences
    diffndv = -9999

    dem1_fn = args.fn1
    dem2_fn = args.fn2

    if dem1_fn == dem2_fn:
        sys.exit('Input filenames are identical')

    fn_list = [dem1_fn, dem2_fn]

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.dirname(os.path.abspath(dem1_fn))

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outprefix = os.path.splitext(os.path.split(dem1_fn)[1])[0]+'_'+os.path.splitext(os.path.split(dem2_fn)[1])[0]

    #Compute dz/dt rate if possible, in m/yr
    if args.rate:
        #Extract basename
        #This was a hack to work with timestamp array filenames that have geoid offset applied
        adj = ''
        if '-adj' in dem1_fn:
            adj = '-adj' 
        dem1_fn_base = re.sub(adj, '', os.path.splitext(dem1_fn)[0]) 
        dem2_fn_base = re.sub(adj, '', os.path.splitext(dem2_fn)[0]) 

        #Attempt to load ordinal timestamp arrays (for mosaics) if present
        t1_fn = dem1_fn_base+'_ts.tif'
        t2_fn = dem2_fn_base+'_ts.tif'
        t_unit = 'day'
        if not os.path.exists(t1_fn) and not os.path.exists(t2_fn):
            #Try to find processed output from dem_mosaic index
            #These are decimal years
            t1_fn = dem1_fn_base+'index_ts.tif'
            t2_fn = dem2_fn_base+'index_ts.tif'
            t_unit = 'year'
        if os.path.exists(t1_fn) and os.path.exists(t2_fn):
            fn_list.extend([t1_fn, t2_fn])
        else:
            #Attempt to extract timestamps from input filenames
            t1 = timelib.fn_getdatetime(dem1_fn)
            t2 = timelib.fn_getdatetime(dem2_fn)
            if t1 is not None and t2 is not None and t1 != t2:  
                dt = t2 - t1
                year = timedelta(days=365.25)
                t_factor = abs(dt.total_seconds()/year.total_seconds()) 
                print("Time differences is %s, dh/%0.3f" % (dt, t_factor))
            else:
                print("Unable to extract timestamps for input images")
                args.rate = False


    print("Warping DEMs to same res/extent/proj")
    #This will check input param for validity, could do beforehand
    ds_list = warplib.memwarp_multi_fn(fn_list, extent=args.te, res=args.tr, t_srs=args.t_srs)
    dem1_ds = ds_list[0]
    dem2_ds = ds_list[1]

    print("Loading input DEMs into masked arrays")
    dem1 = iolib.ds_getma(dem1_ds, 1)
    dem2 = iolib.ds_getma(dem2_ds, 1)

    #Compute relative elevation difference with Eulerian approach 
    print("Computing eulerian elevation difference")
    diff_euler = dem2 - dem1

    #Check to make sure inputs actually intersect
    #if not np.any(~dem1.mask*~dem2.mask):
    if diff_euler.count() == 0:
        sys.exit("No valid overlap between input DEMs")

    if len(fn_list) == 4:
        t1_ds = ds_list[2]
        t2_ds = ds_list[3]
        print("Loading timestamps into masked arrays")
        t1 = iolib.ds_getma(t1_ds)
        t2 = iolib.ds_getma(t2_ds)
        #Compute dt in years 
        t_factor = t2 - t1
        if t_unit == 'day':
            t_factor /= 365.25

    if True:
        print("Eulerian elevation difference stats:")
        diff_euler_stats = malib.print_stats(diff_euler)
        diff_euler_med = diff_euler_stats[5]

    if True:
        print("Writing Eulerian elevation difference map")
        dst_fn = os.path.join(outdir, outprefix+'_dz_eul.tif')
        print(dst_fn)
        iolib.writeGTiff(diff_euler, dst_fn, dem1_ds, ndv=diffndv)
        if args.rate:
            print("Writing Eulerian rate map")
            dst_fn = os.path.join(outdir, outprefix+'_dz_eul_rate.tif')
            print(dst_fn)
            iolib.writeGTiff(diff_euler/t_factor, dst_fn, dem1_ds, ndv=diffndv)
            if len(fn_list) == 4:
                print("Writing time difference map")
                dst_fn = os.path.join(outdir, outprefix+'_dz_eul_dt.tif')
                print(dst_fn)
                iolib.writeGTiff(t_factor, dst_fn, dem1_ds, ndv=diffndv)

    if False:
        print("Writing Eulerian relative elevation difference map")
        diff_euler_rel = diff_euler - diff_euler_med
        dst_fn = os.path.join(outdir, outprefix+'_dz_eul_rel.tif')
        print(dst_fn)
        iolib.writeGTiff(diff_euler_rel, dst_fn, dem1_ds, ndv=diffndv)

    if False:
        print("Writing out DEM2 with median elevation difference removed")
        dst_fn = os.path.splitext(dem2_fn)[0]+'_med'+diff_euler_med+'.tif'
        print(dst_fn)
        iolib.writeGTiff(dem2 - diff_euler_med, dst_fn, dem1_ds, ndv=diffndv)

    if False:
        print("Writing Eulerian elevation difference percentage map")
        diff_euler_perc = 100.0*diff_euler/dem1
        dst_fn = os.path.join(outdir, outprefix+'_dz_eul_perc.tif')
        print(dst_fn)
        iolib.writeGTiff(diff_euler_perc, dst_fn, dem1_ds, ndv=diffndv)

if __name__ == "__main__":
    main()
