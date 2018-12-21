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
    parser.add_argument('-rate', action='store_true', help='Attempt to generate change rate products (dz/dt) in m/yr (default: %(default)s)')
    parser.add_argument('fn1', type=str, help='Raster filename 1')
    parser.add_argument('fn2', type=str, help='Raster filename 2')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    #This is output ndv, avoid using 0 for differences
    diffndv = -9999

    r1_fn = args.fn1
    r2_fn = args.fn2

    if r1_fn == r2_fn:
        sys.exit('Input filenames are identical')

    fn_list = [r1_fn, r2_fn]

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.dirname(os.path.abspath(r1_fn))

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outprefix = os.path.splitext(os.path.split(r1_fn)[1])[0]+'_'+os.path.splitext(os.path.split(r2_fn)[1])[0]

    #Compute dz/dt rate if possible, in m/yr
    if args.rate:
        #Extract basename
        #This was a hack to work with timestamp array filenames that have geoid offset applied
        adj = ''
        if '-adj' in r1_fn:
            adj = '-adj' 
        r1_fn_base = re.sub(adj, '', os.path.splitext(r1_fn)[0]) 
        r2_fn_base = re.sub(adj, '', os.path.splitext(r2_fn)[0]) 

        #Attempt to load ordinal timestamp arrays (for mosaics) if present
        """
        import glob
        t1_fn = glob.glob(r1_fn_base+'*_ts*.tif')
        t2_fn = glob.glob(r1_fn_base+'*_ts*.tif')
        t_unit = 'year'
        """
        t1_fn = r1_fn_base+'_ts.tif'
        t2_fn = r2_fn_base+'_ts.tif'
        t_unit = 'day'
        if not os.path.exists(t1_fn) and not os.path.exists(t2_fn):
            #Try to find processed output from r_mosaic index
            #These are decimal years
            t1_fn = r1_fn_base+'index_ts.tif'
            t2_fn = r2_fn_base+'index_ts.tif'
            t_unit = 'year'
        print(t1_fn, t2_fn)
        if os.path.exists(t1_fn) and os.path.exists(t2_fn):
            fn_list.extend([t1_fn, t2_fn])
        else:
            #Attempt to extract timestamps from input filenames
            t1 = timelib.fn_getdatetime(r1_fn)
            t2 = timelib.fn_getdatetime(r2_fn)
            if t1 is not None and t2 is not None and t1 != t2:  
                dt = t2 - t1
                year = timedelta(days=365.25)
                t_factor = abs(dt.total_seconds()/year.total_seconds()) 
                print("Time differences is %s, dh/%0.3f" % (dt, t_factor))
            else:
                print("Unable to extract timestamps for input images")
                args.rate = False


    print("Warping rasters to same res/extent/proj")
    #This will check input param for validity, could do beforehand
    ds_list = warplib.memwarp_multi_fn(fn_list, extent=args.te, res=args.tr, t_srs=args.t_srs, r='cubic')
    r1_ds = ds_list[0]
    r2_ds = ds_list[1]

    print("Loading input rasters into masked arrays")
    r1 = iolib.ds_getma(r1_ds, 1)
    r2 = iolib.ds_getma(r2_ds, 1)

    #Compute relative difference 
    print("Computing raster difference")
    diff = r2 - r1

    #Check to make sure inputs actually intersect
    if diff.count() == 0:
        sys.exit("No valid overlap between input rasters")

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
        print("Raster difference stats:")
        diff_stats = malib.print_stats(diff)
        diff_med = diff_stats[5]

    if True:
        print("Writing raster difference map")
        dst_fn = os.path.join(outdir, outprefix+'_diff.tif')
        print(dst_fn)
        iolib.writeGTiff(diff, dst_fn, r1_ds, ndv=diffndv)
        if args.rate:
            print("Writing rate map")
            dst_fn = os.path.join(outdir, outprefix+'_diff_rate.tif')
            print(dst_fn)
            iolib.writeGTiff(diff/t_factor, dst_fn, r1_ds, ndv=diffndv)
            if len(fn_list) == 4:
                print("Writing time difference map")
                dst_fn = os.path.join(outdir, outprefix+'_diff_dt.tif')
                print(dst_fn)
                iolib.writeGTiff(t_factor, dst_fn, r1_ds, ndv=diffndv)

    if False:
        print("Writing relative raster difference map")
        diff_rel = diff - diff_med
        dst_fn = os.path.join(outdir, outprefix+'_diff_rel.tif')
        print(dst_fn)
        iolib.writeGTiff(diff_rel, dst_fn, r1_ds, ndv=diffndv)

    if False:
        print("Writing out raster2 with median difference removed")
        dst_fn = os.path.splitext(r2_fn)[0]+'_med'+diff_med+'.tif'
        print(dst_fn)
        iolib.writeGTiff(r2 - diff_med, dst_fn, r1_ds, ndv=diffndv)

    if False:
        print("Writing raster difference percentage map (relative to raster1)")
        diff_perc = 100.0*diff/r1
        dst_fn = os.path.join(outdir, outprefix+'_diff_perc.tif')
        print(dst_fn)
        iolib.writeGTiff(diff_perc, dst_fn, r1_ds, ndv=diffndv)

if __name__ == "__main__":
    main()
