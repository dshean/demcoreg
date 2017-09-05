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

from pygeotools.lib import timelib 
from pygeotools.lib import iolib
from pygeotools.lib import malib
from pygeotools.lib import warplib

def getparser():
    parser = argparse.ArgumentParser(description="Compute difference between two rasters")
    parser.add_argument('-outdir', default=None, help='Output directory')
    parser.add_argument('-tr', default='max', help='Output resolution (default: %(default)s)')
    parser.add_argument('-te', default='intersection', help='Output extent (default: %(default)s)')
    parser.add_argument('-t_srs', default='first', help='Output projection (default: %(default)s)')
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

    print("Warping DEMs to same res/extent/proj")
    #This will check input param for validity, could do beforehand
    dem1_ds, dem2_ds = warplib.memwarp_multi_fn(fn_list, extent=args.te, res=args.tr, t_srs=args.t_srs)

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.split(dem1_fn)[0]
    outprefix = os.path.splitext(os.path.split(dem1_fn)[1])[0]+'_'+os.path.splitext(os.path.split(dem2_fn)[1])[0]

    print("Loading input DEMs into masked arrays")
    dem1 = iolib.ds_getma(dem1_ds, 1)
    dem2 = iolib.ds_getma(dem2_ds, 1)

    #Extract basename
    adj = ''
    if '-adj' in dem1_fn:
        adj = '-adj' 
    dem1_fn_base = re.sub(adj, '', os.path.splitext(dem1_fn)[0]) 
    dem2_fn_base = re.sub(adj, '', os.path.splitext(dem2_fn)[0]) 

    #Compute dz/dt rates if possible, in m/yr
    rates = True 
    if rates:
        #Attempt to load timestamp arrays (for mosaics) if present
        t1_fn = dem1_fn_base+'_ts.tif'
        t2_fn = dem2_fn_base+'_ts.tif'
        if os.path.exists(t1_fn) and os.path.exists(t2_fn):
            print("Preparing timestamp arrays")
            t1_ds, t2_ds = warplib.memwarp_multi_fn([t1_fn, t2_fn], extent=dem1_ds, res=dem1_ds)
            print("Loading timestamps into masked arrays")
            t1 = iolib.ds_getma(t1_ds)
            t2 = iolib.ds_getma(t2_ds)
            #Compute dt in days
            t_factor = t2 - t1
            t_factor /= 365.25
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
                rates = False

    #Check to make sure inputs actually intersect
    #Masked pixels are True
    if not np.any(~dem1.mask*~dem2.mask):
        sys.exit("No valid overlap between input data")

    #Compute common mask
    print("Generating common mask")
    common_mask = malib.common_mask([dem1, dem2])

    #Compute relative elevation difference with Eulerian approach 
    print("Computing eulerian elevation difference")
    diff_euler = np.ma.array(dem2-dem1, mask=common_mask)

    if True:
        print("Eulerian elevation difference stats:")
        diff_euler_stats = malib.print_stats(diff_euler)
        diff_euler_med = diff_euler_stats[5]

    if True:
        print("Writing Eulerian elevation difference map")
        dst_fn = os.path.join(outdir, outprefix+'_dz_eul.tif')
        print(dst_fn)
        iolib.writeGTiff(diff_euler, dst_fn, dem1_ds, ndv=diffndv)
        if rates:
            print("Writing Eulerian rate map")
            dst_fn = os.path.join(outdir, outprefix+'_dz_eul_rate.tif')
            print(dst_fn)
            iolib.writeGTiff(diff_euler/t_factor, dst_fn, dem1_ds, ndv=diffndv)

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
