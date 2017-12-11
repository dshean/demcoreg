#! /usr/bin/env python

"""
Utility to sample input raster for a set of input xyz points and print statistics
"""

import sys
import os
import csv    

from osgeo import gdal
import numpy as np
from pygeotools.lib import geolib, iolib, malib

r_fn = sys.argv[1]
pt_fn = sys.argv[2]

#In the future, allow user to specify point_srs, fields
#For now, assume output is from filter_glas.py, coords are lat/lon
xy_srs = geolib.wgs_srs
#Field numbers of x y z coordinates
xcol = 3
ycol = 2
zcol = 4

print("Loading points: %s" % pt_fn)
#This needs to return a header
pts = iolib.readcsv(pt_fn)

#Assume that all of our points are preselected to fall withing raster extent
#Can add option for spatial filter, see filter_glas.py

r_ds = gdal.Open(r_fn)

print("\nInput raster: %s" % r_fn)
print("Input points: %s" % pt_fn)
print("\nSampling %i points\n" % pts.shape[0])

#This returns median and mad arrays for all values within pad
#Should add option to return number of pixels in sample 
#Use 'glas' here to sample 70 m spot
samp = geolib.sample(r_ds, pts[:,xcol], pts[:,ycol], xy_srs=xy_srs, pad='glas', count=True)
samp_idx = ~(np.ma.getmaskarray(samp[:,0]))
nsamp = samp_idx.nonzero()[0].size

if nsamp == 0:
    sys.exit("No valid samples")
else:
    pts_mask = np.ma.hstack([pts[samp_idx], samp[samp_idx], (pts[samp_idx,zcol] - samp[samp_idx,0])[:,np.newaxis]])

    print("Sample difference (raster - point) statistics:")
    malib.print_stats(pts_mask[:,-1])

    """
    #Update header
    out_hdr_str = None
    if hdr is not None:
        out_hdr = hdr + ['samp_med','samp_nmad','samp_count']
        out_hdr_str = ', '.join(fieldnames)

    #Format for new cols
    #fmt = ['%0.2f', '%0.2f', '%i']
    #fmt = None
    """

    out_pt_fn = os.path.splitext(r_fn)[0]+'_'+os.path.splitext(os.path.split(pt_fn)[-1])[0]+'_sample.csv'
    print("\nWriting out %i points with valid samples:\n%s\n" % (nsamp, out_pt_fn))
    #np.savetxt(out_pt_fn, pts_mask, fmt=fmt, delimiter=',', header=out_hdr_str)
    np.savetxt(out_pt_fn, pts_mask, delimiter=',')

#Should migrate robust_stats functionality to malib
#robust_stats.py -col -1 out_pt_fn
