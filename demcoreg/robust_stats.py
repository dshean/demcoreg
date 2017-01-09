#! /usr/bin/env python

#David Shean
#dshean@gmail.com

import os
import sys
import argparse

import numpy as np

#This takes in the sampled output from pc_align and computes statistics
#Compute and print robust statistics for before or after samples (csv or tif) from pc_align_wrapper.sh

def getparser():
    parser = argparse.ArgumentParser(description="Compute robust stats from dz or pc_align sample.csv")
    parser.add_argument('fn', type=str, help='Input filename')
    parser.add_argument('-outdir', default=None, help='Output directory')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    fn = args.fn
    ext = os.path.splitext(fn)[1]

    #This expects the sample.csv, not the errors.csv
    if 'csv' in ext:
        a = np.loadtxt(fn, delimiter=',', skiprows=1)
        #Signed difference values are in column 5
        dz_m = a[:,4]
        dz_m = dz_m[~np.isnan(dz_m)]
    #If pc_align was run with reference grid, then load the dz raster
    elif 'tif' in ext: 
        from pygeotools.lib import iolib
        a = iolib.fn_getma(fn)
        dz_m = a.compressed()
    else:
        sys.exit('Unsupported input type')

    dz_m_abs = np.abs(dz_m)

    #Extract fn date
    #d = f[0:13]
    #print("Date: %s" % d)

    #print("Filename: %s" % f)

    print("Count: %i" % (dz_m.shape[0] - 1))

    rmse = np.sqrt(np.sum(dz_m**2)/dz_m.size)
    print("RMSE: %0.3f" % rmse)

    mean = np.mean(dz_m)
    print("Mean Error: %0.3f" % mean) 

    std = np.std(dz_m)
    print("Standard Deviation: %0.3f" % std) 

    #thresh = 3 * std
    med = np.median(dz_m)
    print("Median Error: %0.3f" % med) 

    p16, p84 = np.percentile(dz_m, (15.9, 84.2))
    spread = p84 - p16
    print("16th Percentile: %0.3f" % p16)
    print("84th Percentile: %0.3f" % p84)
    print("Spread: %0.3f" % spread)

    absmed = np.median(dz_m_abs)
    print("Absolute Median Error: %0.3f" % absmed) 

    mad = np.median(np.abs(dz_m - med))
    #print("MAD: %0.3f" % mad)
    nmad = 1.4826 * mad
    print("NMAD: %0.3f" % nmad)

    p68, p95 = np.percentile(dz_m_abs, (68.3, 95.0))
    print("68.3th Percentile: %0.3f" % p68)
    print("95th Percentile: %0.3f" % p95)

    #se = np.std(dz_m) / np.sqrt(dz_m.size)
    #print "Standard Error: %0.3f" % se 

if __name__ == "__main__":
    main()
