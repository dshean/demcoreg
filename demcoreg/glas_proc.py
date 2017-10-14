#! /usr/bin/env python

#David Shean
#dshean@gmail.com

#Utility to process ICESat-1 GLAS products, filter and clip to specified bounding box
#Input is HDF5 GLAH14
#https://nsidc.org/data/GLAH14/versions/34
#http://nsidc.org/data/docs/daac/glas_altimetry/data-dictionary-glah14.html

import os, sys
from datetime import datetime, timedelta
import argparse

import h5py
import numpy as np
from osgeo import gdal

from pygeotools.lib import timelib, geolib, iolib, malib, filtlib

#This is needed for LULC products
import dem_mask

#Before running, download all GLAH14 products
#lftp ftp://n5eil01u.ecs.nsidc.org/DP5/GLAS/
#mirror --parallel=16 GLAH14.034

"""
cd GLAH14.034
lfs setstripe -c 32 .
for site in conus hma
do
    parallel --progress --delay 1 -j 32 "~/src/demcoreg/demcoreg/glas_proc.py {} $site" ::: */*.H5
    #Combine output
    for ext in ${site}.csv ${site}_refdemfilt.csv ${site}_refdemfilt_lulcfilt.csv
    do
        first=$(ls */*$ext | head -1)
        head -1 $first > GLAH14_$ext
        cat */*$ext | sort -n | grep -v lat >> GLAH14_$ext
    done
done
"""

#Clip to glacier polygons
#clipsrc=/Volumes/d/hma/rgi/rgi_hma_aea_110kmbuffer_wgs84.shp
#vrt=GLAH14_tllz_hma_lulcfilt_demfilt.vrt
#ogr2ogr -progress -overwrite -clipsrc $clipsrc ${vrt%.*}_clip.shp $vrt

def getparser():
    parser = argparse.ArgumentParser(description="Process and filter ICESat GLAS points")
    parser.add_argument('fn', type=str, help='GLAH14 HDF5 filename')
    site_choices = geolib.site_dict.keys()
    parser.add_argument('sitename', type=str, choices=site_choices, help='Site name')
    #parser.add_argument('--rockfilter', action='store_true', help='Only output points over exposed rock using NLCD or bareground')
    parser.add_argument('-extent', type=str, default=None, help='Specify output spatial extent ("xmin xmax ymin ymax"). Otherwise, use default specified for sitename in pygeotools/lib/geolib')
    parser.add_argument('-refdem_fn', type=str, default=None, help='Specify alternative reference DEM for filtering. Otherwise use NED or SRTM')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    fn = args.fn
    sitename = args.sitename
    #User-specified output extent
    #Note: not checked, untested
    if args.extent is not None:
        extent = (args.extent).split()
    else:
        extent = (geolib.site_dict[sitename]).extent
    if args.refdem_fn is not None:
        refdem_fn = args.refdem_fn
    else:
        refdem_fn = (geolib.site_dict[sitename]).refdem_fn
    
    #Max elevation difference between shot and sampled DEM
    max_z_DEM_diff = 200
    #Max elevation std for sampled DEM values in padded window around shot
    max_DEMhiresArElv_std = 50.0

    f = h5py.File(fn)
    t = f.get('Data_40HZ/Time/d_UTCTime_40')[:]

    #pyt0 = datetime(1, 1, 1, 0, 0)
    #utct0 = datetime(1970, 1, 1, 0, 0)
    #t0 = datetime(2000, 1, 1, 12, 0, 0)
    #offset_s = (t0 - utct0).total_seconds()
    offset_s = 946728000.0
    t += offset_s
    dt = timelib.np_utc2dt(t)
    dt_o = timelib.dt2o(dt)
    #dts = timelib.np_print_dt(dt)
    #dt_decyear = timelib.np_dt2decyear(dt)
    dt_int = np.array([ts.strftime('%Y%m%d') for ts in dt], dtype=long)

    lat = np.ma.masked_equal(f.get('Data_40HZ/Geolocation/d_lat')[:], 1.7976931348623157e+308)
    lon = np.ma.masked_equal(f.get('Data_40HZ/Geolocation/d_lon')[:], 1.7976931348623157e+308)
    lon = geolib.lon360to180(lon)
    z = np.ma.masked_equal(f.get('Data_40HZ/Elevation_Surfaces/d_elev')[:], 1.7976931348623157e+308)

    print('Input: %i' % z.count())

    #Now spatial filter - should do this up front
    x = lon
    y = lat
    xmin, xmax, ymin, ymax = extent
    #This is True if point is within extent
    valid_idx = ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax))

    #Prepare output array
    #out = np.ma.vstack([dt_decyear, dt_o, dt_int, lat, lon, z]).T
    out = np.ma.vstack([dt_o, dt_int, lat, lon, z]).T
    #Create a mask to ensure all four values are valid for each point
    mask = ~(np.any(np.ma.getmaskarray(out), axis=1))
    mask *= valid_idx
    out = out[mask]
    valid_idx = ~(np.any(np.ma.getmaskarray(out), axis=1))

    #Lon and lat indices
    xcol = 3
    ycol = 2
    zcol = 4

    if out.shape[0] == 0:
        sys.exit("No points within specified extent\n")
    else:
        print("Spatial filter: %i" % out.shape[0])

    #out_fmt = ['%0.8f', '%0.8f', '%i', '%0.6f', '%0.6f', '%0.2f'] 
    #out_hdr = ['dt_decyear, dt_ordinal', 'dt_YYYYMMDD', 'lat', 'lon', 'z_WGS84']
    out_fmt = ['%0.8f', '%i', '%0.6f', '%0.6f', '%0.2f'] 
    out_hdr = ['dt_ordinal', 'dt_YYYYMMDD', 'lat', 'lon', 'z_WGS84']

    """
    ICESat-1 filters
    """
    #Saturation Correction Flag
    #These are 0 to 5, not_saturated inconsequential applicable not_computed not_applicable
    sat_corr_flg = f.get('Data_40HZ/Quality/sat_corr_flg')[mask]
    #valid_idx *= (sat_corr_flg < 2)

    #Correction to elevation for saturated waveforms
    #Notes suggest this might not be desirable over land
    satElevCorr = np.ma.masked_equal(f.get('Data_40HZ/Elevation_Corrections/d_satElevCorr')[mask], 1.7976931348623157e+308)
    #z[sat_corr_flg < 3] += satElevCorr.filled(0.0)[sat_corr_flg < 3]
    out[:,zcol] += satElevCorr.filled(0.0)

    #Correction to elevation based on post flight analysis for biases determined for each campaign
    ElevBiasCorr = np.ma.masked_equal(f.get('Data_40HZ/Elevation_Corrections/d_ElevBiasCorr')[mask], 1.7976931348623157e+308)
    out[:,zcol] += ElevBiasCorr.filled(0.0)

    #Surface elevation (T/P ellipsoid) minus surface elevation (WGS84 ellipsoid).
    #Approximately 0.7 m, so WGS is lower; need to subtract from d_elev
    deltaEllip = np.ma.masked_equal(f.get('Data_40HZ/Geophysical/d_deltaEllip')[mask], 1.7976931348623157e+308)
    out[:,zcol] -= deltaEllip

    #These are 1 for valid, 0 for invalid
    valid_idx *= ~(np.ma.getmaskarray(out[:,zcol]))
    print("z corrections: %i" % valid_idx.nonzero()[0].size)

    if False:
        #Reflectivity, not corrected for atmospheric effects
        reflctUC = np.ma.masked_equal(f.get('Data_40HZ/Reflectivity/d_reflctUC')[mask], 1.7976931348623157e+308)
        #This was minimum used for ice sheets
        min_reflctUC = 0.025
        valid_idx *= (reflctUC > min_reflctUC).data
        print("reflctUC: %i" % valid_idx.nonzero()[0].size)

    if False:
        #The Standard deviation of the difference between the functional fit and the received echo \
        #using alternate parameters. It is directly taken from GLA05 parameter d_wfFitSDev_1
        LandVar = np.ma.masked_equal(f.get('Data_40HZ/Elevation_Surfaces/d_LandVar')[mask], 1.7976931348623157e+308)
        #This was max used for ice sheets
        max_LandVar = 0.04
        valid_idx *= (LandVar < max_LandVar).data
        print("LandVar: %i" % valid_idx.nonzero()[0].size)

    if True:
        #Flag indicating whether the elevations on this record should be used.
        #0 = valid, 1 = not valid
        elev_use_flg = f.get('Data_40HZ/Quality/elev_use_flg')[mask].astype('Bool')
        valid_idx *= ~elev_use_flg
        print("elev_use_flg: %i" % valid_idx.nonzero()[0].size)

    if False:
        #Cloud contamination; Indicates if Gain > flag value, indicating probable cloud contamination.
        elv_cloud_flg = f.get('Data_40HZ/Elevation_Flags/elv_cloud_flg')[mask].astype('Bool')
        valid_idx *= ~elv_cloud_flg
        print("elv_cloud_flg: %i" % valid_idx.nonzero()[0].size)

    if False: 
        #Full resolution 1064 Quality Flag; 0 - 12 indicate Cloud detected
        FRir_qa_flg = f.get('Data_40HZ/Atmosphere/FRir_qa_flg')[mask]
        valid_idx *= (FRir_qa_flg == 15).data
        print("FRir_qa_flg: %i" % valid_idx.nonzero()[0].size)

    if False:
        #This is elevation extracted from SRTM30
        DEM_elv = np.ma.masked_equal(f.get('Data_40HZ/Geophysical/d_DEM_elv')[mask], 1.7976931348623157e+308)
        z_DEM_diff = np.abs(out[:,zcol] - DEM_elv)
        valid_idx *= (z_DEM_diff < max_z_DEM_diff).data
        print("z_DEM_diff: %i" % valid_idx.nonzero()[0].size)

        #d_DEMhiresArElv is a 9 element array of high resolution DEM values. The array index corresponds to the position of the DEM value relative to the spot. (5) is the footprint center.
        DEMhiresArElv = np.ma.masked_equal(f.get('Data_40HZ/Geophysical/d_DEMhiresArElv')[mask], 1.7976931348623157e+308)
        DEMhiresArElv_std = np.ma.std(DEMhiresArElv, axis=1)
        valid_idx *= (DEMhiresArElv_std < max_DEMhiresArElv_std).data
        print("max_DEMhiresArElv_std: %i" % valid_idx.nonzero()[0].size)
        #Compute slope

    #Apply cumulative filter to output
    out = out[valid_idx]

    out_fn = os.path.splitext(fn)[0]+'_%s.csv' % sitename
    print("Writing out %i records to: %s\n" % (out.shape[0], out_fn))
    out_fmt_str = ', '.join(out_fmt)
    out_hdr_str = ', '.join(out_hdr)
    np.savetxt(out_fn, out, fmt=out_fmt_str, delimiter=',', header=out_hdr_str)
    iolib.writevrt(out_fn, x='lon', y='lat')

    #Extract our own DEM values - should be better than default GLAS reference DEM stats
    if True:
        print("Loading reference DEM: %s" % refdem_fn)
        dem_ds = gdal.Open(refdem_fn)
        print("Converting coords for DEM")
        dem_mX, dem_mY = geolib.ds_cT(dem_ds, out[:,xcol], out[:,ycol], geolib.wgs_srs)
        print("Sampling")
        dem_samp = geolib.sample(dem_ds, dem_mX, dem_mY, pad='glas')
        abs_dem_z_diff = np.abs(out[:,zcol] - dem_samp[:,0])

        valid_idx *= ~(np.ma.getmaskarray(abs_dem_z_diff))
        print("Valid DEM extract: %i" % valid_idx.nonzero()[0].size)
        valid_idx *= (abs_dem_z_diff < max_z_DEM_diff).data
        print("Valid abs DEM diff: %i" % valid_idx.nonzero()[0].size)
        valid_idx *= (dem_samp[:,1] < max_DEMhiresArElv_std).data
        print("Valid DEM mad: %i" % valid_idx.nonzero()[0].size)

        if valid_idx.nonzero()[0].size == 0:
            sys.exit("No valid points remain")

        out = np.ma.hstack([out, dem_samp])
        out_fmt.extend(['%0.2f', '%0.2f'])
        out_hdr.extend(['z_refdem_med_WGS84', 'z_refdem_nmad'])

        #Apply cumulative filter to output
        out = out[valid_idx]

        out_fn = os.path.splitext(out_fn)[0]+'_refdemfilt.csv'
        print("Writing out %i records to: %s\n" % (out.shape[0], out_fn))
        out_fmt_str = ', '.join(out_fmt)
        out_hdr_str = ', '.join(out_hdr)
        np.savetxt(out_fn, out, fmt=out_fmt_str, delimiter=',', header=out_hdr_str)
        iolib.writevrt(out_fn, x='lon', y='lat')

    #This will sample land-use/land-cover or percent bareground products
    #Can be used to isolate points over exposed rock
    #if args.rockfilter: 
    if True:
        #This should automatically identify appropriate LULC source based on refdem extent
        lulc_source = dem_mask.get_lulc_source(dem_ds)
        #Looks like NED extends beyond NCLD, force use NLCD for conus
        #if sitename == 'conus':
        #    lulc_source = 'nlcd'
        lulc_ds = dem_mask.get_lulc_ds_full(dem_ds, lulc_source)
        print("Converting coords for LULC")
        lulc_mX, lulc_mY = geolib.ds_cT(lulc_ds, out[:,xcol], out[:,ycol], geolib.wgs_srs)
        print("Sampling LULC: %s" % lulc_source)
        #Note: want to make sure we're not interpolating integer values for NLCD
        #Should be safe with pad=0, even with pad>0, should take median, not mean
        lulc_samp = geolib.sample(lulc_ds, lulc_mX, lulc_mY, pad=0)
        l = lulc_samp[:,0].data
        if lulc_source == 'nlcd':
            #This passes rock and ice pixels
            valid_idx = np.logical_or((l==31),(l==12))
        elif lulc_source == 'bareground':
            #This preserves pixels with bareground percentation >85%
            minperc = 85
            valid_idx = (l >= minperc)
        else:
            print("Unknown LULC source")
        print("LULC: %i" % valid_idx.nonzero()[0].size)
        if l.ndim == 1:
            l = l[:,np.newaxis]
        out = np.ma.hstack([out, l])
        out_fmt.append('%i')
        out_hdr.append('lulc')

        #Apply cumulative filter to output
        out = out[valid_idx]

        out_fn = os.path.splitext(out_fn)[0]+'_lulcfilt.csv'
        print("Writing out %i records to: %s\n" % (out.shape[0], out_fn))
        out_fmt_str = ', '.join(out_fmt)
        out_hdr_str = ', '.join(out_hdr)
        np.savetxt(out_fn, out, fmt=out_fmt_str, delimiter=',', header=out_hdr_str)
        iolib.writevrt(out_fn, x='lon', y='lat')

if __name__ == "__main__":
    main()
