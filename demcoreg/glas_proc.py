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

#Download all GLAH14 products
#lftp ftp://n5eil01u.ecs.nsidc.org/DP5/GLAS/
#mirror --parallel=16 GLAH14.034

#Ben's script for processing: index_point_data_h5.m

#cd /nobackupp8/deshean/icesat_glas/GLAH14.034
#lfs setstripe -c 32 .
#parallel --progress --delay 1 -j 32 '~/src/demcoreg/demcoreg/glas_proc.py {}' ::: */*.H5
#cat */*conus_lulcfilt_demfilt.csv | sort -n | grep -v lat > GLAH14_tllz_conus_lulcfilt_demfilt.csv
#cat */*hma_lulcfilt_demfilt.csv | sort -n | grep -v lat > GLAH14_tllz_hma_lulcfilt_demfilt.csv
#clipsrc=/Volumes/d/hma/rgi/rgi_hma_aea_110kmbuffer_wgs84.shp
#vrt=GLAH14_tllz_hma_lulcfilt_demfilt.vrt
#ogr2ogr -progress -overwrite -clipsrc $clipsrc ${vrt%.*}_clip.shp $vrt

def getparser():
    parser = argparse.ArgumentParser(description="Process and filter ICESat GLAS points")
    parser.add_argument('fn', type=str, help='GLAH14 HDF5 filename')
    parser.add_argument('-extent', type=str, default=None, help='Output spatial extent')
    parser.add_argument('-name', type=str, default=None, help='Output file prefix')
    parser.add_argument('-dem_fn', type=str, default=None, help='Output file prefix')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    fn = args.fn

    #Need better args handling, for now, these are hardcoded below
    """
    if args.extent is not None:
        extent = (args.extent).split()
    if args.name is not None:
        name = args.name
    if args.dem_fn is not None:
        dem_fn = args.dem_fn
    """

    #Max elevation difference between shot and sampled DEM
    max_z_DEM_diff = 200
    #Max elevation std for sampled DEM values in padded window around shot
    max_DEMhiresArElv_std = 50.0

    #name = 'hma'
    name = 'conus'

    if name == 'conus':
        #CONUS
        #xmin, xmax, ymin, ymax
        extent = (-125, -104, 31, 50)
        #NED 1/3 arcsec (10 m)
        dem_fn = '/nobackup/deshean/rpcdem/ned13/ned13_tiles_glac24k_115kmbuff.vrt'
        #NED 1 arcsec (30 m)
        #dem_fn = '/nobackup/deshean/rpcdem/ned1/ned1_tiles_glac24k_115kmbuff.vrt'
        #LULC
        lulc_fn = dem_mask.get_nlcd_fn()
    elif name == 'hma':
        #HMA
        extent = (66, 106, 25, 47)
        #SRTM-GL1 1 arcsec (30-m)
        dem_fn = '/nobackup/deshean/rpcdem/hma/srtm1/hma_srtm_gl1.vrt'
        lulc_fn = dem_mask.get_bareground_fn()
    else:
        sys.exit("Other sites not supported at this time")

    f = h5py.File(fn)
    t = f.get('Data_40HZ/Time/d_UTCTime_40')[:]

    #pyt0 = datetime(1, 1, 1, 0, 0)
    #utct0 = datetime(1970, 1, 1, 0, 0)
    #t0 = datetime(2000, 1, 1, 12, 0, 0)
    #offset_s = (t0 - utct0).total_seconds()
    offset_s = 946728000.0
    t += offset_s
    dt = timelib.np_utc2dt(t)
    dto = timelib.dt2o(dt)
    #dts = timelib.np_print_dt(dt)

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

    out = np.ma.vstack([dto, lat, lon, z]).T
    #These is a single mask, where True, all input variables are valid
    mask = ~(np.any(np.ma.getmaskarray(out), axis=1))
    mask *= valid_idx
    out = out[mask]
    valid_idx = ~(np.any(np.ma.getmaskarray(out), axis=1))

    if out.shape[0] == 0:
        sys.exit("No points within specified extent\n")
    else:
        print("Spatial filter: %i" % out.shape[0])

    dto = out[:,0]
    x = out[:,1]
    y = out[:,2]
    z = out[:,3]

    #Saturation Correction Flag
    #These are 0 to 5, not_saturated inconsequential applicable not_computed not_applicable
    sat_corr_flg = f.get('Data_40HZ/Quality/sat_corr_flg')[mask]
    #valid_idx *= (sat_corr_flg < 2)

    #Correction to elevation for saturated waveforms
    #Notes suggest this might not be desirable over land
    satElevCorr = np.ma.masked_equal(f.get('Data_40HZ/Elevation_Corrections/d_satElevCorr')[mask], 1.7976931348623157e+308)
    #z[sat_corr_flg < 3] += satElevCorr.filled(0.0)[sat_corr_flg < 3]
    z += satElevCorr.filled(0.0)

    #Correction to elevation based on post flight analysis for biases determined for each campaign
    ElevBiasCorr = np.ma.masked_equal(f.get('Data_40HZ/Elevation_Corrections/d_ElevBiasCorr')[mask], 1.7976931348623157e+308)
    z += ElevBiasCorr.filled(0.0)

    #Surface elevation (T/P ellipsoid) minus surface elevation (WGS84 ellipsoid).
    #Approximately 0.7 m, so WGS is lower; need to subtract from d_elev
    deltaEllip = np.ma.masked_equal(f.get('Data_40HZ/Geophysical/d_deltaEllip')[mask], 1.7976931348623157e+308)
    z -= deltaEllip

    #These are 1 for valid, 0 for invalid
    valid_idx *= ~(np.ma.getmaskarray(z))
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
        z_DEM_diff = np.abs(z - DEM_elv)
        valid_idx *= (z_DEM_diff < max_z_DEM_diff).data
        print("z_DEM_diff: %i" % valid_idx.nonzero()[0].size)

        #d_DEMhiresArElv is a 9 element array of high resolution DEM values. The array index corresponds to the position of the DEM value relative to the spot. (5) is the footprint center.
        DEMhiresArElv = np.ma.masked_equal(f.get('Data_40HZ/Geophysical/d_DEMhiresArElv')[mask], 1.7976931348623157e+308)
        DEMhiresArElv_std = np.ma.std(DEMhiresArElv, axis=1)
        valid_idx *= (DEMhiresArElv_std < max_DEMhiresArElv_std).data
        print("max_DEMhiresArElv_std: %i" % valid_idx.nonzero()[0].size)
        #Compute slope

    #out = np.ma.array(out, mask=~(valid_idx))
    #valid_idx = ~(np.any(np.ma.getmaskarray(out), axis=1))

    if False:
        lulc_ds = gdal.Open(lulc_fn)
        print("Converting coords for LULC")
        lulc_mX, lulc_mY = geolib.ds_cT(lulc_ds, out[:,2], out[:,1], geolib.wgs_srs)
        print("Sampling LULC")
        lulc_samp = geolib.sample(lulc_ds, lulc_mX, lulc_mY, pad=0)
        l = lulc_samp[:,0].data
        if 'nlcd' in lulc_fn:
            #l = l[:,np.newaxis]
            #This passes rock and ice pixels
            valid_idx *= np.logical_or((l==31),(l==12))
        else:
            minperc = 85
            valid_idx *= (l >= minperc)
        print("LULC: %i" % valid_idx.nonzero()[0].size)
        if l.ndim == 1:
            l = l[:,np.newaxis]


    #Extract our own DEM values
    dem_ds = gdal.Open(dem_fn)
    print("Converting coords for DEM")
    dem_mX, dem_mY = geolib.ds_cT(dem_ds, out[:,2], out[:,1], geolib.wgs_srs)
    print("Sampling DEM")
    dem_samp = geolib.sample(dem_ds, dem_mX, dem_mY, pad='glas')
    abs_dem_z_diff = np.abs(out[:,3] - dem_samp[:,0])

    valid_idx *= ~(np.ma.getmaskarray(abs_dem_z_diff))
    print("Valid DEM extract: %i" % valid_idx.nonzero()[0].size)
    valid_idx *= (abs_dem_z_diff < max_z_DEM_diff).data
    print("Valid abs DEM diff: %i" % valid_idx.nonzero()[0].size)
    valid_idx *= (dem_samp[:,1] < max_DEMhiresArElv_std).data
    print("Valid DEM mad: %i" % valid_idx.nonzero()[0].size)

    if valid_idx.nonzero()[0].size == 0:
        sys.exit("No valid points remain")

    print("Writing out\n")
    dt_dyear = timelib.np_dt2decyear(timelib.o2dt(out[:,0]))[:,np.newaxis]
    if False:
        out = np.ma.hstack([dt_dyear, out, dem_samp, l])
        fmt = '%0.8f, %0.10f, %0.6f, %0.6f, %0.2f, %0.2f, %0.2f, %i'
    else:
        out = np.ma.hstack([dt_dyear, out, dem_samp])
        fmt = '%0.8f, %0.10f, %0.6f, %0.6f, %0.2f, %0.2f, %0.2f'
    out = out[valid_idx]
    #out_fn = os.path.splitext(fn)[0]+'_tllz_%s_lulcfilt_demfilt.csv' % name
    out_fn = os.path.splitext(fn)[0]+'_tllz_%s_demfilt.csv' % name
    #header = 't,lat,lon,elev'
    np.savetxt(out_fn, out, fmt=fmt, delimiter=',')
    iolib.writevrt(out_fn, x='field_4', y='field_3')

if __name__ == "__main__":
    main()
