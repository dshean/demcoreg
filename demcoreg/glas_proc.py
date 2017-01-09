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

#Download all GLAH14 products
#lftp ftp://n5eil01u.ecs.nsidc.org/DP5/GLAS/
#mirror --parallel=16 GLAH14.034

#Ben's script for processing: index_point_data_h5.m

#parallel -j 32 '../glas_proc.py {}' ::: */*.H5
#cat */*conus_demfilt.csv | sort -n | grep -v lat > GLAH14_tllz_conus_demfilt.csv
#cat */*hma_demfilt.csv | sort -n | grep -v lat > GLAH14_tllz_conus_hmafilt.csv

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

    """
    Need better args handling, for now hardcode
    if args.extent is not None:
        extent = (args.extent).split()
    if args.name is not None:
        name = args.name
    if args.dem_fn is not None:
        dem_fn = args.dem_fn
    """

    #CONUS
    #xmin, xmax, ymin, ymax
    #conus_extent = (-125, -104, 31, 50)
    #extent = conus_extent
    #name = 'conus'
    #NED 1/3 arcsec (10 m)
    #dem_fn = '/nobackup/deshean/rpcdem/ned13/ned13_tiles_glac24k_115kmbuff.vrt'
    #NED 1 arcsec (30 m)
    #dem_fn = '/nobackup/deshean/rpcdem/ned1/ned1_tiles_glac24k_115kmbuff.vrt'

    #HMA
    hma_extent = (66, 106, 25, 47)
    extent = hma_extent
    name = 'hma'
    #SRTM-GL1 1 arcsec (30-m)
    dem_fn = '/nobackup/deshean/rpcdem/hma/srtm1/hma_srtm_gl1.vrt'

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
    dts = timelib.np_print_dt(dt)

    lat = np.ma.masked_equal(f.get('Data_40HZ/Geolocation/d_lat')[:], 1.7976931348623157e+308)

    lon = np.ma.masked_equal(f.get('Data_40HZ/Geolocation/d_lon')[:], 1.7976931348623157e+308)
    lon = geolib.lon360to180(lon)

    z = np.ma.masked_equal(f.get('Data_40HZ/Elevation_Surfaces/d_elev')[:], 1.7976931348623157e+308)

    print('Input: %i' % z.count())

    #Saturation Correction Flag
    #These are 0 to 5, not_saturated inconsequential applicable not_computed not_applicable
    sat_corr_flg = f.get('Data_40HZ/Quality/sat_corr_flg')[:]
    #valid_idx *= (sat_corr_flg < 2)

    #Correction to elevation for saturated waveforms
    #Notes suggest this might not be desirable over land
    satElevCorr = np.ma.masked_equal(f.get('Data_40HZ/Elevation_Corrections/d_satElevCorr')[:], 1.7976931348623157e+308)
    #z[sat_corr_flg < 3] += satElevCorr.filled(0.0)[sat_corr_flg < 3]
    z += satElevCorr.filled(0.0)

    #Correction to elevation based on post flight analysis for biases determined for each campaign
    ElevBiasCorr = np.ma.masked_equal(f.get('Data_40HZ/Elevation_Corrections/d_ElevBiasCorr')[:], 1.7976931348623157e+308)
    z += ElevBiasCorr.filled(0.0)

    #Surface elevation (T/P ellipsoid) minus surface elevation (WGS84 ellipsoid).
    #Approximately 0.7 m, so WGS is lower; need to subtract from d_elev
    deltaEllip = np.ma.masked_equal(f.get('Data_40HZ/Geophysical/d_deltaEllip')[:], 1.7976931348623157e+308)
    z -= deltaEllip

    #These are 1 for valid, 0 for invalid
    valid_idx = ~(np.ma.getmaskarray(z))
    print("z corrections: %i" % valid_idx.nonzero()[0].size)

    #Reflectivity, not corrected for atmospheric effects
    reflctUC = np.ma.masked_equal(f.get('Data_40HZ/Reflectivity/d_reflctUC')[:], 1.7976931348623157e+308)
    #This was minimum used for ice sheets
    min_reflctUC = 0.025
    #valid_idx *= (reflctUC > min_reflctUC).data
    print("reflctUC: %i" % valid_idx.nonzero()[0].size)

    #The Standard deviation of the difference between the functional fit and the received echo using alternate parameters. It is directly taken from GLA05 parameter d_wfFitSDev_1
    LandVar = np.ma.masked_equal(f.get('Data_40HZ/Elevation_Surfaces/d_LandVar')[:], 1.7976931348623157e+308)
    #This was max used for ice sheets
    max_LandVar = 0.04
    #valid_idx *= (LandVar < max_LandVar)
    #print("LandVar: %i" % valid_idx.nonzero()[0].size)

    #Flag indicating whether the elevations on this record should be used.
    #0 = valid, 1 = not valid
    elev_use_flg = f.get('Data_40HZ/Quality/elev_use_flg')[:].astype('Bool')
    valid_idx *= ~elev_use_flg
    print("elev_use_flg: %i" % valid_idx.nonzero()[0].size)

    #Cloud contamination; Indicates if Gain > flag value, indicating probable cloud contamination.
    elv_cloud_flg = f.get('Data_40HZ/Elevation_Flags/elv_cloud_flg')[:].astype('Bool')
    #valid_idx *= ~elv_cloud_flg
    print("elv_cloud_flg: %i" % valid_idx.nonzero()[0].size)

    #Full resolution 1064 Quality Flag; 0 - 12 indicate Cloud detected
    FRir_qa_flg = f.get('Data_40HZ/Atmosphere/FRir_qa_flg')[:]
    #valid_idx *= (FRir_qa_flg == 15)

    #This is elevation extracted from SRTM30
    DEM_elv = np.ma.masked_equal(f.get('Data_40HZ/Geophysical/d_DEM_elv')[:], 1.7976931348623157e+308)
    z_DEM_diff = np.abs(z - DEM_elv)
    max_z_DEM_diff = 100
    #valid_idx *= (z_DEM_diff < max_z_DEM_diff).data
    print("z_DEM_diff: %i" % valid_idx.nonzero()[0].size)

    #d_DEMhiresArElv is a 9 element array of high resolution DEM values. The array index corresponds to the position of the DEM value relative to the spot. (5) is the footprint center.
    DEMhiresArElv = np.ma.masked_equal(f.get('Data_40HZ/Geophysical/d_DEMhiresArElv')[:], 1.7976931348623157e+308)
    DEMhiresArElv_std = np.ma.std(DEMhiresArElv, axis=1)
    max_DEMhiresArElv_std = 30.0
    #valid_idx *= (DEMhiresArElv_std < max_DEMhiresArElv_std).data
    print("max_DEMhiresArElv_std: %i" % valid_idx.nonzero()[0].size)
    #Compute slope

    a = np.ma.vstack([dto, lat, lon, z]).T
    mask = ~(np.any(a.mask, axis=1))
    mask *= ~valid_idx
    out = a[mask]

    #out_fn = os.path.splitext(fn)[0]+'_tllz.csv'
    #header = 't,lat,lon,elev'
    #np.savetxt(out_fn, out, fmt='%0.10f,%0.6f,%0.6f,%0.3f', delimiter=',', header=header)

    #Now spatial filter - should do this up front
    x = out[:,2]
    y = out[:,1]
    xmin, xmax, ymin, ymax = extent
    idx = ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax))
    out2 = out[idx]
    if out2.shape[0] == 0:
        sys.exit("No points within specified extent\n")

    out_fn = os.path.splitext(fn)[0]+'_%s_filt.csv' % name
    header = 't,lat,lon,elev'
    #np.savetxt(out_fn, out2, fmt='%0.10f,%0.6f,%0.6f,%0.3f', delimiter=',', header=header)

    #Extract our own DEM values
    dem_ds = gdal.Open(dem_fn)
    dem_shape = (dem_ds.RasterYSize, dem_ds.RasterXSize)
    dem_gt = dem_ds.GetGeoTransform()
    b = dem_ds.GetRasterBand(1)
    b_ndv = iolib.get_ndv_b(b)

    t = out2[:,0]
    lon = out2[:,2]
    lat = out2[:,1]
    z = out2[:,3]

    #Convert lat/lon to projected srs
    mX, mY, mZ = geolib.cT_helper(lon, lat, z, geolib.wgs_srs, geolib.get_ds_srs(dem_ds))
    #Convert to pixel indices
    pX, pY = geolib.mapToPixel(mX, mY, dem_gt)
    #Mask anything outside image dimensions
    pX = np.ma.masked_outside(pX, 0, dem_shape[1]-1)
    pY = np.ma.masked_outside(pY, 0, dem_shape[0]-1)
    common_mask = np.logical_or(np.ma.getmaskarray(pX), np.ma.getmaskarray(pY))

    #This is diameter of ICESat shot
    spotsize=70
    pad=int(np.ceil(((spotsize/dem_gt[1])-1)/2))
    xwin=pad*2+1
    ywin=pad*2+1
    #Create circular mask to simulate spot
    #This only makes sense for for xwin > 3 
    #circ_mask = filtlib.circular_mask(xwin)

    #pX_int = np.around(np.clip(pX-pad, 0, dem_shape[1]-1)).astype(int)
    #pY_int = np.around(np.clip(pY-pad, 0, dem_shape[0]-1)).astype(int)
    t_int = t[~common_mask].data
    z_int = t[~common_mask].data
    pX_int = pX[~common_mask].data
    pY_int = pY[~common_mask].data
    print("Valid DEM extent: %i" % pX_int.size)

    #Round to nearest integer indices
    pX_int = np.around(pX_int).astype(int)
    pY_int = np.around(pY_int).astype(int)
    #Create empty array to hold output
    dem_stats = np.full((pX_int.size, 2), -9999.0, dtype=np.float32)

    for i in range(pX_int.size):
        #Could have float offsets here with GDAL resampling
        #dem = np.ma.masked_equal(b.ReadAsArray(xoff=pX_int[i]-pad, yoff=pY_int[i]-pad, win_xsize=xwin, win_ysize=ywin, resample_alg=gdal.GRA_Cubic), b_ndv)
        dem = np.ma.masked_equal(b.ReadAsArray(xoff=pX_int[i]-pad, yoff=pY_int[i]-pad, win_xsize=xwin, win_ysize=ywin, resample_alg=gdal.GRA_NearestNeighbour), b_ndv)
        #dem = np.ma.array(dem, circ_mask)
        if dem.count() > 5:
            #v = [window_y[0]:window_y[1],window_x[0]:window_x[1]].reshape(m.shape[0], np.ptp(window_x)*np.ptp(window_y))
            vmed = malib.fast_median(dem)
            vmad = malib.mad(dem)
            dem_stats[i][0] = vmed
            dem_stats[i][1] = vmad
            #vals, resid, coef = geolib.ma_fitplane(dem, gt=[0, dem_gt[1], 0, 0, 0, dem_gt[5]], perc=None)
            #Compute slope and aspect from plane
            #rmse = malib.rmse(resid)

    dem_stats = np.ma.masked_equal(dem_stats, -9999)
    out2 = out2[~common_mask]
    abs_dem_z_diff = np.abs(out2[:,3] - dem_stats[:,0])

    dem_valid_idx = ~(np.ma.getmaskarray(abs_dem_z_diff))
    print("Valid DEM extract: %i" % dem_valid_idx.nonzero()[0].size)
    dem_valid_idx *= (abs_dem_z_diff < max_z_DEM_diff)
    print("Valid abs DEM diff: %i" % dem_valid_idx.nonzero()[0].size)
    dem_valid_idx *= (dem_stats[:,1] < max_DEMhiresArElv_std) 
    print("Valid DEM mad: %i\n" % dem_valid_idx.nonzero()[0].size)

    out3 = out2[dem_valid_idx]
    out_fn = os.path.splitext(fn)[0]+'_tllz_hma_demfilt.csv'
    #header = 't,lat,lon,elev'
    np.savetxt(out_fn, out3, fmt='%0.10f,%0.6f,%0.6f,%0.3f', delimiter=',')

if __name__ == "__main__":
    main()
