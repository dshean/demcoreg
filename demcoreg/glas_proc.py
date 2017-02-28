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

#parallel -j 32 '../glas_proc.py {}' ::: */*.H5
#cat */*conus_lulcfilt_demfilt.csv | sort -n | grep -v lat > GLAH14_tllz_conus_lulcfilt_demfilt.csv
#cat */*hma_lulcfilt_demfilt.csv | sort -n | grep -v lat > GLAH14_tllz_hma_lulcfilt_demfilt.csv
#clipsrc=/Volumes/d/hma/rgi/rgi_hma_aea_110kmbuffer_wgs84.shp
#vrt=GLAH14_tllz_hma_lulcfilt_demfilt.vrt
#ogr2ogr -progress -overwrite -clipsrc $clipsrc ${vrt%.*}_clip.shp $vrt

def writevrt(out_csv,srs='EPSG:4326'):
    out_vrt = os.path.splitext(out_csv)[0]+'.vrt'
    out_csv = os.path.split(out_csv)[-1]
    x = 'field_4'
    y = 'field_3'
    f = open(out_vrt, 'w')
    f.write('<OGRVRTDataSource>\n')
    f.write('   <OGRVRTLayer name="%s">\n' % os.path.splitext(out_csv)[0])
    f.write('        <SrcDataSource>%s</SrcDataSource>\n' % out_csv)
    f.write('        <GeometryType>wkbPoint</GeometryType>\n')
    f.write('        <LayerSRS>%s</LayerSRS>\n' % srs)
    f.write('        <GeometryField encoding="PointFromColumns" x="%s" y="%s"/>\n' % (x, y))
    f.write('    </OGRVRTLayer>\n')
    f.write('</OGRVRTDataSource>\n')
    f.close()

def ds_sample_coord(ds, x, y, xy_srs=geolib.wgs_srs):
    """Convert input coordinates to map coordinates of input dataset
    """
    #Convert lat/lon to projected srs
    mX, mY, mZ = geolib.cT_helper(x, y, 0, xy_srs, geolib.get_ds_srs(ds))
    return mX, mY

#This is arbitrary sampling function
#Assumes input map coords are identical to ds srs
def sample(ds, mX, mY, bn=1, pad=0, circ=False):
    """Sample input dataset at given coordinates
    """
    shape = (ds.RasterYSize, ds.RasterXSize)
    gt = ds.GetGeoTransform()
    b = ds.GetRasterBand(bn)
    b_ndv = iolib.get_ndv_b(b)
    b_dtype = b.DataType
    np_dtype = iolib.gdal2np_dtype(b)

    #This will sample an area corresponding to diameter of ICESat shot
    if pad == 'glas':
        spotsize = 70
        pad = int(np.ceil(((spotsize/gt[1])-1)/2))

    mX = np.atleast_1d(mX)
    mY = np.atleast_1d(mY)
    #Convert to pixel indices
    pX, pY = geolib.mapToPixel(mX, mY, gt)
    #Mask anything outside image dimensions
    pX = np.ma.masked_outside(pX, 0, shape[1]-1)
    pY = np.ma.masked_outside(pY, 0, shape[0]-1)
    common_mask = (~(np.logical_or(np.ma.getmaskarray(pX), np.ma.getmaskarray(pY)))).nonzero()[0]

    #Define x and y sample windows
    xwin=pad*2+1
    ywin=pad*2+1
    #This sets the minimum number of valid pixels, default 50%
    min_samp_perc = 50 
    min_samp = int(np.ceil((min_samp_perc/100.)*xwin*ywin))
    #Create circular mask to simulate spot
    #This only makes sense for for xwin > 3 
    if circ:
        circ_mask = filtlib.circular_mask(xwin)
        min_samp = int(np.ceil((min_samp_perc/100.)*circ_mask.nonzero()[0].size))

    pX_int = pX[common_mask].data
    pY_int = pY[common_mask].data
    #Round to nearest integer indices
    pX_int = np.around(pX_int).astype(int)
    pY_int = np.around(pY_int).astype(int)
    #print("Valid extent: %i" % pX_int.size)

    #Create empty array to hold output
    stats = np.full((pX_int.size, 2), b_ndv, dtype=np_dtype)

    r = gdal.GRA_NearestNeighbour
    #r = gdal.GRA_Cubic

    for i in range(pX_int.size):
        #Could have float offsets here with GDAL resampling
        samp = np.ma.masked_equal(b.ReadAsArray(xoff=pX_int[i]-pad, yoff=pY_int[i]-pad, win_xsize=xwin, win_ysize=ywin, resample_alg=r), b_ndv)
        if circ:
            samp = np.ma.array(samp, circ_mask)
        if samp.count() >= min_samp:
            if min_samp > 1:
                #samp_med = samp.mean()
                samp_med = malib.fast_median(samp)
                #samp_mad = samp.std()
                samp_mad = malib.mad(samp)
                stats[i][0] = samp_med 
                stats[i][1] = samp_mad
            else:
                stats[i][0] = samp[0]
                stats[i][1] = 0
            #vals, resid, coef = geolib.ma_fitplane(samp, gt=[0, gt[1], 0, 0, 0, gt[5]], perc=None)
            #Compute slope and aspect from plane
            #rmse = malib.rmse(resid)

    stats = np.ma.masked_equal(stats, b_ndv)
    #Create empty array with as input points
    out = np.full((pX.size, 2), b_ndv, dtype=np_dtype)
    #Populate with valid samples
    out[common_mask, :] = stats
    out = np.ma.masked_equal(out, b_ndv)
    return out

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

    #CONUS
    #xmin, xmax, ymin, ymax
    conus_extent = (-125, -104, 31, 50)
    extent = conus_extent
    name = 'conus'
    #NED 1/3 arcsec (10 m)
    dem_fn = '/nobackup/deshean/rpcdem/ned13/ned13_tiles_glac24k_115kmbuff.vrt'
    #NED 1 arcsec (30 m)
    #dem_fn = '/nobackup/deshean/rpcdem/ned1/ned1_tiles_glac24k_115kmbuff.vrt'
    #LULC
    lulc_fn = dem_mask.get_nlcd()

    #HMA
    hma_extent = (66, 106, 25, 47)
    extent = hma_extent
    name = 'hma'
    #SRTM-GL1 1 arcsec (30-m)
    dem_fn = '/nobackup/deshean/rpcdem/hma/srtm1/hma_srtm_gl1.vrt'
    lulc_fn = dem_mask.get_bareground()

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

    lulc_ds = gdal.Open(lulc_fn)
    print("Converting coords for LULC")
    lulc_mX, lulc_mY = ds_sample_coord(lulc_ds, out[:,2], out[:,1], geolib.wgs_srs)
    print("Sampling LULC")
    lulc_samp = sample(lulc_ds, lulc_mX, lulc_mY, pad=0)
    l = lulc_samp[:,0].data
    if 'nlcd' in lulc_fn:
        #l = l[:,np.newaxis]
        #This passes rock and ice pixels
        valid_idx *= np.logical_or((l==31),(l==12))
    else:
        minperc = 85
        valid_idx *= (l >= minperc)
    print("LULC: %i" % valid_idx.nonzero()[0].size)

    #Extract our own DEM values
    dem_ds = gdal.Open(dem_fn)
    print("Converting coords for DEM")
    dem_mX, dem_mY = ds_sample_coord(dem_ds, out[:,2], out[:,1], geolib.wgs_srs)
    print("Sampling DEM")
    dem_samp = sample(dem_ds, dem_mX, dem_mY, pad='glas')
    abs_dem_z_diff = np.abs(out[:,3] - dem_samp[:,0])

    valid_idx *= ~(np.ma.getmaskarray(abs_dem_z_diff))
    print("Valid DEM extract: %i" % valid_idx.nonzero()[0].size)
    valid_idx *= (abs_dem_z_diff < max_z_DEM_diff).data
    print("Valid abs DEM diff: %i" % valid_idx.nonzero()[0].size)
    valid_idx *= (dem_samp[:,1] < max_DEMhiresArElv_std).data
    print("Valid DEM mad: %i" % valid_idx.nonzero()[0].size)

    if valid_idx.nonzero()[0].size == 0:
        sys.exit("No valid points remain")

    if l.ndim == 1:
        l = l[:,np.newaxis]

    print("Writing out\n")
    dt_dyear = timelib.np_dt2decyear(timelib.o2dt(out[:,0]))[:,np.newaxis]
    out = np.ma.hstack([dt_dyear, out, dem_samp, l])
    out = out[valid_idx]
    out_fn = os.path.splitext(fn)[0]+'_tllz_%s_lulcfilt_demfilt.csv' % name
    #header = 't,lat,lon,elev'
    np.savetxt(out_fn, out, fmt='%0.8f, %0.10f, %0.6f, %0.6f, %0.2f, %0.2f, %0.2f, %i', delimiter=',')
    writevrt(out_fn)

if __name__ == "__main__":
    main()
