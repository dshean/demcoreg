#! /usr/bin/env python
"""
Utility to automate reference surface identification for raster co-registration

Note: Initial run may take a long time to download and process required data (NLCD, global bareground, glacier polygons)

Can control location of these data files with DATADIR environmental variable

export DATADIR=dir

Dependencies: gdal, wget, requests, bs4

"""

#To do: 
#Add minium valid pixel count check - if not met, relax some criteria

import sys
import os
import subprocess
import glob
import argparse
from collections import OrderedDict

from osgeo import gdal, ogr, osr
import numpy as np

from datetime import datetime, timedelta

from pygeotools.lib import iolib
from pygeotools.lib import warplib
from pygeotools.lib import geolib
from pygeotools.lib import timelib

def get_nlcd_fn(datadir=None):
    """Calls external shell script `get_nlcd.sh` to fetch:

    2011 Land Use Land Cover (nlcd) grids, 30 m
    
    http://www.mrlc.gov/nlcd11_leg.php
    """
    if datadir is None:
        datadir = iolib.get_datadir()
    #This is original filename, which requires ~17 GB
    #nlcd_fn = os.path.join(datadir, 'nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img')
    #get_nlcd.sh now creates a compressed GTiff, which is 1.1 GB
    nlcd_fn = os.path.join(datadir, 'nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.tif')
    if not os.path.exists(nlcd_fn):
        cmd = ['get_nlcd.sh',]
        subprocess.call(cmd)
    return nlcd_fn

def get_bareground_fn(datadir=None):
    """Calls external shell script `get_bareground.sh` to fetch:

    ~2010 global bare ground, 30 m

    Note: unzipped file size is 64 GB! Original products are uncompressed, and tiles are available globally (including empty data over ocean)

    The shell script will compress all downloaded tiles using lossless LZW compression.

    http://landcover.usgs.gov/glc/BareGroundDescriptionAndDownloads.php
    """
    if datadir is None:
        datadir = iolib.get_datadir()
    bg_fn = os.path.join(datadir, 'bare2010/bare2010.vrt')
    if not os.path.exists(bg_fn):
        cmd = ['get_bareground.sh',]
        subprocess.call(cmd)
    return bg_fn 

#Download latest global RGI glacier db
def get_glacier_poly(datadir=None):
    """Calls external shell script `get_glacier_poly.sh` to fetch:

    Randolph Glacier Inventory (RGI) glacier outline shapefiles 

    Full RGI database: rgi50.zip is 410 MB

    The shell script will unzip and merge regional shp into single global shp
    
    http://www.glims.org/RGI/
    """
    if datadir is None:
        datadir = iolib.get_datadir()
    rgi_fn = os.path.join(datadir, 'rgi50/regions/rgi50_merge.shp')
    if not os.path.exists(rgi_fn):
        cmd = ['get_glacier_poly.sh',]
        subprocess.call(cmd)
    return rgi_fn 

#Update glacier polygons
def get_icemask(ds, datadir=None, glac_shp_fn=None):
    """Generate glacier polygon raster mask for input Dataset res/extent
    """
    if datadir is None:
        datadir = iolib.get_datadir()
    print("Masking glaciers")
    #Use updated glacier outlines to mask glaciers and perennial snowfields 
    #nlcd has an older glacier mask
    #Downloaded from http://www.glims.org/RGI/rgi50_files/02_rgi50_WesternCanadaUS.zip
    #ogr2ogr -t_srs EPSG:32610 02_rgi50_WesternCanadaUS_32610.shp 02_rgi50_WesternCanadaUS.shp
    #Manual selection over study area in QGIS
    if glac_shp_fn is None:
        glac_shp_dir = os.path.join(datadir, 'rgi50/regions')
        glac_shp_fn = os.path.join(glac_shp_dir, 'rgi50_merge.shp')

    if not os.path.exists(glac_shp_fn):
        print("Unable to locate glacier shp: %s" % glac_shp_fn)
    else:
        print("Found glacier shp: %s" % glac_shp_fn)

    #All of the proj, extent, handling should now occur in shp2array
    icemask = geolib.shp2array(glac_shp_fn, ds)
    return icemask

#Create rockmask from nlcd and remove glaciers
def mask_nlcd(ds, valid='rock+ice+water', datadir=None, mask_glaciers=True, out_fn=None):
    """Generate raster mask for exposed rock in NLCD data
    """
    print("Loading nlcd")
    b = ds.GetRasterBand(1)
    l = b.ReadAsArray()
    print("Isolating rock")
    #Original nlcd products have nan as ndv
        #12 - ice
        #31 - rock
        #11 - open water, includes rivers
        #52 - shrub, <5 m tall, >20%
        #42 - evergreeen forest
    #Should use data dictionary here for general masking
    #Using 'rock+ice+water' preserves the most pixels, although could be problematic over areas with lakes
    if valid == 'rock':
        mask = (l==31)
    elif valid == 'rock+ice':
        mask = np.logical_or((l==31),(l==12))
    elif valid == 'rock+ice+water':
        mask = np.logical_or(np.logical_or((l==31),(l==12)),(l==11))
    elif valid == 'not_forest':
        mask = ~(np.logical_or(np.logical_or((l==41),(l==42)),(l==43)))
    elif valid == 'not_forest+not_water':
        mask = ~(np.logical_or(np.logical_or(np.logical_or((l==41),(l==42)),(l==43)),(l==11)))
    else:
        print("Invalid mask type")
        mask = None
    #Write out original data
    if out_fn is not None:
        print("Writing out %s\n" % out_fn)
        iolib.writeGTiff(l, out_fn, ds)
    l = None
    if mask_glaciers:
        #Need better handling here, check to make sure 
        #Use updated 24k glacier outlines
        if datadir is None:
            datadir = iolib.get_datadir()
        glac_shp_fn = os.path.join(datadir, 'conus_glacierpoly_24k/conus_glacierpoly_24k_32610.shp')
        if not os.path.exists(glac_shp_fn):
            glac_shp_fn = None
        icemask = get_icemask(ds, datadir=datadir, glac_shp_fn=glac_shp_fn)
        if icemask is not None:
            mask *= icemask
    return mask

def mask_bareground(ds, minperc=80, mask_glaciers=True, out_fn=None):
    """Generate raster mask for exposed bare ground from global bareground data
    """
    print("Loading bareground")
    b = ds.GetRasterBand(1)
    l = b.ReadAsArray()
    print("Masking pixels with <%0.1f%% bare ground\n" % minperc)
    if minperc < 0.0 or minperc > 100.0:
        sys.exit("Invalid bare ground percentage")
    mask = (l>minperc)
    #Write out original data
    if out_fn is not None:
        print("Writing out %s\n" % out_fn)
        iolib.writeGTiff(l, out_fn, ds)
    l = None
    if mask_glaciers:
        icemask = get_icemask(ds)
        if icemask is not None:
            mask *= icemask
    return mask

def get_lulc_source(ds):
    """For a given input dataset extent, select the appropriate mask source (NLCD vs. global bareground)
    """
    ds_geom = geolib.ds_geom(ds)
    lulc_fn = get_nlcd_fn()
    lulc_ds = gdal.Open(lulc_fn)
    lulc_geom = geolib.ds_geom(lulc_ds)
    #If the dem geom is within CONUS (nlcd extent), use it
    geolib.geom_transform(ds_geom, t_srs=lulc_geom.GetSpatialReference())

    if lulc_geom.Contains(ds_geom):
        print("Using NLCD 30m data for rockmask")
        lulc_source = 'nlcd'
    else:
        print("Using global 30m bare ground data for rockmask")
        #Otherwise for Global, use 30 m Bare Earth percentage 
        lulc_source = 'bareground'
        #lulc_fn = get_bareground_fn()
        #lulc_ds = gdal.Open(lulc_fn)
    return lulc_source

def get_lulc_ds_full(ds, lulc_source=None):
    if lulc_source is None:
        lulc_source = get_lulc_source(ds)
    if lulc_source == 'nlcd':
        lulc_ds_full = gdal.Open(get_nlcd_fn())
    elif lulc_source == 'bareground':
        lulc_ds_full = gdal.Open(get_bareground_fn())
    return lulc_ds_full

def get_lulc_ds_warp(ds, lulc_source=None):
    if lulc_source is None:
        lulc_source = get_lulc_source(ds)
    if lulc_source == 'nlcd':
        #Note: want to process LULC with nearest to avoid interpolating values
        rs = 'near'
    else: 
        rs = 'cubicspline'
    lulc_ds_full = get_lulc_ds_full(ds, lulc_source)
    lulc_ds_warp = warplib.memwarp_multi([lulc_ds_full,], res=ds, extent=ds, t_srs=ds, r=rs)[0]
    return lulc_ds_warp

#Generate a rockmask
def get_lulc_mask(ds, lulc_source=None, mask_glaciers=True, filter='rock+ice+water', bareground_thresh=80, out_fn=None):
    if lulc_source is None:
        lulc_source = get_lulc_source(ds)
    if out_fn is not None:
        #Write out clipped, warped version of LULC for reference
        out_fn = os.path.splitext(out_fn)[0]+'_%s.tif' % lulc_source
    lulc_ds_warp = get_lulc_ds_warp(ds)
    #Note: these now have RGI 5.0 glacier polygons removed
    if lulc_source == 'nlcd':
        print("Applying NLCD LULC filter, preserving: %s" % filter)
        rockmask = mask_nlcd(lulc_ds_warp, valid=filter, mask_glaciers=mask_glaciers, out_fn=out_fn)
    elif lulc_source == 'bareground':
        print("Applying bareground percent filter (masking values >= %0.1f%%)" % bareground_thresh)
        rockmask = mask_bareground(lulc_ds_warp, minperc=bareground_thresh, mask_glaciers=mask_glaciers, out_fn=out_fn)
    return rockmask

def get_snodas(dem_dt, outdir=None, code=1036):
    """Function to fetch and process SNODAS snow depth products for input datetime

    http://nsidc.org/data/docs/noaa/g02158_snodas_snow_cover_model/index.html

    Product codes:
    1036 is snow depth
    1034 is SWE

    filename format: us_ssmv11036tS__T0001TTNATS2015042205HP001.Hdr

    """
    import tarfile
    import gzip
    snodas_ds = None
    snodas_url_str = None
    #Note: unmasked products (beyond CONUS) are only available from 2010-present
    if dem_dt >= datetime(2003,9,30) and dem_dt < datetime(2010,1,1):
        snodas_url_str = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked/%Y/%m_%b/SNODAS_%Y%m%d.tar'
        tar_subfn_str_fmt = 'us_ssmv1%itS__T0001TTNATS%%Y%%m%%d05HP001.%s.gz'
    elif dem_dt >= datetime(2010,1,1):
        snodas_url_str = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/unmasked/%Y/%m_%b/SNODAS_unmasked_%Y%m%d.tar'
        tar_subfn_str_fmt = './zz_ssmv1%itS__T0001TTNATS%%Y%%m%%d05HP001.%s.gz'
    else:
        print("No SNODAS data available for input date")

    if snodas_url_str is not None:
        snodas_url = dem_dt.strftime(snodas_url_str)
        snodas_tar_fn = iolib.getfile(snodas_url, outdir=outdir)
        print("Unpacking")
        tar = tarfile.open(snodas_tar_fn)
        #gunzip to extract both dat and Hdr files, tar.gz
        for ext in ('dat', 'Hdr'):
            tar_subfn_str = tar_subfn_str_fmt % (code, ext)
            tar_subfn_gz = dem_dt.strftime(tar_subfn_str)
            tar_subfn = os.path.splitext(tar_subfn_gz)[0]
            print(tar_subfn)
            if outdir is not None:
                tar_subfn = os.path.join(outdir, tar_subfn)
            if not os.path.exists(tar_subfn):
                #Should be able to do this without writing intermediate gz to disk
                tar.extract(tar_subfn_gz)
                with gzip.open(tar_subfn_gz, 'rb') as f:
                    outf = open(tar_subfn, 'wb')
                    outf.write(f.read())
                    outf.close()
                os.remove(tar_subfn_gz)

        #Need to delete 'Created by module comment' line from Hdr, can contain too many characters
        bad_str = 'Created by module comment'
        snodas_fn = tar_subfn
        f = open(snodas_fn)
        output = []
        for line in f:
            if not bad_str in line:
                output.append(line)
        f.close()
        f = open(snodas_fn, 'w')
        f.writelines(output)
        f.close()

        #Return GDAL dataset for extracted product
        snodas_ds = gdal.Open(snodas_fn)
    return snodas_ds

def get_modis_tile_list(ds):
    """Helper function to identify MODIS tiles that intersect input geometry

    modis_gird.py contains dictionary of tile boundaries (tile name and WKT polygon ring from bbox)

    See: https://modis-land.gsfc.nasa.gov/MODLAND_grid.html
    """
    from demcoreg import modis_grid
    modis_dict = modis_grid.modis_dict
    for key in modis_dict:
        modis_dict[key] = ogr.CreateGeometryFromWkt(modis_dict[key])
    geom = geolib.ds_geom(ds)
    geom_dup = geolib.geom_dup(geom)
    ct = osr.CoordinateTransformation(geom_dup.GetSpatialReference(), geolib.wgs_srs)
    geom_dup.Transform(ct)
    tile_list = []
    for key, val in list(modis_dict.items()):
        if geom_dup.Intersects(val):
            tile_list.append(key)
    return tile_list

def get_modscag(dem_dt, outdir=None, tile_list=('h08v04', 'h09v04', 'h10v04', 'h08v05', 'h09v05'), pad_days=7):
    """Function to fetch and process MODSCAG fractional snow cover products for input datetime

    Products are tiled in MODIS sinusoidal projection

    example url: https://snow-data.jpl.nasa.gov/modscag-historic/2015/001/MOD09GA.A2015001.h07v03.005.2015006001833.snow_fraction.tif

    """

    #Could also use global MODIS 500 m snowcover grids, 8 day
    #http://nsidc.org/data/docs/daac/modis_v5/mod10a2_modis_terra_snow_8-day_global_500m_grid.gd.html
    #These are HDF4, sinusoidal
    #Should be able to load up with warplib without issue

    import re
    import requests
    from bs4 import BeautifulSoup
    auth = iolib.get_auth()
    pad_days = timedelta(days=pad_days)
    dt_list = timelib.dt_range(dem_dt-pad_days, dem_dt+pad_days+timedelta(1), timedelta(1))
    out_vrt_fn_list = []
    for dt in dt_list:
        out_vrt_fn = os.path.join(outdir, dt.strftime('%Y%m%d_snow_fraction.vrt'))
        #If we already have a vrt and it contains all of the necessary tiles
        if os.path.exists(out_vrt_fn):
            vrt_ds = gdal.Open(out_vrt_fn)
            if np.all([np.any([tile in sub_fn for sub_fn in vrt_ds.GetFileList()]) for tile in tile_list]):
                out_vrt_fn_list.append(out_vrt_fn)
                continue
        #Otherwise, download missing tiles and rebuilt
        #Try to use historic products
        modscag_fn_list = []
        #Note: not all tiles are available for same date ranges in historic vs. real-time
        #Need to repeat search tile-by-tile
        for tile in tile_list:
            modscag_url_str = 'https://snow-data.jpl.nasa.gov/modscag-historic/%Y/%j/' 
            modscag_url_base = dt.strftime(modscag_url_str)
            print("Trying: %s" % modscag_url_base)
            r = requests.get(modscag_url_base, auth=auth)
            modscag_url_fn = []
            if r.ok:
                parsed_html = BeautifulSoup(r.content, "html.parser")
                modscag_url_fn = parsed_html.findAll(text=re.compile('%s.*snow_fraction.tif' % tile))
            if not modscag_url_fn:
                #Couldn't find historic, try to use real-time products
                modscag_url_str = 'https://snow-data.jpl.nasa.gov/modscag/%Y/%j/' 
                modscag_url_base = dt.strftime(modscag_url_str)
                print("Trying: %s" % modscag_url_base)
                r = requests.get(modscag_url_base, auth=auth)
            if r.ok: 
                parsed_html = BeautifulSoup(r.content, "html.parser")
                modscag_url_fn = parsed_html.findAll(text=re.compile('%s.*snow_fraction.tif' % tile))
            if not modscag_url_fn:
                print("Unable to fetch MODSCAG for %s" % dt)
            else:
                #OK, we got
                #Now extract actual tif filenames to fetch from html
                parsed_html = BeautifulSoup(r.content, "html.parser")
                #Fetch all tiles
                modscag_url_fn = parsed_html.findAll(text=re.compile('%s.*snow_fraction.tif' % tile))
                if modscag_url_fn:
                    modscag_url_fn = modscag_url_fn[0]
                    modscag_url = os.path.join(modscag_url_base, modscag_url_fn)
                    print(modscag_url)
                    modscag_fn = os.path.join(outdir, os.path.split(modscag_url_fn)[-1])
                    if not os.path.exists(modscag_fn):
                        iolib.getfile2(modscag_url, auth=auth, outdir=outdir)
                    modscag_fn_list.append(modscag_fn)
        #Mosaic tiles - currently a hack
        if modscag_fn_list:
            cmd = ['gdalbuildvrt', '-vrtnodata', '255', out_vrt_fn]
            cmd.extend(modscag_fn_list)
            print(cmd)
            subprocess.call(cmd, shell=False)
            out_vrt_fn_list.append(out_vrt_fn)
    return out_vrt_fn_list

def proc_modscag(fn_list, extent=None, t_srs=None):
    """Process the MODSCAG products for full date range, create composites and reproject
    """
    #Use cubic spline here for improve upsampling 
    ds_list = warplib.memwarp_multi_fn(fn_list, res='min', extent=extent, t_srs=t_srs, r='cubicspline')
    stack_fn = os.path.splitext(fn_list[0])[0] + '_' + os.path.splitext(os.path.split(fn_list[-1])[1])[0] + '_stack_%i' % len(fn_list) 
    #Create stack here - no need for most of mastack machinery, just make 3D array
    #Mask values greater than 100% (clouds, bad pixels, etc)
    ma_stack = np.ma.array([np.ma.masked_greater(iolib.ds_getma(ds), 100) for ds in np.array(ds_list)], dtype=np.uint8)

    stack_count = np.ma.masked_equal(ma_stack.count(axis=0), 0).astype(np.uint8)
    stack_count.set_fill_value(0)
    stack_min = ma_stack.min(axis=0).astype(np.uint8)
    stack_min.set_fill_value(0)
    stack_max = ma_stack.max(axis=0).astype(np.uint8)
    stack_max.set_fill_value(0)
    stack_med = np.ma.median(ma_stack, axis=0).astype(np.uint8)
    stack_med.set_fill_value(0)

    out_fn = stack_fn + '_count.tif'
    iolib.writeGTiff(stack_count, out_fn, ds_list[0])
    out_fn = stack_fn + '_max.tif'
    iolib.writeGTiff(stack_max, out_fn, ds_list[0])
    out_fn = stack_fn + '_min.tif'
    iolib.writeGTiff(stack_min, out_fn, ds_list[0])
    out_fn = stack_fn + '_med.tif'
    iolib.writeGTiff(stack_med, out_fn, ds_list[0])

    ds = gdal.Open(out_fn)
    return ds

def getparser():
    filter_choices = ['rock', 'rock+ice', 'rock+ice+water', 'not_forest', 'not_forest+not_water']
    parser = argparse.ArgumentParser(description="Identify control surfaces for DEM co-registration") 
    parser.add_argument('dem_fn', type=str, help='DEM filename')
    #parser.add_argument('-outdir', default=None, help='Output directory')
    parser.add_argument('--toa', action='store_true', help='Use top-of-atmosphere reflectance values (requires pregenerated "dem_fn_toa.tif")')
    parser.add_argument('--toa_thresh', type=float, default=0.4, help='Top-of-atmosphere reflectance threshold (default: %(default)s, valid range 0.0-1.0), mask values greater than this value')
    parser.add_argument('--snodas', action='store_true', help='Use SNODAS snow depth products')
    parser.add_argument('--snodas_thresh', type=float, default=0.2, help='SNODAS snow depth threshold (default: %(default)s m), mask values greater than this value')
    parser.add_argument('--modscag', action='store_true', help='Use MODSCAG fractional snow cover products')
    parser.add_argument('--modscag_thresh', type=float, default=50, help='MODSCAG fractional snow cover percent threshold (default: %(default)s%%, valid range 0-100), mask greater than this value')
    parser.add_argument('--bareground_thresh', type=float, default=80, help='Percent bareground threshold (default: %(default)s%%, valid range 0-100), mask greater than this value (only relevant for global bareground data)')
    parser.add_argument('--no_icemask', action='store_true', help="Don't mask glacier polygons")
    parser.add_argument('--filter', type=str, default='rock+ice+water', choices=filter_choices, help='Preserve these LULC pixels (default: %(default)s)') 
    parser.add_argument('--dilate', type=int, default=None, help='Dilate mask with this many iterations (default: %(default)s)')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    #Write out all mask products for the input DEM
    writeall = True

    mask_glaciers = True
    if args.no_icemask:
        mask_glaciers = False

    #Define top-level directory containing DEM
    topdir = os.getcwd()

    #This directory should contain nlcd grid, glacier outlines
    datadir = iolib.get_datadir() 

    dem_fn = args.dem_fn
    dem_ds = gdal.Open(dem_fn)
    print(dem_fn)

    #Extract DEM timestamp
    dem_dt = timelib.fn_getdatetime(dem_fn)

    #This will hold datasets for memwarp and output processing
    ds_dict = OrderedDict() 
    ds_dict['dem'] = dem_ds

    ds_dict['lulc'] = None
    #lulc_source = get_lulc_source(dem_ds)
    #lulc_ds_full = get_lulc_ds_full(dem_ds)
    #ds_dict['lulc'] = lulc_ds_full

    ds_dict['snodas'] = None
    if args.snodas:
        #Get SNODAS snow depth products for DEM timestamp
        snodas_min_dt = datetime(2003,9,30)
        if dem_dt >= snodas_min_dt: 
            snodas_outdir = os.path.join(datadir, 'snodas')
            if not os.path.exists(snodas_outdir):
                os.makedirs(snodas_outdir)
            snodas_ds = get_snodas(dem_dt, snodas_outdir)
            if snodas_ds is not None:
                ds_dict['snodas'] = snodas_ds

    ds_dict['modscag'] = None
    #Get MODSCAG products for DEM timestamp
    #These tiles cover CONUS
    #tile_list=('h08v04', 'h09v04', 'h10v04', 'h08v05', 'h09v05')
    if args.modscag:
        modscag_min_dt = datetime(2000,2,24)
        if dem_dt < modscag_min_dt: 
            print("\nWarning: DEM timestamp (%s) is before earliest MODSCAG timestamp (%s)\nSkipping..." \
                    % (dem_dt, modscag_min_dt))
        else:
            tile_list = get_modis_tile_list(dem_ds)
            print(tile_list)
            pad_days=7
            modscag_outdir = os.path.join(datadir, 'modscag')
            if not os.path.exists(modscag_outdir):
                os.makedirs(modscag_outdir)
            modscag_fn_list = get_modscag(dem_dt, modscag_outdir, tile_list, pad_days)
            if modscag_fn_list:
                modscag_ds = proc_modscag(modscag_fn_list, extent=dem_ds, t_srs=dem_ds)
                ds_dict['modscag'] = modscag_ds

    #TODO: need to clean this up
    #Better error handling
    #Disabled for now
    #Use reflectance values to estimate snowcover
    ds_dict['toa'] = None
    if args.toa:
        #Use top of atmosphere scaled reflectance values (0-1)
        dem_dir = os.path.split(os.path.split(os.path.abspath(dem_fn))[0])[0]
        toa_fn = glob.glob(os.path.join(dem_dir, '*toa.tif'))
        if not toa_fn:
            cmd = ['toa.sh', dem_dir]
            print(cmd)
            subprocess.call(cmd)
            toa_fn = glob.glob(os.path.join(dem_dir, '*toa.tif'))
        toa_fn = toa_fn[0]
        toa_ds = gdal.Open(toa_fn)
        ds_dict['toa'] = toa_ds  

    #Cull all of the None ds from the ds_dict
    for k,v in ds_dict.items():
        if v is None:
            del ds_dict[k]

    #Warp all masks to DEM extent/res
    #Note: use cubicspline here to avoid artifacts with negative values
    if len(ds_dict) > 0:
        ds_list = warplib.memwarp_multi(ds_dict.values(), res=dem_ds, extent=dem_ds, t_srs=dem_ds, r='cubicspline')
        #Update 
        for n, key in enumerate(ds_dict.keys()):
            ds_dict[key] = ds_list[n] 
    
    #lulc_ds_warp = get_lulc_ds_warp(dem_ds)
    #ds_dict['lulc'] = lulc_ds_warp

    print(' ')
    #Need better handling of ds order based on input ds here

    dem = iolib.ds_getma(ds_dict['dem'])
    newmask = ~(np.ma.getmaskarray(dem))

    #Generate a rockmask
    #Note: these now have RGI 5.0 glacier polygons removed
    #if 'lulc' in ds_dict.keys():
    #We are almost always going to want LULC mask
    out_fn_base = os.path.splitext(dem_fn)[0]
    rockmask = get_lulc_mask(dem_ds, mask_glaciers=mask_glaciers, \
            filter=args.filter, bareground_thresh=args.bareground_thresh, out_fn=out_fn_base)
    if writeall:
        out_fn = os.path.splitext(dem_fn)[0]+'_rockmask.tif'
        print("Writing out %s\n" % out_fn)
        iolib.writeGTiff(rockmask, out_fn, src_ds=ds_dict['dem'])
    newmask = np.logical_and(rockmask, newmask)

    if 'snodas' in ds_dict.keys():
        #SNODAS snow depth filter
        snodas_thresh = args.snodas_thresh 
        #snow depth values are mm, convert to meters
        snodas_depth = iolib.ds_getma(ds_dict['snodas'])/1000.
        if snodas_depth.count() > 0:
            print("Applying SNODAS snow depth filter (masking values >= %0.2f m)" % snodas_thresh)
            if writeall:
                out_fn = os.path.splitext(dem_fn)[0]+'_snodas_depth.tif'
                print("Writing out %s" % out_fn)
                iolib.writeGTiff(snodas_depth, out_fn, src_ds=ds_dict['dem'])
            snodas_mask = np.ma.masked_greater(snodas_depth, snodas_thresh)
            #This should be 1 for valid surfaces with no snow, 0 for snowcovered surfaces
            snodas_mask = ~(np.ma.getmaskarray(snodas_mask))
            if writeall:
                out_fn = os.path.splitext(dem_fn)[0]+'_snodas_mask.tif'
                print("Writing out %s\n" % out_fn)
                iolib.writeGTiff(snodas_mask, out_fn, src_ds=ds_dict['dem'])
            newmask = np.logical_and(snodas_mask, newmask)
        else:
            print("SNODAS grid for input location and timestamp is empty!\nSkipping...\n")

    if 'modscag' in ds_dict.keys():
        #MODSCAG percent snowcover
        modscag_thresh = args.modscag_thresh
        print("Applying MODSCAG fractional snow cover percent filter (masking values >= %0.1f%%)" % modscag_thresh)
        modscag_perc = iolib.ds_getma(ds_dict['modscag'])
        if writeall:
            out_fn = os.path.splitext(dem_fn)[0]+'_modscag_perc.tif'
            print("Writing out %s" % out_fn)
            iolib.writeGTiff(modscag_perc, out_fn, src_ds=ds_dict['dem'])
        modscag_mask = (modscag_perc.filled(0) >= modscag_thresh) 
        #This should be 1 for valid surfaces with no snow, 0 for snowcovered surfaces
        modscag_mask = ~(modscag_mask)
        if writeall:
            out_fn = os.path.splitext(dem_fn)[0]+'_modscag_mask.tif'
            print("Writing out %s\n" % out_fn)
            iolib.writeGTiff(modscag_mask, out_fn, src_ds=ds_dict['dem'])
        newmask = np.logical_and(modscag_mask, newmask)

    if 'toa' in ds_dict.keys():
        #TOA reflectance filter
        toa_thresh = args.toa_thresh
        print("Applying TOA filter (masking values >= %0.2f)" % toa_thresh)
        toa = iolib.ds_getma(ds_dict['toa'])
        toa_mask = np.ma.masked_greater(toa, toa_thresh)
        #This should be 1 for valid surfaces, 0 for snowcovered surfaces
        toa_mask = ~(np.ma.getmaskarray(toa_mask))
        if writeall:
            out_fn = os.path.splitext(dem_fn)[0]+'_toamask.tif'
            print("Writing out %s\n" % out_fn)
            iolib.writeGTiff(toa_mask, out_fn, src_ds=ds_dict['dem'])
        newmask = np.logical_and(toa_mask, newmask)

    if False:
        #Filter based on expected snowline
        #Simplest approach uses altitude cutoff
        max_elev = 1500 
        newdem = np.ma.masked_greater(dem, max_elev)
        newmask = np.ma.getmaskarray(newdem)

    print("Generating final mask to use for reference surfaces, and applying to input DEM")
    #Now invert to use to create final masked array
    newmask = ~newmask

    #Dilate the mask
    if args.dilate is not None:
        niter = args.dilate 
        print("Dilating mask with %i iterations" % niter)
        from scipy import ndimage
        newmask = ~(ndimage.morphology.binary_dilation(~newmask, iterations=niter))
    
    #Check that we have enough pixels, good distribution

    #Apply mask to original DEM - use these surfaces for co-registration
    newdem = np.ma.array(dem, mask=newmask)

    min_validpx_count = 100
    min_validpx_std = 10
    validpx_count = newdem.count()
    validpx_std = newdem.std()
    print("\n%i valid pixels in output ref.tif" % validpx_count)
    print("%0.2f m std output ref.tif\n" % validpx_std)
    #if (validpx_count > min_validpx_count) and (validpx_std > min_validpx_std):
    if (validpx_count > min_validpx_count):
        #Write out final mask
        out_fn = os.path.splitext(dem_fn)[0]+'_ref.tif'
        print("Writing out %s\n" % out_fn)
        iolib.writeGTiff(newdem, out_fn, src_ds=ds_dict['dem'])
    else:
        print("Not enough valid pixels!")

if __name__ == "__main__":
    main()
