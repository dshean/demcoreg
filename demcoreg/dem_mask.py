#! /usr/bin/env python
"""
Utility to automate reference surface identification for raster co-registration

Note: Initial run may take a long time to download and process required data (NLCD, global bareground, glacier polygons)

Can control location of these data files with DATADIR environmental variable

export DATADIR=dir

Dependencies: gdal, wget, requests, bs4

"""

#To do: 
#Integrate 1-km LULC data: http://www.landcover.org/data/landcover/
#TODO: need to clean up toa handling

import sys
import os
import subprocess
import glob
import argparse

from osgeo import gdal, ogr, osr
import numpy as np

from datetime import datetime, timedelta

from pygeotools.lib import iolib, warplib, geolib, timelib

datadir = iolib.get_datadir()

def get_nlcd_fn():
    """Calls external shell script `get_nlcd.sh` to fetch:

    2011 Land Use Land Cover (nlcd) grids, 30 m
    
    http://www.mrlc.gov/nlcd11_leg.php
    """
    #This is original filename, which requires ~17 GB
    #nlcd_fn = os.path.join(datadir, 'nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img')
    #get_nlcd.sh now creates a compressed GTiff, which is 1.1 GB
    nlcd_fn = os.path.join(datadir, 'nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.tif')
    if not os.path.exists(nlcd_fn):
        cmd = ['get_nlcd.sh',]
        #subprocess.call(cmd)
        sys.exit("Missing nlcd data source. If already downloaded, specify correct datadir. If not, run `%s` to download" % cmd[0])
    return nlcd_fn

def get_bareground_fn():
    """Calls external shell script `get_bareground.sh` to fetch:

    ~2010 global bare ground, 30 m

    Note: unzipped file size is 64 GB! Original products are uncompressed, and tiles are available globally (including empty data over ocean)

    The shell script will compress all downloaded tiles using lossless LZW compression.

    http://landcover.usgs.gov/glc/BareGroundDescriptionAndDownloads.php
    """
    bg_fn = os.path.join(datadir, 'bare2010/bare2010.vrt')
    if not os.path.exists(bg_fn):
        cmd = ['get_bareground.sh',]
        sys.exit("Missing bareground data source. If already downloaded, specify correct datadir. If not, run `%s` to download" % cmd[0])
        #subprocess.call(cmd)
    return bg_fn 

#Download latest global RGI glacier db
def get_glacier_poly():
    """Calls external shell script `get_rgi.sh` to fetch:

    Randolph Glacier Inventory (RGI) glacier outline shapefiles 

    Full RGI database: rgi50.zip is 410 MB

    The shell script will unzip and merge regional shp into single global shp
    
    http://www.glims.org/RGI/
    """
    #rgi_fn = os.path.join(datadir, 'rgi50/regions/rgi50_merge.shp')
    #Update to rgi60, should have this returned from get_rgi.sh
    rgi_fn = os.path.join(datadir, 'rgi60/regions/rgi60_merge.shp')
    if not os.path.exists(rgi_fn):
        cmd = ['get_rgi.sh',]
        sys.exit("Missing rgi glacier data source. If already downloaded, specify correct datadir. If not, run `%s` to download" % cmd[0])
        #subprocess.call(cmd)
    return rgi_fn 

#Update glacier polygons
def get_icemask(ds, glac_shp_fn=None):
    """Generate glacier polygon raster mask for input Dataset res/extent
    """
    print("Masking glaciers")
    if glac_shp_fn is None:
        glac_shp_fn = get_glacier_poly()

    if not os.path.exists(glac_shp_fn):
        print("Unable to locate glacier shp: %s" % glac_shp_fn)
    else:
        print("Found glacier shp: %s" % glac_shp_fn)

    #All of the proj, extent, handling should now occur in shp2array
    icemask = geolib.shp2array(glac_shp_fn, ds)
    return icemask

#Create nlcd mask 
def get_nlcd_mask(nlcd_ds, filter='not_forest', out_fn=None):
    """Generate raster mask for specified NLCD LULC filter
    """
    print("Loading NLCD LULC")
    b = nlcd_ds.GetRasterBand(1)
    l = b.ReadAsArray()
    print("Filtering NLCD LULC with: %s" % filter)
    #Original nlcd products have nan as ndv
        #12 - ice
        #31 - rock
        #11 - open water, includes rivers
        #52 - shrub, <5 m tall, >20%
        #42 - evergreeen forest
    #Should use data dictionary here for general masking
    #Using 'rock+ice+water' preserves the most pixels, although could be problematic over areas with lakes
    if filter == 'rock':
        mask = (l==31)
    elif filter == 'rock+ice':
        mask = np.logical_or((l==31),(l==12))
    elif filter == 'rock+ice+water':
        mask = np.logical_or(np.logical_or((l==31),(l==12)),(l==11))
    elif filter == 'not_forest':
        mask = ~(np.logical_or(np.logical_or((l==41),(l==42)),(l==43)))
    elif filter == 'not_forest+not_water':
        mask = ~(np.logical_or(np.logical_or(np.logical_or((l==41),(l==42)),(l==43)),(l==11)))
    else:
        print("Invalid mask type")
        mask = None
    #Write out original data
    if out_fn is not None:
        print("Writing out %s" % out_fn)
        iolib.writeGTiff(l, out_fn, nlcd_ds)
    l = None
    return mask

def get_bareground_mask(bareground_ds, bareground_thresh=60, out_fn=None):
    """Generate raster mask for exposed bare ground from global bareground data
    """
    print("Loading bareground")
    b = bareground_ds.GetRasterBand(1)
    l = b.ReadAsArray()
    print("Masking pixels with <%0.1f%% bare ground" % bareground_thresh)
    if bareground_thresh < 0.0 or bareground_thresh > 100.0:
        sys.exit("Invalid bare ground percentage")
    mask = (l>bareground_thresh)
    #Write out original data
    if out_fn is not None:
        print("Writing out %s" % out_fn)
        iolib.writeGTiff(l, out_fn, bareground_ds)
    l = None
    return mask

def get_snodas_ds(dem_dt, code=1036):
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

    outdir = os.path.join(datadir, 'snodas')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

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
    modis_dict = {}
    for key in modis_grid.modis_dict:
        modis_dict[key] = ogr.CreateGeometryFromWkt(modis_grid.modis_dict[key])
    geom = geolib.ds_geom(ds)
    geom_dup = geolib.geom_dup(geom)
    ct = osr.CoordinateTransformation(geom_dup.GetSpatialReference(), geolib.wgs_srs)
    geom_dup.Transform(ct)
    tile_list = []
    for key, val in list(modis_dict.items()):
        if geom_dup.Intersects(val):
            tile_list.append(key)
    return tile_list

def get_modscag_fn_list(dem_dt, tile_list=('h08v04', 'h09v04', 'h10v04', 'h08v05', 'h09v05'), pad_days=7):
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

    outdir = os.path.join(datadir, 'modscag')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    out_vrt_fn_list = []
    for dt in dt_list:
        out_vrt_fn = os.path.join(outdir, dt.strftime('%Y%m%d_snow_fraction.vrt'))
        #If we already have a vrt and it contains all of the necessary tiles
        if os.path.exists(out_vrt_fn):
            vrt_ds = gdal.Open(out_vrt_fn)
            if np.all([np.any([tile in sub_fn for sub_fn in vrt_ds.GetFileList()]) for tile in tile_list]):
                out_vrt_fn_list.append(out_vrt_fn)
                continue
        #Otherwise, download missing tiles and rebuild
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

def get_toa_fn(dem_fn):
    toa_fn = None
    #Original approach, assumes DEM file is in *00/dem_*/*DEM_32m.tif
    #dem_dir = os.path.split(os.path.split(os.path.abspath(dem_fn))[0])[0]
    dem_dir_list = os.path.split(os.path.realpath(dem_fn))[0].split(os.sep)
    import re
    #Get index of the top level pair directory containing toa (WV02_20140514_1030010031114100_1030010030896000)
    r_idx = [i for i, item in enumerate(dem_dir_list) if re.search('(_10)*(_10)*00$', item)]
    if r_idx:
        r_idx = r_idx[0]
        #Reconstruct dir
        dem_dir = (os.sep).join(dem_dir_list[0:r_idx+1])
        #Find toa.tif in top-level dir
        toa_fn = glob.glob(os.path.join(dem_dir, '*toa.tif'))
        if not toa_fn:
            ortho_fn = glob.glob(os.path.join(dem_dir, '*ortho*.tif'))
            if ortho_fn:
                cmd = ['toa.sh', dem_dir]
                print(cmd)
                subprocess.call(cmd)
                toa_fn = glob.glob(os.path.join(dem_dir, '*toa.tif'))
        if toa_fn:
            toa_fn = toa_fn[0]
        else:
            toa_fn = None
    if toa_fn is None:
        sys.exit("Unable to locate TOA dataset")
    return toa_fn

#TOA reflectance filter
def get_toa_mask(toa_ds, toa_thresh=0.4):
    print("Applying TOA filter (masking values >= %0.2f)" % toa_thresh)
    toa = iolib.ds_getma(toa_ds)
    toa_mask = np.ma.masked_greater(toa, toa_thresh)
    #This should be 1 for valid surfaces, 0 for snowcovered surfaces
    toa_mask = ~(np.ma.getmaskarray(toa_mask))
    return toa_mask

def check_mask_list(mask_list):
    temp = []
    for m in mask_list:
        if m not in mask_choices: 
            print("Invalid mask choice: %s" % m)
        else:
            temp.append(m)
    return temp

def get_mask(dem_ds, mask_list, dem_fn=None, writeout=False, outdir=None, args=None):
    mask_list = check_mask_list(mask_list)
    if 'none' in mask_list:
        newmask = False
    else:
        #Basename for output files
        if outdir is not None:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
        else:
            outdir = os.path.split(dem_fn)[0]

        if dem_fn is not None:
            #Extract DEM timestamp
            dem_dt = timelib.fn_getdatetime(dem_fn)
            out_fn_base = os.path.join(outdir, os.path.splitext(dem_fn)[0])
        
        if args is None:
            #Get default values
            parser = getparser()
            args = parser.parse_args(['',])

        newmask = True
        
        if 'glaciers' in mask_list:
            icemask = get_icemask(dem_ds)
            if writeout:
                out_fn = out_fn_base+'_ice_mask.tif'
                print("Writing out %s" % out_fn)
                iolib.writeGTiff(icemask, out_fn, src_ds=dem_ds)
            newmask = np.logical_and(icemask, newmask)

        #Need to process NLCD separately, with nearest neighbor inteprolatin
        if 'nlcd' in mask_list and args.nlcd_filter is not 'none':
            rs = 'near'
            nlcd_ds = gdal.Open(get_nlcd_fn())
            nlcd_ds_warp = warplib.memwarp_multi([nlcd_ds,], res=dem_ds, extent=dem_ds, t_srs=dem_ds, r=rs)[0]
            out_fn = None
            if writeout:
                out_fn = out_fn_base+'_nlcd.tif'
            nlcdmask = get_nlcd_mask(nlcd_ds_warp, filter=args.nlcd_filter, out_fn=out_fn)
            if writeout:
                out_fn = os.path.splitext(out_fn)[0]+'_mask.tif'
                print("Writing out %s" % out_fn)
                iolib.writeGTiff(nlcdmask, out_fn, src_ds=dem_ds)
            newmask = np.logical_and(nlcdmask, newmask)

        if 'bareground' in mask_list and args.bareground_thresh > 0:
            bareground_ds = gdal.Open(get_bareground_fn())
            bareground_ds_warp = warplib.memwarp_multi([bareground_ds,], res=dem_ds, extent=dem_ds, t_srs=dem_ds, r='cubicspline')[0]
            out_fn = None
            if writeout:
                out_fn = out_fn_base+'_bareground.tif'
            baregroundmask = get_bareground_mask(bareground_ds_warp, bareground_thresh=args.bareground_thresh, out_fn=out_fn)
            if writeout:
                out_fn = os.path.splitext(out_fn)[0]+'_mask.tif'
                print("Writing out %s" % out_fn)
                iolib.writeGTiff(baregroundmask, out_fn, src_ds=dem_ds)
            newmask = np.logical_and(baregroundmask, newmask)

        if 'snodas' in mask_list and args.snodas_thresh > 0:
            #Get SNODAS snow depth products for DEM timestamp
            snodas_min_dt = datetime(2003,9,30)
            if dem_dt >= snodas_min_dt: 
                snodas_ds = get_snodas_ds(dem_dt)
                if snodas_ds is not None:
                    snodas_ds_warp = warplib.memwarp_multi([snodas_ds,], res=dem_ds, extent=dem_ds, t_srs=dem_ds, r='cubicspline')[0]
                    #snow depth values are mm, convert to meters
                    snodas_depth = iolib.ds_getma(snodas_ds_warp)/1000.
                    if snodas_depth.count() > 0:
                        print("Applying SNODAS snow depth filter (masking values >= %0.2f m)" % args.snodas_thresh)
                        out_fn = None
                        if writeout:
                            out_fn = out_fn_base+'_snodas_depth.tif'
                            print("Writing out %s" % out_fn)
                            iolib.writeGTiff(snodas_depth, out_fn, src_ds=dem_ds)
                        snodas_mask = np.ma.masked_greater(snodas_depth, args.snodas_thresh)
                        snodas_mask = ~(np.ma.getmaskarray(snodas_mask))
                        if writeout:
                            out_fn = os.path.splitext(out_fn)[0]+'_mask.tif'
                            print("Writing out %s" % out_fn)
                            iolib.writeGTiff(snodas_mask, out_fn, src_ds=dem_ds)
                        newmask = np.logical_and(snodas_mask, newmask)
                    else:
                        print("SNODAS grid for input location and timestamp is empty")

        #These tiles cover CONUS
        #tile_list=('h08v04', 'h09v04', 'h10v04', 'h08v05', 'h09v05')
        if 'modscag' in mask_list and args.modscag_thresh > 0:
            modscag_min_dt = datetime(2000,2,24)
            if dem_dt < modscag_min_dt: 
                print("Warning: DEM timestamp (%s) is before earliest MODSCAG timestamp (%s)" \
                        % (dem_dt, modscag_min_dt))
            else:
                tile_list = get_modis_tile_list(dem_ds)
                print(tile_list)
                pad_days=7
                modscag_fn_list = get_modscag_fn_list(dem_dt, tile_list=tile_list, pad_days=pad_days)
                if modscag_fn_list:
                    modscag_ds = proc_modscag(modscag_fn_list, extent=dem_ds, t_srs=dem_ds)
                    modscag_ds_warp = warplib.memwarp_multi([modscag_ds,], res=dem_ds, extent=dem_ds, t_srs=dem_ds, r='cubicspline')[0]
                    print("Applying MODSCAG fractional snow cover percent filter (masking values >= %0.1f%%)" % args.modscag_thresh)
                    modscag_fsca = iolib.ds_getma(modscag_ds_warp)
                    out_fn = None
                    if writeout:
                        out_fn = out_fn_base+'_modscag_fsca.tif'
                        print("Writing out %s" % out_fn)
                        iolib.writeGTiff(modscag_fsca, out_fn, src_ds=dem_ds)
                    modscag_mask = (modscag_fsca.filled(0) >= args.modscag_thresh) 
                    modscag_mask = ~(modscag_mask)
                    if writeout:
                        out_fn = os.path.splitext(out_fn)[0]+'_mask.tif'
                        print("Writing out %s" % out_fn)
                        iolib.writeGTiff(modscag_mask, out_fn, src_ds=dem_ds)
                    newmask = np.logical_and(modscag_mask, newmask)

        #Use reflectance values to estimate snowcover
        if 'toa' in mask_list:
            #Use top of atmosphere scaled reflectance values (0-1)
            toa_ds = gdal.Open(get_toa_fn(dem_fn))
            toa_mask = get_toa_mask(toa_ds, args.toa_thresh)
            if writeout:
                out_fn = out_fn_base+'_toa_mask.tif'
                print("Writing out %s" % out_fn)
                iolib.writeGTiff(toa_mask, out_fn, src_ds=dem_ds)
            newmask = np.logical_and(toa_mask, newmask)

        if False:
            #Filter based on expected snowline
            #Simplest approach uses altitude cutoff
            max_elev = 1500 
            newdem = np.ma.masked_greater(dem, max_elev)
            newmask = np.ma.getmaskarray(newdem)

        print("Generating final mask to use for reference surfaces, and applying to input DEM")
        #Now invert to use to create final masked array
        #True (1) represents "invalid" pixel to match numpy ma convetion 
        newmask = ~newmask

        #Dilate the mask
        if args.dilate is not None:
            niter = args.dilate 
            print("Dilating mask with %i iterations" % niter)
            from scipy import ndimage
            newmask = ~(ndimage.morphology.binary_dilation(~newmask, iterations=niter))

    return newmask

#Can add "mask_list" argument, instead of specifying individually
mask_choices = ['toa', 'snodas', 'modscag', 'bareground', 'glaciers', 'nlcd', 'none']
def getparser():
    parser = argparse.ArgumentParser(description="Identify control surfaces for DEM co-registration") 
    parser.add_argument('dem_fn', type=str, help='DEM filename')
    parser.add_argument('--outdir', default=None, help='Directory for output products')
    parser.add_argument('--writeout', action='store_true', help='Write out all intermediate products, instead of only final tif')
    #parser.add_argument('-datadir', default=None, help='Data directory containing reference data sources (NLCD, bareground, etc)')
    parser.add_argument('--toa', action='store_true', help='Use top-of-atmosphere reflectance values (requires pregenerated "dem_fn_toa.tif")')
    parser.add_argument('--toa_thresh', type=float, default=0.4, help='Top-of-atmosphere reflectance threshold (default: %(default)s, valid range 0.0-1.0), mask values greater than this value')
    parser.add_argument('--snodas', action='store_true', help='Use SNODAS snow depth products')
    parser.add_argument('--snodas_thresh', type=float, default=0.2, help='SNODAS snow depth threshold (default: %(default)s m), mask values greater than this value')
    parser.add_argument('--modscag', action='store_true', help='Use MODSCAG fractional snow cover products')
    parser.add_argument('--modscag_thresh', type=float, default=50, help='MODSCAG fractional snow cover percent threshold (default: %(default)s%%, valid range 0-100), mask greater than this value')
    parser.add_argument('--bareground', action='store_true', help="Enable bareground filter")
    parser.add_argument('--bareground_thresh', type=float, default=60, help='Percent bareground threshold (default: %(default)s%%, valid range 0-100), mask greater than this value (only relevant for global bareground data)')
    parser.add_argument('--glaciers', action='store_true', help="Mask glacier polygons")
    parser.add_argument('--nlcd', action='store_true', help="Enable NLCD LULC filter (for CONUS)")
    nlcd_filter_choices = ['rock', 'rock+ice', 'rock+ice+water', 'not_forest', 'not_forest+not_water', 'none']
    parser.add_argument('--nlcd_filter', type=str, default='not_forest', choices=nlcd_filter_choices, help='Preserve these NLCD pixels (default: %(default)s)') 
    parser.add_argument('--dilate', type=int, default=None, help='Dilate mask with this many iterations (default: %(default)s)')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    mask_list = []
    if args.toa: mask_list.append('toa') 
    if args.snodas: mask_list.append('snodas') 
    if args.modscag: mask_list.append('modscag') 
    if args.bareground: mask_list.append('bareground') 
    if args.glaciers: mask_list.append('glaciers') 
    if args.nlcd: mask_list.append('nlcd') 

    if not mask_list:
        parser.print_help()
        sys.exit("Must specify at least one mask type")

    #This directory should or will contain the relevant data products
    #if args.datadir is None:
    #    datadir = iolib.get_datadir() 

    dem_fn = args.dem_fn
    dem_ds = gdal.Open(dem_fn)
    print(dem_fn)

    #Get DEM masked array
    dem = iolib.ds_getma(dem_ds)
    print("%i valid pixels in original input tif" % dem.count())

    #Set up cascading mask preparation
    #True (1) represents "valid" unmasked pixel, False (0) represents "invalid" pixel to be masked
    #Initialize the mask
    #newmask = ~(np.ma.getmaskarray(dem))
    
    newmask = get_mask(dem_ds, mask_list, dem_fn=dem_fn, writeout=args.writeout, outdir=args.outdir, args=args)
    
    #Apply mask to original DEM - use these surfaces for co-registration
    newdem = np.ma.array(dem, mask=newmask)

    #Check that we have enough pixels, good distribution
    min_validpx_count = 100
    min_validpx_std = 10
    validpx_count = newdem.count()
    validpx_std = newdem.std()
    print("%i valid pixels in masked output tif to be used as ref" % validpx_count)
    print("%0.2f std in masked output tif to be used as ref" % validpx_std)
    #if (validpx_count > min_validpx_count) and (validpx_std > min_validpx_std):
    if (validpx_count > min_validpx_count):
        out_fn = os.path.join(args.outdir, os.path.splitext(dem_fn)[0]+'_ref.tif')
        print("Writing out %s" % out_fn)
        iolib.writeGTiff(newdem, out_fn, src_ds=dem_ds)
    else:
        print("Not enough valid pixels!")


if __name__ == "__main__":
    main()
