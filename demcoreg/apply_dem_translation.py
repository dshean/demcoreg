#! /usr/bin/env python

#David Shean
#dshean@gmail.com

#Take translation-only transform from ASP pc_align, update projected geotransform and apply z offset to DEM
#This offers significant performance gains over writing out new PC and running point2dem

import sys, os

import numpy as np
from osgeo import gdal, osr

from pygeotools.lib import geolib
from pygeotools.lib import iolib

#This should probably be moved to geolib
def get_proj_shift(src_c, src_shift, s_srs, t_srs, inv_trans=False):
    if s_srs.IsSame(t_srs):
        proj_shift = src_shift
    else:
        src_c_shift = src_c + src_shift
        src2proj = osr.CoordinateTransformation(s_srs, t_srs)
        proj_c = np.array(src2proj.TransformPoint(*src_c))
        proj_c_shift = np.array(src2proj.TransformPoint(*src_c_shift))
        if inv_trans:
            proj_shift = proj_c - proj_c_shift
        else:
            proj_shift = proj_c_shift - proj_c
    #Reduce unnecessary precision
    proj_shift = np.around(proj_shift, 3)
    return proj_shift

def main():
    gdal_opt = ['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=IF_SAFER']

    if len(sys.argv) != 3:
        sys.exit("Usage: %s dem.tif pc_align.log" % os.path.basename(sys.argv[0]))

    src_fn = sys.argv[1]
    log_fn = sys.argv[2]

    #Prepare input image
    src_ds = gdal.Open(src_fn)
    print("Loading input DEM")
    src_a = iolib.ds_getma(src_ds)
    #Overhead for np.ma may be unnecessary here, could just load as np.array
    #src_a = src_ds.GetRasterBand(1).ReadAsArray()
    src_ndv = src_a.fill_value 

    #Extract info from input log
    #These are coordinates of source points
    ecef_centroid_str = 'Centroid of source points (Cartesian, meters):'
    llz_centroid_str = 'Centroid of source points (lat,lon,z):'
    ecef_shift_str = 'Translation vector (Cartesian, meters):'
    ned_shift_str = 'Translation vector (North-East-Down, meters):'
    llz_shift_str = 'Translation vector (lat,lon,z):'

    #Need to extract from log to know how to compute translation
    #if ref is csv and src is dem, want to transform source_center + shift
    inv_trans = False 
    #if ref is dem and src is csv, want to inverse transform ref by shift applied at (source_center - shift)
    #inv_trans_str = '--save-inv-transformed-reference-points'
    inv_trans_str = 'trans_reference: true'

    log = open(log_fn)
    for line in log:
        if inv_trans_str in line:
            inv_trans = True
        if ecef_centroid_str in line:
            ecef_c = np.fromstring(line.split('Vector3')[1][1:-1], sep=',')
        if llz_centroid_str in line:
            llz_c = np.fromstring(line.split('Vector3')[1][1:-1], sep=',')
        if ecef_shift_str in line:
            ecef_shift = np.fromstring(line.split('Vector3')[1][1:-1], sep=',')
        if ned_shift_str in line:
            ned_shift = np.fromstring(line.split('Vector3')[1][1:-1], sep=',')
        if llz_shift_str in line:
            llz_shift = np.fromstring(line.split('Vector3')[1][1:-1], sep=',')
            break
    log.close()

    #Reorder lat,lon,z to lon,lat,z (x,y,z)
    i = [1, 0, 2]
    llz_c = llz_c[i]
    llz_shift = llz_shift[i]

    #s_srs = geolib.wgs_srs
    #src_c = llz_c 
    #src_shift = llz_shift

    s_srs = geolib.ecef_srs
    src_c = ecef_c
    src_shift = ecef_shift
    #src_c_shift = src_c + src_shift

    #Determine shift in original dataset coords
    t_srs = geolib.get_ds_srs(src_ds)
    if t_srs is None:
        sys.exit("Unable to determine src_ds proj")
    proj_shift = get_proj_shift(src_c, src_shift, s_srs, t_srs, inv_trans)

    #Should compute distances in local projection
    #Local orthographic
    ortho_proj4 = '+proj=ortho +ellps=WGS84 +lat_0=%f +lon_0=%f' % (llz_c[1], llz_c[0])
    ortho_srs = osr.SpatialReference()
    ortho_srs.ImportFromProj4(ortho_proj4)
    ortho_shift = get_proj_shift(src_c, src_shift, s_srs, ortho_srs, inv_trans)

    #Local stereographic
    stereo_proj4 = '+proj=stere +ellps=WGS84 +lat_0=%f +lon_0=%f' % (llz_c[1], llz_c[0])
    stereo_srs = osr.SpatialReference()
    stereo_srs.ImportFromProj4(stereo_proj4)
    stereo_shift = get_proj_shift(src_c, src_shift, s_srs, stereo_srs, inv_trans)

    print("ECEF shift", ecef_shift, np.linalg.norm(ecef_shift))
    print("local ned shift", ned_shift, np.linalg.norm(ecef_shift))
    print("local ortho shift", ortho_shift, np.linalg.norm(ortho_shift))
    print("local stereo shift", stereo_shift, np.linalg.norm(stereo_shift))
    print("proj shift", proj_shift, np.linalg.norm(proj_shift))

    #out_fmt = "VRT"
    out_fmt = "TIF"
    out_fn = os.path.splitext(src_fn)[0]+'_trans'
    if out_fmt == "VRT": 
        print("Writing vrt with scaled values")
        dst_fn = out_fn+'.vrt'
        driver = gdal.GetDriverByName("VRT")
        dst_ds = driver.CreateCopy(dst_fn, src_ds, 0)
        dst_ds.GetRasterBand(1).SetOffset(proj_shift[2])
    else:
        dst_fn = out_fn+'.tif'
        driver = gdal.GetDriverByName("GTiff")
        #Create might be faster here 
        print("Copying input dataset")
        dst_ds = driver.CreateCopy(dst_fn, src_ds, 0, options=gdal_opt)
        #Apply vertical shift
        dst_b = dst_ds.GetRasterBand(1)
        print("Writing out z-shifted band")
        dst_b.SetNoDataValue(float(src_ndv))
        dst_b.WriteArray(np.around((src_a + proj_shift[2]).filled(src_ndv), decimals=3))

    dst_gt = list(dst_ds.GetGeoTransform())
    #Apply horizontal shift directly to geotransform
    dst_gt[0] += proj_shift[0]
    dst_gt[3] += proj_shift[1] 
    dst_ds.SetGeoTransform(dst_gt)

    print("Writing out: ", dst_fn)
    dst_ds = None
    src_ds = None

if __name__ == "__main__":
    main()
