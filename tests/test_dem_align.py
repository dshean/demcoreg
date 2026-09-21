"""Planted shift recovery tests for dem_align.py, synthetic terrain, no external data"""
import os
import glob
import json

import numpy as np
import pytest
import matplotlib
matplotlib.use('Agg')
from osgeo import gdal, osr

from demcoreg import dem_align

def make_dem(fn, xres=10., yres=10., size_m=4000., dx=0., dy=0., dz=0.):
    #Analytic terrain evaluated at true map coordinates, so a planted shift is exact
    ns = int(size_m/xres)
    nl = int(size_m/yres)
    ulx, uly = 500000., 4100000.
    x = ulx + (np.arange(ns) + 0.5)*xres
    y = uly - (np.arange(nl) + 0.5)*yres
    xx, yy = np.meshgrid(x - dx, y - dy)
    z = 1500 + 300*np.sin(xx/500.)*np.cos(yy/700.) + 100*np.sin(xx/170. + yy/230.)
    #Ridges and valleys, so there are sharp features for ncc
    z += 150*np.abs(np.sin(xx/410.))*np.abs(np.cos(yy/290.))
    ds = gdal.GetDriverByName('GTiff').Create(fn, ns, nl, 1, gdal.GDT_Float32)
    ds.SetGeoTransform([ulx, xres, 0, uly, 0, -yres])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32610)
    ds.SetProjection(srs.ExportToWkt())
    b = ds.GetRasterBand(1)
    b.SetNoDataValue(-9999)
    b.WriteArray((z + dz).astype(np.float32))
    ds = None

#Surface of src is displaced by (dx, dy, dz) relative to ref
#dem_align reports the shift to apply to src, so expect the opposite sign
planted = (13., -7., 2.5)

def run_dem_align(tmp_path, mode, xres=10., yres=10., extra_args=[]):
    ref_fn = str(tmp_path / 'ref.tif')
    src_fn = str(tmp_path / 'src.tif')
    make_dem(ref_fn)
    make_dem(src_fn, xres=xres, yres=yres, dx=planted[0], dy=planted[1], dz=planted[2])
    outdir = str(tmp_path / 'out')
    dem_align.main(['-mode', mode, '-mask_list', 'none', '-outdir', outdir] + extra_args + [ref_fn, src_fn])
    stats = json.load(open(glob.glob(os.path.join(outdir, '*_align_stats.json'))[0]))
    return stats['shift'], outdir, src_fn

@pytest.mark.parametrize('xres,yres', [(10., 10.), (9., 11.)])
def test_nuth_planted_shift(tmp_path, xres, yres):
    shift, outdir, src_fn = run_dem_align(tmp_path, 'nuth', xres, yres)
    assert shift['dx'] == pytest.approx(-planted[0], abs=0.1)
    assert shift['dy'] == pytest.approx(-planted[1], abs=0.1)
    assert shift['dz'] == pytest.approx(-planted[2], abs=0.05)
    #Shifted output and filtered output share the src grid (issue #69)
    src = gdal.Open(src_fn)
    for suffix in ('*_align.tif', '*_align_filt.tif'):
        ds = gdal.Open(glob.glob(os.path.join(outdir, suffix))[0])
        assert (ds.RasterYSize, ds.RasterXSize) == (src.RasterYSize, src.RasterXSize)

def test_ncc_planted_shift(tmp_path):
    #Parabolic sub-pixel peak, expect agreement to a fraction of a 10 m pixel
    shift, _, _ = run_dem_align(tmp_path, 'ncc', extra_args=['-max_iter', '10'])
    assert shift['dx'] == pytest.approx(-planted[0], abs=1.)
    assert shift['dy'] == pytest.approx(-planted[1], abs=1.)
    assert shift['dz'] == pytest.approx(-planted[2], abs=0.1)

def test_sad_planted_shift(tmp_path):
    #Limit search window, sad evaluates every integer offset
    shift, _, _ = run_dem_align(tmp_path, 'sad', extra_args=['-max_offset', '30'])
    assert shift['dx'] == pytest.approx(-planted[0], abs=1.)
    assert shift['dy'] == pytest.approx(-planted[1], abs=1.)
    assert shift['dz'] == pytest.approx(-planted[2], abs=0.1)

def test_max_offset_is_horizontal(tmp_path):
    #Large vertical offset with small horizontal offset should not exceed max_offset (issue #65)
    ref_fn = str(tmp_path / 'ref.tif')
    src_fn = str(tmp_path / 'src.tif')
    make_dem(ref_fn)
    make_dem(src_fn, dx=5., dy=5., dz=150.)
    outdir = str(tmp_path / 'out')
    dem_align.main(['-mode', 'nuth', '-mask_list', 'none', '-max_dz', '300', '-outdir', outdir, ref_fn, src_fn])
    stats = json.load(open(glob.glob(os.path.join(outdir, '*_align_stats.json'))[0]))
    assert stats['shift']['dz'] == pytest.approx(-150., abs=0.05)

def test_unprojected_input_exits(tmp_path):
    ref_fn = str(tmp_path / 'ref.tif')
    geo_fn = str(tmp_path / 'geo.tif')
    make_dem(ref_fn)
    gdal.Warp(geo_fn, ref_fn, dstSRS='EPSG:4326')
    with pytest.raises(SystemExit) as e:
        dem_align.main(['-mask_list', 'none', '-outdir', str(tmp_path / 'out'), ref_fn, geo_fn])
    assert 'projected CRS' in str(e.value.code)
