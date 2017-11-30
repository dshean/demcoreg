#! /usr/bin/env python

#Filter preprocessed ICESat-1 GLAS points for a given input raster

#First run glas_proc.py - see usage
#Also expects *DEM_32m_ref.tif output from dem_mask.py

#Then run filter_glas.py for a single DEM fn or a list of DEM fn
#parallel --jobs 16 --delay 1 '~/src/demcoreg/demcoreg/filter_glas.py {}' ::: */dem*/*DEM_32m.tif

import sys
import os

import numpy as np
from osgeo import gdal
from pygeotools.lib import geolib, iolib, malib, timelib

import matplotlib.pyplot as plt

from imview.lib import gmtColormap, pltlib
cpt_rainbow = gmtColormap.get_rainbow()

site = 'hma'

#Minimum number of points required to write out _ref.csv
min_pts = 100

#Maximum value of surface slope to use
max_slope = 20.

pt_srs = geolib.wgs_srs
#This is time column in YYYYMMDD
tcol = 0
xcol = 3
ycol = 2
zcol = 4

#Padding in pixels for sample radius
#Since we're likely dealing with 32-m products here, can just use pad=1
pad = 1
#pad = 'glas'

glas_dir = '/nobackupp8/deshean/icesat_glas'
#ext = 'GLAH14_%s_refdemfilt_lulcfilt' % site
ext = 'GLAH14_%s_refdemfilt' % site
glas_npz_fn = os.path.join(glas_dir, ext+'.npz')

if not os.path.exists(glas_npz_fn):
    glas_csv_fn = os.path.splitext(glas_npz_fn)[0]+'.csv'
    print("Loading csv: %s" % glas_csv_fn)
    glas_pts = np.loadtxt(glas_csv_fn, delimiter=',', skiprows=1, dtype=None)
    print("Saving npz: %s" % glas_npz_fn)
    np.savez_compressed(glas_npz_fn, glas_pts)
else:
    #This takes ~5 seconds to load ~9M records with 8 fields
    print("Loading npz: %s" % glas_npz_fn)
    glas_pts = np.load(glas_npz_fn)['arr_0']

dem_fn_list = sys.argv[1:]
for n,dem_fn in enumerate(dem_fn_list):
    print("%i of %i" % (n+1, len(dem_fn_list)))
    #Lat/lon extent filter
    print("Loading DEM: %s" % dem_fn)
    dem_ds = gdal.Open(dem_fn)
    dem_ma = iolib.ds_getma(dem_ds)
    dem_extent_wgs84 = geolib.ds_extent(dem_ds, t_srs=pt_srs)
    xmin, ymin, xmax, ymax = dem_extent_wgs84
    print("Applying spatial filter") 
    x = glas_pts[:,xcol]
    y = glas_pts[:,ycol]
    idx = ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)) 
    if idx.nonzero()[0].size == 0:
        print("No points after spatial filter")
        continue

    print("Sampling DEM at masked point locations") 
    glas_pts_fltr = glas_pts[idx]

    print("Writing out %i points after spatial filter" % glas_pts_fltr.shape[0]) 
    out_csv_fn = os.path.splitext(dem_fn)[0]+'_%s.csv' % ext

    # dt_ordinal, dt_YYYYMMDD, lat, lon, z_WGS84 
    fmt = '%0.8f, %i, %0.6f, %0.6f, %0.2f'
    if glas_pts_fltr.shape[1] == 7:
        # dt_ordinal, dt_YYYYMMDD, lat, lon, z_WGS84, z_refdem_med_WGS84, z_refdem_nmad
        fmt += ', %0.2f, %0.2f'
    elif glas_pts_fltr.shape[1] == 8:
        # dt_ordinal, dt_YYYYMMDD, lat, lon, z_WGS84, z_refdem_med_WGS84, z_refdem_nmad, lulc
        fmt += ', %0.2f, %0.2f, %i'
    np.savetxt(out_csv_fn, glas_pts_fltr, fmt=fmt, delimiter=',')

    x_fltr = glas_pts_fltr[:,xcol]
    y_fltr = glas_pts_fltr[:,ycol]
    z_fltr = glas_pts_fltr[:,zcol]

    dem_mask_fn = os.path.splitext(dem_fn)[0]+'_ref.tif'
    if os.path.exists(dem_mask_fn):
        print("Loading Masked DEM: %s" % dem_mask_fn)
        dem_mask_ds = gdal.Open(dem_mask_fn) 
        dem_mask = iolib.ds_getma(dem_mask_ds) 
    else:
        dem_mask_ds = dem_ds
        dem_mask = dem_ma

    #Convert input xy coordinates to raster coordinates
    mX_fltr, mY_fltr, mZ = geolib.cT_helper(x_fltr, y_fltr, 0, pt_srs, geolib.get_ds_srs(dem_mask_ds))
    pX_fltr, pY_fltr = geolib.mapToPixel(mX_fltr, mY_fltr, dem_mask_ds.GetGeoTransform())
    pX_fltr = np.atleast_1d(pX_fltr)
    pY_fltr = np.atleast_1d(pY_fltr)

    #Sample raster
    #This returns median and mad for ICESat footprint
    samp = geolib.sample(dem_mask_ds, mX_fltr, mY_fltr, pad=pad)
    samp_idx = ~(np.ma.getmaskarray(samp[:,0]))
    npts = samp_idx.nonzero()[0].size
    if npts < min_pts:
        print("Not enough points after sampling valud pixels, post bareground mask (%i < %i)" % (npts, min_pts))
        continue
       
    if True:
        print("Applying slope filter, masking points with slope > %0.1f" % max_slope)
        slope_ds = geolib.gdaldem_mem_ds(dem_mask_ds, processing='slope', returnma=False)
        slope_samp = geolib.sample(slope_ds, mX_fltr, mY_fltr, pad=pad)
        slope_samp_idx = (slope_samp[:,0] <= max_slope).data
        samp_idx = np.logical_and(slope_samp_idx, samp_idx)

    npts = samp_idx.nonzero()[0].size
    if npts < min_pts:
        print("Not enough points after %0.1f deg slope mask (%i < %i)" % (max_slope, npts, min_pts))
        continue

    glas_pts_fltr_mask = glas_pts_fltr[samp_idx]

    if os.path.exists(dem_mask_fn):
        print("Writing out %i points after mask" % glas_pts_fltr_mask.shape[0]) 
        out_csv_fn_mask = os.path.splitext(out_csv_fn)[0]+'_ref.csv'
        #Could add DEM samp columns here
        np.savetxt(out_csv_fn_mask, glas_pts_fltr_mask, fmt=fmt, delimiter=',')

    x_fltr_mask = glas_pts_fltr_mask[:,xcol]
    y_fltr_mask = glas_pts_fltr_mask[:,ycol]
    z_fltr_mask = glas_pts_fltr_mask[:,zcol]
    mX_fltr_mask, mY_fltr_mask, mZ = geolib.cT_helper(x_fltr_mask, y_fltr_mask, 0, pt_srs, geolib.get_ds_srs(dem_mask_ds))
    pX_fltr_mask, pY_fltr_mask = geolib.mapToPixel(mX_fltr_mask, mY_fltr_mask, dem_mask_ds.GetGeoTransform())
    pX_fltr_mask = np.atleast_1d(pX_fltr_mask)
    pY_fltr_mask = np.atleast_1d(pY_fltr_mask)

    dz = z_fltr_mask - samp[samp_idx,0]

    if True:
        print "Creating plot of %i output points" % x_fltr.shape[0]
        fig_kw = {'figsize':(10,7.5)}
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharex=True, sharey=True, **fig_kw)

        #Plot DEM color shaded relief
        hs_ma = geolib.gdaldem_wrapper(dem_fn)
        hs_clim = malib.calcperc(hs_ma, perc=(0.5, 99.5))
        dem_clim = malib.calcperc(dem_ma)
        ax1.imshow(hs_ma, cmap='gray', clim=hs_clim)
        im1 = ax1.imshow(dem_ma, cmap=cpt_rainbow, clim=dem_clim, alpha=0.5)
        cbar = pltlib.add_cbar(ax1, im1, label='DEM Elev. (m WGS84)')
       
        #Plot all color points over shaded relief
        im2 = ax2.imshow(hs_ma, cmap='gray', clim=hs_clim, alpha=0.5)
        #Plot all points in black
        sc2 = ax2.scatter(pX_fltr, pY_fltr, s=0.5, c='k', edgecolors='none')
        #Plot valid in color
        c = z_fltr_mask 
        sc2 = ax2.scatter(pX_fltr_mask, pY_fltr_mask, s=0.5, c=c, cmap=cpt_rainbow, vmin=dem_clim[0], vmax=dem_clim[1], edgecolors='none')
        cbar = pltlib.add_cbar(ax2, sc2, label='Pt Elev. (m WGS84)')

        #Plot time
        c = glas_pts_fltr[:,tcol]
        c_decyear = timelib.np_dt2decyear(timelib.np_o2dt(c))
        c = c_decyear
        #vmin = c.min()
        #vmax = c.max()
        vmin = 2003.14085699
        vmax = 2009.77587047
        #vmin = 20030220
        #vmax = 20091011
        im3 = ax3.imshow(hs_ma, cmap='gray', clim=hs_clim, alpha=0.5)
        sc3 = ax3.scatter(pX_fltr, pY_fltr, s=1, c=c, vmin=vmin, vmax=vmax, edgecolors='none')
        #cbar = pltlib.add_cbar(ax3, sc3, label='Pt Year', cbar_kwargs={'format':'%0.2f'})
        cbar = pltlib.add_cbar(ax3, sc3, label='Pt Year')

        #Plot dz
        c = dz
        vmin, vmax = malib.calcperc(c, perc=(5, 95))
        absmax = np.max(np.abs([vmin, vmax]))
        vmin = -absmax
        vmax = absmax
        im4 = ax4.imshow(hs_ma, cmap='gray', clim=hs_clim, alpha=0.5)
        sc4 = ax4.scatter(pX_fltr_mask, pY_fltr_mask, s=2, c=c, cmap='RdYlBu', vmin=vmin, vmax=vmax, edgecolors='none')
        cbar = pltlib.add_cbar(ax4, sc4, label='GCP - DEM (m)')

        for ax in (ax1, ax2, ax3, ax4):
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            ax.set_aspect('equal', 'box-forced')

        title='%s \n %i valid points (%i initial)' % (os.path.splitext(os.path.split(dem_fn)[1])[0], pX_fltr_mask.shape[0], pX_fltr.shape[0])
        fig.suptitle(title)
        fig.tight_layout()
        #This adjusts subplots to fit suptitle
        plt.subplots_adjust(top=0.92)
        fig_fn = os.path.splitext(out_csv_fn)[0]+'.png'
        print "Saving figure: %s" % fig_fn
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight', pad_inches=0)
        plt.close(fig)
