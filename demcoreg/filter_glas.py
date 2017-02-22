#! /usr/bin/env python

#Can be run for a single DEM fn or a list of DEM fn
#parallel --jobs 16 --delay 1 '~/src/demcoreg/demcoreg/filter_glas.py {}' ::: */dem*/*DEM_32m.tif

import sys
import os

import numpy as np
from osgeo import gdal
from pygeotools.lib import geolib, iolib, malib

import matplotlib.pyplot as plt

from imview.lib import gmtColormap, pltlib
cpt_rainbow = gmtColormap.get_rainbow()
import glas_proc

min_pts = 100
glas_npz_fn = '/nobackupp8/deshean/icesat_glas/GLAH14_tllz_hma_lulcfilt_demfilt.npz'
#glas_csv_fn = '/nobackupp8/deshean/icesat_glas/GLAH14_tllz_hma_lulcfilt_demfilt.csv'
#glas_pts = np.loadtxt(glas_csv_fn, delimiter=',')
#np.savez_compressed(glas_npz_fn, glas_pts)

#This takes ~5 seconds to load ~9M records with 8 fields
print("Loading points: %s" % glas_npz_fn)
glas_pts = np.load(glas_npz_fn)['arr_0']
pt_srs = geolib.wgs_srs

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
    x = glas_pts[:,3]
    y = glas_pts[:,2]
    idx = ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)) 
    if idx.nonzero()[0].size == 0:
        print("No points after spatial filter")
        continue

    dem_mask_fn = os.path.splitext(dem_fn)[0]+'_ref.tif'
    print("Loading Masked DEM: %s" % dem_mask_fn)
    dem_mask_ds = gdal.Open(dem_mask_fn) 
    dem_mask = iolib.ds_getma(dem_mask_ds) 

    print("Sampling DEM at masked point locations") 
    glas_pts_fltr = glas_pts[idx]

    print("Writing out %i points after spatial filter" % glas_pts_fltr.shape[0]) 
    out_csv_fn = os.path.splitext(dem_fn)[0]+'_glas.csv'
    fmt='%0.8f, %0.10f, %0.6f, %0.6f, %0.2f, %0.2f, %0.2f, %i'
    np.savetxt(out_csv_fn, glas_pts_fltr, fmt=fmt, delimiter=',')

    x_fltr = glas_pts_fltr[:,3]
    y_fltr = glas_pts_fltr[:,2]
    z_fltr = glas_pts_fltr[:,4]
    mX_fltr, mY_fltr, mZ = geolib.cT_helper(x_fltr, y_fltr, 0, pt_srs, geolib.get_ds_srs(dem_mask_ds))
    pX_fltr, pY_fltr = geolib.mapToPixel(mX_fltr, mY_fltr, dem_mask_ds.GetGeoTransform())
    pX_fltr = np.atleast_1d(pX_fltr)
    pY_fltr = np.atleast_1d(pY_fltr)
    #This returns median and mad
    samp = glas_proc.sample(dem_mask_ds, mX_fltr, mY_fltr, pad=1)
    samp_idx = ~(np.ma.getmaskarray(samp[:,0]))
    if samp_idx.nonzero()[0].size == 0:
        print("No points after mask")
        continue
        
    glas_pts_fltr_mask = glas_pts_fltr[samp_idx]

    print("Writing out %i points after mask" % glas_pts_fltr_mask.shape[0]) 
    out_csv_fn = os.path.splitext(dem_fn)[0]+'_glas_ref.csv'
    #Could add DEM samp columns here
    #fmt='%0.8f, %0.10f, %0.6f, %0.6f, %0.2f, %0.2f, %0.2f, %i'
    np.savetxt(out_csv_fn, glas_pts_fltr_mask, fmt=fmt, delimiter=',')

    x_fltr_mask = glas_pts_fltr_mask[:,3]
    y_fltr_mask = glas_pts_fltr_mask[:,2]
    z_fltr_mask = glas_pts_fltr_mask[:,4]
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
        c = glas_pts_fltr[:,0]
        #vmin = c.min()
        #vmax = c.max()
        vmin = 2003.14085699
        vmax = 2009.77587047
        im3 = ax3.imshow(hs_ma, cmap='gray', clim=hs_clim, alpha=0.5)
        sc3 = ax3.scatter(pX_fltr, pY_fltr, s=1, c=c, vmin=vmin, vmax=vmax, edgecolors='none')
        cbar = pltlib.add_cbar(ax3, sc3, label='Pt Year', cbar_kwargs={'fmt':'%0.2f'})

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
        print "Saving figure"
        fig_fn = os.path.splitext(dem_fn)[0]+'_glas.png'
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight', pad_inches=0)
        plt.close(fig)
