#! /usr/bin/env python

#Create plot of dem_align results for many input files

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
from pygeotools.lib import geolib, malib

def make_plot3d(x, y, z, title=None, orthogonal_fig=True):
    cmean = np.mean([x,y,z], axis=1)
    cmed = np.median([x,y,z], axis=1)

    ce90 = geolib.CE90(x,y)
    le90 = geolib.LE90(z)
    coefs = [ce90, ce90, le90]
    maxdim = np.ceil(np.max([np.max(np.abs([x, y, z])), ce90, le90]))
    
    if orthogonal_fig:
        from matplotlib.patches import Ellipse
        fig_ortho = plt.figure(figsize=(10,4))
        #fig_ortho = plt.figure()
        title='Co-registration Translation Vector Components\nn=%i, mean: (%0.2f, %0.2f, %0.2f)\nCE90: %0.2f, LE90: %0.2f' % (x.shape[0], cmean[0], cmean[1], cmean[2], ce90, le90)
        plt.suptitle(title) 

        m = '.'

        ax = fig_ortho.add_subplot(131)
        ax.plot(x, y, color='b', linestyle='None', marker=m, label='ICP correction vector')
        ax.plot(cmean[0], cmean[1], color='r', linestyle='None', marker='s', label='Mean')
        #ax.scatter(x, y)
        #ax.scatter(cmean[0], cmean[1], color='r', marker='s')
        ax.set_xlim(-maxdim, maxdim)
        ax.set_ylim(-maxdim, maxdim)
        ax.minorticks_on()
        ax.set_aspect('equal')
        ax.set_xlabel('X offset (m)')
        ax.set_ylabel('Y offset (m)')
        e = Ellipse((0,0), 2*ce90, 2*ce90, linewidth=0, alpha=0.1)
        ax.add_artist(e)
        plt.legend(prop={'size':8}, numpoints=1, loc='upper left')

        ax = fig_ortho.add_subplot(132)
        ax.plot(x, z, color='b', linestyle='None', marker=m, label='ICP correction vector')
        ax.plot(cmean[0], cmean[2], color='r', linestyle='None', marker='s', label='Mean')
        #ax.scatter(x, z)
        #ax.scatter(cmean[0], cmean[2], color='r', marker='s')
        ax.set_xlim(-maxdim, maxdim)
        ax.set_ylim(-maxdim, maxdim)
        ax.minorticks_on()
        ax.set_aspect('equal')
        ax.set_xlabel('X offset (m)')
        ax.set_ylabel('Z offset (m)')
        e = Ellipse((0,0), 2*ce90, 2*le90, linewidth=0, alpha=0.1)
        ax.add_artist(e)

        ax = fig_ortho.add_subplot(133)
        ax.plot(y, z, color='b', linestyle='None', marker=m, label='ICP correction vector')
        ax.plot(cmean[1], cmean[2], color='r', linestyle='None', marker='s', label='Mean')
        #ax.scatter(y, z)
        #ax.scatter(cmean[1], cmean[2], color='r', marker='s')
        ax.set_xlim(-maxdim, maxdim)
        ax.set_ylim(-maxdim, maxdim)
        ax.minorticks_on()
        ax.set_aspect('equal')
        ax.set_xlabel('Y offset (m)')
        ax.set_ylabel('Z offset (m)')
        e = Ellipse((0,0), 2*ce90, 2*le90, linewidth=0, alpha=0.1)
        ax.add_artist(e)
        
        plt.tight_layout()

        #Note: postscript doesn't properly handle tansparency
        fig_fn = 'dem_align_translation_vec_local_meters_orthogonal.pdf'
        plt.savefig(fig_fn, dpi=600, bbox_inches='tight')

def make_map(x, y, z, cx, cy):
    f, axa = plt.subplots(3, sharex=True, sharey=True, figsize=(5,10))
    vmin, vmax = (-15, 15)
    s=5
    cmap='RdYlBu'
    opt={'edgecolor':'k', 'vmin':vmin, 'vmax':vmax, 'cmap':cmap, 's':s, 'lw':0.3}
    sc = axa[0].scatter(cx, cy, c=x, **opt)
    axa[0].set_title("X-offset required to align")
    axa[0].set_aspect('equal')
    axa[1].scatter(cx, cy, c=y, **opt) 
    axa[1].set_title("Y-offset required to align")
    axa[2].scatter(cx, cy, c=z, **opt) 
    axa[2].set_title("Z-offset required to align")
    f.colorbar(sc, ax=axa.ravel().tolist())
    fig_fn = 'dem_align_map.png'
    f.savefig(fig_fn, dpi=300, bbox_inches='tight')

print("Building fn_list")
fn_list = glob.glob('*dem_align/*align.tif')
fn_list = sys.argv[1:]
print("Isolating x, y, z offsets")
xyz = np.array([np.array([a[1:] for a in np.array(fn.split('_'))[-4:-1]], dtype=float) for fn in fn_list])
print("Extracting center coords")
t_srs = geolib.hma_aea_srs
#t_srs = geolib.wgs_srs
ll = np.array([geolib.get_center(gdal.Open(fn), t_srs=t_srs) for fn in fn_list])
cy = ll[:,1]
cx = ll[:,0]
print(xyz.size)
m = np.sqrt(np.sum(np.square(xyz), axis=1))

#Throw out gross outliers
filter=True
if filter:
    stats = malib.print_stats(m)
    print(stats)
    f=3.5
    #thresh = stats[5]+stats[6]*f
    thresh = 90
    idx = (m > thresh)
    bad_fn = np.array(fn_list)[idx]
    np.savetxt('bad_fn.txt', bad_fn, fmt='%s')
    xyz = xyz[~idx]
    cx = cx[~idx]
    cy = cy[~idx]
    print(xyz.size)

x = xyz[:,0]
y = xyz[:,1]
z = xyz[:,2]

print("Creating plot")
make_plot3d(x, y, z)
print("Creating map")
make_map(x,y,z,cx,cy)
