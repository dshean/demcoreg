#! /usr/bin/env python

#Create plot of dem_align results for many input files

import glob
import numpy as np
import matplotlib.pyplot as plt
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
        #fig_fn = 'icp_translation_vec_proj_meters_orthogonal.pdf'
        fig_fn = 'icp_translation_vec_local_meters_orthogonal.pdf'
        plt.savefig(fig_fn, dpi=600, bbox_inches='tight')

print("Building fn_list")
fn_list = glob.glob('*dem_align/*align.tif')
print("Isolating x, y, z offsets")
xyz = np.array([np.array([a[1:] for a in np.array(fn.split('_'))[-4:-1]], dtype=float) for fn in fn_list])
print(xyz.size)
m = np.sqrt(np.sum(np.square(xyz), axis=1))
stats = malib.print_stats(m)
f=3.5
thresh = stats[5]+stats[6]*f
thresh = 90
idx = (m > thresh)
bad_fn = np.array(fn_list)[idx]
np.savetxt('bad_fn.txt', bad_fn, fmt='%s')
xyz = xyz[~idx]
print(xyz.size)

x = xyz[:,0]
y = xyz[:,1]
z = xyz[:,2]

print("Creating plot")
make_plot3d(x, y, z)

