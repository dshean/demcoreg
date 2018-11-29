#! /usr/bin/env python

#Create plot of dem_align results for many input files

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
import pandas as pd
from pygeotools.lib import geolib, malib, iolib

#Run as:
#dem_align_post.py $(ll *tif | grep 'alongtrack/WV03' | awk '{print $9}' | sed 's#.tif#_dem_align/*[0-9]_align.tif#')
#dem_align_post.py $(ll *tif | grep 'alongtrack/WV03' | awk '{print $9}' | sed 's#.tif#_dem_align/*tiltcorr_align.tif#')
#dem_align_post.py $(ll *tif | grep 'alongtrack/WV0[12]' | awk '{print $9}' | sed 's#.tif#_dem_align/*tiltcorr_align.tif#')
#dem_align_post.py $(ll *tif | grep 'QB02' | awk '{print $9}' | sed 's#.tif#_dem_align/*align.tif#')
#cat wv3_at_list.txt | sed 's#.tif#_dem_align/*align.tif#' > wv3_at_list_align.txt

outdir = 'dem_align_aster'
#outdir = 'dem_align_at_wv3'
#outdir = 'dem_align_at_wv12'
#outdir = 'dem_align_qb'
#outdir = 'dem_align_noqb'
#outdir = 'dem_align_at'
#outdir = 'dem_align_ct'

if not os.path.exists(outdir):
    os.makedirs(outdir)
out_fn_prefix = os.path.join(outdir, outdir)

#Throw out gross outliers
filter=True
mv_bad=True
#WV/GE
#outlier_mag_thresh = 20 
#ASTER
outlier_mag_thresh = 90

def make_plot3d(x, y, z, title=None, orthogonal_fig=True):
    cmean = np.mean([x,y,z], axis=1)
    cstd = np.std([x,y,z], axis=1)
    cmed = np.median([x,y,z], axis=1)
    cnmad = malib.mad([x,y,z], axis=1)
    x_corr = x - cmean[0]
    y_corr = y - cmean[1]
    z_corr = z - cmean[2]
   
    ce90 = geolib.CE90(x,y)
    ce90_corr = geolib.CE90(x_corr,y_corr)
    le90 = geolib.LE90(z)
    le90_corr = geolib.LE90(z_corr)

    coefs = [ce90, ce90, le90]
    #maxdim = np.ceil(np.max([np.max(np.abs([x, y, z])), ce90, le90]))
    maxdim = np.ceil(np.max([np.percentile(np.abs([x, y, z]), 99), ce90, le90]))
    
    if orthogonal_fig:
        from matplotlib.patches import Ellipse
        #fig_ortho, axa = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(10,5))
        fig_ortho, axa = plt.subplots(1, 3, figsize=(10,5))
        title = 'Co-registration Translation Vector Components, n=%i\n' % x.shape[0]
        title += 'mean: (%0.2f, %0.2f, %0.2f), std: (%0.2f, %0.2f, %0.2f)\n' % (tuple(cmean) + tuple(cstd))
        title += 'med: (%0.2f, %0.2f, %0.2f), nmad: (%0.2f, %0.2f, %0.2f)\n' % (tuple(cmed) + tuple(cnmad))
        title += 'CE90: %0.2f (Bias-corrected: %0.2f), LE90: %0.2f (Bias-corrected: %0.2f)' % (ce90, ce90_corr, le90, le90_corr)
        plt.suptitle(title) 

        dot_prop={'color':'k', 'linestyle':'None', 'marker':'.', 'ms':3, 'label':'ICP correction vector', 'alpha':0.5}
        mean_prop={'color':'r', 'linestyle':'None', 'marker':'o', 'label':'Mean'}

        for ax in axa:
            ax.set_xlim(-maxdim, maxdim)
            ax.set_ylim(-maxdim, maxdim)
            ax.minorticks_on()
            ax.set_aspect('equal')

        axa[0].plot(x, y, **dot_prop)
        axa[0].plot(cmean[0], cmean[1], **mean_prop)
        axa[0].set_xlabel('X offset (m)')
        axa[0].set_ylabel('Y offset (m)')
        e = Ellipse((0,0), 2*ce90, 2*ce90, linewidth=0, alpha=0.1)
        axa[0].add_artist(e)
        axa[0].legend(prop={'size':8}, numpoints=1, loc='upper left')

        axa[1].plot(x, z, **dot_prop)
        axa[1].plot(cmean[0], cmean[2], **mean_prop)
        axa[1].set_xlabel('X offset (m)')
        axa[1].set_ylabel('Z offset (m)')
        e = Ellipse((0,0), 2*ce90, 2*le90, linewidth=0, alpha=0.1)
        axa[1].add_artist(e)

        axa[2].plot(y, z, **dot_prop)
        axa[2].plot(cmean[1], cmean[2], **mean_prop)
        axa[2].set_xlabel('X offset (m)')
        axa[2].set_ylabel('Z offset (m)')
        e = Ellipse((0,0), 2*ce90, 2*le90, linewidth=0, alpha=0.1)
        axa[2].add_artist(e)
        
        plt.tight_layout()

        #Note: postscript doesn't properly handle tansparency
        fig_fn = '%s_translation_vec_local_meters_orthogonal.pdf' % out_fn_prefix
        plt.savefig(fig_fn, dpi=600, bbox_inches='tight')

def make_map(x, y, z, cx, cy):
    f, axa = plt.subplots(3, sharex=True, sharey=True, figsize=(5,10))
    axa[0].set_aspect('equal')
    maxdim = np.ceil(np.percentile(np.abs([x, y, z]), 99))
    #vmin, vmax = (-15, 15)
    vmin, vmax = (-maxdim, maxdim)
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
    fig_fn = '%s_map.png' % out_fn_prefix
    f.savefig(fig_fn, dpi=300, bbox_inches='tight')

print("Building fn_list")
#fn_list = glob.glob('*dem_align/*align.tif')
fn_list = sys.argv[1:]
fn_list = iolib.fn_list_valid(fn_list)
if not fn_list:
    sys.exit("No valid input files")
print("Isolating x, y, z offsets")
delim='_nuth_'
xyz = np.array([np.array([a[1:] for a in np.array(os.path.split(fn)[-1].split(delim)[-1].split('_'))[0:3]], dtype=float) for fn in fn_list])
print("Extracting center coords")
t_srs = geolib.hma_aea_srs
#t_srs = geolib.conus_aea_srs
#t_srs = geolib.wgs_srs
ll = np.array([geolib.get_center(gdal.Open(fn), t_srs=t_srs) for fn in fn_list])
cy = ll[:,1]
cx = ll[:,0]
m = np.sqrt(np.sum(np.square(xyz), axis=1))

df = pd.DataFrame(xyz, index=fn_list, columns=['x','y','z'])
df['m'] = m
df['cy'] = cy
df['cx'] = cx

df = df.sort_values(by='m', ascending=False)
print(df.shape[0])

if filter:
    print("Correction magnitude")
    stats = malib.print_stats(m)
    print("Outlier magnitude threshold:")
    print(outlier_mag_thresh)
    f=3.5
    print(stats[5]+stats[6]*f)
    idx = (df['m'] > outlier_mag_thresh)
    df[idx].to_csv('%s_bad_fn.txt' % out_fn_prefix, sep=' ', columns='m', index=True, header=False, float_format='%0.2f')

    if mv_bad:
        print("Moving bad solutions")
        import shutil
        baddir = '%s_bad' % out_fn_prefix
        os.makedirs(baddir)
        for i in df[idx].index:
            shutil.move(os.path.split(i)[0], baddir)

    df[~idx].to_csv('%s_good_fn.txt' % out_fn_prefix, sep=' ', columns='m', index=True, header=False)
    df = df[~idx]
    print(df.shape[0])

print("Creating plot")
make_plot3d(df['x'], df['y'], df['z'])
print("Creating map")
make_map(df['x'], df['y'], df['z'], df['cx'], df['cy'])
