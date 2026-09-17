#! /usr/bin/env python
"""
Plotting functions for demcoreg figures

Subset of imview.lib.pltlib, copied here so demcoreg does not require imview
"""

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from pygeotools.lib import malib, geolib

def get_rainbow(rev=False):
    """
    Load GMT rainbow color palette (rainbow.cpt) as matplotlib colormap
    """
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rainbow.cpt')
    z = []
    rgb = []
    with open(fn) as f:
        for line in f:
            p = line.split()
            #Skip comments and background, foreground, nodata colors
            if not p or p[0].startswith('#') or p[0] in ('B', 'F', 'N'):
                continue
            #Each line is z1 r g b z2 r g b, keep the first edge of each segment
            z.append(float(p[0]))
            rgb.append([int(i)/255. for i in p[1:4]])
            last_z = float(p[4])
            last_rgb = [int(i)/255. for i in p[5:8]]
    #Add the second edge of the final segment
    z.append(last_z)
    rgb.append(last_rgb)
    z = np.array(z)
    pos = (z - z[0])/(z[-1] - z[0])
    cmap = LinearSegmentedColormap.from_list('cpt_rainbow', list(zip(pos, rgb)))
    if rev:
        cmap = cmap.reversed(name='cpt_rainbow_r')
    return cmap

#Note: not registered with matplotlib, as imview.lib.pltlib registers the same names
#Use cmap=pltlib.cpt_rainbow
cpt_rainbow = get_rainbow()
cpt_rainbow_r = get_rainbow(rev=True)

#Global imshow keyword arguments
#Note: matpltolib v2.0+ interpolates across masked values, disable interpolation for now
imshow_kwargs = {'interpolation':'none'}

#Global cbar keyword arguments
cbar_kwargs = {'orientation':'vertical'}

def iv(a, ax=None, clim=None, clim_perc=(2,98), cmap=cpt_rainbow, label=None, title=None, \
        ds=None, res=None, hillshade=False, scalebar=True):
    """
    Quick image viewer with standardized display settings
    """
    if ax is None:
        f,ax = plt.subplots()
    ax.set_aspect('equal')
    if clim is None:
        clim = get_clim(a, clim_perc)
    cm = cmap_setndv(cmap, cmap)
    alpha=1.0
    if hillshade:
        if ds is not None:
            hs = geolib.gdaldem_mem_ds(ds, processing='hillshade', computeEdges=True, returnma=True)
            b_cm = cmap_setndv('gray', cmap)
            #Set the overlay bad values to completely transparent, otherwise darkens the bg
            cm.set_bad(alpha=0)
            bg_clim_perc = (2,98)
            bg_clim = get_clim(hs, bg_clim_perc)
            bgplot = ax.imshow(hs, cmap=b_cm, clim=bg_clim, **imshow_kwargs)
            alpha = 0.5
    if scalebar:
        if ds is not None:
            #Get resolution at center of dataset
            ccoord = geolib.get_center(ds, t_srs=geolib.wgs_srs)
            #Compute resolution in local cartesian coordinates at center
            c_srs = geolib.localortho(*ccoord)
            res = geolib.get_res(ds, c_srs)[0]
        if res is not None:
            sb_loc = best_scalebar_location(a)
            add_scalebar(ax, res, location=sb_loc)
    imgplot = ax.imshow(a, cmap=cm, clim=clim, alpha=alpha, **imshow_kwargs)
    cbar = add_cbar(ax, imgplot, label=label, arr=a, clim=clim)
    hide_ticks(ax)
    if title is not None:
        ax.set_title(title)
    plt.tight_layout()
    return ax

def get_clim(a, clim_perc=(2,98)):
    """
    Computer percentile stretch for input array
    """
    clim = malib.calcperc(a, clim_perc)
    if clim[0] == clim[1]:
        if clim[0] > a.fill_value:
            clim = (a.fill_value, clim[0])
        else:
            clim = (clim[0], a.fill_value)
    return clim

def get_cbar_extend(a, clim=None):
    """
    Determine whether we need to add triangles to ends of colorbar
    """
    if clim is None:
        clim = get_clim(a)
    extend = 'both'
    if a.min() >= clim[0] and a.max() <= clim[1]:
        extend = 'neither'
    elif a.min() >= clim[0] and a.max() > clim[1]:
        extend = 'max'
    elif a.min() < clim[0] and a.max() <= clim[1]:
        extend = 'min'
    return extend

def cmap_setndv(cmap1, cmap2=None):
    """
    Set default nodata mapping for colorbar
    """
    #Get copy of matplotlib colormap object, input can be name or object
    cmap = plt.get_cmap(cmap1).copy()
    if cmap2 is None:
        cmap2 = cmap1
    if 'inferno' in plt.get_cmap(cmap2).name:
        #Set nodata to opaque gray
        cmap.set_bad('0.5', alpha=1)
    else:
        #Set nodata to opaque black
        cmap.set_bad('k', alpha=1)
    return cmap

#Turn off ticks and tick labels
def hide_ticks(ax):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

def best_scalebar_location(a, length_pad=0.2, height_pad=0.1):
    """
    Attempt to determine best corner for scalebar based on number of unmasked pixels
    """
    a = malib.checkma(a)
    length = int(a.shape[1]*length_pad)
    height = int(a.shape[0]*height_pad)
    d = {}
    d['upper right'] = a[0:height,-length:].count()
    d['upper left'] = a[0:height,0:length].count()
    d['lower right'] = a[-height:,-length:].count()
    d['lower left'] = a[-height:,0:length].count()
    loc = min(d, key=d.get)
    return loc

def add_scalebar(ax, res, location='lower right', arr=None, kwargs=None):
    from matplotlib_scalebar.scalebar import ScaleBar
    if kwargs is None:
        kwargs = {}
    if arr is not None:
        location=best_scalebar_location(arr)
    sb = ScaleBar(res, location=location, border_pad=0.5, **kwargs)
    ax.add_artist(sb)

def add_cbar(ax, mappable, label=None, arr=None, clim=None, cbar_kwargs=cbar_kwargs, fontsize=10, format=None):
    """
    Add colorbar to axes for previously plotted mappable (output from imshow)
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    fig = ax.get_figure()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad="2%")
    #Don't modify the global defaults
    cbar_kwargs = dict(cbar_kwargs)
    if arr is not None and clim is not None:
        cbar_kwargs['extend'] = get_cbar_extend(arr, clim=clim)
    if format is not None:
        cbar_kwargs['format'] = format
    cbar = fig.colorbar(mappable, cax=cax, **cbar_kwargs)
    if label is not None:
        cbar.set_label(label, size=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    #Set colorbar to be opaque, even if image is transparent
    cbar.set_alpha(1)
    #Need this to update after alpha change
    cbar._draw_all()
    return cbar
