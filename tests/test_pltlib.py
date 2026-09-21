"""Tests for demcoreg.pltlib"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from demcoreg import pltlib

def test_cpt_rainbow_not_registered():
    #imview registers these names unconditionally, so demcoreg must not
    assert 'cpt_rainbow' not in matplotlib.colormaps
    assert pltlib.cpt_rainbow.name == 'cpt_rainbow'
    assert pltlib.cpt_rainbow(0.0)[:3] == (144/255., 0., 111/255.)

def test_add_cbar_does_not_modify_defaults():
    a = np.ma.array(np.random.rand(20, 30))
    f, ax = plt.subplots()
    im = ax.imshow(a)
    pltlib.add_cbar(ax, im, arr=a, clim=(0.2, 0.8), format='%i')
    assert pltlib.cbar_kwargs == {'orientation':'vertical'}

def test_iv():
    a = np.ma.masked_less(np.random.rand(20, 30), 0.1)
    ax = pltlib.iv(a, res=10.)
    assert ax.images[0].get_cmap().name == 'cpt_rainbow'
    ax = pltlib.iv(a, cmap='RdBu', clim=(0, 1), scalebar=False)
    assert ax.images[0].get_clim() == (0, 1)
    ax = pltlib.iv(a, cmap='inferno')
    #Shared colormap object is not modified
    assert tuple(pltlib.cpt_rainbow.get_bad()) == (0., 0., 0., 0.)
