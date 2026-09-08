#!/usr/bin/env python

from setuptools import setup

#To prepare a new release, build from a clean archive and publish:
#git archive --format=tar --prefix=demcoreg/ vX.Y.Z | tar -x -C /tmp && (cd /tmp/demcoreg && uv build && uv publish dist/*)

setup(name='demcoreg',
    version='1.1.3',
    description='Utilities for DEM co-registration',
    author='David Shean',
    author_email='dshean@gmail.com',
    license='MIT',
    url='https://github.com/dshean/demcoreg',
    packages=['demcoreg'],
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    python_requires='>=3.8',
    install_requires=['numpy','gdal','pygeotools>=1.1.2','wget'],
    #Note: this will create local copy of executable scripts
    #scripts=['demcoreg/pc_align_wrapper.sh', 'demcoreg/apply_dem_translation.py', 'demcoreg/compute_diff.py', 'demcoreg/dem_align.py', 'demcoreg/dem_mask.py', 'demcoreg/dem_coreg.sh', 'demcoreg/dem_coreg_all.sh', 'demcoreg/vol_stats.py', 'demcoreg/robust_stats.py', 'demcoreg/glas_proc.py', 'demcoreg/filter_glas.py', 'demcoreg/get_nlcd.sh', 'demcoreg/get_bareground.sh', 'demcoreg/get_rgi.sh']
)

