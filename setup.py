#!/usr/bin/env python

from distutils.core import setup

#To prepare a new release
#python setup.py sdist upload

setup(name='demcoreg',
    version='0.3.0',
    description='Utilities for DEM co-registration',
    author='David Shean',
    author_email='dshean@gmail.com',
    license='MIT',
    url='https://github.com/dshean/demcoreg',
    packages=['demcoreg'],
    long_description=open('README.md').read(),
    install_requires=['numpy','gdal','pygeotools','wget'],
    #Note: this will write to /usr/local/bin
    scripts=['demcoreg/pc_align_wrapper.sh', 'demcoreg/apply_dem_translation.py', \
            'demcoreg/compute_dz.py', 'demcoreg/dem_align.py', \
            'demcoreg/dem_mask.py', 'demcoreg/dem_coreg.sh', 'demcoreg/dem_coreg_all.sh', \
            'demcoreg/vol_stats.py', 'demcoreg/robust_stats.py', \
            'demcoreg/glas_proc.py', 'demcoreg/filter_glas.py', \
            'demcoreg/get_nlcd.sh', 'demcoreg/get_bareground.sh', 'demcoreg/get_glacier_poly.sh']
)

