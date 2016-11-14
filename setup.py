#!/usr/bin/env python

from distutils.core import setup

#To prepare a new release
#python setup.py sdist upload

setup(name='demcoreg',
    version='0.1.2',
    description='Utilities for DEM co-registration',
    author='David Shean',
    author_email='dshean@gmail.com',
    license='MIT',
    url='https://github.com/dshean/demcoreg',
    packages=['demcoreg'],
    long_description=open('README.md').read(),
    install_requires=['numpy','gdal','pygeotools','wget'],
    #Note: this will write to /usr/local/bin
    scripts=['demcoreg/pc_align_wrapper.sh', 'demcoreg/apply_dem_translation.py', 'demcoreg/compute_dz.py', 'demcoreg/dem_align.py', 'demcoreg/robust_stats.py', 'demcoreg/dem_mask.py']
)

