# demcoreg
Python and shell scripts for co-registration of rasters, specifically digital elevation models (DEMs).

## Overview
All DEMs have some horizontal and vertical geolocation error.  It is important to remove relative offsets when differencing multiple DEMs for elevation change analyses.  These tools offer several options to solve this problem.  Most solve for the sub-pixel horizontal shift and vertical offset required to minimize errors over "static" control surfaces.  The ASP pc_align tool can also solve for more complex transformations with rotations and scaling.  

## Features
- Multiple co-registration algorithms (ICP, NCC, SAD, Nuth and Kaab [2011])
- Automatic determination of static control surfaces (i.e., exposed bedrock) for arbitrary DEM timestamp
- Least-squares optimization to correct a group of DEMs

### Some useful command-line utilities (run with no arguments for usage)
- dem_mask.py - generate mask of snow-free rock surfaces using reflectance, LULC, SNODAS, MODSCAG
- pc_align_wrapper.sh - wrapper around NASA Ames Stereo Pipeline pc_align utility for iterative closest point co-registration 
- apply_dem_translation.py - update geotransform and applies vertical offset
- compute_dz.py - simple DEM difference calculation
- robust_stats.py - print out robust statistics for sampled DEM differences before/after co-registration
- ...

- coreglib.py - implementation of various co-registration algorithms: Nuth and Kaab (2011), normalized cross-correlation with sub-pixel refinement, sum of absolute differences

## Examples 

## Documentation

http://demcoreg.readthedocs.io

## Installation

Install the latest release from PyPI:

    pip install demcoreg 

**Note**: by default, this will deploy executable scripts in /usr/local/bin

### Building from source

Clone the repository and install:

    git clone https://github.com/dshean/demcoreg.git
    pip install -e demcoreg

The -e flag ("editable mode", setuptools "develop mode") will allow you to modify source code and immediately see changes.

### Core requirements 
- [GDAL/OGR](http://www.gdal.org/)
- [NumPy](http://www.numpy.org/)
- [pygeotools](https://github.com/dshean/pygeotools)

### Optional requirements (needed for some functionality) 
- [matplotlib](http://matplotlib.org/)
- [SciPy](https://www.scipy.org/)
- [NASA Ames Stereo Pipeline (ASP)](https://ti.arc.nasa.gov/tech/asr/intelligent-robotics/ngt/stereo/)

## License

This project is licensed under the terms of the MIT License.
