# demcoreg

Python and shell scripts for co-registration of rasters, specifically digital elevation models (DEMs).

## Overview

## Features
- Wrapper for ASP pc_align function

### Some useful command-line utilities (run with no arguments for usage)
- apply_dem_translation.py - update geotransform and apply vertical offset
- ...

## Examples 

## Installation

Install the latest release from PyPI:

    pip install demcoreg 

**Note**: by default, this will deploy executable scripts in /usr/local/bin

### Building from source

Clone the repository and install:

    git clone https://github.com/dshean/demcoreg.git
    pip install -e demcoreg/

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
