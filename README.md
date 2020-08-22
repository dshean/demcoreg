[![DOI](https://zenodo.org/badge/72886193.svg)](https://zenodo.org/badge/latestdoi/72886193)

# demcoreg
Python and shell scripts for co-registration of rasters, specifically horizontal and vertical alignment of digital elevation models (DEMs).

## Overview
All DEMs have some horizontal and vertical geolocation error.  It is important to remove relative offsets when differencing DEMs for elevation change analyses.  These tools offer several options to solve this problem.  Most solve for the sub-pixel horizontal shift and vertical offset required to minimize errors over "static" control surfaces.  The ASP `pc_align` tool can also solve for more complex transformations with rotations and scaling. 

## Features
- Multiple co-registration algorithms:
    - Automated, iterative implementation of the algorithm outlined by [Nuth and Kaab (2011)](https://www.the-cryosphere.net/5/271/2011/tc-5-271-2011.html)
    - Wrappers around the [Ames Stereo Pipeline](https://ti.arc.nasa.gov/tech/asr/intelligent-robotics/ngt/stereo/) Iterative Closest Point (ICP) `pc_align` utility (https://stereopipeline.readthedocs.io/en/latest/tools/pc_align.html)
    - Normalized cross-correlation (NCC) with sub-pixel refinment
    - Sum of absolute differences (SAD)
- Command-line utilities for raster differencing with necessary resampling (`compute_diff.py`), raster stats (`robust_stats.py`), and raster sampling at point locations (`sample_raster_at_pts.py`)
- Mask preparation and automatic determination of static control surfaces (i.e., exposed bedrock) for a user-specified combination of:
    - RGI glacier polygons
    - Land-use/Land-cover classification:
        - [National Land Cover Database (NLCD) 30-m products for CONUS](https://www.usgs.gov/centers/eros/science/national-land-cover-database?qt-science_center_objects=0#qt-science_center_objects) with pre-configured filter combinations (e.g., 'not_forest+not_water')
        - [Global Bare Ground 30-m products](https://glad.umd.edu/dataset/global-2010-bare-ground-30-m)
    - Snow mask: 
        - Thresholded MODSCAG fSCA from ~2-week period around DEM timestamp
        - Thresholded SNODAS model
        - Thresholded Top-of-atmosphere reflectance values from corresponding orthoimage (requires pregeneration)

### Some useful command-line utilities (run with `-h` option for complete usage)
- `dem_align.py` - robust raster DEM co-registration (e.g., Nuth and Kaab [2011]) for surfaces with variable slope and aspect (e.g., mountains)
- `dem_mask.py` - pre-generate mask to identify "stable" surfaces to use during co-registration
- `pc_align_wrapper.sh` - wrapper around NASA Ames Stereo Pipeline pc_align utility for iterative closest point co-registration 
- `apply_dem_translation.py` - update raster geotransform and remove vertical offset
- `compute_diff.py` - simple DEM difference calculation with intuitive resampling options
- `robust_stats.py` - print out robust raster statistics (e.g,. for DEM difference map before/after co-registration)

## Sample output 
### dem_align.py 
Sample command: `dem_align.py ref_dem.tif src_dem.tif`
![Sample dem_align](docs/20081123_0446_1735796131_1735796132_40m-DEM_hma_nasadem_hgt_lt5m_err_nuth_x+26.19_y+182.36_z-65.52_align_sm.jpg)
![Nuth and Kaab plot](docs/nuth_sample.jpg)

### dem_mask.py
Sample command: `dem_mask.py --toa --bareground --glaciers src_dem.tif`
![Sample dem_mask](docs/dem_mask_example_sm.jpg)

### filter_glas.py output
![Sample filter_glas](docs/20151227_0803_10200100499B7700_10200100496E3000-DEM_32m_glas_sm.jpg)

## Example applications
#### High-mountain Asia
- Co-registration of ~35000 high-resolution DEMs from multiple sensors
- https://github.com/dshean/hma_mb_paper

#### Pine Island Glacier, Antarctica
- Least-squares optimization to correct for offset and "tilt" of ~800 high-resolution DEMs with limited ground control
- https://github.com/dshean/pig_dem_meltrate

## Installation
We are hoping to clean up the code, remove unnecessary dependencies, and streamline installation using conda. For now, we recommend following the "Building from Latest Source" instructions below, to obtain latest features/bugfixes. 

If unfamiliar with this process, or if you are new to Python, bash, and/or git/github, start with these more detailed instructions and notes: [Beginner's guide for installation and basic usage](./docs/beginners_doc.md)

### Building from Latest Source (recommended)
1. Assuming you have working Python3 install with GDAL and NumPy, install [`pygeotools`](https://github.com/dshean/pygeotools)
1. Clone the `demcoreg` repository: `git clone https://github.com/dshean/demcoreg.git`
1. Perform developer install with pip: `pip install -e demcoreg`
    - *The -e flag ("editable mode", setuptools "develop mode") will allow you to modify source code and immediately see changes. Useful if you need to make minor tweaks or bugfixes (please submit a Pull Request!)*
1. Optionally, append the demcoreg subdirectory containing scripts to your PATH: `export PATH=${PATH}:$PWD/demcoreg/demcoreg` (replacing `$PWD` with the absolute path to the cloned demcoreg repository)
    - *To make this permanent, add that line to your shell config file (e.g., ~/.bashrc).* 

### Simple install with PyPI
`pip install demcoreg`
    
## Documentation
http://demcoreg.readthedocs.io (autogenerated from source code, may be out of date)

## License
This project is licensed under the terms of the MIT License.

### Citation
If you use any of this software for research applications that result in publications, please cite:  
Shean, D. E., O. Alexandrov, Z. Moratto, B. E. Smith, I. R. Joughin, C. C. Porter, Morin, P. J., An automated, open-source pipeline for mass production of digital elevation models (DEMs) from very high-resolution commercial stereo satellite imagery, ISPRS J. Photogramm. Remote Sens, 116, 101-117, doi: [10.1016/j.isprsjprs.2016.03.012](https://doi.org/10.1016/j.isprsjprs.2016.03.012), 2016. [<img src="http://wwwimages.adobe.com/content/dam/acom/en/legal/images/badges/Adobe_PDF_file_icon_24x24.png">](docs/Sheanetal_2016_ISPRS.pdf)
