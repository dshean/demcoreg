# Python/Bash Beginner's guide to demcoreg

The demcoreg README.md provides a basic overview of installation and usage for users with basic Python, bash, and git/github proficiency.  This document is intended to provide additional support for new users of these tools, who often need to co-register two DEMs for a relatively limited application, without all of the advanced options.  

## Operating System
- The intructions below are for a \*nix machine (Linux or Mac OS X).
- For Windows, we reccomend that you install [Ubuntu subsystem for Windows](https://docs.microsoft.com/en-us/windows/wsl/install-win10) and then follow the instructions below. Since you're using a linux subsystem, you want to use miniconda for linux in the first installation step (not Windows). 

## Installation
1. Install conda. Managing python packages is best done using conda. You can download miniconda for your OS/architecture here: https://docs.conda.io/en/latest/miniconda.html
2. Open a terminal/shell
3. Create and activate a Python environment called demcoreg_env with necessary packages installed:
  - Run this command to create an environment: `conda create -c conda-forge -n demcoreg_env python=3.7 gdal=2.4 rasterio geopandas` 
  - After this completes, the new environment must be "activated". You can do this with: `conda activate demcoreg_env`
    - *Note: you will need to activate each time you start a new terminal session*
4. Prepare to install software from github (this needs to be done only once)
  - Create a directory to store code from github repositories. One option is `~/src` which creates a new subdirectory `src` in your home directory, usually `/Users/loginname` (shorthand `~`). In the terminal, run: `mkdir ~/src`
  - Navigate to the new subdirectory: `cd ~/src`
5. Clone github repositories:
  - `git clone https://github.com/dshean/pygeotools.git`
  - `git clone https://github.com/dshean/demcoreg.git`
  - `git clone https://github.com/dshean/imview.git`
6. Install these packages, so you can use them with your conda Python:
  - `pip install -e pygeotools/`
  - `pip install -e demcoreg/`
  - `pip install -e imview/`
7. In addition to Python modules, these packages also contain some command-line scripts.  While you can always run these scripts from the terminal using a full path (e.g., `~/src/pygeotools/pygeotools/warptool.py`), it's convenient to run them using only `warptool.py`. To accomplish this, you can add the directory to the `~/.bashrc` (or `~/.bash_profile`) file in your home directory.
  - To get the full path name `realpath demcoreg/demcoreg`
  - Open `~/.bashrc` (or `~/.bash_profile`) in a text editor of your choice, and add this line to the end of the file: `export PATH="~/src/pygeotools/pygeotools:~/src/demcoreg/demcoreg:~/src/imview/imview:$PATH"`
    - Navigate to your home directory using `cd ~`
    - Type `ls -al` in the command line to view a list of files (including hidden files) in the home directory 
    - Determine if you have `~/.bashrc` (or `~/.bash_profile`) 
      - If neither exists, you can create with `touch ~/.bashrc`
    - Edit the corresponding file with a text editor
      - On OS X, can run `open ~/.bashrc` (or `open ~/.bash_profile`) to open the file with TextEdit.app
      - On most platforms, can edit directly using `nano`, `vim`, `emacs` or other text editor
  - Run `source ~/.bashrc` in your current terminal session (won't need to do this in the future)
8. You may need to reactivate the demcoreg envionment if you ran `source ~/.bashrc`. As in step 3, do so with: `conda activate demcoreg_env`
9. In your terminal, run `dem_align.py -h`.  You should see the usage statement starting with:
```
usage: dem_align.py [-h] [-mode {ncc,sad,nuth,none}]
                    [-mask_list {toa,snodas,modscag,bareground,glaciers,nlcd,none} [{toa,snodas,modscag,bareground,glaciers,nlcd,none} ...]]
                    [-tiltcorr] [-polyorder POLYORDER] [-tol TOL]
                    [-max_offset MAX_OFFSET] [-max_dz MAX_DZ]
                    [-res {min,max,mean,common_scale_factor}]
                    [-slope_lim SLOPE_LIM SLOPE_LIM] [-max_iter MAX_ITER]
                    [-outdir OUTDIR]
                    ref_fn src_fn
```

## Basic command-line usage
- If the above sequence of commands are followed, you have all the command-line tools and python libraries in `demcoreg`, `pygeotools` and `imview` at your disposal (Look at their corresponding readmes for more details).
- `dem_align.py` is the workhorse here which can be used to align two rasters. This is probably the most desirable operation for new users who install this package for aligning/co-registering two DEMs. 
- To align two DEMs named dem1.tif and dem2.tif, run the following command: `dem_align.py -mode nuth dem1.tif dem2.tif`  
*Note: When run for the first time, demcoreg will download glacier polygons from the rgi website which might take some time.*
  - If you have no already downloaded the rgi glacier polygons, you will be prompt to run `get_rgi.sh`
  - In order to run `get_rgi.sh`, you will need to install `wget` if it is not already installed on your computer. You can do so using conda forge or brew: `conda install wget` or `brew install wget`
- This will co-register dem2.tif to dem1.tif, and will store the results in a subdirectory (see the log output). The folders contain the aligned dem2.tif (filename ending with `*_align.tif`), final elevation difference map (filename ending with `*_align_diff.tif`) and the elevation difference map over static surfaces used during alignment (`*_align_diff_filt.tif`). 
- You can open these GeoTiff files in GIS software like QGIS/ArcGIS to analyse the above mentioned files. The output folder also contains a .png file (ending with `*_align.png`). This contains figures of input DEMs, surfaces used for co-registration, elevation difference maps and histograms/stats before and after co-registration (see sample png figure in demcoreg README). 

## Advanced command-line usage
- Run the command as `dem_align.py -h` to get extended usage, explanations of options and default values. 
- By default, `dem_align.py` will exclude DEM pixels over glaciers during co-registration. To use all pixels in the input DEMs, use the option `-mask_list none`. 
- One can limit co-registration to DEM pixels over bare ground or snow-free areas.  Two options for landcover classification are available - the National LandCover Dataset (NLCD) for the United States, and a global bareground dataset. To limit co-registration over these surfaces, use options `-mask_list nlcd` or `-mask_list bareground`. Two or more masking options can be combined with a comma-seperated string: `-mask_list glaciers,nlcd`
- To remove residual elevation differences after co-registration, `dem_align.py` uses polynomial fits of arbitrary order `-tiltcorr`. By default this will use a first-order polynomial (planar fit), though the user can specified higher order fits (e.g., `-tiltcorr 3` for cubic fit)
- Users can also change the slope (`-slope_lim`) and maximum absolute elevation difference (`-max_dz`) filters.  By default, these filters exclude surface slopes outside of 0.1° to 40° and any pixels with absolute elevation difference greater than 100 m. 

## Python API: using demcoreg in a python script/notebook
- Most of the demcoreg code contains modular functions that are called by executables like `dem_align.py`. These can also be called directly from a Python script or interactive session (like Jupyter notebook)
- For example, to access the `get_icemask` function in `dem_mask.py`:
  - Add `from demcoreg import dem_mask` near the top of your Python script or notebook
  - Call the function `dem_mask.get_icemask(GDAL_DataSet)` for your input GDAL DataSet object (e.g., output from `gdal.Open(dem1.tif)`)
  - ([ice_mask_example](https://github.com/dshean/hma_mb_paper/blob/master/notebooks/nogzumpa_dh_dt_error_correlation.ipynb)). 

## Additional help
- The source code contains additional notes and documentation.  Users can file a Github Issue if the encounter errors or bugs.  We hope to improve documentation and simplify the installation process in the future, so could use help if you have the time!
