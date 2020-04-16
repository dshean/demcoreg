# Python/Bash Beginners guide to demcoreg
**Operating System**
- The intructions below are for a *.nix machine (Linux/Mac). There are some reported glitches for Windows.
- If a user wants to install on Windows, please install Ubuntu subsystem for Windows. We will post another set of instructions to help with that.

**Dependent Software**
- This pacakage is dependent on several python packages for full functionality, all of which need to be installed to execute it sucessfully. For new users, this can be daunting, the instructions here are with such users in mind. Powerusers are free to install as per their liking.

**Installations**
- Managing basic python packages is best done using conda, so for new users, we recomended that they download miniconda for their machine (depends on OS and 32 or 64 bit machine)  from here [miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Open a terminal
- Once miniconda is installed, you can create an environment to manage demcoreg related softwares.
- A bare bone environment can be created as `conda create -c conda-forge -n demcoreg python=3.7 gdal=2.4 rasterio geopandas` 
- Once the `demcoreg` environment has been setup, activate it using `conda activate demcoreg`. You will need to do this every time you open a terminal again.
- We will now download and install software from github (this needs to be done only once).
- Use mkdir command to make a directory, where software packages can be downloaded. Command `mkdir software`
- Use command `cd software` to change directory to software.
- Type the following commands (please maintain the order of the commands here).
- git clone https://github.com/dshean/demcoreg.git
- git colne https://github.com/dshean/pygeotools.git
- git clone https://github.com/dshean/imview.git
- pip insatall -e pygeotools/
- pip install -e imview/
- pip install -e demcoreg/
- You will now need to add the path to your .bashrc to use the executables from the command line.
- From the software directory, type `realpath demcoreg/demcoreg`
- Let the output from this command be for example `/softwares/demcoreg/demcoreg`.
- Open your ~/.bashrc in a text-editor of your choice, paste the path provided above as: `export PATH="/softwares/demcoreg/demcoreg:$PATH"`
- Repeat this step with `realpath pygeotools/pygeotools` and `realpath imview/imview`. Also add these paths as specified above in your .bashrc profile.
- Once all these three have been added, run `source ~/.bashrc`.
- Congratulations, You have demcoreg ready to use.

**Basic Usage**
- If the above sequence of commands are followed, you have all the command-line tools and python libraries in `demcoreg`, `pygeotools` and `imview` at your disposal (Look at their corresponding readmes for more details).
- `dem_align.py` is the workhorse here which can be used to align two DEMs. This is probably the most desirable operation for new users who install this package for aligning/co-registering two DEMs. 
- To align two DEMs, in the terminal type `dem_align.py -mode nuth dem1.tif dem2.tif` (*The filenames can be anything*).
*Note: When run for the first time, this will download glacier polygons from the rgi website which might take some time.*
- This will align the dem2.tif to dem1.tif, and will store the results in a folder (see the log output). The folders contain the aligned dem2.tif (filename ending with `*align.tif`), coregistered elevation difference (each pixel representing elevation difference at the pixel, filename ending with `*align_diff.tif`) and elevation difference over static surface (`*align_diff_filt.tif`). 
- You can open them in GIS software like QGIS/ArcGIS to analyse the above mentioned files. The output folder also contains a .png file ending with *align.png. This contains a summary of elevation difference distribution before and after co-registration (like example dem_align.py figure in readme). 

**Additional Usage**
- By default, when invoked, the program uses glacier polygon to define static surfaces. That can be turned off using the option `-mask_list {none}`, if the user needs to co-register DEMs at all pixels. In addition, the tool can use National LandCover Datasets for USA, and global bareground data to define an additional static surface mask. These can be invoked by `-mask_list {nlcd}` or `-mask_list {bareground}`. Two or more masking options can also be applied in combinations, by providing the desired masking operations as a comma seperated list. The command in this case becomes: `-mask_list {glaciers,nlcd}`.
- Solving planar tilts can be invoked using `-tiltcorr`. By default it the tilt is computed using a first order polynomial. That can be changed using the option `-polyorder 4`, if suppose a 4th order polynomail is desired to model the tilt.
- Users can also change the `slope_lim` and `max_dz` filters, which are by default fixed at 0.1 to 40 degrees and 100 m respectively. 
-The above mentioned knobs are the most frequently adjusted options. To get a sense of all parameter tweaks available, run the command as `dem_align.py -h`, which will provide the set of changable options with their default values. 

**Advanced: Using independent package functions in python code/notebook**
- Functions defined in the .py scripts such as `dem_mask.py`,`coreglib.py` or `dem_align.py` might be useful for some independent analysis by users. These can be imported into your notebook/python code. For example, if I wanted to import `get_icemask` function in `dem_mask.py` code from this repo, I would first import `dem_mask` as `from demcoreg import dem_mask`, and then invoke the `get_icemask` function when desired in the code as `dem_mask.get_icemask(gdal_dataset)` ([ice_mask_example](https://github.com/dshean/hma_mb_paper/blob/master/notebooks/nogzumpa_dh_dt_error_correlation.ipynb)). 
- Extensive documentations for all such functions are not available currently, but can be mined from the code block as you explore the package code, if found useful :D.
