#! /usr/bin/env python

import sys
import os
import argparse
import subprocess

import numpy as np
import matplotlib.pyplot as plt

from pygeotools.lib import warplib, iolib, malib, geolib
from imview.lib import pltlib
    
from pathlib import PurePath
import rasterio as rio
from rasterio.enums import Resampling

def load(fn, name, masked=True, squeeze=True, dtype=np.float32, chunks='auto'):
    """Function to open files with rasterio as xarray DataArrays
    
    Parameters
    -------------
    thisDir: str
        directory path to search
    fn_pattern: str
        regex pattern to match files
    verbose: boolean
        print filenames
    
    Returns
    -------------
    fns: list
        list of filenames matched and sorted
    """
    import rioxarray
    if squeeze:
        arr=rioxarray.open_rasterio(fn, masked=masked, default_name=name, chunks=chunks).squeeze(dim='band', drop=True)
    else:
        arr=rioxarray.open_rasterio(fn, masked=masked, default_name=name, chunks=chunks)
    
    if dtype: arr=arr.astype(dtype)
    
    return arr

def coerce_float(a, delimiter='+', verbose=False):
    """Docstring short description
    Parameters
    -------------
    thisDir: str
        directory path to search
    fn_pattern: str
        regex pattern to match files
    verbose: boolean
        print filenames
    
    Returns
    -------------
    fns: list
        list of filenames matched and sorted
    """
    if verbose: print(a)
    try:
        a=float(a)
    except:
        a = a.split('+')[-1]
        if verbose: print(a)
        a=float(a)
    if verbose: print(a)
    return a

def extract_shifts(align_fn, ref_fn, verbose=False):
    """Extract x, y, and z adjustments based on input align fn and coerce to float
    
    Parameters
    -------------
    align_fn: str
        aligned filename
    ref_fn: str
        reference filename
    verbose: boolean
        print statement checks
    
    Returns
    -------------
    x, y, z : tuple
        tuple of float shifts
    """
    x = os.path.basename(align_fn).split(f'{os.path.basename(ref_fn)[:-4]}')[-1].split('_x')[1].split('_')[0]
    y = os.path.basename(align_fn).split(f'{os.path.basename(ref_fn)[:-4]}')[-1].split('_y')[1].split('_')[0]
    z = os.path.basename(align_fn).split(f'{os.path.basename(ref_fn)[:-4]}')[-1].split('_z')[1].split('_')[0]
    if verbose: print(x, y, z)
    x = coerce_float(x)
    y = coerce_float(y)
    z = coerce_float(z)
    if verbose: print(x, y, z)
    if verbose: print(type(x), type(y), type(z))
    return x, y, z

def zshiftstr(zadjust):
    """Pre-pend sign for zshift >=0, add "+" and retain hundredths place

    Parameters
    -------------
    zadjust: float
        dz to stringify
    
    Returns
    -------------
    zstr: str
        string with prepended sign
    """
    if zadjust >= 0: zstr=f'+{zadjust:.2f}'
    else: zstr = f'{zadjust:.2f}'
    print(zstr)
    return zstr

def find_mode(arr, dec=2, verbose=False):
    '''By default, omits nans if present and computes mode over entire array.
    Can return multiple modes.
    If arr dtype is float, will round values to the specified decimals
    e.g., 1 = tenths (0.1), 2 = hundredths (0.01)
    
    Parameters
    -------------
    arr: np.array
        array from which to derive mode
    dec: int
        decimal places to round
    verbose: boolean
        print mode value and count
    
    Returns
    -------------
    modes: list
        list of identified modes
    '''
    # first remove nans
    arr = arr[~np.isnan(arr)]
    
    # then round the values
    arr = arr * 10**dec
    arr = arr.astype(int)
    
    # then compute mode(s) by using a dict and a counter loop approach
    # this could take awhile for large arrays
    # https://www.nobledesktop.com/learn/python/mode-in-python
    vals = {}
    for v in arr:
        if not v in vals:
            vals[v] = 1
        else:
            vals[v]+=1
    
    # this will return all the modes
    modes = [g for g,l in vals.items() if l==max(vals.values())]
    if verbose:
        _ = [print(f'{m / 10**dec}: {vals[m]}') for m in modes]
    modes = [m / 10**dec for m in modes]
    return modes

def getparser():
    parser = argparse.ArgumentParser(description="Perform two-stage co-registration using stable surface z shift", \
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('ref_fn', type=str, help='Reference DEM filename')
    parser.add_argument('stable_fn', type=str, help='Stable surface masked referenced DEM filename to be shifted')
    parser.add_argument('align_fn', type=str, help='Full surface aligned DEM filename to be stable surface median z shifted')
    parser.add_argument('-new_align_fn', default=None, type=str, help='Updated full surface two stage aligned DEM filename to output')
    parser.add_argument('-outdir', default=None, help='Output directory')
    parser.add_argument('-realrun', default=True, type=bool, help='Run calculations and write out files')
    parser.add_argument('-v', default=False, type=bool, help='Increase verbosity')
    return parser

def main(argv=None):
    parser = getparser()
    args = parser.parse_args()
    
    #Should check that files exist
    ref_fn = args.ref_fn
    stable_fn = args.stable_fn
    align_fn = args.align_fn
    updated_align_fn = args.new_align_fn
    realrun = args.realrun
    verbose = args.v
    outdir = args.outdir
    
    # consider adding these as built in args
    overwrite = False
    resampling = Resampling.average
    
    if outdir is None:
        outdir = PurePath(align_fn).parents[0]
    elif not os.path.exists(outdir):
        os.makedirs(outdir)

    outprefix = f'{PurePath(align_fn).stem.split("_nuth_")[0]}'
    
    print(f'\nCurrent directory: {os.getcwd()}')
    print(f'Running: dem_stable_zalign.py')
    
    print(f"\nReference: {ref_fn}")
    print(f"Stable reference: {stable_fn}")
    print(f"Full surface aligned DEM: {align_fn}") 
    
    # Extract shifts from align_fn based on input reference 
    # Note the horizontal shifts (xy) from full shift alignment png
    x, y, z = extract_shifts(align_fn, ref_fn=ref_fn, verbose=verbose)
    # print(x, y, z)
    
    # Load the aligned and the stable surface tif
    align_arr = load(align_fn, PurePath(align_fn).stem)
    stable_arr = load(stable_fn, PurePath(stable_fn).stem)

    # Explicitly write out no data values
    align_arr.rio.write_nodata(np.nan, inplace=True)
    stable_arr.rio.write_nodata(np.nan, inplace=True)

    if verbose:
        print('\nInput array resolutions:')
        print(f'align: {align_arr.rio.resolution()[0]}')
        print(f'stable ref: {stable_arr.rio.resolution()[0]}')
        
    # resample the stable surface reference to the aligned array
    # matching extents and the input aligned resolution
    stable_arr_reproj = stable_arr.rio.reproject_match(align_arr, 
                                                       resampling=resampling, 
                                                       nodata=np.nan)
    if verbose: print(f'stable ref reproject-matched: {stable_arr_reproj.rio.resolution()[0]}')

    # if verbose: print(f'Reproject matched check: {stable_arr_reproj.shape}, {align_arr.shape}')
    
    # create binary mask to remove smooshed non-stable surfaces
    import copy
    stable_binary = copy.deepcopy(stable_arr)

    stable_binaryval = stable_binary.values
    if verbose: print(f'Starting NaNcount: {np.isnan(stable_binaryval).sum()}')
    stable_binaryval[~np.isnan(stable_binaryval)] = 1
    stable_binaryval[np.isnan(stable_binaryval)] = 0
    stable_binary.values = stable_binaryval
    if verbose: print(f'Ending NaNcount after binarization: {np.isnan(stable_binary.values).sum()}')
    
    # Reproject the binary stable arr to match align_arr
    stable_binary_arr_reproj = stable_binary.rio.reproject_match(align_arr, 
                                                                 resampling=resampling, 
                                                                 nodata=np.nan)
    if verbose: print(f'Shapes check: {stable_binary_arr_reproj.shape}, {stable_arr_reproj.shape}, {align_arr.shape}')
    
    # Set pixel values that do not meet threshold to nan 
    strict_thresh = 1
    if verbose: print(f'Starting NaNcount before thresholding: {np.isnan(stable_arr_reproj.values).sum()}')
    stable_arr_reprojval = stable_arr_reproj.values
    stable_binary_arr_reprojval = stable_binary_arr_reproj.values
    stable_arr_reprojval[stable_binary_arr_reprojval < strict_thresh] = np.nan
    stable_arr_reproj.values = stable_arr_reprojval
    if verbose: print(f'Ending NaNcount before thresholding: {np.isnan(stable_arr_reproj.values).sum()}')
    
    # Undo the z shift from the warped full surface align.tif by adding the 
    # opposite zshift to both the horizontal and the full array
    if verbose: print('\nRe-setting z shifts')
    align_warp_xyalone = align_arr + -z
    align_arr = align_arr + -z
    
    # Compute the median difference over stable surface areas using the 
    # matched stable surface and the undone horizontally aligned array
    # this is done to match the align.tif (coarser) array
    print('\nComputing median difference over stable surfaces')
    stable_diff = align_warp_xyalone - stable_arr_reproj
    # if verbose: print(stable_diff.shape)
    
    # Calculate the median difference over stable surfaces to be used for z shift
    zshift = np.nanmedian(stable_diff.values)
    
    print(f'The median diff of (src - ref) is: {zshift:.2f} m') 
    
    calcmodes=True
    if calcmodes:
        # Include this bit to check the diff between med and mode
        modezshift = find_mode(stable_diff.values)
        if len(modezshift) > 1:
            print(f'The mode diffs of (src - ref) is:')
            _ = [print(f'{m:.2f} m') for m in modezshift]
        else:
            print(f'The mode diff of (src - ref) is: {modezshift[0]:.2f} m') 
        
        # sys.exit(1)
    # if this is still nan after all the rigamarole, then exit out on error
    if np.isnan(zshift):
        print('Z shift is nan, check for overlapping stable surfaces')
        print("\nReference: %s" % ref_fn)
        print("Stable reference: %s" % stable_fn)
        print("Full surface aligned DEM: %s" % align_fn)    
        sys.exit(1)

    # Note the z shift here for writing out file
    zstr = f'{-zshift:.2f}'
    print(f'z shift string for source adjustment is: {zshiftstr(float(zstr))}')

    # Apply z shift to the xy-aligned array that moves the median difference to 0
    print('\nApplying z shift')
    align_warp_fullshift = align_warp_xyalone - zshift

    # Apply the z shift to the full aligned array, not the warped horizontal array
    align_arr = align_arr - zshift

    # This diff is also src-ref now
    stable_diff_updated = align_warp_fullshift - stable_arr_reproj

    # check that the median diff is now 0 over stable surfaces
    print(f'New median diff ≈ 0? \
        {np.isclose(np.nanmedian(stable_diff_updated.values), 0, atol=1e-3)}')

    # Specify filename - this assumes nuth and kaab alignment, could leverage
    # the x bit 
    if updated_align_fn is None:
        updated_align_fn = f'{align_fn.split("_nuth_")[0]}_nuth_twostageshift_x{zshiftstr(x)}_y{zshiftstr(y)}_z{zshiftstr(float(zstr))}_align.tif'
    if verbose: print(f'\nNew align.tif: {updated_align_fn}')
    
    # Update attributes and the name of array
    align_arr.name = PurePath(updated_align_fn).stem
    align_arr.attrs = {'stable_fn': stable_fn}

    # and write it out
    if not os.path.exists(updated_align_fn):
        print("\nDNE, writing out now...")
        print(updated_align_fn)
        if realrun: align_arr.rio.to_raster(updated_align_fn)
    else:
        print('Exists')
    
    # Specify diff fn    
    updated_diff_fn = f'{updated_align_fn[:-4]}_diff.tif'
    if not os.path.exists(updated_diff_fn):
        
        # Pull in the full ref_fn
        ref_arr = load(ref_fn, PurePath(ref_fn).stem)
        ref_arr.rio.write_nodata(np.nan, inplace=True)

        # reproject and match the original align arr
        ref_arr_reproj = ref_arr.rio.reproject_match(align_arr, 
                                                    resampling=resampling, 
                                                    nodata=np.nan)
        if verbose: print(ref_arr_reproj.rio.resolution())
        
        # Calculated the updated diff for all surfaces using the updated align_arr
        diff_updated = align_arr - ref_arr_reproj
        
        # and write it out
        print("\nDNE, writing out now...")
        print(updated_diff_fn)
        if realrun: diff_updated.rio.to_raster(updated_diff_fn)
        del diff_updated
    else:
        print('Exists')
        
    del x, y, z, align_arr, stable_arr, stable_arr_reproj, align_warp_xyalone, \
        stable_diff, zshift, zstr, align_warp_fullshift, stable_diff_updated
    
    print('\nMoving to stats and plotting section')
    align_stats_fn = f'{updated_align_fn[:-4]}_stats.json'
    # print(updated_align_fn)
    if os.path.exists(align_stats_fn) and not overwrite:
        print('Stats file exists, overwrite set to False')
        print(align_stats_fn)
    else:
        # consider adding json stats bits and plotting code here
        orig_diff_fn = f'{align_fn.split("_nuth")[0]}_orig_diff.tif'
        ref_base = PurePath(ref_fn).stem
        # dem_fn = f'{PurePath(align_fn).parents[1]}/{PurePath(align_fn.split("_" + ref_base)[0]).stem}.tif' # worked for stacks
        # dem_fn = f'{PurePath(align_fn).parents[0].stem}/{PurePath(align_fn.split("_" + ref_base)[0]).stem}.tif' # for the stacks_clean set up, Pléiades
        dem_fn = f'{align_fn.split("_dem_align")[0]}.tif'
        print(dem_fn, ref_fn, orig_diff_fn)
        src_fn_list = [dem_fn, ref_fn, orig_diff_fn] 
        
        elev_dslist = warplib.memwarp_multi_fn(src_fn_list=src_fn_list, res='mean')
        src_dem_clip_ds = elev_dslist[0]
        ref_dem_clip_ds = elev_dslist[1]
        ref_dem_orig = iolib.ds_getma(ref_dem_clip_ds)
        src_dem_orig = iolib.ds_getma(src_dem_clip_ds)
        
        # add bits for stable surface shifting plot
        _, stable_ds = warplib.memwarp_multi_fn([orig_diff_fn, stable_fn], res='first')
        stable_orig = iolib.ds_getma(stable_ds)
        
        #Needed for plotting
        ref_dem_hs = geolib.gdaldem_mem_ds(ref_dem_clip_ds, processing='hillshade', returnma=True, computeEdges=True)
        src_dem_hs = geolib.gdaldem_mem_ds(src_dem_clip_ds, processing='hillshade', returnma=True, computeEdges=True)
        diff_orig = src_dem_orig - ref_dem_orig 
        diff_orig_stats = malib.get_stats_dict(diff_orig, full=True)    
        # adapt for stable surface stat plotting only

        align_dslist = warplib.memwarp_multi_fn(src_fn_list=[updated_align_fn, ref_fn], res='mean', dst_ndv=-9999)
        src_dem_align = iolib.ds_getma(align_dslist[0])
        ref_dem_align = iolib.ds_getma(align_dslist[1])
        diff_align = src_dem_align - ref_dem_align
        diff_align_stats = malib.get_stats_dict(diff_align[~np.isnan(diff_align)], full=True)
        
        res = geolib.get_res(ref_dem_clip_ds)[0]
    
        
        dx, dy, dz = extract_shifts(updated_align_fn, ref_fn=ref_fn)
        dm_total = np.sqrt(dx**2 + dy**2 + dz**2)
    
        # write out stats jsons
        align_stats_fn = f'{updated_align_fn[:-4]}_stats.json'
        print(f'Pulling stats for {align_stats_fn}')
        align_stats = {}
        align_stats['src_fn'] = dem_fn 
        align_stats['ref_fn'] = ref_fn
        align_stats['stable_fn'] = stable_fn
        align_stats['updated_align_fn'] = updated_align_fn 
        align_stats['shift'] = {'dx':dx, 'dy':dy, 'dz':dz, 'dm':dm_total}
        align_stats['before'] = diff_orig_stats
        align_stats['after'] = diff_align_stats
    
        import json
        with open(align_stats_fn, 'w') as f:
            print('Dumping stats json')
            if realrun: json.dump(align_stats, f)
    
    if realrun:
        fig_fn = f'{updated_align_fn[:-4]}.png'
        
        if os.path.exists(fig_fn) and not overwrite:
            print('Updated figure png exists')
            print(fig_fn)
        else:
            
            f, axa = plt.subplots(2, 4, figsize=(16, 8))
            for ax in axa.ravel()[:-1]:
                ax.set_facecolor('k')
                pltlib.hide_ticks(ax)
            dem_clim = malib.calcperc(ref_dem_orig, (2,98))
            # cmap = plt.cm.jet
            cmap = 'cpt_rainbow' # plt.cm.rainbow

            kwargs = {'interpolation':'none'}
            axa[0,0].imshow(ref_dem_hs, cmap='gray', **kwargs)
            im = axa[0,0].imshow(ref_dem_orig, cmap=cmap, clim=dem_clim, alpha=0.6, **kwargs)
            pltlib.add_cbar(axa[0,0], im, arr=ref_dem_orig, clim=dem_clim, label=None)
            pltlib.add_scalebar(axa[0,0], res=res)
            axa[0,0].set_title('Reference DEM')

            axa[0,1].imshow(src_dem_hs, cmap='gray', **kwargs)
            im = axa[0,1].imshow(src_dem_orig, cmap=cmap, clim=dem_clim, alpha=0.6, **kwargs)
            pltlib.add_cbar(axa[0,1], im, arr=src_dem_orig, clim=dem_clim, label=None)
            axa[0,1].set_title('Source DEM')

            # just use stable surface here
            arr = np.ones_like(stable_orig.data, dtype='uint8')
            arr[stable_orig.mask == True] = 0
            # arr = np.ones_like(ref_dem_orig.data, dtype='uint8')
            # arr[ref_dem_orig.mask == True] = 0
            # arr[src_dem_orig.mask == True] = 0
            axa[0,2].imshow(arr, clim=(0,1), cmap='gray', **kwargs)
            axa[0,2].set_title('Surfaces for median shift')

            #This is empty
            axa[0,3].axis('off')

            dz_clim = malib.calcperc_sym(diff_orig, (5, 95))
            im = axa[1,0].imshow(diff_orig, cmap='RdBu', clim=dz_clim, **kwargs)
            pltlib.add_cbar(axa[1,0], im, arr=diff_orig, clim=dz_clim, label=None)
            axa[1,0].set_title('Elev. Diff. Before (m)')

            im = axa[1,1].imshow(diff_align, cmap='RdBu', clim=dz_clim, **kwargs)
            pltlib.add_cbar(axa[1,1], im, arr=diff_align, clim=dz_clim, label=None)
            axa[1,1].set_title('Elev. Diff. After (m)')

            tight_dz_clim = (-2.0, 2.0)
            im = axa[1,2].imshow(diff_align, cmap='RdBu', clim=tight_dz_clim, **kwargs)
            pltlib.add_cbar(axa[1,2], im, arr=diff_align, clim=tight_dz_clim, label=None)
            axa[1,2].set_title('Elev. Diff. After (m)')

            arr = np.ones_like(ref_dem_orig.data, dtype='uint8')
            arr[ref_dem_orig.data == -9999] = 0
            updated_arr = np.ones_like(ref_dem_align.data, dtype='uint8')
            updated_arr[ref_dem_align.data == -9999] = 0

            bins = np.linspace(dz_clim[0], dz_clim[1], 128)
            diff_orig_compressed = diff_orig[arr==1]
            diff_align_compressed = diff_align[updated_arr==1]
            axa[1,3].hist(diff_orig_compressed, bins, color='g', label='Before', alpha=0.5)
            axa[1,3].hist(diff_align_compressed, bins, color='b', label='After', alpha=0.5)
            axa[1,3].set_xlim(*dz_clim)
            axa[1,3].axvline(0, color='k', linewidth=0.5, linestyle=':')
            axa[1,3].set_xlabel('Elev. Diff. (m)')
            axa[1,3].set_ylabel('Count (px)')
            axa[1,3].set_title("Source - Reference")

            before_str = 'Before\nmed: %0.2f\nnmad: %0.2f' % (diff_orig_stats['med'], diff_orig_stats['nmad'])
            axa[1,3].text(0.05, 0.95, before_str, va='top', color='g', transform=axa[1,3].transAxes, fontsize=8)
            after_str = 'After\nmed: %0.2f\nnmad: %0.2f' % (diff_align_stats['med'], diff_align_stats['nmad'])
            axa[1,3].text(0.65, 0.95, after_str, va='top', color='b', transform=axa[1,3].transAxes, fontsize=8);

            suptitle = '%s\nx: %+0.2fm, y: %+0.2fm, z: %+0.2fm' % (outprefix, dx, dy, dz)
            f.suptitle(suptitle)
            f.tight_layout()
            plt.subplots_adjust(top=0.90)
            
            print(f"Writing out figure: {fig_fn}")
            f.savefig(fig_fn, dpi=300)

if __name__ == "__main__":
    main()
