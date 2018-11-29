#! /usr/bin/env python

"""
Library of functions that can be used for co-registration of raster data

For many situations, ASP pc_align ICP co-registration is superior to these approaches. See pc_align_wrapper.sh
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

from pygeotools.lib import malib, iolib

def apply_xy_shift(ds, dx, dy, createcopy=True):
    """
    Apply horizontal shift to GDAL dataset GeoTransform
    
    Returns:
    GDAL Dataset copy with updated GeoTransform
    """
    print("X shift: ", dx)
    print("Y shift: ", dy)
   
    #Update geotransform
    gt_orig = ds.GetGeoTransform()
    gt_shift = np.copy(gt_orig)
    gt_shift[0] += dx 
    gt_shift[3] += dy

    print("Original geotransform:", gt_orig)
    print("Updated geotransform:", gt_shift)

    #Update ds Geotransform
    if createcopy:
        ds_align = iolib.mem_drv.CreateCopy('', ds, 0)
    else:
        #Update in place, assume ds is opened as GA_Update
        ds_align = ds
    ds_align.SetGeoTransform(gt_shift)
    return ds_align

def apply_z_shift(ds, dz, createcopy=True):
    if isinstance(dz, np.ndarray):
        print("Z shift offset array mean: ", dz.mean())
    else:
        print("Z shift offset: ", dz) 
    if createcopy:
        ds_shift = iolib.mem_drv.CreateCopy('', ds, 0)
    else:
        ds_shift = ds
    b = ds_shift.GetRasterBand(1)
    a = iolib.b_getma(b)
    a += dz
    b.WriteArray(a.filled())
    return ds_shift

#Function for fitting Nuth and Kaab (2011)
def nuth_func(x, a, b, c):
    y = a * np.cos(np.deg2rad(b-x)) + c
    #Can use Phasor addition, but need to change conversion to offset dx and dy
    #https://stackoverflow.com/questions/12397412/i-know-scipy-curve-fit-can-do-better?rq=1
    #y = a * np.cos(np.deg2rad(x)) + b * np.sin(np.deg2rad(x)) + c
    return y

def compute_offset_sad(dem1, dem2, pad=(9,9), plot=False):
    """Compute subpixel horizontal offset between input rasters using sum of absolute differences (SAD) method
    """
    #This defines the search window size
    #Use half-pixel stride?
    #Note: stride is not properly implemented 
    #stride = 1
    #ref = dem1[::stride,::stride]
    #kernel = dem2[pad[0]:-pad[0]:stride, pad[1]:-pad[1]:stride]
    kernel = dem2[pad[0]:-pad[0], pad[1]:-pad[1]]
    #Want to pad evenly on both sides, so add +1 here
    m = np.zeros((pad[0]*2+1, pad[1]*2+1))
   
    #Find integer pixel offset
    i = j = 0
    for i in range(m.shape[0]):
        print(i)
        for j in range(m.shape[1]):
            print(j)
            ref = dem1[i:i+kernel.shape[0], j:j+kernel.shape[1]]
            diff = ref - kernel
            
            #Remove outliers beyond IQR
            diff_iqr = malib.calcperc(diff, (25,75))
            diff = np.ma.masked_outside(diff, *diff_iqr)
            """ 
            diff_med = np.ma.median(diff)
            diff_mad = malib.mad(diff)
            diff_madr = (diff_med - mad, diff_med + mad)
            diff = np.ma.masked_outside(diff, diff_madr)     
            """
            #Masked areas will decrease sum! Normalize by count of valid pixels
            m[i,j] = np.ma.abs(diff).sum()/diff.count()
    
    #Note, we're dealing with min SAD here, so want to provide -m for sub-pixel refinement 
    m = -m  

    int_argmax = np.array(np.unravel_index(m.argmax(), m.shape))
    int_offset = int_argmax - pad
    
    sp_argmax = np.array(find_subpixel_peak_position(m, 'parabolic'))
    sp_offset = sp_argmax - pad

    if plot:
        plt.figure()
        plt.title('Sum of Absolute Differences')
        plt.imshow(m)
        plt.scatter(*sp_argmax[::-1])
        #plt.show()

    return m, int_offset, sp_offset

#This is a decent full-image normalized cross-correlation routine with sub-pixel refinement
def compute_offset_ncc(dem1, dem2, pad=(9,9), prefilter=False, plot=False): 
    """Compute horizontal offset between input rasters using normalized cross-correlation (NCC) method
    """

    #Apply edge detection filter up front - improves results when input DEMs are same resolution
    if prefilter:
        print("Applying LoG edge-detection filter to DEMs")
        sigma = 1
        import scipy.ndimage
        #Note, ndimage alone propagates Nans and greatly reduces valid data area
        #Use the malib.nanfill wrapper to avoid this
        dem1 = malib.nanfill(dem1, scipy.ndimage.filters.gaussian_laplace, sigma) 
        dem2 = malib.nanfill(dem2, scipy.ndimage.filters.gaussian_laplace, sigma) 

    import scipy.signal
    #Compute max offset given dem spatial resolution
    #Should implement arbirary x and y search space
    #xsearch = (20, 41)
    #ysearch = (-10, 1)
    stride = 1
    ref = dem1[::stride,::stride]
    kernel = dem2[pad[0]:-pad[1]:stride, pad[0]:-pad[1]:stride]
    #kernel = dem2[-ysearch[0]:-ysearch[1]:stride, xsearch[0]:-xsearch[1]:stride]

    #Normalize
    ref = (ref - ref.mean()) / ref.std()
    kernel = (kernel - kernel.mean()) / kernel.std()

    #Consider using astropy.convolve here instead of scipy.correlate?

    print("Adding random noise to masked regions")
    #Generate random noise to fill gaps before correlation in frequency domain
    #Normal distribution N(mean, std^2)
    #ref_noise = ref.mask * ref.std() * np.random.rand(*ref.shape) + ref.mean()
    #kernel_noise = kernel.mask * kernel.std() * np.random.rand(*kernel.shape) + kernel.mean()
    #This provides noise in proper range, but noise propagates to m, peak is in different locations!
    #ref_noise = ref.mask * (ref.min() + ref.ptp() * np.random.rand(*ref.shape))
    #kernel_noise = kernel.mask * (kernel.min() + kernel.ptp() * np.random.rand(*kernel.shape))

    #This provides a proper normal distribution with mean=0 and std=1
    ref_noise = ref.mask * (np.random.randn(*ref.shape))
    kernel_noise = kernel.mask * (np.random.randn(*kernel.shape))
    #Add the noise
    ref = ref.filled(0) + ref_noise
    kernel = kernel.filled(0) + kernel_noise

    print("Running 2D correlation with search window (x,y): %i, %i" % (pad[1], pad[0]))
    m = scipy.signal.correlate2d(ref, kernel, 'valid')
    #This has memory issues, but ndimage filters can handle nan
    #m = scipy.ndimage.filters.correlate(ref, kernel)
   
    print("Computing sub-pixel peak")
    int_argmax = np.array(np.unravel_index(m.argmax(), m.shape))
    int_offset = int_argmax*stride - pad
    #int_offset = int_argmax*stride + np.array([ysearch[0], xsearch[0]]) 

    print(m.argmax())
    print(m.shape)
    print(int_argmax)
    print(int_offset)

    #Find sub-pixel peak
    sp_argmax = np.array(find_subpixel_peak_position(m, 'parabolic'))
    #May need to split this into integer and decimal components, multipy stride*int and add decimal
    #sp_offset = int_offset + (sp_argmax - int_argmax)
    sp_offset = sp_argmax - pad
    #sp_offset = sp_argmax + np.array([ysearch[0], xsearch[0]]) 

    print(sp_argmax)
    print(sp_offset)

    if plot: 
        fig, ax = plt.subplots()
        ax.set_title('NCC offset, parabolic SPR')
        ax.imshow(m)
        #plt.scatter(*int_argmax[::-1])
        ax.scatter(*sp_argmax[::-1])
    else:
        fig = None

    return m, int_offset, sp_offset, fig

#This is the Nuth and Kaab (2011) method
def compute_offset_nuth(dh, slope, aspect, min_count=100, remove_outliers=True, plot=True):
    """Compute horizontal offset between input rasters using Nuth and Kaab [2011] (nuth) method
    """
    import scipy.optimize as optimization

    if dh.count() < min_count:
        sys.exit("Not enough dh samples")
    if slope.count() < min_count:
        sys.exit("Not enough slope/aspect samples")

    #mean_dh = dh.mean()
    #mean_slope = slope.mean()
    #c_seed = (mean_dh/np.tan(np.deg2rad(mean_slope))) 
    med_dh = malib.fast_median(dh)
    med_slope = malib.fast_median(slope)
    c_seed = (med_dh/np.tan(np.deg2rad(med_slope))) 

    x0 = np.array([0.0, 0.0, c_seed])
  
    print("Computing common mask")
    common_mask = ~(malib.common_mask([dh, aspect, slope]))

    #Prepare x and y data
    xdata = aspect[common_mask].data
    ydata = (dh[common_mask]/np.tan(np.deg2rad(slope[common_mask]))).data

    print("Initial sample count:")
    print(ydata.size)

    if remove_outliers:
        print("Removing outliers")
        #print("Absolute dz filter: %0.2f" % max_dz)
        #diff = np.ma.masked_greater(diff, max_dz)
        #print(diff.count())

        #Outlier dz filter
        f = 3
        sigma, u = (ydata.std(), ydata.mean())
        #sigma, u = malib.mad(ydata, return_med=True)
        rmin = u - f*sigma
        rmax = u + f*sigma
        print("3-sigma filter: %0.2f - %0.2f" % (rmin, rmax))
        idx = (ydata >= rmin) & (ydata <= rmax)
        xdata = xdata[idx]
        ydata = ydata[idx]
        print(ydata.size)

    #Generate synthetic data to test curve_fit
    #xdata = np.arange(0,360,0.01)
    #ydata = f(xdata, 20.0, 130.0, -3.0) + 20*np.random.normal(size=len(xdata))
    
    #Limit sample size
    #n = 10000
    #idx = random.sample(range(xdata.size), n)
    #xdata = xdata[idx]
    #ydata = ydata[idx]

    #Compute robust statistics for 1-degree bins
    nbins = 360
    bin_range = (0., 360.)
    bin_width = 1.0
    bin_count, bin_edges, bin_centers = malib.bin_stats(xdata, ydata, stat='count', nbins=nbins, bin_range=bin_range)
    bin_med, bin_edges, bin_centers = malib.bin_stats(xdata, ydata, stat='median', nbins=nbins, bin_range=bin_range)
    #Needed to estimate sigma for weighted lsq
    #bin_mad, bin_edges, bin_centers = malib.bin_stats(xdata, ydata, stat=malib.mad, nbins=nbins, bin_range=bin_range)
    #Started implementing this for more generic binning, needs testing
    #bin_count, x_bin_edges, y_bin_edges = malib.get_2dhist(xdata, ydata, \
    #        xlim=bin_range, nbins=(nbins, nbins), stat='count')

    """
    #Mask bins in grid directions, can potentially contain biased stats
    #Especially true for SGM algorithm
    #badbins = [0, 90, 180, 270, 360]
    badbins = [0, 45, 90, 135, 180, 225, 270, 315, 360]
    bin_stat = np.ma.masked_where(np.around(bin_edges[:-1]) % 45 == 0, bin_stat)
    bin_edges = np.ma.masked_where(np.around(bin_edges[:-1]) % 45 == 0, bin_edges)
    """

    #Remove any bins with only a few points
    min_bin_sample_count = 9
    idx = (bin_count.filled(0) >= min_bin_sample_count) 
    bin_count = bin_count[idx].data
    bin_med = bin_med[idx].data
    #bin_mad = bin_mad[idx].data
    bin_centers = bin_centers[idx]

    fit = None
    fit_fig = None

    #Want a good distribution of bins, at least 1/4 to 1/2 of sinusoid, to ensure good fit
    #Need at least 3 valid bins to fit 3 parameters in nuth_func
    #min_bin_count = 3
    min_bin_count = 90 
    
    #Not going to help if we have a step function between two plateaus, but better than nothing
    #Calculate bin aspect spread
    bin_ptp = np.cos(np.radians(bin_centers)).ptp()
    min_bin_ptp = 1.0 

    #Should iterate here, if not enough bins, increase bin width
    if len(bin_med) >= min_bin_count and bin_ptp >= min_bin_ptp:

        print("Computing fit")
        #Unweighted fit
        fit = optimization.curve_fit(nuth_func, bin_centers, bin_med, x0)[0]

        #Weight by observed spread in each bin 
        #sigma = bin_mad
        #fit = optimization.curve_fit(nuth_func, bin_centers, bin_med, x0, sigma, absolute_sigma=True)[0]

        #Weight by bin count
        #sigma = bin_count.max()/bin_count
        #fit = optimization.curve_fit(nuth_func, bin_centers, bin_med, x0, sigma, absolute_sigma=False)[0]

        print(fit)

        if plot:
            print("Generating Nuth and Kaab plot")
            bin_idx = np.digitize(xdata, bin_edges)
            output = []
            for i in np.arange(1, len(bin_edges)):
                output.append(ydata[bin_idx==i])
            #flierprops={'marker':'.'}
            lw = 0.25
            whiskerprops={'linewidth':lw}
            capprops={'linewidth':lw}
            boxprops={'facecolor':'k', 'linewidth':0}
            medianprops={'marker':'o', 'ms':1, 'color':'r'}
            fit_fig, ax = plt.subplots(figsize=(6,6))
            #widths = (bin_width/2.0)
            widths = 2.5*(bin_count/bin_count.max())
            #widths = bin_count/np.percentile(bin_count, 50)
            #Stride
            s=3
            #This is inefficient, but we have list of arrays with different length, need to filter
            #Reduntant with earlier filter, should refactor
            bp = ax.boxplot(np.array(output)[idx][::s], positions=bin_centers[::s], widths=widths[::s], showfliers=False, \
                    patch_artist=True, boxprops=boxprops, whiskerprops=whiskerprops, capprops=capprops, \
                    medianprops=medianprops)
            bin_ticks = [0, 45, 90, 135, 180, 225, 270, 315, 360]
            ax.set_xticks(bin_ticks)
            ax.set_xticklabels(bin_ticks)
            """
            #Can pull out medians from boxplot
            #We are computing multiple times, inefficient
            bp_bin_med = []
            for medline in bp['medians']:
                bp_bin_med.append(medline.get_ydata()[0])
            """

            #Plot the fit
            f_a = nuth_func(bin_centers, fit[0], fit[1], fit[2])
            nuth_func_str = r'$y=%0.2f*cos(%0.2f-x)+%0.2f$' % tuple(fit)
            ax.plot(bin_centers, f_a, 'b', label=nuth_func_str)

            ax.set_xlabel('Aspect (deg)')
            ax.set_ylabel('dh/tan(slope) (m)')
            ax.axhline(color='gray', linewidth=0.5)

            ax.set_xlim(*bin_range)
            ylim = ax.get_ylim()
            abs_ylim = np.max(np.abs(ylim))
            #abs_ylim = np.max(np.abs([ydata.min(), ydata.max()]))
            #pad = 0.2 * abs_ylim 
            pad = 0
            ylim = (-abs_ylim - pad, abs_ylim + pad)
            minylim = (-10,10)
            if ylim[0] > minylim[0]:
                ylim = minylim
            ax.set_ylim(*ylim)
            ax.legend(prop={'size':8})

    return fit, fit_fig

#Attempt to fit polynomial functions to along-track and cross-track signals
#See demtools for existing code
def fit_at_ct():
    #Derive from image corners in projected array
    #Use known orbintal inclinations, project wgs geometry into srs
    img1_inc
    img2_inc
    #Rotate
    #Stats for rows, cols
    #Fit

#Function copied from from openPIV pyprocess
def find_first_peak(corr):
    """
    Find row and column indices of the first correlation peak.
    
    Parameters
    ----------
    corr : np.ndarray
        the correlation map
        
    Returns
    -------
    i : int
        the row index of the correlation peak
        
    j : int
        the column index of the correlation peak    
    
    corr_max1 : int
        the value of the correlation peak
    
    Original code from openPIV pyprocess

    """    
    ind = corr.argmax()
    s = corr.shape[1] 
    
    i = ind // s 
    j = ind %  s
    
    return i, j, corr.max()

#Function copied from from openPIV pyprocess
def find_subpixel_peak_position(corr, subpixel_method='gaussian'):
    """
    Find subpixel approximation of the correlation peak.
    
    This function returns a subpixels approximation of the correlation
    peak by using one of the several methods available. If requested, 
    the function also returns the signal to noise ratio level evaluated 
    from the correlation map.
    
    Parameters
    ----------
    corr : np.ndarray
        the correlation map.
        
    subpixel_method : string
         one of the following methods to estimate subpixel location of the peak: 
         'centroid' [replaces default if correlation map is negative], 
         'gaussian' [default if correlation map is positive], 
         'parabolic'.
         
    Returns
    -------
    subp_peak_position : two elements tuple
        the fractional row and column indices for the sub-pixel
        approximation of the correlation peak.

    Original code from openPIV pyprocess

    """
    # initialization
    default_peak_position = (corr.shape[0]/2,corr.shape[1]/2)

    # the peak locations
    peak1_i, peak1_j, dummy = find_first_peak(corr)
    
    try:
        # the peak and its neighbours: left, right, down, up
        c = corr[peak1_i, peak1_j]
        cl = corr[peak1_i-1, peak1_j]
        cr = corr[peak1_i+1, peak1_j]
        cd = corr[peak1_i, peak1_j-1] 
        cu = corr[peak1_i, peak1_j+1]
        
        # gaussian fit
        if np.any(np.array([c,cl,cr,cd,cu]) < 0) and subpixel_method == 'gaussian':
            subpixel_method = 'centroid'
        
        try: 
            if subpixel_method == 'centroid':
                subp_peak_position = (((peak1_i-1)*cl+peak1_i*c+(peak1_i+1)*cr)/(cl+c+cr),
                                    ((peak1_j-1)*cd+peak1_j*c+(peak1_j+1)*cu)/(cd+c+cu))
        
            elif subpixel_method == 'gaussian':
                subp_peak_position = (peak1_i + ((np.log(cl)-np.log(cr))/(2*np.log(cl) - 4*np.log(c) + 2*np.log(cr))),
                                    peak1_j + ((np.log(cd)-np.log(cu))/( 2*np.log(cd) - 4*np.log(c) + 2*np.log(cu)))) 
        
            elif subpixel_method == 'parabolic':
                subp_peak_position = (peak1_i +  (cl-cr)/(2*cl-4*c+2*cr),
                                        peak1_j +  (cd-cu)/(2*cd-4*c+2*cu)) 
    
        except: 
            subp_peak_position = default_peak_position
            
    except IndexError:
            subp_peak_position = default_peak_position
            
    return subp_peak_position[0], subp_peak_position[1]

