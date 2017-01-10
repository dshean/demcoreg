#! /usr/bin/env python

"""
Library of functions that can be used for horizontal co-registration of raster data

These were written in 2012-2013, and should be cleaned up and tested before use

The ASP pc_align ICP co-registration is usually superior to these approaches

"""

import numpy as np
import matplotlib.pyplot as plt

from pygeotools.lib import malib

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

#Return offest for simple sum of absolute differences minimum
#Note - this should work fine for control surfaces
#Will fail when surface is changing and aspect is uniform
def compute_offset_sad(dem1, dem2, pad=(9,9), plot=False):
    """Compute horizontal offset between input rasters using sum of absolute differences (SAD) method
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
    
    #Note, we're dealing with min here, so provide -m for spr
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
        plt.show()

    return m, int_offset, sp_offset

#This is a decent full-image normalized cross-correlation routine with sub-pixel refinement
def compute_offset_ncc(dem1, dem2, pad=(9,9), plot=False): 
    """Compute horizontal offset between input rasters using normalized cross-correlation (NCC) method
    """
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
    #print "Running 2D correlation with search window xrange (%i, %i) and yrange (%i, %i)" % (xsearch + ysearch)
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

    #This was taken from PIRT
    #http://code.google.com/p/pirt/
    #import fitting
    #patch_hw = 4
    #patch = m[offset[0]-patch_hw:offset[0]+patch_hw, offset[1]-patch_hw:offset[1]+patch_hw]
    #patch_offset = np.array(fitting.fit_lq2(patch))
    #print patch_offset
    #sp_offset = offset + patch_offset
    #print sp_offset

    #Find sub-pixel peak
    sp_argmax = np.array(find_subpixel_peak_position(m, 'parabolic'))
    #May need to split this into integer and decimal components, multipy stride*int and add decimal
    #sp_offset = int_offset + (sp_argmax - int_argmax)
    sp_offset = sp_argmax - pad
    #sp_offset = sp_argmax + np.array([ysearch[0], xsearch[0]]) 

    print(sp_offset)

    if plot: 
        fig = plt.figure()
        plt.title('NCC offset, parabolic SPR')
        plt.imshow(m)
        #plt.scatter(*int_argmax[::-1])
        plt.scatter(*sp_argmax[::-1])
        plt.show()

    return m, int_offset, sp_offset, fig

#Bin y by x.
#Returns the binned "y" values and the left edges of the bins
def bin_by(x, y, nbins=360):
    bins = np.linspace(x.min(), x.max(), nbins+1)
    # To avoid extra bin for the max value
    bins[-1] += 1 

    indicies = np.digitize(x, bins)

    output = []
    outlen = []
    for i in range(1, len(bins)):
        output.append(y[indicies==i])
        outlen.append(len(output[-1]))

    outnan = np.empty((nbins, np.max(outlen)))
    outnan.fill(np.nan)

    for n in range(len(output)):
        outnan[n][0:len(output[n])] = output[n]

    outbin = np.ma.masked_invalid(outnan)

    # Just return the left edges of the bins
    bins = bins[:-1]

    plt.figure()
    plt.boxplot(output, sym='')
    plt.xlim(0,360)
    plt.xticks(np.arange(0,360,30))
    plt.ylabel('dh/tan(slope) (m)')
    plt.xlabel('Aspect (1-deg bins)')

    plt.figure()
    plt.bar(bins,outlen)
    plt.xlim(0,360)
    plt.xticks(np.arange(0,360,30))
    plt.ylabel('Count')
    plt.xlabel('Aspect (1-deg bins)')
    plt.show()

    return outbin, bins

#Function for fitting Nuth and Kaab (2011)
def func(x, a, b, c):
    y = a * np.cos(np.deg2rad(b-x)) + c
    #Per Suki suggestion, can use Phasor addition
    #y = a * np.cos(np.deg2rad(x)) + b * np.sin(np.deg2rad(x)) + c
    return y

#This is the Nuth and Kaab (2011) method
def compute_offset_nuth(dh, slope, aspect):
    """Compute horizontal offset between input rasters using Nuth and Kaab [2011] (nuth) method
    """
    import scipy.optimize as optimization

    #mean_dh = dh.mean()
    #mean_slope = slope.mean()
    #c_seed = (mean_dh/np.tan(np.deg2rad(mean_slope))) 
    
    med_dh = np.ma.median(dh)
    med_slope = np.ma.median(slope)
    c_seed = (med_dh/np.tan(np.deg2rad(med_slope))) 
    x0 = np.array([0.0, 0.0, c_seed])
    
    #n = 10000
    #idx = random.sample(range(dh.compressed().size), n)
    #xdata = np.array(aspect.compressed()[idx], float)
    #ydata = np.array((dh/np.tan(np.deg2rad(slope))).compressed()[idx], float) 
    xdata = np.array(aspect.compressed(), float)
    ydata = np.array((dh/np.tan(np.deg2rad(slope))).compressed(), float) 

    vals, bins = bin_by(xdata, ydata)
    
    med_vals = np.ma.median(vals, axis=1)   
    med_vals = np.ma.masked_where(np.around(bins) % 45 == 0, med_vals)
    bins = np.ma.masked_where(np.around(bins) % 45 == 0, bins)

    #badbins = [0, 45, 90, 180, 225, 270, 315]
    #vals[badbins] = np.nan
    #bins[badbins] = np.nan
    #vals = np.ma.fix_invalid(vals)
    #bins = np.ma.fix_invalid(bins)
 
    fit = optimization.curve_fit(func, xdata, ydata, x0)[0]
    print(fit) 
    genplot(xdata, ydata, fit) 
    
    #hist, xedges, yedges = np.histogram2d(xdata, ydata, bins=[360,1])

    #Compute median absolute difference, better than stddev?
    y_med = np.median(ydata)
    y_mad = malib.mad(ydata)
    mad_factor = 2
    y_perc = [y_med - y_mad*mad_factor, y_med + y_mad*mad_factor]

    y_idx = ((ydata >= y_perc[0]) & (ydata <= y_perc[1]))
    ydata_clip = ydata[y_idx]
    xdata_clip = xdata[y_idx]

    #Generate synthetic data to test curve_fit
    #xdata = np.arange(0,360,0.01)
    #ydata = f(xdata, 20.0, 130.0, -3.0) + 20*np.random.normal(size=len(xdata))
    
    fit = optimization.curve_fit(func, xdata_clip, ydata_clip, x0)[0]
    print(fit)
    genplot(xdata_clip, ydata_clip, fit) 
    
    fit = optimization.curve_fit(func, bins, med_vals, x0)[0]
    print(fit)
    genplot(bins, np.ma.median(vals, axis=1), fit) 
    genplot(xdata, ydata, fit) 
    
    return fit

def genplot(x, y, fit):
    import random
    a = (np.arange(0,360))
    f_a = func(a, fit[0], fit[1], fit[2])
    if x.size > 10000:
        idx = random.sample(list(range(x.size)), 10000)
    else:
        idx = np.arange(x.size)
    plt.figure()
    plt.xlabel('Aspect (deg)')
    plt.ylabel('dh/tan(slope) (m)')
    plt.plot(x[idx], y[idx], 'r.')
    plt.xlim(0,360)
    plt.axhline(color='k')
    plt.plot(a, f_a, 'b')
    #plt.ylim(-200,200)
    plt.show()
