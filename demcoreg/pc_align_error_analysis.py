#! /usr/bin/env python

#David Shean
#dshean@gmail.com

#Perform error analysis for DEM output from pc_align

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from pygeotools.lib import timelib, geolib, malib

#shell
#filenames=$(ls WV*/dem*/*align/*-DEM.tif)
#error_analysis.py $filenames

#ipython
#filenames = !ls *align/*-DEM.tif
#run ~/src/demtools/error_analysis.py $filenames.s

#Compute standard error for all DEMs
#SE=SD/sqrt(n)

def shift_ll2proj(fn, llz):
    from osgeo import gdal, osr
    from pygeotools.lib import geolib
    ds = gdal.Open(fn)
    s_srs = geolib.wgs_srs
    t_srs = geolib.get_ds_srs(ds)
    shift = None
    if t_srs is not None and not s_srs.IsSame(t_srs):
        #center is lon, lat
        #llz is lat, lon
        c = geolib.get_center(ds, t_srs=s_srs)
        c_shift = [c[0]+llz[1], c[1]+llz[0]]
        ct = osr.CoordinateTransformation(s_srs, t_srs)
        c_proj = list(ct.TransformPoint(*c)[0:2])
        c_shift_proj = list(ct.TransformPoint(*c_shift)[0:2])
        shift = list([c_shift_proj[0] - c_proj[0], c_shift_proj[1] - c_proj[1]])
        shift.append(llz[2])
    return shift

def parse_pc_align_log(fn):
    import re
    error_dict = None
    #Determine log filename
    import glob
    log_fn = glob.glob(fn.rsplit('-DEM', 1)[0]+'*.log')
    if not log_fn:
        log_fn = glob.glob(fn.rsplit('-DEM', 1)[0]+'*align/*.log')

    if not log_fn:
        print "Failed to locate align log for %s" % fn
    else:
        log_fn = log_fn[0]
        print(log_fn)
        f = open(log_fn)

        error_dict = {}
        error_dict['File'] = fn
        error_dict['Date'] = timelib.fn_getdatetime(fn)

        #This handles cases where no sampling was performed
        error_dict['Input Sampled 16th Percentile Error'] = np.nan 
        error_dict['Input Sampled Median Error'] = np.nan 
        error_dict['Input Sampled 84th Percentile Error'] = np.nan 
        error_dict['Input Sampled Error Spread'] = np.nan 
        error_dict['Output Sampled 16th Percentile Error'] = np.nan 
        error_dict['Output Sampled Median Error'] = np.nan 
        error_dict['Output Sampled 84th Percentile Error'] = np.nan 
        error_dict['Output Sampled Error Spread'] = np.nan 
        #error_dict['Translation vector (North-East-Down, meters)'] = [np.nan, np.nan, np.nan]

        #Set default reference type to point
        error_dict['Ref type'] = 'point'

        temp = []
        for line in f:
            key = 'Loaded points'
            if key in line:
                temp.append(int(re.split(':', line.rstrip())[1]))
            key = 'Number of errors'
            if key in line:
                error_dict[key] = int(re.split(':', line.rstrip())[1])
            key = 'Input: error percentile'
            if key in line:
                line_a = re.split(': |, ', line.rstrip())
                error_dict['Input 16th Percentile Error'] = float(line_a[3])
                error_dict['Input Median Error'] = float(line_a[5])
                error_dict['Input 84th Percentile Error'] = float(line_a[7])
            """
            key = 'Input: error mean'
            if key in line:
                line_a = re.split(': |, ', line.rstrip())
                error_dict['Input Mean Error'] = float(line_a[2])
                error_dict['Input Std Error'] = float(line_a[4])
            """
            #This pulls the line 
            #Input: mean of smallest errors: 25%: 7.82061, 50%: 9.71931, 75%: 10.9917, 100%: 12.2715
            #Want the final value
            key = 'Input: mean'
            if key in line:
                line_a = re.split(': |, ', line.rstrip())
                error_dict['Input Mean Error'] = float(line_a[-1])
            key = 'Output: error percentile'
            if key in line:
                line_a = re.split(': |, ', line.rstrip())
                error_dict['Output 16th Percentile Error'] = float(line_a[3])
                error_dict['Output Median Error'] = float(line_a[5])
                error_dict['Output 84th Percentile Error'] = float(line_a[7])
            """
            key = 'Output: error mean'
            if key in line:
                line_a = re.split(': |, ', line.rstrip())
                error_dict['Output Mean Error'] = float(line_a[2])
                error_dict['Output Std Error'] = float(line_a[4])
            """
            key = 'Output: mean'
            if key in line:
                line_a = re.split(': |, ', line.rstrip())
                error_dict['Output Mean Error'] = float(line_a[-1])
            key = 'Translation vector (Cartesian, meters)'
            #Previous versions of pc_align output this
            #key = 'Translation vector (meters)'
            if key in line:
                error_dict['Translation vector (Cartesian, meters)'] = list(float(i) for i in re.split('Vector3\(', line.rstrip())[1][:-1].split(',')) 
                #error_dict['Translation vector (meters)'] = list(float(i) for i in re.split('Vector3\(', line.rstrip())[1][:-1].split(',')) 
            key = 'Translation vector (North-East-Down, meters)'
            if key in line:
                error_dict['Translation vector (North-East-Down, meters)'] = list(float(i) for i in re.split('Vector3\(', line.rstrip())[1][:-1].split(',')) 
            key = 'Translation vector magnitude (meters)'
            if key in line:
                error_dict[key] = float(re.split(':', line.rstrip())[1])
            key = 'Translation vector (lat,lon,z)'
            if key in line:
                error_dict[key] = list(float(i) for i in re.split('Vector3\(', line.rstrip())[1][:-1].split(',')) 
                shift_proj = shift_ll2proj(fn, error_dict[key])
                key = 'Translation vector (Proj meters)'
                error_dict[key] = shift_proj

            #This is the output from the point sampling post-alignment
            key = 'Error percentiles'
            if key in line:
                #This is a hack to capture both sampling of input and output
                if 'Output Sampled 16th Percentile Error' in error_dict:
                    error_dict['Input Sampled 16th Percentile Error'] = error_dict['Output Sampled 16th Percentile Error']
                    error_dict['Input Sampled Median Error'] = error_dict['Output Sampled Median Error']
                    error_dict['Input Sampled 84th Percentile Error'] = error_dict['Output Sampled 84th Percentile Error']
                    error_dict['Input Sampled Error Spread'] = error_dict['Output Sampled Error Spread']
                line_a = re.split(': |, ', line.rstrip())
                error_dict['Output Sampled 16th Percentile Error'] = float(line_a[2])
                error_dict['Output Sampled Median Error'] = float(line_a[4])
                error_dict['Output Sampled 84th Percentile Error'] = float(line_a[6])
                error_dict['Output Sampled Error Spread'] = float(line_a[6]) - float(line_a[2])
            #key = 'compute_dh'
            #Note: these are not computed for absolute values by compute_dh
            key = 'count:'
            if key in line:
                error_dict['Ref type'] = 'grid'
                #This is a hack to capture both sampling of input and output
                if 'Output Sampled 16th Percentile Error' in error_dict:
                    error_dict['Input Sampled 16th Percentile Error'] = error_dict['Output Sampled 16th Percentile Error']
                    error_dict['Input Sampled Median Error'] = error_dict['Output Sampled Median Error']
                    error_dict['Input Sampled 84th Percentile Error'] = error_dict['Output Sampled 84th Percentile Error']
                    error_dict['Input Sampled Error Spread'] = error_dict['Output Sampled Error Spread']
                #Assume the following format for stats:
                #count: 349835 min: -51.39 max: 22.00 mean: 0.29 std: 0.49 med: 0.28 mad: 0.37 \
                #q1: 0.04 q2: 0.54 iqr: 0.50 mode: 0.29 p16: -0.07 p84: 0.66 spread: 0.37
                line_a = re.split(': | ', line.rstrip())
                error_dict['Output Sampled 16th Percentile Error'] = float(line_a[23])
                error_dict['Output Sampled Median Error'] = float(line_a[11])
                error_dict['Output Sampled 84th Percentile Error'] = float(line_a[25])
                error_dict['Output Sampled Error Spread'] = float(line_a[25]) - float(line_a[23])
            key = 'Mean error'
            if key in line:
                if 'Output Sampled Mean Error' in error_dict:
                    error_dict['Input Sampled Mean Error'] = error_dict['Output Sampled Mean Error']
                error_dict['Output Sampled Mean Error'] = float(re.split(':', line.rstrip())[1])
            key = 'RMSE'
            if key in line:
                if 'Output Sampled RMSE' in error_dict:
                    error_dict['Input Sampled RMSE'] = error_dict['Output Sampled RMSE']
                error_dict['Output Sampled RMSE'] = float(re.split(':', line.rstrip())[1])
            key = 'Absolute Median Error'
            if key in line:
                if 'Output Absolute Median Error' in error_dict:
                    error_dict['Input Absolute Median Error'] = error_dict['Output Absolute Median Error']
                error_dict['Output Absolute Median Error'] = float(re.split(':', line.rstrip())[1])
                
        error_dict['Source points'] = temp[0] 
        error_dict['Reference points'] = temp[1] 
    
    return error_dict

def make_plot(x,y,yerr=None,c='k',ms=4,label=None,abs=False):
    y_mean = y.mean()
    y_std = y.std()
    y_med = np.ma.median(y)
    y_nmad = malib.mad(y)
    #plt.plot(x, y, label=label, color=c, marker='o', linestyle='None')
    plt.scatter(x, y, label=label, color=c, marker='o', s=ms)
    if yerr is not None:
        plt.errorbar(x, y, yerr=yerr, color=c, linestyle='None', elinewidth=0.5, capsize=np.sqrt(ms), alpha=0.5)
    plt.axhline(y_med, color=c, linestyle='--', alpha=0.5)
    plt.axhline(y_med + y_nmad, color=c, linestyle=':', alpha=0.5)
    plt.axhline(y_med - y_nmad, color=c, linestyle=':', alpha=0.5)
    plt.axhline(0, color='k', linewidth=0.5, linestyle='-', alpha=0.5)
    ax = plt.gca()
    plt.minorticks_on()
    #ax.tick_params(axis='y',which='minor',left='on')
    if abs:
        ax.set_ylim(bottom=0.0)

#Draw ellipsoid in 3D plot
#http://stackoverflow.com/questions/7819498/plotting-ellipsoid-with-matplotlib
def make_plot3d(x, y, z, title=None, orthogonal_fig=True):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')
    ax.set_xlabel('X offset (m)')
    ax.set_ylabel('Y offset (m)')
    ax.set_zlabel('Z offset (m)')
    if title is not None:
        plt.suptitle(title)
    ax.plot(x, y, z, 'o')

    cmean = np.mean([x,y,z], axis=1)
    cmed = np.median([x,y,z], axis=1)

    ax.scatter(cmean[0], cmean[1], cmean[2], color='r', marker='s')

    ce90 = geolib.CE90(x,y)
    le90 = geolib.LE90(z)
    coefs = [ce90, ce90, le90]
    ax.set_title("CE90: %0.2f, LE90: %0.2f, n=%i" % (ce90, le90, x.shape[0])) 
    
    maxdim = np.ceil(np.max([np.max(np.abs([x, y, z])), ce90, le90]))
    ax.set_xlim(-maxdim, maxdim)
    ax.set_ylim(-maxdim, maxdim)
    ax.set_zlim(-maxdim, maxdim)

    rx, ry, rz = coefs
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    ex = rx * np.outer(np.cos(u), np.sin(v))
    ey = ry * np.outer(np.sin(u), np.sin(v))
    ez = rz * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(ex, ey, ez, rstride=2, cstride=2, linewidth=0, color='b', alpha=0.1)
    #max_radius = max(rx, ry, rz)
    #for axis in 'xyz':
    #    getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))
    if orthogonal_fig:
        from matplotlib.patches import Ellipse
        fig_ortho = plt.figure(figsize=(10,4))
        #fig_ortho = plt.figure()
        title='ICP Alignment Translation Vectors\nn=%i, mean: (%0.2f, %0.2f, %0.2f)\nCE90: %0.2f, LE90: %0.2f' % (x.shape[0], cmean[0], cmean[1], cmean[2], ce90, le90)
        plt.suptitle(title) 

        ax = fig_ortho.add_subplot(131)
        ax.plot(x, y, color='b', linestyle='None', marker='o', label='ICP correction vector')
        ax.plot(cmean[0], cmean[1], color='r', linestyle='None', marker='s', label='Mean')
        #ax.scatter(x, y)
        #ax.scatter(cmean[0], cmean[1], color='r', marker='s')
        ax.set_xlim(-maxdim, maxdim)
        ax.set_ylim(-maxdim, maxdim)
        ax.minorticks_on()
        ax.set_aspect('equal')
        ax.set_xlabel('X offset (m)')
        ax.set_ylabel('Y offset (m)')
        e = Ellipse((0,0), 2*ce90, 2*ce90, linewidth=0, alpha=0.1)
        ax.add_artist(e)
        plt.legend(prop={'size':8}, numpoints=1, loc='upper left')

        ax = fig_ortho.add_subplot(132)
        ax.plot(x, z, color='b', linestyle='None', marker='o', label='ICP correction vector')
        ax.plot(cmean[0], cmean[2], color='r', linestyle='None', marker='s', label='Mean')
        #ax.scatter(x, z)
        #ax.scatter(cmean[0], cmean[2], color='r', marker='s')
        ax.set_xlim(-maxdim, maxdim)
        ax.set_ylim(-maxdim, maxdim)
        ax.minorticks_on()
        ax.set_aspect('equal')
        ax.set_xlabel('X offset (m)')
        ax.set_ylabel('Z offset (m)')
        e = Ellipse((0,0), 2*ce90, 2*le90, linewidth=0, alpha=0.1)
        ax.add_artist(e)

        ax = fig_ortho.add_subplot(133)
        ax.plot(y, z, color='b', linestyle='None', marker='o', label='ICP correction vector')
        ax.plot(cmean[1], cmean[2], color='r', linestyle='None', marker='s', label='Mean')
        #ax.scatter(y, z)
        #ax.scatter(cmean[1], cmean[2], color='r', marker='s')
        ax.set_xlim(-maxdim, maxdim)
        ax.set_ylim(-maxdim, maxdim)
        ax.minorticks_on()
        ax.set_aspect('equal')
        ax.set_xlabel('Y offset (m)')
        ax.set_ylabel('Z offset (m)')
        e = Ellipse((0,0), 2*ce90, 2*le90, linewidth=0, alpha=0.1)
        ax.add_artist(e)
        
        plt.tight_layout()

        #Note: postscript doesn't properly handle tansparency
        #fig_fn = 'icp_translation_vec_proj_meters_orthogonal.pdf'
        fig_fn = 'icp_translation_vec_local_meters_orthogonal.pdf'
        plt.savefig(fig_fn, dpi=600, bbox_inches='tight')

    #Set back to original figure
    plt.figure(fig.number)

def main():
    #filenames = !ls *align/*reference-DEM.tif
    #run ~/src/demtools/error_analysis.py $filenames.s

    if len(sys.argv) < 1:
        sys.exit('No input files provided')

    fn_list = sys.argv[1:]
    n_samp = len(fn_list)

    error_dict_list = []
    for fn in fn_list:
        ed = parse_pc_align_log(fn)
        if 'Translation vector (North-East-Down, meters)' in ed.keys(): 
            error_dict_list.append(ed)  

    import matplotlib.dates as mdates
    #This is used for interactive display of x-value in plot window 
    date_str = '%Y-%m-%d %H:%M'
    date_fmt = mdates.DateFormatter(date_str)
    #ax.fmt_xdata = mdates.DateFormatter(date_fmt)
    months = mdates.MonthLocator() 
    months_int = mdates.MonthLocator(interval=6)  # every n months 
    years = mdates.YearLocator()   # every year
    yearsFmt = mdates.DateFormatter('%Y')
    #ax.xaxis.set_major_formatter(yearsFmt)
    #ax.xaxis.set_major_locator(months_int3)

    print
    print "n:", len(error_dict_list) 

    """
    #ECEF translations
    #key = 'Translation vector (ECEF meters)'
    key = 'Translation vector (Cartesian, meters)'
    #key = 'Translation vector (meters)'
    val = np.array([e[key] for e in error_dict_list])
    #make_plot3d(val[:,0], val[:,1], val[:,2], title=key)
    ce90 = geolib.CE90(val[:,0], val[:,1])
    le90 = geolib.LE90(val[:,2])
    print
    print key
    print "CE90:", ce90 
    print "LE90:", le90 
    print

    #Proj translation
    key = 'Translation vector (Proj meters)' 
    val = np.array([e[key] for e in error_dict_list])
    ce90 = geolib.CE90(val[:,0], val[:,1])
    le90 = geolib.LE90(val[:,2])
    print
    print key
    print "CE90:", ce90 
    print "LE90:", le90 
    print
    print 'Centroid (mean) of offsets (Proj meters): ', np.mean(val, axis=0) 
    print 'Centroid (median) of offsets (Proj meters): ', np.median(val, axis=0) 
    """

    #NOTE: changed default to N-E-D on 9/18/15
    #Can have significant differences for local proj vs. polar stereographic proj
    #Should regenerate all previous figures

    #Local translation on ellipsoid
    #This appears to be local stereographic projection on ellipsoid
    key = 'Translation vector (North-East-Down, meters)'
    val = np.array([e[key] for e in error_dict_list])

    #Reformat (n, e, +d) for (x, y, +z) coord sys
    val[:,[0,1]] = val[:,[1,0]]
    val[:,2] *= -1
    ce90 = geolib.CE90(val[:,0], val[:,1])
    le90 = geolib.LE90(val[:,2])
    print
    print key
    print "CE90:", ce90 
    print "LE90:", le90 
    print
    print 'Centroid (mean) of offsets (local ned meters): ', np.mean(val, axis=0) 
    print 'Centroid (median) of offsets (local ned meters): ', np.median(val, axis=0) 

    #Remove vertical bias
    remove_vertbias = False 
    if remove_vertbias:
        print "Removing vertical bias: %0.2f" % np.mean(val, axis=0)[2]
        val[:,2] -= np.mean(val, axis=0)[2]

    remove_outliers = False 
    #Flag outliers
    x_mag = val[:,0]
    y_mag = val[:,1]
    h_mag = np.sqrt(val[:,0]**2 + val[:,1]**2)
    v_mag = val[:,2]
    mag = np.sqrt(val[:,0]**2 + val[:,1]**2 + val[:,2]**2)
    abs_thresh = 10.0
    p = 98.0
    p_thresh = np.percentile(h_mag, p)
    #print "Outliers with horiz error >= %0.2f (%0.1f%%)" % (p_thresh, p)
    print "Outliers:" 
    #idx = (h_mag >= p_thresh).nonzero()[0]
    idx = (h_mag >= ce90).nonzero()[0]
    idx = np.unique(np.hstack([idx, ((np.abs(v_mag) >= le90).nonzero()[0])]))

    #Print all
    #idx = np.arange(h_mag.size)
    #idx_sort = np.argsort(mag[idx])
    #idx = idx[idx_sort]

    print 'name, m, h, v, x, y, z'
    for i in idx:
        print error_dict_list[i]['File'], mag[i], h_mag[i], v_mag[i], val[i,0:3]
        #Delete from list
        if remove_outliers:
            print "Removing from calculation"
            del error_dict_list[i]

    if remove_vertbias or remove_outliers:
        print
        print "Updated values"
        print key
        print "CE90:", geolib.CE90(val[:,0], val[:,1])
        print "LE90:", geolib.LE90(val[:,2])
        print
        print 'Centroid (mean) of offsets (local ned meters): ', np.mean(val, axis=0) 
        print 'Centroid (median) of offsets (local ned meters): ', np.median(val, axis=0) 

    #Extract dates
    date_vec = np.array([e['Date'] for e in error_dict_list])
    x = date_vec

    make_plot3d(val[:,0], val[:,1], val[:,2], title=key)
    #Note: there is a bug in pdf that displayes surface lines
    #fig_fn = 'icp_translation_vec_proj_meters.png'
    fig_fn = 'icp_translation_vec_local_meters.png'
    #plt.savefig(fig_fn, dpi=600, bbox_inches='tight')

    fig, ax = plt.subplots(1)
    key = 'Translation vector (lat,lon,z)'
    plt.title('ICP translation vector (lat,lon,z): Z component')
    val = np.array([e[key] for e in error_dict_list])
    y = val[:,2]
    make_plot(x,y,c='b',label=key, abs=False)
    fig.autofmt_xdate()
    ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_major_locator(months_int)
    #ax.xaxis.set_major_formatter(date_fmt)
    ax.fmt_xdata = date_fmt
    ax.set_ylabel('Z offset (m)')

    fig, ax = plt.subplots(1)
    key = 'Translation vector magnitude (meters)'
    plt.title('ICP Translation vector magnitude (meters)')
    y = np.array([e[key] for e in error_dict_list])
    make_plot(x,y,c='b',label=key, abs=True)
    fig.autofmt_xdate()
    ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_major_locator(months_int)
    #ax.xaxis.set_major_formatter(date_fmt)
    ax.fmt_xdata = date_fmt
    ax.set_ylabel('Offset (m)')

    fig, ax = plt.subplots(1)
    key = 'Number of errors'
    plt.title('Number of error samples')
    nerr = np.array([e[key] for e in error_dict_list])
    make_plot(x,nerr,c='b',label=key)
    fig.autofmt_xdate()
    ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_major_locator(months_int)
    #ax.xaxis.set_major_formatter(date_fmt)
    ax.fmt_xdata = date_fmt
    ax.set_ylabel('N samples')

    """
    fig, ax = plt.subplots(1)
    plt.title('ICP Standard Deviation')
    key = 'Input Std Error'
    in_std = np.array([e[key] for e in error_dict_list])
    make_plot(x,in_std,c='r',label=key)
    key = 'Output Std Error'
    out_std = np.array([e[key] for e in error_dict_list])
    make_plot(x,out_std,c='b',label=key)
    fig.autofmt_xdate()
    ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_major_locator(months_int)
    #ax.xaxis.set_major_formatter(date_fmt)
    ax.fmt_xdata = date_fmt
    plt.legend(scatterpoints=1)
    """

    fig, ax = plt.subplots(1)
    plt.title('ICP Mean Error')
    key = 'Input Mean Error'
    in_mean = np.array([e[key] for e in error_dict_list])
    #make_plot(x,in_mean,c='r',label=key,yerr=in_std)
    make_plot(x,in_mean,c='r',label=key, abs=True)
    key = 'Output Mean Error'
    out_mean = np.array([e[key] for e in error_dict_list])
    #make_plot(x,out_mean,c='b',label=key,yerr=out_std)
    make_plot(x,out_mean,c='b',label=key, abs=True)
    fig.autofmt_xdate()
    ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_major_locator(months_int)
    #ax.xaxis.set_major_formatter(date_fmt)
    ax.fmt_xdata = date_fmt
    ax.set_ylabel('Mean error (m)')
    plt.legend(scatterpoints=1, loc='upper left', prop={'size':8})

    fig, ax = plt.subplots(1)
    plt.title('ICP Median Error')
    key = 'Input 16th Percentile Error'
    in_16p = np.array([e[key] for e in error_dict_list])
    key = 'Input 84th Percentile Error'
    in_84p = np.array([e[key] for e in error_dict_list])
    key = 'Input Median Error'
    in_med = np.array([e[key] for e in error_dict_list])
    make_plot(x,in_med,c='r',label=key,yerr=[in_med - in_16p, in_84p - in_med], abs=True)
    key = 'Output 16th Percentile Error'
    out_16p = np.array([e[key] for e in error_dict_list])
    key = 'Output 84th Percentile Error'
    out_84p = np.array([e[key] for e in error_dict_list])
    key = 'Output Median Error'
    out_med = np.array([e[key] for e in error_dict_list])
    make_plot(x,out_med,c='b',label=key,yerr=[out_med - out_16p, out_84p - out_med], abs=True)
    fig.autofmt_xdate()
    ax.fmt_xdata = mdates.DateFormatter(date_fmt)
    ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_major_locator(months_int)
    #ax.xaxis.set_major_formatter(date_fmt)
    ax.fmt_xdata = date_fmt
    ax.set_ylabel('Median error (m)')
    plt.legend(scatterpoints=1, loc='upper left', prop={'size':8})
    fig_fn = 'icp_median_error.pdf'
    plt.savefig(fig_fn, dpi=600, bbox_inches='tight')

    fig, ax = plt.subplots(1)
    plt.title('Sampled Median Error')
    key = 'Input Sampled 16th Percentile Error'
    in_16p = np.ma.fix_invalid([e[key] for e in error_dict_list])
    if in_16p.count() > 0:
        key = 'Input Sampled 84th Percentile Error'
        in_84p = np.ma.fix_invalid([e[key] for e in error_dict_list])
        key = 'Input Sampled Median Error'
        in_med = np.ma.fix_invalid([e[key] for e in error_dict_list])
        in_spread = in_84p - in_16p
        make_plot(x,in_med,c='r',label=key,yerr=[in_med - in_16p, in_84p - in_med], abs=True)

        key = 'Output Sampled 16th Percentile Error'
        out_16p = np.ma.fix_invalid([e[key] for e in error_dict_list])
        key = 'Output Sampled 84th Percentile Error'
        out_84p = np.ma.fix_invalid([e[key] for e in error_dict_list])
        key = 'Output Sampled Median Error'
        out_med = np.ma.fix_invalid([e[key] for e in error_dict_list])
        out_spread = out_84p - out_16p

        p = 90.0
        out_med_thresh = np.percentile(out_med, p)
        out_spread_thresh = np.percentile(out_spread, p)
        #print "Outliers with horiz error >= %0.2f (%0.1f%%)" % (p_thresh, p)
        print
        print "Sampled Error Outliers:" 
        #idx = (h_mag >= p_thresh).nonzero()[0]
        idx = (out_med >= out_med_thresh).nonzero()[0]
        idx = np.unique(np.hstack([idx, ((out_spread >= out_spread_thresh).nonzero()[0])]))
        #Print all
        idx = np.arange(out_med.size)
        idx_sort = np.argsort(out_med[idx])
        idx = idx[idx_sort]
        print 'name, samp_mederrr, samp_errspread, nerr'
        for i in idx:
            print error_dict_list[i]['File'], out_med[i], out_spread[i], nerr[i]
            #Delete from list
            if remove_outliers:
                print "Removing from calculation"
                del error_dict_list[i]
        print
        print 'Input sampled median error (spread/2): %0.2f (%0.2f)' % (np.median(in_med), np.median(in_spread)/2.)
        print 'Output sampled median error (spread/2): %0.2f (%0.2f)' % (np.median(out_med), np.median(out_spread)/2.)
        print

        make_plot(x,out_med,c='b',label=key,yerr=[out_med - out_16p, out_84p - out_med], abs=True)
        fig.autofmt_xdate()
        ax.set_ylabel('Median error (m)')
        ax.fmt_xdata = mdates.DateFormatter(date_fmt)
        ax.xaxis.set_minor_locator(months)
        #ax.xaxis.set_major_locator(months_int)
        #ax.xaxis.set_major_formatter(date_fmt)
        ax.fmt_xdata = date_fmt
        ax.set_ylabel('Median error (m)')
        plt.legend(scatterpoints=1, loc='upper left', prop={'size':8})
        ax.set_ylim(-15,15)
        fig_fn = 'sampled_median_error.pdf'
        #fig_fn = 'sampled_median_error_2014-2016.pdf'
        #from datetime import datetime
        #ax.set_xlim(datetime(2014,1,1),datetime(2016,7,1))
        plt.savefig(fig_fn, dpi=600, bbox_inches='tight')
        #plt.show()

if __name__ == '__main__':
    main()
