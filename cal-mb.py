#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# cal-mb.py
#
# example usage:
# 
# python3 cal-mb.py 57Fe_calib_raw_data.ws5
# or
# python3 cal-mb.py 57Fe_calib_raw_data.ws5 -s
#
# '-s' to show the plot window 
#
# '-fl' to fold the raw spectrum to the left


##########################################################################################
import sys                                  #sys
import os                                   #os file processing
import datetime                             #calc mod. time of .ws5 from timestamp
import argparse                             #argument parser

import matplotlib.pyplot as plt             #plots
import numpy as np                          #summation and other math

from lmfit.models import LorentzianModel, ConstantModel         #fit
from lmfit import Parameters                                    #fit

from scipy.signal import find_peaks, peak_prominences           #peak finding
from scipy import interpolate               #interpolation of channel intens. for folding
##########################################################################################

#for windows console
sys.stdout.reconfigure(encoding = 'utf-8') #unicode

##########################################################################################
#import data                                                                             #
########################################################################################## 
def op_imp(data_file):
    try:
        #read raw data 
        ws5_raw_data = np.loadtxt(data_file, comments=['#','<'])
        mod_date = datetime.datetime.fromtimestamp(os.path.getmtime(data_file))
    except IOError:
        #file not found -> exit here
        print(f"'{data_file}'" + " not found. Exit.")
        sys.exit(1)
    except ValueError:
        print('Warning! Wrong data format. Exit.')
        sys.exit(1)
    #return array with intensity data from MCA (.ws5), filename and modification 
    #date and time of the .ws5 file 
    return ws5_raw_data, data_file, mod_date.strftime("%d.%m.%Y %H:%M:%S")      

##########################################################################################
#fit data                                                                                #
########################################################################################## 
def do_the_fit(data, mean_stdev_fold_i = 1):
    #fit data with N Lorentz functions; N = number of detected peaks in raw data
    #get the number of channels
    N_chan = len(data)
    #channels from 1 to max channel, step = 1 channel
    chan = np.linspace(1, N_chan, N_chan)
    #find peaks in raw data; has to be '-' since they are minima (not maxima)
    peaks, _ = find_peaks(-data)
    #get prominences = peak separation from noise
    prominences = peak_prominences(-np.array(data), peaks)[0]
    #identify peaks from noise
    #normal max(prominences)/3; found an example where only max(prominences)/4 worked
    peaks, _ = find_peaks(-data, prominence=(max(prominences)/4, max(prominences)))
    #generate a list of Lorentz functions
    lorentzlist = list()
    #initialize parameters for lmfit
    params = Parameters()
    
    #number of Lorentz functions = number of detected peaks 
    for index, peak in enumerate(peaks): 
        #define Lorentz function for lmfit
        lorentz = LorentzianModel(prefix='s'+str(index)+"_") 
        #guess parameters for Lorentz function 
        params.update(lorentz.guess(data, x = chan))
        #update parameters for Lorentz function, x = x from detected peak (channel), 
        #amplitude = y from detected peak (intensity of raw data)
        params.update(lorentz.make_params(center = chan[peak], amplitude = data[peak]))
        #add to list of Lorentz functions
        lorentzlist.append(lorentz)
    
    #constant model for bg or offset (y0)
    bg = ConstantModel()
    #start is max intensity from raw data
    bg_params = bg.make_params(c = {'value':max(data), 'vary':True})
    
    #complete model for lmfit: list of Lorentz functions + offset
    curve = np.sum(lorentzlist) + bg
    #do the fit with weights in case of the folded spectrum
    #weights = 1 in case of the unfolded spectrum
    fit_result = curve.fit(data, params + bg_params, x = chan, \
       weights = 1/mean_stdev_fold_i, scale_covar=True)
    #return fit results
    return fit_result
    
##########################################################################################
#calc FP or v0                                                                           #
##########################################################################################
def calc_FP_v0(fit_result):
    #calc FP or v0 from the mean of the centers of the fitted Lorentz functions
    #a way to estimate the number of Lorentz functions (and centers) from the fit results
    n_center = len([s for s in fit_result.var_names if 'center' in s])
    #generate a list of centers
    centerlist=list()
    for index in range(n_center):
        centerkey = 's'+ str(index) + '_center'
        #get the postion (channel number) of the center of each Lorentz function
        centerlist.append(fit_result.uvars[centerkey])
    #mean of all centers
    FP_v0 = np.mean(centerlist)
    #return FP (folding point) or v0 (channel where the velocity is zero 
    #and the list with centers (for later labeling in the plot)
    return FP_v0, centerlist

##########################################################################################
#fold spectrum                                                                           #
##########################################################################################
def fold_spec(data, FP):
    #fold the spectrum
    #folding_diff = (FP.nominal_value - 256.5)*2
    #get the number of channels
    N_chan = len(data)
    #'(FP - 256.5)*2' for 512 channels, if channel 1 is 1 (and not zero)
    folding_diff = (FP.nominal_value - (int(N_chan/2)+0.5))*2
    #found an example where folding_diff < 0; abs() correct?
    if folding_diff < 0:
        folding_diff = abs(folding_diff)
    #channels from 1 to max channel, step 1 channel 
    chan = np.linspace(1, N_chan, N_chan)
    #interpolate channels, to operate with channel floating point numbers (xxx.xx)
    data_ichan = interpolate.interp1d(chan, data, bounds_error=True, kind = 'linear')
    #lhs (left hand side) channels; note that it goes from high to low
    lhs_chan = np.linspace(int(N_chan/2), 1, int(N_chan/2))
    #rhs (right hand side) channels 
    rhs_chan = np.linspace(int(N_chan/2) + 1, N_chan, int(N_chan/2))
    #folding to the left
    if args.foldleft:
        #lhs (left hand side) channels; note that it goes from low to high
        lhs_chan = np.linspace(1, int(N_chan/2), int(N_chan/2))
        #rhs (right hand side) channels; note that it goes from high to low
        rhs_chan = np.linspace(N_chan, int(N_chan/2) + 1, int(N_chan/2))
    #add the intensities of lhs + folding difference and rhs channels pairwise
    folded_intens = (np.add(data_ichan(lhs_chan+folding_diff), data_ichan(rhs_chan)))
    
    #error from folding
    #assuming that lhs intensity should be equal to rhs intensity
    #since the sum of lhs intensity and rhs intensity is utilized, lhs and rhs intensity
    #are doubled (*2)
    lhs_i_2x = data_ichan(lhs_chan + folding_diff) * 2
    rhs_i_2x = data_ichan(rhs_chan) * 2
    #stddev for error bar plot in ax0
    stdev_fold_i = np.std([lhs_i_2x, rhs_i_2x], axis = 0)
    #mean stddev (one value) from sqrt of the mean of variances of 
    #lhs and rhs intensity * 2
    #the mean stddev is used as weight in the fit and for the residuals
    mean_stdev_fold_i = np.sqrt(np.mean(np.var([lhs_i_2x, rhs_i_2x], axis = 0)))
    
    #return intensities of the folded spectrum ant the mean stddev 
    return folded_intens, mean_stdev_fold_i

##########################################################################################
#calc vmax                                                                               #
##########################################################################################
def calc_vmax(fit_result, data):
    #calc vmax
    #get the number of channels
    N_chan = len(data)
    #multiplicator for velocity per channel (f) to get vmax; it is 127.5 for 256 channels
    #from the manual of mcal from E.B.
    chan_mul = ((N_chan) / 2 - 1) / 2               #127.5 for 265
    #a way to estimate the number of Lorentz functions (and centers) from the fit results
    n_center = len([s for s in fit_result.var_names if 'center' in s])
    #generate a list of centers
    centerlist=list()
    
    for index in range(n_center):
        centerkey = 's'+ str(index) + '_center'
         #get the postion (channel number) of the center of each Lorentz function
        centerlist.append(fit_result.uvars[centerkey])
    
    #from the manual of mcal from E.B.
    #the well known quadrupole splitting values of 57Fe
    #
    #---------------------------------------------------
    #   |       |       |         |       |       |
    #   |       |       |<-1.667->|       |       |
    #   |       |                         |       |
    #   |       |<---------6.167--------->|       |
    #   |                                         |
    #   |<----------------10.657----------------->|
    #   |                  mm/s                   |
    #
    #f is the velocity / channel
    if n_center == 6:
        #folded spectrum with 3 doublets 
        #but here it is the difference of the single Lorentz functions form the outermost
        #to the innermost pair
        #no recalculation with doublets like in mcal from E.B.
        #could not observe a large difference in the final parameters 
        f = 10.657 / (centerlist[0] - centerlist[5]) + \
             6.167 / (centerlist[1] - centerlist[4]) + \
             1.677 / (centerlist[2] - centerlist[3])
        f = (f / 3) # mean f
        #f * chan_mul is vmax
        vmax = f * chan_mul
    elif n_center == 4: 
        #folded spectrum with 2 doublets 
        f = 6.167 / (centerlist[0] - centerlist[3]) + \
            1.677 / (centerlist[1] - centerlist[2])
        f = (f / 2) # mean f
        vmax = f * chan_mul
    elif n_center == 2:
        #folded spectrum with 1 doublet 
        f = 1.667 / (centerlist[0] - centerlist[1]) 
        vmax = f * chan_mul
    else:
        #limited to 6, 4 or 2 lines (3, 2, 1 doublet(s))
        print('The script can only handle folded spectra with 6, 4 or 2 peaks. Exit.')
        sys.exit(1)
    #return vmax and f (velocity per channel)
    return vmax, abs(f)

##########################################################################################
#plot results                                                                            #
##########################################################################################
def plot(ws5_raw_data, folded_intens, unfolded_spec, folded_spec, FP, v0, vmax, f, 
        filename, centerlist_FP, centerlist_v0, mod_date, mean_stdev_fold_i = 1):
    #plot the results
    #two subplots (unfolded & folded)
    fig, (ax0, ax1) = plt.subplots(2, 1)
    #title: filename + modification date of the .ws5 file
    fig.suptitle(filename + ' (' + str(mod_date) + ')')
    #all channel x values 1...512 for example for the upper plot
    x_raw = np.linspace(1,len(ws5_raw_data),len(ws5_raw_data))
    #title sub plot 0
    ax0.set_title('raw spectrum')
    #raw data points
    ax0.plot(x_raw, 
             ws5_raw_data, 
             '.', 
             color = 'steelblue', 
             label = 'raw data')
    #best fit curve
    ax0.plot(x_raw, 
             unfolded_spec.best_fit, 
             '-',
             color ='darkorange', 
             label ='best fit ' 
                    + '('+r'$R^2 =$ ' + '{:.4}'.format(unfolded_spec.rsquared)+')')
    #residuals        
    ax0.plot(x_raw, 
             1 - unfolded_spec.residual + max(ws5_raw_data)*0.01 + max(ws5_raw_data), 
             '.', 
             color = 'lightgrey', 
             label = 'residuals')
             
    #vertical line for FP (folding point) 
    ax0.axvline(x = FP.nominal_value, 
                ls = '--', 
                color = 'black')
    
    #label the FP line
    ax0.annotate('$FP = {:.4f}$'.format(FP.nominal_value),
                 (FP.nominal_value,min(ws5_raw_data)),
                 textcoords = "offset points",
                 xytext=(5,1), 
                 rotation=90)
                 
    #label the centers of the Lorentz functions (used for calculation of FP)
    for index,center in enumerate(centerlist_FP):
        heightkey = 's'+ str(index) + '_height'
        height = max(ws5_raw_data) + unfolded_spec.params[heightkey]
        ax0.annotate('{:.2f}'.format(center.nominal_value),
                     (center.nominal_value, height),
                     textcoords = "offset points",
                     xytext=(8,0), 
                     rotation=90, 
                     size=6)
    
    #label y axis, set ticks and show legend
    ax0.set_ylabel('intensity')
    ax0.set_xticks(np.linspace(1,len(ws5_raw_data),8))
    ax0.set_xlim([1, len(ws5_raw_data)])
    ax0.legend(fancybox = True, shadow = True, loc='upper right', prop={'size': 6})
    
    #R² is wrongly calculated from lmfit in case of weights <> 1
    r_squared = 1 - (folded_spec.residual * mean_stdev_fold_i).var() / \
                np.var(folded_spec.data)
    #all channel x values 1...256 for example for the lower plot
    x_fold = np.linspace(1,len(folded_intens),len(folded_intens))
    #title of sub plot 1
    ax1.set_title('folded spectrum')
    #folded raw data points
    ax1.plot(x_fold, 
             folded_intens, 
             '.',
             color = 'steelblue', 
             label = 'folded raw data')
    #best fit curve for folded data
    ax1.plot(x_fold, 
             folded_spec.best_fit, 
             '-', 
             color='darkorange', 
             label='best fit ' 
                   #+ '('+r'$R^2 =$ ' + '{:.4}'.format(folded_spec.rsquared)+')')
                   #R² is wrongly calculated from lmfit in case of weights <> 1
                   + '('+r'$R^2 =$ ' + '{:.4}'.format(r_squared)+')')
    #residuals for folded data
    #residuals multiplied with mean stddev
    ax1.plot(x_fold, 
             1 - folded_spec.residual * mean_stdev_fold_i + max(folded_intens)*0.01 +\
                 max(folded_intens), 
             '.', 
             color = 'lightgrey', 
             label = 'residuals')
    #vertical line for v0 (channel where the velocity is zero)
    ax1.axvline(x = v0.nominal_value, 
               ls = '--', 
               color = 'black')
    #label the v0 line
    ax1.annotate('$v_0 = {:.4f}$'.format(v0.nominal_value),
                 (v0.nominal_value,min(folded_intens)),
                 textcoords = "offset points",
                 xytext = (5,1), 
                 rotation = 90)
            
    #label the centers of the Lorentz functions (used for calculation of v0)
    for index,center in enumerate(centerlist_v0):
        heightkey = 's'+ str(index) + '_height'
        height = max(folded_intens) + folded_spec.params[heightkey]
        ax1.annotate('{:.2f}'.format(center.nominal_value),
                     (center.nominal_value, height),
                     textcoords = "offset points",
                     xytext = (8,0), 
                     rotation = 90, 
                     size = 6)
                     
    #label x axis           
    ax1.set_xlabel('channel no.')
    #label y axis, set ticks and show legend
    ax1.set_ylabel('intensity')
    ax1.set_xticks(np.linspace(1,len(folded_intens),16))
    ax1.set_xlim([1, len(folded_intens)])
    ax1.legend(fancybox = True, shadow = True, loc='upper right', prop={'size': 6})
    
    #print the github link to source and documentation
    ax1.annotate('https://github.com/radi0sus/cal-mb',
                xy = (1.01, 0.5),
                xycoords = 'axes fraction',
                ha = 'left',
                va = 'center',
                rotation = 270,
                fontsize = 8)
                
    #print the results of the fits, FP, v0, vmax and f
    ax1.annotate(u'$FP (\mathrm{{c}}) = {:.4fP}$  '.format(FP) 
               + u'$v_0 (\mathrm{{c}}) = {:.4fP}$  '.format(v0) 
               + u'$v_\mathrm{{max}} = {:.4fP}$ mm/s  '.format(vmax) 
               + u'$f = {:.4fP}$ mm/s per c'.format(f),
                xy = (-0.05, -0.2),
                xycoords = 'axes fraction',
                ha = 'left',
                va = 'center',
                rotation = 0,
                fontsize = 8)
                
    #secondary x-axis (secax) with velocity on top of ax1
    def chan_2_vel(x_fold, vmax = vmax.nominal_value, v0 = v0.nominal_value, \
                   N_chan = len(ws5_raw_data)):
        return vmax - (vmax + vmax)/(N_chan/2-1)*(x_fold+ (N_chan/2-1)/2 - v0)
    
    start, end = ax1.get_xlim()
    secax = ax1.secondary_xaxis('top', functions=(chan_2_vel,chan_2_vel))
    secax.set_ticks(np.arange(-end, end, 1))
    secax.set_xlabel('velocity $\mathrm{/mm \cdot s^{-1}}$')
    secax.tick_params(axis='x',direction='in')
    #end secondary x-axis (secax) with velocity on top of ax1           
   
    #arrange the plot window and show the plot
    mng = plt.get_current_fig_manager()
    mng.resize(1024,768)
    #(windows) low-res N = 1.2
    #high-res N = 1.5
    N = 1.5
    params = plt.gcf()
    plSize = params.get_size_inches()
    params.set_size_inches((plSize[0]*N*1, plSize[1]*N*1.5))
    plt.tight_layout()
    plt.show()

##########################################################################################
#argument parser                                                                         #
#parse arguments                                                                         #
##########################################################################################
parser = argparse.ArgumentParser(prog = 'cal-mb', 
                          description = 'Easily calibrate Mößbauer (MB) spectra')

#filename is required
parser.add_argument('filename', help = 'file with raw data from 57Fe foil (ws5)')

#show the matplotlib window
parser.add_argument('-s','--show',
    default=0, action='store_true',
    help='show the plot window')

#folding to the left
parser.add_argument('-fl','--foldleft',
    default=0, action='store_true',
    help='fold raw spectrum to the left')
    
#parse arguments
args = parser.parse_args()

##########################################################################################
#main                                                                                    #
##########################################################################################

#import data, filename and modification date of the .ws5 file
ws5_raw_data, filename, mod_date = op_imp(args.filename)

#fit unfolded spectrum
unfolded_spec = do_the_fit(ws5_raw_data)
#calculate FP (folding point)
FP, centerlist_FP = calc_FP_v0(unfolded_spec)
#fold spectrum
folded_intens, mean_stdev_fold_i = fold_spec(ws5_raw_data, FP)
#fit folded spectrum
folded_spec = do_the_fit(folded_intens, mean_stdev_fold_i)
#calculate v0 (channel where velocity is zero)
v0, centerlist_v0 = calc_FP_v0(folded_spec)
#calc vmax (maximum velocity) and f (velocity per channel)
vmax, f = calc_vmax(folded_spec, ws5_raw_data)

#print results, filename and modification date of the .ws5 file
print('======================================')
print('Results for', filename, ':')
print('File modified on', mod_date)
print('--------------------------------------')
print('FP (channel) =', u'{:.4fP}'.format(FP))
print('v₀ (channel) =', u'{:.4fP}'.format(v0))
print('vmax /mm·s⁻¹ = ', u'{: .4fP}'.format(vmax))
print('f /mm·s⁻¹/c  =  ', u'{:.4fP}'.format(f))
#print statistics from folded spectra
print('======================================')
print('Statistics (folded data with weights):')
print('--------------------------------------')
print('data points : ', '{}'.format(folded_spec.ndata))
print('variables   : ', '{}'.format(folded_spec.nvarys))
print('mean σ data : ', '{:.2f}'.format(mean_stdev_fold_i))
print('χ²          : ', '{:.2f}'.format(folded_spec.chisqr))
print('red. χ²     : ', '{:.2f}'.format(folded_spec.redchi))
 #R² is wrongly calculated from lmfit in case of weights <> 1
r_squared = 1 - (folded_spec.residual * mean_stdev_fold_i).var() / \
            np.var(folded_spec.data)
print('R²          : ', '{:.4f}'.format(r_squared))
print('======================================')

if args.show:
    #plot on request
    plot(ws5_raw_data, folded_intens, unfolded_spec, folded_spec, FP, v0, vmax, f, 
    filename, centerlist_FP, centerlist_v0, mod_date, mean_stdev_fold_i )
