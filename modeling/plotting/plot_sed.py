#! /usr/bin/env python
# USAGE:
# $python plot_sed.py [file_with_the_fluxes]
# file_with_the_fluxes should be like:
# survey	band	lambda(nm)	F_nu(Jy)	sigma(Jy)
# FUV	0.1528	0.15	0.0005720	0.0002480
# ...
# NOTE: First line is a header always. The columns are always the same which are separated by the \t delimeter. 
#
# Options:
# -skirt_files [skirt_sed_output1] ... [skirt_sed_outputN]  --- This will plot all SEDs given in the skirt output files with the fluxes (specified in SKIRT as *_sed.dat)
# -spread --- This option for plotting the spread. The [skirt_sed_output] should be re-written: lambda - first collumn, middle (best) model - second collumn, maximal emission - third column
#   minimal emission - fourth column
# -ima_name [name_of_the_pdf_file] --- Specify the output pdf file with the plotted SED
#
# EXAMPLES: python ~/CurrentWork/ImaPrep/IMAN/SKIRT/plot_sed.py fluxes.dat   - This will only plot SED points without the model
#           python ~/CurrentWork/ImaPrep/IMAN/SKIRT/plot_sed.py fluxes.dat -skirt_files sed.dat   - This will plot the SED points together with the model, given in sed.dat
#           python ~/CurrentWork/ImaPrep/IMAN/SKIRT/plot_sed.py fluxes.dat -skirt_files spread.dat model.dat -spread - This will create SEDs for the point, spread and the model

import sys
from os.path import exists
from os import remove
from scipy.ndimage.filters import convolve
from scipy import signal
import numpy as np
from numpy import *
from numpy import max, sum
import pyfits
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
import subprocess
import os
import glob
from matplotlib.pyplot import cm
import matplotlib.patches as mpatches

import matplotlib.gridspec as gridspec

from scipy.interpolate import interp1d

line_styles = ['-','--','-.',':']




def plotseds(data, model_data, plotfile, labels, figsize=(10,6), xlim=None, ylim=None):
    assert plotfile.endswith(".pdf")

    # Setup the figure
    figure = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 1,height_ratios=[4,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1],sharex=ax1)    

    # Set colours and markers for the flux points
    WAVELENGTH,FLUX,FLUX_ERR = data
    color=iter(cm.rainbow(np.linspace(0,1,len(WAVELENGTH))))  
    Labels = list(set(list(labels)))
    filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']
    Markers = filled_markers[:len(Labels)]
    

    used_labels = []
    for k in range(len(WAVELENGTH)):
      c=next(color)
      if labels[k] in Labels:
	MARKER = Markers[Labels.index(labels[k])]
	if labels[k] not in used_labels:
	  error_bar = np.array([[fabs(log10(FLUX[k])-log10(FLUX[k]-FLUX_ERR[k])) ,fabs(log10(FLUX[k])-log10(FLUX[k]+FLUX_ERR[k]))]]).T
	  used_labels.append(labels[k])
	  ax1.errorbar(WAVELENGTH[k], log10(FLUX[k]), yerr=error_bar,fmt=MARKER,markersize=7,color=c,markeredgecolor='black', ecolor=c, capthick=2)
	  ax1.plot(WAVELENGTH[k], log10(FLUX[k]),MARKER,markersize=7,color=c,markeredgecolor='black', markerfacecolor=c, label=labels[k])
	  ax2.errorbar(WAVELENGTH[k], 0., yerr=FLUX_ERR[k]/FLUX[k] * 100.,fmt=MARKER,markersize=7,color=c,markeredgecolor='black', ecolor=c, capthick=2)
	else:
	  error_bar = np.array([[fabs(log10(FLUX[k])-log10(FLUX[k]-FLUX_ERR[k])) ,fabs(log10(FLUX[k])-log10(FLUX[k]+FLUX_ERR[k]))]]).T
	  ax1.errorbar(WAVELENGTH[k], log10(FLUX[k]), yerr=error_bar,fmt=MARKER,markersize=7,color=c,markeredgecolor='black', ecolor=c, capthick=2)	  
	  ax2.errorbar(WAVELENGTH[k], 0., yerr=FLUX_ERR[k]/FLUX[k] * 100.,fmt=MARKER,markersize=7,color=c,markeredgecolor='black', ecolor=c, capthick=2)
	  

    if model_data[0]=='N':
      print 'No models were added!'
    else:
      for k in range(len(model_data)):
	wave_skirt, data_skirt, data_skirt_bottom, data_skirt_top, label_model = model_data[k]
	ax1.plot(wave_skirt,log10(data_skirt),line_styles[k],color='black',label=label_model)
	if len(data_skirt_bottom)>1 and len(data_skirt_top)>1 and k==0:
	  ax1.fill_between(wave_skirt, log10(data_skirt_bottom), log10(data_skirt_top), where=log10(data_skirt_top) <= log10(data_skirt_bottom), facecolor='cyan', edgecolor='cyan', interpolate=True, alpha=0.5)
	  ax1.plot([], [], color='cyan', linewidth=10, label='spread')
	  
	f2 = interp1d(wave_skirt,data_skirt, kind='cubic')     
	ax2.plot(WAVELENGTH,-(FLUX-f2(WAVELENGTH))/FLUX * 100.,line_styles[k],color='black',label='model')

      
    ax2.axhline(y=0.,color='black',ls='-.')
    ax2.set_ylim(-95,95)

    
    # Set axis limits if requested
    if xlim != None: ax1.set_xlim(xlim)
    if ylim != None: ax1.set_ylim(ylim)

    figure.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in figure.axes[:-1]], visible=False)
    
    # Add axis labels and a legend
    ax2.set_xscale('log')
    ax1.set_xscale('log')
    ax2.set_xlabel(r"Wavelength $\lambda\,[\mu \mathrm{m}]$", fontsize='large')
    ax1.set_ylabel(r"Log $F_\nu$$[Jy]$", fontsize='large')
    ax2.set_ylabel(r"Residuals $[\%]$", fontsize='large')
    
    # Add the legend
    ax1.legend(numpoints=1,loc=4,frameon=True,ncol=2,fontsize=11)

    # Save the figure
    plt.savefig(plotfile, bbox_inches='tight', pad_inches=0.25)
    plt.show()
    plt.close()
    return True






file_with_data = str(sys.argv[1])
if '-skirt_files' in sys.argv:
  SKIRT_DATA = []
  for k in range(sys.argv.index('-skirt_files')+1,len(sys.argv)):
    if sys.argv[k]!='-ima_name' and sys.argv[k]!='-spread':
      if '-spread' in sys.argv:
	try:
	  WAVELENGTH_SKIRT,FLUX_SKIRT,FLUX_SKIRT_BOTTOM,FLUX_SKIRT_TOP = loadtxt(sys.argv[k], usecols=[0,1,2,3],dtype=float, unpack=True,skiprows=7,delimiter=' ')
	  SKIRT_DATA.append([WAVELENGTH_SKIRT,FLUX_SKIRT,FLUX_SKIRT_BOTTOM,FLUX_SKIRT_TOP,str(sys.argv[k]).split('/')[-1].split('.')[0]])
	except:
	  WAVELENGTH_SKIRT,FLUX_SKIRT,FLUX_SKIRT_BOTTOM,FLUX_SKIRT_TOP = loadtxt(sys.argv[k], usecols=[0,1,2,3],dtype=float, unpack=True,skiprows=7,delimiter='\t')
	  SKIRT_DATA.append([WAVELENGTH_SKIRT,FLUX_SKIRT,FLUX_SKIRT_BOTTOM,FLUX_SKIRT_TOP,str(sys.argv[k]).split('/')[-1].split('.')[0]]) 
      else:
	WAVELENGTH_SKIRT,FLUX_SKIRT = loadtxt(sys.argv[k], usecols=[0,1],dtype=float, unpack=True,skiprows=7,delimiter=' ')
	SKIRT_DATA.append([WAVELENGTH_SKIRT,FLUX_SKIRT,[0.],[0.],str(sys.argv[k]).split('/')[-1].split('.')[0]])
    else:
      break
else:
  SKIRT_DATA = 'NONE'
  


  
if '-ima_name' in sys.argv:
  plotfile = str(sys.argv[sys.argv.index('-ima_name')+1])
else:
  plotfile = 'sed.pdf'
  
  
# survey	band	lambda(nm)	F_nu(Jy)	sigma(Jy)
SURVEY,BAND,WAVELENGTH,FLUX,FLUX_ERR = loadtxt(file_with_data, usecols=[0,1,2,3,4],dtype=str, unpack=True,skiprows=1,delimiter='\t')
WAVELENGTH = np.array(WAVELENGTH,float)
FLUX = np.array(FLUX,float)
FLUX_ERR = np.array(FLUX_ERR,float)

LABELS = SURVEY



plotseds([WAVELENGTH,FLUX,FLUX_ERR],SKIRT_DATA, plotfile, LABELS, figsize=(10,6), xlim=None, ylim=(-5,2))
print 'Done!'