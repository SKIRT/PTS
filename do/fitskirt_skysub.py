#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.fitskirt_skysub
# This routine reads in a files list, a region list and a polynome order
# the fourth and fifth argument are optional and used to change the
# default values for plotting and linear fitting.
# Apart from the main function, functions have been copied from 
# Sebastien Viaene's skysub.py script.
#
# REQUIREMENTS: 
# - the files list should be free of ".fits" extension
# - the regions file should be made in ds9 and use boxes saved as "image" or soatng
# - integers used as the order of the polynome fit
#
# OPTIONAL:
# - Boolean controlling wether to plot the result or not
#
# EXAMPLE: >pts fitskirt_skysub files.dat sky.reg 2 False
#
# This subtracts the sky from files in files.dat using the regions
# defined in ds9 and saved as soatng in sky.reg.
# The sky is fitted with a 2nd-order polynomial and the results is 
# not shown during the python routine
#
# -----------------------------------------------------------------

# Import the necessary modules
import pyfits
import numpy as np
import itertools
import matplotlib.pyplot as plt
import sys

# -----------------------------------------------------------------

# This function fits a polynomial to a set of background pixels
def subtractBackground(name, array,regfile,order,linear, plot):
  
    # Generate fitting data...
    xaxis = array.shape[1]
    yaxis = array.shape[0]
    x,y,z = getBgPixels(array,regfile)

    # Compyte mean background for plotting purposes
    meanbg=np.mean(z)
    
    # Fit a 3rd order, 2d polynomial
    m = polyfit2d(x,y,z,order,linear)
    
    # Evaluate it on a grid...
    nx, ny = 20, 20
    xx, yy = np.meshgrid(np.linspace(0, xaxis, nx),
                         np.linspace(0, yaxis, ny))
    zz = polyval2d(xx, yy, m)

    # Make a sky array and subtract it from the object frame
    sky = np.copy(array)
    strip = np.arange(1.0,np.float(xaxis+1.0))
    for i in range(0,yaxis):
        sky[i,:] = polyval2d(strip,i,m)
    
    if plot:
		# Plot polynomial
		fig = plt.figure(1)
		plt.subplot(121)
		plt.imshow(sky[::-1,:], vmin=0.95*meanbg, vmax=1.05*meanbg, cmap='hot', extent=(0,xaxis,0,yaxis))
		plt.colorbar(orientation='horizontal')
		plt.scatter(x, y,c='k',marker='+')
		plt.xlim(0,xaxis)
		plt.ylim(0,yaxis)
	
		# Plot image
		plt.subplot(122)
		plt.imshow(array[::-1,:], vmin=0.95*meanbg, vmax=1.05*meanbg, cmap='hot', extent=(0,xaxis,0,yaxis))
		plt.xlim(0,xaxis)
		plt.ylim(0,yaxis)
		plt.colorbar(orientation='horizontal')
		plt.show()
	
    hdu = pyfits.PrimaryHDU(sky)
    hdu.writeto('sky_'+name+'.fits',clobber=True)
    sframe = array - sky
    # Return the sky subtracted image
    return sframe


# Two-dimensional polynomial fit. 
# Based uppon code provided by Joe Kington
# http://stackoverflow.com/questions/7997152/
# python-3d-polynomial-surface-fit-order-dependent/7997925#7997925
def polyfit2d(x, y, z, order, linear):

    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
        if linear & (i != 0.) & (j != 0.):
            G[:, k] = 0
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m
    
# Values to two-dimensional polynomial fit. 
# Based uppon code provided by Joe Kington
def polyval2d(x, y, m):

    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z

# From a set of sky regions, a set of sky points is constructed
def getBgPixels(array, regfile):
    
    # Open and read the region file (boxes saved in SAOtng format)
    reg_stars = open(regfile,'r')
    lines = reg_stars.readlines()
    line = ''
    fluxes = []
    xco = []
    yco = []
    
    # Loop over the boxes
    for l in range (1,len(lines)):
        line = lines[l]
        line = line.replace('(',' ')
        line = line.replace(')',' ')
        line = line.replace(',',' ')
        line = line.split(' ')
        xco  = np.append(xco, round(float(line[1])))
        yco  = np.append(yco, round(float(line[2])))
        xmin = round(float(line[1])) - round(0.5*float(line[3]))
        xmax = round(float(line[1])) + round(0.5*float(line[3]))
        ymin = round(float(line[2])) - round(0.5*float(line[4]))
        ymax = round(float(line[2])) + round(0.5*float(line[4]))
        # append the median of each box to the set of sky points
        fluxes = np.append(fluxes,np.median(array[ymin:ymax,xmin:xmax]))
    
    # return the set of sky points
    return xco,yco,fluxes

# Main function, this part actually runs when routine is called
def main():

	# Reading arguments from command
	files = str(sys.argv[1])
	regions = str(sys.argv[2])
	order = int(sys.argv[3])

	# If given in command, change values for fitting and plotting
	linear, plot = True, False
	if len(sys.argv) > 4:
		plot = str(sys.argv[4]) in ['true', 'True', '1', 't', 'y','yes']
			
	# Read in frames file
	f = open(files, 'r')
	inputnames = f.read().splitlines()
	nframes = len(inputnames)
	
	# Loop over all frames and subtract background
	for i in range(0,nframes):
		inputname = inputnames[i]+".fits"
		outputname = inputnames[i]+"_sub.fits"
		hdulist = pyfits.open(inputname)
		primhdr = pyfits.getheader(inputname, 0)
		frame = hdulist[0].data
		sframe = subtractBackground(inputnames[i],frame,
							regions,order,linear, plot)
		hdu = pyfits.PrimaryHDU(sframe,primhdr)
		hdu.writeto(outputname,clobber=True)		

if __name__ == '__main__':
    main()