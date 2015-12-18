#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.fitskirt.psf
# IMPORTANT: should be run after fitskirt_skysub
#
# This routine reads in a files list and a region list and determines the FWHM.
# The region list should be circles surrounding the stars used to determine the FWHM.
# Apart from the main function, functions have been copied from 
# Sebastien Viaene's skysub.py script.
#
# REQUIREMENTS: 
# - the files list should be free of ".fits" extension
# - the regions file should be made in ds9 and use circles saved as "image" or soatng
#
# OPTIONAL:
# - boolean for plotting the fitted result or not
#
# EXAMPLE: >pts fitskirt_psf files.dat stars.reg False
#
# This determines the FWHM for every star in stars.reg for all
# the images in files.dat and returns the average value.
#
# -----------------------------------------------------------------

# Import the necessary modules
import pyfits
import numpy as np
import itertools
import matplotlib.pyplot as plt
import sys
from scipy import optimize
from scipy.ndimage import interpolation
import matplotlib.pyplot as plt
from scipy.stats import nanmedian
from pylab import *

# -----------------------------------------------------------------

OK = '\033[92m'
ENDC = '\033[0m'

# Returns a 2D Gaussian distribution form the given parameters.
def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

#  Returns (height, x, y, width_x, width_y) the gaussian parameters of a 2D distribution by calculating its moments
def moments(data):
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

# Returns (height, x, y, width_x, width_y) the gaussian parameters of a 2D distribution found by a fit
def fitgaussian(data):
    params = moments(data)
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) - data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

# read the stars from the ds9 region and compute their central peak position through Gaussian fitting
def getStarPositions(frame,starRegion,plotting):
    
    # Read the ds9 region of stars
    reg_stars = open(starRegion,'r')
    lines = reg_stars.readlines()
    line = ''
    xpos = 0
    ypos = 0
    rad  = 0
    stars = np.zeros((len(lines)-1,4))
    # Loop over the lines, extract the region parameters and fit the stellar profiles with a 2D Gaussion
    for l in range (1,len(lines)):
    	line = lines[l]
    	if "circ" in line:
			line = line.replace('(',' ')
			line = line.replace(')',' ')
			line = line.replace(',',' ')
			line = line.split(' ')
			xpos = round(float(line[1]))
			ypos = round(float(line[2]))
			rad  = round(float(line[3]))
			# cut out a square around the star
			square = np.zeros((2*rad,2*rad))
			square = frame[ypos-rad:ypos+rad,xpos-rad:xpos+rad]

			# fit a 2D Gaussian to the brightness distribution
			params = fitgaussian(square)
			fit = gaussian(*params)
			(height, x, y, width_x, width_y) = params

			# Plot the result
			if plotting:
				matshow(square, cmap=cm.CMRmap)
				contour(fit(*indices(square.shape)), cmap=cm.Blues)
				ax = gca()

				text(0.95, 0.05, """
					pxl max: %.1f
					mod max: %.1f
					x : %.1f
					y : %.1f
					width_x : %.1f
					width_y : %.1f""" %(square[round(x),round(y)],height, x, y, width_x, width_y),
						fontsize=16, horizontalalignment='right',color='white',
						verticalalignment='bottom', transform=ax.transAxes)
				show()

			# Add the fitted parameters to the list of stars.
			# NOTE: for some reason, x and y are interchanged by python.
			stars[l-1,0] = xpos-rad+y
			stars[l-1,1] = ypos-rad+x
			stars[l-1,2] = width_y
			stars[l-1,3] = width_x

    return stars

# Main function, this part actually runs when routine is called
def main():

	# Reading arguments from command
	files = str(sys.argv[1])
	regions = str(sys.argv[2])
	plot =  False
	if len(sys.argv) > 3:
		plot = str(sys.argv[3]) in ['true', 'True', '1', 't', 'y','yes']
		
	# Read in frames file
	f = open(files, 'r')
	inputnames = f.read().splitlines()
	nframes = len(inputnames)
	
	# Loop over all frames and determine the FWHM
	for i in range(0,nframes):
		inputname = inputnames[i]+"_sub.fits"
		hdulist = pyfits.open(inputname)
		primhdr = pyfits.getheader(inputname, 0)
		frame = hdulist[0].data
		refstars = getStarPositions(frame,regions,plot)
		aver_x, aver_y = 0.0, 0.0
		for ref in refstars:
			aver_x+=ref[2]
			aver_y+=ref[3]
		aver_x = aver_x/float(len(refstars))
		aver_y = aver_y/float(len(refstars))
		aver = (aver_x + aver_y)/2.0		
		print "\n "+OK+inputname+" has FWHM : "+str("{:5.2f}".format(aver))+ENDC+"  ("+\
			str("{:5.2f}".format(aver_x))+" in x-dir and "+str("{:5.2f}".format(aver_y))+" in y-dir )"

if __name__ == '__main__':
    main()