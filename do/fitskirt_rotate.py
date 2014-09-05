#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.fitskirt_rotate
# IMPORTANT: should be run after fitskirt_skysub (or fitskirt_psf)
#
# This routine reads in a files list and 2 coordinates sets and then rotates
# the frame so the two sets are placed horizontally around the center.
# The first set is the left edge of the galaxy while
# the second set is the right edge of the galaxy.
#
# REQUIREMENTS: 
# - the files list should be free of ".fits" extension
# - the x-coordinate of the left edge
# - the y-coordinate of the left edge
# - the x-coordinate of the right edge
# - the y-coordinate of the right edge
#
# OPTIONAL:
# - boolean asking to flip the final result 180 degrees
#
# EXAMPLE: >pts fitskirt_rotate files.dat 484 881 1009 635 True
#
# This rotates all frames in files.dat to align 484,881 and 1009,635
# around the center (746,758) and flips the final result vertically.
#
# -----------------------------------------------------------------
 
# Import the necessary modules
import pyfits
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt
from scipy import ndimage
from numpy import *
import sys

# -----------------------------------------------------------------

# Main function, this part actually runs when routine is called
def main():

	# Reading arguments from command
	files = str(sys.argv[1])
	left_x = int(sys.argv[2])
	left_y = int(sys.argv[3])
	right_x = int(sys.argv[4])
	right_y = int(sys.argv[5])

	# If given in command, change values for flipping the image
	flip = False
	if len(sys.argv) > 6:
		flip = str(sys.argv[6]) in ['true', 'True', '1', 't', 'y','yes']
	
	# Read in frames file
	f = open(files, 'r')
	inputnames = f.read().splitlines()
	nframes = len(inputnames)
	
	# Loop over all frames and subtract background
	for i in range(0,nframes):
		inputname = inputnames[i]+"_sub.fits"
		outputname = inputnames[i]+"_rot.fits"
		hdulist = pyfits.open(inputname)
		primhdr = pyfits.getheader(inputname, 0)
		inframe = hdulist[0].data
		where_are_NaNs = isnan(inframe)
		inframe[where_are_NaNs] = 0
		print inframe.shape[0], inframe.shape[1]
		shift_x = inframe.shape[0]/2 - (left_x + right_x)/2
		shift_y = inframe.shape[1]/2 - (left_y + right_y)/2
		print "Shifted frame to center by "+str(shift_x)+","+str(shift_y)
		shiftframe = ndimage.interpolation.shift(inframe,(shift_x, shift_y))
		rotangle = math.degrees(math.atan(float(left_y - right_y)/float(left_x - right_x)))
		rotangle+=180.0 if flip else 0.0
		rotframe = ndimage.interpolation.rotate(shiftframe,rotangle)
		print "Rotated frame over "+str(rotangle)+" degrees"
		hdu = pyfits.PrimaryHDU(rotframe,primhdr)
		hdu.writeto(outputname,clobber=True)			

if __name__ == '__main__':
    main()
