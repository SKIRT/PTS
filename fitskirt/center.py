#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.fitskirt_center
# IMPORTANT: should be run after fitskirt_rotate
#
# This routine reads in a files list and one coordinate set and two lower bounds.
# It cuts out the frames respecting the bounds and placing the coordinate in the center.
# Additionally, a rebinningfactor can be chosen.
#
# REQUIREMENTS: 
# - the files list should be free of ".fits" extension
# - the x-coordinate of the center
# - the y-coordinate of the center
# - the x-coordinate of the lower bound on the x-coordinate (left) 
# - the y-coordinate of the lower bound on the y-coordinate (down) 
#
# OPTIONAL:
# - integer n used to rebin the frame by combining (n,n) to (1,1)
#
# EXAMPLE: >pts fitskirt_center files.dat 1007 992 612 920 2
#
# This centers the frame from 612 to 1402 in x (so 1007 is the center)
# and 920 to 1064 in y (so 992 is the center) and rebins by 2,2
#
# -----------------------------------------------------------------

# Import the necessary modules
import pyfits
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt
from scipy import ndimage
import sys

# Main function, this part actually runs when routine is called
def main():

	# Reading arguments from command
	files = str(sys.argv[1])
	center_x = int(sys.argv[2])
	center_y = int(sys.argv[3])
	low_x = int(sys.argv[4])
	low_y = int(sys.argv[5])
	rebinf = int(sys.argv[6]) if len(sys.argv) > 6 else 1
	print "rebinningfactor is = "+str(rebinf)
	up_x = (center_x-low_x)*2+low_x
	up_y = (center_y-low_y)*2+low_y
	print up_x, up_y
	
	# Read in frames file
	f = open(files, 'r')
	inputnames = f.read().splitlines()
	nframes = len(inputnames)
	
	# Loop over all frames and subtract background
	for i in range(0,nframes):
		inputname = inputnames[i]+"_rot.fits"
		outputname = inputnames[i]+"_ext.fits"
		hdulist = pyfits.open(inputname)
		primhdr = pyfits.getheader(inputname, 0)
		inframe = hdulist[0].data
		mod_y = divmod(up_y-low_y, rebinf)
		mod_x = divmod(up_x-low_x, rebinf)
		extract = inframe[low_y:low_y+mod_y[0]*rebinf,low_x:low_x+mod_x[0]*rebinf]
		print "Extracted frame from x: "+str(low_x)+"-"+str(up_x)+" and y: "+str(low_y)+"-"+str(up_y)
		rebin = extract.reshape(mod_y[0], rebinf, mod_x[0], rebinf).mean(axis=3).mean(axis=1)
		hdu = pyfits.PrimaryHDU(rebin,primhdr)
		hdu.writeto(outputname,clobber=True)				
				
if __name__ == '__main__':
    main()
