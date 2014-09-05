#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.fitskirt_mask
# IMPORTANT: should be run after fitskirt_center
#
# this routine reads in a files list and a region list made in ds9
# the files are then all masked with a given value 
# REQUIREMENTS: 
# - the files list should be free of ".fits" extension
# - the region list with polygons or circles
# OPTIONAL:
# - a float value used to put in the mask, default set to 0.0
# - a boolean used to determine the normalisation of the final frames (True = normalise)
# EXAMPLE: 
# >pts fitskirt_mask files.dat mask.reg 0.0
# this masks all the files in files.dat and puts in the value 0.0
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
from matplotlib.path import Path
from itertools import product

# Main function, this part actually runs when routine is called
def main():
	
	# Reading arguments from command
	files = str(sys.argv[1])
	region = str(sys.argv[2])
	value =	float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
	norm = True
	if len(sys.argv) > 4:
		norm = str(sys.argv[4]) in ['true', 'True', '1', 't', 'y','yes']

	# Read in frames file
	f = open(files, 'r')
	inputnames = f.read().splitlines()
	nframes = len(inputnames)
	
	# Loop over all frames and read in mask regions
	for i in range(0,nframes):
		inputname = inputnames[i]+"_ext.fits"
		outputname = inputnames[i]+"_norm.fits"
		hdulist = pyfits.open(inputname)
		primhdr = pyfits.getheader(inputname, 0)
		inframe = hdulist[0].data
		reg = open(region,'r')
		paths=[]
		
		# Loop over all regions 
		for line in reg:
		
			# Detect circle regions and fill them with the mask value
			if 'circle(' in line:
				param = ((line.split('circle(')[1]).split(')')[0]).split(',')
				a, b ,r  = int(float(param[0])), int(float(param[1])), int(float(param[2])) 
				y,x = np.ogrid[-b:inframe.shape[0]-b, -a:inframe.shape[1]-a]
				mask = x*x + y*y <= r*r
				inframe[mask] = value
				
			# Detect polygon regions and add them to the path list
			elif 'polygon(' in line:
				param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
				param2 = [None]*(len(param)/2)
				for i in xrange(0,len(param2)):	param2[i] = (int(param[2*i]),int(param[2*i+1])) 
				param2.append(param2[0])
				codes = []
				codes.append(Path.MOVETO)
				for i in xrange(1,len(param2)-1): codes.append(Path.LINETO)
				codes.append(Path.CLOSEPOLY)
				path = Path(param2, codes)
				paths.append(path)
		
		# Loop over the image and all polygons and fill them with the mask value	
		nx, ny =inframe.shape[1], inframe.shape[0]
		for i, j, path in product(range(0,nx), range(0,ny), paths):
			inframe[j][i]=value if path.contains_point((i,j)) else inframe[j][i]
		
		# Normalise the frame and write it out
		total = inframe.sum()
		print "The total value of the frame is ", str(total)
		if norm:
			for i, j in product(range(0,nx), range(0,ny)):
				inframe[j][i]=inframe[j][i]/total
		hdu = pyfits.PrimaryHDU(inframe,primhdr)
		hdu.writeto(outputname,clobber=True)	

if __name__ == '__main__':
    main()
