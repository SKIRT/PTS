#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.polyfitcomparison Create a comparison sheet for the results of a polychromatic FitSKIRT fit.
#
# This script creates a comparison sheet (in the form of a PNG image) for the results of a polychromatic FitSKIRT fit.
# The in/out filenames and other parameters are hardcoded in the script.
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import matplotlib.colors as col
import matplotlib.cm as cm

# Import the relevant PTS classes and modules
from pts.core.basics.rgbimage import RGBImage

# -----------------------------------------------------------------

# wavelengths and names of reference images used during fit; number of output simulations
residualsname="ranges_Residual_64_"
wavelengthsname="ranges_Best_64_"
refs=["B.fits", "V.fits", "R.fits", "I.fits"]
numSimulations = 1

# colormaps and scales used
frame_cmap = "afmhot"
frame_min = 0.0
frame_max = 0.0015

res=[]
waves=[]
for i in range(0,len(refs)):
    
    print(i)
    res.append(residualsname+str(i)+".fits")
    waves.append(wavelengthsname+str(i)+".fits")

print(res)
print(waves)

# Define new colormap for residuals
def discrete_cmap(N=8):
    # define individual colors as hex values
    cpool = [ '#000000', '#00EE00', '#0000EE', '#00EEEE', '#EE0000',
              '#FFFF00', '#EE00EE', '#FFFFFF']
    cmap_i8 = col.ListedColormap(cpool[0:N], 'i8')
    cm.register_cmap(cmap=cmap_i8)

# -----------------------------------------------------------------

# build a list of residual images
resImages = []
for res in res:
     im = RGBImage(res)
     im.setrange(0,1)
     discrete_cmap();
     im.applycmap("i8")
     resImages.append(im)

# make all residual images the same shape, and remember the final shape of a single image
for i in resImages:
     for j in resImages:
          i.enlargecanvas(j)
imageShape = resImages[0].shape

# stack the residual images on top of each other
fullRes = resImages[0]
for i in range(1,len(resImages)):
     fullRes.addbelow(resImages[i])

fullRes.saveto("residuals.png")

# -----------------------------------------------------------------

# build a list of reference images
refImages = []
for ref in refs:
     im = RGBImage(ref)
     im.setrange(frame_min,frame_max)
     im.applycmap(frame_cmap)
     refImages.append(im)

# make all reference images the same shape, and remember the final shape of a single image
for i in refImages:
     for j in refImages:
          i.enlargecanvas(j)
imageShape = refImages[0].shape

# stack the reference images on top of each other
fullRef = refImages[0]
for i in range(1,len(refImages)):
     fullRef.addbelow(refImages[i])

# create a comparison sheet (for the last simulation result)

# build a list of simulated images
simImages = []
for j in waves:
    im = RGBImage(j)
    im.setrange(frame_min,frame_max)
    im.applycmap(frame_cmap)
    simImages.append(im)

# make all simulated images the same shape as the largest reference image
for j in simImages:
    j.enlargecanvas(imageShape)

# stack the simulated images on top of each other, and scale the values to match the JPEG image
fullSim = simImages[0]
for j in range(1,len(simImages)):
    fullSim.addbelow(simImages[j])

# glue the simulated stack to the right of the reference stack and save
fullRef.addright(fullSim)
fullRef.saveto("compare.png")

#glue the residuals stack to the right of the reference stack and simulated stack and save
fullRef.addright(fullRes)
fullRef.saveto("compare_residuals.png")

# -----------------------------------------------------------------
