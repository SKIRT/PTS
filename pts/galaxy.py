#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from scipy import signal, ndimage, stats

# Import relevant PTS modules
from pts.log import Log

#--------------------------------------------------------------------

## The class Galaxy
#
class Galaxy(object):

    ## The constructor
    def __init__(self, name, log=None):

        # Name
        self._name = name

        # Create a logger or use a pre-existing logger if possible
        self._log = Log() if log is None else log

        # Images
        self._images = []

        # Orientation
        self.orientation = None

        # Position angle
        self._positionangle = None

        # Center
        self._center = None

        # Ellipticity
        self._ellipticity = None

    ## This function
    @property
    def positionangle(self):

        if self._positionangle is None:

            self._log.warning("The orientation has not been determined yet for this galaxy")

        return self._positionangle

    ## This function
    def addimage(self, image):

        # Add this image to the list of images
        self._images.append(image)

    ## This function
    def findorientation(self):

        positionangles = []
        xcenters = []
        ycenters = []
        ellipticities = []

        # For every image of this galaxy
        for image in self._images:

            # Find the parameters from this image
            finder = GalaxyFinder(image.primary.data[::-1,:], quiet=True)

            # Add the results to the appropriate lists
            positionangles.append(finder.theta)
            xcenters.append(finder.xpeak)
            ycenters.append(finder.ypeak)
            ellipticities.append(finder.eps)

        # Set the properties of this galaxy
        self._positionangle = np.mean(positionangles)
        self._center = (np.mean(xcenters), np.mean(ycenters))
        self._ellipticity = np.mean(ellipticities)


#####################################################################
#
# Copyright (C) 1999-2014, Michele Cappellari
# E-mail: cappellari_at_astro.ox.ac.uk
#
# Updated versions of the software are available from my web page
# http://purl.org/cappellari/software
#
# This software is provided as is without any warranty whatsoever.
# Permission to use, for non-commercial purposes is granted.
# Permission to modify for personal or internal use is granted,
# provided this copyright and disclaimer are included unchanged
# at the beginning of the file. All other rights are reserved.
#
#####################################################################

## The class GalaxyFinder
#
class GalaxyFinder(object):

    ## The constructor of the GalaxyFinder class
    #
    #  With nblob=1 find the ellipse of inertia of the largest
    #  connected region in the image, with nblob=2 find the second
    #  in size and so on...
    #
    def __init__(self, img, fraction=0.1, quiet=False, nblob=1, level=None, log=None):

        # Copy the image
        self.img = img

        # Create a logger or use a pre-existing logger if possible
        self._log = Log() if log is None else log

        if len(img.shape) != 2:
            raise ValueError('IMG must be a two-dimensional array')

        a = signal.medfilt(img, 5)

        if level is None:
            level = stats.scoreatpercentile(a, (1 - fraction)*100)

        self.mask = a > level
        labels, nb = ndimage.label(self.mask)   # Get blob indices
        sizes = ndimage.sum(self.mask, labels, np.arange(nb + 1))
        j = np.argsort(sizes)[-nblob]      # find the nblob-th largest blob
        ind = np.flatnonzero(labels == j)

        self.second_moments(img, ind)

        if not quiet:

            self._log.info("Pixels used: " + str(ind.size))
            self._log.info("Peak (x,y): " + str(self.xpeak) + ", " + str(self.ypeak))
            self._log.info("Mean (x,y): " + str(self.xmed)  + ", " + str(self.ymed))
            self._log.info("Theta (deg): " + str(self.theta))
            self._log.info("Eps: " + str(self.eps))
            self._log.info("Sigma along major axis (pix): " + str(self.majoraxis))

    #-------------------------------------------------------------------------

    ## This function
    def plot(self):

        ax = plt.gca()
        ax.imshow(np.log(self.img.clip(self.img[self.xpeak, self.ypeak]/1e4)),
                      cmap='hot', origin='lower', interpolation='nearest')
        ax.imshow(self.mask, cmap='binary', interpolation='nearest',
                      origin='lower', alpha=0.3)
        ax.autoscale(False) # prevents further scaling after imshow()
        mjr = 3.5*self.majoraxis
        yc, xc = self.xmed, self.ymed
        ellipse = patches.Ellipse(xy=(xc, yc), width=2*mjr, fill=False,
                                      height=2*mjr*(1-self.eps), angle=-self.theta,
                                      color='red', linewidth=3)
        ax.add_artist(ellipse)
        ang = np.array([0,np.pi]) - np.radians(self.theta)
        ax.plot(xc - mjr*np.sin(ang), yc + mjr*np.cos(ang), 'g--',
                    xc + mjr*np.cos(ang), yc + mjr*np.sin(ang), 'g-', linewidth=3)
        ax.set_xlabel("pixels")
        ax.set_ylabel("pixels")
        plt.show()

    #-------------------------------------------------------------------------

    ## This function ...
    def second_moments(self, img, ind):

        # Restrict the computation of the first and second moments to
        # the region containing the galaxy, defined by vector IND.
        img1 = img.flat[ind]
        s = img.shape
        x, y = np.unravel_index(ind, s)

        # Compute coefficients of the moment of inertia tensor.
        i = np.sum(img1)
        self.xmed = np.sum(img1*x)/i
        self.ymed = np.sum(img1*y)/i
        x2 = np.sum(img1*x**2)/i - self.xmed**2
        y2 = np.sum(img1*y**2)/i - self.ymed**2
        xy = np.sum(img1*x*y)/i - self.xmed*self.ymed

        # Diagonalize the moment of inertia tensor.
        # theta is the angle, measured counter-clockwise,
        # from the image Y axis to the galaxy major axis.
        self.theta = np.degrees(np.arctan2(2*xy, x2 - y2)/2.) + 90.
        a2 = (x2 + y2)/2. + np.sqrt(((x2 - y2)/2.)**2 + xy**2)
        b2 = (x2 + y2)/2. - np.sqrt(((x2 - y2)/2.)**2 + xy**2)
        self.eps = 1. - np.sqrt(b2/a2)
        self.majoraxis = np.sqrt(a2)

        # If the image has many pixels then compute the coordinates of the
        # highest pixel value inside a 40x40 pixels region centered on the
        # first intensity moments (Xmed,Ymed), otherwise just return the
        # coordinates of the highest pixel value in the whole image.
        n = 20
        xmed1 = int(round(self.xmed))
        ymed1 = int(round(self.ymed))   # Check if subimage fits...
        if n <= xmed1 <= s[0]-n and n <= ymed1 <= s[1]-n:
            img2 = img[xmed1-n:xmed1+n, ymed1-n:ymed1+n]
            ij = np.unravel_index(np.argmax(img2), img2.shape)
            self.xpeak, self.ypeak = ij + np.array([xmed1, ymed1]) - n
        else:            # ...otherwise use full image
            self.xpeak, self.ypeak = np.unravel_index(np.argmax(img), s)

#-------------------------------------------------------------------------
