#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt

# Import Astromagic modules
from .skyobject import SkyObject
from ..tools import analysis
from .vector import Position, Extent
from ..tools import interpolation
from astropy.coordinates import Angle

# Import astronomical modules
import astropy.io.fits as pyfits
import aplpy
from astropy import units as u
from astroquery.ned import Ned
import astroquery.exceptions
from astropy.coordinates import Angle
import astropy.coordinates as coord
from astroquery.vizier import Vizier

# *****************************************************************

class Galaxy(SkyObject):

    """
    This class ...
    """

    def __init__(self, name):

        """
        The constructor ...
        :param ra:
        :param dec:
        :param name:
        :return:
        """

        # Get the properties of this galaxy from online catalogs
        #self._get_properties(name)

        # Obtain more information about this galaxy
        try:

            ned_result = Ned.query_object(name)
            ned_entry = ned_result[0]

            # Get the NGC number
            self.name = ned_entry["Object Name"]

            # Get the redshift
            self.redshift = ned_entry["Redshift"]

            # Get the type (G=galaxy, HII ...)
            self.type = ned_entry["Type"]

        except astroquery.exceptions.RemoteServiceError:

            # Set attributes
            self.name = name
            self.redshift = None
            self.type = None

        # Create a new Vizier object and set the row limit to -1 (unlimited)
        viz = Vizier(keywords=["galaxies", "optical"])
        viz.ROW_LIMIT = -1

        # Query Vizier and obtain the resulting table
        result = viz.query_object(name, catalog=["VII/237"])
        entry = result[0][0]

        # Get the right ascension and the declination
        position = coord.SkyCoord(ra=entry["_RAJ2000"], dec=entry["_DEJ2000"], unit=(u.deg, u.deg), frame='fk5')

        # Get the names given to this galaxy
        self.names = entry["ANames"].split() if entry["ANames"] else None

        # Get the size of the galaxy
        ratio = np.power(10.0, entry["logR25"]) if entry["logR25"] else None
        diameter = np.power(10.0, entry["logD25"])*0.1*u.arcmin if entry["logD25"] else None

        #print("  D25_diameter = ", diameter)

        radial_profiles_result = viz.query_object(name, catalog="J/ApJ/658/1006")

        if len(radial_profiles_result) > 0:

            radial_profiles_entry = radial_profiles_result[0][0]

            self.distance = radial_profiles_entry["Dist"] * u.Unit("Mpc")
            self.inclination = Angle(radial_profiles_entry["i"], u.deg)
            self.d25 = radial_profiles_entry["D25"] * u.arcmin

        else:

            self.distance = None
            self.inclination = None
            self.d25 = None

        # Get the size of major and minor axes
        self.major = diameter
        self.minor = diameter / ratio if diameter is not None and ratio is not None else None

        # Get the position angle of the galaxy
        self.pa = Angle(entry["PA"]-90, u.deg) if entry["PA"] else None

        # Set the principal and companion flags to False initially
        self.principal = False
        self.companion = False

        # Initialize a list for the names of companion galaxies
        self.companions = []
        self.parent = None

        # Call the constructor of the base class
        super(Galaxy, self).__init__(position)

    # *****************************************************************

    def contains(self, position):

        """
        This function ...
        :param star:
        :return:
        """

        # If the position does not lie within the cutout box of the galaxy's source, return False
        if not self.source.cutout.contains(position): return False

        # If it does, check whether the pixel position is masked by the mask of the galaxy's source
        return self.source.mask.masks(self.source.cutout.rel_position(position))

    # *****************************************************************

    def ellipse_parameters(self, wcs, pixelscale, default_radius):

        """
        This function ...
        :param default_radius:
        :return:
        """

        if self.pa is None: angle = Angle(0.0, u.deg)
        else: angle = self.pa

        if self.major is None:

            x_radius = default_radius
            y_radius = default_radius

        elif self.minor is None or angle == 0.0:

            x_radius = 0.5 * self.major.to("arcsec") / pixelscale
            y_radius = x_radius

        else:

            x_radius = 0.5 * self.major.to("arcsec") / pixelscale
            y_radius = 0.5 * self.minor.to("arcsec") / pixelscale

        # Return the parameters
        return self.pixel_position(wcs), Extent(x=x_radius, y=y_radius), angle

    # *****************************************************************

    def fit_model(self, config):

        """
        This function ...
        :param frame:
        :param config:
        :return:
        """

        pass

    # *****************************************************************

    def remove(self, frame, config):

        """
        This function ...
        :param frame:
        :param config:
        :return:
        """

        # If a segment was found that can be identified with a source
        if self.has_source or config.remove_if_undetected:

            # Estimate the background
            self.source.estimate_background(config.remove_method, config.sigma_clip)

            # Replace the frame with the estimated background
            self.source.background.replace(frame, where=self.source.mask)

    # *****************************************************************

    def plot(self, frame):

        """
        This function ...
        :return:
        """

        if self.has_source:

            # Create a HDU from this frame with the image header
            hdu = pyfits.PrimaryHDU(self.source.cutout)

            # Create a figure canvas
            figure = plt.figure(figsize=(15, 15))

            # Create a figure from this frame
            plot = aplpy.FITSFigure(hdu, figure=figure)

            # Plot in color scale
            plot.show_colorscale()

            # Add a color bar if requested
            plot.add_colorbar()

            if self.has_aperture: self.aperture.plot(color='white', lw=1.5, alpha=0.5, ax=plt.gca())

            # Show the plot
            plt.show()

        else:

            # Create a HDU from this frame with the image header
            hdu = pyfits.PrimaryHDU(frame, frame.wcs.to_header())

            # Create a figure canvas
            figure = plt.figure(figsize=(20, 20))

            # Create a figure from this frame
            plot = aplpy.FITSFigure(hdu, figure=figure)

            # Plot in color scale
            plot.show_colorscale()

            # Add a color bar if requested
            plot.add_colorbar()

            if self.has_aperture: self.aperture.plot(color='white', lw=1.5, alpha=0.5, ax=plt.gca())

            # Show the plot
            plt.show()

# *****************************************************************