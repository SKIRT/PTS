#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np

# Import Astromagic modules
from ..tools import analysis
from .vector import Position, Extent
from ..tools import interpolation
from astropy.coordinates import Angle

# Import astronomical modules
from astropy import units as u
from photutils import segment_properties, properties_table
from photutils import EllipticalAperture
from astroquery.ned import Ned
import astroquery.exceptions
from astropy.coordinates import Angle
import astropy.coordinates as coord
from astroquery.vizier import Vizier

# *****************************************************************

class Galaxy(object):

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

        print(self.name)
        print("  type = ", self.type)

        # Create a new Vizier object and set the row limit to -1 (unlimited)
        viz = Vizier(keywords=["galaxies", "optical"])
        viz.ROW_LIMIT = -1

        # Query Vizier and obtain the resulting table
        result = viz.query_object(name, catalog=["VII/237"])
        entry = result[0][0]

        # Get the right ascension and the declination
        self.center = coord.SkyCoord(ra=entry["_RAJ2000"], dec=entry["_DEJ2000"], unit=(u.deg, u.deg), frame='fk5')

        # Get the names given to this galaxy
        self.names = entry["ANames"].split() if entry["ANames"] else None

        # Get the size of the galaxy
        ratio = np.power(10.0, entry["logR25"]) if entry["logR25"] else None
        diameter = np.power(10.0, entry["logD25"])*0.1*u.arcmin if entry["logD25"] else None

        print("  D25_diameter = ", diameter)

        radial_profiles_result = viz.query_object(name, catalog="J/ApJ/658/1006")

        if len(radial_profiles_result) > 0:

            radial_profiles_entry = radial_profiles_result[0][0]

            distance = radial_profiles_entry["Dist"] * u.Unit("Mpc")
            inclination = Angle(radial_profiles_entry["i"], u.deg)
            d25 = radial_profiles_entry["D25"] * u.arcmin

            print("  Distance = ", distance)
            print("  Inclination = ", inclination)
            print("  D25 = ", d25)

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

        # Set the source attribute to None initially
        self.source = None

        # Set the aperture attribute to None initially
        self.aperture = None

    # *****************************************************************

    @property
    def has_source(self):

        """
        This function ...
        :return:
        """

        return self.source is not None

    # *****************************************************************

    @property
    def has_aperture(self):

        """
        This function ...
        :return:
        """

        return self.aperture is not None

    # *****************************************************************

    def ellipse_parameters(self, wcs, pixelscale, default_radius):

        """
        This function ...
        :param default_radius:
        :return:
        """

        # Get the center of the galaxy in pixel coordinates
        x_center, y_center = self.center.to_pixel(wcs, origin=0)

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
        return Position(x=x_center, y=y_center), Extent(x=x_radius, y=y_radius), angle

    # *****************************************************************

    def find_source(self, frame, config):

        """
        This function ...
        :return:
        """

        # Get the parameters of the ellipse
        center, radius, angle = self.ellipse_parameters(frame.wcs, frame.pixelscale, config.initial_radius)

        # Find a source
        self.source = analysis.find_source(frame, center, radius, angle, config)

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
        if self.has_source:

            # Estimate the background
            self.source.estimate_background(config.remove_method, config.sigma_clip)

            #from ..tools import plotting
            #plotting.plot_box(np.ma.masked_array(self.source.background, mask=self.source.mask))

            # Replace the frame with the estimated background
            self.source.background.replace(frame, where=self.source.mask)

    # *****************************************************************

    def find_aperture(self, sigma_level=3.0):

        """
        This function ...
        :return:
        """

        props = segment_properties(self.source.cutout, self.source.mask)
        #tbl = properties_table(props)

        x_shift = self.source.cutout.x_min
        y_shift = self.source.cutout.y_min

        # Since there is only one segment in the self.source.mask (the center segment), the props
        # list contains only one entry (one galaxy)
        properties = props[0]

        # Obtain the position, orientation and extent
        position = (properties.xcentroid.value + x_shift, properties.ycentroid.value + y_shift)
        a = properties.semimajor_axis_sigma.value * sigma_level
        b = properties.semiminor_axis_sigma.value * sigma_level
        theta = properties.orientation.value

        # Create the aperture
        self.aperture = EllipticalAperture(position, a, b, theta=theta)

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