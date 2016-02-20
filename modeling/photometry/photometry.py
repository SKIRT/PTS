#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.photometry.photometry Contains the PhotoMeter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np

# Import astronomical modules
from astropy import units as u
from astropy.coordinates import Angle

# Import the relevant AstroMagic classes and modules
from ...magic.core import Image
from ...magic.tools import headers

# Import the relevant PTS classes and modules
from ..core import ModelingComponent
from ...core.tools import filesystem
from ...core.tools.logging import log
from ..core import ObservedSED

# -----------------------------------------------------------------

class PhotoMeter(ModelingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(PhotoMeter, self).__init__(config)

        # The list of images
        self.images = []

        # The SED
        self.sed = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new PhotoMeter instance
        photometer = cls(arguments.config)

        # Set the input and output path
        photometer.config.path = arguments.path
        photometer.config.input_path = os.path.join(arguments.path, "prep")
        photometer.config.output_path = os.path.join(arguments.path, "phot")

        # A single image can be specified so the photometry is only calculated for that image
        photometer.config.single_image = arguments.image

        # Return the new instance
        return photometer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the prepared images
        self.load_images()

        # 3. Do the photometry
        self.do_photometry()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PhotoMeter, self).setup()

        # Create an observed SED
        self.sed = ObservedSED()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # Loop over all directories in the preparation directory
        for directory_path, directory_name in filesystem.directories_in_path(self.prep_path, names=True):

            # If only a single image has to be processeds, skip the other images
            if self.config.single_image is not None and directory_name != self.config.single_image: continue

            # Determine the filter
            filter_name = directory_name
            if filter_name == "Ha": filter_name = "H alpha"
            filter = headers.get_filter(filter_name)

            # Look for a file called 'result.fits'
            image_path = os.path.join(directory_path, "result.fits")
            if not filesystem.is_file(image_path):

                log.warning("Prepared image could not be found for " + directory_name)
                continue

            # Open the prepared image
            image = Image.from_file(image_path)

            # Set the image name
            image.name = directory_name

            # Set the filter
            image.filter = filter

            # Add the image to the list
            self.images.append(image)

    # -----------------------------------------------------------------

    def do_photometry(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the photometry calculation ...")

        wavelength_column = []
        flux_column = []

        ## GET STRUCTURAL PARAMETERS

        from astroquery.vizier import Vizier
        vizier = Vizier(keywords=["galaxies"])

        # Get parameters from S4G catalog
        result = vizier.query_object(self.galaxy_name, catalog=["J/PASP/122/1397/s4g"])
        table = result[0]

        name = table["Name"][0]

        ra_center = table["_RAJ2000"][0]
        dec_center = table["_DEJ2000"][0]

        major = table["amaj"][0] * u.Unit("arcsec")
        ellipticity = table["ell"][0]
        position_angle = Angle(table["PA"][0] - 90.0, u.Unit("deg"))

        minor = (1.0-ellipticity)*major

        from ...magic.basics import Extent

        radius = Extent(major, minor)


        from ...magic.basics.skygeometry import SkyEllipse, SkyCoord

        center = SkyCoord(ra=ra_center, dec=dec_center, unit=(u.Unit("deg"), u.Unit("deg")), frame='fk5')

        sky_ellipse = SkyEllipse(center, radius, position_angle)

        from ...magic.basics import SkyRegion

        region = SkyRegion()
        region.append(sky_ellipse)
        region_path = self.full_output_path("galaxy.reg")
        region.save(region_path)

        ##

        for image in self.images:

            ellipse = sky_ellipse.to_ellipse(image.wcs)
            from ...magic.basics import Mask

            # Create mask
            mask = Mask.from_shape(ellipse, image.xsize, image.ysize)
            inverted_mask = mask.inverse()

            # Inform the user
            log.debug("Performing photometry for " + image.name + " image ...")

            wavelength_column.append(image.wavelength.to("micron").value)

            # Convert from MJy/sr to Janskys

            # Convert from MJy/sr to Jy/sr
            conversion_factor = 1.0

            # Conversion from Jy to MJy
            conversion_factor *= 1e6

            # Conversion from MJy (per pixel2) to MJy / sr
            pixelscale = image.frames.primary.xy_average_pixelscale
            pixel_factor = (1.0/pixelscale**2).to("pix2/sr").value
            conversion_factor /= pixel_factor

            # Frame in Janskys
            jansky_frame = image.frames.primary * conversion_factor
            jansky_frame[inverted_mask] = 0.0

            # Calculate the total flux in Jansky
            flux = np.sum(jansky_frame)

            jansky_errors_frame = image.frames.errors * conversion_factor
            jansky_errors_frame[inverted_mask] = 0.0
            flux_error = np.sum(jansky_errors_frame)

            # Add this entry to the SED
            self.sed.add_entry(image.filter, flux, flux_error)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write SED table
        self.write_sed()

        # Plot the SED
        self.plot_sed()

    # -----------------------------------------------------------------

    def write_sed(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the output file
        path = self.full_output_path("fluxes.dat")

        # Save the SED
        self.sed.save(path)

    # -----------------------------------------------------------------

    def plot_sed(self):

        """
        This function ...
        :return:
        """

        from ..plotting import SEDPlotter

        plotter = SEDPlotter()

        # Add the SED
        plotter.add_observed_sed(self.sed, "Observation")

        # Determine the full path to the plot file
        path = self.full_output_path("sed.pdf")
        plotter.run(path)

# -----------------------------------------------------------------
