#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.misc.fluxes Contains the ObservedImageMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools.logging import log
from ..tools import filesystem
from ..basics.filter import Filter
from ...magic.core.image import Image
from ...magic.core.frame import Frame
from ..tools.special import remote_filter_convolution

# -----------------------------------------------------------------

class ObservedImageMaker(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(ObservedImageMaker, self).__init__()

        # -- Attributes --

        # The paths to the FITS files produced by SKIRT
        self.fits_paths = None

        # The wavelengths of the simulation
        self.wavelengths = None

        # Filter names
        self.filter_names = ["FUV", "NUV", "u", "g", "r", "i", "z", "H", "J", "Ks", "I1", "I2", "I3", "I4", "W1", "W2",
                             "W3", "W4", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250", "SPIRE 350", "SPIRE 500"]

        # The filters for which the images should be created
        self.filters = None

        # The dictionary containing the images for various SKIRT output datacubes
        self.images = dict()

    # -----------------------------------------------------------------

    def run(self, simulation, output_path=None, filter_names=None, host_id=None):

        """
        This function ...
        :param simulation:
        :param output_path:
        :param filter_names:
        :param host_id:
        :return:
        """

        # Obtain the paths to the 'total' FITS files created by the simulation
        self.fits_paths = simulation.totalfitspaths()

        # Get the list of wavelengths for the simulation
        self.wavelengths = simulation.wavelengths()

        # Set the filter names
        if filter_names is not None: self.filter_names = filter_names

        # Create the filters
        self.create_filters()

        # Make the observed images
        self.make_images(host_id)

        # Write the results
        if output_path is not None: self.write(output_path)

    # -----------------------------------------------------------------

    def create_filters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Constructing the filter objects ...")

        # Initialize the list
        self.filters = []

        # Loop over the different filter names
        for filter_name in self.filter_names:

            # Debugging
            log.debug("Constructing the " + filter_name + " filter ...")

            # Create the filter
            fltr = Filter.from_string(filter_name)

            # Add the filter to the list
            self.filters.append(fltr)

    # -----------------------------------------------------------------

    def make_images(self, host_id=None):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Inform the user
        log.info("Making the observed images (this may take a while) ...")

        # Loop over the different simulated images
        for path in self.fits_paths:

            # Get the name of the datacube (as given by SKIRT)
            datacube_name = filesystem.strip_extension(filesystem.name(path))

            # Debugging
            log.debug("Making the observed images for " + datacube_name + ".fits ...")

            # Create a dictionary to contain the observed images for this FITS file
            images = dict()

            # The filter convolution is performed remotely
            if host_id is not None:

                # Upload the datacube, wavelength grid and filter properties, perform the convolution on the remote and get the resulting image frames back (as a dictionary where the keys are the filter names)
                frames = remote_filter_convolution(host_id, path, self.wavelengths, self.filters)

                # Add the resulting image frames to the dictionary
                for filter_name in frames:
                    # Add the observed image to the dictionary
                    images[filter_name] = frames[filter_name]

            # The calculation is performed locally
            else:

                # Load the simulated image
                datacube = Image.from_file(path)

                # Convert the datacube to a numpy array where wavelength is the third dimension
                fluxdensities = datacube.asarray()

                # densities must be per wavelength instead of per frequency!

                # Loop over the different filters
                for filter in self.filters:

                    # Debugging
                    log.debug("Making the observed image for the " + filter.name + " filter ...")

                    # Calculate the observed image frame
                    data = filter.convolve(self.wavelengths, fluxdensities)
                    frame = Frame(data)

                    # Add the observed image to the dictionary
                    images[filter.name] = frame

            # Add the dictionary of images of the current datacube to the complete images dictionary (with the datacube name as a key)
            self.images[datacube_name] = images

    # -----------------------------------------------------------------

    def write(self, output_path):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Inform the user
        log.info("Writing the images ...")

        # Loop over the different images (self.images is a nested dictionary of dictionaries)
        for datacube_name in self.images:
            for filter_name in self.images[datacube_name]:

                # Determine the path to the output FITS file
                path = filesystem.join(output_path, datacube_name + "__" + filter_name + ".fits")

                # Save the image
                self.images[datacube_name][filter_name].save(path)

# -----------------------------------------------------------------
