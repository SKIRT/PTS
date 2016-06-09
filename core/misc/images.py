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

# Import astronomical modules
from astropy.units import Unit
from astropy import constants

# Import the relevant PTS classes and modules
from ..tools.logging import log
from ..tools import filesystem as fs
from ..basics.filter import Filter
from ...magic.core.image import Image
from ...magic.core.frame import Frame
from ...magic.basics.coordinatesystem import CoordinateSystem
from ..tools.special import remote_filter_convolution, remote_convolution_frame

# -----------------------------------------------------------------

# The speed of light
speed_of_light = constants.c

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

        # The simulation prefix
        self.simulation_prefix = None

        # The paths to the 'total' FITS files produced by SKIRT
        self.fits_paths = None

        # The wavelengths of the simulation
        self.wavelengths = None

        # Filter names
        self.filter_names = ["FUV", "NUV", "u", "g", "r", "i", "z", "H", "J", "Ks", "I1", "I2", "I3", "I4", "W1", "W2",
                             "W3", "W4", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250", "SPIRE 350", "SPIRE 500"]

        # The instrument names
        self.instrument_names = None

        # The filters for which the images should be created
        self.filters = dict()

        # The dictionary containing the images for various SKIRT output datacubes
        self.images = dict()

        # The reference WCS
        self.wcs = None

    # -----------------------------------------------------------------

    def run(self, simulation, output_path=None, filter_names=None, instrument_names=None, wcs_path=None, kernel_paths=None, unit=None, host_id=None):

        """
        This function ...
        :param simulation:
        :param output_path:
        :param filter_names:
        :param instrument_names:
        :param wcs_path:
        :param kernel_paths:
        :param unit:
        :param host_id:
        :return:
        """

        # Obtain the paths to the 'total' FITS files created by the simulation
        self.fits_paths = simulation.totalfitspaths()

        # Get the list of wavelengths for the simulation
        self.wavelengths = simulation.wavelengths()

        # Get the simulation prefix
        self.simulation_prefix = simulation.prefix()

        # Set the filter names
        if filter_names is not None: self.filter_names = filter_names

        # Set the instrument names
        self.instrument_names = instrument_names

        # Create the filters
        self.create_filters()

        # Make the observed images
        self.make_images(host_id)

        # Set the WCS of the created images
        if wcs_path is not None: self.set_wcs(wcs_path)

        # Convolve the image with a given convolution kernel
        if kernel_paths is not None:

            # Check whether the WCS for the image is defined. If not, show a warning and skip the convolution
            if wcs_path is None: log.warning("WCS of the image is not defined, so convolution cannot be performed (the pixelscale is undefined)")
            else: self.convolve(kernel_paths, host_id)

        # Convert the units (WCS has to be loaded!)
        if unit is not None: self.convert_units(unit)

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

        # Loop over the different filter names
        for filter_name in self.filter_names:

            # Debugging
            log.debug("Constructing the " + filter_name + " filter ...")

            # Create the filter
            fltr = Filter.from_string(filter_name)

            # Add the filter to the list
            self.filters[filter_name] = fltr

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

            # Get the name of the instrument
            instr_name = instrument_name(path, self.simulation_prefix)

            # If a list of instruments is defined an this instrument is not in this list, skip it
            if self.instrument_names is not None and instr_name not in self.instrument_names: continue

            # Get the name of the datacube (as given by SKIRT)
            datacube_name = fs.strip_extension(fs.name(path))

            # Debugging
            log.debug("Making the observed images for " + datacube_name + ".fits ...")

            # Create a dictionary to contain the observed images for this FITS file
            images = dict()

            # The filter convolution is performed remotely
            if host_id is not None:

                # Upload the datacube, wavelength grid and filter properties, perform the convolution on the remote and get the resulting image frames back (as a dictionary where the keys are the filter names)
                #frames = remote_filter_convolution(host_id, path, self.wavelengths, self.filters, keep_output=True)
                frames = remote_filter_convolution(host_id, path, self.wavelengths, self.filters)

                # Add the resulting image frames to the dictionary
                for filter_name in frames:
                    # Add the observed image to the dictionary
                    images[filter_name] = frames[filter_name]

            # The calculation is performed locally
            else:

                # Load the simulated image
                datacube = Image.from_file(path, always_call_first_primary=False)

                # Convert the frames from neutral surface brightness to wavelength surface brightness
                for l in range(self.wavelengths):

                    # Get the wavelength
                    wavelength = self.wavelengths[l]

                    # Determine the name of the frame in the datacube
                    frame_name = "frame" + str(l)

                    # Divide this frame by the wavelength in micron
                    datacube.frames[frame_name] /= wavelength

                    # Set the new unit
                    datacube.frames[frame_name].unit = "W / (m2 * arcsec2 * micron)"

                # Convert the datacube to a numpy array where wavelength is the third dimension
                fluxdensities = datacube.asarray()

                # Loop over the different filters
                for filter_name in self.filters:

                    fltr = self.filters[filter_name]

                    # Debugging
                    log.debug("Making the observed image for the " + str(fltr) + " filter ...")

                    # Calculate the observed image frame
                    data = fltr.convolve(self.wavelengths, fluxdensities)
                    frame = Frame(data)

                    # Set the unit of the frame
                    frame.unit = "W/(m2 * arcsec2 * micron)"

                    # Add the observed image to the dictionary
                    images[filter_name] = frame

            # Add the dictionary of images of the current datacube to the complete images dictionary (with the datacube name as a key)
            self.images[datacube_name] = images

    # -----------------------------------------------------------------

    def set_wcs(self, wcs_path):

        """
        This function ...
        :param wcs_path:
        :return:
        """

        # TODO: allow multiple paths (in a dictionary) for the different datacubes (so that for certain instruments the WCS should not be set on the simulated images)

        # Inform the user
        log.info("Setting the WCS of the simulated images ...")

        # Debugging
        log.debug("Loading the coordinate system from '" + wcs_path + "' ...")

        # Load the WCS
        self.wcs = CoordinateSystem.from_file(wcs_path)

        # Loop over the different images and set the WCS
        for datacube_name in self.images:
            for filter_name in self.images[datacube_name]:

                # Debugging
                log.debug("Setting the coordinate system of the " + filter_name + " image of the '" + datacube_name + "' instrument ...")

                # Set the coordinate system for this frame
                self.images[datacube_name][filter_name].wcs = self.wcs

    # -----------------------------------------------------------------

    def convolve(self, kernel_paths, host_id=None):

        """
        This function ...
        :param kernel_paths:
        :param host_id:
        :return:
        """

        # Inform the user
        log.info("Convolving the images ...")

        # If the convolutions must be performed remotely
        if host_id is not None:

            # Loop over the images
            for datacube_name in self.images:
                for filter_name in self.images[datacube_name]:

                    # Check if the name of the image filter is a key in the 'kernel_paths' dictionary. If not, don't convolve.
                    if filter_name not in kernel_paths or kernel_paths[filter_name] is None: continue

                    # Determine the kernel path for this image
                    kernel_path = kernel_paths[filter_name]

                    # Perform the remote convolution
                    self.images[datacube_name][filter_name] = remote_convolution_frame(self.images[datacube_name][filter_name], kernel_path, host_id)

        # The convolution is performed locally
        else:

            # Loop over the images
            for datacube_name in self.images:
                for filter_name in self.images[datacube_name]:

                    # Check if the name of the image filter is a key in the 'kernel_paths' dictionary. If not, don't convolve.
                    if filter_name not in kernel_paths or kernel_paths[filter_name] is None: continue

                    # Load the kernel
                    kernel = Frame.from_file(kernel_paths[filter_name])

                    # Debugging
                    log.debug("Convolving the '" + filter_name + "' image of the '" + datacube_name + "' instrument ...")

                    # Convolve this image frame
                    self.images[datacube_name][filter_name].convolve(kernel)

    # -----------------------------------------------------------------

    def convert_units(self, unit):

        """
        This function ...
        :param self:
        :param unit:
        :return:
        """

        # TODO: right now, this is just an implementation of the conversion from W / (m2 * arcsec2 * micron) to MJy/sr
        # 1 Jy = 1e-26 * W / (m2 * Hz)

        # Inform the user
        log.info("Converting the units of the images to " + str(unit) + " ...")

        # Get the pixelscale
        #pixelscale = self.wcs.xy_average_pixelscale.to("arcsec/pix").value # in arcsec**2 / pixel

        # Loop over the images
        for datacube_name in self.images:
            for filter_name in self.images[datacube_name]:

                # Debugging
                log.debug("Converting the unit of the " + filter_name + " image of the '" + datacube_name + "' instrument ...")

                # Get the pivot wavelength of the filter
                fltr = self.filters[filter_name]
                pivot = fltr.pivotwavelength() * Unit("micron")

                # Determine the conversion factor
                conversion_factor = 1.0

                # From surface brightness to flux density (no)
                #conversion_factor *=

                # From W / (m2 * arcsec2 * micron) to W / (m2 * arcsec2 * Hz)
                conversion_factor *= (pivot ** 2 / speed_of_light).to("micron/Hz").value

                # From W / (m2 * arcsec2 * Hz) to MJy / sr
                conversion_factor *= (Unit("W/(m2 * arcsec2 * Hz)") / Unit("MJy/sr")).to("")

                # Convert
                self.images[datacube_name][filter_name] *= conversion_factor
                self.images[datacube_name][filter_name].unit = "MJy/sr"

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
                path = fs.join(output_path, datacube_name + "__" + filter_name + ".fits")

                # Save the image
                self.images[datacube_name][filter_name].save(path)

# -----------------------------------------------------------------

def instrument_name(datacube_path, prefix):

    """
    This function ...
    :param datacube_path:
    :param prefix:
    :return:
    """

    return fs.name(datacube_path).split("_total.fits")[0].split(prefix + "_")[1]

# -----------------------------------------------------------------
