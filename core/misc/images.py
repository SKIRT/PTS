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
from astropy import constants

# Import the relevant PTS classes and modules
from ..tools.logging import log
from ..tools import filesystem as fs
from ..filter.filter import parse_filter
from ...magic.core.kernel import ConvolutionKernel
from ...magic.core.datacube import DataCube
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.core.remote import RemoteDataCube
from ..simulation.wavelengthgrid import WavelengthGrid
from ..basics.unit import parse_unit as u
from ..basics.configurable import Configurable

# -----------------------------------------------------------------

# The speed of light
speed_of_light = constants.c

# -----------------------------------------------------------------

class ObservedImageMaker(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(ObservedImageMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The simulation prefix
        self.simulation_prefix = None

        # The paths to the 'total' FITS files produced by SKIRT
        self.fits_paths = None

        # The wavelengths of the simulation
        self.wavelengths = None

        # The wavelength grid of the simulation
        self.wavelength_grid = None

        # Filter names
        self.filter_names = ["FUV", "NUV", "u", "g", "r", "i", "z", "H", "J", "Ks", "I1", "I2", "I3", "I4", "W1", "W2",
                             "W3", "W4", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250", "SPIRE 350", "SPIRE 500"]

        # The instrument names
        self.instrument_names = None

        # The filters for which the images should be created
        self.filters = dict()

        # The dictionary containing the different SKIRT output datacubes
        self.datacubes = dict()

        # The dictionary containing the created observation images
        self.images = dict()

        # The reference WCS
        self.wcs = None

        # The kernel paths
        self.kernel_paths = None

        # The target unit
        self.unit = None

        # The host id
        self.host_id = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the wavelength grid
        self.create_wavelength_grid()

        # 3. Create the filters
        self.create_filters()

        # 4. Load the datacubes
        self.load_datacubes()

        # 5. Set the WCS of the datacubes
        if self.wcs is not None: self.set_wcs()

        # 6. Make the observed images
        self.make_images()

        # 7. Do convolutions
        if self.kernel_paths is not None: self.convolve()

        # 8. Do unit conversions
        if self.unit is not None: self.convert_units()

        # 9. Write the results
        if self.output_path is not None: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ObservedImageMaker, self).setup(**kwargs)

        # simulation, output_path=None, filter_names=None, instrument_names=None, wcs_path=None,
        # kernel_paths=None, unit=None, host_id=None
        simulation = kwargs.pop("simulation")
        output_path = kwargs.pop("output_path", None)
        filter_names = kwargs.pop("filter_names", None)
        instrument_names = kwargs.pop("instrument_names", None)
        wcs_path = kwargs.pop("wcs_path", None)
        kernel_paths = kwargs.pop("kernel_paths", None)
        unit = kwargs.pop("unit", None)
        host_id = kwargs.pop("host_id", None)

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

        # Set the output path
        self.config.ouput = output_path

        # If a WCS path is defined (a FITS file)
        if wcs_path is not None:

            # Debugging
            log.debug("Loading the coordinate system from '" + wcs_path + "' ...")

            # Load the WCS
            self.wcs = CoordinateSystem.from_file(wcs_path)

        # Set the kernel paths
        self.kernel_paths = kernel_paths

        # Set the target unit
        self.unit = unit

        # Set the host id
        self.host_id = host_id

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Construct the wavelength grid from the array of wavelengths
        self.wavelength_grid = WavelengthGrid.from_wavelengths(self.wavelengths, "micron")

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
            fltr = parse_filter(filter_name)

            # Add the filter to the list
            self.filters[filter_name] = fltr

    # -----------------------------------------------------------------

    def load_datacubes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SKIRT output datacubes ...")

        # Loop over the different simulated images
        for path in self.fits_paths:

            # Get the name of the instrument
            instr_name = instrument_name(path, self.simulation_prefix)

            # If a list of instruments is defined an this instrument is not in this list, skip it
            if self.instrument_names is not None and instr_name not in self.instrument_names: continue

            # Get the name of the datacube (as given by SKIRT)
            datacube_name = fs.strip_extension(fs.name(path))

            # Debugging
            log.debug("Loading datacube from '" + datacube_name + ".fits' ...")

            ## LOAD AND CONVERT UNITS TO SPECTRAL (WAVELENGTH-DENSITY)

            # Load the datacube (locally or remotely)
            if self.host_id is not None: datacube = RemoteDataCube.from_file(path, self.wavelength_grid, self.host_id)
            else: datacube = DataCube.from_file(path, self.wavelength_grid)

            # Convert the datacube from neutral flux density to wavelength flux density
            datacube.to_wavelength_density("W / (m2 * arcsec2 * micron)", "micron")

            # Add the datacube to the dictionary
            self.datacubes[datacube_name] = datacube

    # -----------------------------------------------------------------

    def set_wcs(self):

        """
        This function ...
        :return:
        """

        # TODO: allow multiple paths (in a dictionary) for the different datacubes
        # (so that for certain instruments the WCS should not be set on the simulated images)

        # Inform the user
        log.info("Setting the WCS of the simulated images ...")

        # Loop over the different datacubes and set the WCS
        for datacube_name in self.datacubes:

            # Debugging
            log.debug("Setting the coordinate system of the " + datacube_name + "' datacube ...")

            # Set the coordinate system for this datacube
            self.datacubes[datacube_name].wcs = self.wcs

    # -----------------------------------------------------------------

    def make_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the observed images (this may take a while) ...")

        # Loop over the datacubes
        for datacube_name in self.datacubes:

            # Debugging
            log.debug("Making the observed images for " + datacube_name + ".fits ...")

            # Create a list of the filter names
            filter_names = self.filters.keys()

            # Create the corresponding list of filters
            filters = self.filters.values()

            # Initialize a dictionary, indexed by the filter names, to contain the images
            images = dict()

            # Create the observed images from the current datacube (the frames get the correct unit, wcs, filter)
            frames = self.datacubes[datacube_name].frames_for_filters(filters, convolve=self.config.spectral_convolution)

            # Add the observed images to the dictionary
            for filter_name, frame in zip(filter_names, frames): images[filter_name] = frame # these frames can be RemoteFrames if the datacube was a RemoteDataCube

            # Add the observed image dictionary for this datacube to the total dictionary (with the datacube name as a key)
            self.images[datacube_name] = images

    # -----------------------------------------------------------------

    def convolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Convolving the images ...")

        # Check whether a WCS was provided. If not, show a warning and skip the convolution
        if self.wcs is None:
            log.warning("WCS of the image is not defined, so convolution cannot be performed (the pixelscale is undefined)")
            return

        # Loop over the images
        for datacube_name in self.images:
            for filter_name in self.images[datacube_name]:

                # Check if the name of the image filter is a key in the 'kernel_paths' dictionary. If not, don't convolve.
                if filter_name not in self.kernel_paths or self.kernel_paths[filter_name] is None: continue

                # Load the kernel
                kernel = ConvolutionKernel.from_file(self.kernel_paths[filter_name])

                # Debugging
                log.debug("Convolving the '" + filter_name + "' image of the '" + datacube_name + "' instrument ...")

                # Convolve this image frame
                self.images[datacube_name][filter_name].convolve(kernel)

    # -----------------------------------------------------------------

    def convert_units(self):

        """
        This function ...
        :return:
        """

        # TODO: right now, this is just an implementation of the conversion from W / (m2 * arcsec2 * micron) to MJy/sr
        # 1 Jy = 1e-26 * W / (m2 * Hz)

        # Inform the user
        log.info("Converting the units of the images to " + str(self.unit) + " ...")

        assert str(self.unit) == "MJy / sr" # TEMPORARY

        # Get the pixelscale
        #pixelscale = self.wcs.average_pixelscale.to("arcsec/pix").value # in arcsec**2 / pixel

        # Loop over the images
        for datacube_name in self.images:
            for filter_name in self.images[datacube_name]:

                # Debugging
                log.debug("Converting the unit of the " + filter_name + " image of the '" + datacube_name + "' instrument ...")

                # Get the pivot wavelength of the filter
                fltr = self.filters[filter_name]
                pivot = fltr.pivot

                # Determine the conversion factor
                conversion_factor = 1.0

                # From surface brightness to flux density (no)
                #conversion_factor *=

                # From W / (m2 * arcsec2 * micron) to W / (m2 * arcsec2 * Hz)
                conversion_factor *= (pivot ** 2 / speed_of_light).to("micron/Hz").value

                # From W / (m2 * arcsec2 * Hz) to MJy / sr
                #conversion_factor *= (Unit("W/(m2 * arcsec2 * Hz)") / Unit("MJy/sr")).to("")
                conversion_factor *= 1e26 * 1e-6 * (u("sr") / u("arcsec2")).to("")

                # Convert
                self.images[datacube_name][filter_name] *= conversion_factor
                self.images[datacube_name][filter_name].unit = "MJy/sr"

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the images ...")

        # Loop over the different images (self.images is a nested dictionary of dictionaries)
        for datacube_name in self.images:
            for filter_name in self.images[datacube_name]:

                # Determine the path to the output FITS file
                path = fs.join(self.output_path, datacube_name + "__" + filter_name + ".fits")

                # Save the image
                self.images[datacube_name][filter_name].saveto(path)

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
