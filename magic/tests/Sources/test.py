#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.utils import lazyproperty
from astropy.modeling.models import Gaussian2D, AiryDisk2D
from photutils.datasets import make_noise_image

# Import the relevant PTS classes and modules
from pts.core.tools.logging import log
from pts.do.commandline import Command
from pts.core.test.implementation import TestImplementation
from pts.magic.core.frame import Frame
from pts.magic.core.mask import Mask
from pts.core.tools import filesystem as fs
from pts.magic.core.dataset import DataSet
from pts.core.basics.unit import parse_unit as u
from pts.modeling.tests.base import m81_data_path
from pts.modeling.basics.properties import GalaxyProperties
from pts.core.filter.filter import parse_filter
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.magic.catalog.fetcher import CatalogFetcher
from pts.magic.basics.stretch import SkyStretch
from pts.magic.region.rectangle import SkyRectangleRegion
from pts.magic.tools import coordinates
from pts.magic.basics.coordinate import SkyCoordinate
from pts.core.tools import stringify
from pts.magic.tools import wavelengths
from pts.magic.catalog.extended import ExtendedSourceCatalog
from pts.magic.catalog.point import PointSourceCatalog
from pts.magic.tools import fitting, statistics
from pts.magic.convolution.kernels import has_variable_fwhm, get_fwhm

# -----------------------------------------------------------------

description = "Test the source extraction"

# -----------------------------------------------------------------

# For M81
fwhms = {"2MASS H": 4.640929858306589 * u("arcsec"),
          "2MASS J": 4.580828087551186 * u("arcsec"),
          "2MASS Ks": 4.662813601376219 * u("arcsec"),
          "SDSS g": 2.015917936060279 * u("arcsec"),
          "SDSS i": 1.85631074608032 * u("arcsec"),
          "SDSS r": 2.026862297071852 * u("arcsec"),
          "SDSS u": 2.327165667182196 * u("arcsec"),
          "SDSS z": 1.841443699129355 * u("arcsec")}

# -----------------------------------------------------------------

# Determine the path to the headers directory
headers_path = fs.join(m81_data_path, "headers")

# -----------------------------------------------------------------

class SourcesTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(SourcesTest, self).__init__(config, interactive)

        # Paths
        self.data_path = None
        self.find_path = None
        self.extract_path = None

        # The galaxy properties
        self.properties = None

        # The coordinate systems
        self.coordinate_systems = dict()

        # The frames
        self.frames = dict()

        # The rotation masks
        self.rotation_masks = dict()

        # The dataset
        self.dataset = None

        # The catalog fetcher
        self.fetcher = CatalogFetcher()

        # Catalogs
        self.extended_source_catalog = None
        self.point_source_catalog = None

        # The real FWHMs
        self.real_fwhms = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function
        self.setup(**kwargs)

        # Load galaxy properties
        self.load_properties()

        # Load the coordinate system
        self.load_coordinate_systems()

        # Initialize the frames
        self.initialize_frames()

        # Set FWHMs
        self.set_fwhms()

        # 2. Make rotation mask
        self.make_rotation_masks()

        # 3. Get catalogs
        self.get_catalogs()

        # 3. Generate the sources
        self.make_sources()

        # 4. Make noise
        self.make_noise()

        # Make sky
        self.make_sky()

        # Create the dataset
        self.create_dataset()

        # FInd
        self.find()

        # Extract
        self.extract()

        # Write
        self.write()

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SourcesTest, self).setup(**kwargs)

        # Set paths
        self.data_path = fs.create_directory_in(self.path, "data")
        self.find_path = fs.create_directory_in(self.path, "find")
        self.extract_path = fs.create_directory_in(self.path, "extract")

    # -----------------------------------------------------------------

    def load_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the galaxy properties ...")

        # Determine the path
        path = fs.join(m81_data_path, "properties.dat")

        # Load
        self.properties = GalaxyProperties.from_file(path)

    # -----------------------------------------------------------------

    def load_coordinate_systems(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the coordinate systems ...")

        nfilters_stars = 0
        nfilters_extra = 0

        # Loop over the header files
        for path, name in fs.files_in_path(headers_path, extension="txt", returns=["path", "name"]):

            # Get the filter and wavelength
            fltr = parse_filter(name)
            wavelength = fltr.effective if fltr.effective is not None else fltr.center

            # SKip Planck
            if fltr.observatory == "Planck": continue

            # Wavelength greater than 25 micron
            if wavelength > wavelengths.ranges.ir.mir.max:

                if nfilters_extra == self.config.nfilters_extra: continue
                else: nfilters_extra += 1

            # Wavelength smaller than 25 micron
            else:

                if nfilters_stars == self.config.nfilters_stars: continue
                else: nfilters_stars += 1

            # Debugging
            log.debug("Loading the coordinate system for the '" + str(fltr) + "' filter ...")

            # Get WCS
            wcs = CoordinateSystem.from_header_file(path)

            # Add the coordinate system
            self.coordinate_systems[fltr] = wcs

            # Break the loop if we have enough
            if nfilters_stars == self.config.nfilters_stars and nfilters_extra == self.config.nfilters_extra: break

    # -----------------------------------------------------------------

    @lazyproperty
    def min_ra(self):

        """
        This function ...
        :return:
        """

        min_ra = None

        # Loop over the coordinate systems
        for fltr in self.coordinate_systems:

            wcs = self.coordinate_systems[fltr]
            wcs_min_ra = wcs.min_ra

            if min_ra is None or min_ra > wcs_min_ra: min_ra = wcs_min_ra

        return min_ra

    # -----------------------------------------------------------------

    @property
    def min_ra_deg(self):

        """
        This function ...
        :return:
        """

        return self.min_ra.to("deg").value

    # -----------------------------------------------------------------

    @lazyproperty
    def max_ra(self):

        """
        This function ...
        :return:
        """

        max_ra = None

        # Loop over the coordinate systems
        for fltr in self.coordinate_systems:

            wcs = self.coordinate_systems[fltr]
            wcs_max_ra = wcs.max_ra

            if max_ra is None or max_ra < wcs_max_ra: max_ra = wcs_max_ra

        return max_ra

    # -----------------------------------------------------------------

    @lazyproperty
    def max_ra_deg(self):

        """
        This function ...
        :return:
        """

        return self.max_ra.to("deg").value

    # -----------------------------------------------------------------

    @lazyproperty
    def ra_center(self):

        """
        This function ...
        :return:
        """

        return 0.5 * (self.min_ra + self.max_ra)

    # -----------------------------------------------------------------

    @lazyproperty
    def ra_center_deg(self):

        """
        This function ...
        :return:
        """

        return self.ra_center.to("deg").value

    # -----------------------------------------------------------------

    @lazyproperty
    def ra_range(self):

        """
        This function ...
        :return:
        """

        the_range = None

        # Loop over the coordinate systems
        for fltr in self.coordinate_systems:

            wcs = self.coordinate_systems[fltr]

            if the_range is None: the_range = wcs.ra_range
            else: the_range.adjust(wcs.ra_range)

        # Return the range
        return the_range

    # -----------------------------------------------------------------

    @lazyproperty
    def min_dec(self):

        """
        This function ...
        :return:
        """

        min_dec = None

        # Loop over the coordinate systems
        for fltr in self.coordinate_systems:

            wcs = self.coordinate_systems[fltr]
            wcs_min_dec = wcs.min_dec

            if min_dec is None or min_dec > wcs_min_dec: min_dec = wcs_min_dec

        # Return
        return min_dec

    # -----------------------------------------------------------------

    @lazyproperty
    def min_dec_deg(self):

        """
        This function ...
        :return:
        """

        return self.min_dec.to("deg").value

    # -----------------------------------------------------------------

    @lazyproperty
    def max_dec(self):

        """
        This function ...
        :return:
        """

        max_dec = None

        # Loop over the coordinate systems
        for fltr in self.coordinate_systems:

            wcs = self.coordinate_systems[fltr]
            wcs_max_dec = wcs.max_dec

            if max_dec is None or max_dec < wcs_max_dec: max_dec = wcs_max_dec

        # Return
        return max_dec

    # -----------------------------------------------------------------

    @lazyproperty
    def dec_center(self):

        """
        This function ...
        :return:
        """

        dec_center = 0.5 * (self.min_dec + self.max_dec)
        return dec_center

    # -----------------------------------------------------------------

    @lazyproperty
    def dec_center_deg(self):

        """
        This function ...
        :return:
        """

        return self.dec_center.to("deg").value

    # -----------------------------------------------------------------

    @lazyproperty
    def max_dec_deg(self):

        """
        This function ...
        :return:
        """

        return self.max_dec.to("deg").value

    # -----------------------------------------------------------------

    @lazyproperty
    def dec_range(self):

        """
        This function ...
        :return:
        """

        the_range = None

        # Loop over the coordinate systems
        for fltr in self.coordinate_systems:

            wcs = self.coordinate_systems[fltr]

            if the_range is None: the_range = wcs.dec_range
            else: the_range.adjust(wcs.dec_range)

        # Return the dec range
        return the_range

    # -----------------------------------------------------------------

    @lazyproperty
    def center(self):

        """
        This function ...
        :return:
        """

        return SkyCoordinate(self.ra_center_deg, self.dec_center_deg, unit="deg")

    # -----------------------------------------------------------------

    @lazyproperty
    def coordinate_range(self):

        """
        This function ...
        :return:
        """

        # Calculate the actual RA and DEC distance in degrees
        ra_distance = abs(coordinates.ra_distance(self.dec_center_deg, self.min_ra_deg, self.max_ra_deg))
        dec_distance = abs(self.max_dec_deg - self.min_dec_deg)

        # Create RA and DEC span as quantities
        ra_span = ra_distance * u("deg")
        dec_span = dec_distance * u("deg")

        # Return the center coordinate and the RA and DEC span
        return self.center, ra_span, dec_span

    # -----------------------------------------------------------------

    @lazyproperty
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Get coordinate range
        center, ra_span, dec_span = self.coordinate_range

        # Create box
        radius = SkyStretch(0.5 * ra_span, 0.5 * dec_span)
        box = SkyRectangleRegion(center, radius)

        # Return the box
        return box

    # -----------------------------------------------------------------

    @lazyproperty
    def min_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the coordinate systems
        for fltr in self.coordinate_systems:

            wcs = self.coordinate_systems[fltr]
            if pixelscale is None or wcs.average_pixelscale < pixelscale: pixelscale = wcs.average_pixelscale

        # Return the minimum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def max_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the coordinate systems
        for fltr in self.coordinate_systems:

            wcs = self.coordinate_systems[fltr]
            if pixelscale is None or wcs.average_pixelscale > pixelscale: pixelscale = wcs.average_pixelscale

        # Return the maximum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    def initialize_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the frames ...")

        # Loop over the filters
        for fltr in self.coordinate_systems:

            # Debugging
            log.debug("Initializing the '" + str(fltr) + "' frame ...")

            # Get the wcs
            wcs = self.coordinate_systems[fltr]

            # Create new frame
            frame = Frame.zeros(wcs.shape)

            # Add the wcs
            frame.wcs = wcs

            # Set the filter
            frame.filter = fltr

            # Set the unit
            frame.unit = "Jy"

            # Add the frame
            self.frames[fltr] = frame

    # -----------------------------------------------------------------

    def set_fwhms(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the FWHMs ...")

        # Loop over the filters
        for fltr in self.frames:

            # Get the fwhm
            if has_variable_fwhm(fltr): fwhm = fwhms[str(fltr)]
            else: fwhm = get_fwhm(fltr)

            # Debugging
            log.debug("The FWHM of the '" + str(fltr) + "' image is " + stringify.stringify(fwhm)[1])

            # Set
            self.real_fwhms[fltr] = fwhm

    # -----------------------------------------------------------------

    def make_rotation_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making rotation masks ...")

        # Loop over the filters
        for fltr in self.frames:

            # Debugging
            log.info("Making rotation mask for the '" + str(fltr) + "' image ...")

            # Get the frame
            frame = self.frames[fltr]

            # Rotate
            if self.config.rotate:

                # Choose a random rotation angle
                angle = Angle(np.random.uniform(-90, 90), "deg")

                # Debugging
                log.debug("The random rotation angle is '" + stringify.stringify(angle)[1] + "'")

                # Create mask
                mask = frame.rotation_mask(angle)

            # Don't rotate
            else: mask = Mask.empty(self.config.shape[1], self.config.shape[0])

            # Set the mask
            self.rotation_masks[fltr] = mask

    # -----------------------------------------------------------------

    def get_catalogs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting catalogs of point and extended sources ...")

        # Adding catalogued sources
        if self.config.add_catalogued_sources: self.fetch_catalogs()
        else: self.initialize_catalogs()

        # Add random sources
        self.create_random_sources()

    # -----------------------------------------------------------------

    def fetch_catalogs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching catalogs ...")

        # Get extended source catalog
        self.extended_source_catalog = self.fetcher.get_extended_source_catalog(self.bounding_box)

        # Fetch
        self.point_source_catalog = self.fetcher.get_point_source_catalog(self.bounding_box, self.min_pixelscale, self.config.point_source_catalogs)

    # -----------------------------------------------------------------

    def initialize_catalogs(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Initializing catalogs ...")

        # Extended sources
        self.extended_source_catalog = ExtendedSourceCatalog()

        # Add principal galaxy

        # Point sources
        self.point_source_catalog = PointSourceCatalog()

    # -----------------------------------------------------------------

    @property
    def star_filters(self):

        """
        This function ...
        :return:
        """

        filters = []

        # Loop over the filters
        for fltr in self.frames:

            # Get wavelength
            wavelength = fltr.effective if fltr.effective is not None else fltr.center

            # Check
            if wavelength > wavelengths.ranges.ir.mir.max: continue
            filters.append(fltr)

        # Return the filters
        return filters

    # -----------------------------------------------------------------

    @property
    def extra_filters(self):

        """
        This function ...
        :return:
        """

        filters = []

        # Loop over the filters
        for fltr in self.frames:

            # Get wavelength
            wavelength = fltr.effective if fltr.effective is not None else fltr.center

            # Check
            if wavelength < wavelengths.ranges.ir.mir.max: continue
            filters.append(fltr)

        # Return the filters
        return filters

    # -----------------------------------------------------------------

    def create_random_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating random point sources ...")

        # Generate random coordinates
        right_ascensions = np.random.uniform(self.min_ra_deg, self.max_ra_deg, size=self.config.nrandom_sources)
        declinations = np.random.uniform(self.min_dec_deg, self.max_dec_deg, size=self.config.nrandom_sources)

        # Loop over the coordinates
        for ra, dec in zip(right_ascensions, declinations):

            # Create a sky coordinate
            coordinate = SkyCoordinate(ra=ra, dec=dec, unit="deg")

            # Add to the point source catalog
            self.point_source_catalog.add_coordinate(coordinate)

    # -----------------------------------------------------------------

    def make_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making sources ...")

        # Extended
        self.make_extended_sources()

        # Point
        self.make_point_sources()

    # -----------------------------------------------------------------

    def make_extended_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making extended sources ...")

    # -----------------------------------------------------------------

    def make_point_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making point sources ...")

        # Loop over the 'star' filters
        for fltr in self.star_filters:

            # Debugging
            log.debug("Making point sources for the '" + str(fltr) + "' image ...")

            # Get y and x
            #y, x = np.indices(self.frames[fltr].shape)

            # Get pixelscale
            pixelscale = self.frames[fltr].average_pixelscale

            # Debugging
            log.debug("The pixelscale of the image is " + stringify.stringify(pixelscale)[1])

            # Determine the FWHM in pixel coordinates
            fwhm_pix = self.real_fwhms[fltr].to("arcsec").value / pixelscale.to("arcsec").value

            counter = 0

            # Loop over the coordinates in the point sources catalog
            for coordinate in self.point_source_catalog.coordinates():

                counter += 1

                # Debugging
                log.debug("Adding source " + str(counter) + " of " + str(len(self.point_source_catalog)) + " ...")

                # Check whether it falls in the frame
                if not self.frames[fltr].contains(coordinate): continue

                # Convert into pixel coordinate
                pixel_coordinate = coordinate.to_pixel(self.frames[fltr].wcs)

                # Check whether not masked
                if self.rotation_masks[fltr][int(pixel_coordinate.y), int(pixel_coordinate.x)]: continue

                # Generate random deviation
                x_deviation = np.random.normal(0.0, 1.)
                y_deviation = np.random.normal(0.0, 1.)

                # Debugging
                log.debug("Random pixel position deviation is (" + str(x_deviation) + ", " + str(y_deviation) + ")")

                # Alter pixel coordinate
                pixel_coordinate.x += x_deviation
                pixel_coordinate.y += y_deviation

                # Generate random deviation from FWHM
                fwhm_deviation = np.random.normal(0.0, 0.05 * fwhm_pix)

                # Debugging
                log.debug("Random FWHM deviation (on a FWHM of " + str(fwhm_pix) + ") is " + str(fwhm_deviation))

                # Add the deviation
                fwhm_pix += fwhm_deviation

                # Generate a random amplitude
                amplitude = float(np.random.uniform(5., 100.))

                # Determine sigma
                sigma = statistics.fwhm_to_sigma * fwhm_pix

                sigma_level = 3.

                # Determine x and y range for evaluation
                min_x = max(int(math.floor(pixel_coordinate.x - sigma_level * sigma)), 0)
                max_x = min(int(math.ceil(pixel_coordinate.x + sigma_level * sigma)), self.frames[fltr].xsize - 1)
                min_y = max(int(math.floor(pixel_coordinate.y - sigma_level * sigma)), 0)
                max_y = min(int(math.ceil(pixel_coordinate.y + sigma_level * sigma)), self.frames[fltr].ysize - 1)

                # Create relative centers
                rel_x = pixel_coordinate.x - min_x
                rel_y = pixel_coordinate.y - min_y

                # 2D GAussian model
                if self.config.psf_model == "gaussian":

                    # Create the model
                    model = Gaussian2D(amplitude=amplitude, x_mean=rel_x, y_mean=rel_y, x_stddev=sigma, y_stddev=sigma)

                # Airy disk model
                elif self.config.psf_model == "airydisk":

                    # Determine the radius
                    radius = fitting.gaussian_sigma_to_airy_radius(sigma)

                    # Create the model
                    model = AiryDisk2D(amplitude=amplitude, x_0=rel_x, y_0=rel_y, radius=radius)

                # Invalid
                else: raise ValueError("Not a valid model")

                # Evaluate the model
                y, x = np.indices((max_y - min_y, max_x - min_x))
                data = model(x, y)

                # Add the data
                self.frames[fltr][min_y:max_y, min_x:max_x] += data

    # -----------------------------------------------------------------

    def make_noise(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding noise ...")

        # Loop over the frames
        for fltr in self.frames:

            # Get the frame
            frame = self.frames[fltr]

            # Make noise and add it
            data = make_noise_image(frame.shape, type='gaussian', mean=0., stddev=self.config.noise_stddev)
            frame += data

            # Mask
            frame[self.rotation_masks[fltr]] = 0.0

    # -----------------------------------------------------------------

    def make_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding sky ...")

        # Loop over the frames
        for fltr in self.frames:

            # Get the frame
            frame = self.frames[fltr]

            # Determine random sky level
            sky_level = np.random.uniform(0.0, 10.)

            # Create sky gradient
            y, x = np.mgrid[:frame.ysize, :frame.xsize]
            sky_gradient =  x * y

            # Normalize so that the maximal sky is 20%
            sky_gradient = sky_gradient / np.max(sky_gradient) * 20

            # Add the sky
            frame += sky_level + sky_gradient

            # Mask
            frame[self.rotation_masks[fltr]] = 0.0

    # -----------------------------------------------------------------

    def make_galaxy(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding galaxy ...")

    # -----------------------------------------------------------------

    def create_dataset(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dataset ...")

        # Initialize
        self.dataset = DataSet()

        # Add the frames
        for fltr in self.frames:

            # Determine name for the image
            name = str(fltr)

            # Determine path for this frame
            path = fs.join(self.data_path, name + ".fits")

            # Save the frame
            self.frames[fltr].saveto(path)

            # Debugging
            log.debug("Adding the '" + name + "' image to the dataset ...")

            # Add the frame to the dataset
            self.dataset.add_path(name, path)

        # Determine database path
        path = fs.join(self.path, "database.dat")

        # Write the dataset
        self.dataset.saveto(path)

    # -----------------------------------------------------------------

    def find(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding the sources ...")

        # Settings
        settings = dict()

        # Input
        input_dict = dict()

        input_dict["dataset"] = self.dataset

        # Construct the command
        command = Command("find_sources", "find sources", settings, input_dict, cwd=None, finish=None)

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    :param temp_path:
    :return:
    """

    pass

# -----------------------------------------------------------------
