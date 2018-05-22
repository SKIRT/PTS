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
import numpy as np

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.modeling.models import Sersic2D

# Import the relevant PTS classes and modules
from pts.core.basics.log import log
from pts.do.commandline import Command
from pts.magic.core.mask import Mask
from pts.core.tools import filesystem as fs
from pts.modeling.tests.base import m81_data_path
from pts.modeling.basics.properties import GalaxyProperties
from pts.core.filter.filter import parse_filter
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.magic.catalog.fetcher import CatalogFetcher
from pts.core.tools import stringify
from pts.magic.tools import wavelengths
from pts.magic.catalog.extended import ExtendedSourceCatalog
from pts.magic.catalog.point import PointSourceCatalog
from pts.magic.tools.catalogs import get_galaxy_info
from pts.magic.region.ellipse import PixelEllipseRegion
from pts.magic.basics.coordinate import PixelCoordinate
from pts.magic.basics.stretch import PixelStretch
from pts.magic.tools import catalogs
from pts.core.remote.remote import Remote
from pts.magic.tests.base import SourcesTestBase

# -----------------------------------------------------------------

description = "Test the source detection and extraction"

# -----------------------------------------------------------------

# Determine the path to the headers directory
headers_path = fs.join(m81_data_path, "headers")

# -----------------------------------------------------------------

class SourcesTest(SourcesTestBase):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SourcesTest, self).__init__(*args, **kwargs)

        # The galaxy properties
        self.properties = None

        # The rotation masks
        self.rotation_masks = dict()

        # The catalog fetcher
        self.fetcher = CatalogFetcher()

        # Catalogs
        self.extended_source_catalog = None

        # The real FWHMs
        self.real_fwhms = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load galaxy properties
        self.load_properties()

        # 3. Load the coordinate system
        self.load_coordinate_systems()

        # 4. Initialize the frames
        self.initialize_frames()

        # 5. Set FWHMs
        self.set_fwhms()

        # 6. Make rotation mask
        self.make_rotation_masks()

        # 7. Get catalogs
        self.get_catalogs()

        # 8. Generate the sources
        self.make_sources()

        # 9. Make noise
        self.make_noise(self.rotation_masks)

        # 10. Make sky
        self.make_sky()

        # 11. Create the dataset
        self.create_dataset(self.rotation_masks)

        # 12. Create the directories
        self.create_directories()

        # 13. Find sources
        self.find()

        # 14. Extract sources
        self.extract()

        # 15. Write
        self.write()

        # 16. Plot
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

        # Set remote
        if self.config.remote is not None:
            self.remote = Remote()
            self.remote.setup(self.config.remote)

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
            self.coordinate_systems.append(wcs, fltr=fltr)

            # Break the loop if we have enough
            if nfilters_stars == self.config.nfilters_stars and nfilters_extra == self.config.nfilters_extra: break

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
            else: mask = Mask.empty_like(frame)

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
        self.extended_source_catalog = self.fetcher.get_extended_source_catalog(self.coordinate_systems.bounding_box)

        # Fetch
        self.point_source_catalog = self.fetcher.get_point_source_catalog(self.coordinate_systems.bounding_box, self.coordinate_systems.min_pixelscale, self.config.point_source_catalogs)

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
        self.add_principal_galaxy()

        # Point sources
        self.point_source_catalog = PointSourceCatalog()

    # -----------------------------------------------------------------

    def add_principal_galaxy(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding principal galaxy ...")

        # Get info
        name, position, redshift, galaxy_type, names, distance, inclination, d25, major, minor, pa = get_galaxy_info("M81", self.properties.center)

        ra = position.ra
        dec = position.dec

        principal = True
        companions = []
        parent = None

        # Add to the catalog
        self.extended_source_catalog.add_entry("M81", ra, dec, redshift, galaxy_type, names, distance, inclination, d25, major, minor, pa, principal, companions, parent)

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

        # Loop over the filters
        for fltr in self.frames:

            # Get the frame
            frame = self.frames[fltr]

            # Loop over the sources
            for index in range(len(self.extended_source_catalog)):

                principal = self.extended_source_catalog["Principal"][index]

                # Determine flux
                if principal: central_flux = 1e2
                else: central_flux = np.random.uniform(10,50)

                # Get pixel position
                coordinate = self.extended_source_catalog.get_position(index).to_pixel(frame.wcs)

                initial_sersic_amplitude = central_flux
                initial_sersic_x_0 = coordinate.x
                initial_sersic_y_0 = coordinate.y

                if principal:

                    s4g_name, pa, ellipticity, n, re, mag = catalogs.get_galaxy_s4g_one_component_info("M81")

                    angle = Angle(pa, "deg")
                    angle_deg = pa

                    effective_radius = re.to("arcsec").value / frame.average_pixelscale.to("arcsec").value

                    initial_sersic_n = n
                    initial_sersic_r_eff = effective_radius
                    initial_sersic_ellip = ellipticity
                    initial_sersic_theta = np.deg2rad(angle_deg)

                    # 1 / axial_ratio = 1 - ellipticity
                    axial_ratio = 1. / (1. - ellipticity)

                else:

                    # Get position angle and axes lengths
                    pa = self.extended_source_catalog["Posangle"][index]
                    major = self.extended_source_catalog["Major"][index]
                    minor = self.extended_source_catalog["Minor"][index]
                    axial_ratio = major / minor

                    angle = Angle(pa, "deg")
                    angle_deg = pa

                    effective_radius = major

                    # Produce guess values
                    initial_sersic_r_eff = effective_radius
                    initial_sersic_n = self.config.galaxy_sersic_index
                    initial_sersic_ellip = (axial_ratio - 1.0) / axial_ratio
                    initial_sersic_theta = np.deg2rad(angle_deg)

                # Produce sersic model from guess parameters, for time trials
                sersic_x, sersic_y = np.meshgrid(np.arange(frame.xsize), np.arange(frame.ysize))
                sersic_model = Sersic2D(amplitude=initial_sersic_amplitude, r_eff=initial_sersic_r_eff,
                                        n=initial_sersic_n, x_0=initial_sersic_x_0,
                                        y_0=initial_sersic_y_0, ellip=initial_sersic_ellip,
                                        theta=initial_sersic_theta)

                sersic_map = sersic_model(sersic_x, sersic_y)

                # Limit galaxy?
                if self.config.limit_galaxy:

                    limit_radius = self.config.galaxy_relative_asymptotic_radius * effective_radius

                    # Create galaxy region
                    galaxy_center = PixelCoordinate(initial_sersic_x_0, initial_sersic_y_0)
                    galaxy_radius = PixelStretch(limit_radius, limit_radius / axial_ratio)
                    galaxy_region = PixelEllipseRegion(galaxy_center, galaxy_radius, angle)

                    galaxy_mask = galaxy_region.to_mask(frame.xsize, frame.ysize)

                    # Set galaxy map zero outside certain radius
                    sersic_map[galaxy_mask.inverse()] = 0.0

                # Add
                frame += sersic_map

            # mask
            frame[self.rotation_masks[fltr]] = 0.0

    # -----------------------------------------------------------------

    def make_point_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making point sources ...")

        # Call the appropriate function
        if self.config.vary_fwhm: self.make_point_sources_variable_fwhm(self.rotation_masks)
        else: self.make_point_sources_fixed_fwhm(self.rotation_masks)

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

    def find(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding the sources ...")

        # Settings
        settings = dict()
        #settings["input"] =
        settings["output"] = self.find_path
        settings["nprocesses"] = self.config.nprocesses

        # Input
        input_dict = dict()
        input_dict["dataset"] = self.dataset
        input_dict["extended_source_catalog"] = self.extended_source_catalog
        input_dict["point_source_catalog"] = self.point_source_catalog
        input_dict["output_paths"] = self.find_paths

        # Construct the command
        command = Command("find_sources", "find sources", settings, input_dict)

        # Run the command
        self.finder = self.run_command(command, remote=self.remote)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
