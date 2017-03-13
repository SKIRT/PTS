#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import inspect
import numpy as np

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.unit import parse_unit as u
from pts.core.basics.configuration import Configuration
from pts.core.simulation.skifile import LabeledSkiFile
from pts.core.basics.range import QuantityRange, RealRange
from pts.core.basics.map import Map
from pts.do.commandline import Command
from pts.core.basics.quantity import parse_quantity
from pts.modeling.basics.instruments import FullInstrument
from pts.modeling.basics.properties import GalaxyProperties
from pts.magic.basics.coordinate import SkyCoordinate
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.magic.basics.stretch import PixelStretch
from pts.magic.basics.coordinate import PixelCoordinate
from pts.magic.basics.pixelscale import Pixelscale
from pts.magic.misc.kernels import AnianoKernels
from pts.core.filter.filter import parse_filter
from pts.core.test.implementation import TestImplementation

# -----------------------------------------------------------------

# Determine path of this directory
this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "fitting a mock spiral galaxy"

# -----------------------------------------------------------------

# Set the filters for which the data can be used for fitting (data that can be trusted well enough)
fitting_filter_names = ["GALEX FUV", "GALEX NUV", "SDSS u", "SDSS g", "SDSS r", "SDSS i", "SDSS z", "WISE W1",
                        "IRAC I1", "IRAC I2", "WISE W2", "IRAC I3", "IRAC I4", "WISE W3", "WISE W4", "MIPS 24mu",
                        "Pacs blue", "Pacs red", "SPIRE PSW", "SPIRE PMW", "SPIRE PLW"]

# -----------------------------------------------------------------

# Determine a name for this galaxy
fake_name = "GALAXY X"

# -----------------------------------------------------------------

class GalaxyTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, path):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(GalaxyTest, self).__init__(path)

        self.free_parameters = Map()

        # Path to the ski file for the reference simulation
        self.reference_ski_path = None

        self.simulation_output_path = None

        # Galaxy properties
        self.properties = None

        # Galaxy properties
        self.wcs = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function
        self.setup(**kwargs)

        # Set paths
        self.set_paths()

        # Set parameters
        self.set_parameters()

        # Set galaxy info
        self.set_galaxy_properties()

        # Create WCS
        self.create_wcs()

        # Create ski file
        self.create_ski()

        # Launch reference simulation
        self.launch_reference()

        # Make observed SED
        self.make_sed()

        # Make observed images
        self.make_images()

        # Setup modelling
        self.setup_modelling()

        # Model
        self.model()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the
        super(GalaxyTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def set_paths(self):

        """
        This function ...
        :return:
        """

        # Reference ski path
        self.reference_ski_path = fs.join(self.path, "galaxy_clumpy.ski")

        # Determine the simulation output path
        self.simulation_output_path = fs.join(self.path, "ref")

    # -----------------------------------------------------------------

    def set_parameters(self):

        # Define the true values for the free parameters
        self.free_parameters.dust_mass = None
        self.free_parameters.fuv_ionizing = None
        self.free_parameters.fuv_young = None

    # -----------------------------------------------------------------

    def set_galaxy_properties(self):

        """
        This function ..
        :return:
        """

        # Galaxy poperties
        galaxy_distance = parse_quantity("3.63 Mpc")
        galaxy_inclination = Angle(59, "deg")
        #galaxy_azimuth = Angle(0, "deg")
        galaxy_pa = Angle(67, "deg")

        # Generate a random coordinate for the center of the galaxy
        ra_random = np.random.rand() * 360.0 * u("deg")
        dec_random = (np.random.rand() * 180.0 - 90.0) * u("deg")
        galaxy_center = SkyCoordinate(ra=ra_random, dec=dec_random)

        # Determine the galaxy size, convert to
        galaxy_size = parse_quantity("100000 lyr")
        galaxy_radius = 0.5 * galaxy_size.to("pc")
        galaxy_radius_arcsec = (galaxy_radius / galaxy_distance).to("arcsec", equivalencies=dimensionless_angles())

        # Determine ellipticity
        ellipticity = 0.5

        # Set the properties
        self.properties = GalaxyProperties(name=fake_name, ngc_name=fake_name, hyperleda_name=fake_name, galaxy_type=None,
                                      center=galaxy_center, major=galaxy_radius, major_arcsec=galaxy_radius_arcsec,
                                      ellipticity=ellipticity, position_angle=galaxy_pa, distance=galaxy_distance, distance_error=None,
                                      inclination=galaxy_inclination, redshift=None, common_name=fake_name)

    # -----------------------------------------------------------------

    def create_wcs(self):

        """
        This function ...
        :return:
        """

        # Create WCS
        size = PixelStretch(1000, 1000)
        center_pixel = PixelCoordinate(500, 500)
        center_sky = self.properties.center
        pixelscale = Pixelscale(parse_quantity("2 arcsec"))
        self.wcs = CoordinateSystem.from_properties(size, center_pixel, center_sky, pixelscale)

    # -----------------------------------------------------------------

    def create_ski(self):

        """
        This function ...
        """

        fraction = 0.5
        count = 1000
        radius = parse_quantity("10 pc")
        cutoff = False
        kernel_type = "uniform"

        # Determine the ski path
        ski_path = fs.join(this_dir_path, "galaxy.ski")

        # Create clumpy ski file
        ski = LabeledSkiFile(ski_path)

        # Add clumpiness to all stellar components
        for component_id in ski.get_stellar_component_ids():
            ski.add_stellar_component_clumpiness(component_id, fraction, count, radius, cutoff, kernel_type)

        # Add clumpiness to all dust components
        for component_id in ski.get_dust_component_ids():
            ski.add_dust_component_clumpiness(component_id, fraction, count, radius, cutoff, kernel_type)

        # Add full instrument that writes out photon counts
        ski.remove_all_instruments()

        # Create a full instrument
        distance = self.properties.distance
        inclination = self.properties.inclination
        #azimuth = self.properties.azimuth
        azimuth = Angle(0, "deg")
        position_angle = self.properties.position_angle
        field_x = parse_quantity("55000 pc")
        field_y = parse_quantity("5500 pc")
        pixels_x = 1000
        pixels_y = 1000
        center_x = parse_quantity("130 pc")
        center_y = parse_quantity("-181 pc")
        scattering_levels = 0
        counts = True  # write photon counts
        instrument = FullInstrument(distance=distance, inclination=inclination, azimuth=azimuth,
                                    position_angle=position_angle,
                                    field_x=field_x, field_y=field_y, pixels_x=pixels_x, pixels_y=pixels_y,
                                    center_x=center_x,
                                    center_y=center_y, scattering_levels=scattering_levels, counts=counts)

        # Add the instrument
        ski.add_instrument("earth", instrument)

        # Save as new ski file
        ski.saveto(self.reference_ski_path)

    # -----------------------------------------------------------------

    def launch_reference(self):

        """
        This function ...
        :return:
        """

        # Settings
        settings_launch = dict()
        settings_launch["ski"] = self.reference_ski_path
        settings_launch["output"] = self.simulation_output_path
        settings_launch["create_output"] = True

        # Input
        input_launch = dict()

        # Launch command
        launch = Command("launch_simulation", "launch the reference simulation", settings_launch, input_launch, cwd=".")
        instance = self.run_command(launch)

    # -----------------------------------------------------------------

    def make_sed(self):

        """
        This function ...
        :return:
        """

        # Settings
        settings_sed = dict()
        settings_sed["spectral_convolution"] = False

        # Input
        input_sed = dict()
        input_sed["simulation_output_path"] = self.simulation_output_path
        input_sed["output_path"] = "."

        # Construct the command
        #create_sed = Command("observed_fluxes", "create the mock SED", settings_sed, input_sed, cwd=".")

        # Add the command
        #commands.append(create_sed)

        # Launch command
        calculate = Command("observed_fluxes", "create the mock SED", settings_sed, input_sed, cwd=".")
        instance = self.run_command(calculate)

        # Determine the path to the mock SED
        mock_sed_path = "spiral_earth_fluxes.dat"

    # -----------------------------------------------------------------

    def make_images(self):

        """
        This function ...
        :return:
        """

        # Create Aniano kernels object
        aniano = AnianoKernels()

        # Set the paths to the kernel for each image
        kernel_paths = dict()
        for filter_name in fitting_filter_names: kernel_paths[filter_name] = aniano.get_psf_path(parse_filter(filter_name))

        # Settings
        settings_images = dict()
        settings_images["spectral_convolution"] = False
        # No output path is specified, so images won't be written out

        # Input
        input_images = dict()
        input_images["simulation_output_path"] = self.simulation_output_path
        input_images["output_path"] = "."
        input_images["filter_names"] = fitting_filter_names
        input_images["instrument_names"] = ["earth"]
        # input_images["wcs_path"] =
        input_images["wcs"] = self.wcs
        input_images["kernel_paths"] = kernel_paths
        input_images["unit"] = "Jy/pix"
        # input_images["host_id"] = "nancy"

        # Construct the command
        #create_images = Command("observed_images", "create the mock images", settings_images, input_images, cwd=".",
        #                        finish=make_data)

        # Add the command
        #commands.append(create_images)

        make = Command("observed_images", "create the mock images", settings_images, input_images, cwd=".")
        image_maker = self.run_command(make)

        # MAKE DATA:

        # Create directory for the images
        ref_path = fs.create_directory_in(self.path, "ref")
        images_path = fs.create_directory_in(ref_path, "images")

        # Determine the name of the datacube
        datacube_names = image_maker.images.keys()
        if len(datacube_names) > 1: raise RuntimeError("Unexpected number of datacubes")
        datacube_name = datacube_names[0]

        # Loop over the images
        for filter_name in image_maker.images[datacube_name]:

            # Get the image
            image = image_maker.images[datacube_name][filter_name]

            # Save the image
            image_path = fs.join(images_path, filter_name + ".fits")
            image.saveto(image_path)

    # -----------------------------------------------------------------

    def setup_modelling(self):

        """
        This function ...
        :return:
        """

        # Settings
        settings_setup = dict()
        settings_setup["type"] = "galaxy"
        settings_setup["name"] = "Galaxy"
        settings_setup["fitting_host_ids"] = None

        # Create input dict for setup
        input_setup = dict()
        input_setup["ngc_name"] = fake_name
        input_setup["hyperleda_name"] = fake_name

        # Construct the command
        #stp = Command("setup", "setup the modeling", settings_setup, input_setup, cwd=".")

        # Add the command
        #commands.append(stp)



    # -----------------------------------------------------------------

    def model(self):

        """
        This function ...
        :return:
        """

        # Settings
        settings_model = dict()
        settings_model["ngenerations"] = 4
        settings_model["nsimulations"] = 20
        settings_model["fitting_settings"] = {"spectral_convolution": False}

        # Input
        input_model = dict()

        # Set galaxy properties
        input_model["properties"] = self.properties

        # Set SEDs
        input_model["seds"] = dict()

        # Set images dictionary
        images = dict()
        for filter_name in fitting_filter_names:
            images[filter_name] = fs.join("../ref/images", filter_name + ".fits")
        input_model["images"] = images

        # Construct the command
        #command = Command("model", "perform the modelling", settings_model, input_model, "./Galaxy")

        # Add the command
        #commands.append(command)

# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    :param temp_path:
    :return:
    """

    return True

# -----------------------------------------------------------------
