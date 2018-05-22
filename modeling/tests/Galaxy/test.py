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
from pts.core.units.parsing import parse_unit as u
from pts.core.simulation.skifile import SkiFile
from pts.core.basics.map import Map
from pts.do.commandline import Command
from pts.core.units.parsing import parse_quantity, parse_angle
from pts.modeling.basics.instruments import FullInstrument
from pts.modeling.basics.properties import GalaxyProperties
from pts.magic.basics.coordinate import SkyCoordinate
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.magic.basics.stretch import PixelStretch
from pts.magic.basics.coordinate import PixelCoordinate
from pts.magic.basics.pixelscale import Pixelscale
from pts.magic.convolution.aniano import AnianoKernels
from pts.core.filter.filter import parse_filter
from pts.core.test.implementation import TestImplementation
from pts.core.basics.log import log
from pts.modeling.basics.models import SersicModel3D, ExponentialDiskModel3D, RingModel3D
from pts.core.prep.wavelengthgrids import WavelengthGridGenerator
from pts.core.prep.dustgrids import DustGridGenerator
from pts.modeling.modeling.galaxy import fitting_filter_names

# -----------------------------------------------------------------

# Determine path of this directory
this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "fitting a mock spiral galaxy"

# -----------------------------------------------------------------

# Determine a name for this galaxy
fake_name = "GALAXY X"

# -----------------------------------------------------------------

# Dictionary with the parameters for the spiral disks for different components
spiral_parameters = dict()

# Young stars
young_parameters = Map()
young_parameters.radial_scale = parse_quantity("2600 pc")
young_parameters.axial_scale = parse_quantity("260 pc")
young_parameters.radial_truncation = parse_quantity("2e4 pc")
young_parameters.axial_truncation = parse_quantity("1e4 pc")
young_parameters.inner_radius = parse_quantity("0 pc")
young_parameters.arms = 2
young_parameters.pitch = parse_angle("18 deg")
young_parameters.radius = parse_quantity("500 pc")
young_parameters.phase = parse_angle("0 deg")
young_parameters.perturbation_weight = 0.3
young_parameters.index = 1
spiral_parameters["young"] = young_parameters

# Ionizing stars
ionizing_parameters = Map()
ionizing_parameters.radial_scale = parse_quantity("2600 pc")
ionizing_parameters.axial_scale = parse_quantity("150 pc")
ionizing_parameters.radial_truncation = parse_quantity("2e4 pc")
ionizing_parameters.axial_truncation = parse_quantity("1e4 pc")
ionizing_parameters.inner_radius = parse_quantity("0 pc")
ionizing_parameters.arms = 2
ionizing_parameters.pitch = parse_angle("18 deg")
ionizing_parameters.radius = parse_quantity("500 pc")
ionizing_parameters.phase = parse_angle("0 deg")
ionizing_parameters.perturbation_weight = 1
ionizing_parameters.index = 1
spiral_parameters["ionizing"] = ionizing_parameters

# Dust
dust_parameters = Map()
dust_parameters.radial_scale = parse_quantity("2600 pc")
dust_parameters.axial_scale = parse_quantity("150 pc")
dust_parameters.radial_truncation = parse_quantity("2e4 pc")
dust_parameters.axial_truncation = parse_quantity("1e4 pc")
dust_parameters.inner_radius = parse_quantity("0 pc")
dust_parameters.arms = 2
dust_parameters.pitch = parse_angle("18 deg")
dust_parameters.radius = parse_quantity("500 pc")
dust_parameters.phase = parse_angle("0 deg")
dust_parameters.perturbation_weight = 1
dust_parameters.index = 1
spiral_parameters["dust"] = dust_parameters

# -----------------------------------------------------------------

titles = dict()
titles["bulge"] = "Evolved stellar bulge"
titles["old"] = "Evolved stellar disk"
titles["young"] = "Young stars"
titles["ionizing"] = "Ionizing stars"
titles["dust"] = "Dust disk"

# -----------------------------------------------------------------

# Bulge
bulge_parameters = Map()
bulge_parameters.flattening = 0.7
bulge_parameters.index = 3.8
bulge_parameters.radius = parse_quantity("1000 pc")
bulge_parameters.metallicity = 0.02
bulge_parameters.age = 10
bulge_parameters.luminosity = 2.2e10 # bolometric

# Young SED
young_sed_parameters = Map()
young_sed_parameters.metallicity = 0.02
young_sed_parameters.age = 6
young_luminosity = 2.4e10 # bolometric # FREE PARAMETER

# Ionizing ring
ionizing_ring_parameters = Map()
ionizing_ring_parameters.radius = parse_quantity("6000 pc")
ionizing_ring_parameters.width = parse_quantity("3000 pc")
ionizing_ring_parameters.height = parse_quantity("150 pc")
ionizing_ring_parameters.metallicity = 0.02
ionizing_ring_parameters.luminosity = 2e9 # bol

# Ionizing SED (disk)
ionizing_sed_parameters = Map()
ionizing_sed_parameters.metallicity = 0.02
ionizing_luminosity = 3e9 # bol # FREE PARAMETER

# Dust mix
dust_mix = Map()
dust_mix.name = "Zubko"
dust_mix.write = False
dust_mix.write_mean = False
dust_mix.write_size = False
dust_mix.graphite_populations = 7
dust_mix.silicate_populations = 7
dust_mix.pah_populations = 5

# Dust ring
dust_ring_parameters = Map()
dust_ring_parameters.radius = parse_quantity("6000 pc")
dust_ring_parameters.width = parse_quantity("3000 pc")
dust_ring_parameters.height = parse_quantity("150 pc")
dust_ring_parameters.mass = parse_quantity("136e5 Msun")
dust_ring_parameters.mix = dust_mix

# Dust disk
dust_mass = parse_quantity("204e5 Msun") # FREE PARAMETER

# -----------------------------------------------------------------

class GalaxyTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(GalaxyTest, self).__init__(*args, **kwargs)

        # Free parameters for fitting
        self.free_parameters = Map()

        # Path to the ski file for the reference simulation
        self.reference_ski_path = None
        self.simulation_output_path = None

        # The reference ski file
        self.ski = None

        # Galaxy properties
        self.properties = None

        # Galaxy properties
        self.wcs = None

        # Configurables
        self.launcher = None
        self.flux_calculator = None
        self.image_maker = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Set parameters
        self.set_parameters()

        # Set galaxy info
        self.set_galaxy_properties()

        # Create WCS
        self.create_wcs()

        # Create galaxy components
        self.create_components()

        # Create instrument

        # Create wavelength grid
        self.create_wavelength_grid()

        # Create dust grid
        self.create_dust_grid()

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

        # Reference ski path
        self.reference_ski_path = fs.join(self.path, "galaxy_clumpy.ski")

        # Determine the simulation output path
        self.simulation_output_path = fs.join(self.path, "ref")

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

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

        # Inform the user
        log.info("Setting galaxy properties ...")

        # Galaxy poperties
        galaxy_distance = parse_quantity("3.63 Mpc")
        galaxy_inclination = Angle(59, "deg")
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

        # Inform the user
        log.info("Creating WCS ...")

        # Create WCS
        size = PixelStretch(1000, 1000)
        center_pixel = PixelCoordinate(500, 500)
        center_sky = self.properties.center
        pixelscale = Pixelscale(parse_quantity("2 arcsec"))
        self.wcs = CoordinateSystem.from_properties(size, center_pixel, center_sky, pixelscale)

    # -----------------------------------------------------------------

    def create_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the galaxy components ...")

        # Create a Sersic model for the bulge
        #self.bulge = SersicModel3D.from_2d(self.components["bulge"], self.properties.inclination, self.disk_pa, azimuth_or_tilt=self.config.bulge_deprojection_method)

        # Create an exponential disk model for the disk
        #self.disk = ExponentialDiskModel3D.from_2d(self.components["disk"], self.properties.inclination, self.disk_pa)

    # -----------------------------------------------------------------

    def create_dust_grid(self):

        """
        This function ...
        :return:
        """

        generator = DustGridGenerator()

        # <CartesianDustGrid writeGrid="true" minX="-2e4 pc" maxX="2e4 pc" minY="-2e4 pc" maxY="2e4 pc" minZ="-500 pc" maxZ="500 pc">
        # <meshX type="MoveableMesh">
        # <SymPowMesh numBins="50" ratio="25"/>
        #</meshX>
        #                <meshY type="MoveableMesh">
        #                    <SymPowMesh numBins="50" ratio="25"/>
        #                </meshY>
        #                <meshZ type="MoveableMesh">
        #                    <SymPowMesh numBins="20" ratio="45"/>
        #                </meshZ>
        #            </CartesianDustGrid>

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        generator = WavelengthGridGenerator()

        # <LogWavelengthGrid writeWavelengths="true" minWavelength="0.1 micron" maxWavelength="1000 micron" points="50"/>

    # -----------------------------------------------------------------

    def create_ski(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Creating ski file ...")

        fraction = 0.5
        count = 1000
        radius = parse_quantity("10 pc")
        cutoff = False
        kernel_type = "uniform"

        # Determine the ski path
        ski_path = fs.join(this_dir_path, "galaxy.ski")

        # Create clumpy ski file
        self.ski = SkiFile(ski_path)

        # Set instrument
        self.set_instrument()

        # Set stellar components
        self.set_stellar_components()

        # Set dust components
        self.set_dust_components()

        # Save as new ski file
        self.ski.saveto(self.reference_ski_path)

    # -----------------------------------------------------------------

    def set_instrument(self):

        """
        This function ...
        :return:
        """

        # Add full instrument that writes out photon counts
        self.ski.remove_all_instruments()

        # Create a full instrument
        distance = self.properties.distance
        inclination = self.properties.inclination
        azimuth = parse_angle("0 deg")
        position_angle = self.properties.position_angle
        field_x = parse_quantity("55000 pc")
        field_y = parse_quantity("5500 pc")
        pixels_x = self.config.nxpixels
        pixels_y = self.config.nypixels
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
        self.ski.add_instrument("earth", instrument)

    # -----------------------------------------------------------------

    def set_stellar_components(self):

        """
        This function ...
        :return:
        """

        # Remove all stellar components
        self.ski.remove_all_stellar_components()

        # self, component_id, radius, perturbation_weight, arms=1, pitch=0.1745329252, , phase=0, index=1
        self.ski.add_stellar_component_clumpiness(component_id, radius, perturbation_weight, arms, pitch, phase, index)

        # Add clumpiness to all stellar components
        for component_id in self.ski.get_stellar_component_ids():
            self.ski.add_stellar_component_clumpiness(component_id, fraction, count, radius, cutoff, kernel_type)

    # -----------------------------------------------------------------

    def set_dust_components(self):

        """
        This function ...
        :return:
        """

        # Remove all dust components
        self.ski.remove_all_dust_components()

        # Add clumpiness to all dust components
        for component_id in self.ski.get_dust_component_ids():
            self.ski.add_dust_component_clumpiness(component_id, fraction, count, radius, cutoff, kernel_type)

    # -----------------------------------------------------------------

    def launch_reference(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the reference simulation ...")

        # Settings
        settings_launch = dict()
        settings_launch["ski"] = self.reference_ski_path
        settings_launch["output"] = self.simulation_output_path
        settings_launch["create_output"] = True

        # Input
        input_launch = dict()

        # Launch command
        launch = Command("launch_simulation", "launch the reference simulation", settings_launch, input_launch, cwd=".")
        self.launcher = self.run_command(launch)

    # -----------------------------------------------------------------

    def make_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the observed mock SED ...")

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
        self.flux_calculator = self.run_command(calculate)

        # Determine the path to the mock SED
        mock_sed_path = "spiral_earth_fluxes.dat"

    # -----------------------------------------------------------------

    def make_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the observed mock images ...")

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
        self.image_maker = self.run_command(make)

        # MAKE DATA:

        # Create directory for the images
        ref_path = fs.create_directory_in(self.path, "ref")
        images_path = fs.create_directory_in(ref_path, "images")

        # Determine the name of the datacube
        datacube_names = self.image_maker.images.keys()
        if len(datacube_names) > 1: raise RuntimeError("Unexpected number of datacubes")
        datacube_name = datacube_names[0]

        # Loop over the images
        for filter_name in self.image_maker.images[datacube_name]:

            # Get the image
            image = self.image_maker.images[datacube_name][filter_name]

            # Save the image
            image_path = fs.join(images_path, filter_name + ".fits")
            image.saveto(image_path)

    # -----------------------------------------------------------------

    def setup_modelling(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting up the modelling ...")

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

        # Inform the user
        log.info("Performing the modelling ...")

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
