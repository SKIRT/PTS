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

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.unit import parse_unit as u
from pts.core.tools.logging import log
from pts.do.commandline import Command
from pts.magic.misc.kernels import AnianoKernels
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.modeling.basics.instruments import FullInstrument
from pts.core.prep.wavelengthgrids import WavelengthGridGenerator
from pts.core.prep.dustgrids import DustGridGenerator
from pts.core.basics.quantity import parse_quantity
from pts.magic.region.list import SkyRegionList
from pts.core.filter.filter import parse_filter
from pts.core.remote.moderator import PlatformModerator
from pts.core.simulation.memory import MemoryRequirement
from pts.core.prep.deploy import Deployer
from pts.core.launch.options import AnalysisOptions
from pts.modeling.tests.base import M81TestBase, m81_data_path, fitting_filter_names

# -----------------------------------------------------------------

# Determine path of this directory
this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "creating the best model of the M81 galaxy based on mock observations"

# -----------------------------------------------------------------

# Rough estimates
# memory requirement: 3.16 GB
serial_memory = 2.0 * u("GB")
parallel_memory = 2.0 * u("GB")

# -----------------------------------------------------------------

class M81Test(M81TestBase):

    """
    This class runs the test on M81
    """

    def __init__(self, config=None, interactive=False):

        """
        This function ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(M81Test, self).__init__(config, interactive)

        # The coordinate systems for the various filters
        self.wcs_paths = None

        # The host ID for remote execution of reference simulation
        self.host_id = None

        # The image maker
        self.image_maker = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the properties
        self.load_properties()

        # 3. Load the components
        self.load_components()

        # 4. Load the wcs
        self.create_wcs()

        # Load the input maps
        self.load_maps()

        # Create the instrument
        self.create_instrument()

        # Create deprojections
        self.create_deprojections()

        # 5. Create wavelength grid
        self.create_wavelength_grid()

        # Create dust grid
        self.create_dust_grid()

        # 6. Create the ski file
        self.create_ski()

        # Write
        self.write()

        # Plot
        self.plot()

        # 7. Launch reference simulation
        self.launch_reference()

        # 10. Setup modelling
        self.setup_modelling()

        # 11. Model
        self.model()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(M81Test, self).setup(**kwargs)

        # Set remote host for reference simulation
        self.set_remote()

        # Deploy SKIRT and PTS
        if self.host_id is not None: self.deploy()

    # -----------------------------------------------------------------

    def set_remote(self):

        """
        This function ...
        :return:
        """

        # Setup the platform moderator
        moderator = PlatformModerator()

        # Launching the reference simulation
        #if self.config.local: moderator.add_local("reference")
        #elif self.config.remotes is not None: moderator.add_single("reference", self.config.remotes)
        #elif self.modeling_config.host_ids is None: moderator.add_local("reference")
        #else: moderator.add_single("other", self.modeling_config.host_ids)
        if self.config.host_ids is None: moderator.add_local("reference")
        else: moderator.add_single("reference", self.config.host_ids)

        # Run the platform moderator
        moderator.run()

        # Set the host ID
        self.host_id = moderator.host_id_for_single("reference")

    # -----------------------------------------------------------------

    def deploy(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deploying SKIRT and PTS ...")

        # Create the deployer
        deployer = Deployer()

        # Set the host ids
        deployer.config.host_ids = [self.host_id]

        # Set the host id on which PTS should be installed (on the host for extra computations and the fitting hosts
        # that have a scheduling system to launch the pts run_queue command)
        #deployer.config.pts_on = self.moderator.all_host_ids

        # Set
        #deployer.config.check = self.config.check_versions

        # Set
        #deployer.config.update_dependencies = self.config.update_dependencies

        # Run the deployer
        deployer.run()

    # -----------------------------------------------------------------

    def create_wcs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the WCS ...")

        min_pixelscale = None

        # Determine the path to the headers directory
        headers_path = fs.join(m81_data_path, "headers")

        # Loop over the header files
        for path, name in fs.files_in_path(headers_path, extension="txt", returns=["path", "name"]):

            # Get the filter
            fltr = parse_filter(name)
            filter_name = str(fltr)

            # Set the path of the file for the filter name
            self.wcs_paths[filter_name] = path

            # Get WCS
            wcs = CoordinateSystem.from_header_file(path)

            # Adjust the pixelscale
            if min_pixelscale is None:
                min_pixelscale = wcs.pixelscale
                self.wcs = wcs
            elif min_pixelscale > wcs.pixelscale:
                min_pixelscale = wcs.pixelscale
                self.wcs = wcs

    # -----------------------------------------------------------------

    def create_wcs_not_working(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the WCS ...")

        # Determine the path to the headers directory
        headers_path = fs.join(m81_data_path, "headers")

        # Loop over the files in the directory
        ra_range = None
        dec_range = None
        min_pixelscale = None
        for path, name in fs.files_in_path(headers_path, extension="txt", returns=["path", "name"]):

            # Get the filter
            fltr = parse_filter(name)

            # Get WCS
            wcs = CoordinateSystem.from_header_file(path)

            # Adjust RA range
            if ra_range is None: ra_range = wcs.ra_range
            else: ra_range.adjust(wcs.ra_range)

            # Adjust DEC range
            if dec_range is None: dec_range = wcs.dec_range
            else: dec_range.adjust(wcs.dec_range)

            # Adjust the pixelscale
            if min_pixelscale is None: min_pixelscale = wcs.pixelscale
            elif min_pixelscale > wcs.pixelscale: min_pixelscale = wcs.pixelscale

        # Create coordinate system
        # size, center_pixel, center_sky, pixelscale
        self.wcs = CoordinateSystem.from_ranges(ra_range, dec_range, min_pixelscale)

    # -----------------------------------------------------------------

    def create_instrument(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instrument ...")

        # Create a full instrument
        # center, distance, inclination, azimuth, position_angle
        azimuth = 0.0
        self.instrument = FullInstrument.from_wcs(self.wcs, self.galaxy_center, self.galaxy_distance, self.galaxy_inclination, azimuth, self.galaxy_position_angle)

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Create the wavelength generator
        generator = WavelengthGridGenerator()

        # Set input
        input_dict = dict()
        input_dict["ngrids"] = 1
        input_dict["npoints"] = self.config.nwavelengths
        input_dict["fixed"] = [self.i1_filter.pivot, self.fuv_filter.pivot]
        input_dict["add_emission_lines"] = True
        input_dict["lines"] = ["Halpha"] # only the H-alpha line is of importance
        input_dict["min_wavelength"] = parse_quantity("0.1 micron")
        input_dict["max_wavelength"] = parse_quantity("1000 micron")
        input_dict["filters"] = [parse_filter(string) for string in fitting_filter_names]

        # Run the generator
        generator.run(**input_dict)

        # Set the wavelength grid
        self.wavelength_grid = generator.single_grid

    # -----------------------------------------------------------------

    def create_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust grid ...")

        # Create the dust grid generator
        generator = DustGridGenerator()

        # Determine truncation ellipse
        disk_ellipse_path = fs.join(m81_data_path, "components", "disk.reg")
        disk_ellipse = SkyRegionList.from_file(disk_ellipse_path)[0]
        truncation_ellipse = 0.82 * disk_ellipse

        # Determine the radius of the galaxy
        semimajor_angular = truncation_ellipse.semimajor  # semimajor axis length of the sky ellipse
        radius_physical = (semimajor_angular * self.galaxy_distance).to("pc", equivalencies=dimensionless_angles())

        # Set properties
        generator.grid_type = "bintree"  # set grid type
        generator.x_radius = radius_physical
        generator.y_radius = radius_physical
        generator.z_radius = 3. * u("kpc")

        # Set input
        input_dict = dict()
        input_dict["ngrids"] = 1
        input_dict["scale"] = self.config.dust_grid_relative_scale * self.deprojections["dust"].pixelscale # in pc
        input_dict["level"] = self.config.dust_grid_min_level
        input_dict["mass_fraction"] = self.config.dust_grid_max_mass_fraction

        # Generate the grid
        generator.run(**input_dict)

        # Set the dust grid
        self.dust_grid = generator.single_grid

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ski file
        self.write_ski()

        # Write the input
        self.write_input()

        # Write the WCS
        self.write_wcs()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the wavelengths
        self.plot_wavelengths()

        # Plot the filters
        self.plot_filters()

    # -----------------------------------------------------------------

    def write_wcs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the reference coordinate system ...")

        # Write the WCS
        self.wcs.saveto(self.reference_wcs_path)

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
        settings_launch["ski"] = self.ski_path
        settings_launch["input"] = self.simulation_input_path
        settings_launch["output"] = self.simulation_output_path
        settings_launch["create_output"] = True
        settings_launch["remote"] = self.host_id
        settings_launch["attached"] = True
        settings_launch["progress_bar"] = True

        # Create the analysis options
        analysis = AnalysisOptions()
        analysis.extraction.path = self.simulation_extract_path
        analysis.plotting.path = self.simulation_plot_path
        analysis.misc.path = self.simulation_misc_path
        analysis.extraction.progress = True
        analysis.extraction.timeline = True
        analysis.extraction.memory = True
        analysis.plotting.progress = True
        analysis.plotting.timeline = True
        analysis.plotting.memory = True
        analysis.plotting.seds = True
        analysis.plotting.grids = True
        seds_path = fs.join(m81_data_path, "seds")
        analysis.plotting.reference_seds = fs.files_in_path(seds_path)
        analysis.misc.fluxes = True
        analysis.misc.images = True
        analysis.misc.observation_filters = fitting_filter_names
        analysis.misc.observation_instruments = ["earth"]
        analysis.misc.make_images_remote = self.host_id
        analysis.misc.images_wcs = self.reference_wcs_path
        analysis.misc.images_unit = "Jy/pix"
        analysis.misc.spectral_convolution = False

        # Create Aniano kernels object
        aniano = AnianoKernels()

        # Set the paths to the kernel for each image
        kernel_paths = dict()
        for filter_name in fitting_filter_names: kernel_paths[filter_name] = aniano.get_psf_path(parse_filter(filter_name))
        analysis.misc.images_kernels = kernel_paths

        # Set the paths to the WCS files for each image
        analysis.misc.rebin_wcs = {"earth": self.wcs_paths}

        # Input
        input_launch = dict()
        #input_launch["memory"] = MemoryRequirement(serial_memory, parallel_memory)

        # Launch command
        launch = Command("launch_simulation", "launch the reference simulation", settings_launch, input_launch, cwd=".")
        self.launcher = self.run_command(launch)

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
        input_setup["ngc_name"] = self.properties.ngc_name
        input_setup["hyperleda_name"] = self.properties.hyperleda_name

        # Construct the command
        stp = Command("setup", "setup the modeling", settings_setup, input_setup, cwd=".")

        # Call the command
        tool = self.run_command(stp)

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
        command = Command("model", "perform the modelling", settings_model, input_model, "./Galaxy")

        # Run the command
        self.modeler = self.run_command(command)

# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    :param temp_path:
    :return:
    """

    pass

# -----------------------------------------------------------------
