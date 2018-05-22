#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.units.parsing import parse_unit as u
from pts.core.basics.log import log
from pts.do.commandline import Command
from pts.magic.convolution.aniano import AnianoKernels
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.modeling.basics.instruments import FullInstrument
from pts.core.filter.filter import parse_filter
from pts.core.remote.moderator import PlatformModerator
from pts.core.prep.deploy import Deployer
from pts.core.launch.options import AnalysisOptions
from pts.modeling.tests.base import M81TestBase, m81_data_path, fitting_filter_names, instrument_name
from pts.modeling.tests.base import seds_path, dustpedia_sed_path
from pts.core.data.sed import ObservedSED
from pts.core.tools import sequences
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.core.tools.utils import lazyproperty

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

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(M81Test, self).__init__(*args, **kwargs)

        # The coordinate systems for the various filters
        self.wcs_paths = None

        # The host ID for remote execution of reference simulation
        self.host_id = None

        # The image maker
        self.image_maker = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the properties
        self.load_properties()

        # 3. Load the components
        self.load_components()

        # 4. Load the wcs
        self.create_wcs()

        # 5. Load the input maps
        self.load_maps()

        # 6. Create the instrument
        self.create_instrument()

        # 7. Create deprojections
        self.create_deprojections()

        # 8. Create wavelength grid
        self.create_wavelength_grid()

        # 9. Create dust grid
        self.create_dust_grid()

        # 10. Create the ski file
        self.create_ski()

        # 11. Write
        self.write()

        # 12. Plot
        self.plot()

        # 13. Launch reference simulation
        self.launch_reference()

        # 14. Get the real parameter values
        self.get_real_parameter_values()

        # 15. Setup modelling
        self.setup_modelling()

        # 16. Model
        self.model()

        # 17. Get best parameter values
        self.get_best_parameter_values()

        # 18. Test
        self.test()

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

    @lazyproperty
    def host(self):

        """
        This function ...
        :return:
        """

        from pts.core.remote.host import load_host
        return load_host(self.host_id)

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
        deployer.config.hosts = [self.host]

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
        settings_launch["ski"] = self.reference_ski_path
        settings_launch["input"] = self.simulation_input_path
        settings_launch["output"] = self.simulation_output_path
        settings_launch["create_output"] = True
        settings_launch["remote"] = self.host_id
        settings_launch["attached"] = True
        settings_launch["show_progress"] = True

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
        analysis.plotting.reference_seds = fs.files_in_path(seds_path)
        analysis.misc.fluxes = True
        analysis.misc.images = True
        analysis.misc.observation_filters = fitting_filter_names
        analysis.misc.observation_instruments = [instrument_name]
        analysis.misc.make_images_remote = self.host_id
        analysis.misc.images_wcs = self.reference_wcs_path
        analysis.misc.images_unit = "Jy/pix"
        analysis.misc.spectral_convolution = self.config.spectral_convolution

        # Set flux error bars
        dustpedia_sed = ObservedSED.from_file(dustpedia_sed_path)
        filter_names = dustpedia_sed.filter_names()
        errors = dustpedia_sed.errors()
        #print([str(error) for error in errors])
        flux_errors = sequences.zip_into_dict(filter_names, [str(error) for error in errors])
        analysis.misc.flux_errors = flux_errors

        # Create Aniano kernels object
        aniano = AnianoKernels()

        # Set the paths to the kernel for each image
        kernel_paths = dict()
        for filter_name in fitting_filter_names: kernel_paths[filter_name] = aniano.get_psf_path(parse_filter(filter_name))
        analysis.misc.images_kernels = kernel_paths

        # Set the paths to the WCS files for each image
        analysis.misc.rebin_wcs = {instrument_name: self.wcs_paths}

        # Input
        input_launch = dict()
        #input_launch["memory"] = MemoryRequirement(serial_memory, parallel_memory)
        input_launch["analysis_options"] = analysis

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
        settings_setup["name"] = self.modeling_name
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
        settings_model["ngenerations"] = self.config.ngenerations
        settings_model["nsimulations"] = self.config.nsimulations
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
            images[filter_name] = fs.join(self.simulation_misc_path, filter_name + ".fits")
        input_model["images"] = images

        # Construct the command
        command = Command("model", "perform the modelling", settings_model, input_model, self.modeling_path)

        # Run the command
        self.modeler = self.run_command(command)

    # -----------------------------------------------------------------

    @property
    def modeling_environment(self):

        """
        This function ...
        :return: 
        """

        return GalaxyModelingEnvironment(self.modeling_path)

    # -----------------------------------------------------------------

    def test(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Testing ...")

        # Check best
        self.check_best()

        # Check database
        self.check_database()

        # Check statistics
        self.check_statistics()

        # Check images
        self.check_images()

    # -----------------------------------------------------------------

    def check_images(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the images ...")

# -----------------------------------------------------------------
