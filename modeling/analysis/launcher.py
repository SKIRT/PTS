#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.launch Contains the BestModelLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools import tables, time
from ...core.tools import filesystem as fs
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.tools.logging import log
from ...core.launch.options import AnalysisOptions
from ...core.launch.options import SchedulingOptions
from ...core.launch.options import LoggingOptions
from ...core.launch.estimate import RuntimeEstimator
from ...core.launch.parallelization import Parallelization
from ...core.simulation.remote import SkirtRemote
from ...magic.misc.kernels import AnianoKernels

# -----------------------------------------------------------------

class AnalysisLauncher(AnalysisComponent):
    
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
        super(AnalysisLauncher, self).__init__(config)

        # -- Attributes --

        # The remote SKIRT environment
        self.remote = SkirtRemote()

        # The name of the analysis run
        self.analysis_run_name = None

        # The path to the analysis run directory
        self.analysis_run_path = None

        # The ski file
        self.ski = None

        # The wavelength grid
        self.wavelength_grid = None

        # The dust grid
        self.dust_grid = None

        # The instruments
        self.instruments = dict()

        # The parallelization scheme
        self.parallelization = None

        # The scheduling options
        self.scheduling_options = None

        # The analysis options
        self.analysis_options = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        self.create_directories()

        # 2. Create the wavelength grid
        self.create_wavelength_grid()

        # 3. Create the dust grid
        self.create_dust_grid()

        # 4. Create the instruments
        self.create_instruments()

        # 5. Adjust the ski file
        self.adjust_ski()

        # 6. Set parallelization
        self.set_parallelization()

        # 7. Estimate the runtime for the simulation
        if self.remote.scheduler: self.estimate_runtime()

        # 9. Set the analysis options
        self.set_analysis_options()

        # 10. Writing
        self.write()

        # 11. Launch the simulation
        self.launch()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisLauncher, self).setup()

        # Setup the remote execution environment
        self.remote.setup(self.config.remote)

        # Generate a name for this analysis run
        self.analysis_run_name = time.unique_name()

        # Create a directory for this analysis run
        self.analysis_run_path = fs.join(self.analysis_path, self.analysis_run_name)

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Create the emission lines instance
        emission_lines = EmissionLines()

        # Fixed wavelengths in the grid
        fixed = [self.i1_filter.pivotwavelength(), self.fuv_filter.pivotwavelength()]

        # Create the grid
        grid, subgrid_npoints, emission_npoints, fixed_npoints = create_one_wavelength_grid(self.config.nwavelengths, emission_lines, fixed)

        # Set the grid
        self.wavelength_grid = grid

    # -----------------------------------------------------------------

    def create_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust grid ...")

        # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
        major_angular = self.truncation_ellipse.major  # major axis length of the sky ellipse
        radius_physical = (major_angular * self.galaxy_properties.distance).to("pc", equivalencies=dimensionless_angles())

        # Get the pixelscale in physical units
        distance = self.galaxy_properties.distance
        pixelscale_angular = self.reference_wcs.average_pixelscale * Unit("pix")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        x_radius = radius_physical
        y_radius = radius_physical
        z_radius = 3. * Unit("kpc")

        x_min = - x_radius
        x_max = x_radius
        y_min = - y_radius
        y_max = y_radius
        z_min = - z_radius
        z_max = z_radius

        x_extent = x_max - x_min
        y_extent = y_max - y_min
        z_extent = z_max - z_min

        # Set the scale
        scale = self.config.dg.rel_scale * pixelscale

        # Create the grid
        grid = create_one_dust_grid(self.config.dg.grid_type, scale, x_extent, y_extent, z_extent, x_min, x_max, y_min, y_max, z_min, z_max, self.config.dg.min_level, self.config.dg.max_mass_fraction)

        # Set the grid
        self.dust_grid = grid

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Loop over the projections
        for projection in ["earth", "faceon", "edgeon"]:

            # Debugging
            log.debug("Creating a full instrument for the " + projection + " projection ...")

            # Create the instrument and add it to the dictionary
            self.instruments[projection] = self.create_instrument("full", projection)

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: self.ski.add_instrument(name, self.instruments[name])

        # Set the number of photon packages
        self.ski.setpackages(self.config.packages)

        # Enable all writing options for analysis
        #self.ski.enable_all_writing_options()

        # Write out the dust grid data
        self.ski.set_write_grid()

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the parallelization scheme ...")

        # If the host uses a scheduling system
        if self.remote.scheduler:

            # Debugging
            log.debug("Remote host (" + self.remote.host_id + ") uses a scheduling system; determining parallelization scheme based on the requested number of nodes (" + str(self.config.nnodes) + ") ...")

            # Create the parallelization scheme from the host configuration and the requested number of nodes
            self.parallelization = Parallelization.for_host(self.remote.host, self.config.nnodes)

        # If the remote host does not use a scheduling system
        else:

            # Debugging
            log.debug("Remote host (" + self.remote.host_id + ") does not use a scheduling system; determining parallelization scheme based on the current load of the system and the requested number of cores per process (" + str(self.config.cores_per_process) + ") ...")

            # Get the amount of (currently) free cores on the remote host
            cores = int(self.remote.free_cores)

            # Determine the number of thread to be used per core
            threads_per_core = self.remote.threads_per_core if self.remote.use_hyperthreading else 1

            # Create the parallelization object
            self.parallelization = Parallelization.from_free_cores(cores, self.config.cores_per_process, threads_per_core)

        # Debugging
        log.debug("Parallelization scheme that will be used: " + str(self.parallelization))

    # -----------------------------------------------------------------

    def estimate_runtime(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the runtime for the simulation based on the timing of previous simulations ...")

        # Create a RuntimeEstimator instance
        estimator = RuntimeEstimator(self.timing_table)

        # Estimate the runtime for the configured number of photon packages and the configured remote host
        runtime = estimator.runtime_for(self.config.remote, self.ski, self.parallelization)

        # Debugging
        log.debug("The estimated runtime for the simulation is " + str(runtime) + " seconds")

        # Create the scheduling options, set the walltime
        self.scheduling_options = SchedulingOptions()
        self.scheduling_options.walltime = runtime

    # -----------------------------------------------------------------

    def set_analysis_options(self):

        """
        This function ...
        :return:
        """

        # Get the names of the filters for which we have photometry
        filter_names = self.get_observed_filter_names()

        # Set the paths to the for each image (except for the SPIRE images)
        kernel_paths = dict()
        aniano = AnianoKernels()
        pacs_red_psf_path = aniano.get_psf_path(self.pacs_red_filter)
        for filter_name in filter_names:
            if "SPIRE" in filter_name: continue
            kernel_paths[filter_name] = pacs_red_psf_path

        # Analysis options
        self.analysis_options = AnalysisOptions()

        # Set options for extraction
        self.analysis_options.extraction.path = self.analysis_extr_path
        self.analysis_options.extraction.progress = True
        self.analysis_options.extraction.timeline = True
        self.analysis_options.extraction.memory = True

        # Set options for plotting
        self.analysis_options.plotting.path = self.analysis_plot_path
        self.analysis_options.plotting.progress = True
        self.analysis_options.plotting.timeline = True
        self.analysis_options.plotting.seds = True
        self.analysis_options.plotting.grids = True
        self.analysis_options.plotting.reference_sed = self.observed_sed_path

        # Set miscellaneous options
        self.analysis_options.misc.path = self.analysis_misc_path
        self.analysis_options.misc.rgb = True
        self.analysis_options.misc.wave = True
        self.analysis_options.misc.fluxes = True
        self.analysis_options.misc.images = True
        self.analysis_options.misc.observation_filters = filter_names
        self.analysis_options.misc.observation_instruments = ["earth"]
        self.analysis_options.misc.make_images_remote = "nancy"
        self.analysis_options.misc.images_wcs = self.reference_wcs
        self.analysis_options.misc.images_kernels = kernel_paths
        self.analysis_options.misc.images_unit = "MJy/sr"

        # Set the paths of the timing and memory table files
        self.analysis_options.timing_table_path = self.timing_table_path
        self.analysis_options.memory_table_path = self.memory_table_path

        # Set the modeling path
        self.analysis_options.modeling_path = self.config.path

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the wavelength grid
        self.write_wavelength_grid()

        # Write the dust grid
        self.write_dust_grid()

        # Write the ski file
        self.write_ski()

    # -----------------------------------------------------------------

    def write_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grid ...")

        # Write the wavelength grid
        self.wavelength_grid.to_skirt_input(self.analysis_wavelengths_path)

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        # Save the ski file
        self.ski.saveto(self.analysis_ski_path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Create the simulation definition
        definition = SingleSimulationDefinition(self.analysis_ski_path, self.analysis_in_path, self.analysis_out_path)

        # Create the logging options
        logging = LoggingOptions(verbose=True, memory=True)

        # Debugging: save the screen output in a text file
        remote_skirt_dir_path = self.remote.skirt_dir
        remote_skirt_run_debug_path = fs.join(remote_skirt_dir_path, "run-debug")
        if not self.remote.is_directory(remote_skirt_run_debug_path): self.remote.create_directory(remote_skirt_run_debug_path)
        screen_output_path = fs.join(remote_skirt_run_debug_path, time.unique_name("screen") + ".txt")

        # Save the script file for manual inspection
        host_id = self.config.remote
        scripts_host_path = fs.join(self.analysis_scripts_path, host_id)
        if not fs.is_directory(scripts_host_path): fs.create_directory(scripts_host_path)
        local_script_path = fs.join(scripts_host_path, time.unique_name() + ".sh")

        # Run the simulation
        simulation = self.remote.run(definition, logging, self.parallelization, scheduling_options=self.scheduling_options,
                                     analysis_options=self.analysis_options, screen_output_path=screen_output_path,
                                     local_script_path=local_script_path)

        # Set the retrieve types
        simulation.retrieve_types = ["log", "sed", "image-total"]

        # Save the simulation file
        simulation.save()

# -----------------------------------------------------------------
