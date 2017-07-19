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
from ...core.tools import filesystem as fs
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.tools.logging import log
from ...core.launch.options import AnalysisOptions
from ...core.launch.options import SchedulingOptions
from ...core.launch.options import LoggingOptions
from ...core.advanced.runtimeestimator import RuntimeEstimator
from ...core.simulation.parallelization import Parallelization
from ...core.simulation.remote import SkirtRemote
from ...magic.convolution.aniano import AnianoKernels
from ...core.advanced.dustgridtool import DustGridTool

# -----------------------------------------------------------------

class AnalysisLauncher(AnalysisComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(AnalysisLauncher, self).__init__(*args, **kwargs)

        # -- Attributes --

        # NEW: THE ANALYSIS RUN
        self.analysis_run = None

        # The remote SKIRT environment
        self.remote = SkirtRemote()

        # The path to the instruments directory
        self.run_instruments_path = None

        # Simulation directories
        self.run_output_path = None
        self.run_extr_path = None
        self.run_plot_path = None
        self.run_misc_path = None

        # Analysis directories
        self.run_attenuation_path = None
        self.run_colours_path = None
        self.run_residuals_path = None
        self.run_heating_path = None

        # The ski file
        self.ski = None

        # # The wavelength grid
        # self.wavelength_grid = None
        #
        # # The dust grid
        # self.dust_grid = None
        #
        # # The instruments
        # self.instruments = dict()

        # The parallelization scheme
        self.parallelization = None

        # The scheduling options
        self.scheduling_options = None

        # The analysis options
        self.analysis_options = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the ski file
        #self.load_ski()
        #self.load_input_paths()

        # 7. Adjust the ski file
        self.adjust_ski()

        # 8. Set parallelization
        self.set_parallelization()

        # 9. Estimate the runtime for the simulation
        if self.remote.scheduler: self.estimate_runtime()

        # 10. Set the analysis options
        self.set_analysis_options()

        # 11. Writing
        self.write()

        # 12. Launch the simulation
        self.launch()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisLauncher, self).setup(**kwargs)

        # NEW: GET THE RUN
        self.analysis_run = self.get_run(self.config.run)

        # Setup the remote execution environment
        self.remote.setup(self.config.remote)

    # -----------------------------------------------------------------

    @property
    def analysis_run_name(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.name

    # -----------------------------------------------------------------

    @property
    def analysis_run_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.path

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # Set the number of photon packages
        self.ski.setpackages(self.config.npackages)

        # Enable dust self-absorption
        self.ski.enable_selfabsorption()

        # Enable transient heating
        self.ski.set_transient_dust_emissivity()

        # Remove the existing instruments
        #self.ski.remove_all_instruments()

        # Add the instruments
        #for name in self.instruments: self.ski.add_instrument(name, self.instruments[name])

        # Enable all writing options for analysis
        #self.ski.enable_all_writing_options()

        # Write out the dust grid data
        #self.ski.set_write_grid()

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
            self.parallelization = Parallelization.for_host(self.remote.host, self.config.nnodes, self.config.data_parallel)

        # If the remote host does not use a scheduling system
        else:

            # Debugging
            log.debug("Remote host (" + self.remote.host_id + ") does not use a scheduling system; determining parallelization scheme based on the current load of the system and the requested number of cores per process (" + str(self.config.cores_per_process) + ") ...")

            # Get the amount of (currently) free cores on the remote host
            cores = int(self.remote.free_cores)

            # Determine the number of thread to be used per core
            threads_per_core = self.remote.threads_per_core if self.remote.use_hyperthreading else 1

            # Create the parallelization object
            self.parallelization = Parallelization.from_free_cores(cores, self.config.cores_per_process, threads_per_core, self.config.data_parallel)

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

        # Determine the number of dust cells by building the tree locally
        ncells = self.estimate_ncells()

        # Estimate the runtime for the configured number of photon packages and the configured remote host
        runtime = estimator.runtime_for(self.ski, self.parallelization, self.remote.host_id, self.remote.cluster_name, self.config.data_parallel, nwavelengths=len(self.wavelength_grid), ncells=ncells)

        # Debugging
        log.debug("The estimated runtime for the simulation is " + str(runtime) + " seconds")

        # Create the scheduling options, set the walltime
        self.scheduling_options = SchedulingOptions()
        self.scheduling_options.walltime = runtime

    # -----------------------------------------------------------------

    def estimate_ncells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the number of dust cells ...")

        # Create simulation directory and output directory
        simulation_path = fs.create_directory_in(self.analysis_run_path, "temp")

        # Initialize dust grid tool
        tool = DustGridTool()

        # Get the dust grid statistics
        statistics = tool.get_statistics(self.ski, simulation_path, self.maps_path, self.galaxy_name)
        return statistics.ncells

    # -----------------------------------------------------------------

    def set_analysis_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the analysis options ...")

        # Set the paths to the for each image (except for the SPIRE images)
        kernel_paths = dict()
        aniano = AnianoKernels()
        pacs_red_psf_path = aniano.get_psf_path(self.pacs_red_filter)
        for filter_name in self.observed_filter_names:
            if "SPIRE" in filter_name: continue
            kernel_paths[filter_name] = pacs_red_psf_path

        # Analysis options
        self.analysis_options = AnalysisOptions()

        # Set options for extraction
        self.analysis_options.extraction.path = self.run_extr_path
        self.analysis_options.extraction.progress = True
        self.analysis_options.extraction.timeline = True
        self.analysis_options.extraction.memory = True

        # Set options for plotting
        self.analysis_options.plotting.path = self.run_plot_path
        self.analysis_options.plotting.progress = True
        self.analysis_options.plotting.timeline = True
        self.analysis_options.plotting.seds = True
        self.analysis_options.plotting.grids = True
        self.analysis_options.plotting.reference_seds = [self.observed_sed_path]

        # Set miscellaneous options
        self.analysis_options.misc.path = self.run_misc_path
        self.analysis_options.misc.rgb = True
        self.analysis_options.misc.wave = True
        self.analysis_options.misc.fluxes = True
        self.analysis_options.misc.images = True
        self.analysis_options.misc.observation_filters = self.observed_filter_names  # the filters for which to create the observations
        self.analysis_options.misc.observation_instruments = ["earth"]
        self.analysis_options.misc.make_images_remote = "nancy"
        # TODO: #self.analysis_options.misc.images_wcs = self.reference_wcs
        self.analysis_options.misc.images_kernels = kernel_paths
        self.analysis_options.misc.images_unit = "MJy/sr"
        self.analysis_options.misc.make_images_remote = self.config.images_remote

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

        # Write the ski file
        self.write_ski()

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file to " + self.ski_file_path + "...")

        # Save the ski file
        self.ski.saveto(self.ski_file_path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulation ...")

        # Create the simulation definition
        definition = SingleSimulationDefinition(self.ski_file_path, self.run_output_path, self.input_paths)

        # Create the logging options
        logging = LoggingOptions(verbose=True, memory=True)

        # Debugging: save the screen output in a text file
        remote_skirt_dir_path = self.remote.skirt_dir
        remote_skirt_run_debug_path = fs.join(remote_skirt_dir_path, "run-debug")
        if not self.remote.is_directory(remote_skirt_run_debug_path): self.remote.create_directory(remote_skirt_run_debug_path)
        screen_output_path = fs.join(remote_skirt_run_debug_path, self.analysis_run_name + ".txt")

        # Determine the path to the launching script file for manual inspection
        local_script_path = fs.join(self.analysis_run_path, self.remote.host_id + ".sh")

        # Run the simulation
        simulation = self.remote.run(definition, logging, self.parallelization, name=self.analysis_run_name,
                                     scheduling_options=self.scheduling_options, analysis_options=self.analysis_options,
                                     screen_output_path=screen_output_path, local_script_path=local_script_path)

        # Set the retrieve types
        simulation.retrieve_types = ["log", "sed", "image-total"]

        # Save the simulation file
        simulation.save()

# -----------------------------------------------------------------
