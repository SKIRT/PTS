#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.explorer Contains the ParameterExplorer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.launch.batchlauncher import BatchLauncher
from .modelgenerators.grid import GridModelGenerator
from .modelgenerators.genetic import GeneticModelGenerator
from .modelgenerators.instinctive import InstinctiveModelGenerator
from ...core.tools import time
from .tables import ParametersTable, ChiSquaredTable, IndividualsTable
from ...core.launch.options import SchedulingOptions
from ...core.advanced.runtimeestimator import RuntimeEstimator
from ...core.tools.stringify import stringify_not_list, stringify
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.advanced.parallelizationtool import ParallelizationTool
from ...core.remote.host import Host
from ...core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter
from .evaluate import get_parameter_values_from_generator, prepare_simulation, generate_simulation_name, get_parameter_values_for_named_individual
from ...core.simulation.input import SimulationInput
from ...core.tools import introspection
from ...core.tools import parallelization as par

# -----------------------------------------------------------------

class GenerationInfo(object):

    """
    This function ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Set info
        self.name = kwargs.pop("name", None)
        self.index = kwargs.pop("index", None)
        self.method = kwargs.pop("method", None)
        self.wavelength_grid_level = kwargs.pop("wavelength_grid_level", None)
        #self.dust_grid_level = kwargs.pop("dust_grid_level", None)
        self.model_representation = kwargs.pop("model_representation", None)
        self.nsimulations = kwargs.pop("nsimulations", None)
        self.npackages = kwargs.pop("npackages", None)
        self.selfabsorption = kwargs.pop("selfabsorption", None)
        self.transient_heating = kwargs.pop("transient_heating", None)

        self.path = kwargs.pop("path", None)
        self.individuals_table_path = kwargs.pop("individuals_table_path", None)
        self.parameters_table_path = kwargs.pop("parameters_table_path", None)
        self.chi_squared_table_path = kwargs.pop("chi_squared_table_path", None)

# -----------------------------------------------------------------

class ParameterExplorer(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(ParameterExplorer, self).__init__(config, interactive)

        # -- Attributes --

        # The fitting run
        self.fitting_run = None

        # The ski template
        self.ski = None

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The generation info
        self.generation = GenerationInfo()

        # The individuals table
        self.individuals_table = None

        # The parameters table
        self.parameters_table = None

        # The chi squared table
        self.chi_squared_table = None

        # The parameter ranges
        self.ranges = dict()

        # The generation index and name
        self.generation_index = None
        self.generation_name = None

        # The model generator
        self.generator = None

        # The simulation input
        self.simulation_input = None

        # A dictionary with the scheduling options for the different remote hosts
        self.scheduling_options = dict()

        # The number of wavelengths used
        self.nwavelengths = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the ski template
        self.load_ski()

        # 3. Set the parameter ranges
        if not self.has_all_ranges: self.set_ranges()

        # 4. Set the generation info
        self.set_info()

        # 5. Create the generation directory
        self.create_generation_directory()

        # 6. Generate the model parameters
        self.generate_models()

        # 7. Set the paths to the input files
        if self.fitting_run.needs_input: self.set_input()

        # 8. Adjust the ski template
        self.adjust_ski()

        # 9. Set the parallelization scheme for local execution if necessary
        if self.only_local: self.set_parallelization_local()

        # 10. Set the parallelization schemes for the different remote hosts
        if self.uses_schedulers: self.set_parallelization_remote()

        # 11. Estimate the runtimes for the different remote hosts
        if self.uses_schedulers: self.estimate_runtimes()

        # 12. Fill the tables for the current generation
        self.fill_tables()

        # 13. Writing
        self.write()

        # 14. Launch the simulations for different parameter values
        self.launch()

    # -----------------------------------------------------------------

    @property
    def uses_remotes(self):

        """
        This function ...
        :return:
        """

        return self.launcher.uses_remotes

    # -----------------------------------------------------------------

    @property
    def only_local(self):

        """
        This function ...
        :return:
        """

        return self.launcher.only_local

    # -----------------------------------------------------------------

    @property
    def has_all_ranges(self):

        """
        This function ...
        :return:
        """

        # Loop over the free parameter labels
        for label in self.fitting_run.free_parameter_labels:

            # If range is already defined
            if label not in self.ranges: return False

        # All ranges are defined
        return True

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ParameterExplorer, self).setup(**kwargs)

        # Load the fitting run
        self.fitting_run = self.load_fitting_run(self.config.name)

        # Get ranges
        if "ranges" in kwargs: self.ranges = kwargs.pop("ranges")

        # Set options for the batch launcher
        self.set_launcher_options()

        # Set the model generator
        self.set_generator()

        # Check whether this is not the first generation so that we can use remotes with a scheduling system
        #if self.ngenerations == 0 and self.uses_schedulers:
        #    raise ValueError("The specified remote hosts cannot be used for the first generation: at least one remote uses a scheduling system")

        # Check whether initialize_fit has been called
        if self.modeling_type == "galaxy":
            if not fs.is_file(self.fitting_run.wavelength_grids_table_path): raise RuntimeError("Call initialize_fit_galaxy before starting the parameter exploration")
            if not fs.is_file(self.fitting_run.dust_grids_table_path): raise RuntimeError("Call initialize_fit_galaxy before starting the parameter exploration")

    # -----------------------------------------------------------------

    def set_launcher_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting options for the batch simulation launcher ...")

        # Set remote host IDs
        remote_host_ids = []
        if self.fitting_run.ngenerations == 0:
            for host_id in self.config.remotes:
                if Host(host_id).scheduler: log.warning("Not using remote host '" + host_id + "' for the initial generation because it uses a scheduling system for launching jobs")
                else: remote_host_ids.append(host_id)
        else: remote_host_ids = self.config.remotes

        # Basic options
        self.launcher.config.shared_input = True                               # The input directories (or files) for the different simulations are shared
        self.launcher.config.remotes = remote_host_ids                         # The remote host(s) on which to run the simulations
        self.launcher.config.attached = self.config.attached                   # Run remote simulations in attached mode
        self.launcher.config.group_simulations = self.config.group             # Group multiple simulations into a single job (because a very large number of simulations will be scheduled) TODO: IMPLEMENT THIS
        self.launcher.config.group_walltime = self.config.walltime             # The preferred walltime for jobs of a group of simulations
        self.launcher.config.timing_table_path = self.fitting_run.timing_table_path        # The path to the timing table file
        self.launcher.config.memory_table_path = self.fitting_run.memory_table_path        # The path to the memory table file
        self.launcher.config.cores_per_process = self.config.cores_per_process # The number of cores per process, for non-schedulers
        self.launcher.config.dry = self.config.dry                             # Dry run (don't actually launch simulations)
        self.launcher.config.progress_bar = True  # show progress bars for local execution

        # Simulation analysis options

        ## General
        self.launcher.config.analysis.relative = True

        ## Logging
        self.launcher.config.logging.verbose = True
        self.launcher.config.logging.memory = True
        #self.launcher.config.logging.allocation = True
        #self.launcher.config.logging.allocation_limit = 1e-5

        ## Extraction
        self.launcher.config.analysis.extraction.path = "extr"    # name of the extraction directory
        self.launcher.config.analysis.extraction.progress = True  # extract progress information
        self.launcher.config.analysis.extraction.timeline = True  # extract the simulation timeline
        self.launcher.config.analysis.extraction.memory = True    # extract memory information

        ## Plotting
        self.launcher.config.analysis.plotting.path = "plot"  # name of the plot directory
        self.launcher.config.analysis.plotting.seds = True    # Plot the output SEDs
        self.launcher.config.analysis.plotting.reference_seds = [self.observed_sed_path]  # the path to the reference SED (for plotting the simulated SED against the reference points)
        self.launcher.config.analysis.plotting.format = "pdf"  # plot format

        ## Miscellaneous
        self.launcher.config.analysis.misc.path = "misc"       # name of the misc output directory
        if self.is_images_modeling:

            self.launcher.config.analysis.misc.fluxes = False
            self.launcher.config.analysis.misc.images = True
            self.launcher.config.analysis.misc.images_wcs = get_images_header_path(self.config.path)
            self.launcher.config.analysis.misc.images_unit = "Jy"
            self.launcher.config.analysis.misc.images_kernels = None
            self.launcher.config.analysis.misc.rebin_wcs = None

        else:

            self.launcher.config.analysis.misc.fluxes = True       # calculate observed fluxes
            self.launcher.config.analysis.misc.images = False
        # observation_filters
        self.launcher.config.analysis.misc.observation_filters = self.observed_filter_names
        # observation_instruments
        self.launcher.config.analysis.misc.observation_instruments = [self.earth_instrument_name]

        # Set spectral convolution flag
        self.launcher.config.analysis.misc.spectral_convolution = self.fitting_run.fitting_configuration.spectral_convolution

        # Set the path to the modeling directory to the simulation object
        self.launcher.config.analysis.modeling_path = self.config.path

        # Set analyser classes
        if self.is_images_modeling: self.launcher.config.analysers = ["pts.modeling.fitting.modelanalyser.ImagesFitModelAnalyser"]
        else: self.launcher.config.analysers = ["pts.modeling.fitting.modelanalyser.SEDFitModelAnalyser"]

    # -----------------------------------------------------------------

    def get_initial_generation_name(self):

        """
        This function ...
        :return: 
        """

        return "initial"

    # -----------------------------------------------------------------

    def get_genetic_generation_name(self, index):

        """
        This function ...
        :param index: 
        :return: 
        """

        return str("Generation" + str(index))

    # -----------------------------------------------------------------

    def set_generator(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the model generator ...")

        # Generate new models based on a simple grid (linear or logarithmic) of parameter values
        if self.config.generation_method == "grid":

            # Set a name for the generation
            self.generation_name = time.unique_name("grid")

            # Create the model generator
            self.generator = GridModelGenerator()

        # Generate new models instinctively based on the current probability distribution of the parameters
        elif self.config.generation_method == "instinctive":

            # Set a name for the generation
            self.generation_name = time.unique_name("instinctive")

            # Create the model generator
            self.generator = InstinctiveModelGenerator()

        # Generate new models using genetic algorithms
        elif self.config.generation_method == "genetic":

            # Not the initial generation
            if self.get_initial_generation_name() in self.fitting_run.generation_names:

                print(self.fitting_run.last_genetic_generation_index)

                # Set index and name
                self.generation_index = self.fitting_run.last_genetic_generation_index + 1
                self.generation_name = self.get_genetic_generation_name(self.generation_index)

            # Initial generation
            else: self.generation_name = self.get_initial_generation_name()

            # Create the generator
            self.generator = GeneticModelGenerator()

        # Invalid generation method
        else: raise ValueError("Invalid generation method: " + str(self.config.generation_method))

        # Debugging
        log.debug("The name of the new generation of parameter exploration is '" + self.generation_name + "'")

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file template ...")

        # Load the labeled ski template file
        self.ski = self.fitting_run.ski_template

    # -----------------------------------------------------------------

    def set_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter ranges ...")

        # Create a definition
        definition = ConfigurationDefinition(write_config=False)

        extra_info = dict()

        # Check if there are any models that have been evaluated
        if self.fitting_run.has_evaluated_models:

            # Inform the user
            # log.info("Determining the parameter ranges based on the current best values and the specified relative ranges ...")

            # Get the best model
            model = self.fitting_run.best_model

            # Debugging
            # log.debug("Using the parameter values of simulation '" + model.simulation_name + "' of generation '" + model.generation_name + "' ...")

            # Get the parameter values of the best model
            parameter_values = model.parameter_values

            # Set info
            for label in parameter_values: extra_info[label] = "parameter value of current best model = " + stringify(parameter_values[label])[1]

        else:

            # Inform the user
            #log.info("Determining the parameter ranges based on the first guess values and the specified relative ranges ...")

            # Get the initial guess values
            parameter_values = self.fitting_run.first_guess_parameter_values

            # Set info
            for label in parameter_values:
                extra_info[label] = "initial parameter value = " + stringify(parameter_values[label])[1]

        # Loop over the free parameters
        for label in self.fitting_run.free_parameter_labels:

            # Skip if range is already defined for this label
            if label in self.ranges: continue

            # Get the default range
            default_range = self.fitting_run.fitting_configuration[label + "_range"]
            ptype, string = stringify_not_list(default_range)

            # Determine description
            description = "the range of " + label
            description += " (" + extra_info[label] + ")"

            # Add the optional range setting for this free parameter
            definition.add_optional(label + "_range", ptype, description, default_range)

        # Get the ranges
        if len(definition) > 0:

            setter = InteractiveConfigurationSetter("ranges", add_cwd=False, add_logging=False)
            config = setter.run(definition)

        # No parameters for which the ranges still have to be specified interactively
        else: config = None

        # Set the ranges
        for label in self.fitting_run.free_parameter_labels:

            # If range is already defined
            if label in self.ranges: continue

            # Set the range
            self.ranges[label] = config[label + "_range"]

    # -----------------------------------------------------------------

    def generate_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the model parameters ...")

        # Set generator options
        self.generator.config.ngenerations = self.config.ngenerations
        self.generator.config.nmodels = self.config.nsimulations

        # Run the model generator
        self.generator.run(fitting_run=self.fitting_run)

    # -----------------------------------------------------------------

    def set_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the generation info ...")

        # Get the previous wavelength grid level
        wavelength_grid_level = self.fitting_run.current_wavelength_grid_level
        #dust_grid_level = self.fitting_run.current_dust_grid_level

        # Determine the wavelength grid level
        #if self.config.refine_wavelengths:
        if self.config.refine_spectral:
            if wavelength_grid_level == self.fitting_run.highest_wavelength_grid_level: log.warning("Cannot refine wavelength grid: highest level reached (" + str(wavelength_grid_level) + ")")
            else: wavelength_grid_level += 1

        # Determine the dust grid level
        #if self.config.refine_dust:
        #    if dust_grid_level == self.fitting_run.highest_dust_grid_level: log.warning("Cannot refine dust grid: highest level reached (" + str(dust_grid_level) + ")")
        #    else: dust_grid_level += 1
        if self.config.refine_spatial: pass

        # Determine the number of photon packages
        if self.config.increase_npackages: npackages = int(self.fitting_run.current_npackages * self.config.npackages_factor)
        else: npackages = self.fitting_run.current_npackages

        # Determine whether selfabsorption should be enabled
        if self.config.selfabsorption is not None: selfabsorption = self.config.selfabsorption
        else: selfabsorption = self.fitting_run.current_selfabsorption

        # Determine whether transient heating should be enabled
        if self.config.transient_heating is not None: transient_heating = self.config.transient_heating
        else: transient_heating = self.fitting_run.current_transient_heating

        # Set the generation info
        self.generation.name = self.generation_name
        self.generation.index = self.generation_index
        self.generation.method = self.config.generation_method
        self.generation.wavelength_grid_level = wavelength_grid_level
        #self.generation.dust_grid_level = dust_grid_level
        self.generation.nsimulations = self.config.nsimulations
        self.generation.npackages = npackages
        self.generation.selfabsorption = selfabsorption
        self.generation.transient_heating = transient_heating

    # -----------------------------------------------------------------

    def set_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the input paths ...")

        # Initialize the SimulationInput object
        self.simulation_input = SimulationInput()

        # Set the paths to the input maps
        for name in self.fitting_run.input_map_paths:
            path = self.fitting_run.input_map_paths[name]
            self.simulation_input.add_file(path, name)

        # Determine and set the path to the appropriate wavelength grid file
        wavelength_grid_path = self.fitting_run.wavelength_grid_path_for_level(self.generation.wavelength_grid_level)
        #self.input_paths.append(wavelength_grid_path)
        self.simulation_input.add_file(wavelength_grid_path)

        # Get the number of wavelengths
        self.nwavelengths = len(WavelengthGrid.from_skirt_input(wavelength_grid_path))

        # Debugging
        log.debug("The wavelength grid for the simulations contains " + str(self.nwavelengths) + " wavelength points")

    # -----------------------------------------------------------------

    def create_generation_directory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the generation directory")

        # Determine the path to the generation directory
        self.generation.path = fs.create_directory_in(self.fitting_run.generations_path, self.generation_name)

        # Determine the path to the individuals table
        self.generation.individuals_table_path = fs.join(self.generation.path, "individuals.dat")

        # Determine the path to the generation parameters table
        self.generation.parameters_table_path = fs.join(self.generation.path, "parameters.dat")

        # Determine the path to the chi squared table
        self.generation.chi_squared_table_path = fs.join(self.generation.path, "chi_squared.dat")

        # Initialize the individuals table
        self.individuals_table = IndividualsTable()

        # Initialize the parameters table
        self.parameters_table = ParametersTable(parameters=self.fitting_run.free_parameter_labels, units=self.fitting_run.parameter_units)

        # Initialize the chi squared table
        self.chi_squared_table = ChiSquaredTable()

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski template for the properties of this generation ...")

        # Set packages
        self.set_npackages()

        # Set self-absoprtion
        self.set_selfabsorption()

        # Set transient heating
        self.set_transient_heating()

        # Set wavelength grid
        if self.fitting_run.has_wavelength_grids: self.set_wavelength_grid()

        # Set model representation
        self.set_representation()

    # -----------------------------------------------------------------

    def set_npackages(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the number of photon packages to " + str(self.generation.npackages) + " ...")

        # Set the number of photon packages per wavelength
        self.ski.setpackages(self.generation.npackages)

    # -----------------------------------------------------------------

    def set_selfabsorption(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Enabling dust self-absorption ..." if self.generation.selfabsorption else "Disabling dust self-absorption ...")

        # Set dust self-absorption
        if self.generation.selfabsorption: self.ski.enable_selfabsorption()
        else: self.ski.disable_selfabsorption()

    # -----------------------------------------------------------------

    def set_transient_heating(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Enabling transient heating ..." if self.generation.transient_heating else "Disabling transient heating ...")

        # Set transient heating
        if self.generation.transient_heating: self.ski.set_transient_dust_emissivity()
        else: self.ski.set_grey_body_dust_emissivity()

    # -----------------------------------------------------------------

    def set_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the name of the wavelengths file to " + fs.name(self.fitting_run.wavelength_grid_path_for_level(self.generation.wavelength_grid_level)) + " (level " + str(self.generation.wavelength_grid_level) + ") ...")

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid(fs.name(self.fitting_run.wavelength_grid_path_for_level(self.generation.wavelength_grid_level)))

    # -----------------------------------------------------------------

    def set_representation(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the model representation ...")

    # -----------------------------------------------------------------

    def set_parallelization_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme for local execution ...")

        # Get properties of the local machine
        nnodes = par.nnodes()
        nsockets = par.sockets_per_node()
        ncores = par.cores_per_socket()
        memory = par.virtual_memory().to("Gbyte")
        threads_per_core = par.nthreads_per_core()
        hyperthreading = threads_per_core > 1
        mpi = introspection.has_mpi()

        # Create the parallelization tool
        tool = ParallelizationTool()

        # Set configuration options
        tool.config.ski = self.ski
        tool.config.input = self.simulation_input

        # Set host properties
        tool.config.nnodes = nnodes
        tool.config.nsockets = nsockets
        tool.config.ncores = ncores
        tool.config.memory = memory

        # MPI available and used
        tool.config.mpi = mpi
        tool.config.hyperthreading = hyperthreading
        tool.config.threads_per_core = threads_per_core

        # Number of dust cells
        tool.config.ncells = None  # number of dust cells (relevant if ski file uses a tree dust grid)

        # Don't show
        tool.config.show = False

        # Run the parallelization tool
        tool.run()

        # Get the parallelization scheme
        parallelization = tool.parallelization

        # Debugging
        log.debug("The parallelization scheme for local execution is " + str(parallelization))

        # Set the parallelization scheme
        self.launcher.set_parallelization_for_local(parallelization)

    # -----------------------------------------------------------------

    def set_parallelization_remote(self):

        """
        This function sets the parallelization scheme for those remote hosts used by the batch launcher that use
        a scheduling system (the parallelization for the other hosts is left up to the batch launcher and will be
        based on the current load of the corresponding system).
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme for the remote host(s) that use a scheduling system ...")

        # Loop over the IDs of the hosts used by the batch launcher that use a scheduling system
        for host in self.launcher.scheduler_hosts:

            # Create the parallelization tool
            tool = ParallelizationTool()

            # Set configuration options
            tool.config.ski = self.ski
            tool.config.input = self.simulation_input

            # Set host properties
            tool.config.nnodes = self.config.nnodes
            tool.config.nsockets = host.cluster.sockets_per_node
            tool.config.ncores = host.cluster.cores_per_sockets
            tool.config.memory = host.cluster.memory

            # MPI available and used
            tool.config.mpi = True
            tool.config.hyperthreading = False # no hyperthreading
            #tool.config.threads_per_core = None

            # Number of dust cells
            tool.config.ncells = None # number of dust cells (relevant if ski file uses a tree dust grid)

            # Run the tool
            tool.run()

            # Get the parallelization
            parallelization = tool.parallelization

            # Debugging
            log.debug("Parallelization scheme for host " + host.id + ": " + str(parallelization))

            # Set the parallelization for this host
            self.launcher.set_parallelization_for_host(host.id, parallelization)

    # -----------------------------------------------------------------

    def estimate_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the runtimes based on the results of previously finished simulations ...")

        # Create a RuntimeEstimator instance
        estimator = RuntimeEstimator.from_file(self.fitting_run.timing_table_path)

        # Initialize a dictionary to contain the estimated walltimes for the different hosts with scheduling system
        walltimes = dict()

        # Loop over the hosts which use a scheduling system and estimate the walltime
        for host in self.launcher.scheduler_hosts:

            # Debugging
            log.debug("Estimating the runtime for host '" + host.id + "' ...")

            # Get the parallelization scheme that we have defined for this remote host
            parallelization = self.launcher.parallelization_for_host(host.id)

            # Visualisation of the distribution of estimated runtimes
            if self.config.visualise: plot_path = fs.join(self.visualisation_path, time.unique_name("explorer_runtimes_" + host.id) + ".pdf")
            else: plot_path = None

            # Estimate the runtime for the current number of photon packages and the current remote host
            runtime = estimator.runtime_for(self.ski, parallelization, host.id, host.cluster_name, self.config.data_parallel, nwavelengths=self.nwavelengths, plot_path=plot_path)

            # Debugging
            log.debug("The estimated runtime for this host is " + str(runtime) + " seconds")

            # Set the estimated walltime
            walltimes[host.id] = runtime

        # Create and set scheduling options for each host that uses a scheduling system
        for host_id in walltimes: self.scheduling_options[host_id] = SchedulingOptions.from_dict({"walltime": walltimes[host_id]})

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")

        # Set the paths to the directories to contain the launch scripts (job scripts) for the different remote hosts
        # Just use the directory created for the generation
        for host_id in self.launcher.host_ids: self.launcher.set_script_path(host_id, self.generation.path)

        # Enable screen output logging for remotes without a scheduling system for jobs
        for host_id in self.launcher.no_scheduler_host_ids: self.launcher.enable_screen_output(host_id)

        # Loop over the simulations, add them to the queue
        for simulation_name in self.simulation_names:

            # Get the parameter values
            parameter_values = self.parameters_table.parameter_values_for_simulation(simulation_name)

            # Prepare simulation directories, ski file, and return the simulation definition
            definition = prepare_simulation(simulation_name, self.ski, parameter_values, self.object_name, self.simulation_input, self.generation.path)

            # Debugging
            log.debug("Adding a simulation to the queue with:")
            log.debug("")
            log.debug(" - name: " + simulation_name)
            log.debug(" - input: " + str(self.simulation_input))
            log.debug(" - ski path: " + definition.ski_path)
            log.debug(" - output path: " + definition.output_path)
            log.debug("")

            # Put the parameters in the queue and get the simulation object
            self.launcher.add_to_queue(definition, simulation_name)

            # Set scheduling options (for the different remote hosts with a scheduling system)
            for host_id in self.scheduling_options: self.launcher.set_scheduling_options(host_id, simulation_name, self.scheduling_options[host_id])

        # Run the launcher, launches the simulations and retrieves and analyses finished simulations
        simulations = self.launcher.run()

        # Check the launched simulations
        self.check_simulations(simulations)

    # -----------------------------------------------------------------

    def check_simulations(self, simulations):

        """
        This function ...
        :param simulations: 
        :return: 
        """

        # Inform the user
        log.info("Checking the simulations ...")

        # Check the number of simulations that were effectively launched
        if self.nmodels == len(simulations):
            log.success("All simulations were scheduled succesfully")
            return

        # No simulations were launched
        if len(simulations) == 0:

            # Show error message
            log.error("No simulations could be launched: removing generation")
            log.error("Try again later")
            log.error("Cleaning up generation and quitting ...")

            # Remove this generation from the generations table
            self.fitting_run.generations_table.remove_entry(self.generation_name)
            self.fitting_run.generations_table.save()

            # Remove the generation directory
            fs.remove_directory(self.generation.path)

            # Quit
            exit()

        # Less simulations were launched
        elif len(simulations) < self.nmodels:

            # Get the names of simulations that were launched
            launched_simulation_names = [simulation.name for simulation in simulations]
            if None in launched_simulation_names: raise RuntimeError("Some or all simulation don't have a name defined")

            # Show error message
            log.error("Launching a simulation for the following models failed:")
            log.error("")

            # Loop over all simulations in the parameters table
            failed_indices = []
            for index, simulation_name in enumerate(self.parameters_table.simulation_names):

                # This simulation is OK
                if simulation_name in launched_simulation_names: continue

                log.error("Model #" + str(index + 1))
                log.error("")
                parameter_values = self.parameters_table.parameter_values_for_simulation(simulation_name)
                for label in parameter_values: log.error(" - " + label + ": " + stringify_not_list(parameter_values[label])[1])
                log.error("")

                failed_indices.append(index)

            # Show error message
            log.error("Removing corresponding entries from the model parameters table ...")

            # Remove rows and save
            self.parameters_table.remove_rows(failed_indices)
            self.parameters_table.save()

        # Unexpected
        else: raise RuntimeError("Unexpected error where nsmulations > nmodels")

    # -----------------------------------------------------------------

    @property
    def model_names(self):

        """
        This function ...
        :return: 
        """

        return self.generator.individual_names

    # -----------------------------------------------------------------

    @property
    def nmodels(self):

        """
        This function ...
        :return:
        """

        return self.generator.nmodels

    # -----------------------------------------------------------------

    @property
    def model_parameters(self):

        """
        This function ...
        :return: 
        """

        return self.generator.parameters

    # -----------------------------------------------------------------

    @property
    def uses_schedulers(self):

        """
        This function ...
        :return:
        """

        return self.launcher.uses_schedulers

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):

        """
        This function ...
        :return: 
        """

        return self.individuals_table.simulation_names

    # -----------------------------------------------------------------

    def fill_tables(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Filling the tables for the current generation ...")

        # Loop over the model names
        counter = 0
        for name in self.model_names:

            # Generate the simulation name
            simulation_name = generate_simulation_name()

            # Debugging
            log.debug("Adding an entry to the individuals table with:")
            log.debug("")
            log.debug(" - Simulation name: " + simulation_name)
            log.debug(" - Individual_name: " + name)
            log.debug("")

            # Add entry
            self.individuals_table.add_entry(simulation_name, name)

            # Get the parameter values
            parameter_values = get_parameter_values_for_named_individual(self.model_parameters, name, self.fitting_run)

            # Debugging
            log.debug("Adding entry to the parameters table with:")
            log.debug("")
            log.debug(" - Simulation name: " + simulation_name)
            for label in parameter_values: log.debug(" - " + label + ": " + stringify_not_list(parameter_values[label])[1])
            log.debug("")

            # Add an entry to the parameters table
            self.parameters_table.add_entry(simulation_name, parameter_values)

            # Increment counter
            counter += 1

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 1. Write the generation info
        self.write_generation()

        # 2. Write the parameters table
        self.write_parameters()

        # 3. Write the (empty) chi squared table
        self.write_chi_squared()

    # -----------------------------------------------------------------

    def write_generation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing generation info ...")

        # Add an entry to the generations table
        self.fitting_run.generations_table.add_entry(self.generation, self.ranges)

        # Save the table
        self.fitting_run.generations_table.save()

    # -----------------------------------------------------------------

    def write_individuals(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the individuals table ...")

        # Save the individuals table
        self.individuals_table.saveto(self.generation.individuals_table_path)

    # -----------------------------------------------------------------

    def write_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the model parameters table ...")

        # Save the parameters table
        self.parameters_table.saveto(self.generation.parameters_table_path)

    # -----------------------------------------------------------------

    def write_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the chi squared table ...")

        # Save the chi squared table
        self.chi_squared_table.saveto(self.generation.chi_squared_table_path)

# -----------------------------------------------------------------
