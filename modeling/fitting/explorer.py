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

# Import standard modules
import traceback
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.launch.batch import BatchLauncher
from .modelgenerators.grid import GridModelGenerator
from .modelgenerators.genetic import GeneticModelGenerator
from ...core.tools import time
from .tables import ParametersTable, ChiSquaredTable, IndividualsTable
from ...core.tools.stringify import stringify_not_list, stringify
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.remote.host import load_host
from ...core.basics.configuration import ConfigurationDefinition, create_configuration_interactive
from .evaluate import prepare_simulation, get_parameter_values_for_named_individual, make_test_definition
from ...core.simulation.input import SimulationInput
from .generation import GenerationInfo, Generation
from ...core.tools.stringify import tostr
from ...core.basics.configuration import prompt_proceed, prompt_integer
from ...core.prep.smile import SKIRTSmileSchema
from ...core.tools.utils import lazyproperty
from ...core.prep.deploy import Deployer
from ...core.basics.range import QuantityRange
from ...core.tools import formatting as fmt
from ...magic.plot.imagegrid import StandardImageGridPlotter

# -----------------------------------------------------------------

class ParameterExplorer(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ParameterExplorer, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The fitting run
        self.fitting_run = None

        # The ski template
        self.ski = None

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The generation info
        self.generation_info = GenerationInfo()

        # The generation object
        self.generation = None

        # The individuals table
        self.individuals_table = None

        # The parameters table
        self.parameters_table = None

        # The chi squared table
        self.chi_squared_table = None

        # The parameter ranges
        self.ranges = dict()

        # Fixed initial parameters
        self.fixed_initial_parameters = None

        # The generation index and name
        self.generation_index = None
        self.generation_name = None

        # The model generator
        self.generator = None

        # The simulation input
        self.simulation_input = None

        # Extra input for the model generator
        self.scales = None
        self.most_sampled_parameters = None
        self.sampling_weights = None

    # -----------------------------------------------------------------

    @property
    def needs_input(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.needs_input

    # -----------------------------------------------------------------

    @property
    def testing(self):

        """
        This function ...
        :return:
        """

        return self.config.test

    # -----------------------------------------------------------------

    @property
    def do_set_ranges(self):

        """
        This function ...
        :return:
        """

        return not self.has_all_ranges

    # -----------------------------------------------------------------

    @property
    def do_create_generation(self):

        """
        This function ...
        :return:
        """

        return not self.testing

    # -----------------------------------------------------------------

    @property
    def do_set_input(self):

        """
        This function ...
        :return:
        """

        return self.needs_input

    # -----------------------------------------------------------------

    @property
    def do_write(self):

        """
        This function ...
        :return:
        """

        return not self.testing

    # -----------------------------------------------------------------

    @property
    def do_plot(self):

        """
        This function ...
        :return:
        """

        return self.config.plot

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the ski template
        self.load_ski()

        # 3. Set the parameter ranges
        if self.do_set_ranges: self.set_ranges()

        # 4. Set the generation info
        self.set_info()

        # 5. Create the generation
        if self.do_create_generation: self.create_generation()

        # 6. Generate the model parameters
        self.generate_models()

        # 7. Set the paths to the input files
        if self.do_set_input: self.set_input()

        # 8. Adjust the ski template
        self.adjust_ski()

        # 9. Fill the tables for the current generation
        self.fill_tables()

        # 10. Writing
        if self.do_write: self.write()

        # 11. Show stuff
        self.show()

        # 12. Plot
        if self.do_plot: self.plot()

        # 13. Launch the simulations for different parameter values
        self.launch_or_finish()

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
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.free_parameter_labels

    # -----------------------------------------------------------------

    @property
    def has_all_ranges(self):

        """
        This function ...
        :return:
        """

        # Loop over the free parameter labels
        for label in self.parameter_labels:

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

        # Run locally?
        if self.config.local: self.config.remotes = []

        # Load the fitting run
        self.fitting_run = self.load_fitting_run(self.config.name)

        # Get ranges
        if "ranges" in kwargs: self.ranges = kwargs.pop("ranges")

        # Get the initial parameter values
        if "fixed_initial_parameters" in kwargs: self.fixed_initial_parameters = kwargs.pop("fixed_initial_parameters")

        # Set options for the batch launcher
        self.set_launcher_options()

        # Check for restarting generations
        if self.config.restart_from_generation is not None: self.clear_for_restart()

        # Set the model generator
        self.set_generator()

        # Check whether this is not the first generation so that we can use remotes with a scheduling system
        #if self.ngenerations == 0 and self.uses_schedulers:
        #    raise ValueError("The specified remote hosts cannot be used for the first generation: at least one remote uses a scheduling system")

        # Check whether the wavelength grids table is present
        if self.is_galaxy_modeling:
            if not fs.is_file(self.fitting_run.wavelength_grids_table_path): raise RuntimeError("Call initialize_fit_galaxy before starting the parameter exploration")

        # Get grid generation settings
        self.scales = kwargs.pop("scales", None)
        self.most_sampled_parameters = kwargs.pop("most_sampled_parameters", None)
        self.sampling_weights = kwargs.pop("sampling_weigths", None)

        # Deploy SKIRT
        if self.has_host_ids and self.config.deploy: self.deploy()

        # Initialize tables
        self.initialize_generation_tables()

    # -----------------------------------------------------------------

    @property
    def nprevious_generations(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.ngenerations

    # -----------------------------------------------------------------

    @property
    def previous_generation_names(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.generation_names

    # -----------------------------------------------------------------

    @property
    def last_previous_generation_name(self):

        """
        This function ...
        :return:
        """

        return self.previous_generation_names[-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def last_previous_generation(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.get_generation(self.last_previous_generation_name)

    # -----------------------------------------------------------------

    @property
    def first_generation(self):

        """
        This function ...
        :return:
        """

        return self.nprevious_generations == 0

    # -----------------------------------------------------------------

    @property
    def nprevious_genetic_generations(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.ngenetic_generations

    # -----------------------------------------------------------------

    @property
    def previous_genetic_generation_names(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.genetic_generations

    # -----------------------------------------------------------------

    @property
    def first_genetic_generation(self):

        """
        This function ...
        :return:
        """

        if not self.genetic_fitting: raise ValueError("Not a genetic generation")
        else: return self.nprevious_generations == 0

    # -----------------------------------------------------------------

    def get_description(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.fitting_run.parameter_descriptions[label]

    # -----------------------------------------------------------------

    @property
    def grid_settings(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.grid_config

    # -----------------------------------------------------------------

    @property
    def genetic_settings(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.genetic_config

    # -----------------------------------------------------------------

    def get_default_npoints(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        if self.grid_fitting: return self.grid_settings[label + "_npoints"]
        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def npoints(self):

        """
        This function ...
        :return:
        """

        # Prompt
        if self.config.prompt_npoints:

            # Check
            if self.config.npoints_all is not None: raise ValueError("Npoints is specified already through 'npoints_all'")

            # Prompt
            npoints_dict = dict()
            for label in self.parameter_labels: npoints_dict[label] = prompt_integer("npoints_" + label, "number of points for the " + self.get_description(label), default=self.get_default_npoints(label))
            return npoints_dict

        # Npoints all
        elif self.config.npoints_all is not None:
            npoints_dict = dict()
            for label in self.parameter_labels: npoints_dict[label] = self.config.npoints_all
            return npoints_dict

        # Npoints dict
        else: return self.config.npoints

    # -----------------------------------------------------------------

    def deploy(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Deploying SKIRT where necessary ...")

        # Create the deployer
        deployer = Deployer()

        # Don't deploy PTS
        deployer.config.skirt = True
        deployer.config.pts = False

        # Don't do anything locally
        deployer.config.local = False

        # Set the host ids
        deployer.config.hosts = self.remote_host_ids

        # Check versions between local and remote
        deployer.config.check = self.config.check_versions

        # Update PTS dependencies
        deployer.config.update_dependencies = self.config.update_dependencies

        # Do clean install
        deployer.config.clean = self.config.deploy_clean

        # Pubkey password
        deployer.config.pubkey_password = self.config.pubkey_password

        # Run the deployer
        deployer.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def remote_host_ids(self):

        """
        This function ...
        :return: 
        """

        # Set remote host IDs
        remote_host_ids = []
        if self.fitting_run.ngenerations == 0:
            for host_id in self.config.remotes:
                if load_host(host_id).scheduler:
                    log.warning("Not using remote host '" + host_id + "' for the initial generation because it uses a scheduling system for launching jobs")
                else: remote_host_ids.append(host_id)
        else: remote_host_ids = self.config.remotes

        # Return the host IDs
        return remote_host_ids

    # -----------------------------------------------------------------

    @lazyproperty
    def remote_hosts(self):

        """
        This function ...
        :return:
        """

        return [load_host(host_id) for host_id in self.remote_host_ids]

    # -----------------------------------------------------------------

    @lazyproperty
    def nhost_ids(self):

        """
        This function ...
        :return:
        """

        return len(self.remote_host_ids)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_host_ids(self):

        """
        This function ...
        :return:
        """

        return self.nhost_ids > 0

    # -----------------------------------------------------------------

    @property
    def modeling_config(self):

        """
        This function ...
        :return:
        """

        return self.environment.modeling_configuration

    # -----------------------------------------------------------------

    @property
    def other_host_ids(self):

        """
        This function ...
        :return:
        """

        if self.config.local_analysis: return []
        elif self.modeling_config.host_ids is None: return []
        else: return self.modeling_config.host_ids

    # -----------------------------------------------------------------

    @property
    def nother_host_ids(self):

        """
        This function ...
        :return:
        """

        return len(self.other_host_ids)

    # -----------------------------------------------------------------

    @property
    def has_other_host_ids(self):

        """
        This function ...
        :return:
        """

        return self.nother_host_ids > 0

    # -----------------------------------------------------------------

    @property
    def other_host_id(self):

        """
        This function ...
        :return:
        """

        if not self.has_other_host_ids: return None
        else: return self.other_host_ids[0]

    # -----------------------------------------------------------------

    @property
    def record_timing(self):

        """
        This function ...
        :return: 
        """

        if self.config.record_timing: return True
        elif len(self.remote_host_ids) > 0:
            log.warning("Record timing will be enabled because remote execution is used")
            return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def record_memory(self):

        """
        This function ...
        :return: 
        """

        if self.config.record_memory: return True
        elif len(self.remote_host_ids) > 0:
            log.warning("Record memory will be enabled because remote execution is used")
            return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def extract_timeline(self):

        """
        This
        :return: 
        """

        return self.record_timing or self.config.extract_timeline

    # -----------------------------------------------------------------

    @property
    def extract_memory(self):

        """
        This function ...
        :return: 
        """

        return self.record_memory or self.config.extract_memory

    # -----------------------------------------------------------------

    @property
    def reference_component_name(self):

        """
        This function ...
        :return:
        """

        return self.representation.reference_deprojection_name

    # -----------------------------------------------------------------

    @property
    def reference_map_path(self):

        """
        This function ...
        :return:
        """

        return self.representation.reference_map_path

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_wcs(self):

        """
        This function ...
        :return:
        """

        from ...magic.basics.coordinatesystem import CoordinateSystem
        if self.reference_map_path is None: return None
        else: return CoordinateSystem.from_file(self.reference_map_path)

    # -----------------------------------------------------------------

    @property
    def reference_wcs_path(self):

        """
        This function ...
        :return:
        """

        return self.reference_map_path

    # -----------------------------------------------------------------

    @property
    def timing_table_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.timing_table_path

    # -----------------------------------------------------------------

    @property
    def memory_table_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.memory_table_path

    # -----------------------------------------------------------------

    def set_launcher_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting options for the batch simulation launcher ...")

        # Basic launcher options
        self.set_basic_launcher_options()

        # Simulation options
        self.set_simulation_options()

        # Analysis options
        self.set_analysis_options()

    # -----------------------------------------------------------------

    def set_basic_launcher_options(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting basic launcher options ...")

        # Write assignments and leftover queues (if not testing)
        self.launcher.config.write = not self.testing

        # Simulations have approximately the same requirements
        self.launcher.config.same_requirements = True

        # Basic options
        self.launcher.config.shared_input = True                               # The input directories (or files) for the different simulations are shared
        self.launcher.config.remotes = self.remote_host_ids                         # The remote host(s) on which to run the simulations
        self.launcher.config.attached = self.config.attached                   # Run remote simulations in attached mode
        self.launcher.config.group_simulations = self.config.group             # Group multiple simulations into a single job (because a very large number of simulations will be scheduled) TODO: IMPLEMENT THIS
        self.launcher.config.group_walltime = self.config.group_walltime             # The preferred walltime for jobs of a group of simulations
        self.launcher.config.cores_per_process = self.config.cores_per_process # The number of cores per process, for non-schedulers
        self.launcher.config.dry = self.config.dry                             # Dry run (don't actually launch simulations, but allow them to be launched manually)
        self.launcher.config.progress_bar = True  # show progress bars for local execution
        self.launcher.config.keep = self.config.keep # keep remote input and output
        self.launcher.config.attached = self.config.attached  # run SKIRT in attached mode
        self.launcher.config.show_progress = self.config.show_progress

        # Memory and timing table paths (for recording and for estimating)
        self.launcher.config.timing_table_path = self.timing_table_path  # The path to the timing table file
        self.launcher.config.memory_table_path = self.memory_table_path  # The path to the memory table file

        # Record timing and memory?
        self.launcher.config.add_timing = self.record_timing
        self.launcher.config.add_memory = self.record_memory

        # Set runtimes plot path
        self.launcher.config.runtimes_plot_path = self.visualisation_path

        # Advanced parallelization options
        self.launcher.config.all_sockets = self.config.all_sockets
        self.launcher.config.nsockets = self.config.nsockets
        self.launcher.config.allow_multisocket_processes = self.config.allow_multisocket_processes

    # -----------------------------------------------------------------

    @property
    def nremotes(self):

        """
        This function ...
        :return:
        """

        return len(self.config.remotes)

    # -----------------------------------------------------------------

    @property
    def has_single_remote(self):

        """
        This function ...
        :return:
        """

        #return self.launcher.has_single_remote
        return self.nremotes == 1

    # -----------------------------------------------------------------

    @property
    def single_host_id(self):

        """
        Thisn function ...
        :return:
        """

        if not self.has_single_remote: raise RuntimeError("Not a single remote host")
        return self.config.remotes[0]

    # -----------------------------------------------------------------

    def set_simulation_options(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting simulation options ...")

        ## General
        self.launcher.config.relative = True

        ## Logging
        self.launcher.config.logging.verbose = True
        self.launcher.config.logging.memory = True

        # Set number of local processes
        self.launcher.set_nprocesses_local(self.config.nprocesses_local)

        # Set number of remote processes
        if self.config.nprocesses_remote is not None:
            if not self.has_single_remote: raise ValueError("Cannot specify number of remote processes when using more than one remote host")
            self.launcher.set_nprocesses_for_host(self.single_host_id, self.config.nprocesses_remote)

        # Set data parallel flag for local execution
        self.launcher.set_data_parallel_local(self.config.data_parallel_local)

        # Set data parallel flag for remote execution
        if self.config.data_parallel_remote is not None:
            if not self.has_single_remote: raise ValueError("Cannot set data parallelization flag when using more than one remote host")
            self.launcher.set_data_parallel_for_host(self.single_host_id, self.config.data_parallel_remote)

    # -----------------------------------------------------------------

    def set_analysis_options(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting analysis options ...")

        # General analysis options
        self.set_general_analysis_options()

        # Extraction options
        self.set_extraction_options()

        # Plotting options
        self.set_plotting_options()

        # Misc options
        self.set_misc_options()

    # -----------------------------------------------------------------

    def set_general_analysis_options(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting general analysis options ...")

        # To create the extr, plot, misc directories relative in the simulation directory
        self.launcher.config.relative_analysis_paths = True

        # Set the path to the modeling directory
        self.launcher.config.analysis.modeling_path = self.config.path

        # Set analyser classes
        if self.is_images_modeling: self.launcher.config.analysers = ["pts.modeling.fitting.modelanalyser.ImagesFitModelAnalyser"]
        else: self.launcher.config.analysers = ["pts.modeling.fitting.modelanalyser.SEDFitModelAnalyser"]

    # -----------------------------------------------------------------

    def set_extraction_options(self):

        """
        This function ....
        :return:
        """

        # Debugging
        log.debug("Setting extraction options ...")

        # Extraction
        self.launcher.config.analysis.extraction.path = "extr"    # name of the extraction directory
        self.launcher.config.analysis.extraction.progress = self.config.extract_progress  # extract progress information
        self.launcher.config.analysis.extraction.timeline = self.extract_timeline # extract the simulation timeline
        self.launcher.config.analysis.extraction.memory = self.extract_memory    # extract memory information

    # -----------------------------------------------------------------

    @property
    def truncated_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.truncated_sed_path

    # -----------------------------------------------------------------

    @property
    def observed_sed_paths(self):

        """
        This function ...
        :return:
        """

        if self.is_galaxy_modeling: return [self.observed_sed_path, self.truncated_sed_path]
        elif self.is_sed_modeling: return [self.observed_sed_path]
        else: raise ValueError("Observed SED not defined for modeling types other than 'galaxy' or 'sed'")

    # -----------------------------------------------------------------

    def set_plotting_options(self):

        """
        Thisn function ...
        :return:
        """

        # Debugging
        log.debug("Setting plotting options ...")

        # Plotting
        self.launcher.config.analysis.plotting.path = "plot"  # name of the plot directory
        self.launcher.config.analysis.plotting.seds = self.config.plot_seds    # Plot the output SEDs
        self.launcher.config.analysis.plotting.reference_seds = self.observed_sed_paths  # the path to the reference SED (for plotting the simulated SED against the reference points)
        self.launcher.config.analysis.plotting.format = self.config.plotting_format  # plot format
        self.launcher.config.analysis.plotting.progress = self.config.plot_progress
        self.launcher.config.analysis.plotting.timeline = self.config.plot_timeline
        self.launcher.config.analysis.plotting.memory = self.config.plot_memory

    # -----------------------------------------------------------------

    @lazyproperty
    def photometry_image_paths_for_fitting_filter_names(self):

        """
        This function ...
        :return:
        """

        # Create new dictionary
        paths = OrderedDict()

        # Loop over the paths
        for filter_name in self.environment.photometry_image_paths_for_filter_names:

            # Skip filters that are not in the fitting filters list
            if filter_name not in self.fitting_filter_names: continue

            # Add the path
            path = self.environment.photometry_image_paths_for_filter_names[filter_name]
            paths[filter_name] = path

        # Return
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def fit_sed_path(self):

        """
        This function ...
        :return:
        """

        if self.config.fit_not_clipped: return self.truncated_sed_path
        else: return self.observed_sed_path

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_sed_paths(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        paths = OrderedDict()

        # Set the SED paths
        paths["Observed clipped fluxes"] = self.observed_sed_path
        paths["Observed truncated fluxes"] = self.truncated_sed_path

        # Return the dictionary
        return paths

    # -----------------------------------------------------------------

    def set_misc_options(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting miscellaneous analysis options ...")

        # Miscellaneous
        self.launcher.config.analysis.misc.path = "misc"       # name of the misc output directory

        # From images
        if self.use_images:

            self.launcher.config.analysis.misc.fluxes_from_images = True # calculate observed fluxes from images
            self.launcher.config.analysis.misc.fluxes = False
            self.launcher.config.analysis.misc.images = False

            # Set instrument and coordinate system path
            self.launcher.config.analysis.misc.fluxes_from_images_instrument = self.earth_instrument_name
            self.launcher.config.analysis.misc.fluxes_from_images_wcs = self.reference_wcs_path

            # Set mask paths
            self.launcher.config.analysis.misc.fluxes_from_images_masks = self.photometry_image_paths_for_fitting_filter_names
            self.launcher.config.analysis.misc.fluxes_from_images_mask_from_nans = True

            # Write the fluxes images
            self.launcher.config.analysis.misc.write_fluxes_images = True

            # Plot fluxes
            self.launcher.config.analysis.misc.plot_fluxes_from_images = True
            self.launcher.config.analysis.misc.plot_fluxes_from_images_reference_seds = self.plot_sed_paths #self.fit_sed_path

            # Set remote for creating images from datacubes
            self.launcher.config.analysis.misc.fluxes_from_images_remote = self.other_host_id
            self.launcher.config.analysis.misc.fluxes_from_images_remote_spectral_convolution = True
            #self.launcher.config.analysis.misc.fluxes_from_images_remote_threshold =
            #self.launcher.config.analysis.misc.fluxes_from_images_remote_npixels_threshold =
            #self.launcher.config.analysis.misc.fluxes_from_images_rebin_remote_threshold =

            # Make a plot of the images
            self.launcher.config.analysis.misc.plot_fluxes_images = True

        # From SEDs
        else:

            self.launcher.config.analysis.misc.fluxes_from_images = False
            self.launcher.config.analysis.misc.fluxes = True       # calculate observed fluxes from SEDs
            self.launcher.config.analysis.misc.images = False

            # Plot fluxes
            self.launcher.config.analysis.misc.plot_fluxes = True
            self.launcher.config.analysis.misc.plot_fluxes_reference_seds = self.plot_sed_paths #self.fit_sed_path

        # Observation filters and observation instruments
        self.launcher.config.analysis.misc.observation_filters = self.fitting_filter_names
        self.launcher.config.analysis.misc.observation_instruments = [self.earth_instrument_name]

        # Set spectral convolution flag
        self.launcher.config.analysis.misc.fluxes_spectral_convolution = self.spectral_convolution
        self.launcher.config.analysis.misc.fluxes_from_images_spectral_convolution = self.spectral_convolution

    # -----------------------------------------------------------------

    def get_initial_generation_name(self):

        """
        This function ...
        :return: 
        """

        return self.fitting_run.get_initial_generation_name()

    # -----------------------------------------------------------------

    def get_genetic_generation_name(self, index):

        """
        This function ...
        :param index: 
        :return: 
        """

        return self.fitting_run.get_genetic_generation_name(index)

    # -----------------------------------------------------------------

    def clear_for_restart(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing things for restarting from generation '" + self.config.restart_from_generation + "' ...")

        # Get the gneration names
        to_clear = self.get_to_clear_generation_names()

        # Get the generations table
        generations_table = self.fitting_run.generations_table
        best_parameters_table = self.fitting_run.best_parameters_table

        # Get the names of the original genetic generations
        #original_genetic_generation_names = generations_table.genetic_generations_with_initial
        original_genetic_generations_with_initial_names_and_indices = generations_table.genetic_generations_with_initial_names_and_indices

        # Keep track of the lowest genetic generation index
        lowest_genetic_generation_index = None
        removed_initial = False

        to_be_removed_paths = []

        # Loop over the generations to be cleared
        for generation_name in to_clear:

            # Prompt to proceed
            if prompt_proceed("Are you absolutely sure all output from generation '" + generation_name + "' can be removed?"):

                # Update the lowest genetic generation index
                if generation_name.startswith("Generation"):
                    index = self.fitting_run.index_for_generation(generation_name)
                    if lowest_genetic_generation_index is None or index < lowest_genetic_generation_index: lowest_genetic_generation_index = index

                if generation_name == "initial": removed_initial = True

                # Remove from generations table
                generations_table.remove_entry(generation_name)

                # Remove from best_parameters table
                best_parameters_table.remove_entry(generation_name)

                # Remove from prob/generations
                prob_generations_path = fs.create_directory_in(self.fitting_run.prob_path, "generations")
                prob_generation_path = fs.join(prob_generations_path, generation_name + ".dat")
                #fs.remove_file(prob_generation_path)
                to_be_removed_paths.append(prob_generation_path)

                # Remove directory from generations/
                generation_directory_path = self.fitting_run.get_generation_path(generation_name)
                #fs.remove_directory(generation_directory_path)
                to_be_removed_paths.append(generation_directory_path)

            # User doesn't want to proceed
            else:

                # Exit with an error
                log.error("Cannot proceed without confirmation")
                exit()

        # IF GENETIC GENERATIONS ARE CLEARED, REPLACE THE MAIN ENGINE, MAIN PRNG AND MAIN OPTIMIZER.CFG
        if removed_initial:

            # Remove
            fs.remove_file(self.fitting_run.main_engine_path)
            fs.remove_file(self.fitting_run.main_prng_path)
            fs.remove_file(self.fitting_run.optimizer_config_path)

        # Some genetic generations are cleared, starting with some lowest genetic generation index
        elif lowest_genetic_generation_index is not None:

            # Search for the last remaining generation
            last_remaining_generation = None

            # Determine name of generation just before this index
            for other_name, other_index in original_genetic_generations_with_initial_names_and_indices:
                if other_index == lowest_genetic_generation_index - 1:
                    last_remaining_generation = other_name
                    break

            if last_remaining_generation: raise RuntimeError("Something went wrong")

            # Determine the path of this generation
            generation = self.fitting_run.get_generation(last_remaining_generation)

            # Determine the paths of the engine, prng and optimizer config
            engine_path = generation.engine_path
            prng_path = generation.prng_path
            optimizer_config_path = generation.optimizer_config_path

            # Replace the main engine, prng and optimizer config
            fs.replace_file(self.fitting_run.main_engine_path, engine_path)
            fs.replace_file(self.fitting_run.main_prng_path, prng_path)
            fs.replace_file(self.fitting_run.optimizer_config_path, optimizer_config_path)

        # Remove everything belonging the cleared generations
        fs.remove_directories_and_files(to_be_removed_paths)

        # Save the generations table
        generations_table.save()

        # Save the best parameters table
        best_parameters_table.save()

    # -----------------------------------------------------------------

    def get_to_clear_generation_names(self):

        """
        This function ...
        :return:
        """

        generation_name = self.config.restart_from_generation

        # Check whether the generation exists
        if generation_name not in self.fitting_run.generation_names: raise ValueError("Generation '" + generation_name + "' does not exist")

        # Generation names to clear
        to_clear = []

        # Grid-type generation
        if "grid" in generation_name:

            # Add to be cleared
            to_clear.append(generation_name)

            # Get the timestamp
            generation_time = time.get_time_from_unique_name(generation_name)

            # Loop over other 'grid' generations
            for other_generation_name in self.fitting_run.grid_generations:

                if other_generation_name == generation_name: continue

                # Get time
                other_generation_time = time.get_time_from_unique_name(other_generation_name)

                # If the time is later, add to generation names to be cleared
                if other_generation_time > generation_time: to_clear.append(generation_name)

        # Initial genetic generation
        elif generation_name == self.get_initial_generation_name():

            # All genetic generations have to be cleared
            to_clear = self.fitting_run.genetic_generations

        # Other genetic generation
        elif generation_name.startswith("Generation"):

            # Add to be cleared
            to_clear.append(generation_name)

            # Get the index of the generation
            index = self.fitting_run.index_for_generation(generation_name)

            # Loop over the other genetic generations
            for other_generation_name in self.fitting_run.genetic_generations:

                if other_generation_name == generation_name: continue

                # Get index of other
                other_index = self.fitting_run.index_for_generation(other_generation_name)

                # If the index is higher, add to be cleared
                if other_index > index: to_clear.append(other_generation_name)

        # Could not understand
        else: raise ValueError("Could not understand the nature of generation '" + generation_name + "'")

        # Retunr the list of generation names to clear
        return to_clear

    # -----------------------------------------------------------------

    @property
    def grid_fitting(self):

        """
        This function ...
        :return:
        """

        return self.config.generation_method == "grid"

    # -----------------------------------------------------------------

    @property
    def genetic_fitting(self):

        """
        This function ...
        :return:
        """

        return self.config.generation_method == "genetic"

    # -----------------------------------------------------------------

    def set_generator(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the model generator ...")

        # Generate new models based on a simple grid (linear or logarithmic) of parameter values
        if self.grid_fitting: self.set_grid_generator()

        # Generate new models using genetic algorithms
        elif self.genetic_fitting: self.set_genetic_generator()

        # Invalid generation method
        else: raise ValueError("Invalid generation method: " + str(self.config.generation_method))

        # Set general options for the model generator
        self.set_generator_options()

        # Debugging
        log.debug("The name of the new generation of parameter exploration is '" + self.generation_name + "'")

    # -----------------------------------------------------------------

    def set_grid_generator(self):

        """
        This function ...
        :param self: 
        :return: 
        """

        # Inform the user
        log.info("Setting grid model generator ...")

        # Set a name for the generation
        #self.generation_name = time.unique_name("grid")

        # Determine grid generation index
        highest_index = self.fitting_run.highest_grid_generation_index
        if highest_index is None: generation_index = 0
        else: generation_index = highest_index + 1

        # Set generation name
        self.generation_name = "grid_" + str(generation_index)

        # Create the model generator
        self.generator = GridModelGenerator()

    # -----------------------------------------------------------------

    def set_genetic_generator(self):

        """
        This function ...
        :param self: 
        :return: 
        """

        # Inform the user
        log.info("Setting genetic model generator ...")

        # Not the initial generation
        if self.get_initial_generation_name() in self.fitting_run.generation_names:

            # Set index and name
            self.generation_index = self.fitting_run.last_genetic_generation_index + 1
            self.generation_name = self.get_genetic_generation_name(self.generation_index)

        # Initial generation
        else: self.generation_name = self.get_initial_generation_name()

        # Create the generator
        self.generator = GeneticModelGenerator()

        # Set recurrence settings
        self.generator.config.check_recurrence = self.config.check_recurrence
        self.generator.config.recurrence_rtol = self.config.recurrence_rtol
        self.generator.config.recurrence_atol = self.config.recurrence_atol

    # -----------------------------------------------------------------

    def set_generator_options(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting general model generator options ...")

        # Set the modeling path for the model generator
        self.generator.config.path = self.config.path

        # Set generator options
        self.generator.config.ngenerations = self.config.ngenerations # only useful for genetic model generator (and then again, cannot be more then 1 yet)
        self.generator.config.nmodels = self.config.nsimulations

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

        # Automatic: based on previous generation(s)
        if self.config.auto_ranges and not self.first_generation: self.determine_ranges()

        # Prompt
        elif self.config.prompt_ranges: self.prompt_ranges()

        # Default
        else: self.set_default_ranges()

    # -----------------------------------------------------------------

    def determine_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the parameter ranges automatically ...")

        # Loop over the free parameters
        for label in self.parameter_labels:

            # If range is already defined
            if label in self.ranges: continue

            # Get the parameter distribution
            if self.fitting_run.has_distribution(label): distribution = self.fitting_run.get_parameter_distribution(label)
            else:
                log.warning("Global parameter distribution for the '" + label + "' parameter not found: using parameter distribution for generation '" + self.last_previous_generation_name + "' ...")
                distribution = self.fitting_run.get_parameter_distribution_for_generation(label, self.last_previous_generation_name)

            # Get the parameter unit
            unit = self.fitting_run.parameter_units[label]

            # Get leading values
            values, min_value, max_value, total_fraction = distribution.get_leading_values_and_edges(self.config.range_probability, logscale=True, return_fraction=True)

            # Add units if necessary
            if not hasattr(min_value, "unit"): min_value = min_value * unit
            if not hasattr(max_value, "unit"): max_value = max_value * unit

            # Debugging
            log.debug("")
            log.debug("Parameter '" + label + "':")
            log.debug(" - Most probable value: " + tostr(distribution.most_frequent))
            log.debug(" - Least probable value: " + tostr(distribution.least_frequent))
            log.debug(" - Number of unique values: " + tostr(distribution.nvalues))
            log.debug(" - Value(s) with (combined) >= " + str(self.config.range_probability * 100) + "% of the probablity: " + tostr(values))
            log.debug(" - Combined probability: " + tostr(total_fraction * 100, round=True, ndigits=5) + "%")
            log.debug(" - Old minimum value: " + tostr(distribution.min_value))
            log.debug(" - Old maximum value: " + tostr(distribution.max_value))
            log.debug(" - New minimum value: " + tostr(min_value))
            log.debug(" - New maximum value: " + tostr(max_value))

            # Set the range for this parameter
            self.ranges[label] = QuantityRange(min_value, max_value)

        # One more line
        log.debug("")

    # -----------------------------------------------------------------

    def prompt_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the parameter ranges ...")

        # Create definition
        definition = self.create_parameter_ranges_definition()

        # Get the ranges
        if len(definition) > 0: config = create_configuration_interactive(definition, "ranges", "parameter ranges", add_cwd=False, add_logging=False, prompt_optional=True)

        # No parameters for which the ranges still have to be specified interactively
        else: config = None

        # Set the ranges
        for label in self.parameter_labels:

            # If range is already defined
            if label in self.ranges: continue

            # Set the range
            self.ranges[label] = config[label + "_range"]

    # -----------------------------------------------------------------

    def set_default_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Using the default parameter ranges for each parameter ...")

        # Loop over the free parameters, add a setting slot for each parameter range
        for label in self.parameter_labels:

            # Skip if range is already defined for this label
            if label in self.ranges: continue

            # Get the default range
            default_range = self.fitting_run.fitting_configuration[label + "_range"]

            # Set the range
            self.ranges[label] = default_range

    # -----------------------------------------------------------------

    def create_parameter_ranges_definition(self):

        """
        This function ...
        :return: 
        """

        # Create a definition
        definition = ConfigurationDefinition(write_config=False)

        # Create info
        extra_info = self.create_parameter_ranges_info()

        # Loop over the free parameters, add a setting slot for each parameter range
        for label in self.parameter_labels:

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

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def create_parameter_ranges_info(self):

        """
        This function ...
        :return: 
        """

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
            for label in parameter_values: extra_info[label] = "initial parameter value = " + stringify(parameter_values[label])[1]

        # Return the info
        return extra_info

    # -----------------------------------------------------------------

    def generate_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the model parameters ...")

        # Run the model generator
        self.generator.run(fitting_run=self.fitting_run, parameter_ranges=self.ranges,
                           fixed_initial_parameters=self.fixed_initial_parameters, generation=self.generation,
                           scales=self.parameter_scales, most_sampled_parameters=self.most_sampled_parameters,
                           sampling_weights=self.sampling_weights, npoints=self.npoints)

        # Set the actual number of simulations for this generation
        self.generation_info.nsimulations = self.nmodels

    # -----------------------------------------------------------------

    @property
    def parameter_scales(self):

        """
        This function ...
        :return:
        """

        scales = dict()

        # Loop over the free parameter
        for label in self.fitting_run.free_parameter_labels:

            # Check whether scales were given as input
            if self.scales is not None and label in self.scales: scales[label] = self.scales[label]
            #elif self.config.scales is not None and label in self.config.scales: scales[label] = self.config.scales[label]
            #else: #raise ValueError("Scale was not set for '" + label + "'")
            #    # Take from grid fitting configuration
            #    scales[label] = self.fitting_run.grid_settings[label + "_scale"]
            else: scales[label] = self.fitting_run.grid_settings[label + "_scale"]

        # Return the scales
        return scales

    # -----------------------------------------------------------------

    @lazyproperty
    def selfabsorption(self):

        """
        This function ...
        :return:
        """

        # Determine whether selfabsorption should be enabled
        if self.config.selfabsorption is not None: return self.config.selfabsorption
        else: return self.fitting_run.current_selfabsorption

    # -----------------------------------------------------------------

    @lazyproperty
    def transient_heating(self):

        """
        This function ...
        :return:
        """

        # Determine whether transient heating should be enabled
        if self.config.transient_heating is not None: return self.config.transient_heating
        else: return self.fitting_run.current_transient_heating

    # -----------------------------------------------------------------

    @lazyproperty
    def spectral_convolution(self):

        """
        This function ...
        :return:
        """

        if self.config.spectral_convolution is not None: return self.config.spectral_convolution
        else: return self.fitting_run.current_spectral_convolution

    # -----------------------------------------------------------------

    @lazyproperty
    def use_images(self):

        """
        This function ...
        :return:
        """

        if self.config.use_images is not None: return self.config.use_images
        else: return self.fitting_run.current_use_images

    # -----------------------------------------------------------------

    @lazyproperty
    def npackages(self):

        """
        This function ...
        :return:
        """

        # Defined?
        if self.config.npackages is not None: npackages = self.config.npackages
        else:

            # Determine the number of photon packages from previous number
            if self.config.increase_npackages: npackages = int(self.fitting_run.current_npackages * self.config.npackages_factor)
            else: npackages = self.fitting_run.current_npackages

        # Check
        if npackages < self.ndust_cells:

            if self.config.adjust_npackages:
                log.debug("Adjusting the number of photon packages from " + str(npackages) + " to the number of dust cells (" + str(self.ndust_cells) + ") ...")
                npackages = int(self.ndust_cells * self.config.ncells_npackages_factor)
            else: log.warning("The number of photon packages (" + str(npackages) + ") is less than the number of dust cells (" + str(self.ndust_cells) + ")")

        # Return the number of photon packages
        return npackages

    # -----------------------------------------------------------------

    @lazyproperty
    def representation(self):

        """
        This function ...
        :return:
        """

        # DETERMINE THE REPRESENTATION
        if self.config.refine_spatial: return self.fitting_run.next_model_representation  # GET NEXT REPRESENTATION (THEY ARE NAMED IN ORDER OF SPATIAL RESOLUTION)

        # Get the previous (current because this generation is just
        else: return self.fitting_run.current_model_representation  # GET LAST REPRESENTATION #self.fitting_run.initial_representation

    # -----------------------------------------------------------------

    @property
    def representation_name(self):

        """
        This function ...
        :return:
        """

        return self.representation.name

    # -----------------------------------------------------------------

    @property
    def ndust_cells(self):

        """
        This function ...
        :return:
        """

        return self.representation.ndust_cells

    # -----------------------------------------------------------------

    @property
    def wavelength_grids_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.wavelength_grids_path

    # -----------------------------------------------------------------

    @property
    def highres_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        if self.config.highres is not None: return self.config.highres
        else: return self.fitting_run.is_highres_current_wavelength_grid

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid_name(self):

        """
        This function ...
        :return:
        """

        # Use high-resolution grids
        if self.highres_wavelength_grid:

            # Target number of wavelengths is defined
            if self.config.nwavelengths is not None: return "highres_" + str(self.config.nwavelengths)

            # Refine from previous high-resolution grid?
            elif self.fitting_run.is_highres_current_wavelength_grid:

                if self.config.refine_spectral: return self.fitting_run.next_wavelength_grid_name
                else: return self.fitting_run.current_wavelength_grid_name

            # Get lowest npoints wavelength grid of high-resolution grids
            else:
                if self.config.refine_spectral: log.warning("Not refining more: previous generation did not use high-resolution wavelength grid")
                return self.fitting_run.lowest_highres_wavelength_grid_name

        # Spectral convolution: use refined grids
        elif self.spectral_convolution:

            # Target number of wavelengths is defined
            if self.config.nwavelengths is not None: return "refined_" + str(self.config.nwavelengths)

            # Refine from previous refined grid?
            if self.fitting_run.is_refined_current_wavelength_grid:

                if self.config.refine_spectral: return self.fitting_run.next_wavelength_grid_name
                else: return self.fitting_run.current_wavelength_grid_name

            # Get lowest npoints wavelength grid of refined grids
            else:
                if self.config.refine_spectral: log.warning("Not refining more: previous generation did not use refined wavelength grid")
                return self.fitting_run.lowest_refined_wavelength_grid_name

        # Basic grids
        else:

            # Target number of wavelengths is defined
            if self.config.nwavelengths is not None: return "basic_" + str(self.config.nwavelengths)

            # Refine from previous basic grid?
            if self.fitting_run.is_basic_current_wavelength_grid:

                if self.config.refine_spectral: return self.fitting_run.next_wavelength_grid_name
                else: return self.fitting_run.current_wavelength_grid_name

            # Get lowest npoints wavelength grid of basic grids
            else:
                if self.config.refine_spectral: raise ValueError("Cannot refine: use 'nwavelengths' to define a specific wavelength grid (and 'highres' to control whether to use high-resolution grid)")
                return self.fitting_run.lowest_basic_wavelength_grid_name

    # -----------------------------------------------------------------

    @property
    def wavelength_grid_filename(self):

        """
        This function ...
        :return:
        """

        return self.wavelength_grid_name + ".dat"

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.wavelength_grids_path, self.wavelength_grid_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return WavelengthGrid.from_skirt_input(self.wavelength_grid_path)

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        return len(self.wavelength_grid)

    # -----------------------------------------------------------------

    @lazyproperty
    def fit_not_clipped(self):

        """
        This function ...
        :return:
        """

        if self.use_images and self.config.fit_not_clipped: raise ValueError("Cannot fit to non-clipped fluxes when using images for calculating observed fluxes (clip masks are used)")
        else: return self.config.fit_not_clipped

    # -----------------------------------------------------------------

    def set_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the generation info ...")

        # Set the generation info
        self.generation_info.name = self.generation_name
        self.generation_info.index = self.generation_index
        self.generation_info.method = self.config.generation_method

        # Wavelength grid and representation
        self.generation_info.wavelength_grid_name = self.wavelength_grid_name
        self.generation_info.model_representation_name = self.representation_name

        # DON'T DO IT HERE YET, GET THE NUMBER OF ACTUAL MODELS SPITTED OUT BY THE MODELGENERATOR (RECURRENCE)
        #self.generation.nsimulations = self.config.nsimulations

        # Set number of photon packages
        self.generation_info.npackages = self.npackages

        # Simulation options
        self.generation_info.selfabsorption = self.selfabsorption
        self.generation_info.transient_heating = self.transient_heating
        self.generation_info.spectral_convolution = self.spectral_convolution
        self.generation_info.use_images = self.use_images

        # Fit options
        self.generation_info.fit_not_clipped = self.fit_not_clipped

    # -----------------------------------------------------------------

    @lazyproperty
    def use_file_tree_dust_grid(self):

        """
        This function ...
        :return:
        """

        smile = SKIRTSmileSchema()
        return smile.supports_file_tree_grids and self.representation.has_dust_grid_tree

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

        # DETERMINE AND SET THE PATH TO THE APPROPRIATE DUST GRID TREE FILE
        if self.use_file_tree_dust_grid: self.simulation_input.add_file(self.representation.dust_grid_tree_path)

        # Determine and set the path to the appropriate wavelength grid file
        self.simulation_input.add_file(self.wavelength_grid_path)

        # Debugging
        log.debug("The wavelength grid for the simulations contains " + str(self.nwavelengths) + " wavelength points")

    # -----------------------------------------------------------------

    def create_generation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the generation directory")

        # Set the path to the generation directory
        self.generation_info.path = fs.create_directory_in(self.fitting_run.generations_path, self.generation_name)

        # Create the generation object
        self.generation = Generation(self.generation_info)

        # Set working directory for the batch launcher
        self.launcher.config.path = self.generation_path

    # -----------------------------------------------------------------

    def initialize_generation_tables(self):

        """
        This function ...
        :return: 
        """

        # Debugging
        log.debug("Initializing generation tables ...")

        # Initialize the individuals table
        self.individuals_table = IndividualsTable()

        # Initialize the parameters table
        self.parameters_table = ParametersTable(parameters=self.parameter_labels, units=self.fitting_run.parameter_units)

        # Initialize the chi squared table
        self.chi_squared_table = ChiSquaredTable()

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_instrument(self):

        """
        This function ...
        :return:
        """

        if self.use_images: return self.representation.simple_instrument
        else: return self.representation.sed_instrument

    # -----------------------------------------------------------------

    @lazyproperty
    def instruments(self):

        """
        This function ...
        :return:
        """

        instrs = dict()
        instrs[self.earth_instrument_name] = self.earth_instrument
        return instrs

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

        # Set instruments
        self.set_instruments()

    # -----------------------------------------------------------------

    def set_npackages(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the number of photon packages to " + str(self.generation_info.npackages) + " ...")

        # Set the number of photon packages per wavelength
        self.ski.setpackages(self.generation_info.npackages)

    # -----------------------------------------------------------------

    def set_selfabsorption(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Enabling dust self-absorption ..." if self.generation_info.selfabsorption else "Disabling dust self-absorption ...")

        # Set dust self-absorption
        if self.generation_info.selfabsorption: self.ski.enable_selfabsorption()
        else: self.ski.disable_selfabsorption()

    # -----------------------------------------------------------------

    def set_transient_heating(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Enabling transient heating ..." if self.generation_info.transient_heating else "Disabling transient heating ...")

        # Set transient heating
        if self.generation_info.transient_heating: self.ski.set_transient_dust_emissivity()
        else: self.ski.set_grey_body_dust_emissivity()

    # -----------------------------------------------------------------

    def set_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the name of the wavelengths file to '" + self.wavelength_grid_filename + "' ...")

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid(self.wavelength_grid_filename)

    # -----------------------------------------------------------------

    def set_representation(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the model representation ...")

        # GET DUST GRID
        if self.use_file_tree_dust_grid:

            # Get the file tree dust grid object
            dust_grid = self.representation.create_file_tree_dust_grid(write=False)

            # Make sure it is only the file name, not a complete path
            dust_grid.filename = fs.name(dust_grid.filename)

        # REGULAR DUST GRID OBJECT
        else: dust_grid = self.representation.dust_grid

        # Set the dust grid
        self.ski.set_dust_grid(dust_grid)

    # -----------------------------------------------------------------

    def set_instruments(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the instruments ...")

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add the instrument
        self.ski.add_instrument(self.earth_instrument_name, self.earth_instrument)

    # -----------------------------------------------------------------

    def launch_or_finish(self):

        """
        This function ...
        :return:
        """

        # Test whether simulations are required, because if the optimizer detects recurrence of earlier models,
        # it is possible that no simulations have to be done

        # Launch simulations
        if self.needs_simulations: self.launch()

        # No simulations need to be launched
        else: self.set_finishing_time()

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")

        # Only if not testing, because otherwise generation path does not exist
        if not self.testing:
            # Set the paths to the directories to contain the launch scripts (job scripts) for the different remote hosts
            # Just use the directory created for the generation
            for host_id in self.launcher.host_ids: self.launcher.set_script_path(host_id, self.generation_info.path)

        # Enable screen output logging for remotes without a scheduling system for jobs
        for host_id in self.launcher.no_scheduler_host_ids: self.launcher.enable_screen_output(host_id)

        # Loop over the simulations, add them to the queue
        for simulation_name in self.simulation_names:

            # Get the parameter values
            parameter_values = self.parameters_table.parameter_values_for_simulation(simulation_name)

            # Prepare simulation directories, ski file, and return the simulation definition
            if not self.testing:
                definition = prepare_simulation(simulation_name, self.ski, parameter_values, self.object_name,
                                                self.simulation_input, self.generation_info.path, scientific=True, fancy=True,
                                                ndigits=self.fitting_run.ndigits_dict)
            else: definition = make_test_definition(simulation_name, self.ski, parameter_values, self.object_name,
                                                    self.simulation_input, scientific=True, fancy=True, ndigits=self.fitting_run.ndigits_dict)

            # Put the parameters in the queue and get the simulation object
            self.launcher.add_to_queue(definition, simulation_name)

        # Set the TEST flag if testing
        self.launcher.config.test = self.testing

        # Run the launcher, launches the simulations and retrieves and analyses finished simulations
        try: self.launcher.run(ncells=self.ndust_cells)
        except Exception as e:

            # Raise the exception again if we are just testing
            if self.testing: raise e

            # Something went wrong launching the simulations, show error message
            log.error("No simulations could be launched: removing generation ...")
            log.error(str(e))
            if log.is_debug: traceback.print_exc()
            log.error("Try again later")
            log.error("Cleaning up generation and quitting ...")

            # Remove this generation from the generations table
            self.fitting_run.generations_table.remove_entry(self.generation_name)
            self.fitting_run.generations_table.save()

            # Remove the generation directory
            fs.remove_directory(self.generation_path)

        # Check the launched simulations
        if not self.testing: self.check_simulations()

        # Save the succesful simulation files in their own directories
        self.save_simulations()

    # -----------------------------------------------------------------

    def set_finishing_time(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting the generation finishing time (there were no simulations for this generation) ...")

        # Set the time and save the table
        self.fitting_run.generations_table.set_finishing_time(self.generation_name, time.timestamp())
        self.fitting_run.generations_table.save()

    # -----------------------------------------------------------------

    @property
    def simulations(self):

        """
        This function ...
        :return:
        """

        return self.launcher.launched_simulations

    # -----------------------------------------------------------------

    def check_simulations(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the simulations ...")

        # Check the number of simulations that were effectively launched
        if self.nmodels == len(self.simulations):
            log.success("All simulations were scheduled succesfully")
            return

        # No simulations were launched
        if len(self.simulations) == 0:

            # Show error message
            log.error("No simulations could be launched: removing generation ...")
            log.error("Try again later")
            log.error("Cleaning up generation and quitting ...")

            # Remove this generation from the generations table
            self.fitting_run.generations_table.remove_entry(self.generation_name)
            self.fitting_run.generations_table.save()

            # Remove the generation directory
            fs.remove_directory(self.generation_path)

            # Quit
            exit()

        # Less simulations were launched
        elif len(self.simulations) < self.nmodels:

            # Get the names of simulations that were launched
            launched_simulation_names = [simulation.name for simulation in self.simulations]
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

    def save_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving the simulation files in the generation's simulation directories ...")

        # Loop over the simulations
        for simulation in self.simulations:

            # Check whether the simulation is succesfully launched
            if simulation.name not in self.parameters_table.simulation_names: continue

            # Determine the filepath
            filepath = fs.join(simulation.base_path, "initial.sim")

            # Save the simulation file
            simulation.saveto(filepath, update_path=False)

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

    @property
    def needs_simulations(self):

        """
        This function ...
        :return: 
        """

        return len(self.simulation_names) > 0

    # -----------------------------------------------------------------

    @property
    def generation_path(self):

        """
        This function ...
        :return: 
        """

        return self.generation_info.path

    # -----------------------------------------------------------------

    @property
    def run_name(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.name

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
            simulation_name = self.object_name + "__" + self.run_name + "__" + self.generation_name + "__" + str(counter)

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
            for label in parameter_values: log.debug(" - " + label + ": " + tostr(parameter_values[label], scientific=True, fancy=True, ndigits=self.fitting_run.ndigits_dict[label]))
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
        self.write_generation_info()

        # 2. Write the generations table
        self.write_generations_table()

        # 2. Write the individuals table
        self.write_individuals()

        # 3. Write the parameters table
        self.write_parameters()

        # 4. Write the (empty) chi squared table
        self.write_chi_squared()

    # -----------------------------------------------------------------

    def write_generation_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the generation info ...")

        # Save as a data file
        self.generation_info.saveto(self.generation.info_path)

    # -----------------------------------------------------------------

    def write_generations_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the generations table ...")

        # Add an entry to the generations table
        self.fitting_run.generations_table.add_entry(self.generation_info, self.ranges, self.parameter_scales)

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

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Show the ranges
        self.show_ranges()

        # Show the number of points
        self.show_npoints()

        # Show the simulation options
        self.show_simulation_options()

        # Show instruments
        self.show_instruments()

        # Show execution options
        self.show_execution_options()

        # Show analysis options
        self.show_analysis_options()

    # -----------------------------------------------------------------

    def show_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the parameter ranges ...")

        print("")
        print(fmt.green + fmt.underlined + "Parameter ranges:" + fmt.reset)
        print("")

        # Loop over the parameters
        for label in self.ranges:

            # Get the range
            parameter_range = self.ranges[label]

            # Show
            print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(parameter_range))

        print("")

    # -----------------------------------------------------------------

    def show_npoints(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the number of points ...")

        print("")
        print(fmt.green + fmt.underlined + "Number of grid points:" + fmt.reset)
        print("")

        # Loop over the parameters
        for label in self.npoints:

            # Get the npoints
            npoints = self.npoints[label]

            # Show
            print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(npoints))

        print("")

    # -----------------------------------------------------------------

    def show_simulation_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing simulation options ...")

        print("")
        print(fmt.green + fmt.underlined + "Simulation options:" + fmt.reset)
        print("")

        print(" - number of wavelengths: " + tostr(self.nwavelengths) + " (" + self.wavelength_grid_name + ")")
        print(" - number of dust cells: " + tostr(self.ndust_cells, scientific_int=False) + " (" + self.representation_name + ")")
        print(" - number of photon packages per wavelength: " + tostr(self.npackages, scientific_int=False))
        print(" - selfabsorption: " + tostr(self.selfabsorption))
        print(" - transient heating: " + tostr(self.transient_heating))
        print(" - dust grid type: " + tostr(self.representation.dust_grid_type))
        if self.representation.has_dust_grid_tree_distribution:
            print(" - dust grid minimum level: " + tostr(self.representation.dust_grid_min_level))
            print(" - dust_grid_maximum_level: " + tostr(self.representation.dust_grid_max_level))

        print("")

    # -----------------------------------------------------------------

    def show_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the instruments ...")

        print("")
        print(fmt.green + fmt.underlined + "Instruments:" + fmt.reset)
        print("")

        # Loop over the instruments
        for name in self.instruments:

            instrument = self.instruments[name]
            instr_class = str(type(instrument).__name__)

            print(" - " + fmt.bold + name + fmt.reset + " (" + instr_class + "):")
            print("")
            print(instrument.to_string(line_prefix="  ", bullet="*", bold=False))

        print("")

    # -----------------------------------------------------------------

    def show_execution_options(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Showing execution options ...")

        print("")
        print(fmt.green + fmt.underlined + "Execution:" + fmt.reset)
        print("")

        print(" - remote hosts: " + tostr(self.remote_host_ids))
        if self.config.cores_per_process is not None: print(" - number of cores per process: " + tostr(self.config.cores_per_process))
        else: print(" - number of cores per process determined automatically")

        if self.uses_remotes:
            if self.config.nprocesses_remote is not None and self.has_single_remote: nprocesses = self.config.nprocesses_remote
            else: nprocesses = None
            if self.config.data_parallel_remote is not None and self.has_single_remote: data_parallel = self.config.data_parallel_remote
            else: data_parallel = None
        else:
            nprocesses = self.config.nprocesses_local
            data_parallel = self.config.data_parallel_local

        if nprocesses is not None: print(" - number of processes: " + tostr(nprocesses))
        else: print(" - number of processes determined automatically")
        if data_parallel is not None: print(" - data parallelization: " + tostr(data_parallel))
        else: print(" - data parallelization enabled or disabled automatically")

        print("")

    # -----------------------------------------------------------------

    def show_analysis_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing analysis options ...")

        print("")
        print(fmt.green + fmt.underlined + "Analysis options:" + fmt.reset)
        print("")

        print(" - spectral convolution: " + tostr(self.spectral_convolution))
        print(" - use images: " + tostr(self.use_images))

        print(" - extract progress: " + tostr(self.config.extract_progress))
        print(" - extract timeline: " + tostr(self.extract_timeline))
        print(" - extract memory: " + tostr(self.extract_memory))
        print(" - plotting format: " + tostr(self.config.plotting_format))

        # From images
        if self.use_images:

            print(" - instrument reference image: " + tostr(self.reference_component_name))
            print(" - reference image xsize: " + tostr(self.reference_wcs.xsize))
            print(" - reference image ysize: " + tostr(self.reference_wcs.ysize))
            print(" - reference image pixelscale: " + tostr(self.reference_wcs.average_pixelscale.to("arcsec")))

            # For plotting:
            #print(" - fluxes_from_images_masks: " + tostr(self.environment.photometry_image_paths_for_filter_names))

        # Observation filters
        print(" - observation_filters: " + tostr(self.observed_filter_names))

        print("")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot reference images
        self.plot_reference_images()

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_images_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.fitting_run.path, "refimages__" + self.generation_name)

    # -----------------------------------------------------------------

    @property
    def fitting_filters(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.fitting_filters

    # -----------------------------------------------------------------

    @property
    def fitting_filter_names(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.fitting_filters]

    # -----------------------------------------------------------------

    def plot_reference_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the reference images ...")

        # Create plotter
        plotter = StandardImageGridPlotter()

        # Set output directory
        plotter.config.output = self.reference_images_path

        # Extra
        plotter.config.normalize = True
        plotter.config.colormap = self.config.reference_images_colormap

        # Write data
        plotter.config.write = self.config.write_reference_images

        # Rebin and crop
        plotter.rebin_to = self.reference_wcs
        plotter.crop_to = self.environment.truncation_box

        # Loop over the filters
        for fltr in self.environment.photometry_image_paths_for_filters:

            # Check whether fitting filter
            #if fltr not in self.fitting_filters: continue

            # Get path
            path = self.environment.photometry_image_paths_for_filters[fltr]

            # Add to plot
            plotter.add_image_from_file(path, masks=False, regions=False)

        # Run the plotter
        plotter.run()

# -----------------------------------------------------------------
