#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeling.base Contains the ModelerBase class, which is the base class for the specific modelers
#  such as the GalaxyModeler, SEDModeler and ImagesModeler.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable, write_input
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ..fitting.explorer import ParameterExplorer
from ..fitting.sedfitting import SEDFitter
from ..component.component import load_modeling_history, get_config_file_path, load_modeling_configuration
from ...core.launch.synchronizer import RemoteSynchronizer
from ...core.prep.deploy import Deployer
from ..fitting.run import get_generations_table, has_unevaluated_generations, get_unevaluated_generations
from ...core.remote.moderator import PlatformModerator
from ...core.tools import stringify
from ...core.tools.loops import repeat_check
from ...core.remote.remote import Remote
from ..fitting.finisher import ExplorationFinisher
from ...core.tools import time
from ..core.steps import commands_after_and_including, output_paths_for_single_command, cached_directory_name_for_single_command
from ...core.basics.configuration import prompt_proceed
from ...core.tools.stringify import tostr
from ...core.tools.utils import lazyproperty
from ..setup import ignore_commands
from ...core.tools.utils import DefaultScope

# -----------------------------------------------------------------

fitting_methods = ["genetic", "grid"]
default_fitting_method = "genetic"

# -----------------------------------------------------------------

class ModelerBase(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ModelerBase, self).__init__(*args, **kwargs)

        # A timestamp
        self.timestamp = None

        # The path to the modeling directory
        self.modeling_path = None

        # The modeling environment
        self.environment = None

        # The modeling configuration
        self.modeling_config = None

        # Platform moderator
        self.moderator = None

        # The modeling history
        self.history = None

        # Fixed names for the fitting run and the model
        self.fitting_run_name = "run_1"
        self.model_name = "model_a"
        self.representation_name = "highres"

        # Parameter ranges
        self.parameter_ranges = None

        # The parameter explorer instance
        self.explorer = None

        # The SED fitter instance
        self.fitter = None

        # The exploration finisher
        self.finisher = None

        # egege
        self.fixed_initial_parameters = None

        # The fitting method (grid, genetic)
        self.fitting_method = None

    # -----------------------------------------------------------------

    @property
    def is_galaxy_modeler(self):

        """
        This function ...
        :return:
        """

        if self.environment is None: raise ValueError("Not yet known")

        from ..core.environment import GalaxyModelingEnvironment
        return isinstance(self.environment, GalaxyModelingEnvironment)

    # -----------------------------------------------------------------

    @property
    def is_sed_modeler(self):

        """
        This function ...
        :return:
        """

        if self.environment is None: raise ValueError("Not yet known")

        from ..core.environment import SEDModelingEnvironment
        return isinstance(self.environment, SEDModelingEnvironment)

    # -----------------------------------------------------------------

    @property
    def is_images_modeler(self):

        """
        This function ...
        :return:
        """

        if self.environment is None: raise ValueError("Not yet known")

        from ..core.environment import ImagesModelingEnvironment
        return isinstance(self.environment, ImagesModelingEnvironment)

    # -----------------------------------------------------------------

    def set_rerun(self):

        """
        This function ...
        :return:
        """

        # Get the commands which have to be removed from the history
        commands = commands_after_and_including(self.config.rerun)

        # Loop over the commands, remove all entries
        for command_name in commands:

            # SKip if not yet present
            if command_name not in self.history: continue

            # Debugging
            log.debug("Removing the '" + command_name + "' command from the modeling history and removing the output ...")

            # Remove
            #self.history.remove_entries_and_save(command_name)

            # Remove the output
            self.remove_output_and_history_for_command(command_name)

    # -----------------------------------------------------------------

    def output_paths_for_command(self, command_name):

        """
        This function ...
        :param command_name:
        :return:
        """

        return output_paths_for_single_command(self.environment, command_name)

    # -----------------------------------------------------------------

    def cached_directory_name_for_command(self, command_name):

        """
        This function ...
        :param command_name:
        :return:
        """

        return cached_directory_name_for_single_command(self.environment, command_name)

    # -----------------------------------------------------------------

    @property
    def cache_host_id(self):

        """
        THis function ..
        :return:
        """

        return self.modeling_config.cache_host_id

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_remote(self):

        """
        This function ...
        :return:
        """

        return Remote(host_id=self.cache_host_id)

    # -----------------------------------------------------------------

    def remove_output_for_command(self, command_name):

        """
        This function ...
        :param command_name:
        :return:
        """

        # Determine output path
        paths = self.output_paths_for_command(command_name)

        # Debugging
        log.debug("Removing the '" + stringify.stringify_paths(paths, base=self.modeling_path)[1] + "' directory ...")

        # Remove
        fs.remove_directories_and_files(paths)

        # Remove cached data remotely
        cached_directory_name = self.cached_directory_name_for_command(command_name)
        if cached_directory_name is not None:
            remote_path = fs.join(self.cache_remote.home_directory, cached_directory_name)
            if self.cache_remote.is_directory(remote_path): self.cache_remote.remove_directory(remote_path)

    # -----------------------------------------------------------------

    def remove_output_and_history_for_command(self, command_name):

        """
        This function ...
        :param command_name:
        :return:
        """

        # Ask to make sure
        if prompt_proceed("are you absolutely sure the output of '" + command_name + "' can be removed?"):

            # Remove the output
            self.remove_output_for_command(command_name)

            # Remove from history file (as if it had never been run before)
            self.history.remove_entries(command_name)

            # Return True, step needs to be performed again
            return True

        # User doesn't want to proceed
        else:

            # Exit with an error message
            log.error("Cannot rerun the " + self.config.rerun + " command when the " + command_name + " cannot be rerun")
            exit()

    # -----------------------------------------------------------------

    def check_needs_step(self, command_name):

        """
        This function ...
        :param command_name:
        :return:
        """

        # Step already finished
        if self.history.is_finished(command_name):

            # Rerun?
            if command_name in self.rerun_commands: return self.remove_output_and_history_for_command(command_name)

            # No rerun
            else: return False

        # Step already started, partial output?
        elif command_name in self.history:

            # Verify that the output can be removed
            if command_name in self.rerun_commands: return self.remove_output_and_history_for_command(command_name)

            # Rerun isn't enabled for this command, but start anyway because it hasn't finished before
            else:

                # Give a warning before starting
                log.warning("The " + command_name + " has been started before but hasn't finished. The commmand will be launched again but partial output may already be present.")
                return True

        # Nothing in history
        else: return True

    # -----------------------------------------------------------------

    def log_path_for_component(self, cls_or_instance):

        """
        This function ...
        :param cls_or_instance:
        :return:
        """

        command_name = cls_or_instance.command_name()
        #print(self.environment.log_path, command_name, self.timestamp)
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        return log_path

    # -----------------------------------------------------------------

    def add_log_path(self, cls_or_instance):

        """
        This function ...
        :param cls_or_instance:
        :return:
        """

        log.add_log_file(self.log_path_for_component(cls_or_instance))

    # -----------------------------------------------------------------

    def write_log(self, cls_or_instance):

        """
        This function ...
        :param cls_or_instance:
        :return:
        """

        return log.log_to_file(self.log_path_for_component(cls_or_instance), filter_level="DEBUG")

    # -----------------------------------------------------------------

    def config_path_for_component(self, cls_or_instance):

        """
        This function ...
        :param cls_or_instance:
        :return:
        """

        command_name = cls_or_instance.command_name()
        config_path = fs.join(self.environment.config_path, command_name + "__" + self.timestamp + ".cfg")
        return config_path

    # -----------------------------------------------------------------

    def register(self, cls_or_instance):

        """
        This function ...
        :param cls_or_instance:
        :return:
        """

        # Get the command name
        command_name = cls_or_instance.command_name()

        # Ignore not-pipeline
        if command_name in ignore_commands: return DefaultScope()

        # Register scope
        return self.history.register(cls_or_instance)

    # -----------------------------------------------------------------

    def write_config(self, instance):

        """
        This function ...
        :param instance:
        :return:
        """

        #return keep_latest_succesful_config(instance, self.config_path_for_component(instance))
        return tag_latest_succesful_config(instance, self.config_path_for_component(instance))

    # -----------------------------------------------------------------

    def input_path_for_component(self, cls_or_instance):

        """
        This function ...
        :param cls_or_instance:
        :return:
        """

        command_name = cls_or_instance.command_name()
        #in_path = fs.create_directory_in(self.environment.in_path, command_name + "_" + self.timestamp)
        in_path = fs.join(self.environment.in_path, command_name)
        return in_path

    # -----------------------------------------------------------------

    def write_input(self, instance, **kwargs):

        """
        This function ...
        :param instance:
        :param kwargs:
        :return:
        """

        return keep_last_succesful_input(instance, self.input_path_for_component(instance), **kwargs)

    # -----------------------------------------------------------------

    @property
    def grid_fitting(self):

        """
        This function ...
        :return: 
        """

        return self.fitting_method == "grid"

    # -----------------------------------------------------------------

    @property
    def genetic_fitting(self):

        """
        This function ...
        :return: 
        """

        return self.fitting_method == "genetic"

    # -----------------------------------------------------------------

    @property
    def configured_fitting_host_ids(self):

        """
        This function ...
        :return:
        """

        if self.modeling_config.fitting_host_ids is None: return []
        else: return self.modeling_config.fitting_host_ids

    # -----------------------------------------------------------------

    @property
    def has_configured_fitting_host_ids(self):

        """
        This function ...
        :return:
        """

        return len(self.configured_fitting_host_ids) > 0

    # -----------------------------------------------------------------

    @property
    def multiple_generations(self):

        """
        This function ...
        :return:
        """

        return self.config.ngenerations > 1

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ModelerBase, self).setup(**kwargs)

        # Create timestamp
        self.timestamp = time.filename_timestamp()

        # Set the path to the modeling directory
        self.modeling_path = self.config.path

        # Check for the presence of the configuration file
        if not fs.is_file(get_config_file_path(self.modeling_path)): raise ValueError("The current working directory (" + self.config.path + ") is not a radiative transfer modeling directory (the configuration file is missing)")
        else: self.modeling_config = load_modeling_configuration(self.modeling_path)

        # Set execution platforms
        self.set_platforms()

        # Check the number of generations
        if self.config.ngenerations > 1 and self.moderator.any_remote: raise ValueError("When remote execution is enabled, the number of generations per run can only be one")

        # Load the modeling history
        self.history = load_modeling_history(self.modeling_path)

        # Clear remotes
        if self.config.clear_remotes: self.moderator.clear_all_hosts()

        # Deploy SKIRT and PTS
        if self.config.deploy: self.deploy()

        # Set the fitting method
        if "fitting_method" in kwargs: self.fitting_method = kwargs.pop("fitting_method")
        else: self.fitting_method = self.config.fitting_method #self.fitting_method = prompt_string("fitting_method", "fitting method", choices=fitting_methods)

    # -----------------------------------------------------------------

    @property
    def fitting_local(self):

        """
        This function ...
        :return: 
        """

        # If local flag is set
        if self.config.fitting_local: return True

        # Fitting remotes have been set
        elif self.config.fitting_remotes is not None: return True

        # If no host ids have been set
        else: return self.modeling_config.fitting_host_ids is None or self.modeling_config.fitting_host_ids == []

    # -----------------------------------------------------------------

    @property
    def other_local(self):

        """
        This function ...
        :return: 
        """

        # If local flag is set
        if self.config.local: return True

        # Remotes have been set
        elif self.config.remotes is not None: return True

        # If no host ids have been set
        else: return self.modeling_config.host_ids is None or self.modeling_config.host_ids == []

    # -----------------------------------------------------------------

    @property
    def fitting_host_ids(self):

        """
        This function ...
        :return: 
        """

        if self.fitting_local: return None
        elif self.config.fitting_remotes is not None: return self.config.fitting_remotes
        elif self.modeling_config.fitting_host_ids is None: return None
        else: return self.modeling_config.fitting_host_ids

    # -----------------------------------------------------------------

    @property
    def other_host_ids(self):

        """
        This function ...
        :return: 
        """

        if self.config.local: return None
        elif self.config.remotes is not None: return self.config.remotes
        elif self.modeling_config.host_ids is None: return None
        else: return self.modeling_config.host_ids

    # -----------------------------------------------------------------

    def set_platforms(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining execution platforms ...")

        # Setup the platform moderator
        self.moderator = PlatformModerator()

        # Set platform(s) for fitting (simulations)
        if self.fitting_local: self.moderator.add_local("fitting")
        else: self.moderator.add_ensemble("fitting", self.fitting_host_ids)

        # Other computations
        if self.other_local: self.moderator.add_local("other")
        else: self.moderator.add_single("other", self.other_host_ids)

        # Run the moderator
        self.moderator.run()

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
        deployer.config.hosts = self.moderator.all_hosts

        # Set the host id on which PTS should be installed (on the host for extra computations and the fitting hosts
        # that have a scheduling system to launch the pts run_queue command)
        deployer.config.pts_on = self.moderator.all_host_ids

        # Set
        deployer.config.check = self.config.check_versions

        # Set
        deployer.config.update_dependencies = self.config.update_dependencies

        # Run the deployer
        deployer.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def rerun_commands(self):

        """
        This function ...
        :return:
        """

        if self.config.rerun is None: return []
        else: return commands_after_and_including(self.config.rerun)

    # -----------------------------------------------------------------

    def fit(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Fitting radiative transfer models ...")

        # Configure the fitting
        if not self.history.has_configured_fit: self.configure_fit()

        # Initialize the fitting
        if not self.history.has_initialized_fit: self.initialize_fit()

        # If we do multiple generations at once
        if self.multiple_generations: self.fit_multiple(**kwargs)

        # We just do one generation now, or finish
        else: self.fit_single(**kwargs)

    # -----------------------------------------------------------------

    def fit_multiple(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return: 
        """

        # Inform the user
        log.info("Fitting multiple generations at once ...")

        # Start: launch the initial generation
        self.start(**kwargs)

        # Advance: launch generations 0 -> (n-1)
        repeat_check(self.advance, self.config.ngenerations, **kwargs)

        # Finish
        self.finish(**kwargs)

    # -----------------------------------------------------------------

    def fit_single(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Inform the user
        log.info("Fitting a single generation ...")

        # Load the generations table
        generations = get_generations_table(self.modeling_path, self.fitting_run_name)

        # If finishing the generation is requested
        if self.config.finish: self.finish(**kwargs)

        # If this is the initial generation
        elif generations.last_generation_name is None: self.start(**kwargs)

        # Advance the fitting with a new generation
        else: self.advance(**kwargs)

    # -----------------------------------------------------------------

    @abstractmethod
    def configure_fit(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def initialize_fit(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def start(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Starting with the initial generation ...")

        # Explore
        self.explore(**kwargs)

    # -----------------------------------------------------------------

    def advance(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Advancing the fitting with a new generation ...")

        # Load the generations table
        generations = get_generations_table(self.modeling_path, self.fitting_run_name)

        # Check whether there is a generation preceeding this one
        if generations.last_generation_name is None: raise RuntimeError("Preceeding generation cannot be found")

        # Debugging
        log.debug("Previous generation: " + generations.last_generation_name)

        # If some generations have not finished, check the status of and retrieve simulations
        if generations.has_unfinished and self.has_configured_fitting_host_ids: self.synchronize()

        # Debugging
        if generations.has_finished: log.debug("There are finished generations: " + stringify.stringify(generations.finished_generations)[1])
        if has_unevaluated_generations(self.modeling_path, self.fitting_run_name): log.debug("There are unevaluated generations: " + stringify.stringify(get_unevaluated_generations(self.modeling_path, self.fitting_run_name))[1])

        # If some generations have finished, fit the SED
        if generations.has_finished and has_unevaluated_generations(self.modeling_path, self.fitting_run_name): self.fit_sed()

        # If all generations have finished, explore new generation of models
        if generations.all_finished:

            # Explore a new generation
            self.explore(**kwargs)
            return True

        # Return False if exploration could not be performed (not all generations had finished)
        else: return False

    # -----------------------------------------------------------------

    def synchronize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Synchronizing with the remotes (retrieving and analysing finished models) ...")

        # Create the remote synchronizer
        synchronizer = RemoteSynchronizer()

        # Set the host IDs
        synchronizer.config.host_ids = self.modeling_config.fitting_host_ids

        # Run the remote synchronizer
        synchronizer.run()

    # -----------------------------------------------------------------

    def fit_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting the SED to the finished generations ...")

        # Configuration settings
        config = dict()
        config["name"] = self.fitting_run_name

        # Create the SED fitter
        self.fitter = SEDFitter(config)

        # Run the fitter
        with self.write_log(self.fitter), self.history.register(self.fitter), self.write_config(self.fitter): self.fitter.run()

    # -----------------------------------------------------------------

    def explore(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Exploring the parameter space ...")

        # Configuration settings
        config = dict()
        config["name"] = self.fitting_run_name

        # Set flags
        if self.moderator.ensemble_is_local("fitting"):
            config["record_timing"] = False
            config["record_memory"] = False
        else:
            config["record_timing"] = True
            config["record_memory"] = True

        # Set recurrence settings
        config["check_recurrence"] = self.config.check_recurrence
        config["recurrence_rtol"] = self.config.recurrence_rtol
        config["recurrence_atol"] = self.config.recurrence_atol

        # Create the parameter explorer
        self.explorer = ParameterExplorer(config, cwd=self.modeling_path)

        # Set the working directory
        self.explorer.config.path = self.modeling_path

        # Set the remote host IDs
        self.explorer.config.remotes = self.moderator.host_ids_for_ensemble("fitting")
        self.explorer.config.attached = self.config.fitting_attached

        # Set the number of generations
        #if self.config.ngenerations is not None: explorer.config.ngenerations = self.config.ngenerations
        # NO: THIS ALWAYS HAVE TO BE ONE: BECAUSE HERE IN THIS CLASS WE ALREADY USE REPEAT(SELF.ADVANCE)
        # IF NGENERATIONS > 1, THE CONTINUOUSOPTIMIZER IS USED INSTEAD OF THE STEPWISEOPTIMIZER
        self.explorer.config.ngenerations = 1

        # Set the number of simulations per generation
        if self.config.nsimulations is not None: self.explorer.config.nsimulations = self.config.nsimulations

        # Set other settings
        self.explorer.config.npackages_factor = self.config.npackages_factor
        self.explorer.config.increase_npackages = self.config.increase_npackages
        #explorer.config.refine_wavelengths = self.config.refine_wavelengths
        self.explorer.config.refine_spectral = self.config.refine_spectral
        #explorer.config.refine_dust = self.config.refine_dust
        self.explorer.config.refine_spatial = self.config.refine_spatial
        self.explorer.config.selfabsorption = self.config.selfabsorption
        self.explorer.config.transient_heating = self.config.transient_heating

        # ADVANCED SETTINGS
        self.explorer.config.restart_from_generation = self.config.restart_from_generation

        # Set the input
        input_dict = dict()
        if self.parameter_ranges is not None: input_dict["ranges"] = self.parameter_ranges

        # Add the fixed parameter values
        if self.fixed_initial_parameters is not None: input_dict["fixed_initial_parameters"] = self.fixed_initial_parameters

        # NEW: Add additional input (such as parameter grid scales)
        input_dict.update(kwargs)

        # Run the parameter explorer
        with self.write_log(self.explorer), self.history.register(self.explorer), self.write_config(self.explorer), self.write_input(self.explorer, **input_dict): self.explorer.run(**input_dict)

    # -----------------------------------------------------------------

    def finish(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Finishing the parameter exploration ...")

        # Configuration settings
        settings = dict()
        settings["name"] = self.fitting_run_name

        # Set the input
        input_dict = dict()

        # Create the exploration finisher
        self.finisher = ExplorationFinisher(settings)

        # NEW: Add additional input (such as parameter grid scales)
        input_dict.update(kwargs)

        # Run the finisher
        with self.write_log(self.finisher), self.history.register(self.finisher), self.write_config(self.finisher), self.write_input(self.finisher, **input_dict): self.finisher.run(**input_dict)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------

class tag_latest_succesful_config(object):

    """
    This class ...
    """

    def __init__(self, instance, config_path, tag="*"):

        """
        Thisn function ...
        :param instance:
        :param config_path:
        :param tag:
        """

        # Set
        self.instance = instance
        self.config_path = config_path
        self.tag = tag

    # -----------------------------------------------------------------

    @property
    def command_name(self):

        """
        This function ...
        :return:
        """

        return self.instance.command_name()

    # -----------------------------------------------------------------

    @property
    def config_directory_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.config_path)

    # -----------------------------------------------------------------

    @property
    def config_filename(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.config_path)

    # -----------------------------------------------------------------

    @property
    def config_name(self):

        """
        This function ...
        :return:
        """

        return fs.strip_extension(self.config_filename)

    # -----------------------------------------------------------------

    def __enter__(self):

        """
        This function ...
        :return:
        """

        # Set the config path to the instance configuration
        self.instance.config.write_config = True
        self.instance.config.config_path = self.config_path

    # -----------------------------------------------------------------

    def __exit__(self, exc_type, exc_val, exc_tb):

        """
        This function ...
        :param exc_type:
        :param exc_val:
        :param exc_tb:
        :return:
        """

        # Error occured
        if exc_type is not None:

            # Warning
            log.error("A " + str(exc_type.__name__) + " occured")

        # No error: tag
        else:

            # Debugging
            log.debug("Tagging configuration file as succesful ...")

            # Rename by adding the tag
            fs.add_suffix(self.config_path, self.tag)

# -----------------------------------------------------------------

class keep_latest_succesful_config(object):

    """
    This class ...
    """

    def __init__(self, instance, config_path):

        """
        This function ...
        :param instance:
        :param config_path:
        """

        # Set
        self.instance = instance
        self.config_path = config_path

    # -----------------------------------------------------------------

    @property
    def command_name(self):

        """
        This function ...
        :return:
        """

        return self.instance.command_name()

    # -----------------------------------------------------------------

    @property
    def config_directory_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.config_path)

    # -----------------------------------------------------------------

    @property
    def config_filename(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.config_path)

    # -----------------------------------------------------------------

    @property
    def config_name(self):

        """
        This function ...
        :return:
        """

        return fs.strip_extension(self.config_filename)

    # -----------------------------------------------------------------

    def __enter__(self):
        
        """
        This function ...
        :return: 
        """

        # Set the config path to the instance configuration
        self.instance.config.write_config = True
        self.instance.config.config_path = self.config_path

    # -----------------------------------------------------------------
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        
        """
        This function ...
        :param exc_type: 
        :param exc_val: 
        :param exc_tb: 
        :return: 
        """

        # Error occured
        if exc_type is not None:

            log.error("A " + str(exc_type.__name__) + " occured")
            log.error("Removing configuration file [" + self.config_path + "] ...")
            fs.remove_file(self.config_path)

        # No error: remove all others
        else:

            paths, names = fs.files_in_path(self.config_directory_path, contains=self.command_name, not_contains=self.config_name, extension="cfg", returns=["path", "name"], unpack=True)
            log.debug("Removing previous configuration files: " + tostr(names) + " ...")
            fs.remove_files(paths)

# -----------------------------------------------------------------

class keep_last_succesful_input(object):

    """
    This function ...
    """

    def __init__(self, instance, input_path, **kwargs):

        """
        This function ...
        :param instance:
        :param input_path:
        """

        # Set
        self.instance = instance
        self.input_path = input_path
        self.input = kwargs

    # -----------------------------------------------------------------

    @property
    def command_name(self):

        """
        This function ...
        :return:
        """

        return self.instance.command_name()

    # -----------------------------------------------------------------

    @property
    def input_directory_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.input_path)

    # -----------------------------------------------------------------

    @property
    def input_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.input_path)

    # -----------------------------------------------------------------

    def __enter__(self):

        """
        This function ...
        :return:
        """

        # Write the input directory
        fs.create_directory(self.input_path)

        # Save the input
        write_input(self.input, self.input_path, light=True)

    # -----------------------------------------------------------------

    def __exit__(self, exc_type, exc_val, exc_tb):

        """
        This function ...
        :param exc_type:
        :param exc_val:
        :param exc_tb:
        :return:
        """

        # Error occured
        if exc_type is not None:

            log.error("A " + str(exc_type.__name__) + " occured")
            log.error("Removing input path [" + self.input_path + "] ...")
            fs.remove_directory(self.input_path)

        # No error: remove all previous directories
        else:

            paths, names = fs.directories_in_path(self.input_directory_path, contains=self.command_name, not_contains=self.input_name, returns=["path", "name"], unpack=True)
            log.debug("Removing previous input directories: " + tostr(names) + " ...")
            fs.remove_directories(paths)

# -----------------------------------------------------------------
