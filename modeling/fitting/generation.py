#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.generation Contains the Generation class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from copy import deepcopy
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.composite import SimplePropertyComposite
from ...core.tools import filesystem as fs
from ..build.component import get_model_definition, get_representation
from .tables import IndividualsTable, ParametersTable, ChiSquaredTable, ParameterProbabilitiesTable, ModelProbabilitiesTable
from ...evolve.optimize.tables import ElitismTable, CrossoverTable, ScoresTable, RecurrenceTable
from ...evolve.optimize.stepwise import load_population
from ...core.basics.configurable import load_input
from ...core.basics.configuration import Configuration
from ...evolve.optimize.components import get_crossover, get_crossover_origins, get_genome_class, get_mutator, get_selector, get_scaling, get_initializator
from ...core.basics.range import IntegerRange, RealRange, QuantityRange
from ...core.tools import sequences
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.launch.tables import SimulationAssignmentTable
from ...core.simulation.remote import get_simulation_for_host, has_simulation_for_host, get_simulation_path_for_host, is_invalid_or_unknown_status
from ...core.remote.host import find_host_ids
from ...core.simulation.screen import ScreenScript
from ...core.simulation.simulation import SkirtSimulation, RemoteSimulation
from ...core.simulation.input import SimulationInput
from ...core.data.sed import ObservedSED, SED
from ...core.simulation.logfile import LogFile
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...magic.core.datacube import DataCube
from ...magic.core.frame import Frame
from ...core.basics.log import log
from ...core.simulation.remote import get_simulations_for_host
from ...core.simulation.remote import finished_name
from ...core.launch.tables import SimulationStatusTable
from ...core.simulation.adapter import adapt_simulation, adapt_analysis
from ...core.basics.configuration import prompt_yn, prompt_string
from ...core.tools.stringify import tostr
from ...core.launch.batch import MissingSimulation
from ...core.simulation.jobscript import SKIRTJobScript, get_host_id_from_jobscript_file
from ...core.simulation.skifile import SkiFile
from ...core.tools.progress import Bar

# -----------------------------------------------------------------

class GenerationInfo(SimplePropertyComposite):

    """
    This function ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(GenerationInfo, self).__init__()

        # Define properties
        self.add_string_property("name", "name of the generation")
        self.add_integer_property("index", "index of the generation")
        self.add_string_property("method", "generation method")
        self.add_string_property("wavelength_grid_name", "wavelength grid name")
        self.add_string_property("model_representation_name", "model representation name")
        self.add_integer_property("nsimulations", "number of simulations")
        self.add_integer_property("npackages", "number of packages")
        self.add_boolean_property("selfabsorption", "dust self-absorption enabled")
        self.add_boolean_property("transient_heating", "transient heating enabled")
        self.add_boolean_property("spectral_convolution", "spectral convolution enabled")
        self.add_boolean_property("use_images", "use images")
        self.add_boolean_property("fit_not_clipped", "fit to not-clipped image fluxes (truncated images)")
        self.add_string_property("path", "generation path")

        # Set properties
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_generation_path(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Determine info path
        info_path = fs.join(path, "info.dat")

        # Check if present
        if not fs.is_file(info_path): raise IOError("The generation info file is not present at '" + info_path + "'")

        # Load the info
        return cls.from_file(info_path)

# -----------------------------------------------------------------

cached_filename = "cached.dat"

# -----------------------------------------------------------------

class Generation(object):
    
    """
    This class...
    """
    
    def __init__(self, info):

        """
        The constructor ...
        :return:
        """

        # Set the generation info
        self.info = info

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, directory_path):

        """
        This function ...
        :param directory_path:
        :return:
        """

        # Load the info
        info = GenerationInfo.from_generation_path(directory_path)

        # Create the generation object and return
        return cls(info)

    # -----------------------------------------------------------------

    @property
    def info_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "info.dat")

    # -----------------------------------------------------------------

    @property
    def name(self):

        """
        This function ...
        :return:
        """

        return self.info.name

    # -----------------------------------------------------------------

    @property
    def path(self):

        """
        This function ..
        :return:
        """

        return self.info.path

    # -----------------------------------------------------------------

    @property
    def generations_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.path)

    # -----------------------------------------------------------------

    @property
    def fitting_run_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.generations_path)

    # -----------------------------------------------------------------

    @property
    def fit_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.fitting_run_path)

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.fit_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def object_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.modeling_path)

    # -----------------------------------------------------------------

    @property
    def has_screen_scripts(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.path, extension="sh")

    # -----------------------------------------------------------------

    @property
    def screen_script_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.path, extension="sh")

    # -----------------------------------------------------------------

    @lazyproperty
    def screen_names(self):

        """
        This function ...
        :return:
        """

        host_ids = find_host_ids(schedulers=False)
        names = []

        # Loop over the screen scripts
        for path in self.screen_script_paths:

            # Get filename
            filename = fs.strip_extension(fs.name(path))

            # Split
            splitted = filename.split("_")
            if len(splitted) == 1: continue

            # Get host ID
            host_id = splitted[-1]
            if host_id not in host_ids: continue

            # Add the screen name
            name = "_".join(splitted[:-1])
            names.append(name)

        # Return the screen names
        return names

    # -----------------------------------------------------------------

    @property
    def nscreens(self):

        """
        This function ...
        :return:
        """

        return len(self.screen_names)

    # -----------------------------------------------------------------

    def get_screen_name(self, script_path):

        """
        This function ...
        :param script_path:
        :return:
        """

        # Get all remote host IDs
        host_ids = find_host_ids(schedulers=False)

        # Get filename
        filename = fs.strip_extension(fs.name(script_path))

        # Split
        splitted = filename.split("_")
        if len(splitted) == 1: return None

        # Get host ID
        host_id = splitted[-1]
        if host_id not in host_ids: return None

        # Add the screen name
        name = "_".join(splitted[:-1])
        return name

    # -----------------------------------------------------------------

    def get_screen_script_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Loop over the screen scripts
        for path in self.screen_script_paths:

            # Check name
            screen_name = self.get_screen_name(path)
            if screen_name == name: return path

        # Error
        raise ValueError("No such screen")

    # -----------------------------------------------------------------

    def get_screen_script(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get path
        path = self.get_screen_script_path(name)

        # Load the screen script
        return ScreenScript.from_file(path)

    # -----------------------------------------------------------------

    @property
    def assignment_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "assignment.dat")

    # -----------------------------------------------------------------

    @property
    def has_assignment_table(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.assignment_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def assignment_table(self):

        """
        This function ...
        :return:
        """

        return SimulationAssignmentTable.from_file(self.assignment_table_path)

    # -----------------------------------------------------------------

    def set_id_for_simulation(self, simulation_name, id):

        """
        This function ...
        :param simulation_name:
        :param id:
        :return:
        """

        self.assignment_table.set_id_for_simulation(simulation_name, id)

    # -----------------------------------------------------------------

    def set_host_for_simulation(self, simulation_name, host_id, cluster_name=None):

        """
        This function ...
        :param simulation_name:
        :param host_id:
        :param cluster_name:
        :return:
        """

        self.assignment_table.set_host_for_simulation(simulation_name, host_id, cluster_name=cluster_name)

    # -----------------------------------------------------------------

    def set_id_and_host_for_simulation(self, simulation_name, id, host_id, cluster_name=None):

        """
        This function ...
        :param simulation_name:
        :param id:
        :param host_id:
        :param cluster_name:
        :return:
        """

        self.assignment_table.set_id_and_host_for_simulation(simulation_name, id, host_id, cluster_name=cluster_name)

    # -----------------------------------------------------------------

    def add_simulation_to_assignment(self, simulation_name, host_id=None, cluster_name=None, simulation_id=None, success=True):

        """
        This function ...
        :param simulation_name:
        :param host_id:
        :param cluster_name:
        :param simulation_id:
        :param success:
        :return:
        """

        self.assignment_table.add_simulation(simulation_name, host_id=host_id, cluster_name=cluster_name, simulation_id=simulation_id, success=success)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_directory_names(self):

        """
        This function ...
        :return:
        """

        # Initialize list for the directory names
        dirnames = []

        # As from now, simulation names all start with the name of the galaxy
        dirnames += fs.directories_in_path(self.path, returns="name", startswith=self.object_name)

        # Before last change, simulation names started with the name of the fitting run
        dirnames += fs.directories_in_path(self.path, returns="name", startswith=self.fitting_run_name)

        # Initially, simulation names were just timestamps
        #if len(dirnames) == 0: dirnames = fs.directories_in_path(self.path, returns="name")
        dirnames += fs.directories_in_path(self.path, returns="name")

        # Return the simulation directory names
        return sequences.unique_values(dirnames)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        #directory_names = self.simulation_directory_names
        #simulation_names = self.individuals_table.simulation_names
        #if not sequences.same_contents(directory_names, simulation_names): raise ValueError("Mismatch between individuals table simulation names and directory names")
        #return simulation_names

        # Get the names from the individuals table
        names = self.individuals_table.simulation_names

        # Check whether there is a directory for each name
        if not sequences.contains_all(self.simulation_directory_names, names):
            raise ValueError("Some simulation directories are missing: " + tostr(sequences.elements_not_in_other(names, self.simulation_directory_names)))

        # Return the names
        return names

    # -----------------------------------------------------------------

    def get_simulation_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.simulation_names:
            if name not in self.simulation_directory_names: raise ValueError("Simulation does not exist in this generation")
            else:
                log.warning("Simulation '" + name + "' is missing from the individuals table ...")
        return fs.join(self.path, name)

    # -----------------------------------------------------------------

    def get_simulation_ski_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.get_simulation_path(name), self.object_name + ".ski")

    # -----------------------------------------------------------------

    def has_ski_file(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_ski_path(name)
        return fs.is_file(path)

    # -----------------------------------------------------------------

    def has_skifile(self, name):
        return self.has_ski_file(name)

    # -----------------------------------------------------------------

    def has_simulation_skifile(self, name):
        return self.has_skifile(name)

    # -----------------------------------------------------------------

    def get_simulation_skifile(self, name):
        return SkiFile(self.get_simulation_ski_path(name))

    # -----------------------------------------------------------------

    def get_simulation_output_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.get_simulation_path(name), "out")

    # -----------------------------------------------------------------

    def has_simulation_output(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_output_path(name)
        return fs.is_directory(path) and not fs.is_empty(path)

    # -----------------------------------------------------------------

    def get_simulation_logfile_path(self, name):

        """
        :param name:
        :return:
        """

        output_path = self.get_simulation_output_path(name)
        return fs.join(output_path, self.object_name + "_log.txt")

    # -----------------------------------------------------------------

    def has_simulation_logfile(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_simulation_logfile_path(name))

    # -----------------------------------------------------------------

    def has_logfile(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.has_simulation_logfile(name)

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation_logfile(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_logfile_path(name)
        return LogFile.from_file(path)

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation_jobscript_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return sequences.get_single(fs.files_in_path(self.get_simulation_path(name), extension="sh", startswith="job"), method="none")

    # -----------------------------------------------------------------

    def has_simulation_jobscript(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_simulation_jobscript_path(name) is not None

    # -----------------------------------------------------------------

    def has_jobscript(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.has_simulation_jobscript(name)

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation_jobscript(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_jobscript_path(name)
        return SKIRTJobScript.from_file(path)

    # -----------------------------------------------------------------

    def get_simulation_sed_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        output_path = self.get_simulation_output_path(name)
        return fs.join(output_path, self.object_name + "_earth_sed.dat")

    # -----------------------------------------------------------------

    def has_simulation_sed(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_simulation_sed_path(name))

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation_sed(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_sed_path(name)
        return SED.from_skirt(path)

    # -----------------------------------------------------------------

    def get_simulation_datacube_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        output_path = self.get_simulation_output_path(name)
        return fs.join(output_path, self.object_name + "_earth_total.fits")

    # -----------------------------------------------------------------

    def has_simulation_datacube(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_simulation_datacube_path(name))

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation_datacube(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_datacube_path(name)
        return DataCube.from_file(path, self.wavelength_grid)

    # -----------------------------------------------------------------

    def get_simulation_extract_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.get_simulation_path(name), "extr")

    # -----------------------------------------------------------------

    def has_extraction_output(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_extract_path(name)
        return fs.is_directory(path) and not fs.is_empty(path)

    # -----------------------------------------------------------------

    def get_simulation_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.get_simulation_path(name), "plot")

    # -----------------------------------------------------------------

    def has_plotting_output(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_plot_path(name)
        return fs.is_directory(path) and not fs.is_empty(path)

    # -----------------------------------------------------------------

    def get_simulation_misc_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.get_simulation_path(name), "misc")

    # -----------------------------------------------------------------

    def get_simulation_sed_plot_path(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        # Determine SED plot path
        plot_path = self.get_simulation_plot_path(name)
        return fs.join(plot_path, "sed.pdf")

    # -----------------------------------------------------------------

    def has_sed_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_sed_plot_path(name)
        return fs.is_file(path)

    # -----------------------------------------------------------------

    def has_misc_output(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_misc_path(name)
        return fs.is_directory(path) and not fs.is_empty(path)

    # -----------------------------------------------------------------

    def get_simulation_misc_fluxes_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.get_simulation_misc_path(name), "fluxes")

    # -----------------------------------------------------------------

    def has_misc_fluxes(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_misc_fluxes_path(name)
        return fs.is_directory(path) and not fs.is_empty(path)

    # -----------------------------------------------------------------

    def get_simulation_misc_image_fluxes_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.get_simulation_misc_path(name), "image fluxes")

    # -----------------------------------------------------------------

    def has_misc_image_fluxes(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_misc_image_fluxes_path(name)
        return fs.is_directory(path) and not fs.is_empty(path)

    # -----------------------------------------------------------------

    def get_simulation_misc_image_fluxes_images_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.get_simulation_misc_image_fluxes_path(name), "images")

    # -----------------------------------------------------------------

    def has_misc_image_fluxes_images(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_misc_image_fluxes_images_path(name)
        return fs.is_directory(path) and not fs.is_empty(path)

    # -----------------------------------------------------------------

    def get_simulation_misc_image_fluxes_images_earth_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.get_simulation_misc_image_fluxes_images_path(name), "earth")

    # -----------------------------------------------------------------

    def has_misc_image_fluxes_images_earth(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_misc_image_fluxes_images_earth_path(name)
        return fs.is_directory(path) and not fs.is_empty(path)

    # -----------------------------------------------------------------

    def get_simulation_misc_differences_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.get_simulation_misc_path(name), "differences.dat")

    # -----------------------------------------------------------------

    def has_misc_differences(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_misc_differences_path(name)
        return fs.is_file(path)

    # -----------------------------------------------------------------

    def get_simulation_misc_differences(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        from .modelanalyser import FluxDifferencesTable
        path = self.get_simulation_misc_differences_path(name)
        return FluxDifferencesTable.from_file(path)

    # -----------------------------------------------------------------

    def get_mock_sed_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Use images
        if self.use_images: return self.get_image_fluxes_path_for_simulation(simulation_name)

        # Regular
        else: return self.get_fluxes_path_for_simulation(simulation_name)

    # -----------------------------------------------------------------

    def has_mock_sed(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Use images
        if self.use_images: return self.has_image_fluxes_for_simulation(simulation_name)

        # Regular
        else: return self.has_fluxes_for_simulation(simulation_name)

    # -----------------------------------------------------------------

    def get_mock_sed(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return ObservedSED.from_file(self.get_mock_sed_path(simulation_name))

    # -----------------------------------------------------------------

    def get_image_fluxes_plot_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).single_images_fluxes_plot

    # -----------------------------------------------------------------

    def get_fluxes_plot_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).single_fluxes_plot

    # -----------------------------------------------------------------

    def has_image_fluxes_plot_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).has_single_images_fluxes_plot

    # -----------------------------------------------------------------

    def has_fluxes_plot_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).has_single_fluxes_plot

    # -----------------------------------------------------------------

    def get_mock_sed_plot_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Use images
        if self.use_images: return self.get_image_fluxes_plot_path_for_simulation(simulation_name)

        # Regular
        else: return self.get_fluxes_plot_path_for_simulation(simulation_name)

    # -----------------------------------------------------------------

    def has_mock_sed_plot(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Use images
        if self.use_images: return self.has_image_fluxes_plot_for_simulation(simulation_name)

        # Regular
        else: return self.has_fluxes_plot_for_simulation(simulation_name)

    # -----------------------------------------------------------------

    def get_simulation_sed_differences_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.get_simulation_misc_path(name), "differences.dat")

    # -----------------------------------------------------------------

    def has_sed_differences(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        path = self.get_simulation_sed_differences_path(name)
        return fs.is_file(path)

    # -----------------------------------------------------------------

    def get_simulation_sed_differences(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        from .modelanalyser import FluxDifferencesTable
        path = self.get_simulation_sed_differences_path(name)
        return FluxDifferencesTable.from_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_paths(self):

        """
        This function ...
        :return:
        """

        # Initialize
        all_cache_paths = set()

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            output_cached_path = self.get_output_cached_path(simulation_name)
            extraction_cached_path = self.get_extraction_cached_path(simulation_name)
            plotting_cached_path = self.get_plotting_cached_path(simulation_name)
            misc_cached_path = self.get_misc_cached_path(simulation_name)

            rel_output_path = fs.relative_to(self.get_simulation_output_path(simulation_name), self.modeling_path)
            rel_extraction_path = fs.relative_to(self.get_simulation_extract_path(simulation_name), self.modeling_path)
            rel_plotting_path = fs.relative_to(self.get_simulation_plot_path(simulation_name), self.modeling_path)
            rel_misc_path = fs.relative_to(self.get_simulation_misc_path(simulation_name), self.modeling_path)

            cache_paths = []
            if output_cached_path is not None: cache_paths.append(output_cached_path.split("/" + rel_output_path)[0])
            if extraction_cached_path is not None: cache_paths.append(extraction_cached_path.split("/" + rel_extraction_path)[0])
            if plotting_cached_path is not None: cache_paths.append(plotting_cached_path.split("/" + rel_plotting_path)[0])
            if misc_cached_path is not None: cache_paths.append(misc_cached_path.split("/" + rel_misc_path)[0])
            if len(cache_paths) == 0: continue
            cache_path = sequences.get_all_equal_value(cache_paths)

            # Add
            all_cache_paths.add(cache_path)

        # Return
        return list(all_cache_paths)

    # -----------------------------------------------------------------

    @property
    def ncache_paths(self):

        """
        This function ...
        :return:
        """

        return len(self.cache_paths)

    # -----------------------------------------------------------------

    @property
    def has_cache_paths(self):

        """
        This function ...
        :return:
        """

        return self.ncache_paths > 0

    # -----------------------------------------------------------------

    @property
    def has_single_cache_path(self):

        """
        This function ...
        :return:
        """

        return self.ncache_paths == 1

    # -----------------------------------------------------------------

    @property
    def single_cache_path(self):

        """
        This function ...
        :return:
        """

        if not self.has_cache_paths: raise ValueError("No cache paths")
        if not self.has_single_cache_path: raise ValueError("Not a single cache path: " + str(self.ncache_paths) + " different paths")
        return self.cache_paths[0]

    # -----------------------------------------------------------------

    def get_output_cached_filepath(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return fs.join(self.get_simulation_output_path(simulation_name), cached_filename)

    # -----------------------------------------------------------------

    def has_output_cached(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return fs.is_file(self.get_output_cached_filepath(simulation_name))

    # -----------------------------------------------------------------

    def get_output_cached_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if not self.has_output_cached(simulation_name): return None
        return fs.get_first_line(self.get_output_cached_filepath(simulation_name))

    # -----------------------------------------------------------------

    def get_extraction_cached_filepath(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return fs.join(self.get_simulation_extract_path(simulation_name), cached_filename)

    # -----------------------------------------------------------------

    def has_extraction_cached(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return fs.is_file(self.get_extraction_cached_filepath(simulation_name))

    # -----------------------------------------------------------------

    def get_extraction_cached_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if not self.has_extraction_cached(simulation_name): return None
        return fs.get_first_line(self.get_extraction_cached_filepath(simulation_name))

    # -----------------------------------------------------------------

    def get_plotting_cached_filepath(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return fs.join(self.get_simulation_plot_path(simulation_name), cached_filename)

    # -----------------------------------------------------------------

    def has_plotting_cached(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return fs.is_file(self.get_plotting_cached_filepath(simulation_name))

    # -----------------------------------------------------------------

    def get_plotting_cached_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if not self.has_plotting_cached(simulation_name): return None
        return fs.get_first_line(self.get_plotting_cached_filepath(simulation_name))

    # -----------------------------------------------------------------

    def get_misc_cached_filepath(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return fs.join(self.get_simulation_misc_path(simulation_name), cached_filename)

    # -----------------------------------------------------------------

    def has_misc_cached(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return fs.is_file(self.get_misc_cached_filepath(simulation_name))

    # -----------------------------------------------------------------

    def get_misc_cached_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if not self.has_misc_cached(simulation_name): return None
        return fs.get_first_line(self.get_misc_cached_filepath(simulation_name))

    # -----------------------------------------------------------------

    def is_retrieved(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Check whether the simulation output is present
        return self.has_simulation_output(name)

    # -----------------------------------------------------------------

    def is_analysed(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Check whether the chi squared is set in the chi squared table
        if self.has_chi_squared_table: return self.chi_squared_table.has_simulation(name)

        # Chi squared is not present!
        else:
            log.warning("Chi squared table for generation " + self.name + " is not present at '" + self.chi_squared_table_path + ": determining whether simulation '" + name + "' is analysed by looking at the miscellaneous output")
            return self.has_misc_output(name)

    # -----------------------------------------------------------------

    @lazyproperty
    def host_ids(self):

        """
        This function ...
        :return:
        """

        return self.assignment_table.unique_host_ids

    # -----------------------------------------------------------------

    def get_simulation_names_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.assignment_table.simulations_for_remote(host_id)

    # -----------------------------------------------------------------

    def get_simulation_ids_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.assignment_table.ids_for_remote(host_id)

    # -----------------------------------------------------------------

    def get_simulation_paths_for_host(self, host_id, as_dict=True, id_or_name="name", not_exist="error"):

        """
        This function ...
        :param host_id:
        :param as_dict:
        :param id_or_name:
        :param not_exist:
        :return:
        """

        paths = OrderedDict()
        simulation_ids = self.assignment_table.ids_for_remote(host_id)
        for simulation_id in simulation_ids:

            # Get the filepath
            filepath = get_simulation_path_for_host(host_id, simulation_id)
            if not fs.is_file(filepath):
                if not_exist == "error": raise ValueError("Simulation file does not exist")
                elif not_exist == "ignore": continue
                elif not_exist == "none": filepath = None
                elif not_exist == "pass": pass
                else: raise ValueError("Invalid input for 'not_exist'")

            # Set key
            if id_or_name == "id": key = simulation_id
            else: key = self.assignment_table.get_simulation_name_for_id(simulation_id)

            # Add
            paths[key] = filepath

        # Return
        if as_dict: return paths
        else: return paths.values()

    # -----------------------------------------------------------------

    def get_simulations_for_host(self, host_id, as_dict=False):

        """
        This function ...
        :param host_id:
        :param as_dict:
        :return:
        """

        simulations = OrderedDict()

        #simulation_ids = self.assignment_table.ids_for_remote(host_id)
        simulation_names, simulation_ids = self.assignment_table.names_and_ids_for_remote(host_id)

        # Loop over the simulation IDs
        for simulation_name, simulation_id in zip(simulation_names, simulation_ids):

            #print(host_id, simulation_id, simulation_name)
            if simulation_id is None:

                if simulation_name is not None:
                    simulation = find_simulation_for_host(simulation_name, host_id)
                    if simulation is not None: simulations[simulation_name] = simulation

            elif has_simulation_for_host(host_id, simulation_id):

                simulation = get_simulation_for_host(host_id, simulation_id)
                simulations[simulation.name] = simulation

            else: log.warning("Simulation with ID '" + str(simulation_id) + "' is not found in the directory remote host '" + host_id + "'")

        # Return the simulations
        if as_dict: return simulations
        else: return simulations.values()

    # -----------------------------------------------------------------

    def in_assignment(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.assignment_table.has_simulation(name)

    # -----------------------------------------------------------------

    def has_simulation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        try: host_id = self.get_host_id(name)
        except MissingSimulation: return False
        if host_id is None: return False

        simulation_id = self.get_simulation_id(name)
        if simulation_id is None: return False
        #print(name, host_id, simulation_id)

        if not has_simulation_for_host(host_id, simulation_id): return False
        else:
            simulation = get_simulation_for_host(host_id, simulation_id)
            return simulation.name == name

    # -----------------------------------------------------------------

    def get_simulation_filepath(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get host ID and simulation ID
        host_id = self.get_host_id(name)
        simulation_id = self.get_simulation_id(name)

        # Return the path
        return get_simulation_path_for_host(host_id, simulation_id)

    # -----------------------------------------------------------------

    def get_simulations_basic_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Initialize list
        simulations = []

        simulation_names = self.assignment_table.simulations_for_remote(host_id)
        for simulation_name in simulation_names:

            # Create the simulation object
            simulation = self.get_simulation_basic(simulation_name, host_id=host_id)

            # Add the simulation
            simulations.append(simulation)

        # Return the simulation objects
        return simulations

    # -----------------------------------------------------------------

    def get_host_id(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        #if not self.has_assignment_table: return None
        #return self.assignment_table.get_host_id_for_simulation(name)

        if self.has_assignment_table: return self._get_host_id_from_assignment(name)
        elif self.has_jobscript(name): return self._get_host_id_from_job_script(name)
        elif self.has_logfile(name): return self._get_host_id_from_logfile(name)
        else: return None

    # -----------------------------------------------------------------

    def _get_host_id_from_assignment(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.assignment_table.get_host_id_for_simulation(name)

    # -----------------------------------------------------------------

    def _get_host_id_from_job_script(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # TOO SLOW
        # Get the jobscript
        #jobscript = self.get_simulation_jobscript(name)
        # Return the host ID
        #return jobscript.host_id

        # Get jobscript path
        filepath = self.get_simulation_jobscript_path(name)
        return get_host_id_from_jobscript_file(filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def all_hosts(self):

        """
        This function ...
        :return:
        """

        from ...core.remote.host import find_hosts
        return find_hosts(as_dict=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def all_host_ids(self):

        """
        Thisf unction ...
        :return:
        """

        return self.all_hosts.keys()

    # -----------------------------------------------------------------

    def _get_host_id_from_logfile(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the logfile
        logfile = self.get_simulation_logfile(name)

        # Return the host ID
        host_name = logfile.host
        if host_name in self.all_host_ids: return host_name
        else:

            for host_id in self.all_host_ids:

                host = self.all_hosts[host_id]
                if not host.has_clusters: continue
                if host_name in host.cluster_names: return host_id

            # Not recognized
            else:
                log.warning("Could not recognized host name '" + host_name + "'")
                return None

    # -----------------------------------------------------------------

    def get_cluster_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        #if not self.has_assignment_table: return None
        #return self.assignment_table.get_cluster_name_for_simulation(name)

        if self.has_assignment_table: return self._get_cluster_name_from_assignment(name)
        elif self.has_logfile(name): return self._get_cluster_name_from_logfile(name)
        else: return None

    # -----------------------------------------------------------------

    def _get_cluster_name_from_assignment(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.assignment_table.get_cluster_name_for_simulation(name)

    # -----------------------------------------------------------------

    def _get_cluster_name_from_logfile(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the logfile
        logfile = self.get_simulation_logfile(name)

        # Return the host ID
        host_name = logfile.host
        if host_name in self.all_host_ids:
            return None
            #host = self.all_hosts[host_name]
            #if host.has_clusters: return None # cluster name not known
            #else: return

        # Must be a cluster name, not recognized as host ID
        else: return host_name

    # -----------------------------------------------------------------

    def get_simulation_id(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        #if not self.has_assignment_table: return None
        #return self.assignment_table.get_simulation_id_for_simulation(name)

        if self.has_assignment_table: return self._get_simulation_id_from_assignment(name)
        else: return None

    # -----------------------------------------------------------------

    def _get_simulation_id_from_assignment(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.assignment_table.get_simulation_id_for_simulation(name)

    # -----------------------------------------------------------------

    def get_simulation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get host ID and simulation ID
        host_id = self.get_host_id(name)
        simulation_id = self.get_simulation_id(name)
        #print(host_id, simulation_id)

        # Create the simulation object (standard or remote)
        if has_simulation_for_host(host_id, simulation_id): simulation = get_simulation_for_host(self.get_host_id(name), self.get_simulation_id(name))
        else: simulation = self.get_simulation_basic(name, host_id=host_id, simulation_id=simulation_id)

        # Load and return the simulation
        #simulation = get_simulation_for_host(self.get_host_id(name), self.get_simulation_id(name))

        # Return the simulation
        if simulation.name != name: raise RuntimeError("Wrong simulation!")
        return simulation

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation_output(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_simulation_or_basic(name).output

    # -----------------------------------------------------------------

    @memoize_method
    def get_extraction_output(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_simulation_or_basic(name).extraction_output

    # -----------------------------------------------------------------

    @memoize_method
    def get_plotting_output(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_simulation_or_basic(name).plotting_output

    # -----------------------------------------------------------------

    @memoize_method
    def get_misc_output(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_simulation_or_basic(name).misc_output

    # -----------------------------------------------------------------

    def has_fluxes_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).has_single_fluxes

    # -----------------------------------------------------------------

    def get_fluxes_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).single_fluxes

    # -----------------------------------------------------------------

    def get_fluxes_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return ObservedSED.from_file(self.get_fluxes_path_for_simulation(simulation_name))

    # -----------------------------------------------------------------

    def has_image_fluxes_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).has_single_image_fluxes

    # -----------------------------------------------------------------

    def get_image_fluxes_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).single_image_fluxes

    # -----------------------------------------------------------------

    def get_image_fluxes_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return ObservedSED.from_file(self.get_image_fluxes_path_for_simulation(simulation_name))

    # -----------------------------------------------------------------

    def has_images_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).has_images_for_fluxes

    # -----------------------------------------------------------------

    def get_image_paths_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).images_for_fluxes

    # -----------------------------------------------------------------

    def get_images_for_simulation(self, simulation_name):

        """
        This function ....
        :param simulation_name:
        :return:
        """

        # Initialize a dictionary for the images
        images = OrderedDict()

        # Loop over the image paths
        filepaths = self.get_image_paths_for_simulation(simulation_name)
        for filepath in filepaths:

            # Get image name
            name = fs.strip_extension(fs.name(filepath))

            # Add to the dictionary
            images[name] = Frame.from_file(filepath)

        # Return the images
        return images

    # -----------------------------------------------------------------

    @property
    def wavelength_grids_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.wavelength_grids_path

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.wavelength_grids_path, self.wavelength_grid_name + ".dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return WavelengthGrid.from_skirt_input(self.wavelength_grid_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def use_file_tree_dust_grid(self):

        """
        This function ...
        :return:
        """

        from ...core.prep.smile import SKIRTSmileSchema
        smile = SKIRTSmileSchema()
        return smile.supports_file_tree_grids and self.model_representation.has_dust_grid_tree

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_input(self):

        """
        This function ...
        :return:
        """

        # Initilize
        simulation_input = SimulationInput()

        # Set the paths to the input maps
        for name in self.fitting_run.input_map_paths:
            path = self.fitting_run.input_map_paths[name]
            simulation_input.add_file(path, name)

        # DETERMINE AND SET THE PATH TO THE APPROPRIATE DUST GRID TREE FILE
        if self.use_file_tree_dust_grid: simulation_input.add_file(self.model_representation.dust_grid_tree_path)

        # Determine and set the path to the appropriate wavelength grid file
        simulation_input.add_file(self.wavelength_grid_path, name="wavelengths.txt")

        # Return the object
        return simulation_input

    # -----------------------------------------------------------------

    def get_simulation_basic(self, name, host_id=None, simulation_id=None):

        """
        This function ...
        :param name:
        :param host_id:
        :param simulation_id:
        :return:
        """

        # Get properties
        ski_path = self.get_simulation_ski_path(name)
        out_path = self.get_simulation_output_path(name)

        # Create the simulation object
        if host_id is not None or simulation_id is not None: simulation = RemoteSimulation(ski_path, self.simulation_input, out_path, host_id=host_id, id=simulation_id, name=name)
        else: simulation = SkirtSimulation(inpath=self.simulation_input, outpath=out_path, ski_path=ski_path, name=name)

        # Set the analysis paths
        simulation.extraction_path = self.get_simulation_extract_path(name)
        simulation.plotting_path = self.get_simulation_plot_path(name)
        simulation.misc_path = self.get_simulation_misc_path(name)

        # Set the modeling path
        simulation.analysis.modeling_path = self.modeling_path

        # Return the simulation
        return simulation

    # -----------------------------------------------------------------

    def get_simulation_or_basic(self, name, host_id=None, simulation_id=None):

        """
        This function ...
        :param name:
        :param host_id:
        :param simulation_id:
        :return:
        """

        if self.has_simulation(name): return self.get_simulation(name)
        else: return self.get_simulation_basic(name, host_id=host_id, simulation_id=simulation_id)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulations(self):

        """
        This function ...
        :return:
        """

        sims = []
        for host_id in self.host_ids:
            simulations = self.get_simulations_for_host(host_id)
            sims.extend(simulations)
        return sims

    # -----------------------------------------------------------------

    @lazyproperty
    def simulations_basic(self):

        """
        This function ...
        :return:
        """

        all_simulations = []

        # DOESN'T WORK WELL WHEN ASSIGNMENT TABLE IS INCORRECT?
        # If there is an assignment table
        if self.has_assignment_table:

             # Loop over the host IDs in the assignment tables
             for host_id in self.host_ids:

                 # Create simulations and add them
                 simulations_host = self.get_simulations_basic_for_host(host_id)
                 all_simulations.extend(simulations_host)

        # No assignment tables
        else:

            # Loop over all simulation names
            for simulation_name in self.simulation_names:

                # Check whether host ID and simulation ID can be determined
                host_id = self.get_host_id(simulation_name)
                simulation_id = self.get_simulation_id(simulation_name)

                # Create simulation and add it
                simulation = self.get_simulation_basic(simulation_name, host_id=host_id, simulation_id=simulation_id)
                all_simulations.append(simulation)

        # Return the simulation objects
        return all_simulations

    # -----------------------------------------------------------------

    @lazyproperty
    def simulations_or_basic(self):

        """
        This function ...
        :return:
        """

        #simulations = []
        simulations = OrderedDict()

        # assignment table present?
        nassignment_simulations = None
        if self.has_assignment_table:

            # Loop over the host IDs in the assignment tables
            for host_id in self.host_ids: simulations.update(self.get_simulations_for_host(host_id, as_dict=True))
                #simulations = self.get_simulations_for_host(host_id)
                #simulations.extend(simulations)

            nassignment_simulations = len(self.assignment_table)

        #else: print("no assignment")

        # Check number of simulations found through the assignment file
        nsimulations = len(simulations)
        #print(nsimulations, self.nsimulations, nassignment_simulations, self.nsimulation_directories)

        if nsimulations == self.nsimulations:

            if self.nsimulation_directories == self.nsimulations: return simulations.values()
            elif self.nsimulation_directories == nassignment_simulations:
                # There are probably more simulations than the info says (maybe because of expansion)
                log.warning("Fixing the number of simulations in the generation info (" + str(nassignment_simulations) + " instead of " + str(nsimulations) + ") ...")
                self.info.nsimulations = nassignment_simulations
                self.info.save()
                self._add_missing_simulations(simulations, from_directory_names=True)
                return simulations.values()

            else: return simulations.values()

        elif nsimulations > self.nsimulations:

            log.warning("The number of simulations in the generation info is lower than expected (" + str(self.nsimulations) + " instead of " + str(nsimulations) + ")")
            return simulations.values()

        else: # less simulations found or none yet (no assignment)

            if nsimulations > 0: log.warning("Only " + str(nsimulations) + " simulation objects were found, " + str(self.nsimulations) + " are expected. Creating additional simulation objects ...")

            # Add missing simulations
            if nassignment_simulations == self.nsimulation_directories: self._add_missing_simulations(simulations, from_directory_names=True)
            else: self._add_missing_simulations(simulations)

            # Return the simulation objects
            return simulations.values()

    # -----------------------------------------------------------------

    def _add_missing_simulations(self, simulations, from_directory_names=False):

        """
        This function ...
        :param simulations:
        :return:
        """

        if from_directory_names: all_names = self.simulation_directory_names
        else: all_names = self.simulation_names

        #print(all_names)

        # Loop over all simulation names, add missing simulations
        for simulation_name in all_names:
            if simulation_name in simulations: continue

            # Check whether host ID and simulation ID can be determined
            host_id = self.get_host_id(simulation_name)
            simulation_id = self.get_simulation_id(simulation_name)

            # Try to find simulation
            simulation = None
            if host_id is None: simulation = find_simulation(simulation_name, self.all_host_ids)
            elif simulation_id is None: simulation = find_simulation_for_host(simulation_name, host_id)
            elif has_simulation_for_host(host_id, simulation_id): simulation = get_simulation_for_host(host_id, simulation_id)

            # Create basic simulation if needed
            if simulation is None: simulation = self.get_simulation_basic(simulation_name, host_id=host_id, simulation_id=simulation_id)

            # Add
            # simulations.append(simulation)
            simulations[simulation_name] = simulation

    # -----------------------------------------------------------------

    @lazyproperty
    def retrieved_simulations(self):

        """
        This function ...
        :return:
        """

        return [simulation for simulation in self.simulations if simulation.retrieved]

    # -----------------------------------------------------------------

    @lazyproperty
    def retrieved_simulations_basic(self):

        """
        This function ...
        :return:
        """

        return [simulation for simulation in self.simulations_basic if self.is_retrieved(simulation.name)]

    # -----------------------------------------------------------------

    @lazyproperty
    def analysed_simulations(self):

        """
        This function ...
        :return:
        """

        return [simulation for simulation in self.simulations if simulation.analysed]

    # -----------------------------------------------------------------

    @lazyproperty
    def analysed_simulations_basic(self):

        """
        This function ...
        :return:
        """

        return [simulation for simulation in self.simulations_basic if self.is_analysed(simulation.name)]

    # -----------------------------------------------------------------

    @lazyproperty
    def nretrieved_simulations(self):

        """
        This function ...
        :return:
        """

        return len(self.retrieved_simulations)

    # -----------------------------------------------------------------

    @lazyproperty
    def nretrieved_simulations_basic(self):

        """
        This function ...
        :return:
        """

        return len(self.retrieved_simulations_basic)

    # -----------------------------------------------------------------

    @lazyproperty
    def nanalysed_simulations(self):

        """
        This function ...
        :return:
        """

        return len(self.analysed_simulations)

    # -----------------------------------------------------------------

    @lazyproperty
    def nanalysed_simulations_basic(self):

        """
        This function ...
        :return:
        """

        return len(self.analysed_simulations_basic)

    # -----------------------------------------------------------------

    @property
    def has_queues(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.path, extension="dat", startswith="queue_")

    # -----------------------------------------------------------------

    @property
    def queue_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.path, extension="dat", startswith="queue_")

    # -----------------------------------------------------------------

    @property
    def local_queue_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "queue_local.dat")

    # -----------------------------------------------------------------

    @property
    def has_local_queue(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.is_file(self.local_queue_path)

    # -----------------------------------------------------------------

    @property
    def has_remote_queues(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.has_files_in_path(self.path, extension="dat", startswith="queue_", exact_not_name="queue_local")

    # -----------------------------------------------------------------

    @property
    def remote_queue_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.path, extension="dat", startswith="queue_", exact_not_name="queue_local")

    # -----------------------------------------------------------------

    @property
    def remote_queue_host_ids(self):

        """
        This function ...
        :return:
        """

        return [name.split("queue_")[1] for name in fs.files_in_path(self.path, extension="dat", startswith="queue_", exact_not_name="queue_local", returns="name")]

    # -----------------------------------------------------------------

    @property
    def index(self):

        """
        This function ...
        :return:
        """

        return self.info.index

    # -----------------------------------------------------------------

    @property
    def method(self):

        """
        This function ...
        :return:
        """

        return self.info.method

    # -----------------------------------------------------------------

    @property
    def wavelength_grid_name(self):

        """
        This function ...
        :return:
        """

        return self.info.wavelength_grid_name

    # -----------------------------------------------------------------

    @property
    def model_representation_name(self):

        """
        This function ...
        :return:
        """

        return self.info.model_representation_name

    # -----------------------------------------------------------------

    @lazyproperty
    def model_representation(self):

        """
        This function ...
        :return:
        """

        name = self.model_representation_name
        return get_representation(self.modeling_path, name)

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        THis function ...
        :return:
        """

        return self.model_representation.model_name

    # -----------------------------------------------------------------

    @lazyproperty
    def model_definition(self):

        """
        THis function ...
        :return:
        """

        name = self.model_name
        return get_model_definition(self.modeling_path, name)

    # -----------------------------------------------------------------

    @property
    def nsimulations(self):

        """
        THis function ...
        :return:
        """

        return self.info.nsimulations

    # -----------------------------------------------------------------

    @property
    def nsimulation_directories(self):

        """
        This function ...
        :return:
        """

        return len(self.simulation_directory_names)

    # -----------------------------------------------------------------

    @property
    def npackages(self):

        """
        This function ...
        :return:
        """

        return self.info.npackages

    # -----------------------------------------------------------------

    @property
    def selfabsorption(self):

        """
        This function ...
        :return:
        """

        return self.info.selfabsorption

    # -----------------------------------------------------------------

    @property
    def transient_heating(self):

        """
        THis function ...
        :return:
        """

        return self.info.transient_heating

    # -----------------------------------------------------------------

    @property
    def spectral_convolution(self):

        """
        This function ...
        :return:
        """

        return self.info.spectral_convolution

    # -----------------------------------------------------------------

    @property
    def use_images(self):

        """
        This function ...
        :return:
        """

        return self.info.use_images

    # -----------------------------------------------------------------

    @property
    def individuals_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "individuals.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def individuals_table(self):

        """
        This function ...
        :return:
        """

        return IndividualsTable.from_file(self.individuals_table_path)

    # -----------------------------------------------------------------

    @property
    def parameters_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "parameters.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def parameters_table(self):

        """
        This function ...
        :return:
        """

        return ParametersTable.from_file(self.parameters_table_path)

    # -----------------------------------------------------------------

    def get_parameter_values_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get a dictionary with the parameter values for this simulation
        return self.parameters_table.parameter_values_for_simulation(simulation_name)

    # -----------------------------------------------------------------

    def get_simulation_name_for_parameter_values(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        return self.parameters_table.simulation_for_parameter_values(values)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_parameter_values(self):
        return self.parameters_table.unique_parameter_values

    # -----------------------------------------------------------------

    @property
    def chi_squared_table_path(self):
        return fs.join(self.path, "chi_squared.dat")

    # -----------------------------------------------------------------

    @property
    def has_chi_squared_table(self):
        return fs.is_file(self.chi_squared_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def chi_squared_table(self):
        return ChiSquaredTable.from_file(self.chi_squared_table_path)

    # -----------------------------------------------------------------

    @property
    def best_simulation_name(self):
        return self.chi_squared_table.best_simulation_name

    # -----------------------------------------------------------------

    @property
    def worst_simulation_name(self):
        return self.chi_squared_table.worst_simulation_name

    # -----------------------------------------------------------------

    def get_chi_squared(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.chi_squared_table.chi_squared_for(simulation_name)

    # -----------------------------------------------------------------

    @property
    def newborns_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "newborns.dat")

    # -----------------------------------------------------------------

    @property
    def newborns(self):

        """
        This function ...
        :return:
        """

        if fs.is_file(self.newborns_path): return load_population(self.newborns_path)
        else: return None

    # -----------------------------------------------------------------

    @property
    def parents_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "parents.dat")

    # -----------------------------------------------------------------

    @property
    def parents(self):

        """
        This function ...
        :return:
        """

        return load_population(self.parents_path)

    # -----------------------------------------------------------------

    @property
    def scores_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "scores")

    # -----------------------------------------------------------------

    @property
    def scores_table(self):

        """
        This function ...
        :return:
        """

        return ScoresTable.from_file(self.scores_table_path)

    # -----------------------------------------------------------------

    @property
    def elitism_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "elitism.dat")

    # -----------------------------------------------------------------

    @property
    def elitism_table(self):

        """
        This function ...
        :return:
        """

        return ElitismTable.from_file(self.elitism_table_path)

    # -----------------------------------------------------------------

    @property
    def recurrence_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "recurrence.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def recurrence_table(self):

        """
        This function ...
        :return:
        """

        return RecurrenceTable.from_file(self.recurrence_path)

    # -----------------------------------------------------------------

    @property
    def crossover_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "crossover.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def crossover_table(self):

        """
        This function ...
        :return:
        """

        return CrossoverTable.from_file(self.crossover_table_path)

    # -----------------------------------------------------------------

    @property
    def optimizer_input_path(self):
        return fs.create_directory_in(self.path, "input")

    # -----------------------------------------------------------------

    @property
    def has_optimizer_input(self):
        return fs.is_directory(self.optimizer_input_path) and not fs.is_empty(self.optimizer_input_path)

    # -----------------------------------------------------------------

    @property
    def optimizer_input(self):
        return load_input(self.optimizer_input_path)

    # -----------------------------------------------------------------

    @property
    def parameter_minima_scalar(self):
        return self.optimizer_input["minima"]

    # -----------------------------------------------------------------

    @property
    def parameter_maxima_scalar(self):
        return self.optimizer_input["maxima"]

    # -----------------------------------------------------------------

    def unit_for_parameter(self, label):
        return self.parameters_table.unit_for(label)

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):
        ## SORTED JUST AS IN FITTINGRUN !!
        return sorted(self.parameters_table.parameter_labels)

    # -----------------------------------------------------------------

    def index_for_parameter(self, label):
        return self.parameter_labels.index(label)

    # -----------------------------------------------------------------

    @property
    def parameter_units(self):

        """
        This function ...
        :return:
        """

        units = dict()

        for label in self.parameter_labels:

            unit = self.unit_for_parameter(label)
            units[label] = unit

        # Return the units
        return units

    # -----------------------------------------------------------------

    @property
    def parameter_minima(self):

        """
        This function ...
        :return:
        """

        minima = dict()

        for minimum, label in zip(self.parameter_minima_scalar, self.parameter_labels):

            unit = self.unit_for_parameter(label)
            if unit is not None: minimum *= unit

            minima[label] = minimum

        # Return
        return minima

    # -----------------------------------------------------------------

    @property
    def parameter_maxima(self):

        """
        This function ...
        :return:
        """

        maxima = dict()

        for maximum, label in zip(self.parameter_maxima_scalar, self.parameter_labels):

            unit = self.unit_for_parameter(label)
            if unit is not None: maximum *= unit

            maxima[label] = maximum

        # Return
        return maxima

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_ranges(self):

        """
        This function ...
        :return:
        """

        # PREVIOUS
        #return self.parameter_ranges_optimizer()

        # Initialize dict
        ranges = OrderedDict()

        # Loop over the parameters
        for label in self.parameter_labels:

            # Get unique values
            unique_values = self.unique_parameter_values[label]

            # Quantity
            unit = self.unit_for_parameter(label)
            if unit is not None: _range = QuantityRange.limits(unique_values)

            # Scalar int
            elif self.parameter_base_types[label] == "integer": _range = IntegerRange.limits(unique_values)

            # Scalar real
            elif self.parameter_base_types[label] == "real": _range = RealRange.limits(unique_values)

            # Invalid
            else: raise ValueError("Unrecognized base type: " + self.parameter_base_types[label])

            # Add the range
            ranges[label] = _range

        # Return
        return ranges

    # -----------------------------------------------------------------

    @property
    def parameter_ranges_optimizer(self):

        """
        This function ...
        :return:
        """

        ranges = dict()

        minima = self.parameter_minima
        maxima = self.parameter_maxima

        for label in self.parameter_labels:

            # Get min and max value
            min_value = minima[label]
            max_value = maxima[label]

            unit = self.unit_for_parameter(label)

            # Quantity
            if unit is not None: range = QuantityRange(min_value, max_value)

            # Scalar value
            else:

                base_type = self.parameter_base_types[label]
                if base_type == "integer": range = IntegerRange(min_value, max_value)
                elif base_type == "real": range = RealRange(min_value, max_value)
                else: raise ValueError("Unrecognized base type: " + base_type)

            # Add the range
            ranges[label] = range

        # Return
        return ranges

    # -----------------------------------------------------------------

    @property
    def parameter_scales(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_input["scales"]

    # -----------------------------------------------------------------

    @property
    def parameter_ndigits(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_input["ndigits"]

    # -----------------------------------------------------------------

    @property
    def ndigits_list(self):

        """
        This function ...
        :return:
        """

        return self.parameter_ndigits

    # -----------------------------------------------------------------

    @property
    def parameter_nbits(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_input["nbits"]

    # -----------------------------------------------------------------

    @property
    def nbits_list(self):

        """
        This function ...
        :return:
        """

        return self.parameter_nbits

    # -----------------------------------------------------------------

    @property
    def total_nbits(self):

        """
        This function ...
        :return:
        """

        return sum(self.nbits_list)

    # -----------------------------------------------------------------

    @property
    def nparameters(self):

        """
        This function ...
        :return:
        """

        return len(self.parameter_labels)

    # -----------------------------------------------------------------

    @property
    def genome_size(self):

        """
        This function ...
        :return:
        """

        # 1D
        if self.is_1d_genome:

            # 1D
            if self.list_genome: return self.nparameters
            elif self.binary_string_genome: return self.total_nbits
            else: raise ValueError("Invalid genome type")

        # 2D
        #elif self.is_2d_genome: raise ValueError("We should not get here")
        raise ValueError("We should not get here")

    # -----------------------------------------------------------------

    @property
    def elitism(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.elitism

    # -----------------------------------------------------------------

    @property
    def check_recurrence(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.check_recurrence

    # -----------------------------------------------------------------

    @property
    def gray_code(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.gray_code

    # -----------------------------------------------------------------

    @property
    def genome_type(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.genome_type

    # -----------------------------------------------------------------

    @property
    def list_genome(self):

        """
        This function ...
        :return:
        """

        return self.genome_type == "list"

    # -----------------------------------------------------------------

    @property
    def binary_string_genome(self):

        """
        This function ...
        :return:
        """

        return self.genome_type == "binary_string"

    # -----------------------------------------------------------------

    @property
    def genome_class(self):

        """
        This function ...
        :return:
        """

        #return genomes[self.genome_dimension][self.genome_type]
        return get_genome_class(self.genome_dimension, self.genome_type)

    # -----------------------------------------------------------------

    ## FIXED FOR MODELLING
    @property
    def genome_dimension(self):

        """
        This function ...
        :return:
        """

        return 1

    # -----------------------------------------------------------------

    ## FIXED FOR MODELLING
    @property
    def is_1d_genome(self):

        """
        This fucntion ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    ## FIXED FOR MODELLING
    @property
    def is_2d_genome(self):

        """
        This function ...
        :return:
        """

        return False

    # -----------------------------------------------------------------

    ## FIXED FOR MODELLING
    @property
    def min_or_max(self):

        """
        This function ...
        :return:
        """

        return "minimize"

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        from .run import FittingRun
        return FittingRun.from_path(self.fitting_run_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def generations_probabilities_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.fitting_run.prob_path, "generations")

    # -----------------------------------------------------------------

    @lazyproperty
    def probabilities_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.generations_probabilities_path, self.name)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_probabilities_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.probabilities_path, "models.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def model_probabilities_table(self):

        """
        This function ...
        :return:
        """

        return ModelProbabilitiesTable.from_file(self.model_probabilities_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def most_probable_model(self):

        """
        This function ...
        :return:
        """

        return self.model_probabilities_table.most_probable_simulation

    # -----------------------------------------------------------------

    def get_parameter_probabilities_table_path(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return fs.join(self.probabilities_path, label + ".dat")

    # -----------------------------------------------------------------

    def has_parameter_probabilities_table(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        path = self.get_parameter_probabilities_table_path(label)
        return fs.is_file(path)

    # -----------------------------------------------------------------

    def get_parameter_probabilities_table(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        path = self.get_parameter_probabilities_table_path(label)
        return ParameterProbabilitiesTable.from_file(path)

    # -----------------------------------------------------------------

    def get_most_probable_parameter_value(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        table = self.get_parameter_probabilities_table(label)
        return table.most_probable_value

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.fitting_run_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_base_types(self):

        """
        This function ...
        :return:
        """

        types = dict()

        ranges = self.parameter_ranges

        for label in self.parameter_labels:

            the_range = ranges[label]

            # Check the type of range
            if isinstance(the_range, IntegerRange): base_type = "integer"
            elif isinstance(the_range, RealRange): base_type = "real"
            elif isinstance(the_range, QuantityRange): base_type = "real" # quantity -> 'real' base type
            else: raise ValueError("Invalid parameter range")

            # Set the type
            types[label] = base_type

        # Return
        return types

    # -----------------------------------------------------------------

    @lazyproperty
    def single_parameter_base_type(self):

        """
        This function ...
        :return:
        """

        assert sequences.all_equal(self.parameter_base_types.values())
        return self.parameter_base_types[self.parameter_base_types.keys()[0]]

    # -----------------------------------------------------------------

    ## FIXED FOR MODELLING
    @property
    def heterogeneous(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def nindividuals_per_generation(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.nindividuals

    # -----------------------------------------------------------------

    @property
    def mutation_rate(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.mutation_rate

    # -----------------------------------------------------------------

    # DEFINED ABOVE ALREADY
    # @property
    # def nparameters(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     return self.optimizer_config.nparameters

    # -----------------------------------------------------------------

    @property
    def crossover_rate(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.crossover_rate

    # -----------------------------------------------------------------

    @property
    def mutation_method(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.mutation_method

    # -----------------------------------------------------------------

    @property
    def binary_mutation_method(self):

        """
        THis function ...
        :return:
        """

        return self.optimizer_config.binary_mutation_method

    # -----------------------------------------------------------------

    @property
    def mutation_function(self):

        """
        This function ...
        :return:
        """

        # List or binary
        if self.list_genome: return get_mutator(self.genome_type, self.genome_dimension, self.mutation_method, parameter_type=self.single_parameter_base_type, hetero=self.heterogeneous)
        elif self.binary_string_genome: return get_mutator(self.genome_type, self.genome_dimension, self.binary_mutation_method)
        else: raise ValueError("Invalid genome type")

    # -----------------------------------------------------------------

    @property
    def scaling_method(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.scaling_method

    # -----------------------------------------------------------------

    @property
    def scaling_function(self):

        """
        This function ...
        :return:
        """

        return get_scaling(self.scaling_method)

    # -----------------------------------------------------------------

    @property
    def selection_method(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.selector_method

    # -----------------------------------------------------------------

    @property
    def selection_function(self):

        """
        This function ...
        :return:
        """

        return get_selector(self.selection_method)

    # -----------------------------------------------------------------

    @property
    def initializator(self):

        """
        This function ...
        :return:
        """

        return get_initializator(self.genome_type, parameter_type=self.single_parameter_base_type, hetero=self.heterogeneous)

    # -----------------------------------------------------------------

    @property
    def crossover_method(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.crossover_method

    # -----------------------------------------------------------------

    @property
    def crossover_function(self):

        """
        This function ...
        :return:
        """

        return get_crossover(self.genome_type, self.genome_dimension, self.crossover_method, genome_size=self.genome_size)

    # -----------------------------------------------------------------

    @property
    def crossover_origins_function(self):

        """
        This function ...
        :return:
        """

        return get_crossover_origins(self.genome_type, self.genome_dimension, self.crossover_method)

    # -----------------------------------------------------------------

    @property
    def engine_path(self):

        """
        THis function ...
        :return:
        """

        return fs.join(self.path, "engine.pickle")

    # -----------------------------------------------------------------

    @property
    def prng_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "prng.pickle")

    # -----------------------------------------------------------------

    @property
    def optimizer_config_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "optimizer.cfg")

    # -----------------------------------------------------------------

    @property
    def optimizer_config(self):

        """
        This function ...
        :return:
        """

        return Configuration.from_file(self.optimizer_config_path)

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

    @property
    def timing_table(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.timing_table

    # -----------------------------------------------------------------

    @property
    def memory_table(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.memory_table

    # -----------------------------------------------------------------

    def has_timing(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.timing_table.has_simulation(simulation_name)

    # -----------------------------------------------------------------

    def has_memory(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.memory_table.has_simulation(simulation_name)

    # -----------------------------------------------------------------

    def find_simulation_for_assignment(self, simulation_name, find_remotes, remotes=None, produce_missing=False,
                                       correct_paths=False, confirm_correction=False):

        """
        This function ...
        :param simulation_name:
        :param find_remotes:
        :param remotes:
        :param produce_missing:
        :param correct_paths:
        :param confirm_correction:
        :return:
        """

        # Give warning
        log.warning("Looking for the simulation '" + simulation_name + "' on any remote host ...")
        simulation = find_simulation(simulation_name, find_remotes)

        # Simulation is not found
        if simulation is None:

            # Produce missing simulation file?
            if produce_missing:

                # Debugging
                log.warning("Simulation was not found, producing simulation file ...")

                # raise RuntimeError("Cannot produce missing simulations when ")
                host_id = prompt_string("host_id", "string", "specify the remote host ID for simulation '" + simulation_name + "'")

                # Produce
                self.produce_missing_simulation_file(simulation_name, host_id=host_id, remotes=remotes,
                                                     correct_paths=correct_paths,
                                                     confirm_correction=confirm_correction)

                # Set flag
                load_simulation = True

            # Don't produce: give error for missing simulation file
            else: raise RuntimeError("Simulation '" + simulation_name + "' was not found")

            changed_assignment = False

        # Simulation is found
        else:

            actual_host_id = simulation.host_id
            actual_id = simulation.id
            cluster_name = simulation.cluster_name
            log.warning("Simulation identified as '" + str(actual_id) + "' for host '" + actual_host_id + "'")

            # Add to assignment scheme
            log.warning("Adding simulation to the assignment table ...")
            self.add_simulation_to_assignment(simulation_name, host_id=actual_host_id, cluster_name=cluster_name,
                                              simulation_id=actual_id)

            changed_assignment = True
            load_simulation = True

        # Return flags
        return changed_assignment, load_simulation

    # -----------------------------------------------------------------

    def get_status(self, remotes=None, lazy=False, find_simulations=False, find_remotes=None, produce_missing=False,
                   retrieve=False, check_paths=False, correct_paths=False, confirm_correction=False, fix_success=True):

        """
        This function gets the status of the simulations
        :param remotes: must be a SKIRT remotes ensemble object, or None (with reduced functionality)
        :param lazy:
        :param find_simulations:
        :param find_remotes:
        :param produce_missing:
        :param retrieve:
        :param check_paths:
        :param correct_paths:
        :param confirm_correction:
        :param fix_success:
        :return:
        """

        # Initialize a list for the status of the different simulations
        status_list = []

        # Initialize flag for when assignment scheme is changed
        changed_assignment = False

        # Get the simulation names
        if self.has_assignment_table: nassignment_simulations = len(self.assignment_table)
        else: nassignment_simulations = None

        # Problem with the number of simulations
        if self.nsimulation_directories != self.nsimulations:
            if self.nsimulation_directories == nassignment_simulations:
                log.warning("Fixing the number of simulations in the generation info (" + str(nassignment_simulations) + " instead of " + str(self.nsimulations) + ") ...")
                self.info.nsimulations = nassignment_simulations
                self.info.save()
                all_simulation_names = self.simulation_directory_names
            else: all_simulation_names = self.simulation_names
        else: all_simulation_names = self.simulation_names

        # Show progress bar of getting status
        log.info("Checking status for each simulation ...")
        nsimulations = len(all_simulation_names)
        with Bar(label='', expected_size=nsimulations, every=1, add_datetime=True) as bar:

            # Loop over the simulations
            for index, simulation_name in enumerate(all_simulation_names):

                # Debugging
                #log.debug("Getting the status of simulation '" + simulation_name + "' ...")

                # Show progress
                bar.show(float(index + 1))

                # Get the status
                simulation_status, changed_assignment_sim = self.get_simulation_status(simulation_name, remotes=remotes,
                                                                                       lazy=lazy, find_simulations=find_simulations,
                                                                                       find_remotes=find_remotes, produce_missing=produce_missing,
                                                                                       retrieve=retrieve, check_paths=check_paths,
                                                                                       correct_paths=correct_paths, confirm_correction=confirm_correction,
                                                                                       fix_success=fix_success)
                if changed_assignment_sim: changed_assignment = True
                #print(simulation_name, simulation_status)

                # Add the status
                status_list.append(simulation_status)

        # Save the assignment table if it has been adapted
        if changed_assignment:
            log.debug("Saving the changed assignment table for the generation ...")
            self.assignment_table.save()

        # Create the simulation status table
        status = SimulationStatusTable.from_columns(self.simulation_names, status_list)

        # Return the status table
        return status

    # -----------------------------------------------------------------

    def get_simulation_status(self, simulation_name, remotes=None, lazy=False, find_simulations=False,
                              find_remotes=None, produce_missing=False, retrieve=False, check_paths=False,
                              correct_paths=False, confirm_correction=False, fix_success=True):

        """
        This function ...
        :param simulation_name:
        :param remotes:
        :param lazy:
        :param find_simulations:
        :param find_remotes:
        :param produce_missing:
        :param retrieve:
        :param check_paths:
        :param correct_paths:
        :param confirm_correction:
        :param fix_success:
        :return:
        """

        # Initialize flag for when assignment scheme is changed
        changed_assignment = False

        load_simulation = True

        # Don't load any simulation
        if lazy: load_simulation = False

        # Simulation not in assignment table
        elif self.has_assignment_table and not self.in_assignment(simulation_name):

            # Give warning
            log.warning("The simulation '" + simulation_name + "' is missing from the assignment table")

            # Find simulation?
            if find_simulations:

                changed_a, load_sim = self.find_simulation_for_assignment(simulation_name, find_remotes, remotes=remotes,
                                                                          produce_missing=produce_missing, correct_paths=correct_paths,
                                                                          confirm_correction=confirm_correction)
                changed_assignment |= changed_a
                load_simulation |= load_sim

            # Don't find: don't load the simulation
            else: load_simulation = False

        # Simulation not found
        elif not self.has_simulation(simulation_name):

            # Give warning
            log.warning("The simulation file for '" + simulation_name + "' is not present anymore")

            # Get host ID and simulation ID
            try: host_id = self.get_host_id(simulation_name)
            except MissingSimulation: host_id = None
            try: simulation_id = self.get_simulation_id(simulation_name)
            except MissingSimulation: simulation_id = None

            # Host ID or simulation ID is not found
            if host_id is None or simulation_id is None:

                if host_id is None: log.warning("Could not determine the host ID for simulation '" + simulation_name + "'")
                if simulation_id is None: log.warning("Could not determine the simulation ID for simulation '" + simulation_name + "'")
                load_simulation = False

            # Find simulation?
            elif find_simulations:

                # Give warning
                if simulation_id is not None: log.warning("Simulation ID was '" + str(simulation_id) + "'")
                if host_id is not None: log.warning("Simulation host ID was host ID is '" + host_id + "'")

                # Give warning
                log.warning("Looking for the simulation '" + simulation_name + "' on other remote hosts ...")
                simulation = find_simulation(simulation_name, find_remotes)

                # Simulation is not found
                if simulation is None:

                    # Produce missing simulation file?
                    if produce_missing:

                        # Warning
                        log.warning("Simulation was not found, producing simulation file ...")

                        # Produce
                        self.produce_missing_simulation_file(simulation_name, host_id, simulation_id=simulation_id,
                                                             remotes=remotes, correct_paths=correct_paths,
                                                             confirm_correction=confirm_correction)

                        # Set flag
                        load_simulation = True

                    # Don't produce: give error for missing simulation file
                    else: raise RuntimeError("Simulation '" + simulation_name + "' was not found")

                # Simulation is found
                else:

                    actual_host_id = simulation.host_id
                    actual_id = simulation.id
                    cluster_name = simulation.cluster_name
                    log.warning("Simulation identified as '" + str(actual_id) + "' for host '" + actual_host_id + "'")
                    log.warning("Fixing assignment table ...")
                    self.set_id_and_host_for_simulation(simulation_name, actual_id, actual_host_id, cluster_name=cluster_name)
                    changed_assignment = True
                    load_simulation = True

            # Don't load the simulation
            else: load_simulation = False

        # Load the simulation
        if load_simulation:

            # Load the simulation
            simulation = self.get_simulation(simulation_name)

            # Get properties
            host_id = simulation.host_id
            #simulation_id = simulation.id

            # Check paths
            if check_paths:
                if remotes is None or host_id not in remotes:
                    #log.warning("Cannot check paths if remote is not passed")
                    try: check_simulation_paths(simulation) # without checking remote paths
                    except RuntimeError as e:
                        if correct_paths:
                            log.warning(str(e))
                            log.warning("Fixing simulation paths ...")
                            correct_simulation_and_analysis_paths(simulation, confirm=confirm_correction)
                        else: raise e
                else:
                    try: check_simulation_paths(simulation, remote=remotes[host_id])
                    except RuntimeError as e:
                        if correct_paths:
                            log.warning(str(e))
                            log.warning("Fixing simulation paths ...")
                            correct_simulation_and_analysis_paths(simulation, confirm=confirm_correction)
                        else: raise e

        # Don't load the simulation
        else:

            # No simulation object
            simulation = None

            # Get host ID and simulation ID from generation assignment table
            try: host_id = self.get_host_id(simulation_name)
            except MissingSimulation: host_id = None
            #simulation_id = self.get_simulation_id(simulation_name)

        # No simulation object
        if simulation is None:

            has_misc = self.has_misc_output(simulation_name)
            has_plotting = self.has_plotting_output(simulation_name)
            has_extraction = self.has_extraction_output(simulation_name)
            has_timing = self.has_timing(simulation_name)
            has_memory = self.has_memory(simulation_name)

            # Has chi squared
            if self.is_analysed(simulation_name): simulation_status = "analysed"

            # Has any analysis output
            elif has_extraction or has_plotting or has_misc:

                analysed = []
                if has_extraction: analysed.append("extraction")
                if has_plotting: analysed.append("plotting")
                if has_misc: analysed.append("misc")
                if has_timing and has_memory: analysed.append("batch")
                simulation_status = "analysed: " + ", ".join(analysed)

            # Has simulation output
            elif self.is_retrieved(simulation_name): simulation_status = "retrieved"

            # No simulation output
            else: simulation_status = "unknown"

        # Simulation object
        else:

            # Already analysed
            if simulation.analysed: simulation_status = "analysed"

            # Partly analysed
            elif simulation.analysed_any:

                analysed = []
                if simulation.analysed_all_extraction: analysed.append("extraction")
                if simulation.analysed_all_plotting: analysed.append("plotting")
                if simulation.analysed_all_misc: analysed.append("misc")
                if simulation.analysed_batch: analysed.append("batch")
                if simulation.analysed_scaling: analysed.append("scaling")
                if simulation.analysed_all_extra: analysed.append("extra")

                if len(analysed) > 0: simulation_status = "analysed: " + ", ".join(analysed)
                else: simulation_status = "analysed: started"

            # Retrieved
            elif simulation.retrieved: simulation_status = "retrieved"

            # Finished
            elif simulation.finished: simulation_status = finished_name

            # Not retrieved
            elif remotes is not None:

                # Not yet retrieved, what is the status?
                screen_states_host = remotes.screens[host_id]
                jobs_status_host = remotes.jobs[host_id]
                if host_id in remotes: simulation_status = remotes[host_id].get_simulation_status(simulation, screen_states=screen_states_host, jobs_status=jobs_status_host)
                else: simulation_status = "unknown"

            # Unknown
            else: simulation_status = "unknown"

            #print(fix_success, self.assignment_table.is_launched(simulation.name), is_invalid_or_unknown_status(simulation_status))
            # Check success flag in assignment
            if fix_success and not self.assignment_table.is_launched(simulation.name) and not is_invalid_or_unknown_status(simulation_status):
                log.warning("Settting the launch of simulation '" + simulation.name + "' as succesful in the assignment table as this was not yet done")
                self.assignment_table.set_success_for_simulation(simulation.name)
                changed_assignment = True

            # Retrieve finished simulations?
            if simulation_status == finished_name and retrieve:
                if remotes is None or host_id not in remotes: log.warning("Cannot retrieve simulations if remotes are not passed")
                else:
                    remotes[host_id].retrieve_simulation(simulation)
                    simulation_status = "retrieved"

        # Return the status
        return simulation_status, changed_assignment

    # -----------------------------------------------------------------

    def produce_missing_simulation_file(self, simulation_name, host_id, simulation_id=None, remotes=None,
                                        correct_paths=False, confirm_correction=False):

        """
        This function ...
        :param simulation_name:
        :param host_id:
        :param simulation_id:
        :param remotes:
        :param correct_paths:
        :param confirm_correction:
        :return:
        """

        # Check
        if host_id is None: raise ValueError("Cannot produce missing simulation when host ID is unknown")

        # Find first other simulation for this host that exists
        other_simulation = None
        for other_simulation_name in self.get_simulation_names_for_host(host_id):

            if not self.has_simulation(other_simulation_name): continue
            other_simulation = self.get_simulation(other_simulation_name)
            break

        if other_simulation is None: raise RuntimeError("Cannot produce missing simulation file: no other existing simulation file for host '" + host_id + "'")

        # Produce the simulation from the other simulation
        cluster_name = self.get_cluster_name(simulation_name)
        other_remote_input_path = other_simulation.remote_input_path
        simulation = produce_simulation_from_other(other_simulation, simulation_name, host_id, simulation_id, cluster_name)

        # Check simulation paths
        if remotes is None or host_id not in remotes: log.warning("Cannot check the simulation paths: remote is not passed")
        else:
            try:
                check_simulation_paths(simulation, remote=remotes[host_id],
                                       other_remote_input_path=other_remote_input_path)
            except RuntimeError as e:
                if correct_paths:
                    log.warning(str(e))
                    log.warning("Fixing simulation paths ...")
                    correct_simulation_and_analysis_paths(simulation, confirm=confirm_correction)
                else: raise e

        # Save the simulation
        simulation.save()

# -----------------------------------------------------------------

# Initialize dictionary
simulations_hosts = dict()

# -----------------------------------------------------------------

def get_simulations(host_id):

    """
    This function gets simulations for a remote host
    :param host_id:
    :return:
    """

    # Get the simulations
    simulations_host = get_simulations_for_host(host_id, as_dict=True)

    # Set the simulations
    simulations_hosts[host_id] = simulations_host

    # Return the simulations
    return simulations_host

# -----------------------------------------------------------------

def find_simulation(simulation_name, host_ids):

    """
    This function ...
    :param simulation_name:
    :param host_ids:
    :return:
    """

    # Loop over the remotes
    for host_id in host_ids:

        simulation = find_simulation_for_host(simulation_name, host_id)
        if simulation is not None: return simulation

    # Return the simulation
    return None

# -----------------------------------------------------------------

def find_simulation_for_host(simulation_name, host_id):

    """
    This function ...
    :param simulation_name:
    :param host_id:
    :return:
    """

    if host_id not in simulations_hosts: simulations_host = get_simulations(host_id)
    else: simulations_host = simulations_hosts[host_id]

    # Check whether the simulation name is in the
    if simulation_name in simulations_host:
        simulation = simulations_host[simulation_name]
        assert simulation.host_id == host_id  # check that the simulation object has the correct host ID as attribute
        return simulation

    else: return None

# -----------------------------------------------------------------

def produce_simulation_from_other(other_simulation, simulation_name, host_id=None, simulation_id=None, cluster_name=None):

    """
    This function ...
    :param other_simulation:
    :param simulation_name:
    :param host_id:
    :param simulation_id:
    :param cluster_name:
    :return:
    """

    # Debugging
    log.debug("Producing simulation '" + simulation_name + "' from simulation '" + other_simulation.name + "' ...")

    # Determine simulations directories
    simulation_path = other_simulation.path
    simulation_dirpath = fs.directory_of(simulation_path)

    # Determine other simulation name
    other_simulation_name = other_simulation.name

    # Copy the simulation
    simulation = deepcopy(other_simulation)

    # Set the new path for the simulation
    simulation.path = fs.join(simulation_dirpath, str(simulation_id) + ".sim")

    # Set the simulation name, ID, host ID and cluster name
    simulation.name = simulation_name
    simulation.id = simulation_id # can be None, shouldn't be the same as original simulation
    if host_id is not None: simulation.host_id = host_id
    if cluster_name is not None: simulation.cluster_name = cluster_name

    # Create config for adapting
    adapt_config = dict()
    adapt_config["contains"] = "path"
    adapt_config["replace_string"] = (other_simulation_name, simulation_name)
    adapt_config["types"] = ["string"]
    adapt_config["only_replacements"] = True

    # Adapt all paths in the simulation file
    adapt_simulation(simulation, adapt_config)
    adapt_analysis(simulation, adapt_config)

    # Return the simulation
    return simulation

# -----------------------------------------------------------------

def check_simulation_paths(simulation, remote=None, other_remote_input_path=None):

    """
    This function ...
    :param simulation:
    :param remote:
    :param other_remote_input_path:
    :return:
    """

    # Get simulation name
    simulation_name = simulation.name

    # Check local paths
    if simulation.has_input_directory and (not fs.is_directory(simulation.input_path) or fs.is_empty(simulation.input_path)): raise IOError("Local input directory '" + simulation.input_path + "' is missing or empty")
    if not fs.is_directory(simulation.output_path): raise IOError("Local output directory '" + simulation.output_path + "' is missing")
    if not fs.is_directory(simulation.base_path): raise IOError("Local simulation directory '" + simulation.base_path + "' is missing")
    if not fs.is_file(simulation.ski_path): raise IOError("Local ski file '" + simulation.ski_path + "' is missing")
    if simulation_name not in simulation.output_path: raise RuntimeError("Something went wrong: wrong output path for simulation '" + simulation_name + "'")
    if simulation_name not in simulation.base_path: raise RuntimeError("Something went wrong: wrong base path for simulation '" + simulation_name + "'")
    if simulation_name not in simulation.ski_path: raise RuntimeError("Something went wrong: wrong skifile path for simulation '" + simulation_name + "'")

    # Check analysis paths
    if simulation.analysis.extraction.path is not None:
        path = simulation.analysis.extraction.path
        if not fs.is_directory(path): raise IOError("Extraction directory '" + path + "' is missing")
        if simulation_name not in path: raise RuntimeError("Something went wrong: wrong extraction path for simulation '" + simulation_name + "'")
    if simulation.analysis.plotting.path is not None:
        path = simulation.analysis.plotting.path
        if not fs.is_directory(path): raise IOError("Plotting directory '" + path + "' is missing")
        if simulation_name not in path: raise RuntimeError("Something went wrong: wrong plotting path for simulation '" + simulation_name + "'")
    if simulation.analysis.misc.path is not None:
        path = simulation.analysis.misc.path
        if not fs.is_directory(path): raise IOError("Misc directory '" + path + "' is missing")
        if simulation_name not in path: raise RuntimeError("Something went wrong: wrong misc path for simulation '" + simulation_name + "'")

    # Check remote paths
    if remote is not None:
        if not remote.is_file(simulation.remote_ski_path): raise IOError("Remote ski file '" + simulation.remote_ski_path + "' is missing")
        if not remote.is_directory(simulation.remote_simulation_path): raise IOError("Remote simulation path '" + simulation.remote_simulation_path + "' is missing")
        if not remote.is_directory(simulation.remote_input_path) or remote.is_empty(simulation.remote_input_path):
            if other_remote_input_path is not None and remote.is_directory(other_remote_input_path) and not remote.is_empty(other_remote_input_path): simulation.remote_input_path = other_remote_input_path
            else: raise IOError("Remote input path for simulation '" + simulation_name + "' is not found")
        if not remote.is_directory(simulation.remote_output_path): raise IOError("Remote output path '" + simulation.remote_output_path + "' is missing")
        if simulation_name not in simulation.remote_ski_path: raise RuntimeError("Something went wrong: wrong remote skifile path for simulation '" + simulation_name + "'")
        if simulation_name not in simulation.remote_simulation_path: raise RuntimeError("Something went wrong: wrong remote base path for simulation '" + simulation_name + "'")
        if simulation_name not in simulation.remote_output_path: raise RuntimeError("Something went wrong: wrong remote output path for simulation '" + simulation_name + "'")

# -----------------------------------------------------------------

def correct_simulation_and_analysis_paths(simulation, confirm=False):

    """
    This function ...
    :param simulation:
    :param confirm:
    :return:
    """

    # Correct simulation paths
    correct_simulation_paths(simulation, confirm=confirm)

    # Correct analysis paths
    correct_analysis_paths(simulation.analysis, simulation.name, confirm=confirm)

    # Save the simulation
    simulation.save()

# -----------------------------------------------------------------

def correct_simulation_paths(simulation, confirm=False):

    """
    This function ...
    :param simulation:
    :param confirm:
    :return:
    """

    # Get simulation name
    simulation_name = simulation.name

    # Check output path
    if simulation_name not in simulation.output_path:

        log.debug("Fixing local output path ...")
        other_simulation_name = fs.name(fs.directory_of(simulation.output_path))

        log.debug("Path was pointing to output path of simulation '" + other_simulation_name + "'")
        correct_path = simulation.output_path.replace(other_simulation_name, simulation_name)

        if not confirm or prompt_yn("replace", "Replace '" + simulation.output_path + "' with '" + correct_path + "'?", default=True):
            log.debug("Replacing with new path: '" + correct_path + "' ...")
            simulation.output_path = correct_path

    # Check ski path
    if simulation_name not in simulation.ski_path:

        log.debug("Fixing local skifile path ...")
        other_simulation_name = fs.name(fs.directory_of(simulation.ski_path))

        log.debug("Path was pointing to skifile path of simulation '" + other_simulation_name + "'")
        correct_path = simulation.ski_path.replace(other_simulation_name, simulation_name)

        if not confirm or prompt_yn("replace", "Replace '" + simulation.ski_path + "' with '" + correct_path + "'?", default=True):
            log.debug("Replacing with new path: '" + correct_path + "' ...")
            simulation.ski_path = correct_path

    # Check base path
    if simulation_name not in simulation.base_path:

        log.debug("Fixing local base path ...")
        other_simulation_name = fs.name(simulation.base_path)

        log.debug("Path was pointing to base path of simulation '" + other_simulation_name + "'")
        #correct_path = simulation.base_path.replace(other_simulation_name, simulation_name)
        #log.debug("Replacing with new path: '" + correct_path + "' ...")
        #simulation.base_path = correct_path
        raise RuntimeError("Simulation base path cannot be fixed")

    # Check remote ski path
    if simulation_name not in simulation.remote_ski_path:

        log.debug("Fixing remote skifile path ...")
        other_simulation_name = fs.name(fs.directory_of(simulation.remote_ski_path))

        log.debug("Path was pointing to remote skifile path of simulation '" + other_simulation_name + "'")
        correct_path = simulation.remote_ski_path.replace(other_simulation_name, simulation_name)

        if not confirm or prompt_yn("replace", "Replace '" + simulation.remote_ski_path + "' with '" + correct_path + "'?", default=True):
            log.debug("Replacing with new path: '" + correct_path + "' ...")
            simulation.remote_ski_path = correct_path

    # Check remote base path
    if simulation_name not in simulation.remote_simulation_path:

        log.debug("Fixing remote base path ...")
        other_simulation_name = fs.name(simulation.remote_simulation_path)

        log.debug("Path was pointing to remote base path of simulation '" + other_simulation_name + "'")
        correct_path = simulation.remote_simulation_path.replace(other_simulation_name, simulation_name)

        if not confirm or prompt_yn("replace", "Replace '" + simulation.remote_simulation_path + "' with '" + correct_path + "'?", default=True):
            log.debug("Replacing with new path: '" + correct_path + "' ...")
            simulation.remote_simulation_path = correct_path

    # Check remote output path
    if simulation_name not in simulation.remote_output_path:

        log.debug("Fixing remote output path ...")
        other_simulation_name = fs.name(fs.directory_of(simulation.remote_output_path))

        log.debug("Path was pointing to remote output path of simulation '" + other_simulation_name + "'")
        correct_path = simulation.remote_output_path.replace(other_simulation_name, simulation_name)

        if not confirm or prompt_yn("replace", "Replace '" + simulation.remote_output_path + "' with '" + correct_path + "'?", default=True):
            log.debug("Replacing with new path: '" + correct_path + "' ...")
            simulation.remote_output_path = correct_path

# -----------------------------------------------------------------

def correct_analysis_paths(analysis, simulation_name, confirm=False):

    """
    This function ...
    :param analysis:
    :param simulation_name:
    :param confirm:
    :return:
    """

    # Check extraction path
    if analysis.extraction.path is not None:

        log.debug("Fixing extraction path ...")
        other_simulation_name = fs.name(fs.directory_of(analysis.extraction.path))

        log.debug("Path was pointing to extraction path of simulation '" + other_simulation_name + "'")
        correct_path = analysis.extraction.path.replace(other_simulation_name, simulation_name)

        if not confirm or prompt_yn("replace", "Replace '" + analysis.extraction.path + "' with '" + correct_path + "'?", default=True):
            log.debug("Replacing with new path: '" + correct_path + "' ...")
            analysis.extraction.path = correct_path

    # Check plotting path
    if analysis.plotting.path is not None:

        log.debug("Fixing plotting path ...")
        other_simulation_name = fs.name(fs.directory_of(analysis.plotting.path))

        log.debug("Path was pointing to plotting path of simulation '" + other_simulation_name + "'")
        correct_path = analysis.plotting.path.replace(other_simulation_name, simulation_name)

        if not confirm or prompt_yn("replace", "Replace '" + analysis.plotting.path + "' with '" + correct_path + "'?", default=True):
            log.debug("Replacing with new path: '" + correct_path + "' ...")
            analysis.plotting.path = correct_path

    # Check misc path
    if analysis.misc.path is not None:

        log.debug("Fixing miscellaneous path ...")
        other_simulation_name = fs.name(fs.directory_of(analysis.misc.path))

        log.debug("Path was pointing to misc path of simulation '" + other_simulation_name + "'")
        correct_path = analysis.misc.path.replace(other_simulation_name, simulation_name)

        if not confirm or prompt_yn("replace", "Replace '" + analysis.misc.path + "' with '" + correct_path + "'?", default=True):
            log.debug("Replacing with new path: '" + correct_path + "' ...")
            analysis.misc.path = correct_path

# -----------------------------------------------------------------
