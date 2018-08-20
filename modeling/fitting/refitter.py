#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.refitter Contains the Refitter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import OrderedDict, defaultdict

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty, memoize_method
from .tables import WeightsTable
from ...core.tools import filesystem as fs
from ...core.tools import tables
from .modelanalyser import FluxDifferencesTable
from ...core.tools import sequences
from ...core.misc.fluxes import ObservedFluxCalculator
from .tables import ChiSquaredTable, BestParametersTable, ModelProbabilitiesTable, ParameterProbabilitiesTable
from ...core.basics.distribution import Distribution
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr
from ...core.tools import nr, numbers
from ...core.plot.distribution import plot_distributions, plot_distribution
from .run import FittingRun
from .tables import GenerationsTable
from ...core.filter.filter import parse_filter
from .generation import GenerationInfo
from .weights import WeightsCalculator
from .backup import FitBackupper
from .modelanalyser import calculate_chi_squared_from_differences

# -----------------------------------------------------------------

earth_instrument_name = "earth"

# -----------------------------------------------------------------

class Refitter(FittingComponent):
    
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
        super(Refitter, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The fitting run
        self.fitting_run = None

        # The table of weights for each band
        self.weights = None

        # The fluxes
        self.fluxes = defaultdict(defaultdict)

        # The flux differences
        self.differences = defaultdict(defaultdict)

        # The chi squared tables
        self.chi_squared_tables = dict()

        # The best parameters table
        self.best_parameters_table = None

        # Model probability tables
        self.model_probabilities = dict()

        # Parameter probability tables
        self.parameter_probabilities = dict()

        # Parameter probabilities for all generations
        self.parameter_probabilities_all = dict()

        # Probability distributions
        self.distributions = dict()

        # Best simulations
        self.best_simulations = dict()
        self.best_simulation_counts = dict()

        # New run, directory path and new fitting configuration
        self.new_run = None
        self.new_fitting_config = None

        # New paths
        self.new_prob_generation_paths = dict()
        self.new_generation_paths = dict()
        self.new_simulation_paths = defaultdict(dict)
        self.new_simulation_misc_paths = defaultdict(dict)

    # -----------------------------------------------------------------

    @property
    def do_create_run(self):

        """
        This function ...
        :return:
        """

        return self.as_run

    # -----------------------------------------------------------------

    @property
    def do_backup_run(self):

        """
        This function ...
        :return:
        """

        return (self.in_place and not self.individual_simulations) and fs.is_empty(self.backup_path)

    # -----------------------------------------------------------------

    @property
    def do_show(self):

        """
        This function ...
        :return:
        """

        return self.config.show

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

        # Create run directory
        if self.do_create_run: self.create_run()

        # Make backups where necessary
        self.backup()

        # Get the weights
        self.get_weights()

        # Get the fluxes
        self.get_fluxes()

        # Get the differences
        self.get_differences()

        # Calculate the chi squared values
        self.calculate_chi_squared()

        # Get the parameters of the best models for each generation
        self.get_best_parameters()

        # Calculate the probabilities
        self.calculate_probabilities()

        # Calculate the probability distributions
        self.create_distributions()

        # Get best simulation
        self.get_best_simulations()

        # Create the configuration
        self.create_config()

        # Writing
        self.write()

        # Show
        if self.do_show: self.show()

        # Plot
        if self.do_plot: self.plot()

    # -----------------------------------------------------------------

    @property
    def refitting_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.refitting_path

    # -----------------------------------------------------------------

    @lazyproperty
    def path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.refitting_path, self.config.name)

    # -----------------------------------------------------------------

    @lazyproperty
    def prob_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.path, "prob")

    # -----------------------------------------------------------------

    @lazyproperty
    def prob_generations_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.prob_path, "generations")

    # -----------------------------------------------------------------

    @lazyproperty
    def prob_generations_paths(self):

        """
        This function ...
        :return:
        """

        paths = dict()
        for generation_name in self.generation_names: paths[generation_name] = fs.create_directory_in(self.prob_generations_path, generation_name)
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def prob_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.prob_path, "parameters")

    # -----------------------------------------------------------------

    @lazyproperty
    def prob_distributions_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.prob_path, "distributions")

    # -----------------------------------------------------------------

    @lazyproperty
    def generations_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.path, "generations")

    # -----------------------------------------------------------------

    @lazyproperty
    def generation_paths(self):

        """
        This function ...
        :return:
        """

        paths = dict()
        for generation_name in self.generation_names: paths[generation_name] = fs.create_directory_in(self.generations_path, generation_name)
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def chi_squared_table_paths(self):

        """
        This function ...
        :return:
        """

        paths = dict()
        for generation_name in self.generation_names: paths[generation_name] = fs.join(self.generation_paths[generation_name], "chi_squared.dat")
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_paths(self):

        """
        This function ...
        :return:
        """

        paths = defaultdict(dict)
        for generation_name in self.generation_names:
            generation = self.generations[generation_name]
            for simulation_name in generation.simulation_names:
                paths[generation_name][simulation_name] = fs.create_directory_in(self.generation_paths[generation_name], simulation_name)
        return paths

    # -----------------------------------------------------------------

    @property
    def as_run(self):

        """
        This function ...
        :return:
        """

        return self.config.as_run is not None

    # -----------------------------------------------------------------

    @property
    def new_run_name(self):

        """
        This function ...
        :return:
        """

        return self.config.as_run

    # -----------------------------------------------------------------

    @property
    def in_place(self):

        """
        This function ...
        :return:
        """

        return self.config.in_place is not None

    # -----------------------------------------------------------------

    @property
    def backup_name(self):

        """
        This function ...
        :return:
        """

        return self.config.in_place

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Refitter, self).setup(**kwargs)

        # Load the fitting run
        if kwargs.get("fitting_run", None) is not None: self.fitting_run = kwargs.pop("fitting_run")
        else: self.fitting_run = self.load_fitting_run(self.config.run)

        # Check options
        if self.as_run and self.in_place: raise ValueError("Cannot refit as new fitting run and refit in-place simultaneously")
        if self.config.name is None and not (self.as_run or self.in_place): raise ValueError("Refitting name must be specified when not refitting as new run or in-place")
        if self.in_place:
            if self.config.generations is not None and self.config.simulations is None: raise ValueError("Cannot specify generations without specifying individual simulations when refitting is done in place")
        if self.as_run and self.individual_simulations: raise ValueError("Cannot refit as new fitting run and specify individual simulation names")
        if self.individual_simulations:
            if self.reweigh: raise ValueError("Cannot reweigh when refitting only individual simulations")

        # Create the table to contain the weights
        self.weights = WeightsTable()

        # Load chi squared tables
        self.load_chi_squared()

        # Initialize best parameters table
        self.best_parameters_table = BestParametersTable(parameters=self.free_parameter_labels, units=self.parameter_units)
        self.best_parameters_table._setup()

    # -----------------------------------------------------------------

    def load_chi_squared(self):

        """
        This fucntion ...
        :return:
        """

        # Initialize chi squared table for each generation
        for generation_name in self.generation_names:

            # Only some simulations of the single generation are going to be refitted
            if self.individual_simulations:

                # No chi squared table found?
                if not self.single_generation.has_chi_squared_table:

                    # Determine backup filepath, look for chi squared backup table
                    backup_filepath = fs.get_backup_filepath(self.single_generation.chi_squared_table_path, prefix=self.backup_name, sep="__")
                    if fs.is_file(backup_filepath):

                        # Load
                        table = ChiSquaredTable.from_file(backup_filepath)
                        table.path = None

                    # Not found?
                    else:

                        # Warn
                        log.warning("Chi squared table for the '" + self.single_generation_name + "' generation could not be found: creating new one, filling in the chi squared based on the flux differences ...")

                        # Create chi squared table
                        table = ChiSquaredTable()  # new one anyway
                        for simulation_name in self.single_generation.simulation_names:
                            if simulation_name in self.simulation_names: continue # don't get chi squared, is going to be refitted anyway
                            if self.single_generation.has_sed_differences(simulation_name): differences = self.single_generation.get_simulation_sed_differences(simulation_name)
                            else:
                                #backup_differences_filepath = fs.get_backup_filepath(self.single_generation.get_simulation_sed_differences_path(simulation_name), prefix=self.backup_name, sep="__")
                                #if not fs.is_file(backup_differences_filepath): raise IOError("Could not find existing differences file for simulation '" + simulation_name + "'")
                                #differences = FluxDifferencesTable.from_file(backup_differences_filepath)
                                log.warning("No flux differences can be found for simulation '" + simulation_name + "': skipping ...")
                                continue
                            chi_squared = calculate_chi_squared_from_differences(differences, self.nfree_parameters)
                            table.add_entry(simulation_name, chi_squared)

                # Load the existing chi squared table
                else: table = self.single_generation.chi_squared_table

            # All simulations of the generation are goging to be refitted
            else: table = ChiSquaredTable()  # new table since all simulations are going to be re-evaluated anyway

            # Set the table for the generation
            self.chi_squared_tables[generation_name] = table

    # -----------------------------------------------------------------

    @lazyproperty
    def backup_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.refitting_path, self.backup_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def backup_generations_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.backup_path, "generations")

    # -----------------------------------------------------------------

    @memoize_method
    def backup_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.create_directory_in(self.backup_generations_path, generation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def backup_path_for_simulation(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return fs.create_directory_in(self.backup_path_for_generation(generation_name), simulation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_prob(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.fitting_run.prob_path, recursive=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_best(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.fitting_run.best_path, recursive=True)

    # -----------------------------------------------------------------

    def backup(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making backups before refitting ...")

        # Backup individual simulations results
        if self.individual_simulations:

            # Backup the chi squared table
            self.backup_chi_squared()

            # Backup simulations
            self.backup_simulations()

        # Make a backup of the current fitting run
        if self.do_backup_run: self.backup_run()

    # -----------------------------------------------------------------

    def backup_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating backup of chi squared table for generation '" + self.single_generation_name + "' ...")

        # Copy the file
        fs.backup_file(self.single_generation.chi_squared_table_path, remove=True, prefix=self.backup_name, sep="__", exists="pass", check_filepath=False, remove_if_exists=True)

    # -----------------------------------------------------------------

    def backup_simulations(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.info("Making a backup of simulations that will be refitted ...")

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check whether the simulation has differences
            if not self.single_generation.has_sed_differences(simulation_name):
                log.warning("No differences table for simulation '" + simulation_name + "' of generation '" + self.single_generation_name + "'")
                continue

            # Debugging
            log.debug("Creating backup of differences table for simulation '" + simulation_name + "' ...")

            # Copy the file
            fs.backup_file(self.single_generation.get_simulation_sed_differences_path(simulation_name), remove=True, prefix=self.backup_name, sep="__", exists="pass", check_filepath=False, remove_if_exists=True)

    # -----------------------------------------------------------------

    def backup_run(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.info("Making a backup of the current fitting results ...")

        # Create the backupper
        backupper = FitBackupper()

        # Set the name for the backup
        backupper.config.name = self.backup_name

        # Set the generation names
        backupper.config.generations = self.generation_names

        # Set modeling path
        backupper.config.path = self.config.path

        # Run the backupper
        backupper.run(fitting_run=self.fitting_run)

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.model_for_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def new_run_path(self):

        """
        This function ...
        :return:
        """

        return self.new_run.path

    # -----------------------------------------------------------------

    @property
    def new_best_path(self):

        """
        This function ...
        :return:
        """

        return self.new_run.best_path

    # -----------------------------------------------------------------

    @property
    def new_generations_path(self):

        """
        This function ...
        :return:
        """

        return self.new_run.generations_path

    # -----------------------------------------------------------------

    @property
    def new_prob_path(self):

        """
        This function ...
        :return:
        """

        return self.new_run.prob_path

    # -----------------------------------------------------------------

    @property
    def new_prob_generations_path(self):

        """
        Thisfunction ...
        :return:
        """

        return self.new_run.prob_generations_path

    # -----------------------------------------------------------------

    @property
    def new_prob_parameters_path(self):

        """
        This function ...
        :return:
        """

        return self.new_run.prob_parameters_path

    # -----------------------------------------------------------------

    @property
    def new_prob_distributions_path(self):

        """
        This function ...
        :return:
        """

        return self.new_run.prob_distributions_path

    # -----------------------------------------------------------------

    @property
    def new_geometries_path(self):

        """
        This function ...
        :return:
        """

        return self.new_run.geometries_path

    # -----------------------------------------------------------------

    @property
    def new_wavelength_grids_path(self):

        """
        This function ...
        :return:
        """

        return self.new_run.wavelength_grids_path

    # -----------------------------------------------------------------

    @property
    def new_generations_table_path(self):

        """
        This function ...
        :return:
        """

        return self.new_run.generations_table_path

    # -----------------------------------------------------------------

    def create_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the new fitting run ...")

        # Clone the fitting run
        self.new_run = clone_fitting_run(self.fitting_run, self.new_run_name, generations=self.generation_names,
                          new_prob_generation_paths=self.new_prob_generation_paths,
                          new_generation_paths=self.new_generation_paths, new_simulation_paths=self.new_simulation_paths,
                          new_simulation_misc_paths=self.new_simulation_misc_paths, clone_configuration=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def new_generations_table(self):

        """
        This function ...
        :return:
        """

        return GenerationsTable.from_file(self.new_generations_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def generation_names(self):

        """
        This function ...
        :return:
        """

        if self.config.generations is not None: return self.config.generations
        else: return self.fitting_run.generation_names

    # -----------------------------------------------------------------

    @property
    def ngenerations(self):
        return len(self.generation_names)

    # -----------------------------------------------------------------

    @property
    def has_single_generation(self):
        return self.ngenerations == 1

    # -----------------------------------------------------------------

    @property
    def single_generation_name(self):
        if not self.has_single_generation: raise ValueError("Not a single generation")
        return self.generation_names[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def generations(self):

        """
        This function ...
        :return:
        """

        gens = OrderedDict()
        for name in self.generation_names: gens[name] = self.fitting_run.get_generation(name)
        return gens

    # -----------------------------------------------------------------

    @property
    def single_generation(self):
        return self.generations[self.single_generation_name]

    # -----------------------------------------------------------------

    @property
    def individual_simulations(self):

        """
        This function ...
        :return:
        """

        # Simulation names are specified
        if self.config.simulations is not None:

            # Check
            if not self.has_single_generation: raise ValueError("Cannot specify individual simulation names when there are multiple generations specified")
            return True

        # All simulations of generation(s)
        else: return False

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        if self.config.simulations is not None: return self.config.simulations
        else: return self.single_generation.simulation_names

    # -----------------------------------------------------------------

    def simulation_names_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        if self.individual_simulations: return self.simulation_names
        else: return self.generations[generation_name].simulation_names

    # -----------------------------------------------------------------

    @lazyproperty
    def not_filters(self):

        """
        This function ...
        :return:
        """

        return map(parse_filter, self.config.not_filters)

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):

        """
        This function ...
        :return:
        """

        if self.config.filters is not None: filters = self.config.filters
        else: filters = self.fitting_run.fitting_filters

        # Remove certain filters?
        if self.config.not_filters is not None: return sequences.removed(filters, self.not_filters)
        else: return filters

    # -----------------------------------------------------------------

    @lazyproperty
    def filter_names(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def different_filters(self):

        """
        This function ...
        :return:
        """

        #if self.config.filters is not None: return sequences.same_contents(self.config.filters, self.fitting_run.fitting_filters)
        #else: return False
        return not sequences.same_contents(self.filters, self.fitting_run.fitting_filters)

    # -----------------------------------------------------------------

    @lazyproperty
    def more_filters(self):

        """
        This function ...
        :return:
        """

        #if self.config.filters is None: return False
        #else: return len(self.config.filters) > len(self.fitting_run.fitting_filters)
        return len(self.filters) > len(self.fitting_run.fitting_filters)

    # -----------------------------------------------------------------

    @lazyproperty
    def additional_filters(self):

        """
        This function ...
        :return:
        """

        #if self.config.filters is None: return False
        #else: return sequences.has_other(self.config.filters, self.fitting_run.fitting_filters)
        answer = sequences.has_other(self.filters, self.fitting_run.fitting_filters)

        # Give warning
        if answer: log.warning("Fluxes are requested in filters additional to the original fitting filters: this will require the mock flux calculation to be performed again")

        # Return
        return answer

    # -----------------------------------------------------------------

    @property
    def reweigh(self):

        """
        This function ...
        :return:
        """

        if self.config.reweigh is not None: return self.config.reweigh
        else: return self.different_filters

    # -----------------------------------------------------------------

    def get_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting weights ...")

        # Recalculate the weights
        if self.reweigh: self.calculate_weights()

        # Load the original weights
        else: self.load_weights()

    # -----------------------------------------------------------------

    def load_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the weights to give to each band ...")

        # Load the weights
        self.weights = self.fitting_run.weights

    # -----------------------------------------------------------------

    def calculate_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the weight to give to each band ...")

        # Calculate
        calculator = WeightsCalculator(self.config.weighing)
        calculator.config.write = False
        calculator.config.show = True
        calculator.run(filters=self.filters)

        # Set weights
        self.weights = calculator.table

    # -----------------------------------------------------------------

    @property
    def reflux(self):

        """
        This function ...
        :return:
        """

        if self.config.reflux is not None: return self.config.reflux
        else: return self.different_spectral_convolution_any_generation or self.different_fluxes_from_images_any_generation or self.additional_filters

    # -----------------------------------------------------------------

    def get_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting mock observed fluxes ...")

        # Recalculate the fluxes
        if self.reflux: self.calculate_fluxes()

        # Load the original fluxes
        else: self.load_fluxes()

    # -----------------------------------------------------------------

    def load_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading mock observed fluxes ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Loading the fluxes for the '" + generation_name + "' generation ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Loop over the simulation names
            for simulation_name in generation.simulation_names:

                # Check whether mock SED is present
                if not generation.has_mock_sed(simulation_name):
                    log.warning("No mock fluxes for '" + simulation_name + "': skipping ...")
                    continue

                # Debugging
                log.debug("Loading mock fluxes for the '" + simulation_name + "' simulation ...")

                # Load the fluxes
                fluxes = generation.get_mock_sed(simulation_name)

                # Set fluxes
                self.fluxes[generation_name][simulation_name] = fluxes

                # Copy the fluxes directory into the new generation directory
                if self.as_run:
                    fluxes_path = fs.join(generation.get_simulation_misc_fluxes_path(simulation_name))
                    fs.copy_directory(fluxes_path, self.new_simulation_misc_paths[generation_name][simulation_name])

    # -----------------------------------------------------------------

    def spectral_convolution_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        if self.config.spectral_convolution is not None: return self.config.spectral_convolution
        else: return self.generations[generation_name].spectral_convolution

    # -----------------------------------------------------------------

    def different_spectral_convolution_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        if self.config.spectral_convolution is None: return False
        else: return self.config.spectral_convolution != self.generations[generation_name].spectral_convolution

    # -----------------------------------------------------------------

    @lazyproperty
    def different_spectral_convolution_any_generation(self):

        """
        This function ...
        :return:
        """

        for generation_name in self.generation_names:
            if self.different_spectral_convolution_for_generation(generation_name): return True
        return False

    # -----------------------------------------------------------------

    def fluxes_from_images_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        if self.config.fluxes_from_images is not None: return self.config.fluxes_from_images
        else: return self.generations[generation_name].use_images

    # -----------------------------------------------------------------

    def different_fluxes_from_images_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        if self.config.fluxes_from_images is None: return False
        else: return self.fluxes_from_images_for_generation(generation_name) != self.generations[generation_name].use_images

    # -----------------------------------------------------------------

    @lazyproperty
    def different_fluxes_from_images_any_generation(self):

        """
        This function ...
        :return:
        """

        for generation_name in self.generation_names:
            if self.different_fluxes_from_images_for_generation(generation_name): return True
        return False

    # -----------------------------------------------------------------

    def calculate_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating mock observed fluxes ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Calculating the mock observed fluxes for the '" + generation_name + "' generation ...")

            # Get the generation
            generation = self.generations[generation_name]
            from_images = self.fluxes_from_images_for_generation(generation_name)
            spectral_convolution = self.spectral_convolution_for_generation(generation_name)

            # Loop over the simulation names
            for simulation_name in self.simulation_names_for_generation(generation_name):

                # Get the simulation
                simulation = generation.get_simulation_or_basic(simulation_name)

                # Debugging
                log.debug("Calculating the mock observed fluxes for the '" + simulation_name + "' simulation ...")

                # Determine fluxes calculation output path
                if self.as_run:
                    misc_path = self.new_simulation_misc_paths[generation_name][simulation_name]
                    fluxes_path = fs.create_directory_in(misc_path, "fluxes")
                    output_path = fluxes_path
                else: output_path = self.simulation_paths[generation_name][simulation_name]

                # Calculate fluxes
                if from_images: fluxes = self.calculate_observed_fluxes_from_images_for_simulation(simulation, output_path, spectral_convolution)
                else: fluxes = self.calculate_observed_fluxes_for_simulation(simulation, output_path, spectral_convolution)

                # Set fluxes
                self.fluxes[generation_name][simulation_name] = fluxes

    # -----------------------------------------------------------------

    def calculate_observed_fluxes_for_simulation(self, simulation, output_path, spectral_convolution):

        """
        This function ...
        :param simulation:
        :param output_path:
        :param spectral_convolution:
        :return:
        """

        # Create the flux calculator
        flux_calculator = ObservedFluxCalculator()

        # Set spectral convolution flag
        flux_calculator.config.spectral_convolution = spectral_convolution

        # Set plot flag
        flux_calculator.config.plot = True

        # Run
        flux_calculator.run(simulation=simulation, output_path=output_path, filter_names=self.filter_names, instrument_names=[earth_instrument_name], reference_sed=self.reference_sed)

        # Get the mock observed SED for the earth instrument
        return flux_calculator.mock_seds[earth_instrument_name]

    # -----------------------------------------------------------------

    def calculate_observed_fluxes_from_images_for_simulation(self, simulation, output_path, spectral_convolution):

        """
        This function ...
        :param simulation:
        :param output_path:
        :param spectral_convolution:
        :return:
        """

        # Create
        # image_flux_calculator = ObservedFluxCalculator()
        #
        # # Set spectral convolution flag
        # image_flux_calculator.config.spectral_convolution = spectral_convolution
        #
        # # Set plot flag
        # image_flux_calculator.config.plot = True
        #
        # # Set from images flag
        # image_flux_calculator.config.from_images = True
        # image_flux_calculator.config.write_images = True
        #
        # # Run
        # image_flux_calculator.run(simulation=simulation,
        #                                filter_names=self.filter_names,
        #                                instrument_names=[earth_instrument_name],
        #                                coordinate_systems=self.fluxes_from_images_coordinate_systems,
        #                                masks=self.fluxes_from_images_masks, reference_sed=self.fluxes_from_images_reference_sed)
        #
        # # Get the mock observed SED for the earth instrument
        # return image_flux_calculator.mock_seds[earth_instrument_name]

        raise NotImplementedError("Not yet implemented")

    # -----------------------------------------------------------------

    @property
    def rediff(self):

        """
        Thisf unction ...
        :return:
        """

        if self.config.rediff is not None:
            if not self.config.rediff and self.config.additional_error is not None: raise ValueError("The differences have to be calculated again if an additional relative error has to be used")
            return self.config.rediff
        else: return self.different_filters or self.reflux or self.reweigh or self.config.additional_error # also after reweighing because chi squared terms are calculated during differences!

    # -----------------------------------------------------------------

    def get_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting differences ...")

        # Recalculate the differences
        if self.rediff: self.calculate_differences()

        # Load the original differences
        else: self.load_differences()

    # -----------------------------------------------------------------

    def load_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the differences between observed and simulated fluxes ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Loading the differences for the '" + generation_name + "' generation ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Loop over the simulation names
            #for simulation_name in generation.simulation_names:
            for simulation_name in self.simulation_names_for_generation(generation_name):

                # Individual simulations: differences files have been backupped and originals have already been removed
                if self.individual_simulations:

                    backup_filepath = fs.prepended_filepath(generation.get_simulation_sed_differences_path(simulation_name), self.backup_name, sep="__")
                    if not fs.is_file(backup_filepath): raise IOError("Differences backup file is not found for simulation '" + simulation_name + "'")
                    differences_filepath = backup_filepath

                else:

                    # Get the differences filepath
                    if not generation.has_sed_differences(simulation_name): raise IOError("Differences file is not found for simulation '" + simulation_name + "'")
                    differences_filepath = generation.get_simulation_sed_differences_path(simulation_name)

                # Debugging
                log.debug("Loading the differences for the '" + simulation_name + "' simulation ...")

                # Load
                #differences = generation.get_simulation_sed_differences(simulation_name)
                differences = FluxDifferencesTable.from_file(differences_filepath)

                # Set table
                self.differences[generation_name][simulation_name] = differences

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_sed(self):

        """
        This function ...
        :return:
        """

        # Get the original observed SED
        sed = self.observed_sed.copy()

        # Add the additional relative errors
        if self.config.additional_error is not None: sed.add_relative_error(self.config.additional_error)

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def get_fluxdensity(self, fltr, unit="Jy"):

        """
        This function ...
        :param fltr:
        :param unit:
        :return:
        """

        # Find the corresponding flux in the SED derived from observation
        #observed_fluxdensity = self.observed_sed.photometry_for_band(instrument, band, unit="Jy").value
        return self.reference_sed.photometry_for_filter(fltr, unit=unit)

    # -----------------------------------------------------------------

    def get_fluxdensity_error(self, fltr, unit="Jy"):

        """
        This function ...
        :param fltr:
        :param unit:
        :return:
        """

        # Find the corresponding flux error in the SED derived from observation
        #observed_fluxdensity_error = self.observed_sed.error_for_band(instrument, band, unit="Jy").average.to("Jy").value
        return self.reference_sed.error_for_filter(fltr, unit=unit)

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the differences between observed and simulated fluxes ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Calculating differences for the '" + generation_name + "' generation ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Loop over the simulations
            for simulation in generation.analysed_simulations_basic:

                # Get simulation name
                simulation_name = simulation.name
                if self.individual_simulations and simulation_name not in self.simulation_names: continue

                # Debugging
                log.debug("Calculating differences for the '" + simulation_name + "' simulation ...")

                # Initialize the differences table
                differences = FluxDifferencesTable()

                # Get mock SED
                if not generation.has_mock_sed(simulation_name): raise IOError("Mock SED is not found for simulation '" + simulation_name + "'")
                mock_sed = generation.get_mock_sed(simulation_name)

                # Loop over the entries in the fluxdensity table (SED) derived from the simulation
                #for i in range(len(mock_sed)):
                for fltr in self.filters:

                    # Find index in the mock SED
                    i = mock_sed.index_for_filter(fltr)

                    # Get instrument, band and flux density
                    instrument = mock_sed["Instrument"][i]
                    band = mock_sed["Band"][i]
                    #fluxdensity = mock_sed["Photometry"][i]
                    fluxdensity = mock_sed.get_photometry(i, unit="Jy")

                    # Find the corresponding flux in the SED derived from observation
                    #observed_fluxdensity = self.observed_sed.photometry_for_band(instrument, band, unit="Jy").value
                    # Find the corresponding flux error in the SED derived from observation
                    #observed_fluxdensity_error = self.observed_sed.error_for_band(instrument, band, unit="Jy").average.to("Jy").value
                    observed_fluxdensity = self.get_fluxdensity(fltr)
                    observed_fluxdensity_error = self.get_fluxdensity_error(fltr)

                    # If no match with (instrument, band) is found in the observed SED
                    if observed_fluxdensity is None:
                        log.warning("The observed flux density could not be found for the " + instrument + " " + band + " band")
                        continue

                    # Calculate the difference
                    #print(fluxdensity, observed_fluxdensity)
                    #observed_fluxdensity_value = observed_fluxdensity.value # is in Jy
                    difference = fluxdensity - observed_fluxdensity
                    relative_difference = float(difference / observed_fluxdensity)

                    # Find the index of the current band in the weights table
                    index = tables.find_index(self.weights, key=[instrument, band], column_name=["Instrument", "Band"])
                    if index is None: continue  # Skip this band if a weight is not found
                    weight = self.weights["Weight"][index]

                    # Calculate the chi squared term
                    difference_value = difference.to("Jy").value
                    error_value = observed_fluxdensity_error.average.to("Jy").value
                    chi_squared_term = weight * difference_value ** 2 / error_value ** 2

                    # Add entry to the table
                    differences.add_entry(instrument, band, difference, relative_difference, chi_squared_term)

                # Set table
                self.differences[generation_name][simulation_name] = differences

            # Show the number of simulations of this generation for which differences have been calculated
            nsimulations = len(self.differences[generation_name])
            log.debug("Flux differences have been calculated for " + str(nsimulations) + " simulations of generation '" + generation_name + "'")

    # -----------------------------------------------------------------

    @property
    def nfree_parameters(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.nfree_parameters

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating chi squared values ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Calculating chi squared values for the generation '" + generation_name + "' ...")

            # Get the chi squared table
            chi_squared_table = self.chi_squared_tables[generation_name]

            # Loop over the simulations
            for simulation_name in self.differences[generation_name]:

                # Get the differences table
                differences = self.differences[generation_name][simulation_name]

                # Debugging
                log.debug("Calculating the chi squared value for simulation '" + simulation_name + "' ...")

                # Calculate
                chi_squared = calculate_chi_squared_from_differences(differences, self.nfree_parameters)

                # Loop over the
                if simulation_name in chi_squared_table.simulation_names:

                    if not self.individual_simulations: raise RuntimeError("Something went wrong: chi squared table can only already be containing to be refitted simulations when refitting only select simulations for one generation")
                    else: chi_squared_table.set_chi_squared(simulation_name, chi_squared)

                # Add entry to the chi squared table
                else: chi_squared_table.add_entry(simulation_name, chi_squared)

    # -----------------------------------------------------------------

    def get_best_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the best parameter values ...")

        # Loop over the finished generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Getting the best parameter values for generation '" + generation_name + "' ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Get the name of the simulation with the lowest chi squared value
            # WHY WAS THIS USING THE ORIGINAL CHI SQUARED TABLE?
            #best_simulation_name = generation.chi_squared_table.best_simulation_name
            #chi_squared = generation.chi_squared_table.chi_squared_for(best_simulation_name)
            best_simulation_name = self.chi_squared_tables[generation_name].best_simulation_name
            chi_squared = self.chi_squared_tables[generation_name].chi_squared_for(best_simulation_name)

            # Get the parameter values
            values = generation.parameters_table.parameter_values_for_simulation(best_simulation_name)

            # Add an entry to the best parameters table file
            self.best_parameters_table.add_entry(generation_name, values, chi_squared)

    # -----------------------------------------------------------------

    def calculate_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the probabilities ...")

        # Get the probability tables
        self.get_model_probabilities()

        # Calculate the parameter probabilities for each generation separately
        self.calculate_parameter_probabilities()

        # Calculate the combined probabilities for all generations
        self.calculate_all_parameter_probabilities()

    # -----------------------------------------------------------------

    def get_simulation_names_parameters_and_chi_squared_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Load the parameter table
        parameter_table = self.fitting_run.parameters_table_for_generation(generation_name)

        # Get the NEW chi squared table
        chi_squared_table = self.chi_squared_tables[generation_name].copy()

        # Sort the table for decreasing chi squared value
        chi_squared_table.sort("Chi squared")
        chi_squared_table.reverse()

        # Get the chi squared values
        chi_squared_values = list(chi_squared_table["Chi squared"])

        # Initialize lists
        simulation_names = []
        parameter_values = []

        # Loop over the simulations
        for i in range(len(chi_squared_table)):

            # Get the simulation name
            simulation_name = chi_squared_table["Simulation name"][i]

            # Get a dictionary with the parameter values for this simulation
            values = parameter_table.parameter_values_for_simulation(simulation_name)

            # Add to lists
            simulation_names.append(simulation_name)
            parameter_values.append(values)

        # Return
        return simulation_names, parameter_values, chi_squared_values

    # -----------------------------------------------------------------

    def get_model_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the model probabilities for the generations ...")

        # Loop over the finished generations
        for generation_name in self.generation_names:

            # Get simulation names with parameter values and chi squared values
            simulation_names, parameter_values, chi_squared_values = self.get_simulation_names_parameters_and_chi_squared_for_generation(generation_name)
            nsimulations = len(simulation_names)

            # Calculate the probability for each model
            probabilities = np.exp(-0.5 * np.asarray(chi_squared_values))

            # Create the probabilities table
            probabilities_table = ModelProbabilitiesTable(parameters=self.free_parameter_labels, units=self.fitting_run.parameter_units)

            # Add the entries to the model probabilities table
            for i in range(nsimulations):

                # Get the simulation name
                simulation_name = simulation_names[i]

                # Get a dictionary with the parameter values for this simulation
                values = parameter_values[i]

                # Get the probability
                probability = probabilities[i]

                # Add an entry to the table
                probabilities_table.add_entry(simulation_name, values, probability)

            # Add to the dictionary
            self.model_probabilities[generation_name] = probabilities_table

    # -----------------------------------------------------------------

    @property
    def free_parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.free_parameter_labels

    # -----------------------------------------------------------------

    def calculate_parameter_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the probabilities of the different parameter values for each generation ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Initialize dictionary for this generation
            self.parameter_probabilities[generation_name] = dict()

            # Loop over the free parameters
            for label in self.free_parameter_labels:

                # Create a set for the unique values
                unique_values = set()

                # Initialize a ParameterProbabilitiesTable instance for this parameter
                table = ParameterProbabilitiesTable()

                # Loop over the values of this parameter for this generation, and expand the set accordingly
                for value in self.model_probabilities[generation_name][label]: unique_values.add(value)

                # Get a (sorted) list of all the unique values for this parameter
                unique_values = sorted(list(unique_values))

                # Add an entry for each unique parameter value that has been encountered
                for value in unique_values:

                    # Get an array of the probabilities of all the models that have this unique value
                    simulation_indices = self.model_probabilities[generation_name][label] == value
                    individual_probabilities = self.model_probabilities[generation_name]["Probability"][simulation_indices]

                    # Combine the individual probabilities
                    combined_probability = np.sum(np.asarray(individual_probabilities))

                    # Add an entry to the table
                    table.add_entry(value, combined_probability)

                # Set the table
                self.parameter_probabilities[generation_name][label] = table

    # -----------------------------------------------------------------

    def calculate_all_parameter_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the probabilities of the different parameter values for all generations ...")

        # Loop over the free parameters
        for label in self.free_parameter_labels:

            # Create a set for the unique values
            unique_values = set()

            # Loop over the generations to extract all unique values for the parameter
            for generation_name in self.model_probabilities:

                # Loop over the values of this parameter for this generation, and expand the set accordingly
                for value in self.model_probabilities[generation_name][label]: unique_values.add(value)

            # Get a (sorted) list of all the unique values for this parameter
            unique_values = sorted(list(unique_values))

            # Initialize a ParameterProbabilitiesTable instance for this parameter
            table = ParameterProbabilitiesTable()

            # Add an entry for each unique parameter value that has been encountered
            for value in unique_values:

                # Add the probabilities from all models that have this value
                individual_probabilities = []

                # Loop over the generations
                for generation_name in self.model_probabilities:
                    simulation_indices = self.model_probabilities[generation_name][label] == value
                    individual_probabilities += list(self.model_probabilities[generation_name]["Probability"][simulation_indices])

                # Combine the individual probabilities
                combined_probability = np.sum(np.array(individual_probabilities))

                # Add an entry to the table
                table.add_entry(value, combined_probability)

            # Set the table
            self.parameter_probabilities_all[label] = table

    # -----------------------------------------------------------------

    def create_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the probability distributions ...")

        # Loop over the free parameters
        for label in self.free_parameter_labels:

            # Debugging
            log.debug("Creating distribution for the '" + label + "' parameter ...")

            # Convert the probability lists into NumPy arrays and normalize them
            normalized_probabilities = np.array(self.parameter_probabilities_all[label]["Probability"]) / sum(self.parameter_probabilities_all[label]["Probability"])

            # Create the probability distributions for the different parameters
            self.distributions[label] = Distribution.from_probabilities(label, normalized_probabilities, self.parameter_probabilities_all[label]["Value"])

    # -----------------------------------------------------------------

    def get_best_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the best simulations for each generation ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Getting best " + str(self.config.nbest) + " simulations of generation '" + generation_name + "' ...")

            # Get the chi squared table
            chi_squared = self.chi_squared_tables[generation_name]
            parameters = self.fitting_run.parameters_table_for_generation(generation_name)

            # Show
            simulation_names, counts = get_best_simulations(self.config.nbest, self.free_parameter_labels, chi_squared, parameters, self.parameter_units)

            # Set
            self.best_simulations[generation_name] = simulation_names
            self.best_simulation_counts[generation_name] = counts

    # -----------------------------------------------------------------

    @lazyproperty
    def best_simulation_counts_distributions(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        counts_distributions = defaultdict(dict)

        # Loop over the generations
        for generation_name in self.generation_names:

            # Get counts for this generation
            counts = self.best_simulation_counts[generation_name]

            # Make counts distributions
            for label in self.free_parameter_labels:

                # Create distribution
                counts_distributions[generation_name][label] = Distribution.from_counts(label, counts[label].values(), counts[label].keys(), sort=True)

        # Return
        return counts_distributions

    # -----------------------------------------------------------------

    def create_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the fitting run configuration ...")

        # Get original fitting configuration
        self.new_fitting_config = self.fitting_run.fitting_configuration

        # Adapt fitting filters?
        if self.different_filters: self.new_fitting_config.filters = self.filters
        #else: log.warning("No diff filters") #raise RuntimeError("No diff filters")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the fitting configuration
        self.write_config()

        # Write the runs table
        if self.as_run: self.write_table()

        # Weights
        self.write_weights()

        # Write the observed SED
        self.write_sed()

        # Fluxes?
        if self.reflux: self.write_fluxes()

        # Differences
        self.write_differences()

        # Chi squared tables
        self.write_chi_squared()

        # Write the model probabilities
        self.write_model_probabilities()

        # Write the parameter probabilities
        self.write_parameter_probabilities()

        # Write all generations parameter probabilities
        self.write_all_parameter_probabilities()

        # Write the ski file of the best simulation
        self.write_best_parameters()

        # Write the probability distributions in table format
        self.write_distributions()

    # -----------------------------------------------------------------

    @property
    def new_fitting_config_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.new_run_path, "configuration.cfg")

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the new fitting configuration ...")

        # Determine the path
        if self.as_run: path = self.new_fitting_config_path
        elif self.in_place:
            if self.different_filters: path = self.fitting_run.fitting_configuration_path
            else: return # not necessary
        else: path = fs.join(self.path, "configuration.cfg")

        # Save the config
        self.new_fitting_config.saveto(path)

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the runs table ...")

        # Open the runs table, add the new run and save
        table = self.runs_table
        table.add_run(self.new_run)
        table.save()

    # -----------------------------------------------------------------

    @property
    def weights_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "weights.dat")

    # -----------------------------------------------------------------

    def write_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with weights ...")

        # Determine the path
        if self.as_run: path = fs.join(self.new_run_path, "weights.dat")
        elif self.in_place:
            if self.reweigh: path = self.fitting_run.weights_table_path
            else: return # not necessary
        else: path = self.weights_table_path

        # Write the table with weights
        self.weights.saveto(path)

    # -----------------------------------------------------------------

    def write_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the observed SED ...")

        # Determine the path
        if self.as_run: path = fs.join(self.new_run_path, "observed_sed.dat")
        elif self.in_place:
            if self.individual_simulations: return # don't write the changed SED (additional errors)
            else: path = fs.join(self.fitting_run.path, "observed_sed.dat")
        else: path = fs.join(self.path, "observed_sed.dat")

        # Write the reference SED
        self.reference_sed.saveto(path)

    # -----------------------------------------------------------------

    def write_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the mock fluxes ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Writing the mock fluxes of generation '" + generation_name + "' ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Determine fluxes directory name
            if generation.use_images: directory_name = "image fluxes"
            else: directory_name = "fluxes"

            # Loop over the simulations
            for simulation_name in self.simulation_names_for_generation(generation_name):

                # Check
                if simulation_name not in self.fluxes[generation_name]:
                    log.warning("No fluxes for simulation '" + simulation_name + "' of generation '" + generation_name + "'")
                    continue

                # Debugging
                log.debug("Writing the mock fluxes of simulation '" + simulation_name + "' ...")

                # Get the fluxes
                fluxes = self.fluxes[generation_name][simulation_name]

                # Determine the path
                if self.as_run: path = fs.join(self.new_simulation_misc_paths[generation_name][simulation_name], directory_name, "fluxes_earth.dat")
                elif self.in_place: path = fs.join(self.fitting_run.get_generation_path(generation_name), simulation_name, "misc", directory_name, "fluxes_earth.dat")
                else: path = fs.join(self.simulation_paths[generation_name][simulation_name], "fluxes_earth.dat")

                # Save
                fluxes.saveto(path)

    # -----------------------------------------------------------------

    def write_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the flux difference tables ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Writing the flux difference tables of generation '" + generation_name + "' ...")

            # Loop over the simulations
            for simulation_name in self.simulation_names_for_generation(generation_name):

                # Check
                if simulation_name not in self.differences[generation_name]:
                    log.warning("No differences for simulation '" + simulation_name + "' of generation '" + generation_name + "'")
                    continue

                # Debugging
                log.debug("Writing the flux difference table of simulation '" + simulation_name + "' ...")

                # Get the differences table
                differences = self.differences[generation_name][simulation_name]

                # Determine the path
                if self.as_run: path = fs.join(self.new_simulation_misc_paths[generation_name][simulation_name], "differences.dat")
                elif self.in_place: path = fs.join(self.fitting_run.get_generation_path(generation_name), simulation_name, "misc", "differences.dat")
                else: path = fs.join(self.simulation_paths[generation_name][simulation_name], "differences.dat")

                # Save
                differences.saveto(path)

    # -----------------------------------------------------------------

    def write_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the chi squared tables ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Writing the chi squared table of generation '" + generation_name + "' ...")

            # Get the chi squared table
            chi_squared_table = self.chi_squared_tables[generation_name]

            # Determine the path
            if self.as_run: path = fs.join(self.new_generation_paths[generation_name], "chi_squared.dat")
            elif self.in_place: path = fs.join(self.fitting_run.get_generation_path(generation_name), "chi_squared.dat")
            else: path = self.chi_squared_table_paths[generation_name]

            # Save the table
            chi_squared_table.saveto(path)

    # -----------------------------------------------------------------

    def get_model_probabilities_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        if self.as_run: return fs.join(self.new_prob_generation_paths[generation_name], "models.dat")
        elif self.in_place: return fs.join(self.fitting_run.prob_generations_path, generation_name, "models.dat")
        else: return fs.join(self.prob_generations_paths[generation_name], "models.dat")

    # -----------------------------------------------------------------

    def write_model_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the model probabilities for each generation ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Determine the path
            path = self.get_model_probabilities_table_path_for_generation(generation_name)

            # Check whether the directory exists for this generation
            dirpath = fs.directory_of(path)
            if not fs.is_directory(dirpath): fs.create_directory(dirpath, recursive=True)

            # Get the table
            table = self.model_probabilities[generation_name]

            # Save the table
            table.saveto(path)

    # -----------------------------------------------------------------

    def get_parameter_probabilities_table_path_for_generation(self, generation_name, parameter_label):

        """
        This function ...
        :param generation_name:
        :param parameter_label:
        :return:
        """

        if self.as_run: return fs.join(self.new_prob_generation_paths[generation_name], parameter_label + ".dat")
        elif self.in_place: return fs.join(self.fitting_run.prob_generations_path, generation_name, parameter_label + ".dat")
        else: return fs.join(self.prob_generations_paths[generation_name], parameter_label + ".dat")

    # -----------------------------------------------------------------

    def write_parameter_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the parameter probabilities for each generation ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Loop over the free parameters
            for label in self.free_parameter_labels:

                # Get the parameter probabilities
                probabilities = self.parameter_probabilities[generation_name][label]

                # Determine the path
                path = self.get_parameter_probabilities_table_path_for_generation(generation_name, label)

                # Save the table
                probabilities.saveto(path)

    # -----------------------------------------------------------------

    def get_description(self, parameter_label):

        """
        This function ...
        :param parameter_label:
        :return:
        """

        return self.fitting_run.parameter_descriptions[parameter_label]

    # -----------------------------------------------------------------

    def get_parameter_probabilities_table_path(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Determine the path for the table
        if self.as_run: return fs.join(self.new_prob_parameters_path, label + ".dat")
        elif self.in_place: return fs.join(self.fitting_run.prob_parameters_path, label + ".dat")
        else: return fs.join(self.prob_parameters_path, label + ".dat")

    # -----------------------------------------------------------------

    def write_all_parameter_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the parameter probabilities from all generations ...")

        # Loop over the probability tables for the different free parameter
        for label in self.fitting_run.free_parameter_labels:

            # Debugging
            log.debug("Writing the parameter probabilities for the " + self.get_description(label) + " ...")

            # Get path
            path = self.get_parameter_probabilities_table_path(label)

            # Save the table
            self.parameter_probabilities_all[label].saveto(path)

    # -----------------------------------------------------------------

    def write_best_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the best model parameters table ...")

        # Determine the path
        if self.as_run: path = fs.join(self.new_run_path, "best_parameters.dat")
        elif self.in_place: path = fs.join(self.fitting_run.best_parameters_table_path)
        else: path = fs.join(self.path, "best_parameters.dat")

        # Save the best parameters table
        self.best_parameters_table.saveto(path)

    # -----------------------------------------------------------------

    def get_parameter_distribution_table_path(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Determine the path for the table
        if self.as_run: return fs.join(self.new_prob_distributions_path, label + ".dat")
        elif self.in_place: return fs.join(self.fitting_run.prob_distributions_path, label + ".dat")
        else: return fs.join(self.prob_distributions_path, label + ".dat")

    # -----------------------------------------------------------------

    def write_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the probability distributions ...")

        # Loop over the entries in the 'probabilities' table
        for label in self.distributions:

            # Debugging
            log.debug("Writing the probability distribution of the " + self.get_description(label) + " ...")

            # Get path
            path = self.get_parameter_distribution_table_path(label)

            # Write the table of probabilities for this parameter
            self.distributions[label].saveto(path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Show the best simulations
        self.show_best()

        # Show statistics
        self.show_statistics()

    # -----------------------------------------------------------------

    @property
    def initial_parameter_values(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.first_guess_parameter_values

    # -----------------------------------------------------------------

    @property
    def parameter_ranges(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.free_parameter_ranges

    # -----------------------------------------------------------------

    @property
    def parameter_units(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.parameter_units

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_scales(self):

        """
        This function ...
        :return:
        """

        # Get the scales
        grid_settings = self.fitting_run.grid_settings
        parameter_scales = dict()
        for label in self.free_parameter_labels:
            key = label + "_scale"
            parameter_scales[label] = grid_settings[key]
        return parameter_scales

    # -----------------------------------------------------------------

    def show_best(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing best simulations of each generation ...")

        # Show initial parameter values
        print("")
        print("Initial parameter values:")
        print("")
        for label in self.free_parameter_labels: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.initial_parameter_values[label]))

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Showing best simulations of generation '" + generation_name + "' ...")

            # Show ranges
            print("")
            print("Parameter ranges:")
            print("")
            for label in self.free_parameter_labels: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.parameter_ranges[label]))
            print("")

            # Get the best simulation names
            simulation_names = self.best_simulations[generation_name]

            # Get the chi squared table
            chi_squared_table = self.chi_squared_tables[generation_name]
            parameters_table = self.fitting_run.parameters_table_for_generation(generation_name)

            # Show the best simulations
            show_best_simulations(simulation_names, chi_squared_table, parameters_table, self.parameter_units, self.parameter_scales, self.initial_parameter_values)
            print("")

            # Loop over the simulations
            for index, simulation_name in enumerate(simulation_names):

                # Determine the plot path
                sed_plot_path = self.generations[generation_name].get_simulation_sed_plot_path(simulation_name)

                # Check whether the file is present
                if not fs.is_file(sed_plot_path):
                    log.warning("SED plot file for simulation '" + simulation_name + "' is not present: skipping ...")
                    continue

                # Copy the SED plot file to the generation path
                if not self.as_run:

                    # Determine filename
                    filename = "best_" + str(index) + "_" + simulation_name + ".pdf"

                    # Copy plot file
                    if self.as_run: fs.copy_file(sed_plot_path, self.new_generation_paths[generation_name], new_name=filename)
                    elif self.in_place: fs.copy_file(sed_plot_path, self.fitting_run.get_generation_path(generation_name), new_name=filename)
                    else: fs.copy_file(sed_plot_path, self.generation_paths[generation_name], new_name=filename)

                # Show the plot
                if self.config.show_best_sed: fs.open_file(sed_plot_path)

            # Show counts
            if self.config.show_counts:

                # Get the counts
                counts = self.best_simulation_counts[generation_name]

                print("Counts in best simulations:")
                for label in self.free_parameter_labels:

                    print("")
                    print(" - " + fmt.bold + label + fmt.reset + ":")
                    print("")
                    for value in sorted(counts[label].keys()):

                        count = counts[label][value]
                        relcount = float(count) / self.config.nbest
                        print("    * " + tostr(value) + ": " + str(counts[label][value]) + " (" + tostr(relcount * 100) + "%)")

                print("")

            # Plot the distributions
            if self.config.plot_counts:

                # Get the distributions
                counts_distributions = self.best_simulation_counts_distributions[generation_name]

                # Determine path
                if self.as_run: path = fs.join(self.new_generation_paths[generation_name], "best_counts.pdf")
                elif self.in_place: path = fs.join(self.fitting_run.get_generation_path(generation_name), "best_counts.pdf")
                else: path = fs.join(self.generation_paths[generation_name], "best_counts.pdf")

                # Plot
                plot_distributions(counts_distributions, logscale=True, panels=True, frequencies=True, path=path)

    # -----------------------------------------------------------------

    def show_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing statistics for each generation ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Showing statistics for generation '" + generation_name + "' ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Get the chi squared table and parameters table
            chi_squared = self.chi_squared_tables[generation_name]
            parameters = self.fitting_run.parameters_table_for_generation(generation_name)

            # Get best simulation name and chi squared
            best_simulation_name, best_chi_squared = chi_squared.best_simulation_name_and_chi_squared

            # Get best simulation parameter values
            best_parameter_values = parameters.parameter_values_for_simulation(best_simulation_name)

            # Most probable model: should be same as simulation with lowest chi squared
            most_probable_simulation_name = self.model_probabilities[generation_name].most_probable_simulation
            assert most_probable_simulation_name == best_simulation_name

            # Show whether or not the best simulation corresponds to the best simulation of the original fit
            if most_probable_simulation_name == generation.most_probable_model: log.success("The best simulation from this refitting is the same as for the original fit (" + most_probable_simulation_name + ")")
            else: log.success("The best simulation from this refitting (" + most_probable_simulation_name + ") is different from the best simulation of the original fit (" + generation.most_probable_model + ")")

            print("")
            print("Statistics:")

            # Get the counts distributions
            counts_distributions = self.best_simulation_counts_distributions[generation_name]

            # Loop over the free parameter labels
            for label in self.free_parameter_labels:

                print("")
                print(" - " + fmt.bold + label + fmt.reset + ":")
                print("")

                # Get most probable parameter value
                most_probable_value = generation.get_most_probable_parameter_value(label)

                print("    * Initial guess value: " + tostr(self.initial_parameter_values[label], ndigits=3))
                print("    * Best simulation value: " + tostr(best_parameter_values[label], ndigits=3))
                print("    * Most probable value: " + tostr(most_probable_value, ndigits=3) + " " + tostr(self.parameter_units[label]))
                if self.config.nbest > 1: print("    * Most counted in " + str(self.config.nbest) + " best simulations: " + tostr(counts_distributions[label].most_frequent, ndigits=3) + " " + tostr(self.parameter_units[label]))
            print("")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Chi squared values
        if self.config.plot_chi_squared: self.plot_chi_squared()

        # Plot the parameter probabilities per generation
        if self.config.plot_probabilities: self.plot_probabilities()

        # Plot the probability distributions as histograms
        if self.config.plot_distributions: self.plot_distributions()

    # -----------------------------------------------------------------

    def plot_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting chi squared tables ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Plotting chi squared table for generation '" + generation_name + "' ...")

            # Set title
            title = "Chi squared values of generation '" + generation_name.replace("_", "\_") + "'"

            # Determine path
            if self.as_run: path = fs.join(self.new_generation_paths[generation_name], "chi_squared.pdf")
            elif self.in_place: path = fs.join(self.fitting_run.get_generation_path(generation_name), "chi_squared.pdf")
            else: path = fs.join(self.generation_paths[generation_name], "chi_squared.pdf")

            # Plot chi squared?
            plot_distribution(self.chi_squared_tables[generation_name].distribution, title=title, logscale=True, path=path)

    # -----------------------------------------------------------------

    def plot_probabilities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the parameter probabilities for each generation ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Determine the path
            if self.as_run: path = fs.join(self.new_prob_generation_paths[generation_name], "probabilities.pdf")
            elif self.in_place: path = fs.join(self.fitting_run.prob_generations_path, generation_name, "probabilities.pdf")
            else: path = fs.join(self.prob_generations_paths[generation_name], "probabilities.pdf")

            # Make distributions
            distributions = dict()
            for label in self.free_parameter_labels:
                distribution = Distribution.from_probabilities(label, self.parameter_probabilities[generation_name][label]["Probability"], self.parameter_probabilities[generation_name][label]["Value"])
                distribution.normalize(method="sum")
                distributions[label] = distribution

            # Plot in different panels
            plot_distributions(distributions, panels=True, extrema=True, frequencies=True, path=path, logscale=True)

    # -----------------------------------------------------------------

    def plot_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the probability distributions ...")

        # Determine the path
        if self.as_run: path = fs.join(self.new_prob_distributions_path, "distributions.pdf")
        elif self.in_place: path = fs.join(self.fitting_run.prob_distributions_path, "distributions.pdf")
        else: path = fs.join(self.prob_distributions_path, "distributions.pdf")

        # Plot in different panels
        plot_distributions(self.distributions, panels=True, extrema=True, frequencies=True, path=path, logscale=True)

# -----------------------------------------------------------------

def get_best_simulations(nsimulations, parameter_labels, chi_squared_table, parameters_table, parameter_units, return_unique_values_scalar=False):

    """
    This function ...
    :param nsimulations:
    :param parameter_labels:
    :param chi_squared_table:
    :param parameters_table:
    :param parameter_units:
    :param return_unique_values_scalar:
    :return:
    """

    # Initialize the counts dictionary
    counts = dict()
    for label in parameter_labels: counts[label] = defaultdict(int)

    # Create list of unique values for each free parameter
    unique_values = parameters_table.unique_parameter_values
    unique_values_scalar = dict()
    for label in unique_values:
        values = list(sorted([value.to(parameter_units[label]).value for value in unique_values[label]]))
        unique_values_scalar[label] = values

    # Fill with zeros
    for label in parameter_labels:
        for value in unique_values_scalar[label]: counts[label][value] = 0

    # Get best simulation names
    simulation_names = chi_squared_table.get_best_simulation_names(nsimulations)

    # Loop over the simulations
    for index, simulation_name in enumerate(simulation_names):

        # Get parameter values
        parameter_values = parameters_table.parameter_values_for_simulation(simulation_name)

        # Get the counts
        for label in parameter_values:

            # Get properties
            unit = parameter_units[label]
            value = parameter_values[label]
            value_scalar = value.to(unit).value

            # Add count for this parameter
            counts[label][value_scalar] += 1

    # Return the simulation names
    if return_unique_values_scalar: return simulation_names, counts, unique_values_scalar
    else: return simulation_names, counts

# -----------------------------------------------------------------

def show_best_simulations(simulation_names, chi_squared_table, parameters_table, parameter_units, parameter_scales, initial_values):

    """
    This function ...
    :param simulation_names:
    :param chi_squared_table:
    :param parameters_table:
    :param parameter_units:
    :param parameter_scales:
    :param initial_values:
    :return:
    """

    # Create list of unique parameter values for each free parameter
    unique_values = parameters_table.unique_parameter_values
    unique_values_scalar = dict()
    for label in unique_values:
        values = list(sorted([value.to(parameter_units[label]).value for value in unique_values[label]]))
        unique_values_scalar[label] = values

    # Create list for chi squared values
    chi_squared_values = []

    # Create list for parameter values
    parameters = []

    # Loop over the simulations
    for index, simulation_name in enumerate(simulation_names):

        # Get the chi squared value
        chisq = chi_squared_table.chi_squared_for(simulation_name)

        # Get parameter values
        parameter_values = parameters_table.parameter_values_for_simulation(simulation_name)

        # Add the values
        chi_squared_values.append(chisq)
        parameters.append(parameter_values)

    # Show
    show_best_simulations_impl(simulation_names, chi_squared_values, parameters, parameter_units, unique_values_scalar, parameter_scales, initial_values)

# -----------------------------------------------------------------

def show_best_simulations_impl(simulation_names, chi_squared_values, parameters, parameter_units, unique_values_scalar, parameter_scales, initial_values):

    """
    This function ...
    :param simulation_names:
    :param chi_squared_values:
    :param parameters:
    :param unique_values_scalar:
    :param parameter_scales:
    :param initial_values:
    :return:
    """

    # Checks
    nsimulations = len(simulation_names)
    if len(chi_squared_values) != nsimulations: raise ValueError("Invalid number of chi squared values")
    if len(parameters) != nsimulations: raise ValueError("Invalid number of parameter dictionaries")

    # Loop over the simulations
    for index, simulation_name, chisq, parameter_values in zip(range(nsimulations), simulation_names, chi_squared_values, parameters):

        # Show
        print("")
        print(" " + str(index + 1) + ". " + fmt.green + fmt.underlined + simulation_name + fmt.reset)
        print("")

        # Show chi squared and parameter values
        print("  - chi squared: " + str(chisq))
        print("  - parameters:")
        for label in parameter_values:

            # Get properties
            unit = parameter_units[label]
            initial_value = initial_values[label]
            initial_value_scalar = initial_value.to(unit).value
            value = parameter_values[label]
            value_scalar = value.to(unit).value
            nunique_values = len(unique_values_scalar[label])

            # Get index of value and index of initial value
            i = unique_values_scalar[label].index(value_scalar)
            j = nr.locate_continuous(unique_values_scalar[label], initial_value_scalar, scale=parameter_scales[label])
            if j is not None:
                j = int(round(j))
                if not numbers.is_integer(j, absolute=False): raise NotImplementedError("Not yet implemented")
                #j = int(round(j))

            # Show
            indicator = "[ "
            for _ in range(nunique_values):
                if j is not None and _ == j: character = "+"
                else: character = "o"
                if _ == i: indicator += fmt.red + character + fmt.reset + " "
                else: indicator += character + " "
            indicator += "]"

            # Show
            print("     * " + fmt.bold + label + fmt.reset + ": " + tostr(value) + "   " + indicator)

# -----------------------------------------------------------------

def get_and_show_best_simulations(nsimulations, parameter_labels, chi_squared_table, parameters_table, parameter_units, parameter_scales, initial_values):

    """
    This function ...
    :param nsimulations:
    :param parameter_labels:
    :param chi_squared_table:
    :param parameters_table:
    :param parameter_units:
    :param parameter_scales:
    :param initial_values:
    :return:
    """

    # Get best simulations and counts
    simulation_names, counts, unique_values_scalar = get_best_simulations(nsimulations, parameter_labels, chi_squared_table, parameters_table, parameter_units, return_unique_values_scalar=True)

    # Create list for chi squared values
    chi_squared_values = []

    # Create list for parameter values
    parameters = []

    # Loop over the simulations
    for index, simulation_name in enumerate(simulation_names):

        # Get the chi squared value
        chisq = chi_squared_table.chi_squared_for(simulation_name)

        # Get parameter values
        parameter_values = parameters_table.parameter_values_for_simulation(simulation_name)

        # Add the values
        chi_squared_values.append(chisq)
        parameters.append(parameter_values)

    # Show best simulations
    show_best_simulations_impl(simulation_names, chi_squared_values, parameters, parameter_units, unique_values_scalar, parameter_scales, initial_values)

    # Return the simulation names
    return simulation_names, counts

# -----------------------------------------------------------------

def clone_fitting_run(fitting_run, new_run_name, generations=None, new_prob_generation_paths=None,
                      new_generation_paths=None, new_simulation_paths=None, new_simulation_misc_paths=None,
                      clone_configuration=True):

    """
    This function ...
    :param fitting_run:
    :param new_run_name:
    :param generations:
    :param new_prob_generation_paths: supposed to be EMPTY dict
    :param new_generation_paths: supposed to be EMPTY dict
    :param new_simulation_paths: supposed to be EMPTY dict
    :param new_simulation_misc_paths: supposed to be EMPTY dict
    :param clone_configuration:
    :return:
    """

    # Get modeling object name
    object_name = fitting_run.object_name

    # Generate new fitting run directory
    new_run_path = fs.create_directory_in(fitting_run.fit_path, new_run_name)

    # Get the model name
    model_name = fitting_run.model_name

    # Create the fitting run
    new_run = FittingRun(new_run_path, new_run_name, model_name, passive=True)

    # Clone the fitting configuration
    if clone_configuration: new_configuration_path = fs.copy_file(fitting_run.fitting_configuration_path, new_run_path)
    else: new_configuration_path = None

    # Create directories
    new_best_path = fs.create_directory_in(new_run_path, "best")
    new_generations_path = fs.create_directory_in(new_run_path, "generations")
    new_prob_path = fs.create_directory_in(new_run_path, "prob")

    # Create prob subdirectories
    new_prob_generations_path = fs.create_directory_in(new_prob_path, "generations")
    new_prob_parameters_path = fs.create_directory_in(new_prob_path, "parameters")
    new_prob_distributions_path = fs.create_directory_in(new_prob_path, "distributions")

    # Create prob generation paths
    if generations is not None:
        for generation_name in generations:
            path = fs.create_directory_in(new_prob_generations_path, generation_name)
            if new_prob_generation_paths is not None: new_prob_generation_paths[generation_name] = path

    # Copy directories
    new_geometries_path = fs.copy_directory(fitting_run.geometries_path, new_run_path)
    new_wavelength_grids_path = fs.copy_directory(fitting_run.wavelength_grids_path, new_run_path)

    # Copy template ski file
    fs.copy_file(fitting_run.template_ski_path, new_run_path)

    # Copy generations table
    new_generations_table_path = fs.copy_file(fitting_run.generations_table_path, new_run_path)

    # Open the generations table
    generations_table = GenerationsTable.from_file(new_generations_table_path)

    # Remove rows
    if generations is not None: generations_table.remove_other_entries(generations)
    else: generations_table.remove_all_rows()

    # Save the table
    generations_table.save()

    # Copy timing and memory tables
    fs.copy_file(fitting_run.timing_table_path, new_run_path)
    fs.copy_file(fitting_run.memory_table_path, new_run_path)

    # Copy input maps file
    fs.copy_file(fitting_run.input_maps_file_path, new_run_path)

    # If generation names are given
    if generations is not None:

        # Copy the generations
        for generation_name in generations:

            # Get original generation path
            generation_path = fitting_run.get_generation_path(generation_name)

            # Create generation directory
            new_generation_path = fs.create_directory_in(new_generations_path, generation_name)
            if new_generation_paths is not None: new_generation_paths[generation_name] = new_generation_path

            # Loop over the simulation names
            for path, simulation_name in fs.directories_in_path(generation_path, returns=["path", "name"]):

                # Make a new simulation directory
                new_path = fs.create_directory_in(new_generation_path, simulation_name)
                if new_simulation_paths is not None: new_simulation_paths[generation_name][simulation_name] = new_path

                # Copy ski file, output, plot and extr
                ski_path = fs.join(path, object_name + ".ski")
                out_path = fs.join(path, "out")
                extr_path = fs.join(path, "extr")
                plot_path = fs.join(path, "plot")
                fs.copy_file(ski_path, new_path)
                fs.copy_directory(out_path, new_path)
                fs.copy_directory(extr_path, new_path)
                fs.copy_directory(plot_path, new_path)

                # Create misc path
                misc_path = fs.create_directory_in(new_path, "misc")
                if new_simulation_misc_paths is not None: new_simulation_misc_paths[generation_name][simulation_name] = misc_path

            # Copy files in the generation directory
            fs.copy_files_from_directory(generation_path, new_generation_path, exact_not_name=["info", "chi_squared"])

            # Create the generation info
            info_path = fs.join(generation_path, "info.dat")
            info = GenerationInfo.from_file(info_path)

            # Set the new generation path
            info.path = new_generation_path

            # Save the new info
            new_info_path = fs.join(new_generation_path, "info.dat")
            info.saveto(new_info_path)

    # Return the new fitting run
    return new_run

# -----------------------------------------------------------------
