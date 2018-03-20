#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.statistics Contains the FittingStatistics class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.configurable import InteractiveConfigurable
from ...core.basics.configuration import ConfigurationDefinition, prompt_settings, parse_arguments
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty
from ..core.environment import load_modeling_environment
from ...core.tools import strings
from ...core.plot.distribution import plot_distribution
from ...core.tools.utils import memoize_method
from .sedfitting import chi_squared_table_to_probabilities
from ...core.basics.distribution import Distribution
from ...core.plot.sed import SEDPlotter
from ...core.tools import sequences
from .refitter import get_and_show_best_simulations

# -----------------------------------------------------------------

_help_command_name = "help"
_history_command_name = "history"
_status_command_name = "status"

_generations_command_name = "generations"
_best_command_name = "best"
_plot_command_name = "plot"

_terms_command_name = "terms"
_ranks_command_name = "ranks"
_chisquared_command_name = "chisquared"
_prob_command_name = "prob"
_seds_command_name = "seds"

# -----------------------------------------------------------------

# Define commands
commands = OrderedDict()

# Standard commands
commands[_help_command_name] = ("show_help", False, "show help", None)
commands[_history_command_name] = ("show_history_command", True, "show history of executed commands", None)
commands[_status_command_name] = ("show_status_command", True, "show generation status", None)

# Other commands
commands[_generations_command_name] = ("show_generations", False, "show generations", None)
commands[_best_command_name] = ("show_best_command", True, "show best models", None)
commands[_plot_command_name] = (None, None, "plotting", None)

# -----------------------------------------------------------------

# Plot commands
plot_commands = OrderedDict()
plot_commands[_terms_command_name] = ("plot_terms_command", True, "plot the chi squared terms", "generation")
plot_commands[_ranks_command_name] = ("plot_ranks_command", True, "plot the chi squared as a function of rank", "generation")
plot_commands[_chisquared_command_name] = ("plot_chi_squared_command", True, "plot the distribution of chi squared values", "generation")
plot_commands[_prob_command_name] = ("plot_probabilities_command", True, "plot the distribution of probabilities", "generation")
plot_commands[_seds_command_name] = ("plot_seds_command", True, "plot the simulation SEDs of a generation", "generation")

# Set subcommands
subcommands = OrderedDict()
subcommands[_plot_command_name] = plot_commands

# -----------------------------------------------------------------

class FittingStatistics(InteractiveConfigurable):
    
    """
    This class...
    """

    _commands = commands
    _subcommands = subcommands

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(FittingStatistics, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def do_commands(self):

        """
        This function ...
        :return:
        """

        return self.config.commands is not None and len(self.config.commands) > 0

    # -----------------------------------------------------------------

    @property
    def do_interactive(self):

        """
        This function ...
        :return:
        """

        return self.config.interactive

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Run commands
        if self.do_commands: self.run_commands()

        # 3. Interactive
        if self.do_interactive: self.interactive()

        # Show
        self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(FittingStatistics, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def get_generation_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("generation", "string", "generation name")

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_generation_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True,
                                interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        if command_definition is not None: parse_command = splitted[index:]
        else: parse_command = splitted[index:index + 1]  # only image name

        # Get the definition
        definition = self.get_generation_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get settings interactively
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False, add_config_path=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get generation name
        generation_name = config.pop("generation")

        # Return
        return splitted, generation_name, config

    # -----------------------------------------------------------------

    def get_generation_name_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, generation_name, config = self.parse_generation_command(command, name=name, interactive=interactive)

        # Return the image name
        return generation_name

    # -----------------------------------------------------------------

    def get_generation_name_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, generation_name, config = self.parse_generation_command(command, command_definition=command_definition, name=name, interactive=interactive)

        # Return image name and config
        return generation_name, config

    # -----------------------------------------------------------------

    @lazyproperty
    def show_status_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def show_status_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.show_status_definition, **kwargs)

        # Show
        self.show_status(extra=config.extra_columns, path=config.path, refresh=config.refresh)

    # -----------------------------------------------------------------

    def show_status(self, extra=None, path=None, refresh=None):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the generation status ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def modeling_environment(self):

        """
        This function ...
        :return:
        """

        return load_modeling_environment(self.config.path)

    # -----------------------------------------------------------------

    @property
    def clipped_sed(self):

        """
        This function ...
        :return:
        """

        return self.modeling_environment.observed_sed

    # -----------------------------------------------------------------

    @property
    def truncated_sed(self):

        """
        This function ...
        :return:
        """

        return self.modeling_environment.truncated_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_runs(self):

        """
        This function ...
        :return:
        """

        return self.modeling_environment.fitting_runs

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.fitting_runs.load(self.config.run)

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

    @property
    def generation_names(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.generation_names

    # -----------------------------------------------------------------

    @memoize_method
    def get_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.fitting_run.get_generation(generation_name)

    # -----------------------------------------------------------------

    def get_simulation_names(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_generation(generation_name).simulation_names

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_generation(generation_name).get_simulation_or_basic(simulation_name)

    # -----------------------------------------------------------------

    def get_simulation_output(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_generation(generation_name).get_simulation_output(simulation_name)

    # -----------------------------------------------------------------

    def get_extraction_output(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_generation(generation_name).get_extraction_output(simulation_name)

    # -----------------------------------------------------------------

    def get_plotting_output(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_generation(generation_name).get_plotting_output(simulation_name)

    # -----------------------------------------------------------------

    def get_misc_output(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_generation(generation_name).get_misc_output(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_flux_differences(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_generation(generation_name).has_misc_differences(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def get_flux_differences(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_generation(generation_name).get_simulation_misc_differences(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def get_flux_differences_wavelengths(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_flux_differences(generation_name, simulation_name).wavelengths(unit="micron", add_unit=False)

    # -----------------------------------------------------------------

    @memoize_method
    def get_flux_differences_terms(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_flux_differences(generation_name, simulation_name).chi_squared_terms()

    # -----------------------------------------------------------------

    def get_chi_squared_terms_distribution(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        # Get the values
        wavelengths = self.get_flux_differences_wavelengths(generation_name, simulation_name)
        terms = self.get_flux_differences_terms(generation_name, simulation_name)

        # Create distribution over the wavelength range
        return Distribution.from_columns("Wavelength", wavelengths, terms, unit="micron", y_name="Weighed squared difference")

    # -----------------------------------------------------------------

    @memoize_method
    def has_sed(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_generation(generation_name).has_simulation_sed(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def get_sed(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_generation(generation_name).get_simulation_sed(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_mock_sed(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_generation(generation_name).has_mock_sed(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def get_mock_sed(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return self.get_generation(generation_name).get_mock_sed(simulation_name)

    # -----------------------------------------------------------------

    def get_chi_squared_table(self, generation_name):

        """
        This fucntion ...
        :param generation_name:
        :return:
        """

        return self.get_generation(generation_name).chi_squared_table

    # -----------------------------------------------------------------

    @memoize_method
    def get_chi_squared_values(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_chi_squared_table(generation_name).chi_squared_values

    # -----------------------------------------------------------------

    @memoize_method
    def get_chi_squared_distribution(self, generation_name, nbins=50, clip_above=None):

        """
        This function ...
        :param generation_name:
        :param nbins:
        :param clip_above:
        :return:
        """

        # Get the values
        values = self.get_chi_squared_values(generation_name)

        # Create and return distribution
        return Distribution.from_values("Chi squared", values, nbins=nbins, clip_above=clip_above, density=False)

    # -----------------------------------------------------------------

    @memoize_method
    def get_chi_squared_ranks_distribution(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Get the values
        values = self.get_chi_squared_values(generation_name)

        # Create and return the distribution
        return Distribution.by_rank("Simulation rank", values, y_name="Chi squared")

    # -----------------------------------------------------------------

    @memoize_method
    def get_probabilities_table(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Get tables
        chi_squared = self.get_chi_squared_table(generation_name)
        parameters = self.get_parameters_table(generation_name)

        # Create the probabilities table and return it
        return chi_squared_table_to_probabilities(chi_squared, parameters)

    # -----------------------------------------------------------------

    @memoize_method
    def get_probabilities(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_probabilities_table(generation_name).probabilities

    # -----------------------------------------------------------------

    @memoize_method
    def get_probability_distribution(self, generation_name, nbins=50, clip_below=None):

        """
        This function ...
        :param generation_name:
        :param nbins:
        :param clip_below:
        :return:
        """

        # Get the values
        values = self.get_probabilities(generation_name)

        # Create and return the distribution
        return Distribution.from_values("Probability", values, nbins=nbins, ignore_value=0, clip_below=clip_below, density=False)

    # -----------------------------------------------------------------

    @memoize_method
    def get_probability_ranks_distribution(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Get the values
        values = self.get_probabilities(generation_name)

        # Create and return the distribution
        return Distribution.by_rank("Simulation rank", values, y_name="Probability")

    # -----------------------------------------------------------------

    def get_probabilities_nzeros(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_probabilities_table(generation_name).nzeros

    # -----------------------------------------------------------------

    def get_assignment_table(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_generation(generation_name).assignment_table

    # -----------------------------------------------------------------

    def get_individuals_table(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_generation(generation_name).individuals_table

    # -----------------------------------------------------------------

    def get_parameters_table(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_generation(generation_name).parameters_table

    # -----------------------------------------------------------------

    def get_scores_table(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_generation(generation_name).scores_table

    # -----------------------------------------------------------------

    def get_elitism_table(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_generation(generation_name).elitism_table

    # -----------------------------------------------------------------

    def get_recurrence_table(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_generation(generation_name).recurrence_table

    # -----------------------------------------------------------------

    def get_crossover_table(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_generation(generation_name).crossover_table

    # -----------------------------------------------------------------

    def get_model_probabilities_table(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.get_generation(generation_name).model_probabilities_table

    # -----------------------------------------------------------------

    def show_generations(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Loop over the generations
        for generation_name in self.generation_names: print(" - " + generation_name)
        print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def show_best_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("nsimulations", "positive_integer", "number of simulations to show", 5)
        return definition

    # -----------------------------------------------------------------

    def show_best_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :return:
        """

        # Get generation name
        generation_name, config = self.get_generation_name_and_config_from_command(command, self.show_best_definition, **kwargs)

        # Show
        self.show_best(generation_name, nsimulations=config.nsimulations)

    # -----------------------------------------------------------------

    def show_best(self, generation_name, nsimulations=5):

        """
        This function ...
        :param generation_name:
        :param nsimulations:
        :return:
        """

        # Debugging
        log.debug("Showing best " + str(nsimulations) + " simulations for generation '" + generation_name + "' ...")

        # Get parameter ranges
        parameter_ranges = fitting_run.free_parameter_ranges
        parameter_labels = parameter_ranges.keys()

        # Get parameter units
        parameter_units = fitting_run.parameter_units

        # Show ranges
        print("")
        print("Parameter ranges:")
        print("")
        for label in parameter_labels: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(parameter_ranges[label]))

        # -----------------------------------------------------------------

        # Get the scales
        grid_settings = fitting_run.grid_settings
        parameter_scales = dict()
        for label in parameter_labels:
            key = label + "_scale"
            parameter_scales[label] = grid_settings[key]

        # -----------------------------------------------------------------

        # Get the chi squared table
        chi_squared = generation.chi_squared_table

        # -----------------------------------------------------------------

        # Load parameters table
        parameters = generation.parameters_table

        # -----------------------------------------------------------------

        initial_parameter_values = fitting_run.first_guess_parameter_values

        # Show
        simulation_names, counts = get_and_show_best_simulations(nsimulations, parameter_labels, chi_squared,
                                                                 parameters, parameter_units, parameter_scales,
                                                                 initial_parameter_values)
        print("")

        # FOR PLOT SED
        # # Show SED?
        # if config.sed:
        #
        #     # Loop over the simulations
        #     for simulation_name in simulation_names:
        #         # Determine the plot path
        #         sed_plot_path = generation.get_simulation_sed_plot_path(simulation_name)
        #
        #         # Show the plot
        #         fs.open_file(sed_plot_path)

    # -----------------------------------------------------------------

    def show_counts(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Make counts distributions
        counts_distributions = dict()
        for label in parameter_labels: counts_distributions[label] = Distribution.from_counts(label,
                                                                                              counts[label].values(),
                                                                                              counts[label].keys(),
                                                                                              sort=True)

        # Show counts
        if config.counts:
            print("Counts in best simulations:")
            for label in parameter_labels:
                print("")
                print(" - " + fmt.bold + label + fmt.reset + ":")
                print("")
                for value in sorted(counts[label].keys()):
                    count = counts[label][value]
                    relcount = float(count) / config.nsimulations
                    print("    * " + tostr(value) + ": " + str(counts[label][value]) + " (" + tostr(
                        relcount * 100) + "%)")

        # Plot the distributions
        if config.plot_counts: plot_distributions(counts_distributions, logscale=True, panels=True, frequencies=True)

    # -----------------------------------------------------------------

    def show_parameters(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Get best simulation name and chi squared
        best_simulation_name, best_chi_squared = chi_squared.best_simulation_name_and_chi_squared

        # Get best simulation parameter values
        best_parameter_values = parameters.parameter_values_for_simulation(best_simulation_name)

        # Most probable model: should be same as simulation with lowest chi squared
        most_probable_simulation_name = generation.most_probable_model
        # print(most_probable_simulation_name)
        assert most_probable_simulation_name == best_simulation_name

        print("")
        print("Statistics:")

        # Loop over the free parameter labels
        for label in parameter_labels:

            print("")
            print(" - " + fmt.bold + label + fmt.reset + ":")
            print("")

            # Get most probable parameter value
            most_probable_value = generation.get_most_probable_parameter_value(label)

            print("    * Initial guess value: " + tostr(initial_parameter_values[label], ndigits=3))
            print("    * Best simulation value: " + tostr(best_parameter_values[label], ndigits=3))
            print("    * Most probable value: " + tostr(most_probable_value, ndigits=3) + " " + tostr(
                parameter_units[label]))
            if config.nsimulations > 1: print(
                "    * Most counted in " + str(config.nsimulations) + " best simulations: " + tostr(
                    counts_distributions[label].most_frequent, ndigits=3) + " " + tostr(parameter_units[label]))
        print("")

    # -----------------------------------------------------------------

    def plot_best(self, generation_name):

        """
        This function ...
        :return:
        """

        # Create the SED plotter
        plotter = SEDPlotter()

        # Inform the user
        log.info("Adding the observed SEDs ...")

        # Add observed SEDs
        plotter.add_sed(environment.observed_sed, "Clipped observed fluxes")
        plotter.add_sed(environment.truncated_sed, "Truncated observed fluxes")

        # Loop over the simulation names
        nseds = config.nsimulations
        for index, simulation_name in enumerate(simulation_names):
            # Debugging
            log.debug("Adding SED of the '" + simulation_name + "' simulation (" + str(index + 1) + " of " + str(
                nseds) + ") ...")

            # Get the simulated SED
            sed = generation.get_simulation_sed(simulation_name)

            # Add the SED
            plotter.add_sed(sed, simulation_name, ghost=True, residuals=False)

        # Run the plotter
        plotter.run(output=config.seds_output)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_terms_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("simulation", "string", "simulation name")
        definition.add_optional("path", "string", "write the plot file")
        return definition

    # -----------------------------------------------------------------

    def plot_terms_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get generation name
        generation_name, config = self.get_generation_name_and_config_from_command(command, self.plot_terms_definition, **kwargs)

        # Plot
        self.plot_terms(generation_name, config.simulation)

    # -----------------------------------------------------------------

    def plot_terms(self, generation_name, simulation_name, path=None):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :param path:
        :return:
        """

        # Get the distribution
        distribution = self.get_chi_squared_terms_distribution(generation_name, simulation_name)

        # Plot
        plot_distribution(distribution, logscale=True, statistics=False, color="xkcd:orange", path=path)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_ranks_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("chisquared_or_prob", "string", "plot chi squared values or probabilities", choices=["chisquared", "prob"])
        definition.add_optional("path", "string", "write the plot file")
        return definition

    # -----------------------------------------------------------------

    def plot_ranks_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get generation name
        generation_name, config = self.get_generation_name_and_config_from_command(command, self.plot_ranks_definition, **kwargs)

        # Chi squared values
        if config.chisquared_or_prob == "chisquared": self.plot_ranks_chisquared(generation_name)

        # Probabilities
        elif config.chisquared_or_prob == "prob": self.plot_ranks_probabilities(generation_name)

    # -----------------------------------------------------------------

    def plot_ranks_chisquared(self, generation_name, path=None):

        """
        This function ...
        :param generation_name:
        :param path:
        :return:
        """

        # Get the rank distribution
        distribution = self.get_chi_squared_ranks_distribution(generation_name)

        # Plot histogram of the chi squared values w.r.t. simulation rank
        plot_distribution(distribution, statistics=False, path=path)

    # -----------------------------------------------------------------

    def plot_ranks_probabilities(self, generation_name, path=None):

        """
        This function ...
        :param generation_name:
        :param path:
        :return:
        """

        # Get the rank distribution
        distribution = self.get_probability_ranks_distribution(generation_name)

        # Plot histogram of the probability values w.r.t. simulation rank
        plot_distribution(distribution, statistics=False, color="red", logfrequency=True, path=path)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_chi_squared_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_optional("nbins", "positive_integer", "number of bins", 50)
        definition.add_optional("max_chisquared", "positive_real", "maximum chi squared to show")
        definition.add_optional("path", "string", "save the plot file")
        return definition

    # -----------------------------------------------------------------

    def plot_chi_squared_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get generation name
        generation_name, config = self.get_generation_name_and_config_from_command(command, self.plot_chi_squared_definition, **kwargs)

        # Plot
        self.plot_chi_squared(generation_name, path=config.path, nbins=config.nbins, max_chisquared=config.max_chisquared)

    # -----------------------------------------------------------------

    def plot_chi_squared(self, generation_name, path=None, nbins=50, max_chisquared=None):

        """
        This function ...
        :param generation_name:
        :param path:
        :param nbins:
        :param max_chisquared:
        :return:
        """

        # Make distribution of chi squared values and plot histogram
        distribution = self.get_chi_squared_distribution(generation_name, nbins=nbins, clip_above=max_chisquared)

        # Plot
        plot_distribution(distribution, path=path)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_probabilities_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_optional("nbins", "positive_integer", "number of bins", 50)
        definition.add_optional("min_probability", "positive_real", "minimum probability to show")
        definition.add_optional("path", "string", "save the plot file")
        return definition

    # -----------------------------------------------------------------

    def plot_probabilities_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get generation name
        generation_name, config = self.get_generation_name_and_config_from_command(command, self.plot_probabilities_definition, **kwargs)

        # Plot
        self.plot_probabilities(generation_name, path=config.path, nbins=config.nbins, min_probability=config.min_probability)

    # -----------------------------------------------------------------

    def plot_probabilities(self, generation_name, path=None, nbins=50, min_probability=None):

        """
        This function ...
        :param generation_name:
        :param path:
        :param nbins:
        :param min_probability:
        :return:
        """

        # Make distribution of probability values and plot histogram
        distribution = self.get_probability_distribution(generation_name, nbins=nbins, clip_below=min_probability)

        # Plot the distribution
        plot_distribution(distribution, color="red", path=path)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_seds_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_optional("path", "string", "save the plot file")
        definition.add_optional("additional_error", "percentage", "additional percentual error for the observed flux points")
        definition.add_optional("random", "positive_integer", "pick a specified number of random simulations to plot")
        return definition

    # -----------------------------------------------------------------

    def plot_seds_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the generation name
        generation_name, config = self.get_generation_name_and_config_from_command(command, self.plot_seds_definition, **kwargs)

        # Plot
        self.plot_seds(generation_name, path=config.path, additional_error=config.additional_error, random=config.random)

    # -----------------------------------------------------------------

    def plot_seds(self, generation_name, path=None, additional_error=None, random=None):

        """
        This function ...
        :param generation_name:
        :param path:
        :param additional_error:
        :param random:
        :return:
        """

        # Debugging
        log.debug("Plotting the model SEDs of generation '" + generation_name + "' ...")

        # Create the SED plotter
        plotter = SEDPlotter()

        # Inform the user
        log.info("Loading the observed SEDs ...")

        # Add relative error
        if additional_error is not None:
            clipped_sed = self.clipped_sed.copy()
            truncated_sed = self.truncated_sed.copy()
            clipped_sed.add_relative_error(additional_error)
            truncated_sed.add_relative_error(additional_error)
        else: clipped_sed, truncated_sed = self.clipped_sed, self.truncated_sed

        # Inform the user
        log.info("Adding the observed SEDs ...")

        # Add observed SEDs
        plotter.add_sed(clipped_sed, "Clipped observed fluxes")
        plotter.add_sed(truncated_sed, "Truncated observed fluxes")

        # Add simulation SEDs
        log.info("Adding the simulated SEDs ...")

        # Get the simulation names
        names = self.get_simulation_names(generation_name)
        nsimulations = len(names)

        # Set number of SEDs
        if random is not None: nseds = random
        else: nseds = nsimulations

        # Get selected simulation names
        if random is not None: names = sequences.random_subset(names, random, avoid_duplication=True)

        # Loop over the simulation names
        for index, simulation_name in enumerate(names):

            # Debugging
            log.debug("Adding SED of the '" + simulation_name + "' simulation (" + str(index + 1) + " of " + str(nseds) + ") ...")

            # Get the simulated SED
            sed = self.get_sed(generation_name, simulation_name)

            # Add the SED
            plotter.add_sed(sed, simulation_name, ghost=True, residuals=False)

        # Run the plotter
        plotter.run(output=path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Show status
        self.show_status()

    # -----------------------------------------------------------------

    @property
    def class_pts_user_path(self):

        """
        This function ...
        :return:
        """

        #return introspection.pts_user_manager_dir
        return None

# -----------------------------------------------------------------
