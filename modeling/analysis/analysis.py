#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.analysis Contains the Analysis class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools.utils import lazyproperty
from .component import AnalysisComponent
from ...core.basics.configurable import InteractiveConfigurable
from ...core.basics.log import log
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr
from ...core.basics.configuration import ConfigurationDefinition
from ...core.plot.wavelengthgrid import plot_wavelength_grid
from .run import contributions

# -----------------------------------------------------------------

# Standard commands
_help_command_name = "help"
_history_command_name = "history"
_status_command_name = "status"

# Other commands
_model_command_name = "model"

# Plot commands
_wavelengths_command_name = "wavelengths"
_dustgrid_command_name = "grid"
_sed_command_name = "sed"

# -----------------------------------------------------------------

# Define commands
commands = OrderedDict()

# Standard commands
commands[_help_command_name] = ("show_help", False, "show help", None)
commands[_history_command_name] = ("show_history_command", True, "show history of executed commands", None)
commands[_status_command_name] = ("show_status_command", True, "show analysis status", None)

# Show stuff
commands[_model_command_name] = ("show_model", False, "show the model properties", None)

# Plot stuff
commands[_wavelengths_command_name] = ("plot_wavelengths_command", True, "plot the wavelength grid", None)
commands[_dustgrid_command_name] = ("plot_grid_command", True, "plot the dust grid", None)
commands[_sed_command_name] = (None, None, "plot SEDs", None)

# -----------------------------------------------------------------

_total_name = "total"
_old_bulge_name = "old_bulge"
_old_disk_name = "old_disk"
_old_name = "old"
_young_name = "young"
_sfr_name = "sfr"
_unevolved_name = "unevolved"

_sfr_stellar_name = "sfr_stellar"
_sfr_dust_name = "sfr_dust"
_dust_name = "dust"

# -----------------------------------------------------------------

# SED subcommands
sed_commands = OrderedDict()
sed_commands[_total_name] = ("plot_total_sed_command", True, "plot the SED of the total simulation", None)
sed_commands[_old_bulge_name] = ("plot_old_bulge_sed_command", True, "plot the SED of the old stellar bulge", None)
sed_commands[_old_disk_name] = ("plot_old_disk_sed_command", True, "plot the SED of the old stellar disk", None)
sed_commands[_old_name] = ("plot_old_sed_command", True, "plot the SED of the old stars", None)
sed_commands[_young_name] = ("plot_young_sed_command", True, "plot the SED of the young stars", None)
sed_commands[_sfr_name] = ("plot_sfr_sed_command", True, "plot the SED of the star formation regions", None)
sed_commands[_sfr_stellar_name] = ("plot_sfr_stellar_sed_command", True, "plot the stellar SED of the star formation regions", None)
sed_commands[_sfr_dust_name] = ("plot_sfr_dust_sed_command", True, "plot the dust SED of the star formation regions", None)
sed_commands[_unevolved_name] = ("plot_unevolved_sed_command", True, "plot the SED of the unevolved stellar population (young + sfr)", None)
sed_commands[_dust_name] = ("plot_dust_sed_command", True, "plot the dust emission SED", None)

# -----------------------------------------------------------------

# Set subcommands
subcommands = OrderedDict()
subcommands[_sed_command_name] = sed_commands

# -----------------------------------------------------------------

earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"
orientations = (earth_name, faceon_name, edgeon_name)

# -----------------------------------------------------------------

observed_name = "observed"
intrinsic_name = "intrinsic"
default_observed_intrinsic = (observed_name, intrinsic_name)
observed_intrinsic_choices = default_observed_intrinsic

# -----------------------------------------------------------------

default_contributions = ("total",)

# -----------------------------------------------------------------

class Analysis(AnalysisComponent, InteractiveConfigurable):

    """
    This class ...
    """

    _commands = commands
    _subcommands = subcommands

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        InteractiveConfigurable.__init__(self, no_config=True)
        AnalysisComponent.__init__(self, *args, **kwargs)

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
        Thisf unction ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Run commands
        if self.do_commands: self.run_commands()

        # 3. Interactive
        if self.do_interactive: self.interactive()

        # 4. Show
        self.show()

        # 5. Write the history
        if self.has_commands: self.write_history()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Analysis, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def analysis_run(self):

        """
        This function ...
        :return:
        """

        return self.get_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.model_name

    # -----------------------------------------------------------------

    @property
    def model(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.model

    # -----------------------------------------------------------------

    @property
    def parameter_values(self):

        """
        This function ...
        :return:
        """

        return self.model.parameter_values

    # -----------------------------------------------------------------

    @property
    def free_parameter_values(self):

        """
        This function ...
        :return:
        """

        return self.model.free_parameter_values

    # -----------------------------------------------------------------

    @property
    def other_parameter_values(self):

        """
        This function ...
        :return:
        """

        return self.model.other_parameter_values

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values

    # -----------------------------------------------------------------

    @property
    def generation_name(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.generation_name

    # -----------------------------------------------------------------

    @property
    def simulation_name(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.simulation_name

    # -----------------------------------------------------------------

    @property
    def chi_squared(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.chi_squared

    # -----------------------------------------------------------------

    @property
    def fitting_run_name(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.fitting_run_name

    # -----------------------------------------------------------------

    @property
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.fitting_run

    # -----------------------------------------------------------------

    @property
    def from_fitting(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.from_fitting

    # -----------------------------------------------------------------

    @property
    def wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.wavelength_grid

    # -----------------------------------------------------------------

    @property
    def dust_grid(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.dust_grid

    # -----------------------------------------------------------------

    @lazyproperty
    def show_status_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Return the definition
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
        self.show_status()

    # -----------------------------------------------------------------

    def show_status(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing status ...")

    # -----------------------------------------------------------------

    def show_model(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Showing the model properties ...")

        # Show the model name
        print("")
        print(fmt.yellow + fmt.underlined + "Model name" + fmt.reset + ": " + self.model_name)
        if self.generation_name is not None: print(fmt.yellow + fmt.underlined + "Generation name" + fmt.reset + ": " + self.generation_name)
        if self.simulation_name is not None: print(fmt.yellow + fmt.underlined + "Simulation name" + fmt.reset + ": " + self.simulation_name)
        if self.chi_squared is not None: print(fmt.yellow + fmt.underlined + "Chi-squared" + fmt.reset + ": " + tostr(self.chi_squared))
        print("")

        # Show the free parameter values
        print(fmt.cyan + fmt.underlined + "Free parameter values:" + fmt.reset)
        print("")
        for label in self.free_parameter_values: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.free_parameter_values[label]))
        print("")

        # Show the other parameter values
        print(fmt.cyan + fmt.underlined + "Other parameter values:" + fmt.reset)
        print("")
        for label in self.other_parameter_values: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.other_parameter_values[label]))
        print("")

        # Derived parameter values
        print(fmt.cyan + fmt.underlined + "Derived parameter values:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values[label]))
        print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_wavelengths_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_wavelengths_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    def plot_wavelengths(self):

        """
        This function ...
        :return:
        """

        # Plot the wavelength grid
        plot_wavelength_grid(self.wavelength_grid, "wavelengths")

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_grid_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_grid_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    def plot_grid(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_total_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("orientation", "string", "instrument orientation", earth_name, orientations)
        return definition

    # -----------------------------------------------------------------

    def plot_total_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_old_bulge_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_old_bulge_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_old_disk_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_old_disk_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_old_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED, or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_old_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_young_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED, or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_young_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_sfr_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED, or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_sfr_sed_command(self, command, **kwargs):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_sfr_stellar_sed_definition(self):

        """
        Thisf unction ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED, or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_sfr_stellar_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_sfr_dust_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED, or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_sfr_dust_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_unevolved_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED, or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_unevolved_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Parse the command
        config = self.get_config_from_command(command, self.plot_unevolved_sed_definition, **kwargs)

        # Plot
        self.plot_unevolved_sed(**config)

    # -----------------------------------------------------------------

    def plot_unevolved_sed(self, observed_intrinsic=observed_intrinsic_choices):

        """
        This function ...
        :param observed_intrinsic:
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_dust_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("contributions", "string_list", "stellar contributions", default_contributions, choices=contributions)
        return definition

    # -----------------------------------------------------------------

    def plot_dust_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the configuration
        config = self.get_config_from_command(command, **kwargs)
        config.pop("_path")

        # Plot
        self.plot_dust_sed(**config)

    # -----------------------------------------------------------------

    def plot_dust_sed(self, contributions=default_contributions):

        """
        This function ...
        :param contributions:
        :return:
        """



    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    @property
    def history_filename(self):

        """
        This function ...
        :return:
        """

        return "analysis"

# -----------------------------------------------------------------
