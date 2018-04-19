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
from ...core.tools.utils import lazyproperty, memoize_method
from .component import AnalysisComponent
from ...core.basics.configurable import InteractiveConfigurable
from ...core.basics.log import log
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr
from ...core.basics.configuration import ConfigurationDefinition
from ...core.plot.wavelengthgrid import plot_wavelength_grid
#from .run import contributions
from ...core.plot.grids import make_grid_plot
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ..core.environment import load_modeling_environment
from ...core.plot.sed import plot_seds, SEDPlotter

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

#_sfr_stellar_name = "sfr_stellar"
#_sfr_dust_name = "sfr_dust"

_stellar_name = "stellar"
_dust_name = "dust"

_contributions_name = "contributions"
_components_name = "components"

# -----------------------------------------------------------------

# SED subcommands
sed_commands = OrderedDict()

## TOTAL
sed_commands[_total_name] = ("plot_total_sed_command", True, "plot the SED of the total simulation", None)
sed_commands[_stellar_name] = ("plot_stellar_sed_command", True, "plot the stellar SED(s)", None)
sed_commands[_dust_name] = ("plot_dust_sed_command", True, "plot the dust SED(s)", None)

## CONTRIBUTIONS
sed_commands[_contributions_name] = ("plot_contribution_seds_command", True, "plot the contributions to the total SED(s)")
sed_commands[_components_name] = ("plot_component_seds_command", True, "plot the SED(s) for different components")

## OLD BULGE
sed_commands[_old_bulge_name] = ("plot_old_bulge_sed_command", True, "plot the SED of the old stellar bulge", None)

## OLD DISK
sed_commands[_old_disk_name] = ("plot_old_disk_sed_command", True, "plot the SED of the old stellar disk", None)

## OLD
sed_commands[_old_name] = ("plot_old_sed_command", True, "plot the SED of the old stars", None)

## YOUNG
sed_commands[_young_name] = ("plot_young_sed_command", True, "plot the SED of the young stars", None)

## SFR
sed_commands[_sfr_name] = ("plot_sfr_sed_command", True, "plot the SED of the star formation regions", None)
#sed_commands[_sfr_stellar_name] = ("plot_sfr_stellar_sed_command", True, "plot the stellar SED of the star formation regions", None)
#sed_commands[_sfr_dust_name] = ("plot_sfr_dust_sed_command", True, "plot the dust SED of the star formation regions", None)

## UNEVOLVED
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
default_orientations = (earth_name,)

# -----------------------------------------------------------------

observed_name = "observed"
intrinsic_name = "intrinsic"
default_observed_intrinsic = (observed_name, intrinsic_name)
observed_intrinsic_choices = default_observed_intrinsic

# -----------------------------------------------------------------

#default_contributions = ("total",)

# -----------------------------------------------------------------

grid_orientations = ["xy", "xz", "yz", "xyz"]

# -----------------------------------------------------------------

clipped_name = "clipped"
truncated_name = "truncated"
asymptotic_sed = "asymptotic"

default_sed_references = (clipped_name, truncated_name)
sed_reference_descriptions = dict()
sed_reference_descriptions[clipped_name] = "Observed clipped fluxes"
sed_reference_descriptions[truncated_name] = "Observed truncated fluxes"
sed_reference_descriptions[asymptotic_sed] = "Observed asymptotic fluxes"

# -----------------------------------------------------------------

bulge = "bulge"
disk = "disk"
old = "old"
young = "young"
sfr = "sfr"
unevolved = "unevolved"
total = "total"

# Make lists
components = [bulge, disk, old, young, sfr, unevolved, total]
default_components = [total, old, young, sfr]

# -----------------------------------------------------------------

from ..core.model import contributions, total_contribution, direct_contribution, scattered_contribution, dust_contribution, transparent_contribution
from ..core.model import dust_direct_contribution, dust_scattered_contribution
#contributions = []
default_contributions = [total_contribution, direct_contribution, scattered_contribution, dust_contribution, transparent_contribution]

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
    def derived_parameter_values_total(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_total

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_bulge(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_bulge

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_disk(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_disk

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_old(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_old

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_young(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_young

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_sfr(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_sfr

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_dust(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_dust

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

        # Derived parameter values of total model
        print(fmt.cyan + fmt.underlined + "Derived parameter values of total model:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_total: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_total[label]))
        print("")

        # Derived parameter values of bulge
        print(fmt.cyan + fmt.underlined + "Derived parameter values of bulge:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_bulge: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_bulge[label]))
        print("")

        # Derived parameter values of disk
        print(fmt.cyan + fmt.underlined + "Derived parameter values of disk:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_disk: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_disk[label]))
        print("")

        # Derived parameter values of old component
        print(fmt.cyan + fmt.underlined + "Derived parameter values of old stellar component:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_old: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_old[label]))
        print("")

        # Derived parameter values of young component
        print(fmt.cyan + fmt.underlined + "Derived parameter values of young stellar component:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_young: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_young[label]))
        print("")

        # Derived parameter values of SF component
        print(fmt.cyan + fmt.underlined + "Derived parameter values of SFR component:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_sfr: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_sfr[label]))
        print("")

        # Derived parameter values of dust component
        print(fmt.cyan + fmt.underlined + "Derived parameter values of dust component:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_dust: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_dust[label]))
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

        # Get the config
        config = self.get_config_from_command(command, self.plot_wavelengths_definition)

        # Plot the wavelengths
        self.plot_wavelengths(**config)

    # -----------------------------------------------------------------

    def plot_wavelengths(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Plot the wavelength grid
        plot_wavelength_grid(self.wavelength_grid, "wavelengths", **kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_grid_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("orientation", "string", "plotting viewpoint", choices=grid_orientations)
        definition.add_optional("path", "string", "path for the plot file")
        return definition

    # -----------------------------------------------------------------

    def plot_grid_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.plot_grid_definition)

        # Get the path
        orientation = config.pop("orientation")
        path = config.pop("path")

        # Plot
        self.plot_grid(orientation, path=path)

    # -----------------------------------------------------------------

    def plot_grid(self, orientation, path=None, show=None):

        """
        This function ...
        :param orientation:
        :param path:
        :param show:
        :return:
        """

        # Debugging
        log.debug("Plotting the dust grid from orientation '" + orientation + "' ...")

        # Determine filepath
        if path is None:
            show = True
            path = fs.join(introspection.pts_temp_dir, "grid_" + orientation + ".pdf")

        # Determine grid filepath
        if orientation == "xy": grid_path = self.model.grid_xy_filepath
        elif orientation == "xz": grid_path = self.model.grid_xz_filepath
        elif orientation == "yz": grid_path = self.model.grid_yz_filepath
        elif orientation == "xyz": grid_path = self.model.grid_xyz_filepath
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Plot
        make_grid_plot(grid_path, path)

        # Open the plot?
        if show: fs.open_file(path)

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

    @property
    def asymptotic_sed(self):

        """
        This function ...
        :return:
        """

        return self.modeling_environment.asymptotic_sed

    # -----------------------------------------------------------------

    @memoize_method
    def get_reference_seds(self, additional_error=None, references=default_sed_references):

        """
        This function ...
        :param additional_error:
        :param references:
        :return:
        """

        # Debugging
        log.debug("Loading the observed SEDs ...")

        # Create dictionary
        seds = OrderedDict()

        # Add observed SEDs
        for label in references:

            # Get sed
            if label == clipped_name: sed = self.clipped_sed
            elif label == truncated_name: sed = self.truncated_sed
            elif label == asymptotic_sed: sed = self.asymptotic_sed
            else: raise ValueError("Invalid reference SED name")

            # Add relative error?
            if additional_error is not None:
                sed = sed.copy()
                sed.add_or_set_relative_error(additional_error)

            # Add
            description = sed_reference_descriptions[label]
            seds[description] = sed

        # Return the seds
        return seds

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_total_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("orientations", "string_list", "instrument orientation", default_orientations, choices=orientations)
        definition.add_flag("add_references", "add reference SEDs", False)
        definition.add_optional("additional_error", "percentage", "additional percentual error for the observed flux points")
        return definition

    # -----------------------------------------------------------------

    def plot_total_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.plot_total_sed_definition, **kwargs)

        # Plot
        self.plot_total_sed(orientations=config.orientations, add_references=config.add_references, additional_error=config.additional_error)

    # -----------------------------------------------------------------

    def plot_total_sed(self, orientations=default_orientations, add_references=False, additional_error=None, path=None,
                       show_file=False, title=None, format=None):

        """
        This function ...
        :param orientations:
        :param add_references:
        :param additional_error:
        :param path:
        :param show_file:
        :param title:
        :param format:
        :return:
        """

        # Debugging
        log.debug("Plotting total SED(s) ...")

        # Add SEDs
        #seds = OrderedDict()
        #if add_references: seds.update(self.get_reference_seds(additional_error=additional_error))
        #for orientation in orientations:
        #    if orientation == earth_name: seds[orientation] = self.model.observed_total_sed
        #    elif orientation == faceon_name: seds[orientation] = self.model.faceon_observed_total_sed
        #    elif orientation == edgeon_name: seds[orientation] = self.model.edgeon_observed_total_sed
        #    else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Plot
        #plot_seds(seds, models_residuals=True)

        # Create SED plotter
        #plotter = SEDPlotter(kwargs)
        plotter = SEDPlotter()

        # Add SEDs
        #for name in seds:
        #    sed = seds[name]
        #    plotter.add_sed(sed, label=name)

        # Add references?
        if add_references: plotter.add_seds(self.get_reference_seds(additional_error=additional_error))

        # Add orientations
        for orientation in orientations:

            if orientation == earth_name: plotter.add_sed(self.model.observed_total_sed, earth_name)
            elif orientation == faceon_name: plotter.add_sed(self.model.faceon_observed_total_sed, faceon_name, residuals=False)
            elif orientation == edgeon_name: plotter.add_sed(self.model.edgeon_observed_total_sed, edgeon_name, residuals=False)
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Set filepath, if plot is to be shown as file
        if path is None and show_file:
            if format is None: raise ValueError("Format has to be specified")
            path = fs.join(introspection.pts_temp_dir, "total_seds." + format)

        # Run the plotter
        plotter.run(title=title, output=path)

        # Show file
        if show_file: fs.open_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_stellar_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_stellar_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
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
        return definition

    # -----------------------------------------------------------------

    def plot_dust_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_contribution_seds_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("contributions", "string_list", "contributions", default_contributions, choices=contributions)
        definition.add_optional("component", "string", "component", total, choices=components)
        return definition

    # -----------------------------------------------------------------

    def plot_contribution_seds_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_contribution_seds_definition, **kwargs)

        # TODO: use 'component' option to do this for different components?

        # Plot
        self.plot_contribution_seds(config.contributions)

    # -----------------------------------------------------------------

    def plot_contribution_seds(self, contributions, path=None, title=None, show_file=False, format=None):

        """
        This function ...
        :param contributions:
        :param path:
        :param title:
        :param show_file:
        :param format:
        :return:
        """

        # Debugging
        log.debug("Plotting contribution SEDs ...")

        # Create SED plotter
        plotter = SEDPlotter()

        # Loop over the contributions
        for contribution in contributions:

            # TODO: different components (simulations)?
            if contribution == total_contribution: plotter.add_sed(self.model.observed_total_sed, contribution)
            elif contribution == direct_contribution: plotter.add_sed(self.model.observed_total_sed_direct, contribution, residuals=False)
            elif contribution == scattered_contribution: plotter.add_sed(self.model.observed_total_sed_scattered, contribution, residuals=False)
            elif contribution == dust_contribution: plotter.add_sed(self.model.observed_total_sed_dust, contribution, residuals=False)
            elif contribution == dust_direct_contribution: plotter.add_sed(self.model.observed_total_sed_dust_direct, contribution, residuals=False)
            elif contribution == dust_scattered_contribution: plotter.add_sed(self.model.observed_total_sed_dust_scattered, contribution, residuals=False)
            elif contribution == transparent_contribution: plotter.add_sed(self.model.observed_total_sed_transparent, contribution, residuals=False)
            else: raise ValueError("Invalid contribution: '" + contribution + "'")

        # Set filepath, if plot is to be shown as file
        if path is None and show_file:
            if format is None: raise ValueError("Format has to be specified")
            path = fs.join(introspection.pts_temp_dir, "contribution_seds." + format)

        # Run the plotter
        plotter.run(title=title, output=path)

        # Show file
        if show_file: fs.open_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_component_seds_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("components", "string_list", "components", default_components, choices=components)
        return definition

    # -----------------------------------------------------------------

    def plot_component_seds_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_component_seds_definition, **kwargs)

        # Plot
        self.plot_component_seds(config.components)

    # -----------------------------------------------------------------

    def plot_component_seds(self, components, path=None, title=None, show_file=False, format=None):

        """
        This function ...
        :param components:
        :param path:
        :param title:
        :param show_file:
        :param format:
        :return:
        """

        # Debugging
        log.debug("Plotting component SEDs ...")

        # Create SED plotter
        plotter = SEDPlotter()

        # Add references?
        #if add_references: plotter.add_seds(self.get_reference_seds(additional_error=additional_error))

        # Loop over the components
        for component in components:

            if component == total: plotter.add_sed(self.model.observed_total_sed, total)
            elif component == bulge: plotter.add_sed(self.model.observed_old_bulge_sed, bulge, residuals=False)
            elif component == disk: plotter.add_sed(self.model.observed_old_disk_sed, disk, residuals=False)
            elif component == old: plotter.add_sed(self.model.observed_old_sed, old, residuals=False)
            elif component == young: plotter.add_sed(self.model.observed_young_sed, young, residuals=False)
            elif component == sfr: plotter.add_sed(self.model.observed_sfr_sed, sfr, residuals=False)
            elif component == unevolved: plotter.add_sed(self.model.observed_unevolved_sed, unevolved, residuals=False)
            else: raise ValueError("Invalid component: '" + component + "'")

        # Set filepath, if plot is to be shown as file
        if path is None and show_file:
            if format is None: raise ValueError("Format has to be specified")
            path = fs.join(introspection.pts_temp_dir, "component_seds." + format)

        # Run the plotter
        plotter.run(title=title, output=path)

        # Show file
        if show_file: fs.open_file(path)

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
        definition.add_positional_optional("components", "string_list", "components", default_components, choices=components)
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

    def plot_dust_sed(self, components=default_components):

        """
        This function ...
        :param components:
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
