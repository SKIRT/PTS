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
from .component import AnalysisRunComponent
from ...core.basics.configurable import InteractiveConfigurable
from ...core.basics.log import log
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr
from ...core.basics.configuration import ConfigurationDefinition
from ...core.plot.wavelengthgrid import plot_wavelength_grid
from ...core.plot.grids import make_grid_plot
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ..core.environment import load_modeling_environment
from ...core.plot.sed import plot_seds, SEDPlotter, plot_sed
from ...core.config.plot_seds import definition as plot_seds_definition
from ..config.evaluate_analysis import definition as evaluate_analysis_definition
from ...core.plot.attenuation import plot_attenuation_curve, plot_attenuation_curves
from ..config.analyse_cell_heating import definition as analyse_cell_heating_definition
from ..config.analyse_projected_heating import definition as analyse_projected_heating_definition
from .heating.cell import CellDustHeatingAnalyser
from .heating.projected import ProjectedDustHeatingAnalyser
from ..config.analyse_properties import definition as analyse_properties_definition
from .properties import PropertiesAnalyser
from ..config.analyse_cell_energy import definition as analyse_cell_energy_definition
from ..config.analyse_projected_energy import definition as analyse_projected_energy_definition
from .energy.cell import CellEnergyAnalyser
from .energy.projected import ProjectedEnergyAnalyser
from ...magic.tools.plotting import plot_frame, plot_frame_contours
from ...core.filter.filter import Filter, parse_filter
from ...core.tools import types
from ...magic.plot.imagegrid import StandardImageGridPlotter, ResidualImageGridPlotter
from .evaluation import AnalysisModelEvaluator
from ...core.tools import sequences

from .properties import bol_map_name, intr_stellar_map_name, obs_stellar_map_name, diffuse_dust_map_name, dust_map_name
from .properties import scattered_map_name, absorbed_diffuse_map_name, fabs_diffuse_map_name, fabs_map_name, stellar_mass_map_name, ssfr_map_name
from .properties import attenuated_map_name, direct_map_name, sfr_map_name, i1_map_name, intr_i1_map_name, fuv_map_name
from .properties import intr_fuv_map_name, dust_mass_map_name, stellar_lum_map_name, intr_dust_map_name
from .properties import diffuse_mass_map_name, mass_map_name, earth_name, faceon_name, edgeon_name

# -----------------------------------------------------------------

# Define names of maps to show
#total_map_names = (bol_map_name, intr_stellar_map_name, obs_stellar_map_name, dust_map_name, dust_with_internal_map_name, scattered_map_name, absorbed_map_name, absorbed_with_internal_map_name, attenuated_map_name, direct_map_name,)
total_map_names = (bol_map_name, intr_stellar_map_name, obs_stellar_map_name, diffuse_dust_map_name, dust_map_name, scattered_map_name, absorbed_diffuse_map_name, fabs_diffuse_map_name, fabs_map_name, attenuated_map_name, direct_map_name, sfr_map_name, stellar_mass_map_name,ssfr_map_name,)
bulge_map_names = (bol_map_name, direct_map_name, i1_map_name, intr_i1_map_name, dust_map_name,)
disk_map_names = (bol_map_name, direct_map_name, i1_map_name, intr_i1_map_name, dust_map_name,)
old_map_names = (bol_map_name, direct_map_name, i1_map_name, intr_i1_map_name, dust_map_name,)
young_map_names = (bol_map_name, direct_map_name, fuv_map_name, intr_fuv_map_name, dust_map_name,)
sfr_map_names = (bol_map_name, direct_map_name, fuv_map_name, intr_fuv_map_name, sfr_map_name, dust_mass_map_name, stellar_lum_map_name, intr_dust_map_name, dust_map_name)
unevolved_map_names = (bol_map_name, direct_map_name, fuv_map_name, intr_fuv_map_name, sfr_map_name, dust_map_name,)
#dust_map_names = (mass_map_name, total_mass_map_name,) #lum_map_name, total_lum_map_name,)
dust_map_names = (diffuse_mass_map_name, mass_map_name,)

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
_attenuation_command_name = "attenuation"
_map_command_name = "map"
_images_command_name = "images"
_cubes_command_name = "cubes"

# Evaluate
_evaluate_command_name = "evaluate"

# Analysis
_properties_command_name = "properties"
_heating_command_name = "heating"
_energy_command_name = "energy"

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
commands[_attenuation_command_name] = (None, None, "plot attenuation curves", None)
commands[_map_command_name] = (None, None, "plot a map", None)
commands[_images_command_name] = ("plot_images_command", True, "plot the simulated images", None)
commands[_cubes_command_name] = ("plot_cubes_command", True, "plot the simulated datacubes", None)

# Evaluate
commands[_evaluate_command_name] = ("evaluate_command", True, "evaluate the analysis model", None)

# Analysis
commands[_properties_command_name] = ("analyse_properties_command", True, "analyse the model properties", None)
commands[_heating_command_name] = (None, None, "analyse dust heating contributions", None)
commands[_energy_command_name] = (None, None, "analyse the energy budget in the galaxy", None)

# -----------------------------------------------------------------

_bulge_name = "bulge"
_disk_name = "disk"

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
sed_commands[_contributions_name] = ("plot_contribution_seds_command", True, "plot the contributions to the total SED(s)", None)
sed_commands[_components_name] = ("plot_component_seds_command", True, "plot the SED(s) for different components", None)

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

#sed_commands[_sfr_stellar_name] = ("plot_sfr_stellar_sed_command", True, "plot the intrinsic/observed stellar SED of the star formation regions", None)
#sed_commands[_sfr_dust_name] = ("plot_sfr_dust_sed_command", True, "plot the intrinsic/observed dust SED of the star formation regions", None)

## UNEVOLVED
sed_commands[_unevolved_name] = ("plot_unevolved_sed_command", True, "plot the SED of the unevolved stellar population (young + sfr)", None)

# -----------------------------------------------------------------

# Attenuation subcommands
attenuation_commands = OrderedDict()

## TOTAL
attenuation_commands[_total_name] = ("plot_total_attenuation_command", True, "plot the attenuation curve of the model", None)

## CONTRIBUTIONS
attenuation_commands[_components_name] = ("plot_component_attenuation_command", True, "plot the attenuation curves of the different components", None)

## OLD BULGE
attenuation_commands[_old_bulge_name] = ("plot_old_bulge_attenuation_command", True, "plot the attenuation curve of the old stellar bulge", None)

## OLD DISK
attenuation_commands[_old_disk_name] = ("plot_old_disk_attenuation_command", True, "plot the attenuation curve of the old stellar disk", None)

## OLD
attenuation_commands[_old_name] = ("plot_old_attenuation_command", True, "plot the attenuation curve of the old stars", None)

## YOUNG
attenuation_commands[_young_name] = ("plot_young_attenuation_command", True, "plot the attenuation curve of the young stars", None)

## SFR
attenuation_commands[_sfr_name] = ("plot_sfr_attenuation_command", True, "plot the attenuation curve of the star formation regions", None)
# BUT WHAT IS THE *INTRINSIC* SFR ATTENUATION CURVE? (by INTERNAL DUST)

## UNEVOLVED
attenuation_commands[_unevolved_name] = ("plot_unevolved_attenuation_command", True, "plot the attenuation curve of the unevolved stellar population (young + sfr)", None)

# -----------------------------------------------------------------

# Map subcommands
map_commands = OrderedDict()

## TOTAL
map_commands[_total_name] = ("show_total_map_command", True, "show a map of the total model", None)

## Bulge
map_commands[_bulge_name] = ("show_bulge_map_command", True, "show a map of the old stellar bulge component", None)

## Disk
map_commands[_disk_name] = ("show_disk_map_command", True, "show a map of the old stellar disk component", None)

## Old
map_commands[_old_name] = ("show_old_map_command", True, "show a map of the old stellar component", None)

## Young
map_commands[_young_name] = ("show_young_map_command", True, "show a map of the young stellar component", None)

## SFR
map_commands[_sfr_name] = ("show_sfr_map_command", True, "show a map of the SFR component", None)

## Unevolved
map_commands[_unevolved_name] = ("show_unevolved_map_command", True, "show a map of the unevolved stellar component", None)

## Dust
map_commands[_dust_name] = ("show_dust_map_command", True, "show a map of the dust component", None)

# -----------------------------------------------------------------

_cell_name = "cell"
_projected_name = "projected"

# -----------------------------------------------------------------

# Heating subcommands
heating_commands = OrderedDict()

# Cell and projected
heating_commands[_cell_name] = ("analyse_cell_heating_command", True, "analyse the cell heating", None)
heating_commands[_projected_name] = ("analyse_projected_heating_command", True, "analyse the projected heating", None)

# -----------------------------------------------------------------

# Energy subcommands
energy_commands = OrderedDict()

# Cell and projected
energy_commands[_cell_name] = ("analyse_cell_energy_command", True, "analyse the cell energy budget", None)
energy_commands[_projected_name] = ("analyse_projected_energy_command", True, "analyse the projected energy budget", None)

# -----------------------------------------------------------------

# Set subcommands
subcommands = OrderedDict()
subcommands[_sed_command_name] = sed_commands
subcommands[_attenuation_command_name] = attenuation_commands
subcommands[_map_command_name] = map_commands
subcommands[_heating_command_name] = heating_commands
subcommands[_energy_command_name] = energy_commands

# -----------------------------------------------------------------

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

class Analysis(AnalysisRunComponent, InteractiveConfigurable):

    """
    This class ...
    """

    _commands = commands
    _subcommands = subcommands
    _log_section = "ANALYSIS"

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        InteractiveConfigurable.__init__(self, no_config=True)
        AnalysisRunComponent.__init__(self, *args, **kwargs)

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

    def _run(self, **kwargs):

        """
        Thisf unction ...
        :param kwargs:
        :return:
        """

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
        #super(Analysis, self).setup(**kwargs)
        AnalysisRunComponent.setup(self, **kwargs)
        InteractiveConfigurable.setup(self, **kwargs)

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
    def derived_parameter_values_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_unevolved

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
        print(fmt.cyan + fmt.underlined + "Derived parameter values of old bulge stellar component:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_bulge: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_bulge[label]))
        print("")

        # Derived parameter values of disk
        print(fmt.cyan + fmt.underlined + "Derived parameter values of old disk stellar component:" + fmt.reset)
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

        # Derived parameter values of unevolved components
        print(fmt.cyan + fmt.underlined + "Derived parameter values of unevolved stars:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_unevolved: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_unevolved[label]))
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
        definition.add_required("observed_intrinsic", "string_tuple", "plot observed stellar SED, intrinsic stellar SED, or both", choices=observed_intrinsic_choices)
        definition.add_positional_optional("components", "string_list", "components", [total], choices=components)
        return definition

    # -----------------------------------------------------------------

    def plot_stellar_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_stellar_sed_definition, **kwargs)

        # Plot
        self.plot_stellar_sed(config.observed_intrinsic, components=config.components)

    # -----------------------------------------------------------------

    def plot_stellar_sed(self, observed_intrinsic, components, path=None, title=None, show_file=False, format=None):

        """
        This function ...
        :param observed_intrinsic:
        :param components:
        :param path:
        :param title:
        :param show_file:
        :param format:
        :return:
        """

        # Debugging
        log.debug("Plotting stellar SED(s) ...")

        # Create SED plotter
        plotter = SEDPlotter()

        # Add references?
        # if add_references: plotter.add_seds(self.get_reference_seds(additional_error=additional_error))

        # Either observed of intrinsic
        if len(observed_intrinsic) == 1:

            oi = observed_intrinsic[0]

            # Loop over the components
            for component in components:

                # Total simulation
                if component == total:

                    if oi == observed_name: plotter.add_sed(self.model.observed_stellar_sed, total)
                    elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_total_sed, total, residuals=False)
                    else: raise ValueError("")

                # Old bulge simulation
                elif component == bulge:

                    if oi == observed_name: plotter.add_sed(self.model.observed_old_bulge_stellar_sed, bulge, residuals=False)
                    elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_old_bulge, bulge, residuals=False)
                    else: raise ValueError("")

                # Old disk simulation
                elif component == disk:

                    if oi == observed_name: plotter.add_sed(self.model.observed_old_disk_stellar_sed, disk, residuals=False)
                    elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_old_disk, disk, residuals=False)
                    else: raise ValueError("")

                # Old simulation
                elif component == old:

                    if oi == observed_name: plotter.add_sed(self.model.observed_old_stellar_sed, old, residuals=False)
                    elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_old, old, residuals=False)
                    else: raise ValueError("")

                # Young simulation
                elif component == young:

                    if oi == observed_name: plotter.add_sed(self.model.observed_young_stellar_sed, young, residuals=False)
                    elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_young, young, residuals=False)
                    else: raise ValueError("")

                # SFR simulation
                elif component == sfr:

                    if oi == observed_name: plotter.add_sed(self.model.observed_sfr_stellar_sed, sfr, residuals=False)
                    elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_sfr, sfr, residuals=False)
                    else: raise ValueError("")

                # Unevolved simulation
                elif component == unevolved:

                    if oi == observed_name: plotter.add_sed(self.model.observed_unevolved_stellar_sed, unevolved, residuals=False)
                    elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_unevolved, unevolved, residuals=False)
                    else: raise ValueError("")

                # Invalid
                else: raise ValueError("")

        else:

            # ALLOW this for multiple components?

            # Loop over the components
            for component in components:

                # Total simulation
                if component == total:

                    for oi in observed_intrinsic:
                        if oi == observed_name: plotter.add_sed(self.model.observed_stellar_sed, total + " " + observed_name)
                        elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_total_sed, total + " " + intrinsic_name, residuals=False)
                        else: raise ValueError("")

                # Old bulge simulation
                elif component == bulge:

                    for oi in observed_intrinsic:
                        if oi == observed_name: plotter.add_sed(self.model.observed_old_bulge_stellar_sed, bulge + " " + observed_name, residuals=False)
                        elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_old_bulge, bulge + " " + intrinsic_name, residuals=False)
                        else: raise ValueError("")

                # Old disk simulation
                elif component == disk:

                    for oi in observed_intrinsic:
                        if oi == observed_name: plotter.add_sed(self.model.observed_old_disk_stellar_sed, disk + " " + observed_name, residuals=False)
                        elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_old_disk, disk + " " + intrinsic_name, residuals=False)
                        else: raise ValueError("")

                # Old simulation
                elif component == old:

                    for oi in observed_intrinsic:
                        if oi == observed_name: plotter.add_sed(self.model.observed_old_stellar_sed, old + " " + observed_name, residuals=False)
                        elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_old, old + " " + intrinsic_name, residuals=False)
                        else: raise ValueError("")

                elif component == young:

                    for oi in observed_intrinsic:
                        if oi == observed_name: plotter.add_sed(self.model.observed_young_stellar_sed, young + " " + observed_name, residuals=False)
                        elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_young, young + " " + intrinsic_name, residuals=False)
                        else: raise ValueError("")

                elif component == sfr:

                    for oi in observed_intrinsic:
                        if oi == observed_name: plotter.add_sed(self.model.observed_sfr_stellar_sed, sfr + " " + observed_name, residuals=False)
                        elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_sfr, sfr + " " + intrinsic_name, residuals=False)
                        else: raise ValueError("")

                elif component == unevolved:

                    for oi in observed_intrinsic:
                        if oi == observed_name: plotter.add_sed(self.model.observed_unevolved_stellar_sed, unevolved + " " + observed_name, residuals=False)
                        elif oi == intrinsic_name: plotter.add_sed(self.model.intrinsic_sed_unevolved, unevolved + " " + intrinsic_name, residuals=False)
                        else: raise ValueError("")

                else: raise ValueError("Invalid component: '" + component + "'")

        # Set filepath, if plot is to be shown as file
        if path is None and show_file:
            if format is None: raise ValueError("Format has to be specified")
            path = fs.join(introspection.pts_temp_dir, "dust_seds." + format)

        # Run the plotter
        plotter.run(title=title, output=path)

        # Show file
        if show_file: fs.open_file(path)

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

        # Get config
        config = self.get_config_from_command(command, self.plot_dust_sed_definition, **kwargs)

        # Plot
        self.plot_dust_sed(config.components)

    # -----------------------------------------------------------------

    def plot_dust_sed(self, components, title=None, path=None, show_file=False, format=None):

        """
        This function ...
        :param components:
        :param title:
        :param path:
        :param show_file:
        :param format:
        :return:
        """

        # Debugging
        log.debug("Plotting dust SED(s) ...")

        # Create SED plotter
        plotter = SEDPlotter()

        # Add references?
        # if add_references: plotter.add_seds(self.get_reference_seds(additional_error=additional_error))

        # Loop over the components
        for component in components:

            if component == total: plotter.add_sed(self.model.dust_sed, total)
            elif component == bulge: plotter.add_sed(self.model.observed_old_bulge_dust_sed, bulge, residuals=False)
            elif component == disk: plotter.add_sed(self.model.observed_old_disk_dust_sed, disk, residuals=False)
            elif component == old: plotter.add_sed(self.model.observed_old_dust_sed, old, residuals=False)
            elif component == young: plotter.add_sed(self.model.observed_young_dust_sed, young, residuals=False)
            elif component == sfr: plotter.add_sed(self.model.observed_sfr_dust_sed, sfr, residuals=False)
            elif component == unevolved: plotter.add_sed(self.model.observed_unevolved_dust_sed, unevolved, residuals=False)
            else: raise ValueError("Invalid component: '" + component + "'")

            # THERE IS ALSO model.diffuse_dust_sed !

        # Set filepath, if plot is to be shown as file
        if path is None and show_file:
            if format is None: raise ValueError("Format has to be specified")
            path = fs.join(introspection.pts_temp_dir, "dust_seds." + format)

        # Run the plotter
        plotter.run(title=title, output=path)

        # Show file
        if show_file: fs.open_file(path)

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
        definition.import_settings(plot_seds_definition)
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

        contributions = config.pop("contributions")

        # Plot
        self.plot_contribution_seds(contributions, **config)

    # -----------------------------------------------------------------

    def plot_contribution_seds(self, contributions, path=None, title=None, show_file=False, format=None, **kwargs):

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
        plotter = SEDPlotter(kwargs) # **kwargs DOESN'T WORK? (e.g. with min_flux)

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
        definition.import_settings(plot_seds_definition)
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

        # Get
        components = config.pop("components")

        # Plot
        self.plot_component_seds(config.components, **config)

    # -----------------------------------------------------------------

    def plot_component_seds(self, components, path=None, title=None, show_file=False, format=None, **kwargs):

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
        plotter = SEDPlotter(kwargs)

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
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED (stellar), or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_old_bulge_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_old_bulge_sed_definition, **kwargs)

        # Either observed or intrinsic
        if config.observed_intrinsic == 1:

            oi = config.observed_intrinsic[0]

            if oi == observed_name: plot_sed(self.model.observed_old_bulge_sed)
            elif oi == intrinsic_name: plot_sed(self.model.intrinsic_sed_old_bulge)
            else: raise ValueError("")

        # Both
        else:

            seds = OrderedDict()
            for oi in config.observed_intrinsic:
                if oi == observed_name: seds[observed_name] = self.model.observed_old_bulge_sed
                elif oi == intrinsic_name: seds[intrinsic_name] = self.model.intrinsic_sed_old_bulge
                else: raise ValueError("")
            plot_seds(seds, residuals=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_old_disk_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED (stellar), or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_old_disk_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_old_disk_sed_definition, **kwargs)

        # Either observed or intrinsic
        if len(config.observed_intrinsic) == 1:

            oi = config.observed_intrinsic[0]

            if oi == observed_name: plot_sed(self.model.observed_old_disk_sed)
            elif oi == intrinsic_name: plot_sed(self.model.intrinsic_sed_old_disk)
            else: raise ValueError("")

        # Both
        else:

            seds = OrderedDict()
            for oi in config.observed_intrinsic:
                if oi == observed_name: seds[observed_name] = self.model.observed_old_disk_sed
                elif oi == intrinsic_name: seds[intrinsic_name] = self.model.intrinsic_sed_old_disk
                else: raise ValueError("")
            plot_seds(seds, residuals=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_old_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED (stellar), or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_old_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_old_sed_definition, **kwargs)

        # Either observed or intrinsic
        if len(config.observed_intrinsic) == 1:

            oi = config.observed_intrinsic[0]

            if oi == observed_name: plot_sed(self.model.observed_old_sed)
            elif oi == intrinsic_name: plot_sed(self.model.intrinsic_sed_old)
            else: raise ValueError("")

        # Both
        else:

            seds = OrderedDict()
            for oi in config.observed_intrinsic:
                if oi == observed_name: seds[observed_name] = self.model.observed_old_sed
                elif oi == intrinsic_name: seds[intrinsic_name] = self.model.intrinsic_sed_old
                else: raise ValueError("")
            plot_seds(seds, residuals=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_young_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED (stellar), or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_young_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.plot_young_sed_definition, **kwargs)

        # Either observed or intrinsic
        if len(config.observed_intrinsic) == 1:

            oi = config.observed_intrinsic[0]

            if oi == observed_name: plot_sed(self.model.observed_young_sed)
            elif oi == intrinsic_name: plot_sed(self.model.intrinsic_sed_young)
            else: raise ValueError("")

        # Both
        else:

            seds = OrderedDict()
            for oi in config.observed_intrinsic:
                if oi == observed_name: seds[observed_name] = self.model.observed_young_sed
                elif oi == intrinsic_name: seds[intrinsic_name] = self.model.intrinsic_sed_young
                else: raise ValueError("")
            plot_seds(seds, residuals=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_sfr_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED (stellar), or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_sfr_sed_command(self, command, **kwargs):

        """
        This function ...
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.plot_sfr_sed_definition, **kwargs)

        # Either observed or intrinsic
        if len(config.observed_intrinsic) == 1:

            oi = config.observed_intrinsic[0]
            if oi == observed_name: plot_sed(self.model.observed_sfr_sed)
            elif oi == intrinsic_name: plot_sed(self.model.intrinsic_sed_sfr)
            else: raise ValueError("")

        # Both
        else:

            seds = OrderedDict()
            for oi in config.observed_intrinsic:
                if oi == observed_name: seds[observed_name] = self.model.observed_sfr_sed
                elif oi == intrinsic_name: seds[intrinsic_name] = self.model.intrinsic_sed_sfr
                else: raise ValueError("")
            plot_seds(seds, residuals=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_unevolved_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_intrinsic", "string_tuple", "plot observed SED, intrinsic SED (stellar), or both", default_observed_intrinsic, choices=observed_intrinsic_choices)
        return definition

    # -----------------------------------------------------------------

    def plot_unevolved_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.plot_unevolved_sed_definition, **kwargs)

        # Either observed or intrinsic
        if len(config.observed_intrinsic) == 1:

            oi = config.observed_intrinsic[0]
            if oi == observed_name: plot_sed(self.model.observed_unevolved_sed)
            elif oi == intrinsic_name: plot_sed(self.model.intrinsic_sed_unevolved)
            else: raise ValueError("")

        # Both
        else:

            seds = OrderedDict()
            for oi in config.observed_intrinsic:
                if oi == observed_name: seds[observed_name] = self.model.observed_unevolved_sed
                elif oi == intrinsic_name: seds[intrinsic_name] = self.model.intrinsic_sed_unevolved
                else: raise ValueError("")
            plot_seds(seds, residuals=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_total_attenuation_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_total_attenuation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_total_attenuation_definition, **kwargs)

        # Plot
        plot_attenuation_curve(self.model.attenuation_curve, total)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_component_attenuation_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("components", "string_list", "components", default_components, choices=components)
        return definition

    # -----------------------------------------------------------------

    def plot_component_attenuation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_component_attenuation_definition, **kwargs)

        # GetAbsorbed bolometric luminosity
        components = config.pop("components")

        # Plot
        self.plot_component_attenuation(components)

    # -----------------------------------------------------------------

    def plot_component_attenuation(self, components):

        """
        This function ...
        :param components:
        :return:
        """

        # Initialize
        curves = OrderedDict()

        # Add components
        for component in components:

            if component == total: curves[total] = self.model.attenuation_curve
            elif component == bulge: curves[bulge] = self.model.attenuation_curve_old_bulge
            elif component == disk: curves[disk] = self.model.attenuation_curve_old_disk
            elif component == old: curves[old] = self.model.attenuation_curve_old
            elif component == young: curves[young] = self.model.attenuation_curve_young
            elif component == sfr: curves[sfr] = self.model.attenuation_curve_sfr
            elif component == unevolved: curves[unevolved] = self.model.attenuation_curve_unevolved
            else: raise ValueError("Invalid component: '" + component + "'")

        # Plot
        plot_attenuation_curves(curves)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_old_bulge_attenuation_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_old_bulge_attenuation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_old_bulge_attenuation_definition, **kwargs)

        # Plot
        plot_attenuation_curve(self.model.attenuation_curve_old_bulge, bulge)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_old_disk_attenuation_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_old_disk_attenuation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_old_disk_attenuation_definition, **kwargs)

        # Plot
        plot_attenuation_curve(self.model.attenuation_curve_old_disk, disk)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_old_attenuation_definition(self):

        """
        Thisf unction ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_old_attenuation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_old_attenuation_definition, **kwargs)

        # Plot
        plot_attenuation_curve(self.model.attenuation_curve_old, old)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_young_attenuation_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_young_attenuation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_young_attenuation_definition, **kwargs)

        # Plot
        plot_attenuation_curve(self.model.attenuation_curve_young, young)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_sfr_attenuation_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_sfr_attenuation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_sfr_attenuation_definition, **kwargs)

        # Plot
        plot_attenuation_curve(self.model.attenuation_curve_sfr, sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_unevolved_attenuation_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_unevolved_attenuation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_unevolved_attenuation_definition, **kwargs)

        # Plot
        plot_attenuation_curve(self.model.attenuation_curve_unevolved, unevolved)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_images_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Add options
        definition.add_positional_optional("orientation", "string", "orientation of the images", earth_name, choices=orientations)
        #definition.add_optional("filters", "lazy_broad_band_filter_list", "filters for which to plot images", default="GALEX,SDSS,IRAC,Mips 24mu,Herschel", convert_default=True)
        definition.add_optional("filters", "lazy_broad_band_filter_list", "filters for which to plot images", default="FUV,NUV,I1,24mu,Pacs160,SPIRE350", convert_default=True)
        definition.add_flag("residuals", "show residuals", True)
        definition.add_flag("distributions", "show residual distributions", True)
        definition.add_flag("from_evaluation", "use the images created in the evaluation step", None)
        definition.add_flag("spectral_convolution", "use spectral convolution to create images", False)
        definition.add_flag("proper", "use the proper mock observed images if present", True)
        definition.add_flag("only_from_evaluation", "only use filters for which an image was made in the evaluation")
        definition.add_flag("sort_filters", "sort the filters on wavelength", True)
        definition.add_optional("path", "string", "path for the plot file")

        # Return
        return definition

    # -----------------------------------------------------------------

    def plot_images_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_images_definition, **kwargs)

        # Get list of filters
        if config.only_from_evaluation:
            config.from_evaluation = True
            if config.proper: filters = sequences.intersection(config.filters, self.analysis_run.evaluation_proper_image_filters)
            else: filters = sequences.intersection(config.filters, self.analysis_run.evaluation_image_filters)
        else: filters = config.filters

        # Sort the filters on wavelength
        if config.sort_filters: filters = sequences.sorted_by_attribute(filters, "wavelength")

        # Earth
        if config.orientation == earth_name: self.plot_earth_images(filters, residuals=config.residuals,
                                                                    distributions=config.distributions, from_evaluation=config.from_evaluation,
                                                                    spectral_convolution=config.spectral_convolution, proper=config.proper, path=config.path)

        # Face-on
        elif config.orientation == faceon_name: self.plot_faceon_images(filters, spectral_convolution=config.spectral_convolution, path=config.path)

        # Edge-on
        elif config.orientation == edgeon_name: self.plot_edgeon_images(filters, spectral_convolution=config.spectral_convolution, path=config.path)

        # Invalid
        else: raise ValueError("Invalid orientation: '" + config.orientation + "'")

    # -----------------------------------------------------------------

    @property
    def earth_cube(self):

        """
        This function ...
        :return:
        """

        return self.model.total_bolometric_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def faceon_cube(self):

        """
        This function ...
        :return:
        """

        return self.model.total_bolometric_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def edgeon_cube(self):

        """
        This function ...
        :return:
        """

        return self.model.total_bolometric_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @memoize_method
    def get_earth_image(self, filter_or_wavelength, from_evaluation=None, spectral_convolution=False, proper=True):

        """
        This function ...
        :param filter_or_wavelength:
        :param from_evaluation:
        :param spectral_convolution:
        :param proper:
        :return:
        """

        # Filter?
        if isinstance(filter_or_wavelength, Filter):

            # From evaluation
            if from_evaluation is None:

                if proper:
                    if self.analysis_run.has_evaluation_proper_image_for_filter(filter_or_wavelength): from_evaluation = True
                    else: from_evaluation = False
                else:
                    if self.analysis_run.has_evaluation_image_for_filter(filter_or_wavelength): from_evaluation = True
                    else: from_evaluation = False

            # From evaluation
            if from_evaluation:

                if proper: return self.analysis_run.get_evaluation_proper_image_for_filter(filter_or_wavelength)
                else: return self.analysis_run.get_evaluation_image_for_filter(filter_or_wavelength)

            # Not from evaluation
            else: return self.earth_cube.frame_for_filter(filter_or_wavelength, convolve=spectral_convolution)

        # Wavelength
        elif types.is_length_quantity(filter_or_wavelength):

            # Checks
            if spectral_convolution: raise ValueError("Spectral convolution cannot be applied when a wavelength is passed")
            if from_evaluation: raise ValueError("Cannot get image for a particular wavelength from evaluation output")

            # Return the frame for this wavelength
            return self.earth_cube.get_frame_for_wavelength(filter_or_wavelength)

        # Invalid
        else: raise ValueError("Invalid argument")

    # -----------------------------------------------------------------

    def has_residual_image_from_evaluation(self, fltr, proper=True):

        """
        Thisf unction ...
        :param fltr:
        :param proper:
        :return:
        """

        if proper: return self.analysis_run.has_evaluation_proper_residuals_for_filter(fltr)
        else: return self.analysis_run.has_evaluation_residuals_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_residual_image_from_evaluation(self, fltr, proper=True):

        """
        This function ...
        :param fltr:
        :param proper:
        :return:
        """

        if proper: return self.analysis_run.get_evaluation_proper_residuals_for_filter(fltr)
        else: return self.analysis_run.get_evaluation_residuals_for_filter(fltr)

    # -----------------------------------------------------------------

    def has_residual_distribution_from_evaluation(self, fltr, proper=True):

        """
        This function ...
        :param fltr:
        :param proper:
        :return:
        """

        if proper: return self.analysis_run.has_evaluation_proper_residuals_distribution_for_filter(fltr)
        else: return self.analysis_run.has_evaluation_residuals_distribution_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_residual_distribution_from_evaluation(self, fltr, proper=True):

        """
        This function ...
        :param fltr:
        :param proper:
        :return:
        """

        if proper: return self.analysis_run.get_evaluation_proper_residuals_distribution_for_filter(fltr)
        else: return self.analysis_run.get_evaluation_residuals_distribution_for_filter(fltr)

    # -----------------------------------------------------------------

    @memoize_method
    def get_faceon_image(self, filter_or_wavelength, spectral_convolution=False):

        """
        This fnuction ...
        :param filter_or_wavelength:
        :param spectral_convolution:
        :return:
        """

        # Filter?
        if isinstance(filter_or_wavelength, Filter): return self.faceon_cube.frame_for_filter(filter_or_wavelength, convolve=spectral_convolution)

        # Wavelength
        elif types.is_length_quantity(filter_or_wavelength):
            if spectral_convolution: raise ValueError("Spectral convolution cannot be applied when a wavelength is passed")
            return self.faceon_cube.get_frame_for_wavelength(filter_or_wavelength)

        # Invalid
        else: raise ValueError("Invalid argument")

    # -----------------------------------------------------------------

    @memoize_method
    def get_edgeon_image(self, filter_or_wavelength, spectral_convolution=False):

        """
        This function ...
        :param filter_or_wavelength:
        :param spectral_convolution:
        :return:
        """

        # Filter?
        if isinstance(filter_or_wavelength, Filter): return self.edgeon_cube.frame_for_filter(filter_or_wavelength, convolve=spectral_convolution)

        # Wavelength
        elif types.is_length_quantity(filter_or_wavelength):
            if spectral_convolution: raise ValueError("Spectral convolution cannot be applied when a wavelength is passed")
            return self.edgeon_cube.get_frame_for_wavelength(filter_or_wavelength)

        # Invalid
        else: raise ValueError("Invalid argument")

    # -----------------------------------------------------------------

    def plot_earth_images(self, filters, residuals=True, distributions=True, from_evaluation=None,
                          spectral_convolution=False, proper=True, path=None):

        """
        Thisf unction ...
        :param filters:
        :param residuals:
        :param distributions:
        :param from_evaluation:
        :param spectral_convolution:
        :param proper:
        :param path:
        :return:
        """

        # Debugging
        log.debug("Plotting the observed and model images in the earth projection ...")

        # Create the plotter
        plotter = ResidualImageGridPlotter()

        # Set options
        plotter.config.distributions = distributions

        # Set the output filepath
        plotter.config.path = path

        # Loop over the filters
        for fltr in filters:

            # Define the image name
            image_name = str(fltr)

            # Get the frame
            observation = self.get_photometry_frame_for_filter(fltr)

            # Replace zeroes and negatives
            observation.replace_zeroes_by_nans()
            observation.replace_negatives_by_nans()

            # Add the frame to the plotter
            plotter.add_observation(image_name, observation)

            # Get modeled frame
            modeled = self.get_earth_image(fltr, from_evaluation=from_evaluation, spectral_convolution=spectral_convolution, proper=proper)

            # Replace zeroes and negatives
            modeled.replace_zeroes_by_nans()
            modeled.replace_negatives_by_nans()

            # Add the mock image to the plotter
            plotter.add_model(image_name, modeled)

            # Add residuals if present
            if residuals and self.has_residual_image_from_evaluation(fltr, proper=proper):
                residuals = self.get_residual_image_from_evaluation(fltr, proper=proper)
                plotter.add_residuals(image_name, residuals)

            # Add residuals distribution if present
            if distributions and self.has_residual_distribution_from_evaluation(fltr, proper=proper):
                distribution = self.get_residual_distribution_from_evaluation(fltr, proper=proper)
                plotter.add_distribution(image_name, distribution)

        # Run the plotter
        plotter.run()

    # -----------------------------------------------------------------

    def plot_faceon_images(self, filters, spectral_convolution=False, path=None):

        """
        This function ...
        :param filters:
        :param spectral_convolution:
        :param path:
        :return:
        """

        # Inform the user
        log.info("Plotting the model images in the face-on projection ...")

        # Create the plotter
        plotter = StandardImageGridPlotter()

        # Set the output filepath
        plotter.config.path = path

        # Loop over the filters
        for fltr in filters:

            # Determine name
            image_name = str(fltr)

            # Get frame
            frame = self.get_faceon_image(fltr)

            # Replace zeroes and negatives
            frame.replace_zeroes_by_nans()
            frame.replace_negatives_by_nans()

            # Add the mock image to the plotter
            plotter.add_image(image_name, frame)

            # Run the plotter
        plotter.run()

    # -----------------------------------------------------------------

    def plot_edgeon_images(self, filters, spectral_convolution=False, path=None):

        """
        This function ...
        :param filters:
        :param spectral_convolution:
        :param path:
        :return:
        """

        # Inform the user
        log.info("Plotting the model images in the edge-on projection ...")

        # Create the plotter
        plotter = StandardImageGridPlotter()

        # Set the output filepath
        plotter.config.path = path

        # Loop over the filters
        for fltr in filters:

            # Determine name
            image_name = str(fltr)

            # Get frame
            frame = self.get_edgeon_image(fltr)

            # Replace zeroes and negatives
            frame.replace_zeroes_by_nans()
            frame.replace_negatives_by_nans()

            # Add the mock image to the plotter
            plotter.add_image(image_name, frame)

            # Run the plotter
        plotter.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_cubes_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Add options
        definition.add_positional_optional("orientation", "string", "orientation of the datacube")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def plot_cubes_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_cubes_definition, **kwargs)

        # Earth
        if config.orientation == earth_name: self.plot_earth_cube()

        # Face-on
        elif config.orientation == faceon_name: self.plot_faceon_cube()

        # Edge-on
        elif config.orientation == edgeon_name: self.plot_edgeon_cube()

        # Invalid
        else: raise ValueError("Invalid orientation: '" + config.orientation + "'")

    # -----------------------------------------------------------------

    def plot_earth_cube(self):

        """
        This function ...
        :return:
        """

        #from ...magic.tools import plotting
        #from ...magic.core.datacube import DataCube

        # Get simulation prefix
        #prefix = self.get_simulation_prefix(simulation_name)

        # Get the wavelength grid
        #wavelength_grid = self.get_wavelength_grid(simulation_name)

        # Load the datacube
        #datacube = DataCube.from_file(path, wavelength_grid)

        # Plot
        plotting.plot_datacube(datacube, title=instr_name, share_normalization=share_normalization, show_axes=False)

    # -----------------------------------------------------------------

    def plot_faceon_cube(self):

        """
        Thisn function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_edgeon_cube(self):

        """
        Thisn function ...
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluate_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.import_settings(evaluate_analysis_definition)
        return definition

    # -----------------------------------------------------------------

    def evaluate_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.evaluate_definition, **kwargs)

        # Evaluate
        self.evaluate(**config)

    # -----------------------------------------------------------------

    def evaluate(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Create evaluator
        evaluator = AnalysisModelEvaluator(**kwargs)

        # Set modeling path
        evaluator.config.path = self.config.path

        # Set analysis run name
        evaluator.config.run = self.config.run

        # Run
        evaluator.run()

    # -----------------------------------------------------------------

    def show_map(self, frame, contours=False, ncontours=5):

        """
        This function ...
        :param frame:
        :param contours:
        :param ncontours:
        :return:
        """

        # With contours
        if contours: plot_frame_contours(frame, nlevels=ncontours)

        # Just frame
        else: plot_frame(frame)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_total_map_definition(self):

        """
        Thisn function ...
        :return:
        """

        # Create
        definition = ConfigurationDefinition(write_config=False)

        # Options
        definition.add_required("which", "string", "which map to plot", choices=total_map_names)
        definition.add_positional_optional("orientation", "string", "orientation of the map", default=earth_name, choices=orientations)
        definition.add_flag("contours", "show contours", False)
        definition.add_optional("ncontours", "positive_integer", "number of contour lines", 5)

        # Return
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def show_total_map_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    def show_total_map_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set kwargs
        kwargs.update(self.show_total_map_kwargs)

        # Get the configuration
        config = self.get_config_from_command(command, definition=self.show_total_map_definition, **kwargs)

        # Show
        self.show_total_map(config.which, orientation=config.orientation)

    # -----------------------------------------------------------------

    @property
    def total_simulations(self):

        """
        This function ...
        :return:
        """

        return self.model.total_simulations

    # -----------------------------------------------------------------

    @property
    def bulge_simulations(self):

        """
        This function ...
        :return:
        """

        return self.model.bulge_simulations

    # -----------------------------------------------------------------

    @property
    def disk_simulations(self):

        """
        This function ...
        :return:
        """

        return self.model.disk_simulations

    # -----------------------------------------------------------------

    @property
    def old_simulations(self):

        """
        Thisn function ...
        :return:
        """

        return self.model.old_simulations

    # -----------------------------------------------------------------

    @property
    def young_simulations(self):

        """
        Thisfunction ...
        :return:
        """

        return self.model.young_simulations

    # -----------------------------------------------------------------

    @property
    def sfr_simulations(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_simulations

    # -----------------------------------------------------------------

    @property
    def unevolved_simulations(self):

        """
        This function ...
        :return:
        """

        return self.model.unevolved_simulations

    # -----------------------------------------------------------------

    def show_total_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Debugging
        log.debug("Showing the " + which + " map of the total model from the " + orientation + " orientation ...")

        # Get the map
        frame = self.get_total_map(which, orientation=orientation)

        # Show
        self.show_map(frame)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_bulge_map_definition(self):

        """
        This function ...
        :return:
        """

        # Create
        definition = ConfigurationDefinition(write_config=False)

        # Options
        definition.add_required("which", "string", "which map to plot", choices=bulge_map_names)
        definition.add_positional_optional("orientation", "string", "orientation of the map", default=earth_name, choices=orientations)
        definition.add_flag("contours", "show contours", False)
        definition.add_optional("ncontours", "positive_integer", "number of contour lines", 5)

        # Return
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def show_bulge_map_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    def show_bulge_map_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set kwargs
        kwargs.update(self.show_bulge_map_kwargs)

        # Get the configuration
        config = self.get_config_from_command(command, self.show_bulge_map_definition, **kwargs)

        # Show
        self.show_bulge_map(config.which, orientation=config.orientation)

    # -----------------------------------------------------------------

    def show_bulge_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Debugging
        log.debug("Showing the " + which + " map of the old stellar bulge component from the " + orientation + " orientation ...")

        # Get the map
        frame = self.get_bulge_map(which, orientation=orientation)

        # Show the map
        self.show_map(frame)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_disk_map_definition(self):

        """
        This function ...
        :return:
        """

        # Create
        definition = ConfigurationDefinition(write_config=False)

        # Options
        definition.add_required("which", "string", "which map to plot", choices=disk_map_names)
        definition.add_positional_optional("orientation", "string", "orientation of the map", default=earth_name, choices=orientations)
        definition.add_flag("contours", "show contours", False)
        definition.add_optional("ncontours", "positive_integer", "number of contour lines", 5)

        # Return
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def show_disk_map_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    def show_disk_map_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set kwargs
        kwargs.update(self.show_disk_map_kwargs)

        # Get the configuration
        config = self.get_config_from_command(command, self.show_disk_map_definition, **kwargs)

        # Show
        self.show_disk_map(config.which, orientation=config.orientation)

    # -----------------------------------------------------------------

    def show_disk_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Debugging
        log.debug("Showing the " + which + " map of the old stellar disk component from the " + orientation + " orientation ...")

        # Get the disk map
        frame = self.get_disk_map(which, orientation=orientation)

        # Show
        self.show_map(frame)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_old_map_definition(self):

        """
        This function ...
        :return:
        """

        # Create
        definition = ConfigurationDefinition(write_config=False)

        # Options
        definition.add_required("which", "string", "which map to plot", choices=old_map_names)
        definition.add_positional_optional("orientation", "string", "orientation of the map", default=earth_name, choices=orientations)
        definition.add_flag("contours", "show contours", False)
        definition.add_optional("ncontours", "positive_integer", "number of contour lines", 5)

        # Return
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def show_old_map_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    def show_old_map_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set kwargs
        kwargs.update(self.show_old_map_kwargs)

        # Get the configuration
        config = self.get_config_from_command(command, self.show_old_map_definition, **kwargs)

        # Show
        self.show_old_map(config.which, orientation=config.orientation)

    # -----------------------------------------------------------------

    def show_old_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Debugging
        log.debug("Showing the " + which + " map of the old stellar component from the " + orientation + " orientation ...")

        # Get the map
        frame = self.get_old_map(which, orientation=orientation)

        # Show
        self.show_map(frame)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_young_map_definition(self):

        """
        This function ...
        :return:
        """

        # Create
        definition = ConfigurationDefinition(write_config=False)

        # Options
        definition.add_required("which", "string", "which map to plot", choices=young_map_names)
        definition.add_positional_optional("orientation", "string", "orientation of the map", default=earth_name, choices=orientations)
        definition.add_flag("contours", "show contours", False)
        definition.add_optional("ncontours", "positive_integer", "number of contour lines", 5)

        # Return
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def show_young_map_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    def show_young_map_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set kwargs
        kwargs.update(self.show_young_map_kwargs)

        # Get the configuration
        config = self.get_config_from_command(command, self.show_young_map_definition, **kwargs)

        # Show
        self.show_young_map(config.which, orientation=config.orientation)

    # -----------------------------------------------------------------

    def show_young_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Debugging
        log.debug("Showing the " + which + " map of the young stellar component from the " + orientation + " orientation ...")

        # Get the map
        frame = self.get_young_map(which, orientation=orientation)

        # Show
        self.show_map(frame)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_sfr_map_definition(self):

        """
        This function ...
        :return:
        """

        # Create
        definition = ConfigurationDefinition(write_config=False)

        # Options
        definition.add_required("which", "string", "which map to plot", choices=sfr_map_names)
        definition.add_positional_optional("orientation", "string", "orientation of the map", default=earth_name, choices=orientations)
        definition.add_flag("contours", "show contours", False)
        definition.add_optional("ncontours", "positive_integer", "number of contour lines", 5)

        # Return
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def show_sfr_map_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    def show_sfr_map_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set kwargs
        kwargs.update(self.show_sfr_map_kwargs)

        # Get the configuration
        config = self.get_config_from_command(command, self.show_sfr_map_definition, **kwargs)

        # Show
        self.show_sfr_map(config.which, orientation=config.orientation)

    # -----------------------------------------------------------------

    def show_sfr_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Debugging
        log.debug("Showing the " + which + " map of the SFR stellar component from the " + orientation + " orientation ...")

        # Get the map
        frame = self.get_sfr_map(which, orientation=orientation)

        # Show
        self.show_map(frame)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_unevolved_map_definition(self):

        """
        This function ...
        :return:
        """

        # Create
        definition = ConfigurationDefinition(write_config=False)

        # Options
        definition.add_required("which", "string", "which map to plot", choices=unevolved_map_names)
        definition.add_positional_optional("orientation", "string", "orientation of the map", default=earth_name, choices=orientations)
        definition.add_flag("contours", "show contours", False)
        definition.add_optional("ncontours", "positive_integer", "number of contour lines", 5)

        # Return
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def show_unevolved_map_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    def show_unevolved_map_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set kwargs
        kwargs.update(self.show_unevolved_map_kwargs)

        # Get the configuration
        config = self.get_config_from_command(command, self.show_unevolved_map_definition, **kwargs)

        # Show
        self.show_unevolved_map(config.which, orientation=config.orientation)

    # -----------------------------------------------------------------

    def show_unevolved_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Debugging
        log.debug("Showing the " + which + " map of the unevolved stellar component from the " + orientation + " orientation ...")

        # Get the map
        frame = self.get_unevolved_map(which, orientation=orientation)

        # Show
        self.show_map(frame)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_dust_map_definition(self):

        """
        This function ...
        :return:
        """

        # Create
        definition = ConfigurationDefinition(write_config=False)

        # Options
        definition.add_required("which", "string", "which map to plot", choices=dust_map_names)
        definition.add_positional_optional("orientation", "string", "orientation of the map", default=earth_name, choices=orientations)
        definition.add_flag("contours", "show contours", False)
        definition.add_optional("ncontours", "positive_integer", "number of contour lines", 5)

        # Return
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def show_dust_map_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    def show_dust_map_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set kwargs
        kwargs.update(self.show_dust_map_kwargs)

        # Get the configuration
        config = self.get_config_from_command(command, self.show_dust_map_definition, **self.show_dust_map_kwargs)

        # Show
        self.show_dust_map(config.which, orientation=config.orientation)

    # -----------------------------------------------------------------

    def show_dust_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Debugging
        log.debug("Showing the " + which + " dust map from the " + orientation + " orientation ...")

        # Get the map
        frame = self.get_dust_map(which, orientation=orientation)

        # Show
        self.show_map(frame)

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_properties_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Add settings
        definition.import_settings(analyse_properties_definition)
        definition.remove_setting("run")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analyse_properties_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_properties_definition, **kwargs)

        # Analyse
        self.analyse_properties(config=config)

    # -----------------------------------------------------------------

    def analyse_properties(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Create the analyser
        analyser = PropertiesAnalyser(config=config)

        # Set the modeling path
        analyser.config.path = self.config.path

        # Set the analysis run
        analyser.config.run = self.config.run

        # Run
        analyser.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_cell_heating_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        #print(analyse_cell_heating_definition.property_names)
        #print(analyse_cell_heating_definition.section_names)

        # Add settings
        definition.import_settings(analyse_cell_heating_definition)
        definition.remove_setting("run")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analyse_cell_heating_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_cell_heating_definition, **kwargs)

        # Analyse
        self.analyse_cell_heating(config=config)

    # -----------------------------------------------------------------

    def analyse_cell_heating(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Create the analyser
        analyser = CellDustHeatingAnalyser(config=config)

        # Set the modeling path
        analyser.config.path = self.config.path

        # Set the analysis run
        analyser.config.run = self.config.run

        # Run
        analyser.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_projected_heating_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        #print(analyse_projected_heating_definition.property_names)
        #print(analyse_projected_heating_definition.section_names)

        # Add settings
        definition.import_settings(analyse_projected_heating_definition)
        definition.remove_setting("run")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analyse_projected_heating_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_projected_heating_definition, **kwargs)

        # Analyse
        self.analyse_projected_heating(config=config)

    # -----------------------------------------------------------------

    def analyse_projected_heating(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Create the analyser
        analyser = ProjectedDustHeatingAnalyser(config=config)

        # Set the modeling path
        analyser.config.path = self.config.path

        # Set the analysis run
        analyser.config.run = self.config.run

        # Run
        analyser.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_cell_energy_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Add settings
        definition.import_settings(analyse_cell_energy_definition)
        definition.remove_setting("run")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analyse_cell_energy_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_cell_energy_definition, **kwargs)

        # Analyse
        self.analyse_cell_energy(config=config)

    # -----------------------------------------------------------------

    def analyse_cell_energy(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Create the analyser
        analyser = CellEnergyAnalyser(config=config)

        # Set the modeling path
        analyser.config.path = self.config.path

        # Set the analysis run
        analyser.config.run = self.config.run

        # Run
        analyser.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_projected_energy_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Add settings
        definition.import_settings(analyse_projected_energy_definition)
        definition.remove_setting("run")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analyse_projected_energy_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_projected_energy_definition, **kwargs)

        # Analyse
        self.analyse_projected_energy(config=config)

    # -----------------------------------------------------------------

    def analyse_projected_energy(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Create the analyser
        analyser = ProjectedEnergyAnalyser(config=config)

        # Set the modeling path
        analyser.config.path = self.config.path

        # Set the analysis run
        analyser.config.run = self.config.run

        # Run
        analyser.run()

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
