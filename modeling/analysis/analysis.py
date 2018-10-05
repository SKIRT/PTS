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
from .absorption import AbsorptionAnalyser
from ..config.analyse_absorption import definition as analyse_absorption_definition
from ...core.plot.sed import plot_seds, SEDPlotter, plot_sed
from ...core.config.plot_seds import definition as plot_seds_definition
from ..config.evaluate_analysis import definition as evaluate_analysis_definition
from ...core.plot.attenuation import plot_attenuation_curve, plot_attenuation_curves
from ..config.analyse_cell_heating import definition as analyse_cell_heating_definition
from ..config.analyse_projected_heating import definition as analyse_projected_heating_definition
from ..config.analyse_spectral_heating import definition as analyse_spectral_heating_definition
from ..config.analyse_images import definition as analyse_images_definition
from ..config.analyse_residuals import definition as analyse_residuals_definition
from .heating.cell import CellDustHeatingAnalyser
from .heating.projected import ProjectedDustHeatingAnalyser
from .heating.spectral import SpectralDustHeatingAnalyser
from .images import ImagesAnalyser
from .residuals import ResidualAnalyser
from ..config.analyse_properties import definition as analyse_properties_definition
from .properties import PropertiesAnalyser
from ..config.analyse_cell_energy import definition as analyse_cell_energy_definition
from ..config.analyse_projected_energy import definition as analyse_projected_energy_definition
from .energy.cell import CellEnergyAnalyser
from .energy.projected import ProjectedEnergyAnalyser
from ...magic.tools.plotting import plot_frame, plot_frame_contours, plot_datacube
from ...core.filter.filter import Filter, parse_filter
from ...core.tools import types
from ...magic.plot.imagegrid import StandardImageGridPlotter, ResidualImageGridPlotter
from .evaluation import AnalysisModelEvaluator
from ...core.tools import sequences
from .correlations import CorrelationsAnalyser
from ..misc.examination import ModelExamination
from ..config.analyse_correlations import definition as analyse_correlations_definition
from ..config.analyse_sfr import definition as analyse_sfr_definition
from .sfr import SFRAnalyser
from ...core.units.parsing import parse_unit as u
from ...core.data.sed import SED
from ...magic.core.dataset import StaticDataSet
from ...core.basics.distribution import Distribution

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
_show_command_name = "show"
_properties_command_name = "properties"
_output_command_name = "output"
_data_command_name = "data"
_model_command_name = "model"

# Plot commands
_plot_command_name = "plot"
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
_absorption_command_name = "absorption"
_heating_command_name = "heating"
_energy_command_name = "energy"
_sfr_command_name = "sfr"
_correlations_command_name = "correlations"
_residuals_command_name = "residuals"

# -----------------------------------------------------------------

# Define commands
commands = OrderedDict()

# Standard commands
commands[_help_command_name] = ("show_help", False, "show help", None)
commands[_history_command_name] = ("show_history_command", True, "show history of executed commands", None)
commands[_status_command_name] = ("show_status_command", True, "show analysis status", None)

# Show stuff
commands[_show_command_name] = (None, None, "show analysis results", None)

# Examine the model
commands[_model_command_name] = ("examine_model", False, "examine the radiative transfer model", None)

# Plot stuff
commands[_sed_command_name] = (None, None, "plot SEDs", None)
commands[_attenuation_command_name] = (None, None, "plot attenuation curves", None)
commands[_map_command_name] = (None, None, "plot a map", None)
commands[_plot_command_name] = (None, None, "plot other stuff", None)

# Evaluate
commands[_evaluate_command_name] = ("evaluate_command", True, "evaluate the analysis model", None)

# Analysis
commands[_properties_command_name] = ("analyse_properties_command", True, "analyse the model properties", None)
commands[_absorption_command_name] = ("analyse_absorption_command", True, "analyse the dust absorption", None)
commands[_heating_command_name] = (None, None, "analyse dust heating contributions", None)
commands[_energy_command_name] = (None, None, "analyse the energy budget in the galaxy", None)
commands[_sfr_command_name] = ("analyse_sfr_command", True, "analyse the star formation rates", None)
commands[_correlations_command_name] = ("analyse_correlations_command", True, "analyse the correlations", None)
commands[_images_command_name] = ("analyse_images_command", True, "analyse the simulation images", None)
commands[_residuals_command_name] = ("analyse_residuals_command", True, "analyse the image residuals", None)

# -----------------------------------------------------------------

_bulge_name = "bulge"
_disk_name = "disk"

_total_name = "total"
_old_bulge_name = "old_bulge"
_old_disk_name = "old_disk"
_old_name = "old"
_young_name = "young"
_sfr_name = "sfr"
_sfr_intrinsic_name = "sfr_intrinsic"
_unevolved_name = "unevolved"

#_sfr_stellar_name = "sfr_stellar"
#_sfr_dust_name = "sfr_dust"

_stellar_name = "stellar"
_dust_name = "dust"

_contributions_name = "contributions"
_components_name = "components"

_absorption_name = "absorption"

# -----------------------------------------------------------------

# Show subcommands
show_commands = OrderedDict()

# Properties
show_commands[_properties_command_name] = ("show_properties", False, "show the model properties", None)

# Simulation output and data
show_commands[_output_command_name] = ("show_output", False, "show the simulation output", None)
show_commands[_data_command_name] = ("show_data", False, "show the simulation data available for the model", None)

# -----------------------------------------------------------------

# Plot subcommands
plot_commands = OrderedDict()
plot_commands[_wavelengths_command_name] = ("plot_wavelengths_command", True, "plot the wavelength grid", None)
plot_commands[_dustgrid_command_name] = ("plot_grid_command", True, "plot the dust grid", None)
plot_commands[_residuals_command_name] = ("plot_residuals_command", True, "plot the observed, modeled and residual images", None)
plot_commands[_images_command_name] = ("plot_images_command", True, "plot the simulated images", None)
plot_commands[_cubes_command_name] = ("plot_cubes_command", True, "plot the simulated datacubes", None)

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

## INTRINSIC SFR
sed_commands[_sfr_intrinsic_name] = ("plot_sfr_intrinsic_sed_command", True, "plot the intrinsic (stellar and dust) SED of the star formation regions", None)

## UNEVOLVED
sed_commands[_unevolved_name] = ("plot_unevolved_sed_command", True, "plot the SED of the unevolved stellar population (young + sfr)", None)

# ABSORPTION
sed_commands[_absorption_name] = ("plot_absorption_sed_command", True, "plot absorption SEDs", None)

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
map_commands[_total_name] = ("show_total_map_command", True, "show a map of the total model", None)
map_commands[_bulge_name] = ("show_bulge_map_command", True, "show a map of the old stellar bulge component", None)
map_commands[_disk_name] = ("show_disk_map_command", True, "show a map of the old stellar disk component", None)
map_commands[_old_name] = ("show_old_map_command", True, "show a map of the old stellar component", None)
map_commands[_young_name] = ("show_young_map_command", True, "show a map of the young stellar component", None)
map_commands[_sfr_name] = ("show_sfr_map_command", True, "show a map of the SFR component", None)
map_commands[_unevolved_name] = ("show_unevolved_map_command", True, "show a map of the unevolved stellar component", None)
map_commands[_dust_name] = ("show_dust_map_command", True, "show a map of the dust component", None)

# -----------------------------------------------------------------

_cell_name = "cell"
_projected_name = "projected"
_spectral_name = "spectral"

# -----------------------------------------------------------------

# Heating subcommands
heating_commands = OrderedDict()

# Cell and projected
heating_commands[_cell_name] = ("analyse_cell_heating_command", True, "analyse the cell heating", None)
heating_commands[_projected_name] = ("analyse_projected_heating_command", True, "analyse the projected heating", None)
heating_commands[_spectral_name] = ("analyse_spectral_heating_command", True, "analyse the spectral heating", None)

# -----------------------------------------------------------------

# Energy subcommands
energy_commands = OrderedDict()

# Cell and projected
energy_commands[_cell_name] = ("analyse_cell_energy_command", True, "analyse the cell energy budget", None)
energy_commands[_projected_name] = ("analyse_projected_energy_command", True, "analyse the projected energy budget", None)

# -----------------------------------------------------------------

# Set subcommands
subcommands = OrderedDict()
subcommands[_show_command_name] = show_commands
subcommands[_plot_command_name] = plot_commands
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
stellar_name = "stellar"
intrinsic_name = "intrinsic"

# -----------------------------------------------------------------

default_observed_intrinsic = (observed_name, intrinsic_name)
observed_intrinsic_choices = default_observed_intrinsic

default_observed_stellar_intrinsic = (observed_name, intrinsic_name)
observed_stellar_intrinsic_choices = [observed_name, stellar_name, intrinsic_name]

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

# Photometric quantity
flux_name = "flux"
luminosity_name = "luminosity"
photometric_quantity_names = [flux_name, luminosity_name]
default_photometric_quantity_name = flux_name

# Spectral style
wavelength_style_name = "wavelength"
frequency_style_name = "frequency"
neutral_style_name = "neutral"
spectral_style_names = [wavelength_style_name, frequency_style_name, neutral_style_name]
default_spectral_style = wavelength_style_name

# -----------------------------------------------------------------

default_plotting_format = "pdf"

# -----------------------------------------------------------------

from ..core.model import contributions, total_contribution, direct_contribution, scattered_contribution, dust_contribution, transparent_contribution
from ..core.model import dust_direct_contribution, dust_scattered_contribution
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
        return self.config.commands is not None and len(self.config.commands) > 0

    # -----------------------------------------------------------------

    @property
    def do_interactive(self):
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
    # PHOTOMETRIC UNITS
    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_lum_unit(self):
        return u("W/micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def frequency_lum_unit(self):
        return u("W/Hz")

    # -----------------------------------------------------------------

    @lazyproperty
    def neutral_lum_unit(self):
        return u("Lsun", density=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_flux_unit(self):
        return u("W/m2/micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def frequency_flux_unit(self):
        return u("Jy")

    # -----------------------------------------------------------------

    @lazyproperty
    def neutral_flux_unit(self):
        return u("W/m2", density=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def luminosity_units(self):
        return {wavelength_style_name: self.wavelength_lum_unit, frequency_style_name: self.frequency_lum_unit, neutral_style_name: self.neutral_lum_unit}

    # -----------------------------------------------------------------

    @lazyproperty
    def flux_units(self):
        return {wavelength_style_name: self.wavelength_flux_unit, frequency_style_name: self.frequency_flux_unit, neutral_style_name: self.neutral_flux_unit}

    # -----------------------------------------------------------------

    @lazyproperty
    def photometric_units(self):
        return {luminosity_name: self.luminosity_units, flux_name: self.flux_units}

    # -----------------------------------------------------------------

    @property
    def simulations(self):
        return self.model.simulations

    # -----------------------------------------------------------------

    @property
    def parameter_values(self):
        return self.model.parameter_values

    # -----------------------------------------------------------------

    @property
    def free_parameter_values(self):
        return self.model.free_parameter_values

    # -----------------------------------------------------------------

    @property
    def other_parameter_values(self):
        return self.model.other_parameter_values

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values(self):
        return self.model.derived_parameter_values

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_total(self):
        return self.model.derived_parameter_values_total

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_bulge(self):
        return self.model.derived_parameter_values_bulge

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_disk(self):
        return self.model.derived_parameter_values_disk

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_old(self):
        return self.model.derived_parameter_values_old

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_young(self):
        return self.model.derived_parameter_values_young

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_sfr(self):
        return self.model.derived_parameter_values_sfr

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_unevolved(self):
        return self.model.derived_parameter_values_unevolved

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_dust(self):
        return self.model.derived_parameter_values_dust

    # -----------------------------------------------------------------

    @property
    def generation_name(self):
        return self.analysis_run.generation_name

    # -----------------------------------------------------------------

    @property
    def simulation_name(self):
        return self.analysis_run.simulation_name

    # -----------------------------------------------------------------

    @property
    def chi_squared(self):
        return self.analysis_run.chi_squared

    # -----------------------------------------------------------------

    @property
    def fitting_run_name(self):
        return self.analysis_run.fitting_run_name

    # -----------------------------------------------------------------

    @property
    def fitting_run(self):
        return self.analysis_run.fitting_run

    # -----------------------------------------------------------------

    @property
    def from_fitting(self):
        return self.analysis_run.from_fitting

    # -----------------------------------------------------------------

    @property
    def wavelength_grid(self):
        return self.analysis_run.wavelength_grid

    # -----------------------------------------------------------------

    @property
    def dust_grid(self):
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

    def show_properties(self, **kwargs):

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

    def show_output(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Showing the simulation output ...")

        # TOTAL
        print(fmt.blue + fmt.underlined + "TOTAL" + fmt.reset + ":")
        print("")
        self.total_output.show(line_prefix="  ", dense=True)
        print("")

        # BULGE
        print(fmt.blue + fmt.underlined + "BULGE" + fmt.reset + ":")
        print("")
        self.bulge_output.show(line_prefix="   ", dense=True)
        print("")

        # DISK
        print(fmt.blue + fmt.underlined + "DISK" + fmt.reset + ":")
        print("")
        self.disk_output.show(line_prefix="   ", dense=True)
        print("")

        # OLD
        print(fmt.blue + fmt.underlined + "OLD" + fmt.reset + ":")
        print("")
        self.old_output.show(line_prefix="   ", dense=True)
        print("")

        # YOUNG
        print(fmt.blue + fmt.underlined + "YOUNG" + fmt.reset + ":")
        print("")
        self.young_output.show(line_prefix="   ", dense=True)
        print("")

        # SFR
        print(fmt.blue + fmt.underlined + "SFR" + fmt.reset + ":")
        print("")
        self.sfr_output.show(line_prefix="   ", dense=True)
        print("")

        # UNEVOLVED
        print(fmt.blue + fmt.underlined + "UNEVOLVED" + fmt.reset + ":")
        print("")
        self.unevolved_output.show(line_prefix="   ", dense=True)
        print("")

    # -----------------------------------------------------------------

    def show_data(self, **kwargs):

        """
        This function ...
        """

        # Debugging
        log.debug("Showing the available model data ...")

        # TOTAL
        print(fmt.blue + fmt.underlined + "TOTAL" + fmt.reset + ":")
        print("")
        self.total_data.show(line_prefix="  ", check_valid=False, dense=True)
        print("")

        # BULGE
        print(fmt.blue + fmt.underlined + "BULGE" + fmt.reset + ":")
        print("")
        self.bulge_data.show(line_prefix="   ", check_valid=False, dense=True)
        print("")

        # DISK
        print(fmt.blue + fmt.underlined + "DISK" + fmt.reset + ":")
        print("")
        self.disk_data.show(line_prefix="   ", check_valid=False, dense=True)
        print("")

        # OLD
        print(fmt.blue + fmt.underlined + "OLD" + fmt.reset + ":")
        print("")
        self.old_data.show(line_prefix="   ", check_valid=False, dense=True)
        print("")

        # YOUNG
        print(fmt.blue + fmt.underlined + "YOUNG" + fmt.reset + ":")
        print("")
        self.young_data.show(line_prefix="   ", check_valid=False, dense=True)
        print("")

        # SFR
        print(fmt.blue + fmt.underlined + "SFR" + fmt.reset + ":")
        print("")
        self.sfr_data.show(line_prefix="   ", check_valid=False, dense=True)
        print("")

        # UNEVOLVED
        print(fmt.blue + fmt.underlined + "UNEVOLVED" + fmt.reset + ":")
        print("")
        self.unevolved_data.show(line_prefix="   ", check_valid=False, dense=True)
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
        return load_modeling_environment(self.config.path)

    # -----------------------------------------------------------------

    @property
    def clipped_sed(self):
        return self.modeling_environment.observed_sed

    # -----------------------------------------------------------------

    @property
    def truncated_sed(self):
        return self.modeling_environment.truncated_sed

    # -----------------------------------------------------------------

    @property
    def asymptotic_sed(self):
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
        definition.add_optional("quantity", "string", "flux or luminosity", default_photometric_quantity_name, choices=photometric_quantity_names)
        definition.add_optional("spectral", "string", "spectral style", default_spectral_style, choices=spectral_style_names)
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

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot
        self.plot_total_sed(orientations=config.orientations, add_references=config.add_references,
                            additional_error=config.additional_error, unit=unit)

    # -----------------------------------------------------------------

    def plot_total_sed(self, orientations=default_orientations, add_references=False, additional_error=None, path=None,
                       show_file=False, title=None, format=default_plotting_format, unit=None):

        """
        This function ...
        :param orientations:
        :param add_references:
        :param additional_error:
        :param path:
        :param show_file:
        :param title:
        :param format:
        :param unit:
        :return:
        """

        # Debugging
        log.debug("Plotting total SED(s) ...")

        # Create SED plotter
        plotter = SEDPlotter()

        # Set unit
        if unit is not None: plotter.config.unit = unit
        plotter.config.distance = self.galaxy_distance

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
        definition.add_optional("quantity", "string", "flux or luminosity", default_photometric_quantity_name, choices=photometric_quantity_names)
        definition.add_optional("spectral", "string", "spectral style", default_spectral_style, choices=spectral_style_names)
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

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot
        self.plot_stellar_sed(config.observed_intrinsic, components=config.components, unit=unit)

    # -----------------------------------------------------------------

    def plot_stellar_sed(self, observed_intrinsic, components, path=None, title=None, show_file=False,
                         format=default_plotting_format, unit=None):

        """
        This function ...
        :param observed_intrinsic:
        :param components:
        :param path:
        :param title:
        :param show_file:
        :param format:
        :param unit:
        :return:
        """

        # Debugging
        log.debug("Plotting stellar SED(s) ...")

        # Create SED plotter
        plotter = SEDPlotter()

        # Set unit
        if unit is not None: plotter.config.unit = unit
        plotter.config.distance = self.galaxy_distance

        # Add references?
        # if add_references: plotter.add_seds(self.get_reference_seds(additional_error=additional_error))

        # Either observed of intrinsic
        if len(observed_intrinsic) == 1:

            oi = observed_intrinsic[0]

            # Loop over the components
            for component in components:

                # Set residuals flag
                residuals = component == total and oi == observed_name
                sed = self.model.get_stellar_sed(component, oi)
                name = component

                # Add SED to plotter
                plotter.add_sed(sed, name, residuals=residuals)

        # Both observed and intrinsic
        else:

            # ALLOW this for multiple components?

            # Loop over the components
            for component in components:
                for oi in observed_intrinsic:

                    # Set residuals flag
                    residuals = component == total and oi == observed_name

                    # Get the SED
                    sed = self.model.get_stellar_sed(component, oi)

                    # Add
                    plotter.add_sed(sed, component + " " + oi, residuals=residuals)

        # Set filepath, if plot is to be shown as file
        if path is None and show_file:
            if format is None: raise ValueError("Format has to be specified")
            path = fs.join(introspection.pts_temp_dir, "stellar_seds." + format)

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
        definition.add_optional("quantity", "string", "flux or luminosity", default_photometric_quantity_name, choices=photometric_quantity_names)
        definition.add_optional("spectral", "string", "spectral style", default_spectral_style, choices=spectral_style_names)
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

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot
        self.plot_dust_sed(config.components, unit=unit)

    # -----------------------------------------------------------------

    def plot_dust_sed(self, components, title=None, path=None, show_file=False, format=default_plotting_format, unit=None):

        """
        This function ...
        :param components:
        :param title:
        :param path:
        :param show_file:
        :param format:
        :param unit:
        :return:
        """

        # Debugging
        log.debug("Plotting dust SED(s) ...")

        # Create SED plotter
        plotter = SEDPlotter()

        # Set unit
        if unit is not None: plotter.config.unit = unit
        plotter.config.distance = self.galaxy_distance

        # Add references?
        # if add_references: plotter.add_seds(self.get_reference_seds(additional_error=additional_error))

        # Loop over the components
        for component in components:

            # Get the SED
            sed = self.model.get_dust_sed(component)

            # Add
            plotter.add_sed(sed, component, residuals=False)

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
        definition.add_optional("quantity", "string", "flux or luminosity", default_photometric_quantity_name, choices=photometric_quantity_names)
        definition.add_optional("spectral", "string", "spectral style", default_spectral_style, choices=spectral_style_names)
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
        contributions = config.pop("contributions")

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot
        self.plot_contribution_seds(contributions, unit=unit)

    # -----------------------------------------------------------------

    def get_sed_contribution(self, contribution, component=total):

        """
        This function ...
        :param contribution:
        :param component:
        :return:
        """

        # Get the simulations
        simulations = self.simulations[component]

        # Return the SED
        if contribution == total_contribution: return simulations.observed_sed
        elif contribution == direct_contribution:
            if simulations.has_full_sed: return simulations.observed_sed_direct
            else: return None
        elif contribution == scattered_contribution:
            if simulations.has_full_sed: return simulations.observed_sed_scattered
            else: return None
        elif contribution == dust_contribution:
            if simulations.has_full_sed: return simulations.observed_sed_dust
            else: return None
        elif contribution == dust_direct_contribution:
            if simulations.has_full_sed: return simulations.observed_sed_dust_direct
            else: return None
        elif contribution == dust_scattered_contribution:
            if simulations.has_full_sed: return simulations.observed_sed_dust_scattered
            else: return None
        elif contribution == transparent_contribution:
            if simulations.has_full_sed: return simulations.observed_sed_transparent
            else: return None
        else: raise ValueError("Invalid contribution: '" + contribution + "'")

    # -----------------------------------------------------------------

    def plot_contribution_seds(self, contributions, path=None, title=None, show_file=False, format=default_plotting_format,
                               component=total, unit=None):

        """
        This function ...
        :param contributions:
        :param path:
        :param title:
        :param show_file:
        :param format:
        :param component:
        :param unit:
        :return:
        """

        # Debugging
        log.debug("Plotting contribution SEDs ...")

        # Create SED plotter
        #plotter = SEDPlotter(kwargs) # **kwargs DOESN'T WORK? (e.g. with min_flux)
        plotter = SEDPlotter()

        # Set unit
        if unit is not None: plotter.config.unit = unit
        plotter.config.distance = self.galaxy_distance

        # Loop over the contributions
        for contribution in contributions:

            # Get the contribution SED
            sed = self.get_sed_contribution(contribution, component=component)
            if sed is None:
                log.warning("No '" + contribution + "' SED can be obtained for the '" + component + "' component: skipping ...")
                continue

            # Add
            residuals = contribution == total_contribution and component == total
            plotter.add_sed(sed, contribution, residuals=residuals)

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
        #definition.import_settings(plot_seds_definition)
        definition.add_optional("quantity", "string", "flux or luminosity", default_photometric_quantity_name, choices=photometric_quantity_names)
        definition.add_optional("spectral", "string", "spectral style", default_spectral_style, choices=spectral_style_names)
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

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot
        self.plot_component_seds(components, unit=unit)

    # -----------------------------------------------------------------

    def get_simulation_sed(self, component):

        """
        This function ...
        :param component:
        :return:
        """

        # Return the SED
        return self.simulations[component].observed_sed

    # -----------------------------------------------------------------

    def get_component_sed(self, component, dust_absorption=True, dust_emission=True):

        """
        This function ...
        :param component:
        :param dust_absorption:
        :param dust_emission:
        :return:
        """

        # Simulation SEDs
        if dust_emission:

            # Simulation SED
            if dust_absorption: return self.get_simulation_sed(component)

            # No dust absorption but dust emission, WEEEIRD
            else: raise NotImplementedError("This SED does not make physically sense")

        # No dust emission
        else:

            # With absorption
            if dust_absorption: return self.model.get_observed_stellar_sed(component)

            # No absorption: intrinsic
            else: return self.model.get_intrinsic_stellar_sed(component)

    # -----------------------------------------------------------------

    def get_observed_stellar_or_intrinsic_sed(self, component, observed_stellar_intrinsic):

        """
        This function ...
        :param component:
        :param observed_stellar_intrinsic:
        :return:
        """

        # Set flags
        if observed_stellar_intrinsic == observed_name: dust_absorption = dust_emission = True
        elif observed_stellar_intrinsic == intrinsic_name: dust_absorption = dust_emission = False
        elif observed_stellar_intrinsic == stellar_name:
            dust_absorption = True
            dust_emission = False
        else: raise ValueError("Invalid option for 'observed_stellar_or_intrinsic'")

        # Return
        return self.get_component_sed(component, dust_absorption=dust_absorption, dust_emission=dust_emission)

    # -----------------------------------------------------------------

    def plot_component_seds(self, components, path=None, title=None, show_file=False, format=default_plotting_format, unit=None):

        """
        This function ...
        :param components:
        :param path:
        :param title:
        :param show_file:
        :param format:
        :param unit:
        :return:
        """

        # Debugging
        log.debug("Plotting component SEDs ...")

        # Create SED plotter
        #plotter = SEDPlotter(kwargs)
        plotter = SEDPlotter()

        # Set unit
        if unit is not None: plotter.config.unit = unit
        plotter.config.distance = self.galaxy_distance

        # Add references?
        #if add_references: plotter.add_seds(self.get_reference_seds(additional_error=additional_error))

        # Loop over the components
        for component in components:

            # Get the SED
            sed = self.get_component_sed(component)

            # Add to plot
            residuals = component == total
            plotter.add_sed(sed, component, residuals=residuals)

        # Set filepath, if plot is to be shown as file
        if path is None and show_file:
            if format is None: raise ValueError("Format has to be specified")
            path = fs.join(introspection.pts_temp_dir, "component_seds." + format)

        # Run the plotter
        plotter.run(title=title, output=path)

        # Show file
        if show_file: fs.open_file(path)

    # -----------------------------------------------------------------

    def plot_component_sed(self, component, observed_stellar_intrinsic, unit=None):

        """
        This function ...
        :param component:
        :param observed_stellar_intrinsic:
        :param unit:
        :return:
        """

        # Either observed or intrinsic
        if len(observed_stellar_intrinsic) == 1:

            # Get the SED
            sed = self.get_observed_stellar_or_intrinsic_sed(component, observed_stellar_intrinsic[0])

            # Plot
            plot_sed(sed, unit=unit, distance=self.galaxy_distance)

        # Both
        else:

            seds = OrderedDict()
            for osi in observed_stellar_intrinsic: seds[osi] = self.get_observed_stellar_or_intrinsic_sed(component, osi)
            plot_seds(seds, residuals=False, unit=unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_old_bulge_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_stellar_intrinsic", "string_tuple", "plot observed SED (simulation), observed SED (stellar), intrinsic SED (stellar), or multiple", default_observed_stellar_intrinsic, choices=observed_stellar_intrinsic_choices)
        definition.add_optional("quantity", "string", "flux or luminosity", default_photometric_quantity_name, choices=photometric_quantity_names)
        definition.add_optional("spectral", "string", "spectral style", default_spectral_style, choices=spectral_style_names)
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

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot component
        self.plot_component_sed(bulge, config.observed_stellar_intrinsic, unit=unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_old_disk_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_stellar_intrinsic", "string_tuple", "plot observed SED (simulation), observed SED (stellar), intrinsic SED (stellar), or multiple", default_observed_stellar_intrinsic, choices=observed_stellar_intrinsic_choices)
        definition.add_optional("quantity", "string", "flux or luminosity", default_photometric_quantity_name, choices=photometric_quantity_names)
        definition.add_optional("spectral", "string", "spectral style", default_spectral_style, choices=spectral_style_names)
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

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot component
        self.plot_component_sed(disk, config.observed_stellar_intrinsic, unit=unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_old_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_stellar_intrinsic", "string_tuple", "plot observed SED (simulation), observed SED (stellar), intrinsic SED (stellar), or multiple", default_observed_stellar_intrinsic, choices=observed_stellar_intrinsic_choices)
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

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot component
        self.plot_component_sed(old, config.observed_stellar_intrinsic, unit=unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_young_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_stellar_intrinsic", "string_tuple", "plot observed SED (simulation), observed SED (stellar), intrinsic SED (stellar), or multiple", default_observed_stellar_intrinsic, choices=observed_stellar_intrinsic_choices)
        definition.add_optional("quantity", "string", "flux or luminosity", default_photometric_quantity_name, choices=photometric_quantity_names)
        definition.add_optional("spectral", "string", "spectral style", default_spectral_style, choices=spectral_style_names)
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

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot
        self.plot_component_sed(young, config.observed_stellar_intrinsic, unit=unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_sfr_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_stellar_intrinsic", "string_tuple", "plot observed SED (simulation), observed SED (stellar), intrinsic SED (stellar), or multiple", default_observed_stellar_intrinsic, choices=observed_stellar_intrinsic_choices)
        definition.add_optional("quantity", "string", "flux or luminosity", default_photometric_quantity_name, choices=photometric_quantity_names)
        definition.add_optional("spectral", "string", "spectral style", default_spectral_style, choices=spectral_style_names)
        return definition

    # -----------------------------------------------------------------

    def plot_sfr_sed_command(self, command, **kwargs):

        """
        This function ...
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.plot_sfr_sed_definition, **kwargs)

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot
        self.plot_component_sed(sfr, config.observed_stellar_intrinsic, unit=unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_sfr_intrinsic_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_optional("quantity", "string", "flux or luminosity", default_photometric_quantity_name, choices=photometric_quantity_names)
        definition.add_optional("spectral", "string", "spectral style", default_spectral_style, choices=spectral_style_names)
        return definition

    # -----------------------------------------------------------------

    def plot_sfr_intrinsic_sed_command(self, command, **kwargs):

        """
        This functino ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.plot_sfr_intrinsic_sed_definition, **kwargs)

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot
        self.plot_sfr_intrinsic_sed(unit=unit)

    # -----------------------------------------------------------------

    def plot_sfr_intrinsic_sed(self, unit=None):

        """
        This function ...
        :param unit:
        :return:
        """

        # Get stellar SEDs
        observed_stellar = self.model.get_stellar_sed(sfr, observed_name)
        intrinsic_stellar = self.model.get_stellar_sed(sfr, intrinsic_name)

        # Get intrinsic SEDs
        transparent_stellar = self.model.intrinsic_transparent_sfr_stellar_sed
        dust = self.model.intrinsic_sfr_dust_sed

        # Plot
        seds = OrderedDict()
        seds["observed stellar"] = observed_stellar
        seds["intrinsic"] = intrinsic_stellar
        seds["intrinsic (transparent) stellar"] = transparent_stellar
        seds["intrinsic dust"] = dust

        #print(observed_stellar)
        #print(intrinsic_stellar)
        #print(transparent_stellar)
        #print(dust)

        # Plot
        plot_seds(seds, unit=unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_unevolved_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("observed_stellar_intrinsic", "string_tuple", "plot observed SED (simulation), observed SED (stellar), intrinsic SED (stellar), or multiple", default_observed_stellar_intrinsic, choices=observed_stellar_intrinsic_choices)
        definition.add_optional("quantity", "string", "flux or luminosity", default_photometric_quantity_name, choices=photometric_quantity_names)
        definition.add_optional("spectral", "string", "spectral style", default_spectral_style, choices=spectral_style_names)
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

        # Get photometric unit
        unit = self.photometric_units[config.quantity][config.spectral]

        # Plot
        self.plot_component_sed(unevolved, config.observed_stellar_intrinsic, unit=unit)

    # -----------------------------------------------------------------

    @property
    def heating_path(self):
        return self.analysis_run.heating_path

    # -----------------------------------------------------------------

    @lazyproperty
    def spectral_heating_path(self):
        return fs.join(self.analysis_run.heating_path, "spectral")

    # -----------------------------------------------------------------

    @lazyproperty
    def spectral_heating_cells_path(self):
        return fs.join(self.spectral_heating_path, "3D")

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_cells_sed_filepath(self):
        total_filename = "total_curve_absorption.dat"
        return fs.get_filepath(self.spectral_heating_cells_path, total_filename, error_message="total spectral absorption SED file is not present: run spectral heating analysis first")

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_cells_sed(self):
        return SED.from_file(self.total_absorption_cells_sed_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_cells_sed_filepath(self):
        unevolved_filename = "unevolved_curve_absorption.dat"
        return fs.get_filepath(self.spectral_heating_cells_path, unevolved_filename, error_message="unevolved spectral absorption SED file is not present: run spectral heating analysis first")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_cells_sed(self):
        return SED.from_file(self.unevolved_absorption_cells_sed_filepath)

    # -----------------------------------------------------------------

    def plot_seds(self, **kwargs):
        plot_seds(kwargs, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_absorption_sed_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("component", "string", "component", total, choices=components)
        return definition

    # -----------------------------------------------------------------

    def plot_absorption_sed_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_absorption_sed_definition, **kwargs)

        # Total?
        if config.component == total: self.plot_absorption_sed_total()

        # Bulge
        elif config.component == bulge: self.plot_absorption_sed_bulge()

        # Disk
        elif config.component == disk: self.plot_absorption_sed_disk()

        # Old
        elif config.component == old: self.plot_absorption_sed_old()

        # Young
        elif config.component == young: self.plot_absorption_sed_young()

        # SFR
        elif config.component == sfr: self.plot_absorption_sed_sfr()

        # Unevolved
        elif config.component == unevolved: self.plot_absorption_sed_unevolved()

        # Invalid
        else: raise ValueError("Invalid component '" + config.component + "'")

    # -----------------------------------------------------------------

    def plot_absorption_sed_total(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting absorption for the total model ...")

    # -----------------------------------------------------------------

    def plot_absorption_sed_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting absorption for old bulge component ...")

    # -----------------------------------------------------------------

    def plot_absorption_sed_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting absorption for old disk component ...")

    # -----------------------------------------------------------------

    def plot_absorption_sed_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting absorption for old stars ...")

    # -----------------------------------------------------------------

    def plot_absorption_sed_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting absorption for young stars ...")

    # -----------------------------------------------------------------

    def plot_absorption_sed_sfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting absorption for star formation regions ...")

    # -----------------------------------------------------------------

    def plot_absorption_sed_unevolved(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Plotting absorption for unevolved stars ...")

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

    def get_component_attenuation_curve(self, component):

        """
        This function ...
        :param component:
        :return:
        """

        # Return
        return self.simulations[component].attenuation_curve

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

            # Get curve
            curve = self.get_component_attenuation_curve(component)

            # Add
            curves[component] = curve

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
    def plot_residuals_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Filters
        definition.add_optional("filters", "lazy_broad_band_filter_list", "filters for which to plot images", default="FUV,NUV,I1,MIPS 24mu,Pacs160,SPIRE350", convert_default=True)

        # Save to path
        definition.add_optional("path", "new_path", "save plot to file")

        # Dark mode
        definition.add_flag("dark", "plot in dark mode")

        # Other options
        definition.add_optional("zoom", "positive_real", "zoom from the normal galaxy truncation", 0.7)
        definition.add_optional("scale_xy_ratio", "positive_real", "scale the xy ratio to make plot panes more or less square", 1.)
        definition.add_optional("scale_xy_exponent", "positive_real", "exponent for the xy ratio to make plot panes more or less square", 0.7)
        definition.add_flag("mask", "mask the model image pixels that are invalid in the observed images", True)

        # Return
        return definition

    # -----------------------------------------------------------------

    def plot_residuals_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.plot_residuals_definition, **kwargs)

        # Plot residuals
        self.plot_residuals(config.filters, path=config.path, dark=config.dark, zoom=config.zoom,
                            scale_xy_ratio=config.scale_xy_ratio, scale_xy_exponent=config.scale_xy_exponent,
                            mask_simulated=config.mask)

    # -----------------------------------------------------------------

    def plot_residuals(self, filters, path=None, dark=False, zoom=1., scale_xy_ratio=1., scale_xy_exponent=1., mask_simulated=False):

        """
        Thisn function ...
        :param filters:
        :param path:
        :param dark:
        :param zoom:
        :param scale_xy_ratio:
        :param scale_xy_exponent:
        :param mask_simulated:
        :return:
        """

        from pts.magic.plot.imagegrid import plot_residuals_aplpy

        # Get images
        observations = self.get_observed_images(filters)
        models = self.get_simulated_images(filters)
        residuals = self.get_residual_images(filters)
        distributions = self.get_residual_distributions(filters)

        # Get center and radius
        center = self.galaxy_center
        radius = self.truncation_radius * zoom
        xy_ratio = (self.truncation_box_axial_ratio * scale_xy_ratio)**scale_xy_exponent

        #print(xy_ratio)

        # Plot
        plot_residuals_aplpy(observations, models, residuals, center=center, radius=radius, filepath=path, dark=dark,
                             xy_ratio=xy_ratio, distance=self.galaxy_distance, mask_simulated=mask_simulated)

    # -----------------------------------------------------------------
    # OBSERVED
    # -----------------------------------------------------------------

    def get_observed_images(self, filters):
        return self.static_photometry_dataset.get_frames_for_filters(filters)

    # -----------------------------------------------------------------
    # SIMULATED (MOCK)
    # -----------------------------------------------------------------

    @property
    def simulated_images_path(self):
        return self.analysis_run.images_path

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_images_dataset(self):
        return StaticDataSet.from_directory(self.simulated_images_path)

    # -----------------------------------------------------------------

    def get_simulated_images(self, filters):
        return self.simulated_images_dataset.get_frames_for_filters(filters)

    # -----------------------------------------------------------------
    # RESIDUALS
    # -----------------------------------------------------------------

    @property
    def residual_images_path(self):
        return fs.join(self.analysis_run.residuals_path, "maps")

    # -----------------------------------------------------------------

    @lazyproperty
    def residual_images_dataset(self):
        return StaticDataSet.from_directory(self.residual_images_path)

    # -----------------------------------------------------------------

    def get_residual_images(self, filters):
        return self.residual_images_dataset.get_frames_for_filters(filters)

    # -----------------------------------------------------------------
    # DISTRIBUTIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def residual_distributions_path(self):
        return fs.join(self.analysis_run.residuals_path, "distributions")

    # -----------------------------------------------------------------

    def get_residual_distributions(self, filters):
        distributions = []
        for fltr in filters:
            filepath = fs.join(self.residual_distributions_path, str(fltr) + ".dat")
            if not fs.is_file(filepath): distribution = None
            else: distribution = Distribution.from_file(filepath)
            distributions.append(distribution)
        return distributions

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
        definition.add_optional("filters", "lazy_broad_band_filter_list", "filters for which to plot images", default="FUV,NUV,I1,MIPS 24mu,Pacs160,SPIRE350", convert_default=True)
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
        return self.model.total_bolometric_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def faceon_cube(self):
        return self.model.total_bolometric_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def edgeon_cube(self):
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
        plot_datacube(datacube, title=instr_name, share_normalization=share_normalization, show_axes=False)

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
    # TOTAL SIMULATION
    # -----------------------------------------------------------------

    @property
    def total_simulations(self):
        return self.model.total_simulations

    # -----------------------------------------------------------------

    @property
    def total_simulation(self):
        return self.model.total_simulation

    # -----------------------------------------------------------------

    @property
    def total_output(self):
        return self.model.total_simulation_output

    # -----------------------------------------------------------------

    @property
    def total_data(self):
        return self.model.total_simulation_data

    # -----------------------------------------------------------------
    # BULGE SIMULATION
    # -----------------------------------------------------------------

    @property
    def bulge_simulations(self):
        return self.model.bulge_simulations

    # -----------------------------------------------------------------

    @property
    def bulge_simulation(self):
        return self.model.bulge_simulation

    # -----------------------------------------------------------------

    @property
    def bulge_output(self):
        return self.model.bulge_simulation_output

    # -----------------------------------------------------------------

    @property
    def bulge_data(self):
        return self.model.bulge_simulation_data

    # -----------------------------------------------------------------
    # DISK SIMULATION
    # -----------------------------------------------------------------

    @property
    def disk_simulations(self):
        return self.model.disk_simulations

    # -----------------------------------------------------------------

    @property
    def disk_simulation(self):
        return self.model.disk_simulation

    # -----------------------------------------------------------------

    @property
    def disk_output(self):
        return self.model.disk_simulation_output

    # -----------------------------------------------------------------

    @property
    def disk_data(self):
        return self.model.disk_simulation_data

    # -----------------------------------------------------------------
    # OLD SIMULATION
    # -----------------------------------------------------------------

    @property
    def old_simulations(self):
        return self.model.old_simulations

    # -----------------------------------------------------------------

    @property
    def old_simulation(self):
        return self.model.old_simulation

    # -----------------------------------------------------------------

    @property
    def old_output(self):
        return self.model.old_simulation_output

    # -----------------------------------------------------------------

    @property
    def old_data(self):
        return self.model.old_simulation_data

    # -----------------------------------------------------------------
    # YOUNG SIMULATION
    # -----------------------------------------------------------------

    @property
    def young_simulations(self):
        return self.model.young_simulations

    # -----------------------------------------------------------------

    @property
    def young_simulation(self):
        return self.model.young_simulation

    # -----------------------------------------------------------------

    @property
    def young_output(self):
        return self.model.young_simulation_output

    # -----------------------------------------------------------------

    @property
    def young_data(self):
        return self.model.young_simulation_data

    # -----------------------------------------------------------------
    # SFR SIMULATION
    # -----------------------------------------------------------------

    @property
    def sfr_simulations(self):
        return self.model.sfr_simulations

    # -----------------------------------------------------------------

    @property
    def sfr_simulation(self):
        return self.model.sfr_simulation

    # -----------------------------------------------------------------

    @property
    def sfr_output(self):
        return self.model.sfr_simulation_output

    # -----------------------------------------------------------------

    @property
    def sfr_data(self):
        return self.model.sfr_simulation_data

    # -----------------------------------------------------------------
    # UNEVOLVED SIMULATION
    # -----------------------------------------------------------------

    @property
    def unevolved_simulations(self):
        return self.model.unevolved_simulations

    # -----------------------------------------------------------------

    @property
    def unevolved_simulation(self):
        return self.model.unevolved_simulation

    # -----------------------------------------------------------------

    @property
    def unevolved_output(self):
        return self.model.unevolved_simulation_output

    # -----------------------------------------------------------------

    @property
    def unevolved_data(self):
        return self.model.unevolved_simulation_data

    # -----------------------------------------------------------------
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

        # Inform the user
        log.info("Analysing the model properties ...")

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
    def analyse_absorption_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Add settings
        definition.import_settings(analyse_absorption_definition)
        definition.remove_setting("run")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analyse_absorption_command(self, command, **kwargs):
        
        """
        This function ...
        :param command: 
        :param kwargs: 
        :return: 
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_absorption_definition, **kwargs)

        # Analyse
        self.analyse_absorption(config=config)
        
    # -----------------------------------------------------------------

    def analyse_absorption(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Create the analyser
        analyser = AbsorptionAnalyser(config=config)

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

        # Inform the user
        log.info("Analysing the dust cell heating ...")

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

        # Inform the user
        log.info("Analysing the projected heating ...")

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
    def analyse_spectral_heating_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Change settings
        definition.import_settings(analyse_spectral_heating_definition)
        definition.remove_setting("run")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analyse_spectral_heating_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_spectral_heating_definition, **kwargs)

        # Analyse
        self.analyse_spectral_heating(config=config)

    # -----------------------------------------------------------------

    def analyse_spectral_heating(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Inform the user
        log.info("Analysing the spectral heating ...")

        # Create the analyser
        analyser = SpectralDustHeatingAnalyser(config=config)

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

        # Inform the user
        log.info("Analysing the cell energy balance ...")

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

        # Inform the user
        log.info("Analysing the projected energy balance ...")

        # Create the analyser
        analyser = ProjectedEnergyAnalyser(config=config)

        # Set the modeling path
        analyser.config.path = self.config.path

        # Set the analysis run
        analyser.config.run = self.config.run

        # Run
        analyser.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_sfr_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Add settings
        definition.import_settings(analyse_sfr_definition)
        definition.remove_setting("run")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analyse_sfr_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_sfr_definition, **kwargs)

        # Analyse
        self.analyse_sfr(config=config)

    # -----------------------------------------------------------------

    def analyse_sfr(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Inform the user
        log.info("Analysing the star formation rates ...")

        # Create the analyser
        analyser = SFRAnalyser(config=config)

        # Set the modeling path
        analyser.config.path = self.config.path

        # Set the analysis run
        analyser.config.run = self.config.run

        # Run
        analyser.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_correlations_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Add settings
        definition.import_settings(analyse_correlations_definition)
        definition.remove_setting("run")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analyse_correlations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_correlations_definition, **kwargs)

        # Analyse
        self.analyse_correlations(config=config)

    # -----------------------------------------------------------------

    def analyse_correlations(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Inform the user
        log.info("Analysing the correlations ...")

        # Create the analyser
        analyser = CorrelationsAnalyser(config=config)

        # Set the modeling path
        analyser.config.path = self.config.path

        # Set the analysis run
        analyser.config.run = self.config.run

        # Run
        analyser.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_images_definition(self):

        """
        This fnction ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Add settings
        definition.import_settings(analyse_images_definition)
        definition.remove_setting("run")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analyse_images_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_images_definition, **kwargs)

        # Analyse
        self.analyse_images(config=config)

    # -----------------------------------------------------------------

    def analyse_images(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Inform the user
        log.info("Analysing the mock images ...")

        # Create the analyser
        analyser = ImagesAnalyser(config=config)

        # Set the modeling path
        analyser.config.path = self.config.path

        # Set the analysis run
        analyser.config.run = self.config.run

        # Run
        analyser.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_residuals_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Add settings
        definition.import_settings(analyse_residuals_definition)
        definition.remove_setting("run")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analyse_residuals_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_residuals_definition, **kwargs)

        # Analyse
        self.analyse_residuals(config=config)

    # -----------------------------------------------------------------

    def analyse_residuals(self, config=None):

        """
        Thi function ...
        :param config:
        :return:
        """

        # Inform the user
        log.info("Analysing the image residuals ...")

        # Create the analyser
        analyser = ResidualAnalyser(config=config)

        # Set the modeling path
        analyser.config.path = self.config.path

        # Set the analysis run
        analyser.config.run = self.config.run

        # Run
        analyser.run()

    # -----------------------------------------------------------------

    def examine_model(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Initialize
        examination = ModelExamination()

        # Run
        examination.run(model=self.model)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    # FROM ANALYSISPLOTTER:

    # def load_wavelength_grid(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Loading the wavelength grid ...")
    #
    #     # Determine the path to the wavelength grid file
    #     path = fs.join(self.analysis_path, "in", "wavelengths.txt")
    #
    #     # Load the wavelength grid
    #     if fs.is_file(path): self.wavelength_grid = WavelengthGrid.from_skirt_input(path)
    #
    # # -----------------------------------------------------------------
    #
    # def load_transmission_curves(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Loading the transmission curves ...")
    #
    #     # Load the observed SED
    #     sed = ObservedSED.from_file(self.observed_sed_path)
    #
    #     # Loop over all filters for the points in the SED
    #     for fltr in sed.filters():
    #
    #         # Create the transmission curve
    #         transmission = TransmissionCurve.from_filter(fltr)
    #
    #         # Normalize the transmission curve
    #         transmission.normalize(value=1.0, method="max")
    #
    #         # Add the transmission curve to the dictionary
    #         self.transmission_curves[str(fltr)] = transmission
    #
    # # -----------------------------------------------------------------
    #
    # def plot_wavelengths(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Plotting the wavelength grid ...")
    #
    #     # Create the transmission plotter
    #     plotter = TransmissionPlotter()
    #
    #     plotter.title = "Wavelengths used for analysis"
    #     plotter.transparent = True
    #
    #     # Add the transmission curves
    #     for label in self.transmission_curves: plotter.add_transmission_curve(self.transmission_curves[label], label)
    #
    #     # Add the wavelength points
    #     for wavelength in self.wavelength_grid.wavelengths(): plotter.add_wavelength(wavelength)
    #
    #     # Determine the path to the plot file
    #     path = fs.join(self.plot_analysis_path, "wavelengths.pdf")
    #
    #     # Run the plotter
    #     plotter.run(path, min_wavelength=self.wavelength_grid.min_wavelength, max_wavelength=self.wavelength_grid.max_wavelength, min_transmission=0.0, max_transmission=1.05)

    # -----------------------------------------------------------------

    @property
    def history_filename(self):

        """
        This function ...
        :return:
        """

        return "analysis"

# -----------------------------------------------------------------
