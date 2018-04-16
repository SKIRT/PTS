#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.model Contains the RTModel class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..core.mappings import Mappings
from ...core.filter.filter import parse_filter
from ...core.tools.utils import lazyproperty
from ...core.tools import sequences
from ..config.parameters import distance_name, ionizing_scaleheight_name, sfr_compactness_name, fuv_young_name
from ..config.parameters import old_scaleheight_name, position_angle_name, dust_mass_name, fuv_ionizing_name
from ..config.parameters import metallicity_name, young_scaleheight_name, sfr_covering_name, dust_scaleheight_name
from ..config.parameters import i1_old_name, sfr_pressure_name, inclination_name
from ..config.parameters import modeling_parameter_labels
from ...core.basics.containers import create_subdict
from ...magic.tools import extinction
from ..basics.instruments import SEDInstrument
from ...core.tools import filesystem as fs
from ...core.simulation.skifile import SkiFile
from ...core.prep.smile import get_oligochromatic_template, get_panchromatic_template
from ..build.construct import add_new_stellar_component
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.simulation.simulation import createsimulations
from ...core.simulation.execute import run_simulation
from ...core.basics.log import log
from ...core.data.sed import SED
from ...core.simulation.output import SimulationOutput
from ...core.simulation.data import SimulationData
from ...core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

stellar_dust_sed_split_wavelength = 5. * u("micron")

# -----------------------------------------------------------------

# Names of derived model properties

## Total
obs_total_bol_lum_name = "Observed total bolometric luminosity" # 1
intr_total_bol_lum_name = "Intrinsic total bolometric luminosity" # 2 ; should be same as 1
obs_total_stellar_bol_lum_name = "Observed total stellar bolometric luminosity" # 3
intr_total_stellar_bol_lum_name = "Intrinsic total stellar bolometric luminosity" # 4 ; should be same as 2
bol_attenuation_name = "Total bolometric attenuation" # 5
total_dust_mass_name = "Total dust mass" # 6 with SFR dust mass

## Old bulge
obs_bulge_spec_lum_name = "Observed old stellar bulge specific luminosity" # 7
intr_bulge_spec_lum_name = "Intrinsic old stellar bulge specific luminosity" # 8 part of parameter set
obs_bulge_bol_lum_name = "Observed old stellar bulge bolometric luminosity" # 9
intr_bulge_bol_lum_name = "Intrinsic old stellar bulge bolometric luminosity" # 10
bulge_spec_attenuation_name = "Old stellar bulge specific attenuation" # 11
bulge_bol_attenuation_name = "Old stellar bulge bolometric attenuation" # 12
obs_bulge_total_lum_name = "Old stellar bulge bolometric total luminosity" # 13 stellar + dust ; should be the same as 10 (intrinsic lum)
obs_bulge_dust_lum_name = "Old stellar bulge bolometric dust luminosity" # 14

## Old disk
obs_disk_spec_lum_name = "Observed old stellar disk specific luminosity" #
intr_disk_spec_lum_name = "Intrinsic old stellar disk specific luminosity" # part of parameter set
obs_disk_bol_lum_name = "Observed old stellar disk bolometric luminosity" #
intr_disk_bol_lum_name = "Intrinsic old stellar disk bolometric luminosity" #
disk_spec_attenuation_name = "Old stellar disk specific attenuation" #
disk_bol_attenuation_name = "Old stellar disk bolometric attenuation" #
obs_disk_total_lum_name = "Old stellar disk bolometric total luminosity" # stellar + dust ; should be the same as x (intrinsic lum)
obs_disk_dust_lum_name = "Old stellar disk bolometric dust luminosity" #

## Old (evolved)
obs_old_spec_lum_name = "Observed old stellar specific luminosity" #
intr_old_spec_lum_name = "Intrinsic old stellar specific luminosity" #
obs_old_bol_lum_name = "Observed old stellar bolometric luminosity" #
intr_old_bol_lum_name = "Intrinsic old stellar bolometric luminosity" #
old_spec_attenuation_name = "Old stellar specific attenuation" #
old_bol_attenuation_name = "Old stellar bolometric attenuation" #
obs_old_total_lum_name = "Old stellar bolometric total luminosity" # stellar + dust ; should be the same as x (intrinsic lum)
obs_old_dust_lum_name = "Old stellar bolometric dust luminosity" #

## Young stars
obs_young_spec_lum_name = "Observed young stellar specific luminosity" #
intr_young_spec_lum_name = "Intrinsic young stellar specific luminosity" # part of (free) parameter set
obs_young_bol_lum_name = "Observed young stellar bolometric luminosity" #
intr_young_bol_lum_name = "Intrinsic young stellar bolometric luminosity" #
young_spec_attenuation_name = "Young stellar specific attenuation" #
young_bol_attenuation_name = "Young stellar bolometric attenuation" #
obs_young_total_lum_name = "Young stellar bolometric total luminosity" # stellar + dust ; should be the same as x (intrinsic lum)
obs_young_dust_lum_name = "Young stellar bolometric dust luminosity" #

## Ionizing stars (SFR)
sfr_name = "Star formation rate"
obs_sfr_spec_lum_name = "Observed SFR specific luminosity" #
intr_sfr_spec_lum_name = "Intrinsic SFR specific luminosity" #
obs_sfr_bol_lum_name = "Observed SFR bolometric luminosity" #
intr_sfr_bol_lum_name = "Intrinsic SFR bolometric luminosity" #
sfr_spec_attenuation_name = "SFR specific attenuation" #
sfr_bol_attenuation_name = "SFR bolometric attenuation" #
obs_sfr_stellar_bol_lum_name = "Observed SFR stellar bolometric luminosity" #
intr_sfr_stellar_bol_lum_name = "Intrinsic SFR stellar bolometric luminosity" #
sfr_stellar_mass_name = "SFR stellar mass" # internal stars
sfr_dust_mass_name = "SFR internal dust mass" # INTERNAL SO ONLY THE DUST IN MAPPINGS
sfr_dust_lum_name = "SFR internal bolometric dust luminosity" # INTERNAL SO ONLY THE DUST IN MAPPINGS
obs_sfr_total_lum_name = "SFR bolometric total luminosity" # stellar(SFR with MAPPINGS) + dust ; should be the same as x (intrinsic lum)
obs_sfr_dust_lum_name = "SFR bolometric dust luminosity" # INTERNAL DUST + DUST DISK

## Young + ionizing (unevolved)
obs_unevolved_spec_lum_name = "Observed unevolved stellar specific luminosity" #
intr_unevolved_spec_lum_name = "Intrinsic unevolved stellar specific luminosity" #
obs_unevolved_bol_lum_name = "Observed unevolved stellar bolometric luminosity" #
intr_unevolved_bol_lum_name = "Intrinsic unevolved stellar bolometric luminosity" #
unevolved_spec_attenuation_name = "Unevolved stellar specific attenuation" #
unevolved_bol_attenuation_name = "Unevolved stellar bolometric attenuation" #
obs_unevolved_total_lum_name = "Unevolved stellar bolometric total luminosity" # stellar + dust ; should be the same as x (intrinsic lum)
obs_unevolved_dust_lum_name = "Unevolved stellar bolometric dust luminosity" #

## Dust
dust_lum_name = "Bolometric dust luminosity" #
diffuse_dust_lum_name = "Bolometric diffuse dust luminosity"

# -----------------------------------------------------------------

sed_dirname = "sed"

# -----------------------------------------------------------------

old_bulge_name = "old_bulge"
old_disk_name = "old_disk"
young_name = "young"

# -----------------------------------------------------------------

default_npackages = 1e5
wavelengths_filename = "wavelengths.txt"
map_filename = "map.fits"

# Instruments/orientations
earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

# Contributions
total_contribution = "total"
direct_contribution = "direct"
scattered_contribution = "scattered"
dust_contribution = "dust"
dust_scattered_contribution = "dustscattered"
transparent_contribution = "transparent"
contributions = [total_contribution, direct_contribution, scattered_contribution, dust_contribution, dust_scattered_contribution, transparent_contribution]

# -----------------------------------------------------------------

class RTModel(object):
    
    """
    Objects of this class describe a full radiative transfer model.
    """

    def __init__(self, definition, wavelength_grid=None, simulation_name=None, chi_squared=None,
                 free_parameter_labels=None, free_parameter_values=None, observed_total_output_path=None,
                 observed_bulge_output_path=None, observed_disk_output_path=None, observed_old_output_path=None,
                 observed_young_output_path=None, observed_sfr_output_path=None, observed_unevolved_output_path=None,
                 parameters=None):

        """
        The constructor ...
        :param definition: model definition
        :param wavelength_grid: wavelength grid of the simulation
        :param simulation_name: name for the simulation
        :param chi_squared: reduced chi squared value
        :param free_parameter_labels: labels of free parameters
        :param free_parameter_values: values of free parameters (as dict)
        :param observed_total_output_path: output path of simulation with total stellar contribution + dust
        :param observed_bulge_output_path: output path of simulation with old stellar bulge contribution + dust
        :param observed_disk_output_path: output path of simulation with old stellar disk contribution + dust
        :param observed_old_output_path: output path of simulation with old stellar contribution + dust
        :param observed_young_output_path: output path of simulation with young stellar contribution + dust
        :param observed_sfr_output_path: output path of simulation with SFR contribution + dust
        :param observed_unevolved_output_path: output path of simulation with young stellar and SFR contribution + dust
        :param parameters:
        :return:
        """

        # Set wavelength grid
        self.wavelength_grid = wavelength_grid

        # Other attributes
        self.simulation_name = simulation_name
        self.chi_squared = chi_squared
        self.free_parameter_labels = free_parameter_labels if free_parameter_labels is not None else []

        # The model definition describing the components
        self.definition = definition

        # Set the free parameter values explicitely
        if free_parameter_values is not None:

            # Check
            if free_parameter_labels is not None: raise ValueError("Cannot pass both free parameter labels as free parameter values")
            self.free_parameter_labels = free_parameter_values.keys()

            # Set values
            self.free_parameter_values = free_parameter_values

            # Set a flag that indicates that the definition may not be containing the same parameter values as the explicitly set free parameter values
            self._unsynchronized = True

        # Parameters are determined directly from the definition components
        else: self._unsynchronized = False

        # Save the wavelength grid in SKIRT input file
        if self.has_wavelength_grid: self.wavelength_grid.to_skirt_input(self.wavelength_grid_path)

        # No wavelength grid passed, load the wavelength grid if necessary, and if present
        elif self.has_wavelengths_directory and fs.is_file(self.wavelength_grid_path): self.wavelength_grid = WavelengthGrid.from_skirt_input(self.wavelength_grid_path)

        # Simulation output paths
        self.observed_total_output_path = observed_total_output_path
        self.observed_bulge_output_path = observed_bulge_output_path
        self.observed_disk_output_path = observed_disk_output_path
        self.observed_old_output_path = observed_old_output_path
        self.observed_young_output_path = observed_young_output_path
        self.observed_sfr_output_path = observed_sfr_output_path
        self.observed_unevolved_output_path = observed_unevolved_output_path

        # Set parameters
        if parameters is not None: self.set_parameters(**parameters)

    # -----------------------------------------------------------------

    def set_parameters(self, **values):

        """
        This function ...
        :param values:
        :return:
        """

        # Loop over the values
        for name in values:

            # Set
            value = values[name]
            setattr(self, name, value)

    # -----------------------------------------------------------------

    @property
    def has_observed_total_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_output_path is not None and fs.is_directory(self.observed_total_output_path) and not fs.is_empty(self.observed_total_output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_total_output(self):

        """
        This function ...
        :return:
        """

        if not self.has_observed_total_output: raise ValueError("No output")
        return SimulationOutput.from_directory(self.observed_total_output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_total_data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.observed_total_output)

    # -----------------------------------------------------------------

    @property
    def has_observed_total_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_total_output and self.observed_total_data.has_seds

    # -----------------------------------------------------------------

    @property
    def observed_total_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_data.sed_paths_instruments[earth_name]

    # -----------------------------------------------------------------

    @property
    def has_other_observed_total_sed_orientations(self):

        """
        This function ...
        :return:
        """

        return len(self.observed_total_data.sed_paths_instruments) > 1

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_total_sed_orientation(self):

        """
        This function ...
        :return:
        """

        return faceon_name in self.observed_total_data.sed_paths_instruments

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_total_sed_orientation(self):

        """
        This function ...
        :return:
        """

        return edgeon_name in self.observed_total_data.sed_paths_instruments

    # -----------------------------------------------------------------

    @property
    def faceon_observed_total_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_data.sed_paths_instruments[faceon_name]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_total_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_data.sed_paths_instruments[edgeon_name]

    # -----------------------------------------------------------------

    @property
    def observed_total_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_data.seds[earth_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_sed_total(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_sed

    # -----------------------------------------------------------------

    @property
    def faceon_observed_total_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_data.seds[faceon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_total_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_data.seds[edgeon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def has_other_observed_total_sed_contributions(self):

        """
        This function ...
        :return:
        """

        return len(self.observed_total_data.seds[earth_name]) > 1

    # -----------------------------------------------------------------

    @property
    def observed_total_sed_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_data.seds[earth_name][direct_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_total_sed_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_data.seds[earth_name][scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_total_sed_dust(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_data.seds[earth_name][dust_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_total_sed_dust_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_data.seds[earth_name][dust_scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_total_sed_transparent(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_data.seds[earth_name][transparent_contribution]

    # -----------------------------------------------------------------

    @property
    def intrinsic_total_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_sed_transparent

    # -----------------------------------------------------------------

    @property
    def faceon_intrinsic_total_sed(self):

        """
        This function ...
        :return:
        """

        raise RuntimeError("Contributions are not simulated for faceon instrument (no full instrument)?")

    # -----------------------------------------------------------------

    @property
    def edgeon_intrinsic_total_sed(self):

        """
        This function ...
        :return:
        """

        raise RuntimeError("Contributions are not simulated for edgeon instrument (no full instrument)?")

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_total_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_total_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_total_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_total_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_total_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_total_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_other_observed_total_sed_contributions

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_bulge_output_path is not None and fs.is_directory(self.observed_bulge_output_path) and not fs.is_empty(self.observed_bulge_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_old_bulge_output(self):

        """
        This function ...
        :return:
        """

        if not self.has_observed_old_bulge_output: raise ValueError("No output")
        return SimulationOutput.from_directory(self.observed_bulge_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_old_bulge_data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.observed_old_bulge_output)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_data(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_observed_old_bulge_output and self.observed_old_bulge_data.has_seds

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_disk_output_path is not None and fs.is_directory(self.observed_disk_output_path) and not fs.is_empty(self.observed_disk_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_old_disk_output(self):

        """
        This function ...
        :return:
        """

        if not self.has_observed_old_disk_output: raise ValueError("No output")
        return SimulationOutput.from_directory(self.observed_disk_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_old_disk_data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.observed_old_disk_output)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_data(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_old_disk_output and self.observed_old_disk_data.has_seds

    # -----------------------------------------------------------------

    @property
    def has_observed_old_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_old_output_path is not None and fs.is_directory(self.observed_old_output_path) and not fs.is_empty(self.observed_old_output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_output(self):

        """
        This function ...
        :return:
        """

        if not self.has_observed_old_output: raise ValueError("No output")
        return SimulationOutput.from_directory(self.observed_old_output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.observed_old_output)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_data(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_old_output and self.observed_old_data.has_seds

    # -----------------------------------------------------------------

    @property
    def has_observed_young_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_young_output_path is not None and fs.is_directory(self.observed_young_output_path) and not fs.is_empty(self.observed_young_output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_young_output(self):

        """
        This function ...
        :return:
        """

        if not self.has_observed_young_output: raise ValueError("No output")
        return SimulationOutput.from_directory(self.observed_young_output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_young_data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.observed_young_output)

    # -----------------------------------------------------------------

    @property
    def has_observed_young_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_young_output and self.observed_young_data.has_seds

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_sfr_output_path is not None and fs.is_directory(self.observed_sfr_output_path) and not fs.is_empty(self.observed_sfr_output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sfr_output(self):

        """
        This function ...
        :return:
        """

        if not self.has_observed_sfr_output: raise ValueError("No output")
        return SimulationOutput.from_directory(self.observed_sfr_output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sfr_data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.observed_sfr_output)

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sfr_output and self.observed_sfr_data.has_seds

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_unevolved_output_path is not None and fs.is_directory(self.observed_unevolved_output_path) and not fs.is_empty(self.observed_unevolved_output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_unevolved_output(self):

        """
        This function ...
        :return:
        """

        if not self.has_observed_unevolved_output: raise ValueError("No output")
        return SimulationOutput.from_directory(self.observed_unevolved_output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_unevolved_data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.observed_unevolved_output)

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_unevolved_output and self.observed_unevolved_data.has_seds

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_filter(self):

        """
        This function ...
        :return:
        """

        return parse_filter("IRAC I1")

    # -----------------------------------------------------------------

    @property
    def i1_wavelength(self):

        """
        This function ...
        :return:
        """

        return self.i1_filter.wavelength

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_filter(self):

        """
        This function ...
        :return:
        """

        return parse_filter("GALEX FUV")

    # -----------------------------------------------------------------

    @property
    def fuv_wavelength(self):

        """
        This function ...
        :return:
        """

        return self.fuv_filter.wavelength

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        return modeling_parameter_labels

    # -----------------------------------------------------------------

    @lazyproperty
    def other_parameter_labels(self):

        """
        This function ...
        :return:
        """

        return sequences.elements_not_in_other(self.parameter_labels, self.free_parameter_labels)

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_values(self):

        """
        This function ...
        :return:
        """

        from ..fitting.configuration import get_definition_value

        # Initialize dictionary
        values = OrderedDict()

        # Loop over all parameter labels
        for label in self.parameter_labels:

            # Get the value
            value = get_definition_value(self.definition, label)

            # Add the value
            values[label] = value

        # Return the dictionary of values
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def free_parameter_values(self):

        """
        This function ...
        :return:
        """

        return create_subdict(self.parameter_values, self.free_parameter_labels)

    # -----------------------------------------------------------------

    @lazyproperty
    def other_parameter_values(self):

        """
        This function ...
        :return:
        """

        return create_subdict(self.parameter_values, self.other_parameter_labels)

    # -----------------------------------------------------------------

    @property
    def metallicity(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values[metallicity_name]

    # -----------------------------------------------------------------

    @property
    def sfr_compactness(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values[sfr_compactness_name]

    # -----------------------------------------------------------------

    @property
    def sfr_pressure(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values[sfr_pressure_name]

    # -----------------------------------------------------------------

    @property
    def sfr_covering_factor(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values[sfr_covering_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr(self):

        """
        This function derives the SFR (in Msun / year) from the FUV luminosity of the model and the intrinsic MAPPINGS SED
        :return:
        """

        # Get the SFR
        return Mappings.sfr_for_luminosity(self.metallicity, self.sfr_compactness, self.sfr_pressure, self.sfr_covering_factor, self.intrinsic_fuv_luminosity_sfr, self.fuv_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_sfr(self):

        """
        This function ...
        :return:
        """

        return self.has_mappings

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings(self):

        """
        This function ...
        :return:
        """

        # Create the MAPPINGS template and return it
        return Mappings(self.metallicity, self.sfr_compactness, self.sfr_pressure, self.sfr_covering_factor, self.sfr)

    # -----------------------------------------------------------------

    @property
    def has_mappings(self):

        """
        This function ...
        :return:
        """

        return True # should always be able to be created

    # -----------------------------------------------------------------

    @lazyproperty
    def normalized_mappings(self):

        """
        This function ...
        :return:
        """

        # Create the MAPPINGS template
        return Mappings(self.metallicity, self.sfr_compactness, self.sfr_pressure, self.sfr_covering_factor)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the I1 wavelength
        return self.observed_sed_old_bulge.photometry_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_i1_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.definition.bulge_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_i1_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return True # should be defined in definition

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_attenuation_old_bulge(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_i1_luminosity_old_bulge.to("W/micron", wavelength=self.i1_wavelength, distance=self.distance).value
        intrinsic = self.intrinsic_i1_luminosity_old_bulge.to("W/micron", wavelength=self.i1_wavelength, distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_i1_attenuation_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_i1_luminosity_old_bulge and self.has_intrinsic_i1_luminosity_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bolometric_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.observed_old_bulge_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_old_bulge_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_bolometric_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed_old_bulge.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_attenuation_old_bulge(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_bolometric_luminosity_old_bulge.to("W", distance=self.distance).value
        intrinsic = self.intrinsic_bolometric_luminosity_old_bulge.to("W", distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_bolometric_luminosity_old_bulge and self.has_intrinsic_bolometric_luminosity_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_bulge_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_old_bulge.splice(x_max=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_bulge_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_old_bulge.splice(x_min=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_bulge_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_old_bulge.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_bulge_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_old_bulge_dust_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_old_bulge_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the I1 wavelength
        return self.observed_sed_old_disk.photometry_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_i1_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.definition.old_stars_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_i1_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return True # should be defined in definition

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_attenuation_old_disk(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_i1_luminosity_old_disk.to("W/micron", wavelength=self.i1_wavelength, distance=self.distance).value
        intrinsic = self.intrinsic_i1_luminosity_old_disk.to("W/micron", wavelength=self.i1_wavelength, distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_i1_attenuation_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_i1_luminosity_old_disk and self.has_intrinsic_i1_luminosity_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bolometric_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.observed_old_disk_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_old_disk_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_bolometric_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed_old_disk.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_attenuation_old_disk(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_bolometric_luminosity_old_disk.to("W", distance=self.distance).value
        intrinsic = self.intrinsic_bolometric_luminosity_old_disk.to("W", distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_bolometric_luminosity_old_disk and self.has_intrinsic_bolometric_luminosity_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_disk_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_old_disk.splice(x_max=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_disk_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_old_disk.splice(x_min=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_disk_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_old_disk.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_disk_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_old_disk_dust_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_old_disk_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_old(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the I1 wavelength
        return self.observed_sed_old.photometry_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_old

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_i1_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_i1_luminosity_old_bulge + self.intrinsic_i1_luminosity_old_disk

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_i1_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_i1_luminosity_old_bulge and self.has_intrinsic_i1_luminosity_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_attenuation_old(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_i1_luminosity_old.to("W/micron", wavelength=self.i1_wavelength, distance=self.distance).value
        intrinsic = self.intrinsic_i1_luminosity_old.to("W/micron", wavelength=self.i1_wavelength, distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_i1_attenuation_old(self):

        """
        This fnuction ...
        :return:
        """

        return self.has_observed_i1_luminosity_old and self.has_intrinsic_i1_luminosity_old

    # -----------------------------------------------------------------

    @property
    def attenuation_i1_old(self):

        """
        This function ...
        :return:
        """

        return self.i1_attenuation_old

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bolometric_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.observed_old_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_old_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_bolometric_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed_old.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_old(self):

        """
        This function ...
        :return:
        """
        
        return self.has_intrinsic_sed_old

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_attenuation_old(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_bolometric_luminosity_old.to("W", distance=self.distance).value
        intrinsic = self.intrinsic_bolometric_luminosity_old.to("W", distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_old(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_bolometric_luminosity_old and self.has_intrinsic_bolometric_luminosity_old

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_old.splice(x_max=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_old

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_old.splice(x_min=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_dust_sed(self):
        
        """
        This function ...
        :return: 
        """
        
        return self.has_observed_sed_old

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_old.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_old_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_old

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_old_dust_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_old_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_old_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_young(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the FUV wavelength
        return self.observed_sed_young.photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_young

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values[fuv_young_name]

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return True # part of free parameters

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_attenuation_young(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_fuv_luminosity_young.to("W/micron", wavelength=self.fuv_wavelength, distance=self.distance).value
        intrinsic = self.intrinsic_fuv_luminosity_young.to("W/micron", wavelength=self.fuv_wavelength, distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_fuv_attenuation_young(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_fuv_luminosity_young and self.has_intrinsic_fuv_luminosity_young

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bolometric_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.observed_young_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_young_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_bolometric_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed_young.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed_young

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_attenuation_young(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_bolometric_luminosity_young.to("W", distance=self.distance).value
        intrinsic = self.intrinsic_bolometric_luminosity_young.to("W", distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_young(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_bolometric_luminosity_young and self.has_intrinsic_bolometric_luminosity_young

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_young_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_young.splice(x_max=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_young_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_young

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_young_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_young.splice(x_min=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_young_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_young

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_young_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_young.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_young_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_young

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_young_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_young_dust_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_young_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_young_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the FUV wavelength
        return self.observed_sed_sfr.photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values[fuv_ionizing_name]

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return True # free parameter

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_attenuation_sfr(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_fuv_luminosity_sfr.to("W/micron", wavelength=self.fuv_wavelength, distance=self.distance).value
        intrinsic = self.intrinsic_fuv_luminosity_sfr.to("W/micron", wavelength=self.fuv_wavelength, distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_fuv_attenuation_sfr(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_fuv_luminosity_sfr and self.has_intrinsic_fuv_luminosity_sfr

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.observed_sfr_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sfr_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sfr_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_attenuation_sfr(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_bolometric_luminosity_sfr.to("W", distance=self.distance).value
        intrinsic = self.intrinsic_bolometric_luminosity_sfr.to("W", distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_sfr(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_bolometric_luminosity_sfr and self.has_intrinsic_bolometric_luminosity_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values[dust_mass_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.mappings.dust_mass

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.has_mappings

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.dust_mass + self.sfr_dust_mass

    # -----------------------------------------------------------------

    @property
    def has_total_dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr_dust_mass

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_stellar_mass(self):

        """
        This function ...
        :return:
        """

        return self.mappings.stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_sfr_stellar_mass(self):

        """
        This function ...
        :return:
        """

        #return self.has_mappings
        return False # returns NotImplementedError in Mappings: we don't know the conversion yet between Mappings parameters and stellar mass!

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the FUV wavelength
        return self.observed_sed_unevolved.photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_unevolved

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_fuv_luminosity_young + self.intrinsic_fuv_luminosity_sfr

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_fuv_luminosity_young and self.has_intrinsic_fuv_luminosity_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bolometric_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.observed_unevolved_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_unevolved_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_bolometric_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed_unevolved.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed_unevolved

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_attenuation_unevolved(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_fuv_luminosity_unevolved.to("W/micron", wavelength=self.fuv_wavelength, distance=self.distance).value
        intrinsic = self.intrinsic_fuv_luminosity_unevolved.to("W/micron", wavelength=self.fuv_wavelength, distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_fuv_attenuation_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_fuv_luminosity_unevolved and self.has_intrinsic_fuv_luminosity_unevolved

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_attenuation_unevolved(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_bolometric_luminosity_unevolved.to("W", distance=self.distance).value
        intrinsic = self.intrinsic_bolometric_luminosity_unevolved.to("W", distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_bolometric_luminosity_unevolved and self.has_intrinsic_bolometric_luminosity_unevolved

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_unevolved_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_unevolved.splice(x_max=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_unevolved

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_unevolved_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_unevolved.splice(x_min=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_unevolved

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_unevolved_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_unevolved.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_unevolved

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_unevolved_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_unevolved_dust_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_unevolved_dust_sed

    # -----------------------------------------------------------------

    @property
    def distance(self):

        """
        This function ...
        :return:
        """

        return self.definition.distance

    # -----------------------------------------------------------------

    @property
    def inclination(self):

        """
        This function ...
        :return:
        """

        return self.definition.inclination

    # -----------------------------------------------------------------

    @property
    def position_angle(self):

        """
        This function ...
        :return:
        """

        return self.definition.position_angle

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_instrument(self):

        """
        This function ...
        :return:
        """

        return SEDInstrument.from_properties(self.distance, self.inclination, self.position_angle)

    # -----------------------------------------------------------------

    @property
    def definition_path(self):

        """
        This function ...
        :return:
        """

        return self.definition.path

    # -----------------------------------------------------------------

    @property
    def old_bulge_path(self):

        """
        This function ...
        :return:
        """

        return self.definition.bulge_component_path

    # -----------------------------------------------------------------

    @property
    def old_disk_path(self):

        """
        This function ..
        :return:
        """

        return self.definition.old_stars_component_path

    # -----------------------------------------------------------------

    @property
    def young_path(self):

        """
        This function ...
        :return:
        """

        return self.definition.young_stars_component_path

    # -----------------------------------------------------------------

    @property
    def sfr_path(self):

        """
        This function ...
        :return:
        """

        return self.definition.ionizing_stars_component_path

    # -----------------------------------------------------------------

    @property
    def dust_path(self):

        """
        This function ...
        :return:
        """

        return self.definition.dust_component_path

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_sed_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.old_bulge_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_sed_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.old_disk_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_sed_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.young_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_sed_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.sfr_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_path(self):

        """
        Thisf unction ..
        :return:
        """

        return fs.create_directory_in(self.dust_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_ski_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.old_bulge_sed_path, old_bulge_name + ".ski")

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_ski_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.old_disk_sed_path, old_disk_name + ".ski")

    # -----------------------------------------------------------------

    @lazyproperty
    def young_ski_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.young_sed_path, young_name + ".ski")

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_ski_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.sfr_sed_path, sfr_name + ".ski")

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_skifile(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_bulge_ski_path)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_skifile(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_disk_ski_path)

    # -----------------------------------------------------------------

    @property
    def has_young_skifile(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_ski_path)

    # -----------------------------------------------------------------

    @property
    def has_sfr_skifile(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_ski_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_component(self):

        """
        This function ...
        :return:
        """

        return self.definition.load_bulge_component()

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_component(self):

        """
        This function ...
        :return:
        """

        return self.definition.load_old_stars_component()

    # -----------------------------------------------------------------

    @lazyproperty
    def young_component(self):

        """
        This function ...
        :return:
        """

        return self.definition.load_young_stars_component()

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_component(self):

        """
        This function ...
        :return:
        """

        return self.definition.load_ionizing_stars_component()

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_component(self):

        """
        This function ...
        :return:
        """

        return self.definition.load_dust_disk_component()

    # -----------------------------------------------------------------

    @property
    def has_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.wavelength_grid is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths(self):

        """
        This function ...
        :return:
        """

        return self.wavelength_grid.wavelengths(unit="micron", add_unit=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths_micron(self):

        """
        This function ...
        :return:
        """

        return self.wavelength_grid.wavelengths(unit="micron", asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_deltas(self):

        """
        This function ...
        :return:
        """

        return self.wavelength_grid.deltas(unit="micron", add_unit=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_deltas_micron(self):

        """
        This function ...
        :return:
        """

        return self.wavelength_grid.deltas(unit="micron", asarray=True)

    # -----------------------------------------------------------------

    @property
    def has_wavelengths_directory(self):

        """
        This function ...
        :return:
        """

        return fs.contains_directory(self.definition_path, "wavelengths")

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.definition_path, "wavelengths")

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.wavelengths_path, "grid.txt")

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_skifile(self):

        """
        This function ...
        :return:
        """

        # Load the ski file if it already exists
        if self.has_old_bulge_skifile: return SkiFile(self.old_bulge_ski_path)

        # Create
        else: return self.create_old_bulge_skifile()

    # -----------------------------------------------------------------

    def create_old_bulge_skifile(self):

        """
        This function ...
        :return:
        """

        # Check whether the wavelength grid is defined
        if not self.has_wavelength_grid: raise ValueError("Wavelength grid path must be set")

        # Create a ski template
        ski = get_panchromatic_template()

        # Add the old stellar bulge component
        # print(self.old_bulge_component.parameters)
        add_new_stellar_component(ski, old_bulge_name, self.old_bulge_component)

        # Add the instrument
        ski.add_instrument(earth_name, self.sed_instrument)

        # Set the wavelength grid
        ski.set_file_wavelength_grid(wavelengths_filename)

        # Set the number of photon packages
        ski.setpackages(default_npackages)

        # Remove the dust system
        ski.remove_dust_system()

        # Save the skifile
        ski.saveto(self.old_bulge_ski_path, fix=True)

        # Return the skifile
        return ski

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_skifile(self):

        """
        This fnuction ...
        :return:
        """

        # Load the ski file if it already exists
        if self.has_old_disk_skifile: return SkiFile(self.old_disk_ski_path)

        # Create
        else: return self.create_old_disk_skifile()

    # -----------------------------------------------------------------

    def create_old_disk_skifile(self):

        """
        This function ...
        :return:
        """

        # Check whether the wavelength grid is defined
        if not self.has_wavelength_grid: raise ValueError("Wavelength grid path must be set")

        # Create a ski template
        ski = get_panchromatic_template()

        # Add the old stellar disk component
        # print(self.old_disk_component.parameters)
        add_new_stellar_component(ski, old_disk_name, self.old_disk_component)

        # Add the instrument
        ski.add_instrument(earth_name, self.sed_instrument)

        # Set the wavelength grid
        ski.set_file_wavelength_grid(wavelengths_filename)

        # Set the number of photon packages
        ski.setpackages(default_npackages)

        # Remove the dust system
        ski.remove_dust_system()

        # Save the skifile
        ski.saveto(self.old_disk_ski_path, fix=True)

        # Return the skifile
        return ski

    # -----------------------------------------------------------------

    @lazyproperty
    def young_skifile(self):

        """
        This function ...
        :return:
        """

        # Load the ski file if it already exists
        if self.has_young_skifile: return SkiFile(self.young_ski_path)

        # Create
        else: return self.create_young_skifile()

    # -----------------------------------------------------------------

    def create_young_skifile(self):

        """
        This function ...
        :return:
        """

        # Check whether the wavelength grid is defined
        if not self.has_wavelength_grid: raise ValueError("Wavelength grid path must be set")

        # Create a ski template
        ski = get_panchromatic_template()

        # Add the young stellar component
        # print(self.young_component.parameters)
        add_new_stellar_component(ski, young_name, self.young_component)

        # Add the instrument
        ski.add_instrument(earth_name, self.sed_instrument)

        # Set the wavelength grid
        ski.set_file_wavelength_grid(wavelengths_filename)

        # Set the number of photon packages
        ski.setpackages(default_npackages)

        # Remove the dust system
        ski.remove_dust_system()

        # Save the skifile
        ski.saveto(self.young_ski_path, fix=True)

        # Return the skifile
        return ski

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_skifile(self):

        """
        This function ...
        :return:
        """

        # Load the ski file if it already exists
        if self.has_sfr_skifile: return SkiFile(self.sfr_ski_path)

        # Create
        else: return self.create_sfr_skifile()

    # -----------------------------------------------------------------

    def create_sfr_skifile(self):

        """
        This function ...
        :return:
        """

        # Check whether the wavelength grid is defined
        if not self.has_wavelength_grid: raise ValueError("Wavelength grid path must be set")

        # Create a ski template
        ski = get_panchromatic_template()

        # Add the sfr component
        # print(self.sfr_component.parameters)
        add_new_stellar_component(ski, sfr_name, self.sfr_component)

        # Add the instrument
        ski.add_instrument(earth_name, self.sed_instrument)

        # Set the wavelength grid
        ski.set_file_wavelength_grid(wavelengths_filename)

        # Set the number of photon packages
        ski.setpackages(default_npackages)

        # Remove the dust system
        ski.remove_dust_system()

        # Save the skifile
        ski.saveto(self.sfr_ski_path, fix=True)

        # Return the skifile
        return ski

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_sed_out_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.old_bulge_sed_path, "out")

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_sed_out_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.old_disk_sed_path, "out")

    # -----------------------------------------------------------------

    @lazyproperty
    def young_sed_out_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.young_sed_path, "out")

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_sed_out_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.sfr_sed_path, "out")

    # -----------------------------------------------------------------

    @property
    def old_disk_map_path(self):

        """
        This function ...
        :return:
        """

        return self.definition.old_stars_map_path

    # -----------------------------------------------------------------

    @property
    def young_map_path(self):

        """
        This function ...
        :return:
        """

        return self.definition.young_stars_map_path

    # -----------------------------------------------------------------

    @property
    def sfr_map_path(self):

        """
        This function ...
        :return:
        """

        return self.definition.ionizing_stars_map_path

    # -----------------------------------------------------------------

    @property
    def dust_map_path(self):

        """
        This function ...
        :return:
        """

        return self.definition.dust_map_path

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_input_filepaths(self):

        """
        This function ...
        :return:
        """

        paths = OrderedDict()
        paths[wavelengths_filename] = self.wavelength_grid_path
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_input_filepaths(self):

        """
        This function ...
        :return:
        """

        paths = OrderedDict()
        paths[wavelengths_filename] = self.wavelength_grid_path
        paths[map_filename] = self.old_disk_map_path
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def young_input_filepaths(self):

        """
        This function ...
        :return:
        """

        paths = OrderedDict()
        paths[wavelengths_filename] = self.wavelength_grid_path
        paths[map_filename] = self.young_map_path
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_input_filepaths(self):

        """
        This function ...
        :return:
        """

        paths = OrderedDict()
        paths[wavelengths_filename] = self.wavelength_grid_path
        paths[map_filename] = self.sfr_map_path
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_definition(self):

        """
        This function ...
        :return:
        """

        # Create the skifile if necessary
        if not self.has_old_bulge_skifile: self.create_old_bulge_skifile()

        # Create the definition and return
        return SingleSimulationDefinition(self.old_bulge_ski_path, self.old_bulge_sed_out_path, input_path=self.old_bulge_input_filepaths, name=old_bulge_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_definition(self):

        """
        This function ...
        :return:
        """

        # Create the skifile if necessary
        if not self.has_old_disk_skifile: self.create_old_disk_skifile()

        # Create the definition and return
        return SingleSimulationDefinition(self.old_disk_ski_path, self.old_disk_sed_out_path, input_path=self.old_disk_input_filepaths, name=old_disk_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_definition(self):

        """
        This function ...
        :return:
        """

        # Create the skifile if necessary
        if not self.has_young_skifile: self.create_young_skifile()

        # Create the definition and return
        return SingleSimulationDefinition(self.young_ski_path, self.young_sed_out_path, input_path=self.young_input_filepaths, name=young_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_definition(self):

        """
        This function ...
        :return:
        """

        # Create the skifile if necessary
        if not self.has_sfr_skifile: self.create_sfr_skifile()

        # Create the definition and return
        return SingleSimulationDefinition(self.sfr_ski_path, self.sfr_sed_out_path, input_path=self.sfr_input_filepaths, name=sfr_name)

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_output(self):

        """
        This function ...
        :return:
        """

        return not fs.is_empty(self.old_bulge_sed_out_path)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_output(self):

        """
        This function ...
        :return:
        """

        return not fs.is_empty(self.old_disk_sed_out_path)

    # -----------------------------------------------------------------

    @property
    def has_young_output(self):

        """
        This function ...
        :return:
        """

        return not fs.is_empty(self.young_sed_out_path)

    # -----------------------------------------------------------------

    @property
    def has_sfr_output(self):

        """
        This function ...
        :return:
        """

        return not fs.is_empty(self.sfr_sed_out_path)

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_sed(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.old_bulge_sed_out_path, extension="dat", endswith="_sed")

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_simulation(self):

        """
        This function ...
        :return:
        """

        # Simulation already performed?
        #if self.has_old_bulge_output: return createsimulations(self.old_bulge_sed_out_path, single=True)
        if self.has_old_bulge_sed: return createsimulations(self.old_bulge_sed_out_path, single=True)

        # Show message
        log.info("Running SKIRT for the old stellar bulge component ...")

        # Run simulation
        return run_simulation(self.old_bulge_definition, show_progress=True, debug_output=True)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_sed(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.old_disk_sed_out_path, extension="dat", endswith="_sed")

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_simulation(self):

        """
        This function ...
        :return:
        """

        # Simulation already performed?
        #if self.has_old_disk_output: return createsimulations(self.old_disk_sed_out_path, single=True)
        if self.has_old_disk_sed: return createsimulations(self.old_disk_sed_out_path, single=True)

        # Show message
        log.info("Running SKIRT for the old stellar disk component ...")

        # Run
        return run_simulation(self.old_disk_definition, show_progress=True, debug_output=True)

    # -----------------------------------------------------------------

    @property
    def has_young_sed(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.young_sed_out_path, extension="dat", endswith="_sed")

    # -----------------------------------------------------------------

    @lazyproperty
    def young_simulation(self):

        """
        This function ...
        :return:
        """

        # Simulation already performed?
        #if self.has_young_output: return createsimulations(self.young_sed_out_path, single=True)
        if self.has_young_sed: return createsimulations(self.young_sed_out_path, single=True)

        # Show message
        log.info("Running SKIRT simulation for the young stellar component ...")

        # Run
        return run_simulation(self.young_definition, show_progress=True, debug_output=True)

    # -----------------------------------------------------------------

    @property
    def has_sfr_sed(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.sfr_sed_out_path, extension="dat", endswith="_sed")

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_simulation(self):

        """
        This function ...
        :return:
        """

        # Simulation already performed?
        #if self.has_sfr_output: return createsimulations(self.sfr_sed_out_path, single=True)
        if self.has_sfr_sed: return createsimulations(self.sfr_sed_out_path, single=True)

        # Show message
        log.info("Running SKIRT simulation for the SFR component ...")

        # Run
        return run_simulation(self.sfr_definition, show_progress=True, debug_output=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_output(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_simulation.output

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_sed_filepath(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_output.single_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_output(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_simulation.output

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_sed_filepath(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_output.single_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def young_output(self):

        """
        This fucntion ...
        :return:
        """

        return self.young_simulation.output

    # -----------------------------------------------------------------

    @lazyproperty
    def young_sed_filepath(self):

        """
        This function ...
        :return:
        """

        return self.young_output.single_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_output(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulation.output

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_sed_filepath(self):

        """
        This function ...
        :return:
        """

        return self.sfr_output.single_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sed_old_bulge(self):

        """
        This function ...
        :return:
        """

        return SED.from_skirt(self.old_bulge_sed_filepath)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_old_bulge(self):

        """
        This function ...
        :return:
        """

        return True # should always be able to be created from a new simulation

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sed_old_disk(self):

        """
        This function ...
        :return:
        """

        return SED.from_skirt(self.old_disk_sed_filepath)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_old_disk(self):

        """
        This function ...
        :return:
        """

        return True # should always be able to be created from a new simulation

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sed_old(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed_old_bulge + self.intrinsic_sed_old_disk

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_old(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed_old_bulge and self.has_intrinsic_sed_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sed_young(self):

        """
        This function ...
        :return:
        """

        return SED.from_skirt(self.young_sed_filepath)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_young(self):

        """
        This function ...
        :return:
        """

        return True # should always be able to be created from a new simulation

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sed_sfr(self):

        """
        This function ...
        :return:
        """

        return SED.from_skirt(self.sfr_sed_filepath)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_sfr(self):

        """
        This function ...
        :return:
        """

        return True # should always be able to be created from a new simulation

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sed_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed_young + self.intrinsic_sed_sfr

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed_young and self.has_intrinsic_sed_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sed_sfr_stellar(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed_sfr.splice(x_max=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_sfr_stellar(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sed_sfr_dust(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed_sfr.splice(x_min=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_sfr_dust(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_old_bulge(self):

        """
        This function ...
        :return:
        """

        # WE NEED A SEPARATE SIMULATION FOR JUST THE OLD BULGE FOR THIS
        return None

    # -----------------------------------------------------------------

    @property
    def has_observed_sed_old_bulge(self):

        """
        This function ...
        :return:
        """

        return False

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_old_disk(self):

        """
        This function ...
        :return:
        """

        # WE NEED A SEPARATE SIMULATION FOR JUST THE OLD DISK FOR THIS
        return None

    # -----------------------------------------------------------------

    @property
    def has_observed_sed_old_disk(self):

        """
        This function ...
        :return:
        """

        return False

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_old(self):

        """
        This function ...
        :return:
        """

        #return self.observed_sed_old_bulge + self.observed_sed_old_disk
        return self.observed_old_data.seds[earth_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def has_observed_sed_old(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_old_output and self.observed_old_data.has_seds

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_young(self):

        """
        This function ...
        :return:
        """

        return self.observed_young_data.seds[earth_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def has_observed_sed_young(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_young_output and self.observed_young_data.has_seds

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_sfr(self):

        """
        This function ...
        :return:
        """

        return self.observed_sfr_data.seds[earth_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def has_observed_sed_sfr(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sfr_output and self.observed_sfr_data.has_seds

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.observed_unevolved_data.seds[earth_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def has_observed_sed_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_unevolved_output and self.observed_unevolved_data.has_seds

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_total.splice(x_max=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_total.splice(x_min=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_total_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        # based on self.observed_stellar_sed
        return self.observed_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_total_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_bolometric_luminosity_old + self.intrinsic_bolometric_luminosity_young + self.intrinsic_bolometric_luminosity_sfr

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_old and self.has_intrinsic_bolometric_luminosity_young and self.has_intrinsic_bolometric_luminosity_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_attenuation(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_stellar_luminosity.to("W", distance=self.distance).value
        intrinsic = self.intrinsic_stellar_luminosity.to("W", distance=self.distance).value
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_stellar_luminosity and self.has_intrinsic_stellar_luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity(self):

        """
        This function ...
        :return:
        """

        # based on self.dust_sed
        return self.dust_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.dust_sed - self.observed_sfr_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_sed and self.has_observed_sfr_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.diffuse_dust_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_diffuse_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sfr_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed_sfr.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sfr_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        # From intrinsic_sed_sfr_stellar
        return self.intrinsic_sed_sfr_stellar.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sfr_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sfr_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        # From intrinsic_sed_sfr_dust
        return self.intrinsic_sed_sfr_dust.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sfr_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed_sfr_dust

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sfr_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_sfr.splice(x_max=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sfr_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sfr_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_sfr.splice(x_min=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sfr_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sfr_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_sfr.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sfr_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        # From observed_sfr_stellar_sed
        return self.observed_sfr_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sfr_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        # From observed_sfr_dust_sed
        return self.observed_sfr_dust_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sfr_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sfr_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_sfr.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_total_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_total(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Total
        if self.has_observed_total_bolometric_luminosity: values[obs_total_bol_lum_name] = self.observed_total_bolometric_luminosity
        if self.has_intrinsic_total_bolometric_luminosity: values[intr_total_bol_lum_name] = self.intrinsic_total_bolometric_luminosity
        if self.has_observed_stellar_luminosity: values[obs_total_stellar_bol_lum_name] = self.observed_stellar_luminosity
        if self.has_intrinsic_stellar_luminosity: values[intr_total_stellar_bol_lum_name] = self.intrinsic_stellar_luminosity
        if self.has_bolometric_attenuation: values[bol_attenuation_name] = self.bolometric_attenuation

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_bulge(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Old bulge
        if self.has_observed_i1_luminosity_old_bulge: values[obs_bulge_spec_lum_name] = self.observed_i1_luminosity_old_bulge
        if self.has_intrinsic_i1_luminosity_old_bulge: values[intr_bulge_spec_lum_name] = self.intrinsic_i1_luminosity_old_bulge # part of parameter set
        if self.has_observed_bolometric_luminosity_old_bulge: values[obs_bulge_bol_lum_name] = self.observed_bolometric_luminosity_old_bulge
        if self.has_intrinsic_bolometric_luminosity_old_bulge: values[intr_bulge_bol_lum_name] = self.intrinsic_bolometric_luminosity_old_bulge
        if self.has_i1_attenuation_old_bulge: values[bulge_spec_attenuation_name] = self.i1_attenuation_old_bulge
        if self.has_bolometric_attenuation_old_bulge: values[bulge_bol_attenuation_name] = self.bolometric_attenuation_old_bulge

        if self.has_observed_old_bulge_total_luminosity: values[obs_bulge_total_lum_name] = self.observed_old_bulge_total_luminosity
        if self.has_observed_old_bulge_dust_luminosity: values[obs_bulge_dust_lum_name] = self.observed_old_bulge_dust_luminosity

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_disk(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Old disk
        if self.has_observed_i1_luminosity_old_disk: values[obs_disk_spec_lum_name] = self.observed_i1_luminosity_old_disk
        if self.has_intrinsic_i1_luminosity_old_disk: values[intr_disk_spec_lum_name] = self.intrinsic_i1_luminosity_old_disk # part of parameter set
        if self.has_observed_bolometric_luminosity_old_disk: values[obs_disk_bol_lum_name] = self.observed_bolometric_luminosity_old_disk
        if self.has_intrinsic_bolometric_luminosity_old_disk: values[intr_disk_bol_lum_name] = self.intrinsic_bolometric_luminosity_old_disk
        if self.has_i1_attenuation_old_disk: values[disk_spec_attenuation_name] = self.i1_attenuation_old_disk
        if self.has_bolometric_attenuation_old_disk: values[disk_bol_attenuation_name] = self.bolometric_attenuation_old_disk

        if self.has_observed_old_disk_total_luminosity: values[obs_disk_total_lum_name] = self.observed_old_disk_total_luminosity
        if self.has_observed_old_disk_dust_luminosity: values[obs_disk_dust_lum_name] = self.observed_old_disk_dust_luminosity

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_old(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Old (evolved)
        if self.has_observed_i1_luminosity_old: values[obs_old_spec_lum_name] = self.observed_i1_luminosity_old
        if self.has_intrinsic_i1_luminosity_old: values[intr_old_spec_lum_name] = self.intrinsic_i1_luminosity_old
        if self.has_observed_bolometric_luminosity_old: values[obs_old_bol_lum_name] = self.observed_bolometric_luminosity_old
        if self.has_intrinsic_bolometric_luminosity_old: values[intr_old_bol_lum_name] = self.intrinsic_bolometric_luminosity_old
        if self.has_i1_attenuation_old: values[old_spec_attenuation_name] = self.i1_attenuation_old
        if self.has_bolometric_attenuation_old: values[old_bol_attenuation_name] = self.bolometric_attenuation_old

        if self.has_observed_old_total_luminosity: values[obs_old_total_lum_name] = self.observed_old_total_luminosity
        if self.has_observed_old_dust_luminosity: values[obs_old_dust_lum_name] = self.observed_old_dust_luminosity

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_young(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Young stars
        if self.has_observed_fuv_luminosity_young: values[obs_young_spec_lum_name] = self.observed_fuv_luminosity_young
        if self.has_intrinsic_fuv_luminosity_young: values[intr_young_spec_lum_name] = self.intrinsic_fuv_luminosity_young # part of (free) parameter set
        if self.has_observed_bolometric_luminosity_young: values[obs_young_bol_lum_name] = self.observed_bolometric_luminosity_young
        if self.has_intrinsic_bolometric_luminosity_young: values[intr_young_bol_lum_name] = self.intrinsic_bolometric_luminosity_young
        if self.has_fuv_attenuation_young: values[young_spec_attenuation_name] = self.fuv_attenuation_young
        if self.has_bolometric_attenuation_young: values[young_bol_attenuation_name] = self.bolometric_attenuation_young

        if self.has_observed_young_total_luminosity: values[obs_young_total_lum_name] = self.observed_young_total_luminosity
        if self.has_observed_young_dust_luminosity: values[obs_young_dust_lum_name] = self.observed_young_dust_luminosity

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_sfr(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Ionizing stars (SFR)
        if self.has_sfr: values[sfr_name] = self.sfr
        if self.has_observed_fuv_luminosity_sfr: values[obs_sfr_spec_lum_name] = self.observed_fuv_luminosity_sfr
        if self.has_intrinsic_fuv_luminosity_sfr: values[intr_sfr_spec_lum_name] = self.intrinsic_fuv_luminosity_sfr # part of the (free) parameter set
        if self.has_observed_bolometric_luminosity_sfr: values[obs_sfr_bol_lum_name] = self.observed_bolometric_luminosity_sfr # SEE: observed_sfr_bolometric_luminosity
        if self.has_intrinsic_bolometric_luminosity_sfr: values[intr_sfr_bol_lum_name] = self.intrinsic_bolometric_luminosity_sfr # SEE: intrinsic_sfr_bolometric_luminosity
        if self.has_fuv_attenuation_sfr: values[sfr_spec_attenuation_name] = self.fuv_attenuation_sfr
        if self.has_bolometric_attenuation_sfr: values[sfr_bol_attenuation_name] = self.bolometric_attenuation_sfr
        if self.has_sfr_stellar_mass: values[sfr_stellar_mass_name] = self.sfr_stellar_mass
        if self.has_observed_sfr_stellar_luminosity: values[obs_sfr_stellar_bol_lum_name] = self.observed_sfr_stellar_luminosity
        if self.has_intrinsic_sfr_stellar_luminosity: values[intr_sfr_stellar_bol_lum_name] = self.intrinsic_sfr_stellar_luminosity
        if self.has_sfr_dust_mass: values[sfr_dust_mass_name] = self.sfr_dust_mass
        if self.has_intrinsic_sfr_dust_luminosity: values[sfr_dust_lum_name] = self.intrinsic_sfr_dust_luminosity # intrinsic so only the dust IN MAPPINGS

        if self.has_observed_sfr_total_luminosity: values[obs_sfr_total_lum_name] = self.observed_sfr_total_luminosity
        if self.has_observed_sfr_dust_luminosity: values[obs_sfr_dust_lum_name] = self.observed_sfr_dust_luminosity

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_unevolved(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Young + ionizing (unevolved)
        if self.has_observed_fuv_luminosity_unevolved: values[obs_unevolved_spec_lum_name] = self.observed_fuv_luminosity_unevolved
        if self.has_intrinsic_fuv_luminosity_unevolved: values[intr_unevolved_spec_lum_name] = self.intrinsic_fuv_luminosity_unevolved
        if self.has_observed_bolometric_luminosity_unevolved: values[obs_unevolved_bol_lum_name] = self.observed_bolometric_luminosity_unevolved
        if self.has_intrinsic_bolometric_luminosity_unevolved: values[intr_unevolved_bol_lum_name] = self.intrinsic_bolometric_luminosity_unevolved
        if self.has_fuv_attenuation_unevolved: values[unevolved_spec_attenuation_name] = self.fuv_attenuation_unevolved
        if self.has_bolometric_attenuation_unevolved: values[unevolved_bol_attenuation_name] = self.bolometric_attenuation_unevolved

        if self.has_observed_unevolved_total_luminosity: values[obs_unevolved_total_lum_name] = self.observed_unevolved_total_luminosity
        if self.has_observed_unevolved_dust_luminosity: values[obs_unevolved_dust_lum_name] = self.observed_unevolved_dust_luminosity

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_dust(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Dust
        if self.has_dust_luminosity: values[dust_lum_name] = self.dust_luminosity
        if self.has_diffuse_dust_luminosity: values[diffuse_dust_lum_name] = self.diffuse_dust_luminosity
        if self.has_total_dust_mass: values[total_dust_mass_name] = self.total_dust_mass # with SFR dust mass

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Add values
        values.update(self.derived_parameter_values_total)
        values.update(self.derived_parameter_values_bulge)
        values.update(self.derived_parameter_values_disk)
        values.update(self.derived_parameter_values_old)
        values.update(self.derived_parameter_values_young)
        values.update(self.derived_parameter_values_sfr)
        values.update(self.derived_parameter_values_unevolved)

        # Return
        return values

# -----------------------------------------------------------------
