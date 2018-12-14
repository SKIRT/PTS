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
import warnings
import numpy as np
from copy import deepcopy
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
from ..basics.instruments import SEDInstrument
from ...core.tools import filesystem as fs
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ..simulation.single import SingleComponentSimulations
from ..simulation.multi import MultiComponentSimulations
from ..simulation.simulation import ObservedComponentSimulation
from ...magic.core.frame import Frame
from ...magic.core.list import convolve_and_rebin, convolve_rebin_and_convert
from ...magic.basics.coordinatesystem import CoordinateSystem
from ..basics.projection import GalaxyProjection
from ..basics.instruments import FullSEDInstrument
from ..projection.model import ComponentProjections
from ..simulation.sed import ComponentSED
from ...core.basics.log import log
from ...core.simulation.output import SimulationOutput
from .sfr import salim_fuv_to_sfr, kennicutt_fuv_to_sfr, kennicutt_evans_fuv_to_sfr, kennicutt_tir_to_sfr, calzetti_24um_to_sfr
from .stellar_mass import oliver_stellar_mass, hubble_stage_to_type
from ...core.data.sed import SED
from ...core.tools.serialization import write_dict, load_dict
from ...core.units.parsing import parse_quantity
from ...core.data.attenuation import MappingsAttenuationCurve, AttenuationCurve

# -----------------------------------------------------------------

observed_name = "observed"
stellar_name = "stellar"
intrinsic_name = "intrinsic"

# -----------------------------------------------------------------

total_suffix = " [TOTAL]"
bulge_suffix = " [BULGE]"
disk_suffix = " [DISK]"
old_suffix = " [OLD]"
young_suffix = " [YOUNG]"
sfr_suffix = " [SFR]"
unevolved_suffix = " [UNEVOLVED]"
extra_suffix = " [EXTRA]"
dust_suffix = " [DUST]"

# -----------------------------------------------------------------

# Names of derived model properties

## Total
obs_total_bol_lum_name = "Observed total bolometric luminosity" # 1
intr_total_bol_lum_name = "Intrinsic total bolometric luminosity" # 2 ; should be same as 1
obs_total_stellar_bol_lum_name = "Observed total stellar bolometric luminosity" # 3
intr_total_stellar_bol_lum_name = "Intrinsic total stellar bolometric luminosity" # 4 ; should be same as 2
dust_lum_name = "Bolometric dust luminosity" #
diffuse_dust_lum_name = "Bolometric diffuse dust luminosity"
diffuse_abs_stellar_lum_name = "Absorbed stellar luminosity by diffuse dust"
diffuse_fabs_name = "Fraction of absorbed stellar luminosity by diffuse dust"
fabs_name = "Fraction of absorbed stellar luminosity"
bol_attenuation_name = "Total bolometric attenuation" # 5
direct_stellar_lum_name = "Direct stellar luminosity"

# Special: some total, some young, some sfr, some unevolved, some multiple
sfr_salim_name = "Star formation rate (Salim)"
sfr_ke_name = "Star formation rate (Kennicutt&Evans)"
sfr_tir_name = "Star formation rate (TIR)"
sfr_24um_name = "Star formation rate (24um, Calzetti)"
sfr_mappings_name = "Star formation rate (MAPPINGS)"
sfr_mappings_ke_name = "Star formation rate (MAPPINGS+K&E)"
stellar_mass_name = "Stellar mass"
stellar_mass_intrinsic_name = "Stellar mass (from intrinsic luminosity)"
ssfr_salim_name = "Specific star formation rate (Salim)"
ssfr_ke_name = "Specific star formation rate (Kennicutt&Evans)"
ssfr_tir_name = "Specific star formation rate (TIR)"
ssfr_24um_name = "Specific star formation rate (24um, Calzetti)"
ssfr_mappings_name = "Specific star formation rate (MAPPINGS)"
ssfr_mappings_ke_name = "Specific star formation rate (MAPPINGS+K&E)"

# (specific) star formation rates for total
#sfr_salim_total_name = "Star formation rate (Salim, total)"
#sfr_ke_total_name = "Star formation rate (Kennicutt&Evans, total)"
#ssfr_salim_total_name = "Specific star formation rate (Salim, total)"
#ssfr_ke_total_name = "Specific star formation rate (Kennicutt&Evans, total)"

obs_fuv_spec_lum_name = "Observed FUV specific luminosity"
intr_fuv_spec_lum_name = "Intrinsic FUV specific luminosity"

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

## Extra component
obs_extra_spec_lum_name = "Observed extra stellar specific luminosity" #
intr_extra_spec_lum_name = "Intrinsic extra stellar specific luminosity" # part of (free) parameter set
obs_extra_bol_lum_name = "Observed extra stellar bolometric luminosity" #
intr_extra_bol_lum_name = "Intrinsic extra stellar bolometric luminosity" #
extra_spec_attenuation_name = "Extra stellar specific attenuation" #
extra_bol_attenuation_name = "Extra stellar bolometric attenuation" #
obs_extra_total_lum_name = "Extra stellar bolometric total luminosity" # stellar + dust ; should be the same as x (intrinsic lum)
obs_extra_dust_lum_name = "Extra stellar bolometric dust luminosity" #

## Dust
total_dust_mass_name = "Total dust mass" # 6 with SFR dust mass
diffuse_dust_mass_name = "Diffuse dust mass"

# -----------------------------------------------------------------

sed_dirname = "sed"
projections_dirname = "projections"
mappings_dirname = "mappings"
transparent_dirname = "transparent"

# -----------------------------------------------------------------

total_simulation_name = "total"
bulge_simulation_name = "bulge"
disk_simulation_name = "disk"
old_simulation_name = "old"
young_simulation_name = "young"
sfr_simulation_name = "sfr"
unevolved_simulation_name = "unevolved"
extra_simulation_name = "extra"
dust_simulation_name = "dust"

bulge_component_name = "Evolved stellar bulge"
disk_component_name = "Evolved stellar disk"
young_component_name = "Young stars"
ionizing_component_name = "Ionizing stars"
transparent_ionizing_component_name = "Ionizing stars (transparent)"
extra_component_name = "Extra component"
dust_component_name = "Dust"

bulge_component_description = "old bulge stellar component"
disk_component_description = "old disk stellar component"
young_component_description = "young stellar component"
ionizing_component_description = "ionizing stellar component"
transparent_ionizing_component_description = "ionizing stellar component (transparent)"
extra_component_description = "extra stellar/agn component"
dust_component_description = "dust component"

evolved_component_name = "Evolved stars"
unevolved_component_name = "Unevolved stars"

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
dust_direct_contribution = "dust_direct"
dust_scattered_contribution = "dust_scattered"
transparent_contribution = "transparent"
contributions = [total_contribution, direct_contribution, scattered_contribution, dust_contribution, dust_direct_contribution, dust_scattered_contribution, transparent_contribution]

# -----------------------------------------------------------------

components_name = "components"
output_filename = "output.txt"

# -----------------------------------------------------------------

# MAPPINGS parameters
sfr_parameter = "sfr"
metallicity_parameter = "metallicity"
compactness_parameter = "compactness"
pressure_parameter = "pressure"
covering_factor_parameter = "covering_factor"

# -----------------------------------------------------------------

class RTModel(object):

    """
    Objects of this class describe a full radiative transfer model.
    """

    def __init__(self, definition, wavelength_grid=None, simulation_name=None, chi_squared=None,
                 free_parameter_labels=None, free_parameter_values=None, observed_total_output_path=None,
                 observed_bulge_output_path=None, observed_disk_output_path=None, observed_old_output_path=None,
                 observed_young_output_path=None, observed_sfr_output_path=None, observed_unevolved_output_path=None,
                 observed_extra_output_path=None, parameters=None, center=None, galaxy_name=None, hubble_stage=None,
                 redshift=None, earth_wcs=None, truncation_ellipse=None):

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
        :param observed_extra_output_path: output path of simulation with extra stellar component
        :param parameters:
        :param center: the galaxy center as a sky coordinate
        :param galaxy_name: the name of the galaxy
        :param hubble_stage: the Hubble stage classification of the galaxy
        :param redshift: the redshift of the galaxy
        :param earth_wcs: the celestial coordinate system of all simulated earth datacubes
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
        self.observed_extra_output_path = observed_extra_output_path

        # Set parameters
        if parameters is not None: self.set_parameters(**parameters)

        # Set the center and other galaxy properties
        self.center = center
        self.galaxy_name = galaxy_name
        self.hubble_stage = hubble_stage
        self.redshift = redshift
        self.truncation_ellipse = truncation_ellipse

        # Set the earth coordinate system
        self.earth_wcs = earth_wcs

    # -----------------------------------------------------------------

    @property
    def has_extra_component(self):
        return self.observed_extra_output_path is not None

    # -----------------------------------------------------------------

    @property
    def has_earth_wcs(self):
        return self.earth_wcs is not None

    # -----------------------------------------------------------------

    @property
    def has_hubble_stage(self):
        return self.hubble_stage is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def hubble_type(self):
        if not self.has_hubble_stage: return None
        return hubble_stage_to_type(self.hubble_stage)

    # -----------------------------------------------------------------

    @lazyproperty
    def hubble_subtype(self):
        if not self.has_hubble_stage: return None
        return hubble_stage_to_type(self.hubble_stage, add_subtype=True)[1]

    # -----------------------------------------------------------------

    @property
    def has_center(self):
        return self.center is not None

    # -----------------------------------------------------------------

    @property
    def has_truncation_ellipse(self):
        return self.truncation_ellipse is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_radius(self):
        if not self.has_truncation_ellipse: return None
        return self.truncation_ellipse.semimajor

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_box(self):
        if not self.has_truncation_ellipse: return None
        return self.truncation_ellipse.bounding_box

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_box_axial_ratio(self):
        if not self.has_truncation_ellipse: return None
        return self.truncation_box.axial_ratio

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
    # MEAN STELLAR AGES
    # -----------------------------------------------------------------

    @property
    def old_mean_stellar_age(self):
        return self.definition.old_stars_age

    # -----------------------------------------------------------------

    @property
    def young_mean_stellar_age(self):
        return self.definition.young_stars_age

    # -----------------------------------------------------------------

    @property
    def sfr_mean_stellar_age(self):
        return self.definition.ionizing_stars_age

    # -----------------------------------------------------------------
    # TOTAL SIMULATION INTRINSIC SEDS
    # -----------------------------------------------------------------

    @property
    def intrinsic_sed_path_old_bulge(self):
        return self.old_bulge_sed_filepath

    # -----------------------------------------------------------------

    @property
    def intrinsic_sed_path_old_disk(self):
        return self.old_disk_sed_filepath

    # -----------------------------------------------------------------

    @property
    def intrinsic_sed_path_young(self):
        return self.young_sed_filepath

    # -----------------------------------------------------------------

    @property
    def intrinsic_sed_path_sfr(self):
        return self.sfr_sed_filepath

    # -----------------------------------------------------------------

    @property
    def intrinsic_sed_path_extra(self):
        return self.extra_sed_filepath

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_component_sed_paths(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the intrinsic SEDs
        seds = OrderedDict()

        # Add
        seds[bulge_component_name] = self.intrinsic_sed_path_old_bulge
        seds[disk_component_name] = self.intrinsic_sed_path_old_disk
        seds[young_component_name] = self.intrinsic_sed_path_young
        seds[ionizing_component_name] = self.intrinsic_sed_path_sfr
        if self.has_extra_component: seds[extra_component_name] = self.intrinsic_sed_path_extra

        # Return the seds
        return seds

    # -----------------------------------------------------------------
    # TOTAL SIMULATION INTRINSIC CUBES
    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_evolved_component_cubes(self):
        cubes = OrderedDict()
        if self.bulge_simulations.has_intrinsic_cube and self.disk_simulations.has_intrinsic_cube:
            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_earth
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_earth
        elif self.old_simulations.has_intrinsic_cube: cubes[evolved_component_name] = self.old_intrinsic_stellar_luminosity_cube_earth # old
        else: warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the total simulation, the transparent earth cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_unevolved_component_cubes(self):
        cubes = OrderedDict()
        if self.young_simulations.has_intrinsic_cube and self.sfr_simulations.has_intrinsic_cube:
            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_earth
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_earth
        elif self.unevolved_simulations.has_intrinsic_cube: cubes[unevolved_component_name] = self.unevolved_intrinsic_stellar_luminosity_cube_earth
        else: warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the total simulation, the transparent earth cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_component_cubes(self):
        cubes =  OrderedDict(self.total_simulation_evolved_component_cubes.items() + self.total_simulation_unevolved_component_cubes.items())
        if self.has_extra_component:
            if not self.extra_simulations.has_intrinsic_cube: warnings.warn("No intrinsic cube for the extra component")
            else: cubes[extra_component_name] = self.extra_intrinsic_stellar_luminosity_cube_earth
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_evolved_component_cubes_faceon(self):
        cubes = OrderedDict()
        if self.bulge_simulations.has_intrinsic_cube_faceon and self.disk_simulations.has_intrinsic_cube_faceon:
            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_faceon
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_faceon
        elif self.old_simulations.has_intrinsic_cube_faceon: cubes[evolved_component_name] = self.old_intrinsic_stellar_luminosity_cube_faceon
        else: warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the total simulation, the transparent faceon cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_unevolved_component_cubes_faceon(self):
        cubes = OrderedDict()
        if self.young_simulations.has_intrinsic_cube_faceon and self.sfr_simulations.has_intrinsic_cube_faceon:
            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_faceon
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_faceon
        elif self.unevolved_simulations.has_intrinsic_cube_faceon: cubes[unevolved_component_name] = self.unevolved_intrinsic_stellar_luminosity_cube_faceon
        else: warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the total simulation, the transparent faceon cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_component_cubes_faceon(self):
        cubes = OrderedDict(self.total_simulation_evolved_component_cubes_faceon.items() + self.total_simulation_unevolved_component_cubes_faceon.items())
        if self.has_extra_component:
            if not self.extra_simulations.has_intrinsic_cube_faceon: warnings.warn("No intrinsic faceon cube for the extra component")
            else: cubes[extra_component_name] = self.extra_intrinsic_stellar_luminosity_cube_faceon
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_evolved_component_cubes_edgeon(self):
        cubes = OrderedDict()
        if self.bulge_simulations.has_intrinsic_cube_edgeon and self.disk_simulations.has_intrinsic_cube_edgeon:
            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_edgeon
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_edgeon
        elif self.old_simulations.has_intrinsic_cube_edgeon: cubes[evolved_component_name] = self.old_intrinsic_stellar_luminosity_cube_edgeon
        else: warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the total simulation, the transparent edgeon cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_unevolved_component_cubes_edgeon(self):
        cubes = OrderedDict()
        if self.young_simulations.has_intrinsic_cube_edgeon and self.sfr_simulations.has_intrinsic_cube_edgeon:
            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_edgeon
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_edgeon
        elif self.unevolved_simulations.has_intrinsic_cube_edgeon: cubes[unevolved_component_name] = self.unevolved_intrinsic_stellar_luminosity_cube_edgeon
        else: warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the total simulation, the transparent edgeon cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_component_cubes_edgeon(self):
        cubes = OrderedDict(self.total_simulation_evolved_component_cubes_edgeon.items() + self.total_simulation_unevolved_component_cubes_edgeon.items())
        if self.has_extra_component:
            if not self.extra_simulations.has_intrinsic_cube_edgeon: warnings.warn("No intrinsic edgeon cube for the extra component")
            else: cubes[extra_component_name] = self.extra_intrinsic_stellar_luminosity_cube_edgeon
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_components_path(self):
        return fs.create_directory_in(self.observed_total_simulation_path, components_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_components_path_earth(self):
        return fs.create_directory_in(self.total_simulation_components_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_components_path_faceon(self):
        return fs.create_directory_in(self.total_simulation_components_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_components_path_edgeon(self):
        return fs.create_directory_in(self.total_simulation_components_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_component_cube_paths(self):
        if fs.is_empty(self.total_simulation_components_path_earth): self.create_total_simulation_component_cubes_earth()
        return fs.files_in_path(self.total_simulation_components_path_earth, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_component_cube_paths_faceon(self):
        if fs.is_empty(self.total_simulation_components_path_faceon): self.create_total_simulation_component_cubes_faceon()
        return fs.files_in_path(self.total_simulation_components_path_faceon, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_component_cube_paths_edgeon(self):
        if fs.is_empty(self.total_simulation_components_path_edgeon): self.create_total_simulation_component_cubes_edgeon()
        return fs.files_in_path(self.total_simulation_components_path_edgeon, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    def create_total_simulation_component_cubes_earth(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the cubes of the components of the total simulation from the earth projection (this is a one-time process) ...")

        # Loop over the components
        for component_name in self.total_simulation_component_cubes:
            log.debug("Creating the " + component_name + " component cube ...")
            filepath = fs.join(self.total_simulation_components_path_earth, component_name + ".fits")
            cube = self.total_simulation_component_cubes[component_name]
            cube.saveto(filepath)

    # -----------------------------------------------------------------

    def create_total_simulation_component_cubes_faceon(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the cubes of the components of the total simulation from the faceon projection (this is a one-time process) ...")

        # Loop over the components
        for component_name in self.total_simulation_component_cubes_faceon:
            log.debug("Creating the " + component_name + " component cube ...")
            filepath = fs.join(self.total_simulation_components_path_faceon, component_name + ".fits")
            cube = self.total_simulation_component_cubes_faceon[component_name]
            cube.saveto(filepath)

    # -----------------------------------------------------------------

    def create_total_simulation_component_cubes_edgeon(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the cubes of the components of the total simulation from the edgeon projection (this is a one-time process) ...")

        # Loop over the components
        for component_name in self.total_simulation_component_cubes_edgeon:
            log.debug("Creating the " + component_name + " component cube ...")
            filepath = fs.join(self.total_simulation_components_path_edgeon, component_name + ".fits")
            cube = self.total_simulation_component_cubes_edgeon[component_name]
            cube.saveto(filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_evolved_component_simulations(self):

        """
        This function ...
        :return:
        """

        # Initialize dict
        simulations = OrderedDict()

        # If both bulge and disk simulations have native intrinsic cubes (i.e. full instrument cubes for each of the three instruments)
        if self.bulge_simulations.has_native_intrinsic_cubes and self.disk_simulations.has_native_intrinsic_cubes:

            simulations[bulge_component_name] = self.bulge_simulations
            simulations[disk_component_name] = self.disk_simulations

        # If the old simulation has native intrinsic cubes?
        else: 

            # Check whether old simulation has native intrinsic cubes (transparent cubes)
            if not self.old_simulations.has_native_intrinsic_cubes: warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the total simulation, the transparent cubes will not be available")

            # Add old as component regardless (to define all components of the total simulation)
            simulations[evolved_component_name] = self.old_simulations

        # Return the simulations
        return simulations

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_unevolved_component_simulations(self):

        """
        This function ...
        :return:
        """

        # Initialize dict
        simulations = OrderedDict()

        # If both young and sfr simulations have native intrinsic cubes (i.e. full instrument cubes for each of the three instruments)
        if self.young_simulations.has_native_intrinsic_cubes and self.sfr_simulations.has_native_intrinsic_cubes:

            simulations[young_component_name] = self.young_simulations
            simulations[ionizing_component_name] = self.sfr_simulations

        # If the unevolved simulation has native intrinsic cubes?
        else:

            # Check whether the unevolved simulation has native intrinsic cubes (transparent cubes)
            if not self.unevolved_simulations.has_native_intrinsic_cubes: warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the total simulation, the transparent cubes will not be available")

            # Add unevolved as component regardless (to define all components of the total simulation)
            simulations[unevolved_component_name] = self.unevolved_simulations

        # Return the simulations
        return simulations

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_component_simulations(self):
        simulations = OrderedDict(self.total_simulation_evolved_component_simulations.items() + self.total_simulation_unevolved_component_simulations.items())
        if self.has_extra_component: simulations[extra_component_name] = self.extra_simulations
        return simulations

    # -----------------------------------------------------------------
    # TOTAL SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_total_simulation_path(self):
        return fs.directory_of(self.observed_total_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_total_simulation_output_filepath(self):
        return fs.join(self.observed_total_simulation_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_observed_total_simulation_output_file(self):
        return fs.is_file(self.observed_total_simulation_output_filepath)

    # -----------------------------------------------------------------

    def create_observed_total_simulation_output_file(self):
        output = SimulationOutput.from_directory(self.observed_total_output_path)
        output.saveto(self.observed_total_simulation_output_filepath)
        return output

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_total_simulation_output(self):
        if not self.has_observed_total_simulation_output_file: return self.create_observed_total_simulation_output_file()
        else: return SimulationOutput.from_file(self.observed_total_simulation_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulations(self):

        """
        This function ...
        :return:
        """

        simulations = OrderedDict()
        simulations[total_simulation_name] = self.total_simulations
        simulations[bulge_simulation_name] = self.bulge_simulations
        simulations[disk_simulation_name] = self.disk_simulations
        simulations[old_simulation_name] = self.old_simulations
        simulations[young_simulation_name] = self.young_simulations
        simulations[sfr_simulation_name] = self.sfr_simulations
        simulations[unevolved_simulation_name] = self.unevolved_simulations
        if self.has_extra_component: simulations[extra_simulation_name] = self.extra_simulations
        return simulations

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulations(self):

        """
        This function ...
        :return:
        """

        # Load and return
        return MultiComponentSimulations.from_output(total_simulation_name, self.observed_total_simulation_output, self.total_simulation_component_simulations,
                                                          intrinsic_sed_paths=self.total_simulation_component_sed_paths, distance=self.distance, earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation(self):
        return self.total_simulations.observed

    # -----------------------------------------------------------------

    @property
    def total_simulation_output(self):
        return self.total_simulation.output

    # -----------------------------------------------------------------

    @property
    def total_simulation_data(self):
        return self.total_simulation.data

    # -----------------------------------------------------------------
    # BULGE SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bulge_simulation_path(self):
        return fs.directory_of(self.observed_bulge_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_bulge_simulation_output_filepath(self):
        return fs.join(self.observed_bulge_simulation_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_observed_bulge_simulation_output_file(self):
        return fs.is_file(self.observed_bulge_simulation_output_filepath)

    # -----------------------------------------------------------------

    def create_observed_bulge_simulation_output_file(self):
        output = SimulationOutput.from_directory(self.observed_bulge_output_path)
        output.saveto(self.observed_bulge_simulation_output_filepath)
        return output

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bulge_simulation_output(self):
        if not self.has_observed_bulge_simulation_output_file: return self.create_observed_bulge_simulation_output_file()
        else: return SimulationOutput.from_file(self.observed_bulge_simulation_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_simulations(self):

        """
        This function ...
        :return:
        """

        # To run the intrinsic SED simulation
        sed = self.old_bulge_component_sed

        # Load and return
        return SingleComponentSimulations.from_output(bulge_simulation_name, self.observed_bulge_simulation_output,
                                                      intrinsic_output=sed.output, distance=self.distance,
                                                      map_earth_path=self.old_bulge_map_earth_path, map_faceon_path=self.old_bulge_map_faceon_path, map_edgeon_path=self.old_bulge_map_edgeon_path,
                                                      earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_simulation(self):
        return self.bulge_simulations.observed

    # -----------------------------------------------------------------

    @property
    def bulge_simulation_output(self):
        return self.bulge_simulation.output

    # -----------------------------------------------------------------

    @property
    def bulge_simulation_data(self):
        return self.bulge_simulation.data

    # -----------------------------------------------------------------
    # DISK SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_disk_simulation_path(self):
        return fs.directory_of(self.observed_disk_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_disk_simulation_output_filepath(self):
        return fs.join(self.observed_disk_simulation_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_observed_disk_simulation_output_file(self):
        return fs.is_file(self.observed_disk_simulation_output_filepath)

    # -----------------------------------------------------------------

    def create_observed_disk_simulation_output_file(self):
        output = SimulationOutput.from_directory(self.observed_disk_output_path)
        output.saveto(self.observed_disk_simulation_output_filepath)
        return output

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_disk_simulation_output(self):
        if not self.has_observed_disk_simulation_output_file: return self.create_observed_disk_simulation_output_file()
        else: return SimulationOutput.from_file(self.observed_disk_simulation_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_simulations(self):

        """
        This function ...
        :return:
        """

        # To run the intrinsic SED simulation
        sed = self.old_disk_component_sed

        # Load and return
        return SingleComponentSimulations.from_output(disk_simulation_name, self.observed_disk_simulation_output,
                                                      intrinsic_output=sed.output, distance=self.distance,
                                                      map_earth_path=self.old_disk_map_earth_path, map_faceon_path=self.old_disk_map_faceon_path, map_edgeon_path=self.old_disk_map_edgeon_path,
                                                      earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @property
    def disk_simulation(self):
        return self.disk_simulations.observed

    # -----------------------------------------------------------------

    @property
    def disk_simulation_output(self):
        return self.disk_simulation.output

    # -----------------------------------------------------------------

    @property
    def disk_simulation_data(self):
        return self.disk_simulation.data

    # -----------------------------------------------------------------
    # OLD SIMULATION INTRINSIC SEDs
    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_component_sed_paths(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the intrinsic SEDs
        seds = OrderedDict()

        # Add
        seds[bulge_component_name] = self.intrinsic_sed_path_old_bulge
        seds[disk_component_name] = self.intrinsic_sed_path_old_disk

        # Return
        return seds

    # -----------------------------------------------------------------
    # OLD SIMULATION INTRINSIC CUBES
    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_component_cubes(self):
        cubes = OrderedDict()
        if self.bulge_simulations.has_intrinsic_cube and self.disk_simulations.has_intrinsic_cube:
            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_earth
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_earth
        else: warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the old simulation, the transparent cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_component_cubes_faceon(self):
        cubes = OrderedDict()
        if self.bulge_simulations.has_intrinsic_cube_faceon and self.disk_simulations.has_intrinsic_cube_faceon:
            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_faceon
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_faceon
        else: warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the old simulation, the transparent faceon cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_component_cubes_edgeon(self):
        cubes = OrderedDict()
        if self.bulge_simulations.has_intrinsic_cube_edgeon and self.disk_simulations.has_intrinsic_cube_edgeon:
            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_edgeon
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_edgeon
        else: warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the old simulation, the transparent edgeon cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_components_path(self):
        return fs.create_directory_in(self.observed_old_simulation_path, components_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_components_path_earth(self):
        return fs.create_directory_in(self.old_simulation_components_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_components_path_faceon(self):
        return fs.create_directory_in(self.old_simulation_components_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_components_path_edgeon(self):
        return fs.create_directory_in(self.old_simulation_components_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_component_cube_paths(self):
        if fs.is_empty(self.old_simulation_components_path_earth): self.create_old_simulation_component_cubes_earth()
        return fs.files_in_path(self.old_simulation_components_path_earth, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_component_cube_paths_faceon(self):
        if fs.is_empty(self.old_simulation_components_path_faceon): self.create_old_simulation_component_cubes_faceon()
        return fs.files_in_path(self.old_simulation_components_path_faceon, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_component_cube_paths_edgeon(self):
        if fs.is_empty(self.old_simulation_components_path_edgeon): self.create_old_simulation_component_cubes_edgeon()
        return fs.files_in_path(self.old_simulation_components_path_edgeon, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    def create_old_simulation_component_cubes_earth(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the cubes of the components of the old simulation from the earth projection (this is a one-time process) ...")

        # Loop over the components
        for component_name in self.old_simulation_component_cubes:
            log.debug("Creating the " + component_name + " component cube ...")
            filepath = fs.join(self.old_simulation_components_path_earth, component_name + ".fits")
            cube = self.old_simulation_component_cubes[component_name]
            cube.saveto(filepath)

    # -----------------------------------------------------------------

    def create_old_simulation_component_cubes_faceon(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the cubes of the components of the old simulation from the faceon projection (this is a one-time process) ...")

        # Loop over the components
        for component_name in self.old_simulation_component_cubes_faceon:
            log.debug("Creating the " + component_name + " component cube ...")
            filepath = fs.join(self.old_simulation_components_path_faceon, component_name + ".fits")
            cube = self.old_simulation_component_cubes_faceon[component_name]
            cube.saveto(filepath)

    # -----------------------------------------------------------------

    def create_old_simulation_component_cubes_edgeon(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the cubes of the components of the old simulation from the edgeon projection (this is a one-time process) ...")

        # Loop over the components
        for component_name in self.old_simulation_component_cubes_edgeon:
            log.debug("Creating the " + component_name + " component cube ...")
            filepath = fs.join(self.old_simulation_components_path_edgeon, component_name + ".fits")
            cube = self.old_simulation_component_cubes_edgeon[component_name]
            cube.saveto(filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_component_simulations(self):

        """
        This function ...
        :return:
        """

        # Initialize dict
        simulations = OrderedDict()

        # If both bulge and disk simulations have native intrinsic cubes (i.e. full instrument cubes for each of the three instruments)
        if not (self.bulge_simulations.has_native_intrinsic_cubes and self.disk_simulations.has_native_intrinsic_cubes): warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the old simulation, the transparent cubes will not be available")

        # Add simulations regardless
        simulations[bulge_component_name] = self.bulge_simulations
        simulations[disk_component_name] = self.disk_simulations

        # Return the simulations
        return simulations

    # -----------------------------------------------------------------
    # OLD SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_simulation_path(self):
        return fs.directory_of(self.observed_old_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_old_simulation_output_filepath(self):
        return fs.join(self.observed_old_simulation_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_simulation_output_file(self):
        return fs.is_file(self.observed_old_simulation_output_filepath)

    # -----------------------------------------------------------------

    def create_observed_old_simulation_output_file(self):
        output = SimulationOutput.from_directory(self.observed_old_output_path)
        output.saveto(self.observed_old_simulation_output_filepath)
        return output

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_simulation_output(self):
        if not self.has_observed_old_simulation_output_file: return self.create_observed_old_simulation_output_file()
        else: return SimulationOutput.from_file(self.observed_old_simulation_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulations(self):

        """
        This function ...
        :return:
        """

        # Load and return
        return MultiComponentSimulations.from_output(old_simulation_name, self.observed_old_simulation_output, self.old_simulation_component_simulations,
                                                      intrinsic_sed_paths=self.old_simulation_component_sed_paths,
                                                      distance=self.distance, earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation(self):
        return self.old_simulations.observed

    # -----------------------------------------------------------------

    @property
    def old_simulation_output(self):
        return self.old_simulation.output

    # -----------------------------------------------------------------

    @property
    def old_simulation_data(self):
        return self.old_simulation.data

    # -----------------------------------------------------------------
    # YOUNG SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_young_simulation_path(self):
        return fs.directory_of(self.observed_young_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_young_simulation_output_filepath(self):
        return fs.join(self.observed_young_simulation_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_observed_young_simulation_output_file(self):
        return fs.is_file(self.observed_young_simulation_output_filepath)

    # -----------------------------------------------------------------

    def create_observed_young_simulation_output_file(self):
        output = SimulationOutput.from_directory(self.observed_young_output_path)
        output.saveto(self.observed_young_simulation_output_filepath)
        return output

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_young_simulation_output(self):
        if not self.has_observed_young_simulation_output_file: return self.create_observed_young_simulation_output_file()
        else: return SimulationOutput.from_file(self.observed_young_simulation_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_simulations(self):

        """
        This function ...
        :return:
        """

        # To run the simulation
        sed = self.young_component_sed

        # Load and return
        return SingleComponentSimulations.from_output(young_simulation_name, self.observed_young_simulation_output,
                                                      intrinsic_output=sed.output, distance=self.distance,
                                                      map_earth_path=self.young_map_earth_path, map_faceon_path=self.young_map_faceon_path, map_edgeon_path=self.young_map_edgeon_path,
                                                      earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_simulation(self):
        return self.young_simulations.observed

    # -----------------------------------------------------------------

    @property
    def young_simulation_output(self):
        return self.young_simulation.output

    # -----------------------------------------------------------------

    @property
    def young_simulation_data(self):
        return self.young_simulation.data

    # -----------------------------------------------------------------
    # SFR SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sfr_simulation_path(self):
        return fs.directory_of(self.observed_sfr_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_sfr_simulation_output_filepath(self):
        return fs.join(self.observed_sfr_simulation_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_simulation_output_file(self):
        return fs.is_file(self.observed_sfr_simulation_output_filepath)

    # -----------------------------------------------------------------

    def create_observed_sfr_simulation_output_file(self):
        output = SimulationOutput.from_directory(self.observed_sfr_output_path)
        output.saveto(self.observed_sfr_simulation_output_filepath)
        return output

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sfr_simulation_output(self):
        if not self.has_observed_sfr_simulation_output_file: return self.create_observed_sfr_simulation_output_file()
        else: return SimulationOutput.from_file(self.observed_sfr_simulation_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_simulations(self):

        """
        This function ...
        :return:
        """

        # Load and return
        return SingleComponentSimulations.from_output(sfr_simulation_name, self.observed_sfr_simulation_output,
                                                      intrinsic_output=self.sfr_component_sed.output, distance=self.distance,
                                                      map_earth_path=self.sfr_map_earth_path, map_faceon_path=self.sfr_map_faceon_path, map_edgeon_path=self.sfr_map_edgeon_path,
                                                      earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_simulation(self):
        return self.sfr_simulations.observed

    # -----------------------------------------------------------------

    @property
    def sfr_simulation_output(self):
        return self.sfr_simulation.output

    # -----------------------------------------------------------------

    @property
    def sfr_simulation_data(self):
        return self.sfr_simulation.data

    # -----------------------------------------------------------------
    # EXTRA SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_extra_simulation_path(self):
        return fs.directory_of(self.observed_extra_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_extra_simulation_output_filepath(self):
        return fs.join(self.observed_extra_simulation_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_observed_extra_simulation_output_file(self):
        return fs.is_file(self.observed_extra_simulation_output_filepath)

    # -----------------------------------------------------------------

    def create_observed_extra_simulation_output_file(self):
        output = SimulationOutput.from_directory(self.observed_extra_output_path)
        output.saveto(self.observed_extra_simulation_output_filepath)
        return output

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_extra_simulation_output(self):
        if not self.has_observed_extra_simulation_output_file: return self.create_observed_extra_simulation_output_file()
        else: return SimulationOutput.from_file(self.observed_extra_simulation_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_simulations(self):

        """
        This function ...
        :return:
        """

        # To run the simulation
        sed = self.extra_component_sed

        # Load and return
        return SingleComponentSimulations.from_output(extra_simulation_name, self.observed_extra_simulation_output,
                                                      intrinsic_output=sed.output, distance=self.distance,
                                                      map_earth_path=self.extra_map_earth_path, map_faceon_path=self.extra_map_faceon_path, map_edgeon_path=self.extra_map_edgeon_path,
                                                      earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_simulation(self):
        return self.extra_simulations.observed

    # -----------------------------------------------------------------

    @property
    def extra_simulation_output(self):
        return self.extra_simulation.output

    # -----------------------------------------------------------------

    @property
    def extra_simulation_data(self):
        return self.extra_simulation.data

    # -----------------------------------------------------------------
    # UNEVOLVED SIMULATION INTRINSIC SEDs
    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_sed_paths(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the intrinsic SEDs
        seds = OrderedDict()

        # Add
        seds[young_component_name] = self.intrinsic_sed_path_young
        seds[ionizing_component_name] = self.intrinsic_sed_path_sfr

        # Return
        return seds

    # -----------------------------------------------------------------
    # UNEVOLVED SIMULATION INTRINSIC CUBES
    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_cubes(self):
        cubes = OrderedDict()
        if self.young_simulations.has_intrinsic_cube and self.sfr_simulations.has_intrinsic_cube:
            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_earth
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_earth
        else: warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the unevolved simulation, the transparent cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_cubes_faceon(self):
        cubes = OrderedDict()
        if self.young_simulations.has_intrinsic_cube_faceon and self.sfr_simulations.has_intrinsic_cube_faceon:
            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_faceon
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_faceon
        else: warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the unevolved simulation, the transparent faceon cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_cubes_edgeon(self):
        cubes = OrderedDict()
        if self.young_simulations.has_intrinsic_cube_edgeon and self.sfr_simulations.has_intrinsic_cube_edgeon:
            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_edgeon
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_edgeon
        else: warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the unevolved simulation, the transparent edgeon cube will not be available")
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_components_path(self):
        return fs.create_directory_in(self.observed_unevolved_simulation_path, components_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_components_path_earth(self):
        return fs.create_directory_in(self.unevolved_simulation_components_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_components_path_faceon(self):
        return fs.create_directory_in(self.unevolved_simulation_components_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_components_path_edgeon(self):
        return fs.create_directory_in(self.unevolved_simulation_components_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_cube_paths(self):
        if fs.is_empty(self.unevolved_simulation_components_path_earth): self.create_unevolved_simulation_component_cubes_earth()
        return fs.files_in_path(self.unevolved_simulation_components_path_earth, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_cube_paths_faceon(self):
        if fs.is_empty(self.unevolved_simulation_components_path_faceon): self.create_unevolved_simulation_component_cubes_faceon()
        return fs.files_in_path(self.unevolved_simulation_components_path_faceon, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_cube_paths_edgeon(self):
        if fs.is_empty(self.unevolved_simulation_components_path_edgeon): self.create_unevolved_simulation_component_cubes_edgeon()
        return fs.files_in_path(self.unevolved_simulation_components_path_edgeon, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_cube_paths(self):
        if fs.is_empty(self.unevolved_simulation_components_path_earth): self.create_unevolved_simulation_component_cubes_earth()
        return fs.files_in_path(self.unevolved_simulation_components_path_earth, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_cube_paths_faceon(self):
        if fs.is_empty(self.unevolved_simulation_components_path_faceon): self.create_unevolved_simulation_component_cubes_faceon()
        return fs.files_in_path(self.unevolved_simulation_components_path_faceon, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_cube_paths_edgeon(self):
        if fs.is_empty(self.unevolved_simulation_components_path_edgeon): self.create_unevolved_simulation_component_cubes_edgeon()
        return fs.files_in_path(self.unevolved_simulation_components_path_edgeon, extension="fits", returns="dict")

    # -----------------------------------------------------------------

    def create_unevolved_simulation_component_cubes_earth(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the cubes of the components of the unevolved simulation from the earth projection (this is a one-time process) ...")

        # Loop over the components
        for component_name in self.unevolved_simulation_component_cubes:
            log.debug("Creating the " + component_name + " component cube ...")
            filepath = fs.join(self.unevolved_simulation_components_path_earth, component_name + ".fits")
            cube = self.unevolved_simulation_component_cubes[component_name]
            cube.saveto(filepath)

    # -----------------------------------------------------------------

    def create_unevolved_simulation_component_cubes_faceon(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the cubes of the components of the unevolved simulation from the earth projection (this is a one-time process) ...")

        # Loop over the components
        for component_name in self.unevolved_simulation_component_cubes_faceon:
            log.debug("Creating the " + component_name + " component cube ...")
            filepath = fs.join(self.unevolved_simulation_components_path_faceon, component_name + ".fits")
            cube = self.unevolved_simulation_component_cubes_faceon[component_name]
            cube.saveto(filepath)

    # -----------------------------------------------------------------

    def create_unevolved_simulation_component_cubes_edgeon(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the cubes of the components of the unevolved simulation from the earth projection (this is a one-time process) ...")

        # Loop over the components
        for component_name in self.unevolved_simulation_component_cubes_edgeon:
            log.debug("Creating the " + component_name + " component cube ...")
            filepath = fs.join(self.unevolved_simulation_components_path_edgeon, component_name + ".fits")
            cube = self.unevolved_simulation_component_cubes_edgeon[component_name]
            cube.saveto(filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_simulations(self):

        """
        Thisn function ...
        :return:
        """

        # Initialize dict
        simulations = OrderedDict()

        # If both young and sfr simulations have native intrinsic cubes (i.e. full instrument cubes for each of the three instruments)
        if not (self.young_simulations.has_native_intrinsic_cubes and self.sfr_simulations.has_native_intrinsic_cubes): warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the unevolved simulation, the transparent cubes will not be available")

        # Add simulations regardless (to define which are the unevolved components)
        simulations[young_component_name] = self.young_simulations
        simulations[ionizing_component_name] = self.sfr_simulations

        # Return the simulations
        return simulations

    # -----------------------------------------------------------------
    # UNEVOLVED SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_unevolved_simulation_path(self):
        return fs.directory_of(self.observed_unevolved_output_path)

    # -----------------------------------------------------------------

    @property
    def observed_unevolved_simulation_output_filepath(self):
        return fs.join(self.observed_unevolved_simulation_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_simulation_output_file(self):
        return fs.is_file(self.observed_unevolved_simulation_output_filepath)

    # -----------------------------------------------------------------

    def create_observed_unevolved_simulation_output_file(self):
        output = SimulationOutput.from_directory(self.observed_unevolved_output_path)
        output.saveto(self.observed_unevolved_simulation_output_filepath)
        return output

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_unevolved_simulation_output(self):
        if not self.has_observed_unevolved_simulation_output_file: return self.create_observed_unevolved_simulation_output_file()
        else: return SimulationOutput.from_file(self.observed_unevolved_simulation_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulations(self):

        """
        This function ...
        :return:
        """

        # Load and return
        return MultiComponentSimulations.from_output(unevolved_simulation_name, self.observed_unevolved_simulation_output, self.unevolved_simulation_component_simulations,
                                                      intrinsic_sed_paths=self.unevolved_simulation_component_sed_paths,
                                                      distance=self.distance, earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation(self):
        return self.unevolved_simulations.observed

    # -----------------------------------------------------------------

    @property
    def unevolved_simulation_output(self):
        return self.unevolved_simulation.output

    # -----------------------------------------------------------------

    @property
    def unevolved_simulation_data(self):
        return self.unevolved_simulation.data

    # -----------------------------------------------------------------
    # HAS OUTPUT
    # -----------------------------------------------------------------

    @property
    def has_observed_total_output(self):
        return self.observed_total_output_path is not None and fs.is_directory(self.observed_total_output_path) and not fs.is_empty(self.observed_total_output_path)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_output(self):
        return self.observed_bulge_output_path is not None and fs.is_directory(self.observed_bulge_output_path) and not fs.is_empty(self.observed_bulge_output_path)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_output(self):
        return self.observed_disk_output_path is not None and fs.is_directory(self.observed_disk_output_path) and not fs.is_empty(self.observed_disk_output_path)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_output(self):
        return self.observed_old_output_path is not None and fs.is_directory(self.observed_old_output_path) and not fs.is_empty(self.observed_old_output_path)

    # -----------------------------------------------------------------

    @property
    def has_observed_young_output(self):
        return self.observed_young_output_path is not None and fs.is_directory(self.observed_young_output_path) and not fs.is_empty(self.observed_young_output_path)

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_output(self):
        return self.observed_sfr_output_path is not None and fs.is_directory(self.observed_sfr_output_path) and not fs.is_empty(self.observed_sfr_output_path)

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_output(self):
        return self.observed_unevolved_output_path is not None and fs.is_directory(self.observed_unevolved_output_path) and not fs.is_empty(self.observed_unevolved_output_path)

    # -----------------------------------------------------------------

    @property
    def has_observed_extra_output(self):
        return self.observed_extra_output_path is not None and fs.is_directory(self.observed_extra_output_path) and not fs.is_empty(self.observed_extra_output_path)

    # -----------------------------------------------------------------
    # FILTERS & WAVELENGTHS
    # -----------------------------------------------------------------

    @lazyproperty
    def i1_filter(self):
        return parse_filter("IRAC I1")

    # -----------------------------------------------------------------

    @property
    def i1_wavelength(self):
        return self.i1_filter.wavelength

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_filter(self):
        return parse_filter("GALEX FUV")

    # -----------------------------------------------------------------

    @property
    def fuv_wavelength(self):
        return self.fuv_filter.wavelength

    # -----------------------------------------------------------------

    @lazyproperty
    def mips24_filter(self):
        return parse_filter("MIPS 24")

    # -----------------------------------------------------------------

    @property
    def mips24_wavelength(self):
        return self.mips24_filter.wavelength

    # -----------------------------------------------------------------
    # PARAMETERS
    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):
        return modeling_parameter_labels

    # -----------------------------------------------------------------

    @lazyproperty
    def other_parameter_labels(self):
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
        return create_subdict(self.parameter_values, self.free_parameter_labels)

    # -----------------------------------------------------------------

    @lazyproperty
    def other_parameter_values(self):
        return create_subdict(self.parameter_values, self.other_parameter_labels)

    # -----------------------------------------------------------------

    @property
    def metallicity(self):
        return self.parameter_values[metallicity_name]

    # -----------------------------------------------------------------

    @property
    def sfr_compactness(self):
        return self.parameter_values[sfr_compactness_name]

    # -----------------------------------------------------------------

    @property
    def sfr_pressure(self):
        return self.parameter_values[sfr_pressure_name]

    # -----------------------------------------------------------------

    @property
    def sfr_covering_factor(self):
        return self.parameter_values[sfr_covering_name]

    # -----------------------------------------------------------------

    @property
    def mappings_parameters_path(self):
        return fs.join(self.sfr_mappings_path, "parameters.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_parameters(self):
        return fs.is_file(self.mappings_parameters_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings_parameters(self):
        return load_dict(self.mappings_parameters_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr(self):

        """
        This function derives the SFR (in Msun / year) from the FUV luminosity of the model and the intrinsic MAPPINGS SED
        :return:
        """

        # Has MAPPINGS parameters file?
        if self.has_mappings_parameters: return self.mappings_parameters[sfr_parameter]

        # No parameters file yet
        else:

            # Calculate the SFR
            sfr = Mappings.sfr_for_luminosity(self.metallicity, self.sfr_compactness, self.sfr_pressure, self.sfr_covering_factor, self.intrinsic_fuv_luminosity_sfr, self.fuv_wavelength)

            # Create parameters dict
            parameters = OrderedDict()
            parameters[metallicity_parameter] = self.metallicity
            parameters[compactness_parameter] = self.sfr_compactness
            parameters[pressure_parameter] = self.sfr_pressure
            parameters[covering_factor_parameter] = self.sfr_covering_factor
            parameters[sfr_parameter] = sfr

            # Write the dictionary
            write_dict(parameters, self.mappings_parameters_path)

            # Return the star formation rate
            return sfr

    # -----------------------------------------------------------------

    @property
    def has_sfr(self):
        return self.has_mappings

    # -----------------------------------------------------------------

    @property
    def has_mappings(self):
        return True  # should always be able to be created

    # -----------------------------------------------------------------

    @property
    def mappings_sed_path(self):
        return fs.join(self.sfr_mappings_path, "sed.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_sed(self):
        return fs.is_file(self.mappings_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings_sed(self):
        return SED.from_file(self.mappings_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings(self):
        mappings = Mappings(self.metallicity, self.sfr_compactness, self.sfr_pressure, self.sfr_covering_factor, self.sfr)
        if self.has_mappings_sed: mappings.sed = self.mappings_sed
        else: mappings.sed.saveto(self.mappings_sed_path)
        return mappings

    # -----------------------------------------------------------------

    @property
    def mappings_transparent_sed_path(self):
        return fs.join(self.sfr_mappings_path, "sed_transparent.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_transparent_sed(self):
        return fs.is_file(self.mappings_transparent_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings_transparent_sed(self):
        return SED.from_file(self.mappings_transparent_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings_transparent(self):
        # No metallicity: no dust (?) -> No, doesn't change the curve in the right direction
        # Minimal compactness: minimal absorption
        # Minimal covering factor: minimal absorption
        mappings = Mappings(self.metallicity, 0., self.sfr_pressure, 0., self.sfr)
        if self.has_mappings_transparent_sed: mappings.sed = self.mappings_transparent_sed
        else: mappings.sed.saveto(self.mappings_transparent_sed_path)
        return mappings

    # -----------------------------------------------------------------

    @property
    def mappings_transparent_parameters_path(self):
        return fs.join(self.sfr_mappings_path, "parameters_transparent.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_transparent_parameters(self):
        return fs.is_file(self.mappings_transparent_parameters_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings_transparent_parameters(self):
        
        """
        This function ...
        :return: 
        """

        # Has parameters file?
        if self.has_mappings_transparent_parameters: return load_dict(self.mappings_transparent_parameters_path, ordered=True)

        # No file yet
        else:

            # Set neutral luminosity and flux density
            luminosity = self.mappings_transparent.sed.photometry_at(self.fuv_wavelength, unit=self.intrinsic_fuv_luminosity_sfr.unit, interpolate=True)
            neutral_luminosity = luminosity.to("Lsun", density=True, density_strict=True, wavelength=self.fuv_wavelength)
            fluxdensity = luminosity.to("Jy", wavelength=self.fuv_wavelength, distance=self.distance)

            # Create parameters dict
            parameters = OrderedDict()
            parameters[metallicity_parameter] = self.metallicity
            parameters[compactness_parameter] = 0.
            parameters[pressure_parameter] = self.sfr_pressure
            parameters[covering_factor_parameter] = 0.
            parameters[sfr_parameter] = self.sfr

            parameters["luminosity"] = luminosity
            parameters["neutral_luminosity"] = neutral_luminosity
            parameters["fluxdensity"] = fluxdensity

            # Write the dictionary
            write_dict(parameters, self.mappings_transparent_parameters_path)

            # Return the parameters
            return parameters

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_transparent_sfr(self):
        return self.mappings_transparent_parameters["luminosity"]

    # -----------------------------------------------------------------

    def find_internal_mappings_attenuation(self, min_fuv_attenuation_fraction=0.1, fuv_attenuation_fraction_step=0.01):

        """
        # NORMALIZE THE ATTENUATION CURVE AT FUV TO ONLY ATTENUATE 10% OF THE LIGHT: AN UNDERESTIMATION,
        # SO WE CAN DRIVE UP THE ATTENUATION AND KEEP DOING THAT UNTIL WE FIND A MATCH WITH THE EMITTED DUST ENERGY
        This function ...
        :param min_fuv_attenuation_fraction:
        :param fuv_attenuation_fraction_step:
        :return:
        """

        # Bolometric luminosities
        bolometric_luminosity_unit = parse_quantity("W")

        # 'Observed': after internal attenuation
        observed_sed = self.intrinsic_sfr_stellar_sed
        stellar_sed = self.sfr_simulations.intrinsic_sed
        stellar_luminosity = stellar_sed.integrate().to(bolometric_luminosity_unit, distance=self.distance)

        # Dust
        internal_dust_sed = self.intrinsic_sfr_dust_sed
        dust_luminosity_internal = internal_dust_sed.integrate().to(bolometric_luminosity_unit, distance=self.distance)
        dust_fraction_internal = dust_luminosity_internal.value / stellar_luminosity.value

        # Get attenuation curve
        curve = MappingsAttenuationCurve(from_astropy=False)

        # Loop over the values of attenuation fraction
        att_fraction = min_fuv_attenuation_fraction
        while True:

            # Debugging
            log.debug("Testing a FUV attenuation fraction of " + str(att_fraction) + " ...")

            # Get observed over transparent FUV luminosity
            obs_over_trans = 1. / (att_fraction + 1.)
            attenuation = -2.5 * np.log10(obs_over_trans)

            # Normalize the attenuation curve
            curve.normalize_at(self.fuv_wavelength, attenuation)
            #transmissions = 10 ** (-0.4 * curve.y_array)

            # Correct for attenuation
            transparent_sed = correct_sed_for_attenuation(observed_sed, curve)

            # Get absorbed luminosity under this attenuation curve
            absorbed_sed = transparent_sed - observed_sed
            absorbed_lum = absorbed_sed.integrate().to(bolometric_luminosity_unit, distance=self.distance)
            #print(absorbed_lum)
            absorption_fraction = absorbed_lum.value / stellar_luminosity.value

            # Debugging
            log.debug("The bolometric absorption fraction is " + str(absorption_fraction) + " (dust emission fraction is " + str(dust_fraction_internal) + ")")

            # Return the FUV attenuation
            if absorbed_lum >= dust_luminosity_internal: return attenuation

            # Increase the attenuation fraction
            att_fraction += fuv_attenuation_fraction_step

    # -----------------------------------------------------------------

    @property
    def attenuation_curve_sfr_internal_path(self):
        return fs.join(self.sfr_mappings_path, "internal_attenuation.dat")

    # -----------------------------------------------------------------

    @property
    def has_attenuation_curve_sfr_internal(self):
        return fs.is_file(self.attenuation_curve_sfr_internal_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def attenuation_curve_sfr_internal(self):
        if not self.has_attenuation_curve_sfr_internal:
            fuv_attenuation = self.find_internal_mappings_attenuation()
            attenuation_curve = MappingsAttenuationCurve(wavelength=self.fuv_wavelength, attenuation=fuv_attenuation)
            attenuation_curve.saveto(self.attenuation_curve_sfr_internal_path)
            return attenuation_curve
        else: return AttenuationCurve.from_file(self.attenuation_curve_sfr_internal_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_attenuation_sfr_internal(self):
        return self.attenuation_curve_sfr_internal.attenuation_at(self.fuv_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def unattenuated_mappings_sed(self):
        return correct_sed_for_attenuation(self.intrinsic_sfr_stellar_sed, self.attenuation_curve_sfr_internal)
        
    # -----------------------------------------------------------------

    @lazyproperty
    def attenuation_sed_sfr_internal(self):
        return self.unattenuated_mappings_sed - self.intrinsic_sfr_stellar_sed

    # -----------------------------------------------------------------

    @property
    def normalized_mappings_sed_path(self):
        return fs.join(self.sfr_mappings_path, "normalized_sed.dat")

    # -----------------------------------------------------------------

    @property
    def has_normalized_mappings_sed(self):
        return fs.is_file(self.normalized_mappings_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def normalized_mappings_sed(self):
        return SED.from_file(self.normalized_mappings_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def normalized_mappings(self):
        mappings = Mappings(self.metallicity, self.sfr_compactness, self.sfr_pressure, self.sfr_covering_factor)
        if self.has_normalized_mappings_sed: mappings.sed = self.normalized_mappings_sed
        else: mappings.sed.saveto(self.normalized_mappings_sed_path)
        return mappings

    # -----------------------------------------------------------------
    # TOTAL SIMULATIONS
    # -----------------------------------------------------------------

    @property
    def has_observed_total_bolometric_luminosity(self):
        return self.total_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_total_bolometric_luminosity(self):
        return self.total_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_total_bolometric_luminosity(self):
        return self.total_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_total_bolometric_luminosity(self):
        return self.total_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_stellar_luminosity(self):
        return self.total_simulations.observed_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_luminosity(self):
        return self.total_simulations.has_observed_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_luminosity(self):
        #return self.intrinsic_bolometric_luminosity_old + self.intrinsic_bolometric_luminosity_young + self.intrinsic_bolometric_luminosity_sfr
        return self.total_simulations.intrinsic_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_luminosity(self):
        #return self.has_intrinsic_bolometric_luminosity_old and self.has_intrinsic_bolometric_luminosity_young and self.has_intrinsic_bolometric_luminosity_sfr
        return self.total_simulations.has_intrinsic_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation(self):
        return self.total_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation(self):
        return self.total_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_total_sed(self):
        return self.total_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_total_sed(self):
        return self.total_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_total_sed(self):
        return self.total_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_total_sed(self):
        return self.total_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_total_stellar_sed(self):

        """
        This function ...
        :return:
        """

        # Get the SED
        sed = self.total_simulations.observed_stellar_sed

        # We need to add the internal dust emission part of the the SFR
        #if self.total_simulations.observed_stellar_sed_needs_reprocessed_internal_part: return sed + self.intrinsic_dust_sed_sfr
        #else: return sed

        return sed

    # -----------------------------------------------------------------

    @property
    def has_observed_total_stellar_sed(self):
        return self.total_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_total_dust_sed(self):
        return self.observed_total_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_observed_total_dust_sed(self):
        return self.has_observed_total_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def observed_total_dust_sed_earth(self):
        return self.total_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_total_dust_sed_earth(self):
        return self.total_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_total_dust_sed_faceon(self):
        return self.total_simulations.faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_total_dust_sed_faceon(self):
        return self.total_simulations.has_faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_total_dust_sed_edgeon(self):
        return self.total_simulations.edgeon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_total_dust_sed_edgeon(self):
        return self.total_simulations.has_edgeon_observed_dust_sed

    # -----------------------------------------------------------------
    # OLD BULGE SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_old_bulge(self):
        # don't interpolate, wavelength grid is expected to contain the I1 wavelength
        return self.bulge_simulations.observed_photometry_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_old_bulge(self):
        return self.bulge_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_i1_luminosity_old_bulge(self):
        return self.definition.bulge_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_i1_luminosity_old_bulge(self):
        return True # should be defined in definition

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_old_bulge(self):
        return self.bulge_simulations.observed_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_old_bulge(self):
        return self.bulge_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_old_bulge(self):
        return self.bulge_simulations.intrinsic_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_old_bulge(self):
        return self.bulge_simulations.has_intrinsic_photometry

    # -----------------------------------------------------------------

    @property
    def attenuation_curve_old_bulge(self):
        return self.bulge_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_attenuation_old_bulge(self):
        return self.bulge_simulations.attenuation_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_i1_attenuation_old_bulge(self):
        return self.bulge_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_old_bulge(self):
        return self.bulge_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_old_bulge(self):
        return self.bulge_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_old_bulge(self):
        return self.bulge_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_old_bulge(self):
        return self.bulge_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_old_bulge(self):
        return self.bulge_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_old_bulge(self):
        return self.bulge_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_old_bulge(self):
        return self.bulge_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_old_bulge(self):
        return self.bulge_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_sed(self):
        return self.bulge_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_bulge_sed(self):
        return self.bulge_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_old_bulge_sed(self):
        return self.bulge_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_old_bulge_sed(self):
        return self.bulge_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_bulge_stellar_sed(self):
        return self.bulge_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_stellar_sed(self):
        return self.bulge_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_bulge_dust_sed(self):
        return self.bulge_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_dust_sed(self):
        return self.bulge_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------
    # OLD DISK SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_old_disk(self):
        # don't interpolate, wavelength grid is expected to contain the I1 wavelength
        return self.disk_simulations.observed_photometry_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_old_disk(self):
        return self.disk_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_i1_luminosity_old_disk(self):
        return self.definition.old_stars_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_i1_luminosity_old_disk(self):
        return True # should be defined in definition

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_old_disk(self):
        return self.disk_simulations.observed_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_old_disk(self):
        return self.disk_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_old_disk(self):
        return self.disk_simulations.intrinsic_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_old_disk(self):
        return self.disk_simulations.has_intrinsic_photometry

    # -----------------------------------------------------------------

    @property
    def attenuation_curve_old_disk(self):
        return self.disk_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_attenuation_old_disk(self):
        return self.disk_simulations.attenuation_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_i1_attenuation_old_disk(self):
        return self.disk_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_old_disk(self):
        return self.disk_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_old_disk(self):
        return self.disk_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_old_disk(self):
        return self.disk_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_old_disk(self):
        return self.disk_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_old_disk(self):
        return self.disk_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_old_disk(self):
        return self.disk_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_old_disk(self):
        return self.disk_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_old_disk(self):
        return self.disk_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_sed(self):
        return self.disk_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_disk_sed(self):
        return self.disk_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_old_disk_sed(self):
        return self.disk_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_old_disk_sed(self):
        return self.disk_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_disk_stellar_sed(self):
        return self.disk_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_stellar_sed(self):
        return self.disk_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_disk_dust_sed(self):
        return self.disk_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_dust_sed(self):
        return self.disk_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------
    # OLD SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_old(self):
        # don't interpolate, wavelength grid is expected to contain the I1 wavelength
        return self.old_simulations.observed_photometry_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_old(self):
        return self.old_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_i1_luminosity_old(self):
        return self.old_simulations.intrinsic_photometry_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_i1_luminosity_old(self):
        return self.old_simulations.has_intrinsic_photometry

    # -----------------------------------------------------------------

    @property
    def attenuation_curve_old(self):
        return self.old_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_attenuation_old(self):
        return self.old_simulations.attenuation_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_i1_attenuation_old(self):
        return self.old_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def attenuation_i1_old(self):
        return self.i1_attenuation_old

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_old(self):
        return self.old_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_old(self):
        return self.old_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_old(self):
        return self.old_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_old(self):
        return self.old_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_old(self):
        return self.old_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_old(self):
        return self.old_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_old(self):
        return self.old_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_old(self):
        return self.old_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_old_sed(self):
        return self.old_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_sed(self):
        return self.old_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_old_sed(self):
        return self.old_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_old_sed(self):
        return self.old_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_stellar_sed(self):
        return self.old_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_stellar_sed(self):
        return self.old_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_dust_sed(self):
        return self.observed_old_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_observed_old_dust_sed(self):
        return self.has_observed_old_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def observed_old_dust_sed_earth(self):
        return self.old_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_dust_sed_earth(self):
        return self.old_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_dust_sed_faceon(self):
        return self.old_simulations.faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_dust_sed_faceon(self):
        return self.old_simulations.has_faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_dust_sed_edgeon(self):
        return self.old_simulations.edgeon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_dust_sed_edgeon(self):
        return self.old_simulations.has_edgeon_observed_dust_sed

    # -----------------------------------------------------------------
    # YOUNG SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_young(self):
        # don't interpolate, wavelength grid is expected to contain the FUV wavelength
        return self.young_simulations.observed_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_young(self):
        return self.young_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @property
    def intrinsic_fuv_luminosity_young(self):
        return self.parameter_values[fuv_young_name]

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_young(self):
        return True # part of free parameters

    # -----------------------------------------------------------------

    @property
    def attenuation_curve_young(self):
        return self.young_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_attenuation_young(self):
        return self.young_simulations.attenuation_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_fuv_attenuation_young(self):
        return self.young_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_young(self):
        return self.young_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_young(self):
        return self.young_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_young(self):
        return self.young_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_young(self):
        return self.young_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_young(self):
        return self.young_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_young(self):
        return self.young_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_young(self):
        return self.young_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_young(self):
        return self.young_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_young_sed(self):
        return self.young_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_young_sed(self):
        return self.young_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_young_sed(self):
        return self.young_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_young_sed(self):
        return self.young_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_young_stellar_sed(self):
        return self.young_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_young_stellar_sed(self):
        return self.young_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_young_dust_sed(self):
        return self.observed_young_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_observed_young_dust_sed(self):
        return self.has_observed_young_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def observed_young_dust_sed_earth(self):
        return self.young_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_young_dust_sed_earth(self):
        return self.young_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_young_dust_sed_faceon(self):
        return self.young_simulations.faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_young_dust_sed_faceon(self):
        return self.young_simulations.has_faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_young_dust_sed_edgeon(self):
        return self.young_simulations.edgeon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_young_dust_sed_edgeon(self):
        return self.young_simulations.has_edgeon_observed_dust_sed


    # -----------------------------------------------------------------
    # EXTRA SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_extra(self):
        # don't interpolate, wavelength grid is expected to contain the FUV wavelength
        return self.extra_simulations.observed_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_extra(self):
        return self.extra_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @property
    def attenuation_curve_extra(self):
        return self.extra_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_extra(self):
        return self.extra_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_extra(self):
        return self.extra_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_extra(self):
        return self.extra_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_extra(self):
        return self.extra_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_extra(self):
        return self.extra_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_extra(self):
        return self.extra_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_extra(self):
        return self.extra_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_extra(self):
        return self.extra_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_extra_sed(self):
        return self.extra_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_extra_sed(self):
        return self.extra_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_extra_sed(self):
        return self.extra_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_extra_sed(self):
        return self.extra_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_extra_stellar_sed(self):
        return self.extra_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_extra_stellar_sed(self):
        return self.extra_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_extra_dust_sed(self):
        return self.observed_extra_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_observed_extra_dust_sed(self):
        return self.has_observed_extra_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def observed_extra_dust_sed_earth(self):
        return self.extra_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_extra_dust_sed_earth(self):
        return self.extra_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_extra_dust_sed_faceon(self):
        return self.extra_simulations.faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_extra_dust_sed_faceon(self):
        return self.extra_simulations.has_faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_extra_dust_sed_edgeon(self):
        return self.extra_simulations.edgeon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_extra_dust_sed_edgeon(self):
        return self.extra_simulations.has_edgeon_observed_dust_sed

    # -----------------------------------------------------------------
    # SFR
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_sfr(self):
        # don't interpolate, wavelength grid is expected to contain the FUV wavelength
        return self.sfr_simulations.observed_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_sfr(self):
        return self.sfr_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @property
    def intrinsic_fuv_luminosity_sfr(self):
        return self.parameter_values[fuv_ionizing_name]

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_sfr(self):
        return True # free parameter

    # -----------------------------------------------------------------

    @property
    def attenuation_curve_sfr(self):
        return self.sfr_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_attenuation_sfr(self):
        return self.sfr_simulations.attenuation_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_fuv_attenuation_sfr(self):
        return self.sfr_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_sfr(self):
        return self.sfr_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_sfr(self):
        return self.sfr_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_sfr(self):
        return self.sfr_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_sfr(self):
        return self.sfr_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_sfr(self):
        return self.sfr_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_sfr(self):
        return self.sfr_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    # ONLY FOR SFR
    @property
    def intrinsic_dust_luminosity_sfr(self):
        return self.sfr_simulations.intrinsic_dust_luminosity

    # -----------------------------------------------------------------

    # ONLY FOR SFR
    @property
    def has_intrinsic_dust_luminosity_sfr(self):
        return self.sfr_simulations.has_intrinsic_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_sfr(self):
        return self.sfr_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_sfr(self):
        return self.sfr_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_sed(self):
        return self.sfr_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_sfr_sed(self):
        return self.sfr_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sfr_sed(self):
        return self.sfr_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_sfr_sed(self):
        return self.sfr_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def fit_mappings_from_wavelength(self):
        return parse_quantity("0.1 micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def extrapolate_mappings_from_wavelength(self):
        return parse_quantity("5 micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_transparent_sfr_stellar_sed(self):
        return self.transparent_sfr_sed.extrapolated_from(self.extrapolate_mappings_from_wavelength, self.fit_mappings_from_wavelength, xlog=True, ylog=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sfr_stellar_sed(self):
        return self.intrinsic_sfr_sed.extrapolated_from(self.extrapolate_mappings_from_wavelength, self.fit_mappings_from_wavelength, xlog=True, ylog=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sfr_dust_sed(self):
        #from pts.core.plot.sed import plot_seds
        #plot_seds({"intrinsic": self.intrinsic_sfr_sed, "intrinsic_stellar": self.intrinsic_sfr_stellar_sed})
        return self.intrinsic_sfr_sed - self.intrinsic_sfr_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_sfr_stellar_sed(self):
        return self.sfr_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_stellar_sed(self):
        return self.sfr_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_sfr_dust_sed(self):
        return self.observed_sfr_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_dust_sed(self):
        return self.has_observed_sfr_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def observed_sfr_dust_sed_earth(self):
        return self.sfr_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_dust_sed_earth(self):
        return self.sfr_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_sfr_dust_sed_faceon(self):
        return self.sfr_simulations.faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_dust_sed_faceon(self):
        return self.sfr_simulations.has_faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_sfr_dust_sed_edgeon(self):
        return self.sfr_simulations.edgeon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_dust_sed_edgeon(self):
        return self.sfr_simulations.has_edgeon_observed_sed

    # -----------------------------------------------------------------

    @property
    def diffuse_dust_mass(self):
        return self.parameter_values[dust_mass_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_mass(self):
        return self.mappings.dust_mass

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass(self):
        return self.has_mappings

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_mass(self):
        return self.diffuse_dust_mass + self.sfr_dust_mass

    # -----------------------------------------------------------------

    @property
    def total_dust_mass(self):
        return self.dust_mass

    # -----------------------------------------------------------------

    @property
    def has_total_dust_mass(self):
        return self.has_sfr_dust_mass

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_stellar_mass(self):
        return self.mappings.stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_sfr_stellar_mass(self):
        #return self.has_mappings
        return False # returns NotImplementedError in Mappings: we don't know the conversion yet between Mappings parameters and stellar mass!

    # -----------------------------------------------------------------

    @property
    def observed_stellar_luminosity_sfr(self):
        return self.sfr_simulations.observed_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_luminosity_sfr(self):
        return self.sfr_simulations.has_observed_stellar_luminosity

    # -----------------------------------------------------------------

    # ONLY FOR SFR (NORMALLY INTRINSIC STELLAR = INTRINSIC BOLOMETRIC)
    @property
    def intrinsic_stellar_luminosity_sfr(self):
        return self.sfr_simulations.intrinsic_stellar_luminosity

    # -----------------------------------------------------------------

    # ONLY FOR SFR (NORMALLY INTRINSIC STELLAR = INTRINSIC BOLOMETRIC)
    @property
    def has_intrinsic_stellar_luminosity_sfr(self):
        return self.sfr_simulations.has_intrinsic_stellar_luminosity

    # -----------------------------------------------------------------
    # UNEVOLVED
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_unevolved(self):
        # don't interpolate, wavelength grid is expected to contain the FUV wavelength
        return self.unevolved_simulations.observed_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_unevolved(self):
        return self.unevolved_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_unevolved(self):
        #self.intrinsic_fuv_luminosity_young + self.intrinsic_fuv_luminosity_sfr
        return self.unevolved_simulations.intrinsic_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_unevolved(self):
        #return self.has_intrinsic_fuv_luminosity_young and self.has_intrinsic_fuv_luminosity_sfr
        return self.unevolved_simulations.has_intrinsic_photometry

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_unevolved(self):
        return self.unevolved_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_unevolved(self):
        return self.unevolved_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_unevolved(self):
        return self.unevolved_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_unevolved(self):
        return self.unevolved_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def attenuation_curve_unevolved(self):
        return self.unevolved_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_attenuation_unevolved(self):
        return self.unevolved_simulations.attenuation_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_fuv_attenuation_unevolved(self):
        return self.unevolved_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_unevolved(self):
        return self.unevolved_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_unevolved(self):
        return self.unevolved_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_unevolved(self):
        return self.unevolved_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_unevolved(self):
        return self.unevolved_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_sed(self):
        return self.unevolved_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_unevolved_sed(self):
        return self.unevolved_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_unevolved_sed(self):
        return self.unevolved_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_unevolved_sed(self):
        return self.unevolved_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_unevolved_stellar_sed(self):
        return self.unevolved_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_stellar_sed(self):
        return self.unevolved_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_unevolved_dust_sed(self):
        return self.observed_unevolved_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_dust_sed(self):
        return self.has_observed_unevolved_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def observed_unevolved_dust_sed_earth(self):
        return self.unevolved_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_dust_sed_earth(self):
        return self.unevolved_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_unevolved_dust_sed_faceon(self):
        return self.unevolved_simulations.faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_dust_sed_faceon(self):
        return self.unevolved_simulations.has_faceon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def observed_unevolved_dust_sed_edgeon(self):
        return self.unevolved_simulations.edgeon_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_dust_sed_edgeon(self):
        return self.unevolved_simulations.has_edgeon_observed_dust_sed

    # -----------------------------------------------------------------
    # PROJECTION PROPERTIES
    # -----------------------------------------------------------------

    @property
    def distance(self):
        return self.definition.distance

    # -----------------------------------------------------------------

    @property
    def inclination(self):
        return self.definition.inclination

    # -----------------------------------------------------------------

    @property
    def position_angle(self):
        return self.definition.position_angle

    # -----------------------------------------------------------------
    # INSTRUMENTS
    # -----------------------------------------------------------------

    @lazyproperty
    def sed_instrument(self):
        return SEDInstrument.from_properties(self.distance, self.inclination, self.position_angle)

    # -----------------------------------------------------------------

    @lazyproperty
    def full_sed_instrument(self):
        return FullSEDInstrument.from_properties(self.distance, self.inclination, self.position_angle)

    # -----------------------------------------------------------------
    # MODEL DEFINITION
    # -----------------------------------------------------------------

    @property
    def definition_path(self):
        return self.definition.path

    # -----------------------------------------------------------------

    @property
    def old_bulge_path(self):
        return self.definition.bulge_component_path

    # -----------------------------------------------------------------

    @property
    def old_disk_path(self):
        return self.definition.old_stars_component_path

    # -----------------------------------------------------------------

    @property
    def young_path(self):
        return self.definition.young_stars_component_path

    # -----------------------------------------------------------------

    @property
    def sfr_path(self):
        return self.definition.ionizing_stars_component_path

    # -----------------------------------------------------------------

    @property
    def extra_path(self):
        return self.definition.stellar_paths["extra"]

    # -----------------------------------------------------------------

    @property
    def dust_path(self):
        return self.definition.dust_component_path

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_sed_path(self):
        return fs.create_directory_in(self.old_bulge_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_projections_path(self):
        return fs.create_directory_in(self.old_bulge_path, projections_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_sed_path(self):
        return fs.create_directory_in(self.old_disk_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_projections_path(self):
        return fs.create_directory_in(self.old_disk_path, projections_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_sed_path(self):
        return fs.create_directory_in(self.young_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_projections_path(self):
        return fs.create_directory_in(self.young_path, projections_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_sed_path(self):
        return fs.create_directory_in(self.extra_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_projections_path(self):
        return fs.create_directory_in(self.extra_path, projections_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_sed_path(self):
        return fs.create_directory_in(self.sfr_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_projections_path(self):
        return fs.create_directory_in(self.sfr_path, projections_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_mappings_path(self):
        return fs.create_directory_in(self.sfr_path, mappings_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def transparent_sfr_sed_path(self):
        return fs.create_directory_in(self.sfr_path, transparent_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_path(self):
        return fs.create_directory_in(self.dust_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_projections_path(self):
        return fs.create_directory_in(self.dust_path, projections_dirname)

    # -----------------------------------------------------------------
    # INTRINSIC SED SKI FILES
    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_sed_ski_path(self):
        return fs.join(self.old_disk_sed_path, disk_simulation_name + ".ski")

    # -----------------------------------------------------------------

    @lazyproperty
    def young_sed_ski_path(self):
        return fs.join(self.young_sed_path, young_simulation_name + ".ski")

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_sed_ski_path(self):
        return fs.join(self.sfr_sed_path, sfr_simulation_name + ".ski")

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_sed_ski_path(self):
        return fs.join(self.extra_sed_path, extra_simulation_name + ".ski")

    # -----------------------------------------------------------------

    @property
    def has_old_disk_sed_skifile(self):
        return fs.is_file(self.old_disk_sed_ski_path)

    # -----------------------------------------------------------------

    @property
    def has_young_sed_skifile(self):
        return fs.is_file(self.young_sed_ski_path)

    # -----------------------------------------------------------------

    @property
    def has_sfr_sed_skifile(self):
        return fs.is_file(self.sfr_sed_ski_path)

    # -----------------------------------------------------------------

    @property
    def has_extra_sed_skifile(self):
        return fs.is_file(self.extra_sed_ski_path)

    # -----------------------------------------------------------------
    # COMPONENT MODELS
    #   OLD BULGE
    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_component(self):
        return self.definition.load_bulge_component()

    # -----------------------------------------------------------------

    @property
    def old_bulge_model(self):
        return self.old_bulge_component.model

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_radial_effective_radius(self):
        return self.old_bulge_model.effective_radius

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_vertical_effective_radius(self):
        return self.old_bulge_model.effective_radius * self.old_bulge_model.z_flattening

    # -----------------------------------------------------------------

    @property
    def old_bulge_scaleheight(self):
        return self.old_bulge_vertical_effective_radius

    # -----------------------------------------------------------------
    #   OLD DISK
    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_component(self):
        return self.definition.load_old_stars_component()

    # -----------------------------------------------------------------

    @property
    def old_disk_deprojection(self):
        return self.old_disk_component.deprojection

    # -----------------------------------------------------------------

    @property
    def old_disk_scaleheight(self):
        return self.old_disk_deprojection.scale_height

    # -----------------------------------------------------------------
    #   YOUNG DISK
    # -----------------------------------------------------------------

    @lazyproperty
    def young_component(self):
        return self.definition.load_young_stars_component()

    # -----------------------------------------------------------------

    @property
    def young_deprojection(self):
        return self.young_component.deprojection

    # -----------------------------------------------------------------

    @property
    def young_scaleheight(self):
        return self.young_deprojection.scale_height

    # -----------------------------------------------------------------
    #   SFR DISK
    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_component(self):
        return self.definition.load_ionizing_stars_component()

    # -----------------------------------------------------------------

    @property
    def sfr_deprojection(self):
        return self.sfr_component.deprojection

    # -----------------------------------------------------------------

    @property
    def sfr_scaleheight(self):
        return self.sfr_deprojection.scale_height

    # -----------------------------------------------------------------

    @lazyproperty
    def transparent_sfr_component(self):
        component = deepcopy(self.sfr_component)
        component.parameters_path = None
        component.parameters.fluxdensity = self.mappings_transparent_parameters["fluxdensity"]
        component.parameters.luminosity = self.mappings_transparent_parameters["luminosity"]
        component.parameters.title = "Ionizing stars (transparent)"
        component.parameters.pressure = self.mappings_transparent_parameters["pressure"]
        component.parameters.metallicity = self.mappings_transparent_parameters["metallicity"]
        component.parameters.covering_factor = self.mappings_transparent_parameters["covering_factor"]
        component.parameters.compactness = self.mappings_transparent_parameters["compactness"]
        component.parameters.neutral_luminosity = self.mappings_transparent_parameters["neutral_luminosity"]
        return component

    # -----------------------------------------------------------------
    #   EXTRA STELLAR COMPONENT
    # -----------------------------------------------------------------

    @lazyproperty
    def extra_component(self):
        return self.definition.load_stellar_component("extra")

    # -----------------------------------------------------------------

    @property
    def extra_model(self):
        return self.extra_component.model

    # -----------------------------------------------------------------
    #   DUST DISK
    # -----------------------------------------------------------------

    @lazyproperty
    def dust_component(self):
        return self.definition.load_dust_disk_component()

    # -----------------------------------------------------------------

    @property
    def dust_deprojection(self):
        return self.dust_component.deprojection

    # -----------------------------------------------------------------

    @property
    def dust_scaleheight(self):
        return self.dust_deprojection.scale_height

    # -----------------------------------------------------------------
    # WAVELENGTH GRID
    # -----------------------------------------------------------------

    @property
    def has_wavelength_grid(self):
        return self.wavelength_grid is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths(self):
        return self.wavelength_grid.wavelengths(unit="micron", add_unit=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths_micron(self):
        return self.wavelength_grid.wavelengths(unit="micron", asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_deltas(self):
        return self.wavelength_grid.deltas(unit="micron", add_unit=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_deltas_micron(self):
        return self.wavelength_grid.deltas(unit="micron", asarray=True)

    # -----------------------------------------------------------------

    @property
    def has_wavelengths_directory(self):
        return fs.contains_directory(self.definition_path, "wavelengths")

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths_path(self):
        return fs.create_directory_in(self.definition_path, "wavelengths")

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid_path(self):
        return fs.join(self.wavelengths_path, "grid.txt")

    # -----------------------------------------------------------------
    # TOTAL CUBES
    # -----------------------------------------------------------------

    @property
    def total_bolometric_luminosity_cube_earth(self):
        return self.total_simulations.observed_cube

    # -----------------------------------------------------------------

    @property
    def has_total_bolometric_luminosity_cube_earth(self):
        return self.total_simulations.has_observed_cube

    # -----------------------------------------------------------------

    @property
    def total_bolometric_luminosity_cube_faceon(self):
        return self.total_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_total_bolometric_luminosity_cube_faceon(self):
        return self.total_simulations.has_faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def total_bolometric_luminosity_cube_edgeon(self):
        return self.total_simulations.edgeon_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_total_bolometric_luminosity_cube_edgeon(self):
        return self.total_simulations.has_edgeon_observed_cube

    # -----------------------------------------------------------------

    @property
    def total_intrinsic_stellar_luminosity_cube_earth(self):
        return self.total_simulations.intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def has_total_intrinsic_stellar_luminosity_cube_earth(self):
        return self.total_simulations.has_intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def total_intrinsic_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.faceon_intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def has_total_intrinsic_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.has_faceon_intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def total_intrinsic_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.edgeon_intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def has_total_intrinsic_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.has_edgeon_intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def total_observed_stellar_luminosity_cube_earth(self):
        return self.total_simulations.observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_total_observed_stellar_luminosity_cube_earth(self):
        return self.total_simulations.has_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def total_observed_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_total_observed_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.has_faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def total_observed_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_total_observed_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.has_edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def total_diffuse_dust_luminosity_cube_earth(self):
        return self.total_simulations.observed_diffuse_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_total_diffuse_dust_luminosity_cube_earth(self):
        return self.total_simulations.has_observed_diffuse_dust_cube

    # -----------------------------------------------------------------

    @property
    def total_diffuse_dust_luminosity_cube_faceon(self):
        return self.total_simulations.faceon_observed_diffuse_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_total_diffuse_dust_luminosity_cube_faceon(self):
        return self.total_simulations.has_faceon_observed_diffuse_dust_cube

    # -----------------------------------------------------------------

    @property
    def total_diffuse_dust_luminosity_cube_edgeon(self):
        return self.total_simulations.edgeon_observed_diffuse_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_total_diffuse_dust_luminosity_cube_edgeon(self):
        return self.total_simulations.has_edgeon_observed_diffuse_dust_cube

    # -----------------------------------------------------------------

    @property
    def total_dust_luminosity_cube_earth(self):
        return self.total_simulations.observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_total_dust_luminosity_cube_earth(self):
        return self.total_simulations.has_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def total_dust_luminosity_cube_faceon(self):
        return self.total_simulations.faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_total_dust_luminosity_cube_faceon(self):
        return self.total_simulations.has_faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def total_dust_luminosity_cube_edgeon(self):
        return self.total_simulations.edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_total_dust_luminosity_cube_edgeon(self):
        return self.total_simulations.has_edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def total_scattered_stellar_luminosity_cube_earth(self):
        return self.total_simulations.observed_cube_scattered

    # -----------------------------------------------------------------

    @property
    def has_total_scattered_stellar_luminosity_cube_earth(self):
        return self.total_simulations.has_full_cube

    # -----------------------------------------------------------------

    @property
    def total_scattered_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.faceon_observed_cube_scattered

    # -----------------------------------------------------------------

    @property
    def has_total_scattered_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.has_full_cube_faceon

    # -----------------------------------------------------------------

    @property
    def total_scattered_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.edgeon_observed_cube_scattered

    # -----------------------------------------------------------------

    @property
    def has_total_scattered_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.has_full_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def total_absorbed_diffuse_stellar_sed_earth(self):
        return self.total_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def total_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.total_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_sed_earth(self):
        return self.total_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.total_simulations.has_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def total_absorbed_diffuse_stellar_sed_faceon(self):
        return self.total_simulations.faceon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def total_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_sed_faceon(self):
        return self.total_simulations.has_faceon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def total_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.total_simulations.edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def total_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.total_simulations.has_edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.has_edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def total_attenuated_stellar_luminosity_cube_earth(self):
        return self.total_simulations.observed_cube_attenuated

    # -----------------------------------------------------------------

    @property
    def has_total_attenuated_stellar_luminosity_cube_earth(self):
        return self.total_simulations.has_observed_cube_attenuated

    # -----------------------------------------------------------------

    @property
    def total_attenuated_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.faceon_observed_cube_attenuated

    # -----------------------------------------------------------------

    @property
    def has_total_attenuated_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.has_faceon_observed_cube_attenuated

    # -----------------------------------------------------------------

    @property
    def total_attenuated_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.edgeon_observed_cube_attenuated

    # -----------------------------------------------------------------

    @property
    def has_total_attenuated_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.has_edgeon_observed_cube_attenuated

    # -----------------------------------------------------------------

    @property
    def total_direct_stellar_luminosity_cube_earth(self):
        return self.total_simulations.observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_total_direct_stellar_luminosity_cube_earth(self):
        return self.total_simulations.has_full_cube

    # -----------------------------------------------------------------

    @property
    def total_direct_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.faceon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_total_direct_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.has_full_cube_faceon

    # -----------------------------------------------------------------

    @property
    def total_direct_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.edgeon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_total_direct_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.has_full_cube_edgeon

    # -----------------------------------------------------------------
    # BULGE CUBES
    # -----------------------------------------------------------------

    @property
    def bulge_bolometric_luminosity_cube_earth(self):
        return self.bulge_simulations.observed_cube

    # -----------------------------------------------------------------

    @property
    def bulge_bolometric_luminosity_cube_faceon(self):
        return self.bulge_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def bulge_bolometric_luminosity_cube_edgeon(self):
        return self.bulge_simulations.edgeon_observed_cube

    # -----------------------------------------------------------------

    @property
    def bulge_intrinsic_stellar_luminosity_cube_earth(self):
        return self.bulge_simulations.intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def bulge_intrinsic_stellar_luminosity_cube_faceon(self):
        return self.bulge_simulations.faceon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def bulge_intrinsic_stellar_luminosity_cube_edgeon(self):
        return self.bulge_simulations.edgeon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def bulge_observed_stellar_luminosity_cube_earth(self):
        return self.bulge_simulations.observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_bulge_observed_stellar_luminosity_cube_earth(self):
        return self.bulge_simulations.has_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def bulge_observed_stellar_luminosity_cube_faceon(self):
        return self.bulge_simulations.faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_bulge_observed_stellar_luminosity_cube_faceon(self):
        return self.bulge_simulations.has_faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def bulge_observed_stellar_luminosity_cube_edgeon(self):
        return self.bulge_simulations.edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_bulge_observed_stellar_luminosity_cube_edgeon(self):
        return self.bulge_simulations.has_edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def bulge_direct_stellar_luminosity_cube_earth(self):
        return self.bulge_simulations.observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_bulge_direct_stellar_luminosity_cube_earth(self):
        return self.bulge_simulations.has_full_cube

    # -----------------------------------------------------------------

    @property
    def bulge_direct_stellar_luminosity_cube_faceon(self):
        return self.bulge_simulations.faceon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_bulge_direct_stellar_luminosity_cube_faceon(self):
        return self.bulge_simulations.has_full_cube_faceon

    # -----------------------------------------------------------------

    @property
    def bulge_direct_stellar_luminosity_cube_edgeon(self):
        return self.bulge_simulations.edgeon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_bulge_direct_stellar_luminosity_cube_edgeon(self):
        return self.bulge_simulations.has_full_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def bulge_dust_luminosity_cube_earth(self):
        return self.bulge_simulations.observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_bulge_dust_luminosity_cube(self):
        return self.bulge_simulations.has_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def bulge_dust_luminosity_cube_faceon(self):
        return self.bulge_simulations.faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_bulge_dust_luminosity_cube_faceon(self):
        return self.bulge_simulations.has_faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def bulge_dust_luminosity_cube_edgeon(self):
        return self.bulge_simulations.edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_bulge_dust_luminosity_cube_edgeon(self):
        return self.bulge_simulations.has_edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def bulge_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.bulge_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def bulge_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.bulge_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def bulge_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.bulge_simulations.edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------
    # DISK CUBES
    # -----------------------------------------------------------------

    @property
    def disk_bolometric_luminosity_cube_earth(self):
        return self.disk_simulations.observed_cube

    # -----------------------------------------------------------------

    @property
    def disk_bolometric_luminosity_cube_faceon(self):
        return self.disk_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def disk_bolometric_luminosity_cube_edgeon(self):
        return self.disk_simulations.edgeon_observed_cube

    # -----------------------------------------------------------------

    @property
    def disk_intrinsic_stellar_luminosity_cube_earth(self):
        return self.disk_simulations.intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def disk_intrinsic_stellar_luminosity_cube_faceon(self):
        return self.disk_simulations.faceon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def disk_intrinsic_stellar_luminosity_cube_edgeon(self):
        return self.disk_simulations.edgeon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def disk_observed_stellar_luminosity_cube_earth(self):
        return self.disk_simulations.observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_disk_observed_stellar_luminosity_cube_earth(self):
        return self.disk_simulations.has_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def disk_observed_stellar_luminosity_cube_faceon(self):
        return self.disk_simulations.faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_disk_observed_stellar_luminosity_cube_faceon(self):
        return self.disk_simulations.has_faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def disk_observed_stellar_luminosity_cube_edgeon(self):
        return self.disk_simulations.edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_disk_observed_stellar_luminosity_cube_edgeon(self):
        return self.disk_simulations.has_edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def disk_direct_stellar_luminosity_cube_earth(self):
        return self.disk_simulations.observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_disk_direct_stellar_luminosity_cube_earth(self):
        return self.disk_simulations.has_full_cube

    # -----------------------------------------------------------------

    @property
    def disk_direct_stellar_luminosity_cube_faceon(self):
        return self.disk_simulations.faceon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_disk_direct_stellar_luminosity_cube_faceon(self):
        return self.disk_simulations.has_full_cube_faceon

    # -----------------------------------------------------------------

    @property
    def disk_direct_stellar_luminosity_cube_edgeon(self):
        return self.disk_simulations.edgeon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_disk_direct_stellar_luminosity_cube_edgeon(self):
        return self.disk_simulations.has_full_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def old_disk_direct_stellar_sed(self):
        return self.old_disk_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_old_disk_direct_stellar_sed(self):
        return self.has_old_disk_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def old_disk_direct_stellar_sed_earth(self):
        return self.disk_simulations.observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_old_disk_direct_stellar_sed_earth(self):
        return self.disk_simulations.has_direct_sed

    # -----------------------------------------------------------------

    @property
    def old_disk_direct_stellar_sed_faceon(self):
        return self.disk_simulations.faceon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_old_disk_direct_stellar_sed_faceon(self):
        return self.disk_simulations.has_direct_sed_faceon

    # -----------------------------------------------------------------

    @property
    def old_disk_direct_stellar_sed_edgeon(self):
        return self.disk_simulations.edgeon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_old_disk_direct_stellar_sed_edgeon(self):
        return self.disk_simulations.has_direct_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def disk_dust_luminosity_cube_earth(self):
        return self.disk_simulations.observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_disk_dust_luminosity_cube_earth(self):
        return self.disk_simulations.has_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def disk_dust_luminosity_cube_faceon(self):
        return self.disk_simulations.faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_disk_dust_luminosity_cube_faceon(self):
        return self.disk_simulations.has_faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def disk_dust_luminosity_cube_edgeon(self):
        return self.disk_simulations.edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_disk_dust_luminosity_cube_edgeon(self):
        return self.disk_simulations.has_edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def disk_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.disk_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def disk_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.disk_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def disk_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.disk_simulations.edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------
    # OLD CUBES
    # -----------------------------------------------------------------

    @property
    def old_bolometric_luminosity_cube_earth(self):
        return self.old_simulations.observed_cube

    # -----------------------------------------------------------------

    @property
    def has_old_bolometric_luminosity_cube_earth(self):
        return self.old_simulations.has_observed_cube

    # -----------------------------------------------------------------

    @property
    def old_bolometric_luminosity_cube_faceon(self):
        return self.old_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_old_bolometric_luminosity_cube_faceon(self):
        return self.old_simulations.has_faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def old_bolometric_luminosity_cube_edgeon(self):
        return self.old_simulations.edgeon_oberved_cube

    # -----------------------------------------------------------------

    @property
    def has_old_bolometric_luminosity_cube_edgeon(self):
        return self.old_simulations.has_edgeon_observed_cube

    # -----------------------------------------------------------------

    @property
    def old_intrinsic_stellar_luminosity_cube_earth(self):
        return self.old_simulations.intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def old_intrinsic_stellar_luminosity_cube_faceon(self):
        return self.old_simulations.faceon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def old_intrinsic_stellar_luminosity_cube_edgeon(self):
        return self.old_simulations.edgeon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def old_observed_stellar_luminosity_cube_earth(self):
        return self.old_simulations.observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def old_observed_stellar_luminosity_cube_faceon(self):
        return self.old_simulations.faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def old_observed_stellar_luminosity_cube_edgeon(self):
        return self.old_simulations.edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def old_direct_stellar_luminosity_cube_earth(self):
        return self.old_simulations.observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def old_direct_stellar_luminosity_cube_faceon(self):
        return self.old_simulations.faceon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def old_direct_stellar_luminosity_cube_edgeon(self):
        return self.old_simulations.edgeon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def old_dust_luminosity_cube_earth(self):
        return self.old_simulations.observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_old_dust_luminosity_cube_earth(self):
        return self.old_simulations.has_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def old_dust_luminosity_cube_faceon(self):
        return self.old_simulations.faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_old_dust_luminosity_cube_faceon(self):
        return self.old_simulations.has_faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def old_dust_luminosity_cube_edgeon(self):
        return self.old_simulations.edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_old_dust_luminosity_cube_edgeon(self):
        return self.old_simulations.has_edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def old_absorbed_diffuse_stellar_sed_earth(self):
        return self.old_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def old_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.old_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_old_absorbed_diffuse_stellar_sed_earth(self):
        return self.old_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_old_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.old_simulations.has_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def old_absorbed_diffuse_stellar_sed_faceon(self):
        return self.old_simulations.faceon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def old_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.old_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_old_absorbed_diffuse_stellar_sed_faceon(self):
        return self.old_simulations.has_faceon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_old_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.old_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def old_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.old_simulations.edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def old_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.old_simulations.edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_old_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.old_simulations.has_edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_old_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.old_simulations.has_edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------
    # YOUNG CUBES
    # -----------------------------------------------------------------

    @property
    def young_bolometric_luminosity_cube_earth(self):
        return self.young_simulations.observed_cube

    # -----------------------------------------------------------------

    @property
    def young_bolometric_luminosity_cube_faceon(self):
        return self.young_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def young_bolometric_luminosity_cube_egeon(self):
        return self.young_simulations.edgeon_observed_cube

    # -----------------------------------------------------------------

    @property
    def young_intrinsic_stellar_luminosity_cube_earth(self):
        return self.young_simulations.intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def young_intrinsic_stellar_luminosity_cube_faceon(self):
        return self.young_simulations.faceon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def young_intrinsic_stellar_luminosity_cube_edgeon(self):
        return self.young_simulations.edgeon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def young_observed_stellar_luminosity_cube_earth(self):
        return self.young_simulations.observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_young_observed_stellar_luminosity_cube_earth(self):
        return self.young_simulations.has_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def young_observed_stellar_luminosity_cube_faceon(self):
        return self.young_simulations.faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_young_observed_stellar_luminosity_cube_faceon(self):
        return self.young_simulations.has_faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def young_observed_stellar_luminosity_cube_edgeon(self):
        return self.young_simulations.edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_young_observed_stellar_luminosity_cube_edgeon(self):
        return self.young_simulations.has_edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def young_direct_stellar_luminosity_cube_earth(self):
        return self.young_simulations.observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_young_direct_stellar_luminosity_cube_earth(self):
        return self.young_simulations.has_full_cube

    # -----------------------------------------------------------------

    @property
    def young_direct_stellar_luminosity_cube_faceon(self):
        return self.young_simulations.faceon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_young_direct_stellar_luminosity_cube_faceon(self):
        return self.young_simulations.has_full_cube_faceon

    # -----------------------------------------------------------------

    @property
    def young_direct_stellar_luminosity_cube_edgeon(self):
        return self.young_simulations.edgeon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_young_direct_stellar_luminosity_cube_edgeon(self):
        return self.young_simulations.has_full_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def young_direct_stellar_sed(self):
        return self.young_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_young_direct_stellar_sed(self):
        return self.has_young_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def young_direct_stellar_sed_earth(self):
        return self.young_simulations.observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_young_direct_stellar_sed_earth(self):
        return self.young_simulations.has_direct_sed

    # -----------------------------------------------------------------

    @property
    def young_direct_stellar_sed_faceon(self):
        return self.young_simulations.faceon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_young_direct_stellar_sed_faceon(self):
        return self.young_simulations.has_direct_sed_faceon

    # -----------------------------------------------------------------

    @property
    def young_direct_stellar_sed_edgeon(self):
        return self.young_simulations.edgeon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_young_direct_stellar_sed_edgeon(self):
        return self.young_simulations.has_direct_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def young_dust_luminosity_cube_earth(self):
        return self.young_simulations.observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_young_dust_luminosity_cube_earth(self):
        return self.young_simulations.has_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def young_dust_luminosity_cube_faceon(self):
        return self.young_simulations.faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_young_dust_luminosity_cube_faceon(self):
        return self.young_simulations.has_faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def young_dust_luminosity_cube_edgeon(self):
        return self.young_simulations.edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_young_dust_luminosity_cube_edgeon(self):
        return self.young_simulations.has_edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def young_absorbed_diffuse_stellar_sed_earth(self):
        return self.young_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def young_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.young_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_young_absorbed_diffuse_stellar_sed_earth(self):
        return self.young_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_young_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.young_simulations.has_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def young_absorbed_diffuse_stellar_sed_faceon(self):
        return self.young_simulations.faceon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def young_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.young_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_young_absorbed_diffuse_stellar_sed_faceon(self):
        return self.young_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_young_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.young_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def young_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.young_simulations.edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def young_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.young_simulations.edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_young_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.young_simulations.has_edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_young_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.young_simulations.has_edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------
    # SFR CUBES
    # -----------------------------------------------------------------

    @property
    def sfr_bolometric_luminosity_cube_earth(self):
        return self.sfr_simulations.observed_cube

    # -----------------------------------------------------------------

    @property
    def sfr_bolometric_luminosity_cube_faceon(self):
        return self.sfr_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def sfr_bolometric_luminosity_cube_edgeon(self):
        return self.sfr_simulations.edgeon_observed_cube

    # -----------------------------------------------------------------

    @property
    def sfr_intrinsic_stellar_luminosity_cube_earth(self):
        return self.sfr_simulations.intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def sfr_intrinsic_stellar_luminosity_cube_faceon(self):
        return self.sfr_simulations.faceon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def sfr_intrinsic_stellar_luminosity_cube_edgeon(self):
        return self.sfr_simulations.edgeon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def sfr_observed_stellar_luminosity_cube_earth(self):
        return self.sfr_simulations.observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_sfr_observed_stellar_luminosity_cube_earth(self):
        return self.sfr_simulations.has_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def sfr_observed_stellar_luminosity_cube_faceon(self):
        return self.sfr_simulations.faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_sfr_observed_stellar_luminosity_cube_faceon(self):
        return self.sfr_simulations.has_faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def sfr_observed_stellar_luminosity_cube_edgeon(self):
        return self.sfr_simulations.edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_sfr_observed_stellar_luminosity_cube_edgeon(self):
        return self.sfr_simulations.has_edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def sfr_direct_stellar_luminosity_cube_earth(self):
        return self.sfr_simulations.observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_sfr_direct_stellar_luminosity_cube_earth(self):
        return self.sfr_simulations.has_full_cube

    # -----------------------------------------------------------------

    @property
    def sfr_direct_stellar_luminosity_cube_faceon(self):
        return self.sfr_simulations.faceon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_sfr_direct_stellar_luminosity_cube_faceon(self):
        return self.sfr_simulations.has_full_cube_faceon

    # -----------------------------------------------------------------

    @property
    def sfr_direct_stellar_luminosity_cube_edgeon(self):
        return self.sfr_simulations.edgeon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_sfr_direct_stellar_luminosity_cube_edgeon(self):
        return self.sfr_simulations.has_full_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_direct_stellar_sed(self):
        return self.sfr_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_direct_stellar_sed(self):
        return self.has_sfr_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def sfr_direct_stellar_sed_earth(self):
        return self.sfr_simulations.observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_sfr_direct_stellar_sed_earth(self):
        return self.sfr_simulations.has_direct_sed

    # -----------------------------------------------------------------

    @property
    def sfr_direct_stellar_sed_faceon(self):
        return self.sfr_simulations.faceon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_sfr_direct_stellar_sed_faceon(self):
        return self.sfr_simulations.has_direct_sed_faceon

    # -----------------------------------------------------------------

    @property
    def sfr_direct_stellar_sed_edgeon(self):
        return self.sfr_simulations.edgeon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_sfr_direct_stellar_sed_edgeon(self):
        return self.sfr_simulations.has_direct_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_dust_luminosity_cube_earth(self):
        return self.sfr_simulations.observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_luminosity_cube_earth(self):
        return self.sfr_simulations.has_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def sfr_dust_luminosity_cube_faceon(self):
        return self.sfr_simulations.faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_luminosity_cube_faceon(self):
        return self.sfr_simulations.has_faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def sfr_dust_luminosity_cube_edgeon(self):
        return self.sfr_simulations.edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_luminosity_cube_edgeon(self):
        return self.sfr_simulations.has_edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def sfr_absorbed_diffuse_stellar_sed_earth(self):
        return self.sfr_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def sfr_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.sfr_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_sfr_absorbed_diffuse_stellar_sed_earth(self):
        return self.sfr_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_sfr_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.sfr_simulations.has_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def sfr_absorbed_diffuse_stellar_sed_faceon(self):
        return self.sfr_simulations.faceon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def sfr_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.sfr_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_sfr_absorbed_diffuse_stellar_sed_faceon(self):
        return self.sfr_simulations.has_faceon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_sfr_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.sfr_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def sfr_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.sfr_simulations.edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def sfr_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.sfr_simulations.edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_sfr_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.sfr_simulations.has_edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_sfr_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.sfr_simulations.has_edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------
    # UNEVOLVED CUBES
    # -----------------------------------------------------------------

    @property
    def unevolved_bolometric_luminosity_cube_earth(self):
        return self.unevolved_simulations.observed_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_bolometric_luminosity_cube_faceon(self):
        return self.unevolved_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_bolometric_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.edgeon_observed_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_intrinsic_stellar_luminosity_cube_earth(self):
        return self.unevolved_simulations.intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_intrinsic_stellar_luminosity_cube_faceon(self):
        return self.unevolved_simulations.faceon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_intrinsic_stellar_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.edgeon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_observed_stellar_luminosity_cube_earth(self):
        return self.unevolved_simulations.observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_unevolved_observed_stellar_luminosity_cube_earth(self):
        return self.unevolved_simulations.has_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_observed_stellar_luminosity_cube_faceon(self):
        return self.unevolved_simulations.faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_unevolved_observed_stellar_luminosity_cube_faceon(self):
        return self.unevolved_simulations.has_faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_observed_stellar_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_unevolved_observed_stellar_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.has_edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_direct_stellar_luminosity_cube_earth(self):
        return self.unevolved_simulations.observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_unevolved_direct_stellar_luminosity_cube_earth(self):
        return self.unevolved_simulations.has_full_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_direct_stellar_luminosity_cube_faceon(self):
        return self.unevolved_simulations.faceon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_unevolved_direct_stellar_luminosity_cube_faceon(self):
        return self.unevolved_simulations.has_full_cube_faceon

    # -----------------------------------------------------------------

    @property
    def unevolved_direct_stellar_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.edgeon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_unevolved_direct_stellar_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.has_full_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_direct_stellar_sed(self):
        return self.unevolved_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_direct_stellar_sed(self):
        return self.has_unevolved_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def unevolved_direct_stellar_sed_earth(self):
        return self.unevolved_simulations.observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_unevolved_direct_stellar_sed_earth(self):
        return self.unevolved_simulations.has_direct_sed

    # -----------------------------------------------------------------

    @property
    def unevolved_direct_stellar_sed_faceon(self):
        return self.unevolved_simulations.faceon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_unevolved_direct_stellar_sed_faceon(self):
        return self.unevolved_simulations.has_direct_sed_faceon

    # -----------------------------------------------------------------

    @property
    def unevolved_direct_stellar_sed_edgeon(self):
        return self.unevolved_simulations.edgeon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_unevolved_direct_stellar_sed_edgeon(self):
        return self.unevolved_simulations.has_direct_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_dust_luminosity_cube_earth(self):
        return self.unevolved_simulations.observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_luminosity_cube_earth(self):
        return self.unevolved_simulations.has_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_dust_luminosity_cube_faceon(self):
        return self.unevolved_simulations.faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_luminosity_cube_faceon(self):
        return self.unevolved_simulations.has_faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_dust_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.has_edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def unevolved_absorbed_diffuse_stellar_sed_earth(self):
        return self.unevolved_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def unevolved_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.unevolved_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorbed_diffuse_stellar_sed_earth(self):
        return self.unevolved_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.unevolved_simulations.has_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def unevolved_absorbed_diffuse_stellar_sed_faceon(self):
        return self.unevolved_simulations.faceon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def unevolved_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.unevolved_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorbed_diffuse_stellar_sed_faceon(self):
        return self.unevolved_simulations.has_faceon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.unevolved_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def unevolved_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.unevolved_simulations.edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def unevolved_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.unevolved_simulations.has_edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.has_edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------
    # EXTRA CUBES
    # -----------------------------------------------------------------

    @property
    def extra_bolometric_luminosity_cube_earth(self):
        return self.extra_simulations.observed_cube

    # -----------------------------------------------------------------

    @property
    def extra_bolometric_luminosity_cube_faceon(self):
        return self.extra_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def extra_bolometric_luminosity_cube_egeon(self):
        return self.extra_simulations.edgeon_observed_cube

    # -----------------------------------------------------------------

    @property
    def extra_intrinsic_stellar_luminosity_cube_earth(self):
        return self.extra_simulations.intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def extra_intrinsic_stellar_luminosity_cube_faceon(self):
        return self.extra_simulations.faceon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def extra_intrinsic_stellar_luminosity_cube_edgeon(self):
        return self.extra_simulations.edgeon_intrinsic_stellar_cube

    # -----------------------------------------------------------------

    @property
    def extra_observed_stellar_luminosity_cube_earth(self):
        return self.extra_simulations.observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_extra_observed_stellar_luminosity_cube_earth(self):
        return self.extra_simulations.has_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def extra_observed_stellar_luminosity_cube_faceon(self):
        return self.extra_simulations.faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_extra_observed_stellar_luminosity_cube_faceon(self):
        return self.extra_simulations.has_faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def extra_observed_stellar_luminosity_cube_edgeon(self):
        return self.extra_simulations.edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_extra_observed_stellar_luminosity_cube_edgeon(self):
        return self.extra_simulations.has_edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def extra_direct_stellar_luminosity_cube_earth(self):
        return self.extra_simulations.observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_extra_direct_stellar_luminosity_cube_earth(self):
        return self.extra_simulations.has_full_cube

    # -----------------------------------------------------------------

    @property
    def extra_direct_stellar_luminosity_cube_faceon(self):
        return self.extra_simulations.faceon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_extra_direct_stellar_luminosity_cube_faceon(self):
        return self.extra_simulations.has_full_cube_faceon

    # -----------------------------------------------------------------

    @property
    def extra_direct_stellar_luminosity_cube_edgeon(self):
        return self.extra_simulations.edgeon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def has_extra_direct_stellar_luminosity_cube_edgeon(self):
        return self.extra_simulations.has_full_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def extra_direct_stellar_sed(self):
        return self.extra_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_extra_direct_stellar_sed(self):
        return self.has_extra_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def extra_direct_stellar_sed_earth(self):
        return self.extra_simulations.observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_extra_direct_stellar_sed_earth(self):
        return self.extra_simulations.has_direct_sed

    # -----------------------------------------------------------------

    @property
    def extra_direct_stellar_sed_faceon(self):
        return self.extra_simulations.faceon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_extra_direct_stellar_sed_faceon(self):
        return self.extra_simulations.has_direct_sed_faceon

    # -----------------------------------------------------------------

    @property
    def extra_direct_stellar_sed_edgeon(self):
        return self.extra_simulations.edgeon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_extra_direct_stellar_sed_edgeon(self):
        return self.extra_simulations.has_direct_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def extra_dust_luminosity_cube_earth(self):
        return self.extra_simulations.observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_extra_dust_luminosity_cube_earth(self):
        return self.extra_simulations.has_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def extra_dust_luminosity_cube_faceon(self):
        return self.extra_simulations.faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_extra_dust_luminosity_cube_faceon(self):
        return self.extra_simulations.has_faceon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def extra_dust_luminosity_cube_edgeon(self):
        return self.extra_simulations.edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_extra_dust_luminosity_cube_edgeon(self):
        return self.extra_simulations.has_edgeon_observed_dust_cube

    # -----------------------------------------------------------------

    @property
    def extra_absorbed_diffuse_stellar_sed_earth(self):
        return self.extra_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def extra_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.extra_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_extra_absorbed_diffuse_stellar_sed_earth(self):
        return self.extra_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_extra_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.extra_simulations.has_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def extra_absorbed_diffuse_stellar_sed_faceon(self):
        return self.extra_simulations.faceon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def extra_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.extra_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_extra_absorbed_diffuse_stellar_sed_faceon(self):
        return self.extra_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_extra_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.extra_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def extra_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.extra_simulations.edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def extra_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.extra_simulations.edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_extra_absorbed_diffuse_stellar_sed_edgeon(self):
        return self.extra_simulations.has_edgeon_observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_extra_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.extra_simulations.has_edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------
    # TOTAL MAPS
    # -----------------------------------------------------------------
    # 1. TOTAL BOLOMETRIC LUMINOSITY
    # -----------------------------------------------------------------

    @property
    def total_bolometric_luminosity_map(self):
        return self.total_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_bolometric_luminosity_map(self):
        return self.has_total_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_bolometric_luminosity_map_earth(self):
        return self.total_bolometric_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_bolometric_luminosity_map_earth(self):
        return self.has_total_bolometric_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_bolometric_luminosity_map_faceon(self):
        return self.total_bolometric_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_bolometric_luminosity_map_faceon(self):
        return self.has_total_bolometric_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_bolometric_luminosity_map_edgeon(self):
        return self.total_bolometric_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_bolometric_luminosity_map_edgeon(self):
        return self.has_total_bolometric_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # 2. INTRINSIC STELLAR LUMINOSITY
    # -----------------------------------------------------------------

    @property
    def total_intrinsic_stellar_luminosity_map(self):
        return self.total_intrinsic_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_intrinsic_stellar_luminosity_map(self):
        return self.has_total_intrinsic_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_intrinsic_stellar_luminosity_map_earth(self):
        return self.total_intrinsic_stellar_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_intrinsic_stellar_luminosity_map_earth(self):
        return self.has_total_intrinsic_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_intrinsic_stellar_luminosity_map_faceon(self):
        return self.total_intrinsic_stellar_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_intrinsic_stellar_luminosity_map_faceon(self):
        return self.has_total_intrinsic_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_intrinsic_stellar_luminosity_map_edgeon(self):
        return self.total_intrinsic_stellar_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_intrinsic_stellar_luminosity_map_edgeon(self):
        return self.has_total_intrinsic_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # 3. OBSERVED STELLAR LUMINOSITY
    # -----------------------------------------------------------------

    @property
    def total_observed_stellar_luminosity_map(self):
        return self.total_observed_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_observed_stellar_luminosity_map(self):
        return self.has_total_observed_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_observed_stellar_luminosity_map_earth(self):
        return self.total_observed_stellar_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_observed_stellar_luminosity_map_earth(self):
        return self.has_total_observed_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_observed_stellar_luminosity_map_faceon(self):
        return self.total_observed_stellar_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_observed_stellar_luminosity_map_faceon(self):
        return self.has_total_observed_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_observed_stellar_luminosity_map_edgeon(self):
        return self.total_observed_stellar_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_observed_stellar_luminosity_map_edgeon(self):
        return self.has_total_observed_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # 4. DIFFUSE DUST LUMINOSITY
    # -----------------------------------------------------------------

    @property
    def total_diffuse_dust_luminosity_map(self):
        return self.total_diffuse_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_diffuse_dust_luminosity_map(self):
        return self.has_total_diffuse_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_diffuse_dust_luminosity_map_earth(self):
        return self.total_diffuse_dust_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_diffuse_dust_luminosity_map_earth(self):
        return self.has_total_diffuse_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_diffuse_dust_luminosity_map_faceon(self):
        return self.total_diffuse_dust_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_diffuse_dust_luminosity_map_faceon(self):
        return self.has_total_diffuse_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_diffuse_dust_luminosity_map_edgeon(self):
        return self.total_diffuse_dust_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_diffuse_dust_luminosity_map_edgeon(self):
        return self.has_total_diffuse_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # 5. DUST LUMINOSITY
    # -----------------------------------------------------------------

    @property
    def total_dust_luminosity_map(self):
        return self.total_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_dust_luminosity_map(self):
        return self.has_total_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_luminosity_map_earth(self):
        return self.total_dust_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_dust_luminosity_map_earth(self):
        return self.has_total_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_luminosity_map_faceon(self):
        return self.total_dust_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_dust_luminosity_map_faceon(self):
        return self.has_total_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_luminosity_map_edgeon(self):
        return self.total_dust_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_dust_luminosity_map_edgeon(self):
        return self.has_total_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # 6. SCATTERED STELLAR LUMINOSITY
    # -----------------------------------------------------------------

    @property
    def total_scattered_stellar_luminosity_map(self):
        return self.total_scattered_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_scattered_stellar_luminosity_map(self):
        return self.has_total_scattered_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_scattered_stellar_luminosity_map_earth(self):
        return self.total_scattered_stellar_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_scattered_stellar_luminosity_map_earth(self):
        return self.has_total_scattered_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_scattered_stellar_luminosity_map_faceon(self):
        return self.total_scattered_stellar_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_scattered_stellar_luminosity_map_faceon(self):
        return self.has_total_scattered_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_scattered_stellar_luminosity_map_edgeon(self):
        return self.total_scattered_stellar_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_scattered_stellar_luminosity_map_edgeon(self):
        return self.has_total_scattered_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # 7. ABSORBED STELLAR LUMINOSITY BY DIFFUSE DUST
    # -----------------------------------------------------------------

    @property
    def total_absorbed_diffuse_stellar_luminosity_map(self):
        return self.total_absorbed_diffuse_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_luminosity_map(self):
        return self.has_total_absorbed_diffuse_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorbed_diffuse_stellar_luminosity_map_earth(self):
        return self.total_absorbed_diffuse_stellar_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_luminosity_map_earth(self):
        return self.has_total_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorbed_diffuse_stellar_luminosity_map_faceon(self):
        return self.total_absorbed_diffuse_stellar_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_luminosity_map_faceon(self):
        return self.has_total_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorbed_diffuse_stellar_luminosity_map_edgeon(self):
        return self.total_absorbed_diffuse_stellar_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_luminosity_map_edgeon(self):
        return self.has_total_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # X. ABSORBED STELLAR LUMINOSITY WITH INTERNAL -> THIS INFORMATION (MAPPINGS ABSORPTION) IS NOT AVAILABLE
    # -----------------------------------------------------------------
    # 8. FABS DIFFUSE
    # -----------------------------------------------------------------

    @property
    def total_fabs_diffuse_map(self):
        return self.total_fabs_diffuse_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_fabs_diffuse_map(self):
        return self.has_total_fabs_diffuse_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_fabs_diffuse_map_earth(self):
        return self.total_diffuse_dust_luminosity_map_earth / self.total_intrinsic_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_fabs_diffuse_map_earth(self):
        return self.has_total_diffuse_dust_luminosity_map_earth and self.has_total_intrinsic_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_fabs_diffuse_map_faceon(self):
        return self.total_diffuse_dust_luminosity_map_faceon / self.total_intrinsic_stellar_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_fabs_diffuse_map_faceon(self):
        return self.has_total_diffuse_dust_luminosity_map_faceon and self.has_total_intrinsic_stellar_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_fabs_diffuse_map_edgeon(self):
        return self.total_diffuse_dust_luminosity_map_edgeon / self.total_intrinsic_stellar_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_fabs_diffuse_map_edgeon(self):
        return self.has_total_diffuse_dust_luminosity_map_edgeon and self.has_total_intrinsic_stellar_luminosity_map_edgeon

    # -----------------------------------------------------------------
    # 9. FABS
    # -----------------------------------------------------------------

    @property
    def total_fabs_map(self):
        return self.total_fabs_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_fabs_map(self):
        return self.has_total_fabs_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_fabs_map_earth(self):
        return self.total_dust_luminosity_map_earth / self.total_intrinsic_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_fabs_map_earth(self):
        return self.has_total_dust_luminosity_map_earth and self.has_total_intrinsic_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_fabs_map_faceon(self):
        return self.total_dust_luminosity_map_faceon / self.total_intrinsic_stellar_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_fabs_map_faceon(self):
        return self.has_total_dust_luminosity_map_faceon and self.has_total_intrinsic_stellar_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_fabs_map_edgeon(self):
        return self.total_dust_luminosity_map_edgeon / self.total_intrinsic_stellar_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_fabs_map_edgeon(self):
        return self.has_total_dust_luminosity_map_edgeon and self.has_total_intrinsic_stellar_luminosity_map_edgeon

    # -----------------------------------------------------------------
    # 10. ATTENUATED
    # -----------------------------------------------------------------

    @property
    def total_attenuated_stellar_luminosity_map(self):
        return self.total_attenuated_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_attenuated_stellar_luminosity_map(self):
        return self.has_total_attenuated_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_attenuated_stellar_luminosity_map_earth(self):
        return self.total_attenuated_stellar_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_attenuated_stellar_luminosity_map_earth(self):
        return self.has_total_attenuated_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_attenuated_stellar_luminosity_map_faceon(self):
        return self.total_attenuated_stellar_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_attenuated_stellar_luminosity_map_faceon(self):
        return self.has_total_attenuated_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_attenuated_stellar_luminosity_map_edgeon(self):
        return self.total_attenuated_stellar_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_attenuated_stellar_luminosity_map_edgeon(self):
        return self.has_total_attenuated_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # 11. DIRECT
    # -----------------------------------------------------------------

    @property
    def total_direct_stellar_luminosity_map(self):
        return self.total_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_direct_stellar_luminosity_map(self):
        return self.has_total_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_direct_stellar_luminosity_map_earth(self):
        return self.total_direct_stellar_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_direct_stellar_luminosity_map_earth(self):
        return self.has_total_direct_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_direct_stellar_luminosity_map_faceon(self):
        return self.total_direct_stellar_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_direct_stellar_luminosity_map_faceon(self):
        return self.has_total_direct_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_direct_stellar_luminosity_map_edgeon(self):
        return self.total_direct_stellar_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_total_direct_stellar_luminosity_map_edgeon(self):
        return self.has_total_direct_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def total_direct_stellar_sed(self):
        return self.total_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_total_direct_stellar_sed(self):
        return self.has_total_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def total_direct_stellar_sed_earth(self):
        return self.total_simulations.observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_total_direct_stellar_sed_earth(self):
        return self.total_simulations.has_direct_sed

    # -----------------------------------------------------------------

    @property
    def total_direct_stellar_sed_faceon(self):
        return self.total_simulations.faceon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_total_direct_stellar_sed_faceon(self):
        return self.total_simulations.has_direct_sed_faceon

    # -----------------------------------------------------------------

    @property
    def total_direct_stellar_sed_edgeon(self):
        return self.total_simulations.edgeon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_total_direct_stellar_sed_edgeon(self):
        return self.total_simulations.has_direct_sed_edgeon

    # -----------------------------------------------------------------
    # 12. INTRINSIC FUV LUMINOSITY
    # -----------------------------------------------------------------

    @property
    def intrinsic_fuv_luminosity_map(self):
        return self.intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_map(self):
        return self.has_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_map_earth(self):
        return self.total_intrinsic_stellar_luminosity_cube_earth.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_map_earth(self):
        return self.has_total_intrinsic_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_map_faceon(self):
        return self.total_intrinsic_stellar_luminosity_cube_faceon.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_map_faceon(self):
        return self.has_total_intrinsic_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_map_edgeon(self):
        return self.total_intrinsic_stellar_luminosity_cube_edgeon.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_map_edgeon(self):
        return self.has_total_intrinsic_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # X. 24 MICRON LUMINOSITY
    # -----------------------------------------------------------------

    @property
    def total_observed_luminosity_cube_earth(self):
        return self.total_simulations.observed_cube

    # -----------------------------------------------------------------

    @property
    def has_total_observed_luminosity_cube_earth(self):
        return self.total_simulations.has_observed_cube

    # -----------------------------------------------------------------

    @property
    def total_observed_luminosity_cube_faceon(self):
        return self.total_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_total_observed_luminosity_cube_faceon(self):
        return self.total_simulations.has_faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def total_observed_luminosity_cube_edgeon(self):
        return self.total_simulations.edgeon_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_total_observed_luminosity_cube_edgeon(self):
        return self.total_simulations.has_edgeon_observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def total_24um_luminosity_map_earth(self):
        return self.total_observed_luminosity_cube_earth.get_frame_for_wavelength(self.mips24_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_total_24um_luminosity_map_earth(self):
        return self.has_total_observed_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_24um_luminosity_map_faceon(self):
        return self.total_observed_luminosity_cube_faceon.get_frame_for_wavelength(self.mips24_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_total_24um_luminosity_map_faceon(self):
        return self.has_total_observed_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_24um_luminosity_map_edgeon(self):
        return self.total_observed_luminosity_cube_edgeon.get_frame_for_wavelength(self.mips24_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_total_24um_luminosity_map_edgeon(self):
        return self.has_total_observed_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # 13. STAR FORMATION RATE
    #   SALIM
    #     EARTH
    # -----------------------------------------------------------------

    @property
    def total_star_formation_rate_map_salim(self):
        return self.total_star_formation_rate_map_earth_salim

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_salim(self):
        return self.has_total_star_formation_rate_map_earth_salim

    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_earth_salim(self):
        return salim_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_earth, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_earth_salim(self):
        return self.has_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_faceon_salim(self):
        return salim_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_faceon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_faceon_salim(self):
        return self.has_intrinsic_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_edgeon_salim(self):
        return salim_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_edgeon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_edgeon_salim(self):
        return self.has_intrinsic_fuv_luminosity_map_edgeon

    # -----------------------------------------------------------------
    #   KENNICUTT
    #     EARTH
    # -----------------------------------------------------------------

    @property
    def total_star_formation_rate_map_kennicutt(self):
        return self.total_star_formation_rate_map_earth_kennicutt

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_kennicutt(self):
        return self.has_total_star_formation_rate_map_earth_kennicutt

    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_earth_kennicutt(self):
        return kennicutt_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_earth, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_earth_kennicutt(self):
        return self.has_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_faceon_kennicutt(self):
        return kennicutt_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_faceon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_faceon_kennicutt(self):
        return self.has_intrinsic_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_edgeon_kennicutt(self):
        return kennicutt_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_edgeon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_edgeon_kennicutt(self):
        return self.has_intrinsic_fuv_luminosity_map_edgeon

    # -----------------------------------------------------------------
    #   KENNICUTT & EVANS
    #     EARTH
    # -----------------------------------------------------------------

    @property
    def total_star_formation_rate_map_ke(self):
        return self.total_star_formation_rate_map_earth_ke

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_ke(self):
        return self.has_total_star_formation_rate_map_earth_ke

    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_earth_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_earth, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_earth_ke(self):
        return self.has_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_faceon_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_faceon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_faceon_ke(self):
        return self.has_intrinsic_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_edgeon_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_edgeon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_edgeon_ke(self):
        return self.has_intrinsic_fuv_luminosity_map_edgeon

    # -----------------------------------------------------------------
    #   TIR
    #     EARTH
    # -----------------------------------------------------------------

    @property
    def total_star_formation_rate_map_tir(self):
        return self.total_star_formation_rate_map_earth_tir

    # -----------------------------------------------------------------

    @lazyproperty
    def has_total_star_formation_rate_map_tir(self):
        return self.has_total_star_formation_rate_map_earth_tir

    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_earth_tir(self):
        return kennicutt_tir_to_sfr(self.total_dust_luminosity_map_earth, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_earth_tir(self):
        return self.has_total_dust_luminosity_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_faceon_tir(self):
        return kennicutt_tir_to_sfr(self.total_dust_luminosity_map_faceon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_faceon_tir(self):
        return self.has_total_dust_luminosity_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_edgeon_tir(self):
        return kennicutt_tir_to_sfr(self.total_dust_luminosity_map_edgeon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_edgeon_tir(self):
        return self.has_total_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------
    #   24 micron
    #     EARTH
    # -----------------------------------------------------------------

    @property
    def total_star_formation_rate_map_24um(self):
        return self.total_star_formation_rate_map_earth_24um

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_24um(self):
        return self.has_total_star_formation_rate_map_earth_24um

    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_earth_24um(self):
        return calzetti_24um_to_sfr(self.total_24um_luminosity_map_earth, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_earth_24um(self):
        return self.has_total_24um_luminosity_map_earth

    # -----------------------------------------------------------------
    #    FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_faceon_24um(self):
        return calzetti_24um_to_sfr(self.total_24um_luminosity_map_faceon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_faceon_24um(self):
        return self.has_total_24um_luminosity_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_edgeon_24um(self):
        return calzetti_24um_to_sfr(self.total_24um_luminosity_map_edgeon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_edgeon_24um(self):
        return self.has_total_24um_luminosity_map_edgeon

    # -----------------------------------------------------------------
    # 14. I1 LUMINOSITY
    # -----------------------------------------------------------------

    @property
    def observed_i1_luminosity_map(self):
        return self.observed_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_map(self):
        return self.has_observed_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_map_earth(self):
        return self.total_bolometric_luminosity_cube_earth.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_map_earth(self):
        return self.has_total_bolometric_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_map_faceon(self):
        return self.total_bolometric_luminosity_cube_faceon.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_map_faceon(self):
        return self.has_total_bolometric_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_map_edgeon(self):
        return self.total_bolometric_luminosity_cube_edgeon.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_map_edgeon(self):
        return self.has_total_bolometric_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_i1_luminosity_map_earth(self):
        return self.old_bolometric_luminosity_cube_earth.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_i1_luminosity_map_earth(self):
        return self.has_old_bolometric_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_i1_luminosity_map_faceon(self):
        return self.old_bolometric_luminosity_cube_faceon.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_i1_luminosity_map_faceon(self):
        return self.has_old_bolometric_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_old_i1_luminosity_map_edgeon(self):
        return self.old_bolometric_luminosity_cube_edgeon.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_i1_luminosity_map_edgeon(self):
        return self.has_old_bolometric_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # 15. STELLAR MASS
    #   Total
    # -----------------------------------------------------------------

    @property
    def total_stellar_mass_map(self):
        return self.total_stellar_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_stellar_mass_map(self):
        return self.has_total_stellar_mass_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_stellar_mass_map_earth(self):
        return oliver_stellar_mass(self.observed_i1_luminosity_map_earth, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_stellar_mass_map_earth(self):
        return self.has_observed_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_stellar_mass_map_faceon(self):
        return oliver_stellar_mass(self.observed_i1_luminosity_map_faceon, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_stellar_mass_map_faceon(self):
        return self.has_observed_i1_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_stellar_mass_map_edgeon(self):
        return oliver_stellar_mass(self.observed_i1_luminosity_map_edgeon, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_stellar_mass_map_edgeon(self):
        return self.has_observed_i1_luminosity_map_edgeon

    # -----------------------------------------------------------------
    #   Old
    # -----------------------------------------------------------------

    @property
    def old_stellar_mass_map(self):
        return self.old_stellar_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_stellar_mass_map(self):
        return self.has_old_stellar_mass_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stellar_mass_map_earth(self):
        return oliver_stellar_mass(self.observed_old_i1_luminosity_map_earth, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_old_stellar_mass_map_earth(self):
        return self.has_observed_old_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stellar_mass_map_faceon(self):
        return oliver_stellar_mass(self.observed_old_i1_luminosity_map_faceon, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_old_stellar_mass_map_faceon(self):
        return self.has_observed_old_i1_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stellar_mass_map_edgeon(self):
        return oliver_stellar_mass(self.observed_old_i1_luminosity_map_edgeon, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_old_stellar_mass_map_edgeon(self):
        return self.has_observed_old_i1_luminosity_map_edgeon

    # -----------------------------------------------------------------
    # 16. SPECIFIC STAR FORMATION RATE
    #   SALIM
    #     EARTH
    # -----------------------------------------------------------------

    @property
    def total_ssfr_map_salim(self):
        return self.total_ssfr_map_earth_salim

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_salim(self):
        return self.has_total_ssfr_map_earth_salim

    # -----------------------------------------------------------------

    @lazyproperty
    def total_ssfr_map_earth_salim(self):
        return self.total_star_formation_rate_map_earth_salim / self.total_stellar_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_earth_salim(self):
        return self.has_total_star_formation_rate_map_earth_salim and self.has_total_stellar_mass_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_ssfr_map_faceon_salim(self):
        return self.total_star_formation_rate_map_faceon_salim / self.total_stellar_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_faceon_salim(self):
        return self.has_total_star_formation_rate_map_faceon_salim and self.has_total_stellar_mass_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_ssfr_map_edgeon_salim(self):
        return self.total_star_formation_rate_map_edgeon_salim / self.total_stellar_mass_map_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_edgeon(self):
        return self.has_total_star_formation_rate_map_edgeon_salim and self.has_total_stellar_mass_map_edgeon

    # -----------------------------------------------------------------
    #   KENNICUTT
    #     EARTH
    # -----------------------------------------------------------------

    @property
    def total_ssfr_map_kennicutt(self):
        return self.total_ssfr_map_earth_kennicutt

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_kennicutt(self):
        return self.has_total_ssfr_map_earth_kennicutt

    # -----------------------------------------------------------------

    @lazyproperty
    def total_ssfr_map_earth_kennicutt(self):
        return self.total_star_formation_rate_map_earth_kennicutt / self.total_stellar_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_earth_kennicutt(self):
        return self.has_total_star_formation_rate_map_earth_kennicutt and self.has_total_stellar_mass_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_ssfr_map_faceon_kennicutt(self):
        return self.total_star_formation_rate_map_faceon_kennicutt / self.total_stellar_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_faceon_kennicutt(self):
        return self.has_total_star_formation_rate_map_faceon_kennicutt and self.has_total_stellar_mass_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_ssfr_map_edgeon_kennicutt(self):
        return self.total_star_formation_rate_map_edgeon_kennicutt / self.total_stellar_mass_map_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_edgeon_kennicutt(self):
        return self.has_total_star_formation_rate_map_edgeon_kennicutt and self.has_total_stellar_mass_map_edgeon

    # -----------------------------------------------------------------
    #   KENNICUTT & EVANS
    #     EARTH
    # -----------------------------------------------------------------

    @property
    def total_ssfr_map_ke(self):
        return self.total_ssfr_map_earth_ke

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_ke(self):
        return self.has_total_ssfr_map_earth_ke

    # -----------------------------------------------------------------

    @lazyproperty
    def total_ssfr_map_earth_ke(self):
        return self.total_star_formation_rate_map_earth_ke / self.total_stellar_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_earth_ke(self):
        return self.has_total_star_formation_rate_map_earth_ke and self.has_total_stellar_mass_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_ssfr_map_faceon_ke(self):
        return self.total_star_formation_rate_map_faceon_ke / self.total_stellar_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_faceon_ke(self):
        return self.has_total_star_formation_rate_map_faceon_ke and self.has_total_stellar_mass_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_ssfr_map_edgeon_ke(self):
        return self.total_star_formation_rate_map_edgeon_ke / self.total_stellar_mass_map_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_map_edgeon_ke(self):
        return self.has_total_star_formation_rate_map_edgeon_ke and self.has_total_stellar_mass_map_edgeon

    # -----------------------------------------------------------------
    # BULGE MAPS
    # -----------------------------------------------------------------

    @property
    def old_bulge_intrinsic_i1_luminosity_map(self):
        return self.old_bulge_intrinsic_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_intrinsic_i1_luminosity_map_earth(self):
        return self.old_bulge_map_earth.normalized(to=self.intrinsic_i1_luminosity_old_bulge)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_intrinsic_i1_luminosity_map_faceon(self):
        return self.old_bulge_map_faceon.normalized(to=self.intrinsic_i1_luminosity_old_bulge)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_intrinsic_i1_luminosity_map_edgeon(self):
        return self.old_bulge_map_edgeon.normalized(to=self.intrinsic_i1_luminosity_old_bulge)

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_intrinsic_i1_luminosity_map(self):
        return self.has_old_bulge_intrinsic_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_intrinsic_i1_luminosity_map_earth(self):
        return self.has_intrinsic_i1_luminosity_old_bulge and self.has_old_bulge_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_intrinsic_i1_luminosity_map_faceon(self):
        return self.has_intrinsic_i1_luminosity_old_bulge and self.has_old_bulge_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_intrinsic_i1_luminosity_map_edgeon(self):
        return self.has_intrinsic_i1_luminosity_old_bulge and self.has_old_bulge_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_bulge_bolometric_luminosity_map(self):
        return self.old_bulge_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_bolometric_luminosity_map_earth(self):
        return self.old_bulge_map_earth.normalized(to=self.intrinsic_i1_luminosity_old_bulge)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_bolometric_luminosity_map_faceon(self):
        return self.old_bulge_map_faceon.normalized(to=self.intrinsic_i1_luminosity_old_bulge)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_bolometric_luminosity_map_edgeon(self):
        return self.old_bulge_map_edgeon.normalized(to=self.intrinsic_i1_luminosity_old_bulge)

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_bolometric_luminosity_map(self):
        return self.has_old_bulge_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_bolometric_luminosity_map_earth(self):
        return self.has_intrinsic_bolometric_luminosity_old_bulge and self.has_old_bulge_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_bolometric_luminosity_map_faceon(self):
        return self.has_intrinsic_bolometric_luminosity_old_bulge and self.has_old_bulge_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_bolometric_luminosity_map_edgeon(self):
        return self.has_intrinsic_bolometric_luminosity_old_bulge and self.has_old_bulge_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_bulge_direct_stellar_luminosity_map(self):
        return self.old_bulge_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_direct_stellar_luminosity_map(self):
        return self.has_old_bulge_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_direct_stellar_luminosity_map_earth(self):
        return self.bulge_direct_stellar_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_direct_stellar_luminosity_map_earth(self):
        return self.has_bulge_direct_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_direct_stellar_luminosity_map_faceon(self):
        return self.bulge_direct_stellar_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_direct_stellar_luminosity_map_faceon(self):
        return self.has_bulge_direct_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_direct_stellar_luminosity_map_edgeon(self):
        return self.bulge_direct_stellar_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_direct_stellar_luminosity_map_edgeon(self):
        return self.has_bulge_direct_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def old_bulge_direct_stellar_sed(self):
        return self.old_bulge_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_direct_stellar_sed(self):
        return self.has_old_bulge_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def old_bulge_direct_stellar_sed_earth(self):
        return self.bulge_simulations.observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_direct_stellar_sed_earth(self):
        return self.bulge_simulations.has_direct_sed

    # -----------------------------------------------------------------

    @property
    def old_bulge_direct_stellar_sed_faceon(self):
        return self.bulge_simulations.faceon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_direct_stellar_sed_faceon(self):
        return self.bulge_simulations.has_direct_sed_faceon

    # -----------------------------------------------------------------

    @property
    def old_bulge_direct_stellar_sed_edgeon(self):
        return self.bulge_simulations.edgeon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_direct_stellar_sed_edgeon(self):
        return self.bulge_simulations.has_direct_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def old_bulge_i1_luminosity_map(self):
        return self.old_bulge_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_i1_luminosity_map(self):
        return self.has_old_bulge_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_i1_luminosity_map_earth(self):
        return self.bulge_observed_stellar_luminosity_cube_earth.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_i1_luminosity_map_earth(self):
        return self.has_bulge_observed_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_i1_luminosity_map_faceon(self):
        return self.bulge_observed_stellar_luminosity_cube_faceon.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_i1_luminosity_map_faceon(self):
        return self.has_bulge_observed_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_i1_luminosity_map_edgeon(self):
        return self.bulge_observed_stellar_luminosity_cube_edgeon.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_i1_luminosity_map_edgeon(self):
        return self.has_bulge_observed_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def old_bulge_dust_luminosity_map(self):
        return self.old_bulge_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_dust_luminosity_map(self):
        return self.has_old_bulge_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_dust_luminosity_map_earth(self):
        return self.bulge_dust_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_dust_luminosity_map_earth(self):
        return self.has_bulge_dust_luminosity_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_dust_luminosity_map_faceon(self):
        return self.bulge_dust_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_dust_luminosity_map_faceon(self):
        return self.has_bulge_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_dust_luminosity_map_edgeon(self):
        return self.bulge_dust_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_dust_luminosity_map_edgeon(self):
        return self.has_bulge_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # DISK MAPS
    # -----------------------------------------------------------------

    @property
    def old_disk_map_path(self):
        return self.definition.old_stars_map_path

    # -----------------------------------------------------------------

    @property
    def has_old_disk_map(self):
        return fs.is_file(self.old_disk_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_map(self):
        return Frame.from_file(self.old_disk_map_path)

    # -----------------------------------------------------------------

    @property
    def old_disk_map_psf_filter(self):
        return self.old_disk_map.psf_filter

    # -----------------------------------------------------------------

    @property
    def old_disk_map_fwhm(self):
        return self.old_disk_map.fwhm

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_map_wcs(self):
        return CoordinateSystem.from_file(self.old_disk_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_map_projection(self):

        """
        This function ...
        :return:
        """

        azimuth = 0.0
        if not self.has_center: raise ValueError("Galaxy center coordinate is not defined")
        return GalaxyProjection.from_wcs(self.old_disk_map_wcs, self.center, self.distance, self.inclination, azimuth, self.position_angle)

    # -----------------------------------------------------------------

    @property
    def old_disk_map_shape(self):
        return self.old_disk_map_wcs.shape

    # -----------------------------------------------------------------

    @property
    def old_disk_map_pixelscale(self):
        return self.old_disk_map_wcs.pixelscale

    # -----------------------------------------------------------------

    @property
    def old_disk_intrinsic_i1_luminosity_map(self):
        return self.old_disk_intrinsic_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_intrinsic_i1_luminosity_map_earth(self):
        return self.old_disk_map_earth.normalized(to=self.intrinsic_i1_luminosity_old_disk)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_intrinsic_i1_luminosity_map_faceon(self):
        return self.old_disk_map_faceon.normalized(to=self.intrinsic_i1_luminosity_old_disk)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_intrinsic_i1_luminosity_map_edgeon(self):
        return self.old_disk_map_edgeon.normalized(to=self.intrinsic_i1_luminosity_old_disk)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_intrinsic_i1_luminosity_map(self):
        return self.has_old_disk_intrinsic_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_disk_intrinsic_i1_luminosity_map_earth(self):
        return self.has_intrinsic_i1_luminosity_old_disk and self.has_old_disk_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_disk_intrinsic_i1_luminosity_map_faceon(self):
        return self.has_intrinsic_i1_luminosity_old_disk and self.has_old_disk_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_disk_intrinsic_i1_luminosity_map_edgeon(self):
        return self.has_intrinsic_i1_luminosity_old_disk and self.has_old_disk_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_disk_bolometric_luminosity_map(self):
        return self.old_disk_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_disk_bolometric_luminosity_map(self):
        return self.has_old_disk_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_bolometric_luminosity_map_earth(self):
        return self.old_disk_map_earth.normalized(to=self.intrinsic_bolometric_luminosity_old_disk)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_bolometric_luminosity_map_earth(self):
        return self.has_intrinsic_bolometric_luminosity_old_disk and self.has_old_disk_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_bolometric_luminosity_map_faceon(self):
        return self.old_disk_map_faceon.normalized(to=self.intrinsic_bolometric_luminosity_old_disk)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_bolometric_luminosity_map_faceon(self):
        return self.has_intrinsic_bolometric_luminosity_old_disk and self.has_old_disk_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_bolometric_luminosity_map_edgeon(self):
        return self.old_disk_map_edgeon.normalized(to=self.intrinsic_bolometric_luminosity_old_disk)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_bolometric_luminosity_map_edgeon(self):
        return self.has_intrinsic_bolometric_luminosity_old_disk and self.has_old_disk_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_disk_direct_stellar_luminosity_map(self):
        return self.old_disk_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_disk_direct_stellar_luminosity_map(self):
        return self.has_old_disk_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_direct_stellar_luminosity_map_earth(self):
        return self.disk_direct_stellar_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_disk_direct_stellar_luminosity_map_earth(self):
        return self.has_disk_direct_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_direct_stellar_luminosity_map_faceon(self):
        return self.disk_direct_stellar_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_disk_direct_stellar_luminosity_map_faceon(self):
        return self.has_disk_direct_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_direct_stellar_luminosity_map_edgeon(self):
        return self.disk_direct_stellar_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_disk_direct_stellar_luminosity_map_edgeon(self):
        return self.has_disk_direct_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def old_disk_i1_luminosity_map(self):
        return self.old_disk_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_disk_i1_luminosity_map(self):
        return self.has_old_disk_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_i1_luminosity_map_earth(self):
        return self.disk_observed_stellar_luminosity_cube_earth.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_i1_luminosity_map_earth(self):
        return self.has_disk_observed_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_i1_luminosity_map_faceon(self):
        return self.disk_observed_stellar_luminosity_cube_faceon.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_i1_luminosity_map_faceon(self):
        return self.has_disk_observed_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_i1_luminosity_map_edgeon(self):
        return self.disk_observed_stellar_luminosity_cube_edgeon.get_frame_for_wavelength(self.i1_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_i1_luminosity_map_edgeon(self):
        return self.has_disk_observed_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def old_disk_dust_luminosity_map(self):
        return self.old_disk_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_disk_dust_luminosity_map(self):
        return self.has_old_disk_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_dust_luminosity_map_earth(self):
        return self.disk_dust_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_disk_dust_luminosity_map_earth(self):
        return self.has_disk_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_dust_luminosity_map_faceon(self):
        return self.disk_dust_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_disk_dust_luminosity_map_faceon(self):
        return self.has_disk_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_dust_luminosity_map_edgeon(self):
        return self.disk_dust_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_disk_dust_luminosity_map_edgeon(self):
        return self.has_disk_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # OLD MAPS
    # -----------------------------------------------------------------

    @property
    def old_bolometric_luminosity_map(self):
        return self.old_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bolometric_luminosity_map(self):
        return self.has_old_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bolometric_luminosity_map_earth(self):
        return self.old_bulge_bolometric_luminosity_map_earth + self.old_disk_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bolometric_luminosity_map_earth(self):
        return self.has_old_bulge_bolometric_luminosity_map_earth and self.has_old_disk_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bolometric_luminosity_map_faceon(self):
        return self.old_bulge_bolometric_luminosity_map_faceon + self.old_disk_bolometric_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_bolometric_luminosity_map_faceon(self):
        return self.has_old_bulge_bolometric_luminosity_map_faceon and self.has_old_disk_bolometric_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bolometric_luminosity_map_edgeon(self):
        return self.old_bulge_bolometric_luminosity_map_edgeon + self.old_disk_bolometric_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_bolometric_luminosity_map_edgeon(self):
        return self.has_old_bulge_bolometric_luminosity_map_edgeon and self.has_old_disk_bolometric_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_direct_stellar_luminosity_map(self):
        return self.old_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_direct_stellar_luminosity_map(self):
        return self.has_old_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_direct_stellar_luminosity_map_earth(self):
        return self.old_bulge_direct_stellar_luminosity_map_earth + self.old_disk_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_direct_stellar_luminosity_map_earth(self):
        return self.has_old_bulge_direct_stellar_luminosity_map_earth and self.has_old_disk_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_direct_stellar_luminosity_map_faceon(self):
        return self.old_bulge_direct_stellar_luminosity_map_faceon + self.old_disk_direct_stellar_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_direct_stellar_luminosity_map_faceon(self):
        return self.has_old_bulge_direct_stellar_luminosity_map_faceon and self.has_old_disk_direct_stellar_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_direct_stellar_luminosity_map_edgeon(self):
        return self.old_bulge_direct_stellar_luminosity_map_edgeon + self.old_disk_direct_stellar_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_direct_stellar_luminosity_map_edgeon(self):
        return self.has_old_bulge_direct_stellar_luminosity_map_edgeon and self.has_old_disk_direct_stellar_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_direct_stellar_sed(self):
        return self.old_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_old_direct_stellar_sed(self):
        return self.has_old_direct_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def old_direct_stellar_sed_earth(self):
        return self.old_simulations.observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_old_direct_stellar_sed_earth(self):
        return self.old_simulations.has_direct_sed

    # -----------------------------------------------------------------

    @property
    def old_direct_stellar_sed_faceon(self):
        return self.old_simulations.faceon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_old_direct_stellar_sed_faceon(self):
        return self.old_simulations.has_direct_sed_faceon

    # -----------------------------------------------------------------

    @property
    def old_direct_stellar_sed_edgeon(self):
        return self.old_simulations.edgeon_observed_sed_direct

    # -----------------------------------------------------------------

    @property
    def has_old_direct_stellar_sed_edgeon(self):
        return self.old_simulations.has_direct_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def old_i1_luminosity_map(self):
        return self.old_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_i1_luminosity_map(self):
        return self.has_old_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_i1_luminosity_map_earth(self):
        return self.old_bulge_i1_luminosity_map_earth + self.old_disk_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_i1_luminosity_map_earth(self):
        return self.has_old_bulge_i1_luminosity_map_earth and self.has_old_disk_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_i1_luminosity_map_faceon(self):
        return self.old_bulge_i1_luminosity_map_faceon + self.old_disk_i1_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_i1_luminosity_map_faceon(self):
        return self.has_old_bulge_i1_luminosity_map_faceon and self.has_old_disk_i1_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_i1_luminosity_map_edgeon(self):
        return self.old_bulge_i1_luminosity_map_edgeon + self.old_disk_i1_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_i1_luminosity_map_edgeon(self):
        return self.has_old_bulge_i1_luminosity_map_edgeon and self.has_old_disk_i1_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_intrinsic_i1_luminosity_map(self):
        return self.old_intrinsic_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_intrinsic_i1_luminosity_map(self):
        return self.has_old_intrinsic_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_intrinsic_i1_luminosity_map_earth(self):
        return self.old_bulge_intrinsic_i1_luminosity_map_earth + self.old_disk_intrinsic_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_intrinsic_i1_luminosity_map_earth(self):
        return self.has_old_bulge_intrinsic_i1_luminosity_map_earth and self.has_old_disk_intrinsic_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_intrinsic_i1_luminosity_map_faceon(self):
        return self.old_bulge_intrinsic_i1_luminosity_map_faceon + self.old_disk_intrinsic_i1_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_intrinsic_i1_luminosity_map_faceon(self):
        return self.has_old_bulge_intrinsic_i1_luminosity_map_faceon and self.has_old_disk_intrinsic_i1_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_intrinsic_i1_luminosity_map_edgeon(self):
        return self.old_bulge_intrinsic_i1_luminosity_map_edgeon + self.old_disk_intrinsic_i1_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_intrinsic_i1_luminosity_map_edgeon(self):
        return self.has_old_bulge_intrinsic_i1_luminosity_map_edgeon and self.has_old_disk_intrinsic_i1_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_dust_luminosity_map(self):
        return self.old_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_dust_luminosity_map(self):
        return self.has_old_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_dust_luminosity_map_earth(self):
        return self.old_dust_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_dust_luminosity_map_earth(self):
        return self.has_old_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_dust_luminosity_map_faceon(self):
        return self.old_dust_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_dust_luminosity_map_faceon(self):
        return self.has_old_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_dust_luminosity_map_edgeon(self):
        return self.old_dust_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_old_dust_luminosity_map_edgeon(self):
        return self.has_old_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # YOUNG MAPS
    # -----------------------------------------------------------------

    @property
    def young_map_path(self):
        return self.definition.young_stars_map_path

    # -----------------------------------------------------------------

    @property
    def has_young_map(self):
        return fs.is_file(self.young_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map(self):
        return Frame.from_file(self.young_map_path)

    # -----------------------------------------------------------------

    @property
    def young_map_psf_filter(self):
        return self.young_map.psf_filter

    # -----------------------------------------------------------------

    @property
    def young_map_fwhm(self):
        return self.young_map.fwhm

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map_wcs(self):
        return CoordinateSystem.from_file(self.young_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map_projection(self):

        """
        This function ...
        :return:
        """

        azimuth = 0.0
        if not self.has_center: raise ValueError("Galaxy center coordinate is not defined")
        return GalaxyProjection.from_wcs(self.young_map_wcs, self.center, self.distance, self.inclination, azimuth, self.position_angle)

    # -----------------------------------------------------------------

    @property
    def young_map_shape(self):
        return self.young_map_wcs.shape

    # -----------------------------------------------------------------

    @property
    def young_map_pixelscale(self):
        return self.young_map_wcs.pixelscale

    # -----------------------------------------------------------------
    # INTRINSIC FUV
    # -----------------------------------------------------------------

    @property
    def young_intrinsic_fuv_luminosity_map(self):
        return self.young_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_fuv_luminosity_map_earth(self):
        return self.young_map_earth.normalized(to=self.intrinsic_fuv_luminosity_young)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_fuv_luminosity_map_faceon(self):
        return self.young_map_faceon.normalized(to=self.intrinsic_fuv_luminosity_young)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_fuv_luminosity_map_edgeon(self):
        return self.young_map_edgeon.normalized(to=self.intrinsic_fuv_luminosity_young)

    # -----------------------------------------------------------------

    @property
    def has_young_intrinsic_fuv_luminosity_map(self):
        return self.has_young_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_young_intrinsic_fuv_luminosity_map_earth(self):
        return self.has_intrinsic_fuv_luminosity_young and self.has_young_map_earth

    # -----------------------------------------------------------------

    @property
    def has_young_intrinsic_fuv_luminosity_map_faceon(self):
        return self.has_intrinsic_fuv_luminosity_young and self.has_young_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_intrinsic_fuv_luminosity_map_edgeon(self):
        return self.has_intrinsic_fuv_luminosity_young and self.has_young_map_edgeon

    # -----------------------------------------------------------------
    # BOLOMETRIC LUMINOSITY
    # -----------------------------------------------------------------

    @property
    def young_bolometric_luminosity_map(self):
        return self.young_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def young_bolometric_luminosity_map_earth(self):
        return self.young_map_earth.normalized(to=self.intrinsic_bolometric_luminosity_young)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_bolometric_luminosity_map_faceon(self):
        return self.young_map_faceon.normalized(to=self.intrinsic_bolometric_luminosity_young)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_bolometric_luminosity_map_edgeon(self):
        return self.young_map_edgeon.normalized(to=self.intrinsic_bolometric_luminosity_young)

    # -----------------------------------------------------------------

    @property
    def has_young_bolometric_luminosity_map(self):
        return self.has_young_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_young_bolometric_luminosity_map_earth(self):
        return self.has_intrinsic_bolometric_luminosity_young and self.has_young_map_earth

    # -----------------------------------------------------------------

    @property
    def has_young_bolometric_luminosity_map_faceon(self):
        return self.has_intrinsic_bolometric_luminosity_young and self.has_young_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_bolometric_luminosity_map_edgeon(self):
        return self.has_intrinsic_bolometric_luminosity_young and self.has_young_map_edgeon

    # -----------------------------------------------------------------
    # DIRECT
    # -----------------------------------------------------------------

    @property
    def young_direct_stellar_luminosity_map(self):
        return self.young_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_young_direct_stellar_luminosity_map(self):
        return self.has_young_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def young_direct_stellar_luminosity_map_earth(self):
        return self.young_direct_stellar_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_young_direct_stellar_luminosity_map_earth(self):
        return self.has_young_direct_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def young_direct_stellar_luminosity_map_faceon(self):
        return self.young_direct_stellar_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_young_direct_stellar_luminosity_map_faceon(self):
        return self.has_young_direct_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def young_direct_stellar_luminosity_map_edgeon(self):
        return self.young_direct_stellar_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_young_direct_stellar_luminosity_map_edgeon(self):
        return self.has_young_direct_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # FUV LUMINOSITY
    # -----------------------------------------------------------------

    @property
    def young_fuv_luminosity_map(self):
        return self.young_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_young_fuv_luminosity_map(self):
        return self.has_young_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def young_fuv_luminosity_map_earth(self):
        return self.young_observed_stellar_luminosity_cube_earth.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_young_fuv_luminosity_map_earth(self):
        return self.has_young_observed_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def young_fuv_luminosity_map_faceon(self):
        return self.young_observed_stellar_luminosity_cube_faceon.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_young_fuv_luminosity_map_faceon(self):
        return self.has_young_observed_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def young_fuv_luminosity_map_edgeon(self):
        return self.young_observed_stellar_luminosity_cube_edgeon.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_young_fuv_luminosity_map_edgeon(self):
        return self.has_young_observed_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def young_dust_luminosity_map(self):
        return self.young_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_young_dust_luminosity_map(self):
        return self.has_young_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def young_dust_luminosity_map_earth(self):
        return self.young_dust_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_young_dust_luminosity_map_earth(self):
        return self.has_young_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def young_dust_luminosity_map_faceon(self):
        return self.young_dust_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_young_dust_luminosity_map_faceon(self):
        return self.has_young_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def young_dust_luminosity_map_edgeon(self):
        return self.young_dust_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_young_dust_luminosity_map_edgeon(self):
        return self.has_young_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # STAR FORMATION RATE
    #   SALIM
    #     EARTH
    # -----------------------------------------------------------------

    @property
    def young_star_formation_rate_map_salim(self):
        return self.young_star_formation_rate_map_earth_salim

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_salim(self):
        return self.has_young_star_formation_rate_map_earth_salim

    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_map_earth_salim(self):
        return salim_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_earth, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_earth_salim(self):
        return self.has_young_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_map_faceon_salim(self):
        return salim_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_faceon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_faceon_salim(self):
        return self.has_young_intrinsic_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_map_edgeon_salim(self):
        return salim_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_edgeon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_edgeon_salim(self):
        return self.has_young_intrinsic_fuv_luminosity_map_edgeon

    # -----------------------------------------------------------------
    #   KENNICUTT & EVANS
    #     EARTH
    # -----------------------------------------------------------------

    @property
    def young_star_formation_rate_map_ke(self):
        return self.young_star_formation_rate_map_earth_ke

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_ke(self):
        return self.has_young_star_formation_rate_map_earth_ke

    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_map_earth_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_earth, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_earth_ke(self):
        return self.has_young_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_map_faceon_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_faceon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_faceon_ke(self):
        return fs.is_file(self.young_star_formation_rate_map_faceon_ke)

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_map_edgeon_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_edgeon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_edgeon_ke(self):
        return self.has_young_intrinsic_fuv_luminosity_map_edgeon

    # -----------------------------------------------------------------
    # SFR MAPS
    # -----------------------------------------------------------------

    @property
    def sfr_map_path(self):
        return self.definition.ionizing_stars_map_path

    # -----------------------------------------------------------------

    @property
    def has_sfr_map(self):
        return fs.is_file(self.sfr_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_map(self):
        return Frame.from_file(self.sfr_map_path)

    # -----------------------------------------------------------------

    @property
    def sfr_map_psf_filter(self):
        return self.sfr_map.psf_filter

    # -----------------------------------------------------------------

    @property
    def sfr_map_fwhm(self):
        return self.sfr_map.fwhm

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_map_wcs(self):
        return CoordinateSystem.from_file(self.sfr_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_map_projection(self):

        """
        This function ...
        :return:
        """

        azimuth = 0.0
        if not self.has_center: raise ValueError("Galaxy center coordinate is not defined")
        return GalaxyProjection.from_wcs(self.sfr_map_wcs, self.center, self.distance, self.inclination, azimuth, self.position_angle)

    # -----------------------------------------------------------------

    @property
    def sfr_map_shape(self):
        return self.sfr_map_wcs.shape

    # -----------------------------------------------------------------

    @property
    def sfr_map_pixelscale(self):
        return self.sfr_map_wcs.pixelscale

    # -----------------------------------------------------------------

    @property
    def sfr_intrinsic_fuv_luminosity_map(self):
        return self.sfr_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_fuv_luminosity_map_earth(self):
        return self.sfr_map_earth.normalized(to=self.intrinsic_fuv_luminosity_sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_fuv_luminosity_map_faceon(self):
        return self.sfr_map_faceon.normalized(to=self.intrinsic_fuv_luminosity_sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsinc_fuv_luminosity_map_edgeon(self):
        return self.sfr_map_edgeon.normalized(to=self.intrinsic_fuv_luminosity_sfr)

    # -----------------------------------------------------------------

    @property
    def has_sfr_intrinsic_fuv_luminosity_map(self):
        return self.has_sfr_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_intrinsic_fuv_luminosity_map_earth(self):
        return self.has_intrinsic_fuv_luminosity_sfr and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_intrinsic_fuv_luminosity_map_faceon(self):
        return self.has_intrinsic_fuv_luminosity_sfr and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_sfr_intrinsic_fuv_luminosity_map_edgeon(self):
        return self.has_intrinsic_fuv_luminosity_sfr and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_bolometric_luminosity_map(self):
        return self.sfr_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_bolometric_luminosity_map_earth(self):
        return self.sfr_map_earth.normalized(to=self.intrinsic_bolometric_luminosity_sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_bolometric_luminosity_map_faceon(self):
        return self.sfr_map_faceon.normalized(to=self.intrinsic_bolometric_luminosity_sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_bolometric_luminosity_map_edgeon(self):
        return self.sfr_map_edgeon.normalized(to=self.intrinsic_bolometric_luminosity_sfr)

    # -----------------------------------------------------------------

    @property
    def has_sfr_bolometric_luminosity_map(self):
        return self.has_sfr_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_bolometric_luminosity_map_earth(self):
        return self.has_intrinsic_bolometric_luminosity_sfr and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_bolometric_luminosity_map_faceon(self):
        return self.has_intrinsic_bolometric_luminosity_sfr and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_sfr_bolometric_luminosity_map_edgeon(self):
        return self.has_intrinsic_bolometric_luminosity_sfr and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_direct_stellar_luminosity_map(self):
        return self.sfr_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_direct_stellar_luminosity_map(self):
        return self.has_sfr_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_direct_stellar_luminosity_map_earth(self):
        return self.sfr_direct_stellar_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_sfr_direct_stellar_luminosity_map_earth(self):
        return self.has_sfr_direct_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_direct_stellar_luminosity_map_faceon(self):
        return self.sfr_direct_stellar_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_sfr_direct_stellar_luminosity_map_faceon(self):
        return self.has_sfr_direct_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_direct_stellar_luminosity_map_edgeon(self):
        return self.sfr_direct_stellar_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_sfr_direct_stellar_luminosity_map_edgeon(self):
        return self.has_sfr_direct_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_fuv_luminosity_map(self):
        return self.sfr_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_fuv_luminosity_map(self):
        return self.has_sfr_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_fuv_luminosity_map_earth(self):
        return self.sfr_observed_stellar_luminosity_cube_earth.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_sfr_fuv_luminosity_map_earth(self):
        return self.has_sfr_observed_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_fuv_luminosity_map_faceon(self):
        return self.sfr_observed_stellar_luminosity_cube_faceon.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_sfr_fuv_luminosity_map_faceon(self):
        return self.has_sfr_observed_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_fuv_luminosity_map_edgeon(self):
        return self.sfr_observed_stellar_luminosity_cube_edgeon.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_sfr_fuv_luminosity_map_edgeon(self):
        return self.has_sfr_observed_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def star_formation_rate_map(self):
        return self.star_formation_rate_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def star_formation_rate_map_earth(self):
        return self.sfr_map_earth.normalized(to=self.sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def star_formation_rate_map_faceon(self):
        return self.sfr_map_faceon.normalized(to=self.sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def star_formation_rate_map_edgeon(self):
        return self.sfr_map_edgeon.normalized(to=self.sfr)

    # -----------------------------------------------------------------

    @property
    def has_star_formation_rate_map(self):
        return self.has_star_formation_rate_map_earth

    # -----------------------------------------------------------------

    @property
    def has_star_formation_rate_map_earth(self):
        return self.has_sfr and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_star_formation_rate_map_faceon(self):
        return self.has_sfr and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_star_formation_rate_map_edgeon(self):
        return self.has_sfr and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_dust_mass_map(self):
        return self.sfr_dust_mass_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_mass_map_earth(self):
        return self.sfr_map_earth.normalized(to=self.sfr_dust_mass)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_mass_map_faceon(self):
        return self.sfr_map_faceon.normalized(to=self.sfr_dust_mass)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_mass_map_edgeon(self):
        return self.sfr_map_edgeon.normalized(to=self.sfr_dust_mass)

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass_map(self):
        return self.has_sfr_dust_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass_map_earth(self):
        return self.has_sfr_dust_mass and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass_map_faceon(self):
        return self.has_sfr_dust_mass and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass_map_edgeon(self):
        return self.has_sfr_dust_mass and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_stellar_luminosity_map(self):
        return self.sfr_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_stellar_luminosity_map_earth(self):
        return self.sfr_map_earth.normalized(to=self.intrinsic_stellar_luminosity_sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_stellar_luminosity_map_faceon(self):
        return self.sfr_map_faceon.normalized(to=self.intrinsic_stellar_luminosity_sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_stellar_luminosity_map_edgeon(self):
        return self.sfr_map_edgeon.normalized(to=self.intrinsic_stellar_luminosity_sfr)

    # -----------------------------------------------------------------

    @property
    def has_sfr_stellar_luminosity_map(self):
        return self.has_sfr_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_stellar_luminosity_map_earth(self):
        return self.has_intrinsic_stellar_luminosity_sfr and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_stellar_luminosity_map_faceon(self):
        return self.has_intrinsic_stellar_luminosity_sfr and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_sfr_stellar_luminosity_map_edgeon(self):
        return self.has_intrinsic_stellar_luminosity_sfr and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_intrinsic_dust_luminosity_map(self):
        return self.sfr_intrinsic_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_dust_luminosity_map_earth(self):
        return self.sfr_map_earth.normalized(to=self.intrinsic_dust_luminosity_sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_dust_luminosity_map_faceon(self):
        return self.sfr_map_faceon.normalized(to=self.intrinsic_dust_luminosity_sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_dust_luminosity_map_edgeon(self):
        return self.sfr_map_edgeon.normalized(to=self.intrinsic_dust_luminosity_sfr)

    # -----------------------------------------------------------------

    @property
    def has_sfr_intrinsic_dust_luminosity_map(self):
        return self.has_sfr_intrinsic_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_intrinsic_dust_luminosity_map_earth(self):
        return self.has_intrinsic_dust_luminosity_sfr and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_intrinsic_dust_luminosity_map_faceon(self):
        return self.has_intrinsic_dust_luminosity_sfr and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_sfr_intrinsic_dust_luminosity_map_edgeon(self):
        return self.has_intrinsic_dust_luminosity_sfr and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_dust_luminosity_map(self):
        return self.sfr_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_luminosity_map(self):
        return self.has_sfr_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_luminosity_map_earth(self):
        return self.sfr_dust_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_luminosity_map_earth(self):
        return self.has_sfr_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_luminosity_map_faceon(self):
        return self.sfr_dust_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_luminosity_map_faceon(self):
        return self.has_sfr_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_luminosity_map_edgeon(self):
        return self.sfr_dust_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_luminosity_map_edgeon(self):
        return self.has_sfr_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # UNEVOLVED MAPS
    # -----------------------------------------------------------------

    @property
    def unevolved_intrinsic_fuv_luminosity_map(self):
        return self.unevolved_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_intrinsic_fuv_luminosity_map(self):
        return self.has_unevolved_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_intrinsic_fuv_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        # Get contributions
        young = self.young_intrinsic_fuv_luminosity_map_earth
        sfr = self.sfr_intrinsic_fuv_luminosity_map_earth

        # Uniformize
        #young, sfr = convolve_and_rebin(young, sfr)
        young, sfr = convolve_rebin_and_convert(young, sfr)

        # Sum the contributions
        return young + sfr

    # -----------------------------------------------------------------

    @property
    def has_unevolved_intrinsic_fuv_luminosity_map_earth(self):
        return self.has_young_intrinsic_fuv_luminosity_map_earth and self.has_sfr_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_intrinsic_fuv_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get contributions
        young = self.young_intrinsic_fuv_luminosity_map_faceon
        sfr = self.sfr_intrinsic_fuv_luminosity_map_faceon

        # Combine (no WCS): regrid and recenter??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_intrinsic_fuv_luminosity_map_faceon(self):
        return self.has_young_intrinsic_fuv_luminosity_map_faceon and self.has_sfr_intrinsic_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_intrinsic_fuv_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get contributions
        young = self.young_intrinsic_fuv_luminosity_map_edgeon
        sfr = self.young_intrinsic_fuv_luminosity_map_edgeon

        # Combine (no WCS): regrid and recenter??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_intrinsic_fuv_luminosity_map_edgeon(self):
        return self.has_unevolved_intrinsic_fuv_luminosity_map_edgeon and self.has_sfr_intrinsic_fuv_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_star_formation_rate_map(self):
        return self.unevolved_star_formation_rate_map_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_star_formation_rate_map(self):
        return self.has_unevolved_star_formation_rate_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_star_formation_rate_map_earth(self):
        return salim_fuv_to_sfr(self.unevolved_intrinsic_fuv_luminosity_map_earth, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_star_formation_rate_map_earth(self):
        return self.has_unevolved_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_star_formation_rate_map_faceon(self):
        return salim_fuv_to_sfr(self.unevolved_intrinsic_fuv_luminosity_map_faceon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_star_formation_rate_map_faceon(self):
        return self.has_unevolved_intrinsic_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_star_formation_rate_map_edgeon(self):
        return salim_fuv_to_sfr(self.unevolved_intrinsic_fuv_luminosity_map_edgeon, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_star_formation_rate_map_edgeon(self):
        return self.has_unevolved_intrinsic_fuv_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_fuv_luminosity_map(self):
        return self.has_unevolved_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def unevolved_bolometric_luminosity_map(self):
        return self.unevolved_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_bolometric_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        # Get contributions
        young = self.young_bolometric_luminosity_map_earth
        sfr = self.sfr_bolometric_luminosity_map_earth

        # Unformize
        #young, sfr = convolve_and_rebin(young, sfr)
        young, sfr = convolve_rebin_and_convert(young, sfr)

        # Sum the contributions
        return young + sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_bolometric_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get contributions
        young = self.young_bolometric_luminosity_map_faceon
        sfr = self.sfr_bolometric_luminosity_map_faceon

        # Combine (no WCS): regrid and recenter??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_bolometric_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get contributions
        young = self.young_bolometric_luminosity_map_edgeon
        sfr = self.sfr_bolometric_luminosity_map_edgeon

        # Combine (no WCS): regrid and recenter??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_bolometric_luminosity_map(self):
        return self.has_unevolved_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_bolometric_luminosity_map_earth(self):
        return self.has_young_bolometric_luminosity_map_earth and self.has_sfr_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_bolometric_luminosity_map_faceon(self):
        return self.has_young_bolometric_luminosity_map_faceon and self.has_sfr_bolometric_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_bolometric_luminosity_map_edgeon(self):
        return self.has_young_bolometric_luminosity_map_edgeon and self.has_sfr_bolometric_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_direct_stellar_luminosity_map(self):
        return self.unevolved_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_direct_stellar_luminosity_map(self):
        return self.has_unevolved_direct_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_direct_stellar_luminosity_map_earth(self):
        return self.unevolved_direct_stellar_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_unevolved_direct_stellar_luminosity_map_earth(self):
        return self.has_unevolved_direct_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_direct_stellar_luminosity_map_faceon(self):
        return self.unevolved_direct_stellar_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_unevolved_direct_stellar_luminosity_map_faceon(self):
        return self.has_unevolved_direct_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_direct_stellar_luminosity_map_edgeon(self):
        return self.unevolved_direct_stellar_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_unevolved_direct_stellar_luminosity_map_edgeon(self):
        return self.has_unevolved_direct_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_fuv_luminosity_map(self):
        return self.unevolved_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_fuv_luminosity_map_earth(self):
        return self.unevolved_observed_stellar_luminosity_cube_earth.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_fuv_luminosity_map_earth(self):
        return self.has_unevolved_observed_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_fuv_luminosity_map_faceon(self):
        return self.unevolved_observed_stellar_luminosity_cube_faceon.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_fuv_luminosity_map_faceon(self):
        return self.has_unevolved_observed_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_fuv_luminosity_map_edgeon(self):
        return self.unevolved_observed_stellar_luminosity_cube_edgeon.get_frame_for_wavelength(self.fuv_wavelength, copy=True)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_fuv_luminosity_map_edgeon(self):
        return self.has_unevolved_observed_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_dust_luminosity_map(self):
        return self.unevolved_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_luminosity_map(self):
        return self.has_unevolved_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_luminosity_map_earth(self):
        return self.unevolved_dust_luminosity_cube_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_luminosity_map_earth(self):
        return self.has_unevolved_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_luminosity_map_faceon(self):
        return self.unevolved_dust_luminosity_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_luminosity_map_faceon(self):
        return self.has_unevolved_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_luminosity_map_edgeon(self):
        return self.unevolved_dust_luminosity_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_luminosity_map_edgeon(self):
        return self.has_unevolved_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------
    # DUST MAPS
    # -----------------------------------------------------------------

    @property
    def dust_map_path(self):
        return self.definition.dust_map_path

    # -----------------------------------------------------------------

    @property
    def has_dust_map(self):
        return fs.is_file(self.dust_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map(self):
        return Frame.from_file(self.dust_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map_wcs(self):
        return CoordinateSystem.from_file(self.dust_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map_projection(self):
        azimuth = 0.0
        if not self.has_center: raise ValueError("Galaxy center coordinate is not defined")
        return GalaxyProjection.from_wcs(self.dust_map_wcs, self.center, self.distance, self.inclination, azimuth, self.position_angle)

    # -----------------------------------------------------------------

    @property
    def dust_map_shape(self):
        return self.dust_map_wcs.shape

    # -----------------------------------------------------------------

    @property
    def dust_map_pixelscale(self):
        return self.dust_map_wcs.pixelscale

    # -----------------------------------------------------------------

    @property
    def diffuse_dust_mass_map(self):
        return self.diffuse_dust_mass_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_mass_map_earth(self):
        return self.dust_map_earth.normalized(to=self.diffuse_dust_mass)

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_mass_map_faceon(self):
        return self.dust_map_faceon.normalized(to=self.diffuse_dust_mass)

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_mass_map_edgeon(self):
        return self.dust_map_edgeon.normalized(to=self.diffuse_dust_mass)

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_mass(self):
        return True

    # -----------------------------------------------------------------

    @property
    def has_dust_mass(self):
        return self.has_diffuse_dust_mass and self.has_sfr_dust_mass

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_mass_map(self):
        return self.has_diffuse_dust_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_mass_map_earth(self):
        return self.has_dust_map_earth and self.has_diffuse_dust_mass

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_mass_map_faceon(self):
        return self.has_dust_map_edgeon and self.has_diffuse_dust_mass

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_mass_map_edgeon(self):
        return self.has_dust_map_faceon and self.has_diffuse_dust_mass

    # -----------------------------------------------------------------

    @property
    def dust_mass_map(self):
        return self.dust_mass_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_mass_map_earth(self):

        """
        This function ...
        :return:
        """

        # Get the contributions
        diffuse = self.diffuse_dust_mass_map_earth
        sfr = self.sfr_dust_mass_map_earth

        # Uniformize
        #diffuse, sfr = convolve_and_rebin(diffuse, sfr)
        diffuse, sfr = convolve_rebin_and_convert(diffuse, sfr)

        # Sum the contributions
        return diffuse + sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the contributions
        diffuse = self.diffuse_dust_mass_map_faceon
        sfr = self.sfr_dust_mass_map_faceon

        # Combine (no WCS): regrid and center??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the contributions
        diffuse = self.diffuse_dust_mass_map_edgeon
        sfr = self.sfr_dust_mass_map_edgeon

        # Combine (no WCS): regrid and center??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def has_dust_mass_map(self):
        return self.has_dust_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_dust_mass_map_earth(self):
        return self.has_diffuse_dust_mass_map_earth and self.has_sfr_dust_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_dust_mass_map_faceon(self):
        return self.has_diffuse_dust_mass_map_faceon and self.has_sfr_dust_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_mass_map_edgeon(self):
        return self.has_diffuse_dust_mass_map_edgeon and self.has_sfr_dust_mass_map_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_input_filepaths(self):
        paths = OrderedDict()
        paths[wavelengths_filename] = self.wavelength_grid_path
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_input_filepaths(self):
        paths = OrderedDict()
        paths[wavelengths_filename] = self.wavelength_grid_path
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_input_filepaths(self):
        paths = OrderedDict()
        paths[wavelengths_filename] = self.wavelength_grid_path
        paths[map_filename] = self.old_disk_map_path
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def young_input_filepaths(self):
        paths = OrderedDict()
        paths[wavelengths_filename] = self.wavelength_grid_path
        paths[map_filename] = self.young_map_path
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_input_filepaths(self):
        paths = OrderedDict()
        paths[wavelengths_filename] = self.wavelength_grid_path
        paths[map_filename] = self.sfr_map_path
        return paths

    # -----------------------------------------------------------------

    @property
    def old_bulge_earth_projection_path(self):
        return fs.create_directory_in(self.old_bulge_projections_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_faceon_projection_path(self):
        return fs.create_directory_in(self.old_bulge_projections_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_edgeon_projection_path(self):
        return fs.create_directory_in(self.old_bulge_projections_path, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def old_bulge_earth_projection_output_filepath(self):
        return fs.join(self.old_bulge_earth_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_earth_projection_output_file(self):
        return fs.is_file(self.observed_total_simulation_output_filepath)

    # -----------------------------------------------------------------

    @property
    def old_bulge_faceon_projection_output_filepath(self):
        return fs.join(self.old_bulge_faceon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_faceon_projection_output_file(self):
        return fs.is_file(self.old_bulge_faceon_projection_output_filepath)

    # -----------------------------------------------------------------

    @property
    def old_bulge_edgeon_projection_output_filepath(self):
        return fs.join(self.old_bulge_edgeon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_edgeon_projection_output_file(self):
        return fs.join(self.old_bulge_edgeon_projection_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_projections(self):

        """
        This function ...
        :return:
        """

        return ComponentProjections(bulge_simulation_name, self.old_bulge_model, path=self.old_bulge_projections_path,
                                    description=bulge_component_description,
                                    projection=self.old_disk_projections.projection_earth,
                                    projection_faceon=self.old_disk_projections.projection_faceon,
                                    projection_edgeon=self.old_disk_projections.projection_edgeon, center=self.center,
                                    earth_wcs=self.old_disk_map_wcs, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def old_bulge_map(self):
        return self.old_bulge_map_earth

    # -----------------------------------------------------------------

    @property
    def old_bulge_map_earth(self):
        return self.old_bulge_projections.earth

    # -----------------------------------------------------------------

    @property
    def old_bulge_map_earth_path(self):
        return self.old_bulge_projections.earth_map_path

    # -----------------------------------------------------------------

    @property
    def old_bulge_map_faceon(self):
        return self.old_bulge_projections.faceon

    # -----------------------------------------------------------------

    @property
    def old_bulge_map_faceon_path(self):
        return self.old_bulge_projections.faceon_map_path

    # -----------------------------------------------------------------

    @property
    def old_bulge_map_edgeon(self):
        return self.old_bulge_projections.edgeon

    # -----------------------------------------------------------------

    @property
    def old_bulge_map_edgeon_path(self):
        return self.old_bulge_projections.edgeon_map_path

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_map(self):
        return self.has_old_bulge_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_map_earth(self):
        return True

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_map_faceon(self):
        return True

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_map_edgeon(self):
        return True

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_faceon_projection_path(self):
        return fs.create_directory_in(self.old_disk_projections_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_edgeon_projection_path(self):
        return fs.create_directory_in(self.old_disk_projections_path, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def old_disk_faceon_projection_output_filepath(self):
        return fs.join(self.old_disk_faceon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_faceon_projection_output_file(self):
        return fs.is_file(self.old_disk_faceon_projection_output_filepath)

    # -----------------------------------------------------------------

    @property
    def old_disk_edgeon_projection_output_filepath(self):
        return fs.join(self.old_disk_edgeon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_old_disk_edgeon_projection_output_file(self):
        return fs.is_file(self.old_disk_edgeon_projection_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_projections(self):

        """
        This function ...
        :return:
        """

        return ComponentProjections(disk_simulation_name, self.old_disk_deprojection, path=self.old_disk_projections_path,
                                    earth=False, description=disk_component_description,
                                    input_filepaths=[self.old_disk_map_path], center=self.center,
                                    earth_wcs=self.old_disk_map_wcs, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def old_disk_map_earth(self):
        return self.old_disk_map # just the input map

    # -----------------------------------------------------------------

    @property
    def old_disk_map_earth_path(self):
        return self.old_disk_map_path

    # -----------------------------------------------------------------

    @property
    def old_disk_map_faceon(self):
        return self.old_disk_projections.faceon

    # -----------------------------------------------------------------

    @property
    def old_disk_map_faceon_path(self):
        return self.old_disk_projections.faceon_map_path

    # -----------------------------------------------------------------

    @property
    def old_disk_map_edgeon(self):
        return self.old_disk_projections.edgeon

    # -----------------------------------------------------------------

    @property
    def old_disk_map_edgeon_path(self):
        return self.old_disk_projections.edgeon_map_path

    # -----------------------------------------------------------------

    @property
    def has_old_disk_map_earth(self):
        return self.has_old_disk_map

    # -----------------------------------------------------------------

    @property
    def has_old_disk_map_faceon(self):
        return self.has_old_disk_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------

    @property
    def has_old_disk_map_edgeon(self):
        return self.has_old_disk_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------

    @lazyproperty
    def young_faceon_projection_path(self):
        return fs.create_directory_in(self.young_projections_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_edgeon_projection_path(self):
        return fs.create_directory_in(self.young_projections_path, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def young_faceon_projection_output_filepath(self):
        return fs.join(self.young_faceon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_young_faceon_projection_output_file(self):
        return fs.is_file(self.young_faceon_projection_output_filepath)

    # -----------------------------------------------------------------

    @property
    def young_edgeon_projection_output_filepath(self):
        return fs.join(self.young_edgeon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_young_edgeon_projection_output_file(self):
        return fs.is_file(self.young_edgeon_projection_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_projections(self):

        """
        This function ...
        :return:
        """

        return ComponentProjections(young_simulation_name, self.young_deprojection, path=self.young_projections_path,
                                    earth=False, description=young_component_description,
                                    input_filepaths=[self.young_map_path], center=self.center,
                                    earth_wcs=self.young_map_wcs, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def young_map_earth(self):
        return self.young_map # just the input map

    # -----------------------------------------------------------------

    @property
    def young_map_earth_path(self):
        return self.young_map_path

    # -----------------------------------------------------------------

    @property
    def young_map_faceon(self):
        return self.young_projections.faceon

    # -----------------------------------------------------------------

    @property
    def young_map_faceon_path(self):
        return self.young_projections.faceon_map_path

    # -----------------------------------------------------------------

    @property
    def young_map_edgeon(self):
        return self.young_projections.edgeon

    # -----------------------------------------------------------------

    @property
    def young_map_edgeon_path(self):
        return self.young_projections.edgeon_map_path

    # -----------------------------------------------------------------

    @property
    def has_young_map_earth(self):
        return self.has_young_map

    # -----------------------------------------------------------------

    @property
    def has_young_map_faceon(self):
        return self.has_young_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------

    @property
    def has_young_map_edgeon(self):
        return self.has_young_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_faceon_projection_path(self):
        return fs.create_directory_in(self.sfr_projections_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_edgeon_projection_path(self):
        return fs.create_directory_in(self.sfr_projections_path, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def sfr_faceon_projection_output_filepath(self):
        return fs.join(self.sfr_faceon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_sfr_faceon_projection_output_file(self):
        return fs.is_file(self.sfr_faceon_projection_output_filepath)

    # -----------------------------------------------------------------

    @property
    def sfr_edgeon_projection_output_filepath(self):
        return fs.join(self.sfr_edgeon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_sfr_edgeon_projection_output_file(self):
        return fs.is_file(self.sfr_edgeon_projection_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_projections(self):

        """
        This function ...
        :return:
        """

        return ComponentProjections(sfr_simulation_name, self.sfr_deprojection, path=self.sfr_projections_path,
                                    earth=False, description=ionizing_component_description, input_filepaths=[self.sfr_map_path],
                                    center=self.center, earth_wcs=self.sfr_map_wcs, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def sfr_map_earth(self):
        return self.sfr_map # just the input map

    # -----------------------------------------------------------------

    @property
    def sfr_map_earth_path(self):
        return self.sfr_map_path

    # -----------------------------------------------------------------

    @property
    def sfr_map_faceon(self):
        return self.sfr_projections.faceon

    # -----------------------------------------------------------------

    @property
    def sfr_map_faceon_path(self):
        return self.sfr_projections.faceon_map_path

    # -----------------------------------------------------------------

    @property
    def sfr_map_edgeon(self):
        return self.sfr_projections.edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_map_edgeon_path(self):
        return self.sfr_projections.edgeon_map_path

    # -----------------------------------------------------------------

    @property
    def has_sfr_map_earth(self):
        return self.has_sfr_map

    # -----------------------------------------------------------------

    @property
    def has_sfr_map_faceon(self):
        return self.has_sfr_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------

    @property
    def has_sfr_map_edgeon(self):
        return self.has_sfr_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def extra_earth_projection_path(self):
        return fs.create_directory_in(self.extra_projections_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_faceon_projection_path(self):
        return fs.create_directory_in(self.extra_projections_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_edgeon_projection_path(self):
        return fs.create_directory_in(self.extra_projections_path, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def extra_earth_projection_output_filepath(self):
        return fs.join(self.extra_earth_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_extra_earth_projection_output_file(self):
        return fs.is_file(self.observed_total_simulation_output_filepath)

    # -----------------------------------------------------------------

    @property
    def extra_faceon_projection_output_filepath(self):
        return fs.join(self.extra_faceon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_extra_faceon_projection_output_file(self):
        return fs.is_file(self.extra_faceon_projection_output_filepath)

    # -----------------------------------------------------------------

    @property
    def extra_edgeon_projection_output_filepath(self):
        return fs.join(self.extra_edgeon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_extra_edgeon_projection_output_file(self):
        return fs.join(self.extra_edgeon_projection_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_projections(self):

        """
        This function ...
        :return:
        """

        return ComponentProjections(bulge_simulation_name, self.extra_model, path=self.extra_projections_path,
                                    description=bulge_component_description,
                                    projection=self.old_disk_projections.projection_earth,
                                    projection_faceon=self.old_disk_projections.projection_faceon,
                                    projection_edgeon=self.old_disk_projections.projection_edgeon, center=self.center,
                                    earth_wcs=self.old_disk_map_wcs, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def extra_map(self):
        return self.extra_map_earth

    # -----------------------------------------------------------------

    @property
    def extra_map_earth(self):
        return self.extra_projections.earth

    # -----------------------------------------------------------------

    @property
    def extra_map_earth_path(self):
        return self.extra_projections.earth_map_path

    # -----------------------------------------------------------------

    @property
    def extra_map_faceon(self):
        return self.extra_projections.faceon

    # -----------------------------------------------------------------

    @property
    def extra_map_faceon_path(self):
        return self.extra_projections.faceon_map_path

    # -----------------------------------------------------------------

    @property
    def extra_map_edgeon(self):
        return self.extra_projections.edgeon

    # -----------------------------------------------------------------

    @property
    def extra_map_edgeon_path(self):
        return self.extra_projections.edgeon_map_path

    # -----------------------------------------------------------------

    @property
    def has_extra_map(self):
        return self.has_extra_map_earth

    # -----------------------------------------------------------------

    @property
    def has_extra_map_earth(self):
        return True

    # -----------------------------------------------------------------

    @property
    def has_extra_map_faceon(self):
        return True

    # -----------------------------------------------------------------

    @property
    def has_extra_map_edgeon(self):
        return True

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def dust_faceon_projection_path(self):
        return fs.create_directory_in(self.dust_projections_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_edgeon_projection_path(self):
        return fs.create_directory_in(self.dust_projections_path, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def dust_faceon_projection_output_filepath(self):
        return fs.join(self.dust_faceon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_dust_faceon_projection_output_file(self):
        return fs.is_file(self.dust_faceon_projection_output_filepath)

    # -----------------------------------------------------------------

    @property
    def dust_edgeon_projection_output_filepath(self):
        return fs.join(self.dust_edgeon_projection_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_dust_edgeon_projection_output_file(self):
        return fs.is_file(self.dust_edgeon_projection_output_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_projections(self):

        """
        This function ...
        :return:
        """

        return ComponentProjections(dust_simulation_name, self.dust_deprojection, path=self.dust_projections_path,
                                    earth=False, description=dust_component_description, input_filepaths=[self.dust_map_path],
                                    distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def dust_map_earth(self):
        return self.dust_map # just the input map

    # -----------------------------------------------------------------

    @property
    def dust_map_faceon(self):
        return self.dust_projections.faceon

    # -----------------------------------------------------------------

    @property
    def dust_map_edgeon(self):
        return self.dust_projections.edgeon

    # -----------------------------------------------------------------

    @property
    def has_dust_map_earth(self):
        return self.has_dust_map

    # -----------------------------------------------------------------

    @property
    def has_dust_map_faceon(self):
        return self.has_dust_map # if there is an input map, we can deproject it

    # -----------------------------------------------------------------

    @property
    def has_dust_map_edgeon(self):
        return self.has_dust_map # if there is an input map, we can deproject it

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_component_sed(self):

        """
        This function ...
        :return:
        """

        return ComponentSED(bulge_component_name, self.old_bulge_component, description="old stellar bulge component",
                            path=self.old_bulge_sed_path, input_filepaths=self.old_bulge_input_filepaths,
                            distance=self.distance, inclination=self.inclination, position_angle=self.position_angle,
                            wavelengths_filename=wavelengths_filename)

    # -----------------------------------------------------------------

    @property
    def old_bulge_sed(self):
        return self.old_bulge_component_sed.sed

    # -----------------------------------------------------------------

    @property
    def old_bulge_sed_filepath(self):
        return self.old_bulge_component_sed.sed_filepath

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_component_sed(self):

        """
        This function ...
        :return:
        """

        return ComponentSED(disk_component_name, self.old_disk_component, description="old stellar disk component",
                            path=self.old_disk_sed_path, input_filepaths=self.old_disk_input_filepaths,
                            distance=self.distance, inclination=self.inclination, position_angle=self.position_angle,
                            wavelengths_filename=wavelengths_filename)

    # -----------------------------------------------------------------

    @property
    def old_disk_sed(self):
        return self.old_disk_component_sed.sed

    # -----------------------------------------------------------------

    @property
    def old_disk_sed_filepath(self):
        return self.old_disk_component_sed.sed_filepath

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def young_component_sed(self):

        """
        This function ...
        :return:
        """

        return ComponentSED(young_component_name, self.young_component, description="young stellar component",
                            path=self.young_sed_path, input_filepaths=self.young_input_filepaths,
                            distance=self.distance, inclination=self.inclination, position_angle=self.position_angle,
                            wavelengths_filename=wavelengths_filename)

    # -----------------------------------------------------------------

    @property
    def young_sed(self):
        return self.young_component_sed.sed

    # -----------------------------------------------------------------

    @property
    def young_sed_filepath(self):
        return self.young_component_sed.sed_filepath

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_component_sed(self):

        """
        This function ...
        :return:
        """

        return ComponentSED(ionizing_component_name, self.sfr_component, description=ionizing_component_description,
                            path=self.sfr_sed_path, input_filepaths=self.sfr_input_filepaths,
                            distance=self.distance, inclination=self.inclination, position_angle=self.position_angle,
                            wavelengths_filename=wavelengths_filename)

    # -----------------------------------------------------------------

    @property
    def sfr_sed(self):
        return self.sfr_component_sed.sed

    # -----------------------------------------------------------------

    @property
    def sfr_sed_filepath(self):
        return self.sfr_component_sed.sed_filepath

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def transparent_sfr_component_sed(self):

        """
        This function ...
        :return:
        """

        return ComponentSED(transparent_ionizing_component_name, self.transparent_sfr_component, description=transparent_ionizing_component_description,
                            path=self.transparent_sfr_sed_path, input_filepaths=self.sfr_input_filepaths,
                            distance=self.distance, inclination=self.inclination, position_angle=self.position_angle,
                            wavelengths_filename=wavelengths_filename)

    # -----------------------------------------------------------------

    @property
    def transparent_sfr_sed(self):
        return self.transparent_sfr_component_sed.sed

    # -----------------------------------------------------------------

    @property
    def transparent_sfr_sed_filepath(self):
        return self.transparent_sfr_component_sed.sed_filepath

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def extra_component_sed(self):

        """
        This function ...
        :return:
        """
        return ComponentSED(extra_component_name, self.extra_component, description="extra stellar component",
                            path=self.extra_sed_path, input_filepaths=self.extra_input_filepaths,
                            distance=self.distance, inclination=self.inclination, position_angle=self.position_angle,
                            wavelengths_filename=wavelengths_filename)

    # -----------------------------------------------------------------

    @property
    def extra_sed(self):
        return self.extra_component_sed.sed

    # -----------------------------------------------------------------

    @property
    def extra_sed_filepath(self):
        return self.extra_component_sed.sed_filepath

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def dust_sed(self):
        return self.total_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_dust_sed(self):
        return self.total_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def dust_luminosity(self):
        return self.total_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_dust_luminosity(self):
        return self.total_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_dust_sed_sfr(self):
        return self.sfr_simulations.intrinsic_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_dust_sed_sfr(self):
        return self.sfr_simulations.has_intrinsic_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_sed(self):
        return self.dust_sed - self.intrinsic_dust_sed_sfr

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_sed(self):
        return self.has_dust_sed and self.has_intrinsic_dust_sed_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_luminosity(self):
        return self.diffuse_dust_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_luminosity(self):
        return self.has_diffuse_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def absorbed_diffuse_stellar_luminosity(self):
        return self.total_absorbed_diffuse_stellar_sed_earth.integrate()

    # -----------------------------------------------------------------

    @property
    def has_absorbed_diffuse_stellar_luminosity(self):
        return self.has_total_absorbed_diffuse_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def diffuse_fabs(self):
        return self.total_fabs_diffuse_map.average()

    # -----------------------------------------------------------------

    @property
    def has_diffuse_fabs(self):
        return self.has_total_fabs_diffuse_map

    # -----------------------------------------------------------------

    @property
    def fabs(self):
        return self.total_fabs_map.average()

    # -----------------------------------------------------------------

    @property
    def has_fabs(self):
        return self.has_total_fabs_map

    # -----------------------------------------------------------------

    @lazyproperty
    def direct_stellar_luminosity(self):
        #return self.total_direct_stellar_luminosity_map.sum()
        return self.total_direct_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity(self):
        #return self.has_total_direct_stellar_luminosity_map
        return self.has_total_direct_stellar_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_fuv_luminosity(self):
        #return self.intrinsic_fuv_luminosity_map.sum(add_unit=True, per_area="error")
        return self.intrinsic_total_sed.photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity(self):
        return self.has_intrinsic_fuv_luminosity_map

    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_salim(self):
        return salim_fuv_to_sfr(self.intrinsic_fuv_luminosity, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_salim(self):
        return self.has_intrinsic_fuv_luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.intrinsic_fuv_luminosity, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_ke(self):
        return self.has_intrinsic_fuv_luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_tir(self):
        return kennicutt_tir_to_sfr(self.dust_luminosity, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_tir(self):
        return self.has_dust_luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_24um_luminosity(self):
        # don't interpolate, wavelength grid is expected to contain the FUV wavelength
        return self.total_simulations.observed_photometry_at(self.mips24_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_24um_luminosity(self):
        return self.total_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_24um(self):
        return calzetti_24um_to_sfr(self.observed_24um_luminosity, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_24um(self):
        return self.has_observed_24um_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_i1_luminosity(self):
        return self.observed_total_sed.photometry_at(self.i1_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity(self):
        return self.has_observed_total_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_i1_luminosity(self):
        return self.intrinsic_total_sed.photometry_at(self.i1_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_i1_luminosity(self):
        return self.has_intrinsic_total_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def total_stellar_mass(self):
        return oliver_stellar_mass(self.observed_i1_luminosity, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_stellar_mass(self):
        return self.has_observed_i1_luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def total_stellar_mass_intrinsic(self):
        return oliver_stellar_mass(self.intrinsic_i1_luminosity, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_total_stellar_mass_intrinsic(self):
        return self.has_intrinsic_i1_luminosity

    # -----------------------------------------------------------------

    @property
    def total_ssfr_salim(self):
        #return self.total_ssfr_map_salim.average(add_unit=True) # INCORRECT!
        return self.total_star_formation_rate_salim / self.total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_salim(self):
        #return self.has_total_ssfr_map_salim
        return self.has_total_star_formation_rate_salim and self.has_total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def total_ssfr_ke(self):
        #return self.total_ssfr_map_ke.average(add_unit=True) # INCORRECT!
        return self.total_star_formation_rate_ke / self.total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_ke(self):
        #return self.has_total_ssfr_map_ke
        return self.has_total_star_formation_rate_ke and self.has_total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def total_ssfr_tir(self):
        return self.total_star_formation_rate_tir / self.total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_tir(self):
        return self.has_total_star_formation_rate_tir and self.has_total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def total_ssfr_24um(self):
        return self.total_star_formation_rate_24um / self.total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr_24um(self):
        return self.has_total_star_formation_rate_24um and self.has_total_stellar_mass

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_total(self):

        """
        Total
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Bolometric luminosity (2 values should be the same)
        if self.has_observed_total_bolometric_luminosity: values[obs_total_bol_lum_name] = self.observed_total_bolometric_luminosity
        if self.has_intrinsic_total_bolometric_luminosity: values[intr_total_bol_lum_name] = self.intrinsic_total_bolometric_luminosity

        # Stellar luminosity
        if self.has_intrinsic_stellar_luminosity: values[intr_total_stellar_bol_lum_name] = self.intrinsic_stellar_luminosity
        if self.has_observed_stellar_luminosity: values[obs_total_stellar_bol_lum_name] = self.observed_stellar_luminosity

        # Dust luminosity
        if self.has_diffuse_dust_luminosity: values[diffuse_dust_lum_name] = self.diffuse_dust_luminosity
        if self.has_dust_luminosity: values[dust_lum_name] = self.dust_luminosity

        # Absorbed luminosity
        if self.has_absorbed_diffuse_stellar_luminosity: values[diffuse_abs_stellar_lum_name] = self.absorbed_diffuse_stellar_luminosity

        # Fabs
        if self.has_diffuse_fabs: values[diffuse_fabs_name] = self.diffuse_fabs
        if self.has_fabs: values[fabs_name] = self.fabs

        # Total attenuation
        if self.has_bolometric_attenuation: values[bol_attenuation_name] = self.bolometric_attenuation

        # Direct stellar luminosity
        if self.has_direct_stellar_luminosity: values[direct_stellar_lum_name] = self.direct_stellar_luminosity

        # Star formation rate
        if self.has_total_star_formation_rate_salim: values[sfr_salim_name] = self.total_star_formation_rate_salim
        if self.has_total_star_formation_rate_ke: values[sfr_ke_name] = self.total_star_formation_rate_ke
        if self.has_total_star_formation_rate_tir: values[sfr_tir_name] = self.total_star_formation_rate_tir
        if self.has_total_star_formation_rate_24um: values[sfr_24um_name] = self.total_star_formation_rate_24um

        # Stellar mass
        if self.has_total_stellar_mass: values[stellar_mass_name] = self.total_stellar_mass
        if self.has_total_stellar_mass_intrinsic: values[stellar_mass_intrinsic_name] = self.total_stellar_mass_intrinsic

        # Specific star formation rate
        if self.has_total_ssfr_salim: values[ssfr_salim_name] = self.total_ssfr_salim
        if self.has_total_ssfr_ke: values[ssfr_ke_name] = self.total_ssfr_ke
        if self.has_total_ssfr_tir: values[ssfr_tir_name] = self.total_ssfr_tir
        if self.has_total_ssfr_24um: values[ssfr_24um_name] = self.total_ssfr_24um

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def direct_stellar_luminosity_old_bulge(self):
        #return self.old_bulge_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")
        return self.old_bulge_direct_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_old_bulge(self):
        #return self.has_old_bulge_direct_stellar_luminosity_map
        return self.has_old_bulge_direct_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_stellar_mass(self):
        return oliver_stellar_mass(self.observed_i1_luminosity_old_bulge, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_bulge_stellar_mass(self):
        return self.has_observed_i1_luminosity_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_stellar_mass_intrinsic(self):
        return oliver_stellar_mass(self.intrinsic_i1_luminosity_old_bulge, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_bulge_stellar_mass_intrinsic(self):
        return self.has_intrinsic_i1_luminosity_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_bulge(self):

        """
        Old bulge
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Bolometric luminosity (2 values should be the same)
        if self.has_observed_bolometric_luminosity_old_bulge: values[obs_bulge_bol_lum_name] = self.observed_bolometric_luminosity_old_bulge
        if self.has_intrinsic_bolometric_luminosity_old_bulge: values[intr_bulge_bol_lum_name] = self.intrinsic_bolometric_luminosity_old_bulge

        # Direct stellar luminosity
        if self.has_direct_stellar_luminosity_old_bulge: values[direct_stellar_lum_name] = self.direct_stellar_luminosity_old_bulge

        # I1 specific luminosity
        if self.has_observed_i1_luminosity_old_bulge: values[obs_bulge_spec_lum_name] = self.observed_i1_luminosity_old_bulge
        if self.has_intrinsic_i1_luminosity_old_bulge: values[intr_bulge_spec_lum_name] = self.intrinsic_i1_luminosity_old_bulge  # part of parameter set

        # FUV specific luminosity
        if self.has_observed_fuv_luminosity_old_bulge: values[obs_fuv_spec_lum_name] = self.observed_fuv_luminosity_old_bulge
        if self.has_intrinsic_fuv_luminosity_old_bulge: values[intr_fuv_spec_lum_name] = self.intrinsic_fuv_luminosity_old_bulge

        # Stellar mass
        if self.has_bulge_stellar_mass: values[stellar_mass_name] = self.bulge_stellar_mass
        if self.has_bulge_stellar_mass_intrinsic: values[stellar_mass_intrinsic_name] = self.bulge_stellar_mass_intrinsic

        # Attenuation
        if self.has_i1_attenuation_old_bulge: values[bulge_spec_attenuation_name] = self.i1_attenuation_old_bulge
        if self.has_bolometric_attenuation_old_bulge: values[bulge_bol_attenuation_name] = self.bolometric_attenuation_old_bulge

        # Dust
        if self.has_observed_dust_luminosity_old_bulge: values[obs_bulge_dust_lum_name] = self.observed_dust_luminosity_old_bulge

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def direct_stellar_luminosity_old_disk(self):
        #return self.old_disk_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")
        return self.old_disk_direct_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_old_disk(self):
        #return self.has_old_disk_direct_stellar_luminosity_map
        return self.has_old_disk_direct_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_stellar_mass(self):
        return oliver_stellar_mass(self.observed_i1_luminosity_old_disk, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_disk_stellar_mass(self):
        return self.has_observed_i1_luminosity_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_stellar_mass_intrinsic(self):
        return oliver_stellar_mass(self.intrinsic_i1_luminosity_old_disk, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_disk_stellar_mass_intrinsic(self):
        return self.has_intrinsic_i1_luminosity_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_disk(self):

        """
        Old disk
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Bolometric luminosity (2 values should be the same)
        if self.has_observed_bolometric_luminosity_old_disk: values[obs_disk_bol_lum_name] = self.observed_bolometric_luminosity_old_disk
        if self.has_intrinsic_bolometric_luminosity_old_disk: values[intr_disk_bol_lum_name] = self.intrinsic_bolometric_luminosity_old_disk

        # Direct stellar luminosity
        if self.has_direct_stellar_luminosity_old_disk: values[direct_stellar_lum_name] = self.direct_stellar_luminosity_old_disk

        # I1 specific luminosity
        if self.has_observed_i1_luminosity_old_disk: values[obs_disk_spec_lum_name] = self.observed_i1_luminosity_old_disk
        if self.has_intrinsic_i1_luminosity_old_disk: values[intr_disk_spec_lum_name] = self.intrinsic_i1_luminosity_old_disk # part of parameter set

        # FUV specific luminosity
        if self.has_observed_fuv_luminosity_old_disk: values[obs_fuv_spec_lum_name] = self.observed_fuv_luminosity_old_disk
        if self.has_intrinsic_fuv_luminosity_old_disk: values[intr_fuv_spec_lum_name] = self.intrinsic_fuv_luminosity_old_disk

        # Stellar mass
        if self.has_disk_stellar_mass: values[stellar_mass_name] = self.disk_stellar_mass
        if self.has_disk_stellar_mass_intrinsic: values[stellar_mass_intrinsic_name] = self.disk_stellar_mass_intrinsic

        # Attenuation
        if self.has_i1_attenuation_old_disk: values[disk_spec_attenuation_name] = self.i1_attenuation_old_disk
        if self.has_bolometric_attenuation_old_disk: values[disk_bol_attenuation_name] = self.bolometric_attenuation_old_disk

        # Dust
        if self.has_observed_dust_luminosity_old_disk: values[obs_disk_dust_lum_name] = self.observed_dust_luminosity_old_disk

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def direct_stellar_luminosity_old(self):
        #return self.old_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")
        return self.old_direct_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_old(self):
        #return self.has_old_direct_stellar_luminosity_map
        return self.has_old_direct_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stellar_mass(self):
        return oliver_stellar_mass(self.observed_i1_luminosity_old, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_old_stellar_mass(self):
        return self.has_observed_i1_luminosity_old

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stellar_mass_intrinsic(self):
        return oliver_stellar_mass(self.intrinsic_i1_luminosity_old, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_old_stellar_mass_intrinsic(self):
        return self.has_intrinsic_i1_luminosity_old

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_old(self):

        """
        Old (evolved)
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Bolometric luminosity (2 values should be the same)
        if self.has_observed_bolometric_luminosity_old: values[obs_old_bol_lum_name] = self.observed_bolometric_luminosity_old
        if self.has_intrinsic_bolometric_luminosity_old: values[intr_old_bol_lum_name] = self.intrinsic_bolometric_luminosity_old

        # Direct stellar luminosity
        if self.has_direct_stellar_luminosity_old: values[direct_stellar_lum_name] = self.direct_stellar_luminosity_old

        # I1 specific luminosity
        if self.has_observed_i1_luminosity_old: values[obs_old_spec_lum_name] = self.observed_i1_luminosity_old
        if self.has_intrinsic_i1_luminosity_old: values[intr_old_spec_lum_name] = self.intrinsic_i1_luminosity_old

        # Stellar mass
        if self.has_old_stellar_mass: values[stellar_mass_name] = self.old_stellar_mass
        if self.has_old_stellar_mass_intrinsic: values[stellar_mass_intrinsic_name] = self.old_stellar_mass_intrinsic

        # Attenuation
        if self.has_i1_attenuation_old: values[old_spec_attenuation_name] = self.i1_attenuation_old
        if self.has_bolometric_attenuation_old: values[old_bol_attenuation_name] = self.bolometric_attenuation_old

        # Dust
        if self.has_observed_dust_luminosity_old: values[obs_old_dust_lum_name] = self.observed_dust_luminosity_old

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def direct_stellar_luminosity_young(self):
        #return self.young_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")
        return self.young_direct_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_young(self):
        #return self.has_young_direct_stellar_luminosity_map
        return self.has_young_direct_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_salim(self):
        return salim_fuv_to_sfr(self.intrinsic_fuv_luminosity_young, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_salim(self):
        return self.has_intrinsic_fuv_luminosity_young

    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.intrinsic_fuv_luminosity_young, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_ke(self):
        return self.has_intrinsic_fuv_luminosity_young

    # -----------------------------------------------------------------

    @property
    def young_ssfr_salim(self):
        return self.young_star_formation_rate_salim / self.total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_young_ssfr_salim(self):
        return self.has_young_star_formation_rate_salim and self.has_total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def young_ssfr_ke(self):
        return self.young_star_formation_rate_ke / self.total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_young_ssfr_ke(self):
        return self.has_young_star_formation_rate_ke and self.has_total_stellar_mass

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_young(self):

        """
        Young stars
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Bolometric luminosity (2 values should be the same)
        if self.has_observed_bolometric_luminosity_young: values[obs_young_bol_lum_name] = self.observed_bolometric_luminosity_young
        if self.has_intrinsic_bolometric_luminosity_young: values[intr_young_bol_lum_name] = self.intrinsic_bolometric_luminosity_young

        # Direct stellar luminosity
        if self.has_direct_stellar_luminosity_young: values[direct_stellar_lum_name] = self.direct_stellar_luminosity_young

        # FUV specific luminosity
        if self.has_observed_fuv_luminosity_young: values[obs_young_spec_lum_name] = self.observed_fuv_luminosity_young
        if self.has_intrinsic_fuv_luminosity_young: values[intr_young_spec_lum_name] = self.intrinsic_fuv_luminosity_young # part of (free) parameter set

        # Star formation rate
        if self.has_young_star_formation_rate_salim: values[sfr_salim_name] = self.young_star_formation_rate_salim
        if self.has_young_star_formation_rate_ke: values[sfr_ke_name] = self.young_star_formation_rate_ke

        # Specific star formation rate
        if self.has_young_ssfr_salim: values[ssfr_salim_name] = self.young_ssfr_salim
        if self.has_young_ssfr_ke: values[ssfr_ke_name] = self.young_ssfr_ke

        # Attenuation
        if self.has_fuv_attenuation_young: values[young_spec_attenuation_name] = self.fuv_attenuation_young
        if self.has_bolometric_attenuation_young: values[young_bol_attenuation_name] = self.bolometric_attenuation_young

        # Dust
        if self.has_observed_dust_luminosity_young: values[obs_young_dust_lum_name] = self.observed_dust_luminosity_young

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def direct_stellar_luminosity_sfr(self):
        #return self.sfr_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")
        return self.sfr_direct_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_sfr(self):
        #return self.has_sfr_direct_stellar_luminosity_map
        return self.has_sfr_direct_stellar_sed

    # -----------------------------------------------------------------

    @property
    def sfr_star_formation_rate_mappings(self):
        return self.sfr

    # -----------------------------------------------------------------

    @property
    def has_sfr_star_formation_rate_mappings(self):
        return self.has_sfr

    # -----------------------------------------------------------------

    @property
    def sfr_ssfr_mappings(self):
        return self.sfr_star_formation_rate_mappings / self.total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_sfr_ssfr_mappings(self):
        return self.has_sfr_star_formation_rate_mappings and self.has_total_stellar_mass

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_sfr(self):

        """
        Ionizing stars (SFR)
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Bolometric luminosity (2 values should be the same)
        if self.has_observed_bolometric_luminosity_sfr: values[obs_sfr_bol_lum_name] = self.observed_bolometric_luminosity_sfr
        if self.has_intrinsic_bolometric_luminosity_sfr: values[intr_sfr_bol_lum_name] = self.intrinsic_bolometric_luminosity_sfr

        # Direct stellar luminosity
        if self.has_direct_stellar_luminosity_sfr: values[direct_stellar_lum_name] = self.direct_stellar_luminosity_sfr

        # FUV specific luminosity
        if self.has_observed_fuv_luminosity_sfr: values[obs_sfr_spec_lum_name] = self.observed_fuv_luminosity_sfr
        if self.has_intrinsic_fuv_luminosity_sfr: values[intr_sfr_spec_lum_name] = self.intrinsic_fuv_luminosity_sfr # part of the (free) parameter set

        # SFR
        if self.has_sfr_star_formation_rate_mappings: values[sfr_mappings_name] = self.sfr_star_formation_rate_mappings

        # Specific star formation rate
        if self.has_sfr_ssfr_mappings: values[ssfr_mappings_name] = self.sfr_ssfr_mappings

        # Attenuation
        if self.has_fuv_attenuation_sfr: values[sfr_spec_attenuation_name] = self.fuv_attenuation_sfr
        if self.has_bolometric_attenuation_sfr: values[sfr_bol_attenuation_name] = self.bolometric_attenuation_sfr

        # Stellar
        if self.has_sfr_stellar_mass: values[sfr_stellar_mass_name] = self.sfr_stellar_mass
        if self.has_observed_stellar_luminosity_sfr: values[obs_sfr_stellar_bol_lum_name] = self.observed_stellar_luminosity_sfr
        if self.has_intrinsic_stellar_luminosity_sfr: values[intr_sfr_stellar_bol_lum_name] = self.intrinsic_stellar_luminosity_sfr

        # Dust
        if self.has_sfr_dust_mass: values[sfr_dust_mass_name] = self.sfr_dust_mass
        if self.has_observed_dust_luminosity_sfr: values[obs_sfr_dust_lum_name] = self.observed_dust_luminosity_sfr
        if self.has_intrinsic_dust_luminosity_sfr: values[sfr_dust_lum_name] = self.intrinsic_dust_luminosity_sfr # intrinsic so only the dust in MAPPINGS

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def direct_stellar_luminosity_unevolved(self):
        #return self.unevolved_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")
        return self.unevolved_direct_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_unevolved(self):
        #return self.has_unevolved_direct_stellar_luminosity_map
        return self.has_unevolved_direct_stellar_sed

    # -----------------------------------------------------------------

    @property
    def unevolved_intrinsic_fuv_luminosity(self):
        return self.intrinsic_fuv_luminosity_unevolved

    # -----------------------------------------------------------------

    @property
    def has_unevolved_intrinsic_fuv_luminosity(self):
        return self.has_intrinsic_fuv_luminosity_unevolved

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_star_formation_rate_salim(self):
        return salim_fuv_to_sfr(self.unevolved_intrinsic_fuv_luminosity, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_star_formation_rate_salim(self):
        return self.has_unevolved_intrinsic_fuv_luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_star_formation_rate_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.unevolved_intrinsic_fuv_luminosity, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_star_formation_rate_ke(self):
        return self.has_unevolved_intrinsic_fuv_luminosity

    # -----------------------------------------------------------------

    @property
    def unevolved_star_formation_rate_mappings_ke(self):
        return self.sfr_star_formation_rate_mappings + self.young_star_formation_rate_ke

    # -----------------------------------------------------------------

    @property
    def has_unevolved_star_formation_rate_mappings_ke(self):
        return self.has_sfr_star_formation_rate_mappings and self.has_young_star_formation_rate_ke

    # -----------------------------------------------------------------

    @property
    def unevolved_ssfr_salim(self):
        return self.unevolved_star_formation_rate_salim / self.total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_unevolved_ssfr_salim(self):
        return self.has_unevolved_star_formation_rate_salim and self.has_total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def unevolved_ssfr_ke(self):
        return self.unevolved_star_formation_rate_ke / self.total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_unevolved_ssfr_ke(self):
        return self.has_unevolved_star_formation_rate_ke and self.has_total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def unevolved_ssfr_mappings_ke(self):
        return self.unevolved_star_formation_rate_mappings_ke / self.total_stellar_mass

    # -----------------------------------------------------------------

    @property
    def has_unevolved_ssfr_mappings_ke(self):
        return self.has_unevolved_star_formation_rate_mappings_ke and self.has_total_stellar_mass

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_unevolved(self):

        """
        Young + ionizing (unevolved)
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Bolometric luminosity
        if self.has_observed_bolometric_luminosity_unevolved: values[obs_unevolved_bol_lum_name] = self.observed_bolometric_luminosity_unevolved
        if self.has_intrinsic_bolometric_luminosity_unevolved: values[intr_unevolved_bol_lum_name] = self.intrinsic_bolometric_luminosity_unevolved

        # Direct
        if self.has_direct_stellar_luminosity_unevolved: values[direct_stellar_lum_name] = self.direct_stellar_luminosity_unevolved

        # FUV specific luminosity
        if self.has_observed_fuv_luminosity_unevolved: values[obs_unevolved_spec_lum_name] = self.observed_fuv_luminosity_unevolved
        if self.has_intrinsic_fuv_luminosity_unevolved: values[intr_unevolved_spec_lum_name] = self.intrinsic_fuv_luminosity_unevolved

        # Star formation rate
        if self.has_unevolved_star_formation_rate_salim: values[sfr_salim_name] = self.unevolved_star_formation_rate_salim
        if self.has_unevolved_star_formation_rate_ke: values[sfr_ke_name] = self.unevolved_star_formation_rate_ke
        if self.has_unevolved_star_formation_rate_mappings_ke: values[sfr_mappings_ke_name] = self.unevolved_star_formation_rate_mappings_ke

        # Specific star formation rate
        if self.has_unevolved_ssfr_salim: values[ssfr_salim_name] = self.unevolved_ssfr_salim
        if self.has_unevolved_ssfr_ke: values[ssfr_ke_name] = self.unevolved_ssfr_ke
        if self.has_unevolved_ssfr_mappings_ke: values[ssfr_mappings_ke_name] = self.unevolved_ssfr_mappings_ke

        # Attenuation
        if self.has_fuv_attenuation_unevolved: values[unevolved_spec_attenuation_name] = self.fuv_attenuation_unevolved
        if self.has_bolometric_attenuation_unevolved: values[unevolved_bol_attenuation_name] = self.bolometric_attenuation_unevolved

        # Dust
        if self.has_observed_dust_luminosity_unevolved: values[obs_unevolved_dust_lum_name] = self.observed_dust_luminosity_unevolved

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

        # Mass
        if self.has_total_dust_mass: values[total_dust_mass_name] = self.total_dust_mass # with SFR dust mass

        # Diffuse mass
        if self.has_diffuse_dust_mass: values[diffuse_dust_mass_name] = self.diffuse_dust_mass # ACTUALLY ONE OF THE INTRINSIC PARAMETERS

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
        for key in self.derived_parameter_values_total: values[key + total_suffix] = self.derived_parameter_values_total[key]
        for key in self.derived_parameter_values_bulge: values[key + bulge_suffix] = self.derived_parameter_values_bulge[key]
        for key in self.derived_parameter_values_disk: values[key + disk_suffix] = self.derived_parameter_values_disk[key]
        for key in self.derived_parameter_values_old: values[key + old_suffix] = self.derived_parameter_values_old[key]
        for key in self.derived_parameter_values_young: values[key + young_suffix] = self.derived_parameter_values_young[key]
        for key in self.derived_parameter_values_sfr: values[key + sfr_suffix] = self.derived_parameter_values_sfr[key]
        for key in self.derived_parameter_values_unevolved: values[key + unevolved_suffix] = self.derived_parameter_values_unevolved[key]
        for key in self.derived_parameter_values_dust: values[key + dust_suffix] = self.derived_parameter_values_dust[key]

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Total
        if self.has_total_star_formation_rate_salim: values[sfr_salim_name + total_suffix] = self.total_star_formation_rate_salim
        if self.has_total_star_formation_rate_ke: values[sfr_ke_name + total_suffix] = self.total_star_formation_rate_ke
        if self.has_total_star_formation_rate_tir: values[sfr_tir_name + total_suffix] = self.total_star_formation_rate_tir
        if self.has_total_star_formation_rate_24um: values[sfr_24um_name + total_suffix] = self.total_star_formation_rate_24um

        # Young
        if self.has_young_star_formation_rate_salim: values[sfr_salim_name + young_suffix] = self.young_star_formation_rate_salim
        if self.has_young_star_formation_rate_ke: values[sfr_ke_name + young_suffix] = self.young_star_formation_rate_ke

        # Sfr
        if self.has_sfr_star_formation_rate_mappings: values[sfr_mappings_name + sfr_suffix] = self.sfr_star_formation_rate_mappings

        # Unevolved
        if self.has_unevolved_star_formation_rate_salim: values[sfr_salim_name + unevolved_suffix] = self.unevolved_star_formation_rate_salim
        if self.has_unevolved_star_formation_rate_ke: values[sfr_ke_name + unevolved_suffix] = self.unevolved_star_formation_rate_ke
        if self.has_unevolved_star_formation_rate_mappings_ke: values[sfr_mappings_ke_name + unevolved_suffix] = self.unevolved_star_formation_rate_mappings_ke

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Total
        if self.has_total_ssfr_salim: values[ssfr_salim_name + total_suffix] = self.total_ssfr_salim
        if self.has_total_ssfr_ke: values[ssfr_ke_name + total_suffix] = self.total_ssfr_ke
        if self.has_total_ssfr_tir: values[ssfr_tir_name + total_suffix] = self.total_ssfr_tir
        if self.has_total_ssfr_24um: values[ssfr_24um_name + total_suffix] = self.total_ssfr_24um

        # Young
        if self.has_young_ssfr_salim: values[ssfr_salim_name + young_suffix] = self.young_ssfr_salim
        if self.has_young_ssfr_ke: values[ssfr_ke_name + young_suffix] = self.young_ssfr_ke

        # Sfr
        if self.has_sfr_ssfr_mappings: values[ssfr_mappings_name + sfr_suffix] = self.sfr_ssfr_mappings

        # Unevolved
        if self.has_unevolved_ssfr_salim: values[ssfr_salim_name + unevolved_suffix] = self.unevolved_ssfr_salim
        if self.has_unevolved_ssfr_ke: values[ssfr_ke_name + unevolved_suffix] = self.unevolved_ssfr_ke
        if self.has_unevolved_ssfr_mappings_ke: values[ssfr_mappings_ke_name + unevolved_suffix] = self.unevolved_ssfr_mappings_ke

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_mass_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Total
        if self.has_total_stellar_mass: values[stellar_mass_name + total_suffix] = self.total_stellar_mass
        if self.has_total_stellar_mass_intrinsic: values[stellar_mass_intrinsic_name + total_suffix] = self.total_stellar_mass_intrinsic

        # Bulge
        if self.has_bulge_stellar_mass: values[stellar_mass_name + bulge_suffix] = self.bulge_stellar_mass
        if self.has_bulge_stellar_mass_intrinsic: values[stellar_mass_intrinsic_name + bulge_suffix] = self.bulge_stellar_mass_intrinsic

        # Disk
        if self.has_disk_stellar_mass: values[stellar_mass_name + disk_suffix] = self.disk_stellar_mass
        if self.has_disk_stellar_mass_intrinsic: values[stellar_mass_intrinsic_name + disk_suffix] = self.disk_stellar_mass_intrinsic

        # Old
        if self.has_old_stellar_mass: values[stellar_mass_name + old_suffix] = self.old_stellar_mass
        if self.has_old_stellar_mass_intrinsic: values[stellar_mass_intrinsic_name + old_suffix] = self.old_stellar_mass_intrinsic

        # Return
        return values

    # -----------------------------------------------------------------
    # DUST GRID
    # -----------------------------------------------------------------

    @property
    def has_old_bulge_cell_stellar_density(self):
        return self.bulge_simulations.has_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def old_bulge_cell_stellar_density(self):
        return self.bulge_simulations.cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def has_old_disk_cell_stellar_density(self):
        return self.disk_simulations.has_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def old_disk_cell_stellar_density(self):
        return self.disk_simulations.cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def has_young_cell_stellar_density(self):
        return self.young_simulations.has_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def young_cell_stellar_density(self):
        return self.young_simulations.cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def has_sfr_cell_stellar_density(self):
        return self.sfr_simulations.has_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def sfr_cell_stellar_density(self):
        return self.sfr_simulations.cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def has_cell_properties(self):
        return self.total_simulations.has_cell_properties

    # -----------------------------------------------------------------

    @property
    def cell_volumes(self):
        return self.total_simulations.cell_volumes

    # -----------------------------------------------------------------

    @property
    def cell_volume_unit(self):
        return self.total_simulations.cell_volumes_unit

    # -----------------------------------------------------------------

    @property
    def cell_dust_densities(self):
        return self.total_simulations.cell_dust_densities

    # -----------------------------------------------------------------

    @property
    def cell_mass_fractions(self):
        return self.total_simulations.cell_mass_fractions

    # -----------------------------------------------------------------

    @property
    def cell_optical_depths(self):
        return self.total_simulations.cell_optical_depths

    # -----------------------------------------------------------------

    @property
    def cell_masses(self):
        return self.total_simulations.cell_masses

    # -----------------------------------------------------------------

    @property
    def cell_mass_unit(self):
        return self.total_simulations.cell_mass_unit

    # -----------------------------------------------------------------

    @property
    def cell_temperatures(self):
        return self.total_simulations.cell_temperatures

    # -----------------------------------------------------------------

    @property
    def cell_temperature_unit(self):
        return self.total_simulations.cell_temperature_unit

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_x_coordinates(self):
        return self.total_simulations.cell_x_coordinates

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_y_coordinates(self):
        return self.total_simulations.cell_y_coordinates

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_z_coordinates(self):
        return self.total_simulations.cell_z_coordinates

    # -----------------------------------------------------------------

    @property
    def has_grid_files(self):
        return self.total_simulations.has_grid_files

    # -----------------------------------------------------------------

    @property
    def grid_filepaths(self):
        return self.total_simulations.grid_filepaths

    # -----------------------------------------------------------------

    @property
    def grid_xy_filepath(self):
        return self.total_simulations.grid_xy_filepath

    # -----------------------------------------------------------------

    @property
    def grid_xz_filepath(self):
        return self.total_simulations.grid_xz_filepath

    # -----------------------------------------------------------------

    @property
    def grid_yz_filepath(self):
        return self.total_simulations.grid_yz_filepath

    # -----------------------------------------------------------------

    @property
    def grid_xyz_filepath(self):
        return self.total_simulations.grid_xyz_filepath

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def get_observed_stellar_sed(self, component):

        """
        This function ...
        :param component:
        :return:
        """

        # Get the simulations
        simulations = self.simulations[component]

        # Return the SED
        if simulations.has_observed_stellar_sed: return simulations.observed_stellar_sed
        else: return None

    # -----------------------------------------------------------------

    def get_intrinsic_stellar_sed(self, component):

        """
        This function ...
        :param component:
        :return:
        """

        # Get the simulations
        simulations = self.simulations[component]

        # Return the SED
        if simulations.has_intrinsic_sed: return simulations.intrinsic_sed
        else: return None

    # -----------------------------------------------------------------

    def get_stellar_sed(self, component, observed_intrinsic):

        """
        This function ...
        :param component:
        :param observed_intrinsic:
        :return:
        """

        # Observed
        if observed_intrinsic == observed_name: return self.get_observed_stellar_sed(component)

        # Intrinsic
        elif observed_intrinsic == intrinsic_name: return self.get_intrinsic_stellar_sed(component)

        # Invalid
        else: raise ValueError("Invalid option for 'observed_intrinsic': '" + observed_intrinsic + "'")

    # -----------------------------------------------------------------

    def get_dust_sed(self, component):

        """
        This function ...
        :param component:
        :return:
        """

        # Return
        return self.simulations[component].observed_dust_sed

        #if component == total: return self.dust_sed
        #elif component == bulge: return self.observed_old_bulge_dust_sed
        #elif component == disk: return self.observed_old_disk_dust_sed
        #elif component == old: return self.observed_old_dust_sed
        #elif component == young: return self.observed_young_dust_sed
        #elif component == sfr: return self.observed_sfr_dust_sed
        #elif component == unevolved: return self.observed_unevolved_dust_sed
        #else: raise ValueError("Invalid component: '" + component + "'")

        # THERE IS ALSO model.diffuse_dust_sed !

    # -----------------------------------------------------------------

    @lazyproperty
    def component_maps(self):
        return self.component_maps_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def component_maps_earth(self):
        maps = OrderedDict()
        maps[bulge_simulation_name] = self.old_bulge_map_earth
        maps[disk_simulation_name] = self.old_disk_map_earth
        maps[young_simulation_name] = self.young_map_earth
        maps[sfr_simulation_name] = self.sfr_map_earth
        maps[dust_simulation_name] = self.dust_map_earth
        return maps

    # -----------------------------------------------------------------

    @lazyproperty
    def component_maps_faceon(self):
        maps = OrderedDict()
        maps[bulge_simulation_name] = self.old_bulge_map_faceon
        maps[disk_simulation_name] = self.old_disk_map_faceon
        maps[young_simulation_name] = self.young_map_faceon
        maps[sfr_simulation_name] = self.sfr_map_faceon
        maps[dust_simulation_name] = self.dust_map_faceon
        return maps

    # -----------------------------------------------------------------

    @lazyproperty
    def component_maps_edgeon(self):
        maps = OrderedDict()
        maps[bulge_simulation_name] = self.old_bulge_map_edgeon
        maps[disk_simulation_name] = self.old_disk_map_edgeon
        maps[young_simulation_name] = self.young_map_edgeon
        maps[sfr_simulation_name] = self.sfr_map_edgeon
        maps[dust_simulation_name] = self.dust_map_edgeon
        return maps

    # -----------------------------------------------------------------

    @lazyproperty
    def component_names(self):
        names = OrderedDict()
        names[bulge_simulation_name] = bulge_component_name
        names[disk_simulation_name] = disk_component_name
        names[young_simulation_name] = young_component_name
        names[sfr_simulation_name] = ionizing_component_name
        names[dust_simulation_name] = dust_component_name
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def component_descriptions(self):
        descriptions = OrderedDict()
        descriptions[bulge_simulation_name] = bulge_component_description
        descriptions[disk_simulation_name] = disk_component_description
        descriptions[young_simulation_name] = young_component_description
        descriptions[sfr_simulation_name] = ionizing_component_description
        descriptions[dust_simulation_name] = dust_component_description
        return descriptions

# -----------------------------------------------------------------

def correct_sed_for_attenuation(sed, attenuation_curve):

    """
    Thisf unction ...
    :param sed:
    :param attenuation_curve:
    :return:
    """

    from scipy.interpolate import interp1d

    # Get attenuation data
    wavelengths = attenuation_curve.wavelengths(unit="micron", asarray=True)
    attenuations = attenuation_curve.attenuations(asarray=True)

    # Interpolate attenuation data on the SED wavelengths
    interpolated = interp1d(wavelengths, attenuations)
    sed_wavelengths = sed.wavelengths(unit="micron", asarray=True)
    sed_attenuations = interpolated(sed_wavelengths)

    # Get transparent photometry
    photometry = sed.photometry(unit=sed.unit, asarray=True)
    exponent = sed_attenuations / 2.5
    transparent_photometry = photometry * 10 ** exponent

    # Create SED
    transparent_sed = SED.from_arrays(sed_wavelengths, transparent_photometry, "micron", sed.unit)

    # Return
    return transparent_sed

# -----------------------------------------------------------------
