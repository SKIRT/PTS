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
from ...core.units.parsing import parse_unit as u
from ...core.tools import types

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
sfr_name = "Star formation rate"
stellar_mass_name = "Stellar mass"
ssfr_name = "Specific star formation rate"

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

## Dust
total_dust_mass_name = "Total dust mass" # 6 with SFR dust mass
diffuse_dust_mass_name = "Diffuse dust mass"

# -----------------------------------------------------------------

sed_dirname = "sed"
projections_dirname = "projections"

# -----------------------------------------------------------------

total_simulation_name = "total"
bulge_simulation_name = "bulge"
disk_simulation_name = "disk"
old_simulation_name = "old"
young_simulation_name = "young"
sfr_simulation_name = "sfr"
unevolved_simulation_name = "unevolved"
dust_simulation_name = "dust"

bulge_component_name = "Evolved stellar bulge"
disk_component_name = "Evolved stellar disk"
young_component_name = "Young stars"
ionizing_component_name = "Ionizing stars"

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

class RTModel(object):

    """
    Objects of this class describe a full radiative transfer model.
    """

    def __init__(self, definition, wavelength_grid=None, simulation_name=None, chi_squared=None,
                 free_parameter_labels=None, free_parameter_values=None, observed_total_output_path=None,
                 observed_bulge_output_path=None, observed_disk_output_path=None, observed_old_output_path=None,
                 observed_young_output_path=None, observed_sfr_output_path=None, observed_unevolved_output_path=None,
                 parameters=None, center=None, galaxy_name=None, hubble_stage=None, redshift=None, earth_wcs=None):

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

        # Set parameters
        if parameters is not None: self.set_parameters(**parameters)

        # Set the center and other galaxy properties
        self.center = center
        self.galaxy_name = galaxy_name
        self.hubble_stage = hubble_stage
        self.redshift = redshift

        # Set the earth coordinate system
        self.earth_wcs = earth_wcs

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

        # Return the seds
        return seds

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_component_cubes(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        cubes = OrderedDict()

        ## EVOLVED

        # Bulge & disk?
        if self.bulge_simulations.has_intrinsic_cube and self.disk_simulations.has_intrinsic_cube:

            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_earth
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_earth

        # Old?
        elif self.old_simulations.has_intrinsic_cube: cubes[evolved_component_name] = self.old_intrinsic_stellar_luminosity_cube_earth

        # Not enough data
        else: #raise ValueError("Not enough simulation data")
            warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the total simulation, the transparent earth cube will not be available")
            return None

        ## UNEVOLVED

        # Young & ionizing?
        if self.young_simulations.has_intrinsic_cube and self.sfr_simulations.has_intrinsic_cube:

            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_earth
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_earth

        # Unevolved?
        elif self.unevolved_simulations.has_intrinsic_cube: cubes[unevolved_component_name] = self.unevolved_intrinsic_stellar_luminosity_cube_earth

        # Not enough data
        else: #raise ValueError("Not enough simulation data")
            warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the total simulation, the transparent earth cube will not be available")
            return None

        # Return
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_component_cubes_faceon(self):

        """
        This function ...
        :return:
        """

        # Initialie dictionary
        cubes = OrderedDict()

        ## EVOLVED

        # Bulge @ disk?
        if self.bulge_simulations.has_intrinsic_cube_faceon and self.disk_simulations.has_intrinsic_cube_faceon:

            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_faceon
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_faceon

        # Old?
        elif self.old_simulations.has_intrinsic_cube_faceon: cubes[evolved_component_name] = self.old_intrinsic_stellar_luminosity_cube_faceon

        # Not enough data
        else:
            warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the total simulation, the transparent faceon cube will not be available")
            return None

        ## UNEVOLVED

        # Young & ionizing?
        if self.young_simulations.has_intrinsic_cube_faceon and self.sfr_simulations.has_intrinsic_cube_faceon:

            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_faceon
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_faceon

        # Unevolved?
        elif self.unevolved_simulations.has_intrinsic_cube_faceon: cubes[unevolved_component_name] = self.unevolved_intrinsic_stellar_luminosity_cube_faceon

        # Not enough data
        else:
            warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the total simulation, the transparent faceon cube will not be available")
            return None

        # Return
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation_component_cubes_edgeon(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        cubes = OrderedDict()

        ## EVOLVED

        # Bulge & disk?
        if self.bulge_simulations.has_intrinsic_cube_edgeon and self.disk_simulations.has_intrinsic_cube_edgeon:

            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_edgeon
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_edgeon

        # Old?
        elif self.old_simulations.has_intrinsic_cube_edgeon: cubes[evolved_component_name] = self.old_intrinsic_stellar_luminosity_cube_edgeon

        # Not enough data
        else:
            warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the total simulation, the transparent edgeon cube will not be available")
            return None

        ## UNEVOLVED

        # Young & ionizing?
        if self.young_simulations.has_intrinsic_cube_edgeon and self.sfr_simulations.has_intrinsic_cube_edgeon:

            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_edgeon
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_edgeon

        # Unevolved?
        elif self.unevolved_simulations.has_intrinsic_cube_edgeon: cubes[unevolved_component_name] = self.unevolved_intrinsic_stellar_luminosity_cube_edgeon

        # Not enough data
        else:
            warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the total simulation, the transparent edgeon cube will not be available")
            return None

        # Return
        return cubes

    # -----------------------------------------------------------------
    # TOTAL SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulations(self):

        """
        This function ...
        :return:
        """

        # Load and return
        return MultiComponentSimulations.from_output_path(total_simulation_name, self.observed_total_output_path,
                                                          intrinsic_sed_paths=self.total_simulation_component_sed_paths,
                                                          distance=self.distance, intrinsic_cubes=self.total_simulation_component_cubes,
                                                          intrinsic_cubes_faceon=self.total_simulation_component_cubes_faceon,
                                                          intrinsic_cubes_edgeon=self.total_simulation_component_cubes_edgeon,
                                                          earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_simulation(self):
        #return self.total_simulations.observed
        return ObservedComponentSimulation.from_output_path(self.observed_total_output_path, total_simulation_name, earth_wcs=self.earth_wcs)

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
    def bulge_simulations(self):

        """
        This function ...
        :return:
        """

        # To run the simulation
        sed = self.old_bulge_component_sed

        # Load and return
        return SingleComponentSimulations.from_output_paths(bulge_simulation_name, observed=self.observed_bulge_output_path,
                                                            intrinsic=sed.out_path, distance=self.distance,
                                                            map_earth=self.old_bulge_map_earth, map_faceon=self.old_bulge_map_faceon,
                                                            map_edgeon=self.old_bulge_map_edgeon, earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_simulation(self):
        #return self.bulge_simulations.observed
        return ObservedComponentSimulation.from_output_path(self.observed_bulge_output_path, bulge_simulation_name, earth_wcs=self.earth_wcs)

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
    def disk_simulations(self):

        """
        This function ...
        :return:
        """

        # To run the simulation
        sed = self.old_disk_component_sed

        # Load and return
        return SingleComponentSimulations.from_output_paths(disk_simulation_name, observed=self.observed_disk_output_path,
                                                            intrinsic=sed.out_path, distance=self.distance,
                                                            map_earth=self.old_disk_map_earth, map_faceon=self.old_disk_map_faceon,
                                                            map_edgeon=self.old_disk_map_edgeon, earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @property
    def disk_simulation(self):
        #return self.disk_simulations.observed
        return ObservedComponentSimulation.from_output_path(self.observed_disk_output_path, disk_simulation_name, earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @property
    def disk_simulation_output(self):
        return self.disk_simulation.output

    # -----------------------------------------------------------------

    @property
    def disk_simulation_data(self):
        return self.disk_simulation.data

    # -----------------------------------------------------------------
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

    @lazyproperty
    def old_simulation_component_cubes(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the cubes
        cubes = OrderedDict()

        ## EVOLVED

        # Bulge & disk?
        if self.bulge_simulations.has_intrinsic_cube and self.disk_simulations.has_intrinsic_cube:

            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_earth
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_earth

        # Add
        #cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_earth
        #cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_earth

        # Return
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_component_cubes_faceon(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        cubes = OrderedDict()

        ## EVOLVED

        # Bulge & disk?
        if self.bulge_simulations.has_intrinsic_cube_faceon and self.disk_simulations.has_intrinsic_cube_faceon:

            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_faceon
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_faceon

        # Old?
        # elif self.old_simulations.has_intrinsic_cube_edgeon: cubes[evolved_component_name] = self.old_intrinsic_stellar_luminosity_cube_edgeon

        # Not enough data
        else:
            warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the old simulation, the transparent faceon cube will not be available")
            return None

        # Initialize dictionary for the cubes
        #cubes = OrderedDict()

        # Add
        #cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_faceon
        #cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_faceon

        # Return
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation_component_cubes_edgeon(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        cubes = OrderedDict()

        ## EVOLVED

        # Bulge & disk?
        if self.bulge_simulations.has_intrinsic_cube_edgeon and self.disk_simulations.has_intrinsic_cube_edgeon:

            cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_edgeon
            cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_edgeon

        # Old?
        #elif self.old_simulations.has_intrinsic_cube_edgeon: cubes[evolved_component_name] = self.old_intrinsic_stellar_luminosity_cube_edgeon

        # Not enough data
        else:
            warnings.warn("Not enough simulation data from the evolved intrinsic components. If no full cubes are available for the old simulation, the transparent edgeon cube will not be available")
            return None

        # Initialize dictionary for the cubes
        #cubes = OrderedDict()

        # Add
        #cubes[bulge_component_name] = self.bulge_intrinsic_stellar_luminosity_cube_edgeon
        #cubes[disk_component_name] = self.disk_intrinsic_stellar_luminosity_cube_edgeon

        # Return
        return cubes

    # -----------------------------------------------------------------
    # OLD SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulations(self):

        """
        This function ...
        :return:
        """

        # Load and return
        return MultiComponentSimulations.from_output_path(old_simulation_name, self.observed_old_output_path,
                                                          intrinsic_sed_paths=self.old_simulation_component_sed_paths,
                                                          distance=self.distance, intrinsic_cubes=self.old_simulation_component_cubes,
                                                          intrinsic_cubes_faceon=self.old_simulation_component_cubes_faceon,
                                                          intrinsic_cubes_edgeon=self.old_simulation_component_cubes_edgeon,
                                                          earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_simulation(self):
        #return self.old_simulations.observed
        return ObservedComponentSimulation.from_output_path(self.observed_old_output_path, old_simulation_name, earth_wcs=self.earth_wcs)

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
    def young_simulations(self):

        """
        This function ...
        :return:
        """

        # To run the simulation
        sed = self.young_component_sed

        # Load and return
        return SingleComponentSimulations.from_output_paths(young_simulation_name, observed=self.observed_young_output_path,
                                                            intrinsic=sed.out_path, distance=self.distance,
                                                            map_earth=self.young_map_earth, map_faceon=self.young_map_faceon,
                                                            map_edgeon=self.young_map_edgeon, earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_simulation(self):
        #return self.young_simulations.observed
        return ObservedComponentSimulation.from_output_path(self.observed_young_output_path, young_simulation_name, earth_wcs=self.earth_wcs)

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
    def sfr_simulations(self):

        """
        This function ...
        :return:
        """

        # To run the simulation
        sed = self.sfr_component_sed

        # Load and return
        return SingleComponentSimulations.from_output_paths(sfr_simulation_name, observed=self.observed_sfr_output_path,
                                                            intrinsic=sed.out_path, distance=self.distance,
                                                            map_earth=self.sfr_map_earth, map_faceon=self.sfr_map_faceon,
                                                            map_edgeon=self.sfr_map_edgeon, earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_simulation(self):
        #return self.sfr_simulations.observed
        return ObservedComponentSimulation.from_output_path(self.observed_sfr_output_path, sfr_simulation_name, earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @property
    def sfr_simulation_output(self):
        return self.sfr_simulation.output

    # -----------------------------------------------------------------

    @property
    def sfr_simulation_data(self):
        return self.sfr_simulation.data

    # -----------------------------------------------------------------
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

    @lazyproperty
    def unevolved_simulation_component_cubes(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        cubes = OrderedDict()

        # Add
        #cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_earth
        #cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_earth

        ## UNEVOLVED

        # Young & ionizing?
        if self.young_simulations.has_intrinsic_cube and self.sfr_simulations.has_intrinsic_cube:

            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_earth
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_earth

        # Not enough data
        else:
            warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the unevolved simulation, the transparent cube will not be available")
            return None

        # Return
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_cubes_faceon(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        cubes = OrderedDict()

        # Add
        #cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_faceon
        #cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_faceon

        ## UNEVOLVED

        # Young & ionizing?
        if self.young_simulations.has_intrinsic_cube_faceon and self.sfr_simulations.has_intrinsic_cube_faceon:

            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_faceon
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_faceon

        # Not enough data
        else:
            warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the unevolved simulation, the transparent faceon cube will not be available")
            return None

        # Return
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation_component_cubes_edgeon(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        cubes = OrderedDict()

        # Add
        #cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_edgeon
        #cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_edgeon

        ## UNEVOLVED

        # Young & ionizing?
        if self.young_simulations.has_intrinsic_cube_edgeon and self.sfr_simulations.has_intrinsic_cube_edgeon:

            cubes[young_component_name] = self.young_intrinsic_stellar_luminosity_cube_edgeon
            cubes[ionizing_component_name] = self.sfr_intrinsic_stellar_luminosity_cube_edgeon

        # Not enough data
        else:
            warnings.warn("Not enough simulation data from the unevolved intrinsic components. If no full cubes are available for the unevolved simulation, the transparent edgeon cube will not be available")
            return None

        # Return
        return cubes

    # -----------------------------------------------------------------
    # UNEVOLVED SIMULATIONS
    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulations(self):

        """
        This function ...
        :return:
        """

        # Load and return
        return MultiComponentSimulations.from_output_path(unevolved_simulation_name, self.observed_unevolved_output_path,
                                                          intrinsic_sed_paths=self.unevolved_simulation_component_sed_paths,
                                                          distance=self.distance, intrinsic_cubes=self.unevolved_simulation_component_cubes,
                                                          intrinsic_cubes_faceon=self.unevolved_simulation_component_cubes_faceon,
                                                          intrinsic_cubes_edgeon=self.unevolved_simulation_component_cubes_edgeon,
                                                          earth_wcs=self.earth_wcs)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_simulation(self):
        #return self.unevolved_simulations.observed
        return ObservedComponentSimulation.from_output_path(self.observed_unevolved_output_path, unevolved_simulation_name, earth_wcs=self.earth_wcs)

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
        return self.has_mappings

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings(self):
        # Create the MAPPINGS template and return it
        return Mappings(self.metallicity, self.sfr_compactness, self.sfr_pressure, self.sfr_covering_factor, self.sfr)

    # -----------------------------------------------------------------

    @property
    def has_mappings(self):
        return True # should always be able to be created

    # -----------------------------------------------------------------

    @lazyproperty
    def normalized_mappings(self):
        # Create the MAPPINGS template
        return Mappings(self.metallicity, self.sfr_compactness, self.sfr_pressure, self.sfr_covering_factor)

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
        return self.total_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_total_dust_sed(self):
        return self.total_simulations.has_observed_dust_sed

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
        return self.old_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_dust_sed(self):
        return self.old_simulations.has_observed_dust_sed

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
        return self.young_simulations.observed_bolometric_luminosity

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
        return self.young_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_young_dust_sed(self):
        return self.young_simulations.has_observed_dust_sed

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
        return self.sfr_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_dust_sed(self):
        return self.sfr_simulations.has_observed_dust_sed

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
        return self.unevolved_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_dust_sed(self):
        return self.unevolved_simulations.has_observed_dust_sed

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
    def sfr_sed_path(self):
        return fs.create_directory_in(self.sfr_path, sed_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_projections_path(self):
        return fs.create_directory_in(self.sfr_path, projections_dirname)

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
        return fs.join(self.sfr_sed_path, sfr_name + ".ski")

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
    def total_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.total_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.total_simulations.has_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def total_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_total_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.total_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def total_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.total_simulations.edgeon_observed_cube_absorbed

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
    def old_bolometric_luminosity_cube_faceon(self):
        return self.old_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def old_bolometric_luminosity_cube_edgeon(self):
        return self.old_simulations.edgeon_oberved_cube

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
    def old_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.old_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_old_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.old_simulations.has_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def old_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.old_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_old_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.old_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def old_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.old_simulations.edgeon_observed_cube_absorbed

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
    def young_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.young_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_young_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.young_simulations.has_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def young_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.young_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_young_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.young_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def young_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.young_simulations.edgeon_observed_cube_absorbed

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
    def sfr_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.sfr_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_sfr_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.sfr_simulations.has_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def sfr_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.sfr_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_sfr_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.sfr_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def sfr_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.sfr_simulations.edgeon_observed_cube_absorbed

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
    def unevolved_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.unevolved_simulations.observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorbed_diffuse_stellar_luminosity_cube_earth(self):
        return self.unevolved_simulations.has_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def unevolved_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.unevolved_simulations.faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorbed_diffuse_stellar_luminosity_cube_faceon(self):
        return self.unevolved_simulations.has_faceon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def unevolved_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.edgeon_observed_cube_absorbed

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorbed_diffuse_stellar_luminosity_cube_edgeon(self):
        return self.unevolved_simulations.has_edgeon_observed_cube_absorbed

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
        return salim_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_earth)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_earth_salim(self):
        return self.has_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_faceon_salim(self):
        return salim_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_faceon)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_faceon_salim(self):
        return self.has_intrinsic_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_edgeon_salim(self):
        return salim_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_edgeon)

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
        return kennicutt_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_earth)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_earth_kennicutt(self):
        return self.has_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_faceon_kennicutt(self):
        return kennicutt_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_faceon)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_faceon_kennicutt(self):
        return self.has_intrinsic_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_edgeon_kennicutt(self):
        return kennicutt_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_edgeon)

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
        return kennicutt_evans_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_earth)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_earth_ke(self):
        return self.has_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_faceon_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_faceon)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_faceon_ke(self):
        return self.has_intrinsic_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def total_star_formation_rate_map_edgeon_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_edgeon)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate_map_edgeon_ke(self):
        return self.has_intrinsic_fuv_luminosity_map_edgeon

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
    # 15. STELLAR MASS
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
        return oliver_stellar_mass(self.observed_i1_luminosity_map_earth, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype)

    # -----------------------------------------------------------------

    @property
    def has_total_stellar_mass_map_earth(self):
        return self.has_observed_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_stellar_mass_map_faceon(self):
        return oliver_stellar_mass(self.observed_i1_luminosity_map_faceon, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype)

    # -----------------------------------------------------------------

    @property
    def has_total_stellar_mass_map_faceon(self):
        return self.has_observed_i1_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_stellar_mass_map_edgeon(self):
        return oliver_stellar_mass(self.observed_i1_luminosity_map_edgeon, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype)

    # -----------------------------------------------------------------

    @property
    def has_total_stellar_mass_map_edgeon(self):
        return self.has_observed_i1_luminosity_map_edgeon

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

        """
        This function ...
        :return:
        """

        # Get the old bulge map
        frame = self.old_bulge_map_earth.copy()

        # Normalize to the I1 specific luminosity
        frame.normalize(to=self.intrinsic_i1_luminosity_old_bulge)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_intrinsic_i1_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the old bulge map
        frame = self.old_bulge_map_faceon.copy()

        # Normalize to the I1 specific luminosity
        frame.normalize(to=self.intrinsic_i1_luminosity_old_bulge)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_intrinsic_i1_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the old bulge map
        frame = self.old_bulge_map_edgeon.copy()

        # Normalize to the I1 specific luminosity
        frame.normalize(to=self.intrinsic_i1_luminosity_old_bulge)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_intrinsic_i1_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.has_old_bulge_intrinsic_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_intrinsic_i1_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_i1_luminosity_old_bulge and self.has_old_bulge_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_intrinsic_i1_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_i1_luminosity_old_bulge and self.has_old_bulge_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_intrinsic_i1_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_i1_luminosity_old_bulge and self.has_old_bulge_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_bulge_bolometric_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_bolometric_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        # Get the old bulge map
        frame = self.old_bulge_map_earth.copy()

        # Normalize to the I1 specific luminosity
        frame.normalize(to=self.intrinsic_i1_luminosity_old_bulge)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_bolometric_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the old bulge map
        frame = self.old_bulge_map_faceon.copy()

        # Normalize to the I1 specific luminosity
        frame.normalize(to=self.intrinsic_i1_luminosity_old_bulge)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_bolometric_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the old bulge map
        frame = self.old_bulge_map_edgeon.copy()

        # Normalize to the I1 specific luminosity
        frame.normalize(to=self.intrinsic_i1_luminosity_old_bulge)

        # Return the frame
        return frame

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

        """
        This function ...
        :return:
        """

        # Get the old disk map
        frame = self.old_disk_map_earth.copy()

        # Normalize to the I1 specific luminosity
        frame.normalize(to=self.intrinsic_i1_luminosity_old_disk)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_intrinsic_i1_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the old disk map
        frame = self.old_disk_map_faceon.copy()

        # Normalize to the I1 specific luminosity
        frame.normalize(to=self.intrinsic_i1_luminosity_old_disk)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_intrinsic_i1_luminosity_map_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        # Get the old disk map
        frame = self.old_disk_map_edgeon.copy()

        # Normalize to the I1 specific luminosity
        frame.normalize(to=self.intrinsic_i1_luminosity_old_disk)

        # Return the frame
        return frame

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

        """
        This function ...
        :return:
        """

        # Get the old disk map
        frame = self.old_disk_map_earth.copy()

        # Normalize to the bolometric luminosity
        frame.normalize(to=self.intrinsic_bolometric_luminosity_old_disk)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @property
    def has_old_disk_bolometric_luminosity_map_earth(self):
        return self.has_intrinsic_bolometric_luminosity_old_disk and self.has_old_disk_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_bolometric_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the old disk map
        frame = self.old_disk_map_faceon.copy()

        # Normalize to the bolometric luminosity
        frame.normalize(to=self.intrinsic_bolometric_luminosity_old_disk)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @property
    def has_old_disk_bolometric_luminosity_map_faceon(self):
        return self.has_intrinsic_bolometric_luminosity_old_disk and self.has_old_disk_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_bolometric_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the old disk map
        frame = self.old_disk_map_edgeon.copy()

        # Normalize to the bolometric luminosity
        frame.normalize(to=self.intrinsic_bolometric_luminosity_old_disk)

        # Return the frame
        return frame

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

        """
        This function ...
        :return:
        """

        # Get the young stellar map
        frame = self.young_map_earth.copy()

        # Normalize to the FUV luminosity
        frame.normalize(to=self.intrinsic_fuv_luminosity_young)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_fuv_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the young stellar map
        frame = self.young_map_faceon.copy()

        # Normalize to the FUV luminosity
        frame.normalize(to=self.intrinsic_fuv_luminosity_young)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_fuv_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the young stellar map
        frame = self.young_map_edgeon.copy()

        # Normalize to the FUV luminosity
        frame.normalize(to=self.intrinsic_fuv_luminosity_young)

        # Return the frame
        return frame

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

        """
        This function ...
        :return:
        """

        # Get the young stellar map
        frame = self.young_map_earth.copy()

        # Normalize to the bolometric luminosity
        frame.normalize(to=self.intrinsic_bolometric_luminosity_young)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def young_bolometric_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the young stellar map
        frame = self.young_map_faceon.copy()

        # Normalize to the bolometric luminosity
        frame.normalize(to=self.intrinsic_bolometric_luminosity_young)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def young_bolometric_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the young stellar map
        frame = self.young_map_edgeon.copy()

        # Normalize to the bolometric luminosity
        frame.normalize(to=self.intrinsic_bolometric_luminosity_young)

        # Return the frame
        return frame

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
        return salim_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_earth)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_earth_salim(self):
        return self.has_young_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_map_faceon_salim(self):
        return salim_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_faceon)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_faceon_salim(self):
        return self.has_young_intrinsic_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_map_edgeon_salim(self):
        return salim_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_edgeon)

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
        return kennicutt_evans_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_earth)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_earth_ke(self):
        return self.has_young_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------
    #     FACEON
    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_map_faceon_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_faceon)

    # -----------------------------------------------------------------

    @property
    def has_young_star_formation_rate_map_faceon_ke(self):
        return fs.is_file(self.young_star_formation_rate_map_faceon_ke)

    # -----------------------------------------------------------------
    #     EDGEON
    # -----------------------------------------------------------------

    @lazyproperty
    def young_star_formation_rate_map_edgeon_ke(self):
        return kennicutt_evans_fuv_to_sfr(self.young_intrinsic_fuv_luminosity_map_edgeon)

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

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_earth.copy()

        # Normalize to the FUV luminosity
        frame.normalize(to=self.intrinsic_fuv_luminosity_sfr)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_fuv_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_faceon.copy()

        # Normalize to the FUV luminosity
        frame.normalize(to=self.intrinsic_fuv_luminosity_sfr)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsinc_fuv_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_edgeon.copy()

        # Normalize to the FUV luminosity
        frame.normalize(to=self.intrinsic_fuv_luminosity_sfr)

        # Return the frame
        return frame

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

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_earth.copy()

        # Normalize to the bolometric luminosity
        frame.normalize(to=self.intrinsic_bolometric_luminosity_sfr)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_bolometric_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_faceon.copy()

        # Normalize to the bolometric luminosity
        frame.normalize(to=self.intrinsic_bolometric_luminosity_sfr)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_bolometric_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_edgeon.copy()

        # Normalize to the bolometric luminosity
        frame.normalize(to=self.intrinsic_bolometric_luminosity_sfr)

        # Return the frame
        return frame

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

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_earth.copy()

        # Normalize to the star formation rate
        frame.normalize(to=self.sfr)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def star_formation_rate_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_faceon.copy()

        # Normalize to the star formation rate
        frame.normalize(to=self.sfr)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def star_formation_rate_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_edgeon.copy()

        # Normalize to the star formation rate
        frame.normalize(to=self.sfr)

        # Return the frame
        return frame

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

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_earth.copy()

        # Normalize to the SF dust mass
        frame.normalize(to=self.sfr_dust_mass)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_faceon.copy()

        # Normalize to the SF dust mass
        frame.normalize(to=self.sfr_dust_mass)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_edgeon.copy()

        # Normalize to the SF dust mass
        frame.normalize(to=self.sfr_dust_mass)

        # Return the frame
        return frame

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

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_earth.copy()

        # Normalize
        frame.normalize(to=self.intrinsic_stellar_luminosity_sfr)

        # Return
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_stellar_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_faceon.copy()

        # Normalize
        frame.normalize(to=self.intrinsic_stellar_luminosity_sfr)

        # Return
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_stellar_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_edgeon.copy()

        # Normalize
        frame.normalize(to=self.intrinsic_stellar_luminosity_sfr)

        # Return
        return frame

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

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_earth.copy()

        # Normalize
        frame.normalize(to=self.intrinsic_dust_luminosity_sfr)

        # Return
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_dust_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_faceon.copy()

        # Normalize
        frame.normalize(to=self.intrinsic_dust_luminosity_sfr)

        # Return
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_dust_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the SF map
        frame = self.sfr_map_edgeon.copy()

        # Normalize
        frame.normalize(to=self.intrinsic_dust_luminosity_sfr)

        # Return
        return frame

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
        return salim_fuv_to_sfr(self.unevolved_intrinsic_fuv_luminosity_map_earth)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_star_formation_rate_map_earth(self):
        return self.has_unevolved_intrinsic_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_star_formation_rate_map_faceon(self):
        return salim_fuv_to_sfr(self.unevolved_intrinsic_fuv_luminosity_map_faceon)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_star_formation_rate_map_faceon(self):
        return self.has_unevolved_intrinsic_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_star_formation_rate_map_edgeon(self):
        return salim_fuv_to_sfr(self.unevolved_intrinsic_fuv_luminosity_map_edgeon)

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

        """
        This function ...
        :return:
        """

        return self.definition.dust_map_path

    # -----------------------------------------------------------------

    @property
    def has_dust_map(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.is_file(self.dust_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.dust_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map_wcs(self):

        """
        This function ...
        :return:
        """

        return CoordinateSystem.from_file(self.dust_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map_projection(self):

        """
        This function ...
        :return:
        """

        azimuth = 0.0
        if not self.has_center: raise ValueError("Galaxy center coordinate is not defined")
        return GalaxyProjection.from_wcs(self.dust_map_wcs, self.center, self.distance, self.inclination, azimuth, self.position_angle)

    # -----------------------------------------------------------------

    @property
    def dust_map_shape(self):

        """
        This function ...
        :return:
        """

        return self.dust_map_wcs.shape

    # -----------------------------------------------------------------

    @property
    def dust_map_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.dust_map_wcs.pixelscale

    # -----------------------------------------------------------------

    @property
    def diffuse_dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return self.diffuse_dust_mass_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_mass_map_earth(self):

        """
        This function ...
        :return:
        """

        # Get the dust map
        frame = self.dust_map_earth.copy()

        # Normalize to the dust mass
        frame.normalize(to=self.diffuse_dust_mass)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the dust map
        frame = self.dust_map_faceon.copy()

        # Normalize to the dust mass
        frame.normalize(to=self.diffuse_dust_mass)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the dust map
        frame = self.dust_map_edgeon.copy()

        # Normalize to the dust mass
        frame.normalize(to=self.diffuse_dust_mass)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_mass(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def has_dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.has_diffuse_dust_mass and self.has_sfr_dust_mass

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return self.has_diffuse_dust_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_mass_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_map_earth and self.has_diffuse_dust_mass

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_map_edgeon and self.has_diffuse_dust_mass

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_map_faceon and self.has_diffuse_dust_mass

    # -----------------------------------------------------------------

    @property
    def dust_mass_map(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.has_dust_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_dust_mass_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_diffuse_dust_mass_map_earth and self.has_sfr_dust_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_diffuse_dust_mass_map_faceon and self.has_sfr_dust_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_diffuse_dust_mass_map_edgeon and self.has_sfr_dust_mass_map_edgeon

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
    def old_bulge_projections(self):

        """
        This function ...
        :return:
        """

        return ComponentProjections(bulge_simulation_name, self.old_bulge_model, path=self.old_bulge_projections_path,
                                    description="old bulge stellar component",
                                    projection=self.old_disk_projections.projection_earth,
                                    projection_faceon=self.old_disk_projections.projection_faceon,
                                    projection_edgeon=self.old_disk_projections.projection_edgeon, center=self.center,
                                    earth_wcs=self.old_disk_map_wcs, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def old_bulge_map(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_map_earth

    # -----------------------------------------------------------------

    @property
    def old_bulge_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_projections.earth

    # -----------------------------------------------------------------

    @property
    def old_bulge_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_projections.faceon

    # -----------------------------------------------------------------

    @property
    def old_bulge_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_projections.edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_map(self):

        """
        This function ...
        :return:
        """

        return self.has_old_bulge_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_map_earth(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_map_faceon(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_projections(self):

        """
        This function ...
        :return:
        """

        return ComponentProjections(disk_simulation_name, self.old_disk_deprojection, path=self.old_disk_projections_path,
                                    earth=False, description="old disk stellar component",
                                    input_filepaths=[self.old_disk_map_path], center=self.center,
                                    earth_wcs=self.old_disk_map_wcs, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def old_disk_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_map # just the input map

    # -----------------------------------------------------------------

    @property
    def old_disk_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_projections.faceon

    # -----------------------------------------------------------------

    @property
    def old_disk_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_projections.edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_disk_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_old_disk_map

    # -----------------------------------------------------------------

    @property
    def has_old_disk_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_old_disk_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------

    @property
    def has_old_disk_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_old_disk_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------

    @lazyproperty
    def young_projections(self):

        """
        This function ...
        :return:
        """

        return ComponentProjections(young_simulation_name, self.young_deprojection, path=self.young_projections_path,
                                    earth=False, description="young stellar component",
                                    input_filepaths=[self.young_map_path], center=self.center,
                                    earth_wcs=self.young_map_wcs, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def young_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.young_map # just the input map

    # -----------------------------------------------------------------

    @property
    def young_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.young_projections.faceon

    # -----------------------------------------------------------------

    @property
    def young_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.young_projections.edgeon

    # -----------------------------------------------------------------

    @property
    def has_young_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_young_map

    # -----------------------------------------------------------------

    @property
    def has_young_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_young_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------

    @property
    def has_young_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_young_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_projections(self):

        """
        This function ...
        :return:
        """

        return ComponentProjections(sfr_simulation_name, self.sfr_deprojection, path=self.sfr_projections_path,
                                    earth=False, description="SFR component", input_filepaths=[self.sfr_map_path],
                                    center=self.center, earth_wcs=self.sfr_map_wcs, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def sfr_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.sfr_map # just the input map

    # -----------------------------------------------------------------

    @property
    def sfr_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.sfr_projections.faceon

    # -----------------------------------------------------------------

    @property
    def sfr_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.sfr_projections.edgeon

    # -----------------------------------------------------------------

    @property
    def has_sfr_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr_map

    # -----------------------------------------------------------------

    @property
    def has_sfr_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------

    @property
    def has_sfr_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr_map # if we have the input map, we can deproject it

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_projections(self):

        """
        This function ...
        :return:
        """

        return ComponentProjections(dust_simulation_name, self.dust_deprojection, path=self.dust_projections_path,
                                    earth=False, description="dust component", input_filepaths=[self.dust_map_path],
                                    distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def dust_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.dust_map # just the input map

    # -----------------------------------------------------------------

    @property
    def dust_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.dust_projections.faceon

    # -----------------------------------------------------------------

    @property
    def dust_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.dust_projections.edgeon

    # -----------------------------------------------------------------

    @property
    def has_dust_map_earth(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_dust_map

    # -----------------------------------------------------------------

    @property
    def has_dust_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_map # if there is an input map, we can deproject it

    # -----------------------------------------------------------------

    @property
    def has_dust_map_edgeon(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.old_bulge_component_sed.sed

    # -----------------------------------------------------------------

    @property
    def old_bulge_sed_filepath(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.old_disk_component_sed.sed

    # -----------------------------------------------------------------

    @property
    def old_disk_sed_filepath(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.young_component_sed.sed

    # -----------------------------------------------------------------

    @property
    def young_sed_filepath(self):

        """
        This function ...
        :return:
        """

        return self.young_component_sed.sed_filepath

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_component_sed(self):

        """
        This function ...
        :return:
        """

        return ComponentSED(ionizing_component_name, self.sfr_component, description="SFR component",
                            path=self.sfr_sed_path, input_filepaths=self.sfr_input_filepaths,
                            distance=self.distance, inclination=self.inclination, position_angle=self.position_angle,
                            wavelengths_filename=wavelengths_filename)

    # -----------------------------------------------------------------

    @property
    def sfr_sed(self):

        """
        This function ...
        :return:
        """

        return self.sfr_component_sed.sed

    # -----------------------------------------------------------------

    @property
    def sfr_sed_filepath(self):

        """
        This function ...
        :return:
        """

        return self.sfr_component_sed.sed_filepath

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_dust_sed_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.intrinsic_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_dust_sed_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_intrinsic_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.dust_sed - self.intrinsic_dust_sed_sfr

    # -----------------------------------------------------------------

    @property
    def has_diffuse_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_sed and self.has_intrinsic_dust_sed_sfr

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

    @property
    def absorbed_diffuse_stellar_luminosity(self):
        return self.total_absorbed_diffuse_stellar_luminosity_map.sum(add_unit=True, per_area="error")

    # -----------------------------------------------------------------

    @property
    def has_absorbed_diffuse_stellar_luminosity(self):
        return self.has_total_absorbed_diffuse_stellar_luminosity_map

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

    @property
    def direct_stellar_luminosity(self):
        return self.total_direct_stellar_luminosity_map.sum()

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity(self):
        return self.has_total_direct_stellar_luminosity_map

    # -----------------------------------------------------------------

    @property
    def intrinsic_fuv_luminosity(self):
        return self.intrinsic_fuv_luminosity_map.sum(add_unit=True, per_area="error")

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity(self):
        return self.has_intrinsic_fuv_luminosity_map

    # -----------------------------------------------------------------

    @property
    def total_star_formation_rate(self):
        return salim_fuv_to_sfr(self.intrinsic_fuv_luminosity)

    # -----------------------------------------------------------------

    @property
    def has_total_star_formation_rate(self):
        return self.has_intrinsic_fuv_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_i1_luminosity(self):
        return self.observed_total_sed.photometry_at(self.i1_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity(self):
        return True

    # -----------------------------------------------------------------

    @property
    def total_stellar_mass(self):
        return oliver_stellar_mass(self.observed_i1_luminosity, hubble_type=self.hubble_type, hubble_subtype=self.hubble_subtype)

    # -----------------------------------------------------------------

    @property
    def has_total_stellar_mass(self):
        return True

    # -----------------------------------------------------------------

    @property
    def total_ssfr(self):
        return self.total_ssfr_map.average(add_unit=True)

    # -----------------------------------------------------------------

    @property
    def has_total_ssfr(self):
        return self.has_total_ssfr_map

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
        if self.has_total_star_formation_rate: values[sfr_name] = self.total_star_formation_rate

        # Stellar mass
        if self.has_total_stellar_mass: values[stellar_mass_name] = self.total_stellar_mass

        # Specific star formation rate
        if self.has_total_ssfr: values[ssfr_name] = self.total_ssfr

        # Return
        return values

    # -----------------------------------------------------------------

    @property
    def direct_stellar_luminosity_old_bulge(self):
        return self.old_bulge_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_old_bulge(self):
        return self.has_old_bulge_direct_stellar_luminosity_map

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

        # Attenuation
        if self.has_i1_attenuation_old_bulge: values[bulge_spec_attenuation_name] = self.i1_attenuation_old_bulge
        if self.has_bolometric_attenuation_old_bulge: values[bulge_bol_attenuation_name] = self.bolometric_attenuation_old_bulge

        # Dust
        if self.has_observed_dust_luminosity_old_bulge: values[obs_bulge_dust_lum_name] = self.observed_dust_luminosity_old_bulge

        # Return
        return values

    # -----------------------------------------------------------------

    @property
    def direct_stellar_luminosity_old_disk(self):
        return self.old_disk_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_old_disk(self):
        return self.has_old_disk_direct_stellar_luminosity_map

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

        # Attenuation
        if self.has_i1_attenuation_old_disk: values[disk_spec_attenuation_name] = self.i1_attenuation_old_disk
        if self.has_bolometric_attenuation_old_disk: values[disk_bol_attenuation_name] = self.bolometric_attenuation_old_disk

        # Dust
        if self.has_observed_dust_luminosity_old_disk: values[obs_disk_dust_lum_name] = self.observed_dust_luminosity_old_disk

        # Return
        return values

    # -----------------------------------------------------------------

    @property
    def direct_stellar_luminosity_old(self):
        return self.old_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_old(self):
        return self.has_old_direct_stellar_luminosity_map

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

        # Attenuation
        if self.has_i1_attenuation_old: values[old_spec_attenuation_name] = self.i1_attenuation_old
        if self.has_bolometric_attenuation_old: values[old_bol_attenuation_name] = self.bolometric_attenuation_old

        # Dust
        if self.has_observed_dust_luminosity_old: values[obs_old_dust_lum_name] = self.observed_dust_luminosity_old

        # Return
        return values

    # -----------------------------------------------------------------

    @property
    def direct_stellar_luminosity_young(self):
        return self.young_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_young(self):
        return self.has_young_direct_stellar_luminosity_map

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

        # Attenuation
        if self.has_fuv_attenuation_young: values[young_spec_attenuation_name] = self.fuv_attenuation_young
        if self.has_bolometric_attenuation_young: values[young_bol_attenuation_name] = self.bolometric_attenuation_young

        # Dust
        if self.has_observed_dust_luminosity_young: values[obs_young_dust_lum_name] = self.observed_dust_luminosity_young

        # Return
        return values

    # -----------------------------------------------------------------

    @property
    def direct_stellar_luminosity_sfr(self):
        return self.sfr_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_sfr(self):
        return self.has_sfr_direct_stellar_luminosity_map

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
        if self.has_sfr: values[sfr_name] = self.sfr

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

    @property
    def direct_stellar_luminosity_unevolved(self):
        return self.unevolved_direct_stellar_luminosity_map.sum(add_unit=True, per_area="error")

    # -----------------------------------------------------------------

    @property
    def has_direct_stellar_luminosity_unevolved(self):
        return self.has_unevolved_direct_stellar_luminosity_map

    # -----------------------------------------------------------------

    @property
    def unevolved_intrinsic_fuv_luminosity(self):
        return self.intrinsic_fuv_luminosity_unevolved

    # -----------------------------------------------------------------

    @property
    def has_unevolved_intrinsic_fuv_luminosity(self):
        return self.has_intrinsic_fuv_luminosity_unevolved

    # -----------------------------------------------------------------

    @property
    def unevolved_star_formation_rate(self):
        return salim_fuv_to_sfr(self.unevolved_intrinsic_fuv_luminosity)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_star_formation_rate(self):
        return self.has_unevolved_intrinsic_fuv_luminosity

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
        if self.has_unevolved_star_formation_rate: values[sfr_name] = self.unevolved_star_formation_rate

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

kennicutt_evans_logc = 43.35 # 2012
kennicutt = 1.4e-28 # 1998
salim = 1.08e-28 # 200

# -----------------------------------------------------------------

def kennicutt_evans_fuv_to_sfr(fuv_luminosity, unit=None):

    """
    This function ...
    :param fuv_luminosity:
    :param unit:
    :return:
    """

    # Calculate factor
    calibration = 1./ 10**kennicutt_evans_logc

    # Get the FUV wavelength
    fuv_wavelength = parse_filter("GALEX FUV").wavelength

    from ..core.data import Data3D

    # Frame
    if isinstance(fuv_luminosity, Frame):

        if unit is not None: raise ValueError("Cannot specify unit")

        converted = fuv_luminosity.converted_to("erg/s", density=True, wavelength=fuv_wavelength)
        converted *= calibration
        converted.unit = "Msun/yr"
        return converted

    # 3D data
    elif isinstance(fuv_luminosity, Data3D):

        if unit is not None: raise ValueError("Cannot specify unit")

        factor = fuv_luminosity.unit.conversion_factor("erg/s", density=True, wavelength=fuv_wavelength)
        factor *= calibration
        return fuv_luminosity.converted_by_factor(factor, "Msun/yr", new_name="SFR", new_description="star formation rate (Kennicutt & Evans)")

    # Photometric quantity
    elif types.is_quantity(fuv_luminosity):
        if unit is not None: raise ValueError("Cannot specify unit")
        return fuv_luminosity.to("erg/s", density=True).value * calibration * u("Msun/yr")

    # Array
    elif types.is_array_like(fuv_luminosity):

        if unit is None: raise ValueError("Unit is not specified")
        factor = unit.conversion_factor("erg/s", density=True, wavelength=fuv_wavelength)
        factor *= calibration
        return fuv_luminosity * factor

    # Invalid
    else: raise ValueError("Invalid type for 'fuv_luminosity'")

# -----------------------------------------------------------------

def kennicutt_fuv_to_sfr(fuv_luminosity, unit=None):

    """
    This function ...
    :param fuv_luminosity:
    :param unit:
    :return:
    """

    # Get the FUV wavelength
    fuv_wavelength = parse_filter("GALEX FUV").wavelength

    from ..core.data import Data3D

    # Frame
    if isinstance(fuv_luminosity, Frame):

        if unit is not None: raise ValueError("Cannot specify unit")

        converted = fuv_luminosity.converted_to("erg/s/Hz", wavelength=fuv_wavelength)
        converted *= kennicutt
        converted.unit = "Msun/yr"
        return converted

    # 3D data
    elif isinstance(fuv_luminosity, Data3D):

        if unit is not None: raise ValueError("Cannot specify unit")

        factor = fuv_luminosity.unit.conversion_factor("erg/s/Hz", wavelength=fuv_wavelength)
        factor *= kennicutt
        return fuv_luminosity.converted_by_factor(factor, "Msun/yr", new_name="SFR", new_description="star formation rate (Kennicutt)")

    # Photometric quantity
    elif types.is_quantity(fuv_luminosity):
        if unit is not None: raise ValueError("Cannot specify unit")
        return fuv_luminosity.to("erg/s/Hz").value * kennicutt * u("Msun/yr")

    # Array
    elif types.is_array_like(fuv_luminosity):
        if unit is None: raise ValueError("Unit is not specified")
        factor = unit.conversion_factor("erg/s/Hz", wavelength=fuv_wavelength)
        factor *= kennicutt
        return fuv_luminosity * factor

    # Invalid
    else: raise ValueError("Invalid type for 'fuv_luminosity'")

# -----------------------------------------------------------------

def salim_fuv_to_sfr(fuv_luminosity, unit=None):

    """
    This function ...
    :param fuv_luminosity:
    :param unit:
    :return:
    """

    # Get the FUV wavelength
    fuv_wavelength = parse_filter("GALEX FUV").wavelength

    from ..core.data import Data3D

    # Frame
    if isinstance(fuv_luminosity, Frame):

        if unit is not None: raise ValueError("Cannot specify unit")

        converted = fuv_luminosity.converted_to("erg/s/Hz", wavelength=fuv_wavelength)
        converted *= salim
        converted.unit = "Msun/yr"
        return converted

    # 3D data
    elif isinstance(fuv_luminosity, Data3D):

        if unit is not None: raise ValueError("Cannot specify unit")

        factor = fuv_luminosity.unit.conversion_factor("erg/s/Hz", wavelength=fuv_wavelength)
        factor *= salim
        return fuv_luminosity.converted_by_factor(factor, "Msun/yr", new_name="SFR", new_description="star formation rate (Salim)")

    # Photometric quantity
    elif types.is_quantity(fuv_luminosity):
        if unit is not None: raise ValueError("Cannot specify unit")
        return fuv_luminosity.to("erg/s/Hz") * salim * u("Msun/yr")

    # Array
    elif types.is_array_like(fuv_luminosity):
        if unit is None: raise ValueError("Unit is not specified")
        factor = unit.conversion_factor("erg/s/Hz", wavelength=fuv_wavelength)
        factor *= salim
        return fuv_luminosity * factor

    # Invalid
    else: raise ValueError("Invalid type for 'fuv_luminosity'")

# -----------------------------------------------------------------

def hubble_stage_to_type(stage, add_subtype=False):

    """
    # SOURCE: https://en.wikipedia.org/wiki/Galaxy_morphological_classification
    This function ...
    :param stage:
    :param add_subtype:
    :return:
    """

    if stage < -3.5:

        if add_subtype:

            if stage < -5.5: subtype = "E-"
            elif stage < -4.5: subtype = "E"
            else: subtype = "E+"
            return "E", subtype

        else: return "E"

    elif stage < -0.5:

        if add_subtype:

            if stage < -2.5: subtype = "S0-"
            elif stage < -1.5: subtype = "S00"
            else: subtype = "S0+"
            return "S0", subtype

        else: return "S0"

    elif stage < 0.5:

        if add_subtype: return "S0/a", None
        else: return "S0/a"

    elif stage < 1.5:

        if add_subtype: return "Sa", None
        else: return "Sa"

    elif stage < 2.5:

        if add_subtype: return "Sab", None
        else: return "Sab"

    elif stage < 3.5:

        if add_subtype: return "Sb", None
        else: return "Sb"

    elif stage < 4.5:

        if add_subtype: return "Sbc", None
        else: return "Sbc"

    elif stage < 7.5:

        if add_subtype:
            if stage < 5.5: subtype = "Sc"
            elif stage < 6.5: subtype = "Scd"
            else: subtype = "Sd"
            return "Sc", subtype
        else: return "Sc"

    elif stage < 8.5:

        if add_subtype: return "Sc/Irr", None
        else: return "Sc/Irr"

    else:

        if add_subtype:
            if stage < 9.5: subtype = "Sm"
            elif stage < 10.5: subtype = "Im"
            else: subtype = None
            return "Irr", subtype
        else: return "Irr"

# -----------------------------------------------------------------

def hubble_type_to_stage(hubble_class):

    """
    This function ...
    # SOURCE: https://en.wikipedia.org/wiki/Galaxy_morphological_classification
    :param hubble_class:
    :return:
    """

    # E
    if hubble_class == "E": return -5

    # SUBDIVISION (VAUCOULEURS)
    elif hubble_class == "cE": return -6
    elif hubble_class == "E+": return -4

    # S0
    elif hubble_class == "S0": return -2

    # SUBDIVISION (VAUCOULEURS)
    elif hubble_class == "S0-": return -3
    elif hubble_class == "S00": return -2
    elif hubble_class == "S0+": return -1

    # S0/a
    elif hubble_class == "S0/a": return 0

    # Sa
    elif hubble_class == "Sa": return 1

    # Sab
    elif hubble_class == "Sab": return 2
    # synonym
    elif hubble_class == "Sa-b": return 2

    # Sb
    elif hubble_class == "Sb": return 3

    # Sbc
    elif hubble_class == "Sbc": return 4
    elif hubble_class == "Sb-c": return 4

    # Sc
    elif hubble_class == "Sc": return 5.5
    elif hubble_class == "Scd": return 6
    elif hubble_class == "Sc-d": return 6
    elif hubble_class == "Sd": return 7

    # Sc/Irr
    elif hubble_class == "Sc/Irr": return 8
    elif hubble_class == "Sc-Irr": return 8
    elif hubble_class == "Sdm": return 8

    # Higher
    elif hubble_class == "Sm": return 9
    elif hubble_class == "Irr": return 9.5
    elif hubble_class == "Im": return 10
    elif hubble_class == "I": return 9.5

    # Invalid
    else: raise ValueError("Unknown hubble class: '" + hubble_class + "'")

# -----------------------------------------------------------------

# OLIVER 2010
# (M∗/M⊙)/[νLν(3.6)/L⊙] to be 38.4, 40.8, 27.6, 35.3, 18.7 and 26.7, for types E, Sab, Sbc, Scd, Sdm and sb
# measuring the 3.6 μm monochromatic luminosity in total solar units, not in units of the Sun’s monochromatic 3.6 μm lu- minosity

oliver_stellar_mass_factors = OrderedDict()
oliver_stellar_mass_factors["E"] = 38.4
oliver_stellar_mass_factors["Sab"] = 40.8
oliver_stellar_mass_factors["Sb"] = 26.7
oliver_stellar_mass_factors["Sbc"] = 27.6
oliver_stellar_mass_factors["Scd"] = 35.3
oliver_stellar_mass_factors["Sdm"] = 18.7

# -----------------------------------------------------------------

def get_oliver_stellar_mass_factor(hubble_type, hubble_subtype=None):

    """
    Thisf unction ...
    :param hubble_type:
    :param hubble_subtype:
    :return:
    """

    # Get the factor
    if hubble_type not in oliver_stellar_mass_factors:
        if hubble_subtype is not None:
            if hubble_subtype not in oliver_stellar_mass_factors: raise ValueError("Hubble type '" + hubble_type + "' or '" + hubble_subtype + "' not supported")
            else: factor = oliver_stellar_mass_factors[hubble_subtype]
        else: raise ValueError("Hubble type '" + hubble_type + "' not supported")
    else: factor = oliver_stellar_mass_factors[hubble_type]

    # Return the factor
    return factor

# -----------------------------------------------------------------

def oliver_stellar_mass(i1_luminosity, hubble_type, hubble_subtype=None):

    """
    This function ...
    :param i1_luminosity:
    :param hubble_type:
    :param hubble_subtype:
    :return:
    """

    from ..core.data import Data3D

    # Get the I1 wavelength
    i1_wavelength = parse_filter("IRAC I1").wavelength

    # Get the factor
    oliver_factor = get_oliver_stellar_mass_factor(hubble_type, hubble_subtype=hubble_subtype)

    # Frame
    if isinstance(i1_luminosity, Frame):

        converted = i1_luminosity.converted_to("Lsun", density=True, wavelength=i1_wavelength)
        converted *= oliver_factor
        converted.unit = "Msun"
        return converted

    # 3D data
    elif isinstance(i1_luminosity, Data3D):

        factor = i1_luminosity.unit.conversion_factor("Lsun", density=True, wavelength=i1_wavelength)
        factor *= oliver_factor
        return i1_luminosity.converted_by_factor(factor, "Msun", new_name="Mstar", new_description="Stellar mass (Oliver)")

    # Photometric quantity
    else: return i1_luminosity.to("Lsun", density=True, wavelength=i1_wavelength).value * oliver_factor * u("Msun")

# -----------------------------------------------------------------
