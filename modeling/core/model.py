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
from ..basics.instruments import SEDInstrument
from ...core.tools import filesystem as fs
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ..simulation.single import SingleComponentSimulations
from ..simulation.multi import MultiComponentSimulations
from ...magic.core.frame import Frame
from ...magic.core.list import convolve_and_rebin
from ...magic.basics.coordinatesystem import CoordinateSystem
from ..basics.projection import GalaxyProjection
from ..basics.instruments import FrameInstrument, FullSEDInstrument
from ..simulation.projections import ComponentProjections
from ..simulation.sed import ComponentSED

# -----------------------------------------------------------------

default_scale_heights = 15. # number of times to take the old stellar scale height as the vertical radius of the model

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
                 parameters=None, center=None):

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

        # Set the center
        self.center = center

    # -----------------------------------------------------------------

    @property
    def has_center(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.old_bulge_sed_filepath

    # -----------------------------------------------------------------

    @property
    def intrinsic_sed_path_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_sed_filepath

    # -----------------------------------------------------------------

    @property
    def intrinsic_sed_path_young(self):

        """
        This function ...
        :return:
        """

        return self.young_sed_filepath

    # -----------------------------------------------------------------

    @property
    def intrinsic_sed_path_sfr(self):

        """
        This function ...
        :return:
        """

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
    def total_simulations(self):

        """
        This function ...
        :return:
        """

        # Load and return
        return MultiComponentSimulations.from_output_path(total_simulation_name, self.observed_total_output_path,
                                                          intrinsic_sed_paths=self.total_simulation_component_sed_paths,
                                                          distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_simulations(self):

        """
        This function ...
        :return:
        """

        # Run simulation?
        if not self.has_old_bulge_sed_simulation: self.run_old_bulge_sed_simulation()

        # Load and return
        return SingleComponentSimulations.from_output_paths(bulge_simulation_name, observed=self.observed_bulge_output_path,
                                                            intrinsic=self.old_bulge_sed_out_path, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_simulations(self):

        """
        This function ...
        :return:
        """

        # Run simulation?
        if not self.has_old_disk_sed_simulation: self.run_old_disk_sed_simulation()

        # Load and return
        return SingleComponentSimulations.from_output_paths(disk_simulation_name, observed=self.observed_disk_output_path,
                                                            intrinsic=self.old_disk_sed_out_path, distance=self.distance)

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
    def old_simulations(self):

        """
        This function ...
        :return:
        """

        # Load and return
        return MultiComponentSimulations.from_output_path(old_simulation_name, self.observed_old_output_path,
                                                          intrinsic_sed_paths=self.old_simulation_component_sed_paths,
                                                          distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_simulations(self):

        """
        This function ...
        :return:
        """

        # Run simulation?
        if not self.has_young_sed_simulation: self.run_young_sed_simulation()

        # Load and return
        return SingleComponentSimulations.from_output_paths(young_simulation_name, observed=self.observed_young_output_path,
                                                            distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_simulations(self):

        """
        This function ...
        :return:
        """

        # Run simulation?
        if not self.has_sfr_sed_simulation: self.run_sfr_sed_simulation()

        # Load and return
        return SingleComponentSimulations.from_output_paths(sfr_simulation_name, observed=self.observed_sfr_output_path,
                                                            distance=self.distance)

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
    def unevolved_simulations(self):

        """
        This function ...
        :return:
        """

        # Load and return
        return MultiComponentSimulations.from_output_path(unevolved_simulation_name, self.observed_unevolved_output_path,
                                                          intrinsic_sed_paths=self.unevolved_simulation_component_sed_paths,
                                                          distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_observed_total_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_total_output_path is not None and fs.is_directory(self.observed_total_output_path) and not fs.is_empty(self.observed_total_output_path)

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
    def has_observed_old_disk_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_disk_output_path is not None and fs.is_directory(self.observed_disk_output_path) and not fs.is_empty(self.observed_disk_output_path)

    # -----------------------------------------------------------------

    @property
    def has_observed_old_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_old_output_path is not None and fs.is_directory(self.observed_old_output_path) and not fs.is_empty(self.observed_old_output_path)

    # -----------------------------------------------------------------

    @property
    def has_observed_young_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_young_output_path is not None and fs.is_directory(self.observed_young_output_path) and not fs.is_empty(self.observed_young_output_path)

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_output(self):

        """
        This function ...
        :return:
        """

        return self.observed_sfr_output_path is not None and fs.is_directory(self.observed_sfr_output_path) and not fs.is_empty(self.observed_sfr_output_path)

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
    # -----------------------------------------------------------------

    # TOTAL

    # -----------------------------------------------------------------

    @property
    def has_observed_total_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_total_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_total_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_total_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.observed_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_observed_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        #return self.intrinsic_bolometric_luminosity_old + self.intrinsic_bolometric_luminosity_young + self.intrinsic_bolometric_luminosity_sfr
        return self.total_simulations.intrinsic_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        #return self.has_intrinsic_bolometric_luminosity_old and self.has_intrinsic_bolometric_luminosity_young and self.has_intrinsic_bolometric_luminosity_sfr
        return self.total_simulations.has_intrinsic_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_total_sed(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_total_sed(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_total_sed(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_total_sed(self):

        """
        This function ...
        :return:
        """

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
        if self.total_simulations.observed_stellar_sed_needs_reprocessed_internal_part: return sed + self.intrinsic_dust_sed_sfr
        else: return sed

    # -----------------------------------------------------------------

    @property
    def has_observed_total_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_total_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_total_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    # OLD BULGE

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the I1 wavelength
        return self.bulge_simulations.observed_photometry_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.has_observed_photometry

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

    @property
    def attenuation_curve_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_attenuation_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.attenuation_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_i1_attenuation_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_old_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_sed(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_bulge_sed(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_old_bulge_sed(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_old_bulge_sed(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_bulge_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_bulge_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_bulge_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    # OLD DISK

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the I1 wavelength
        return self.disk_simulations.observed_photometry_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.has_observed_photometry

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

    @property
    def attenuation_curve_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_attenuation_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.attenuation_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_i1_attenuation_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_old_disk(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_sed(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_disk_sed(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_old_disk_sed(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_old_disk_sed(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_disk_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_disk_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_disk_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    # OLD

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_i1_luminosity_old(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the I1 wavelength
        return self.old_simulations.observed_photometry_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_i1_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_i1_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.intrinsic_photometry_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_i1_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.has_intrinsic_photometry

    # -----------------------------------------------------------------

    @property
    def attenuation_curve_old(self):

        """
        This fuction ...
        :return:
        """

        return self.old_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_attenuation_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.attenuation_at(self.i1_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_i1_attenuation_old(self):

        """
        This fnuction ...
        :return:
        """

        return self.old_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def attenuation_i1_old(self):

        """
        This function ...
        :return:
        """

        return self.i1_attenuation_old

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_old(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_old_sed(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_sed(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_old_sed(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_old_sed(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_old_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_old_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.old_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    # YOUNG

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_young(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the FUV wavelength
        return self.young_simulations.observed_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @property
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

    @property
    def attenuation_curve_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_attenuation_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.attenuation_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_fuv_attenuation_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_young_sed(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_young_sed(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_young_sed(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_young_sed(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_young_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_young_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_young_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_young_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    # SFR

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the FUV wavelength
        return self.sfr_simulations.observed_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @property
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

    @property
    def attenuation_curve_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_attenuation_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.attenuation_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_fuv_attenuation_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    # ONLY FOR SFR
    @property
    def intrinsic_dust_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.intrinsic_dust_luminosity

    # -----------------------------------------------------------------

    # ONLY FOR SFR
    @property
    def has_intrinsic_dust_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_intrinsic_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_sed(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_sfr_sed(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sfr_sed(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_sfr_sed(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_sfr_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_sfr_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_sfr_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------

    @property
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

    @property
    def observed_stellar_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.observed_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_observed_stellar_luminosity

    # -----------------------------------------------------------------

    # ONLY FOR SFR (NORMALLY INTRINSIC STELLAR = INTRINSIC BOLOMETRIC)
    @property
    def intrinsic_stellar_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.intrinsic_stellar_luminosity

    # -----------------------------------------------------------------

    # ONLY FOR SFR (NORMALLY INTRINSIC STELLAR = INTRINSIC BOLOMETRIC)
    @property
    def has_intrinsic_stellar_luminosity_sfr(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_intrinsic_stellar_luminosity

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    # UNEVOLVED

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fuv_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        # don't interpolate, wavelength grid is expected to contain the FUV wavelength
        return self.unevolved_simulations.observed_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_observed_fuv_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.has_observed_photometry

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        #self.intrinsic_fuv_luminosity_young + self.intrinsic_fuv_luminosity_sfr
        return self.unevolved_simulations.intrinsic_photometry_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_fuv_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        #return self.has_intrinsic_fuv_luminosity_young and self.has_intrinsic_fuv_luminosity_sfr
        return self.unevolved_simulations.has_intrinsic_photometry

    # -----------------------------------------------------------------

    @property
    def observed_bolometric_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.has_observed_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def intrinsic_bolometric_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.has_intrinsic_bolometric_luminosity

    # -----------------------------------------------------------------

    @property
    def attenuation_curve_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.attenuation_curve

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_attenuation_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.attenuation_at(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def has_fuv_attenuation_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.has_attenuation

    # -----------------------------------------------------------------

    @property
    def observed_dust_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.has_observed_dust_luminosity

    # -----------------------------------------------------------------

    @property
    def bolometric_attenuation_unevolved(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_unevolved(self):

        """
        Thisf unction ...
        :return:
        """

        return self.unevolved_simulations.has_bolometric_attenuation

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_sed(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def observed_unevolved_sed(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.observed_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_unevolved_sed(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_unevolved_sed(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def observed_unevolved_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.has_observed_stellar_sed

    # -----------------------------------------------------------------

    @property
    def observed_unevolved_dust_sed(self):

        """
        Thisn function ...
        :return:
        """

        return self.unevolved_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_unevolved_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_simulations.has_observed_dust_sed

    # -----------------------------------------------------------------
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

    @lazyproperty
    def full_sed_instrument(self):

        """
        This function ...
        :return:
        """

        return FullSEDInstrument.from_properties(self.distance, self.inclination, self.position_angle)

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
    def old_bulge_projections_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.old_bulge_path, projections_dirname)

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
    def old_disk_projections_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.old_disk_path, projections_dirname)

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
    def young_projections_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.young_path, projections_dirname)

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
    def sfr_projections_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.sfr_path, projections_dirname)

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
    def dust_projections_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.dust_path, projections_dirname)

    # -----------------------------------------------------------------

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_sed_ski_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.old_disk_sed_path, disk_simulation_name + ".ski")

    # -----------------------------------------------------------------

    @lazyproperty
    def young_sed_ski_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.young_sed_path, young_simulation_name + ".ski")

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_sed_ski_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.sfr_sed_path, sfr_name + ".ski")

    # -----------------------------------------------------------------

    @property
    def has_old_disk_sed_skifile(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_disk_sed_ski_path)

    # -----------------------------------------------------------------

    @property
    def has_young_sed_skifile(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_sed_ski_path)

    # -----------------------------------------------------------------

    @property
    def has_sfr_sed_skifile(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_sed_ski_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_component(self):

        """
        This function ...
        :return:
        """

        return self.definition.load_bulge_component()

    # -----------------------------------------------------------------

    @property
    def old_bulge_model(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_component.model

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_radial_effective_radius(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_model.effective_radius

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_vertical_effective_radius(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_model.effective_radius * self.old_bulge_model.z_flattening

    # -----------------------------------------------------------------

    @property
    def old_bulge_scaleheight(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_vertical_effective_radius

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_component(self):

        """
        This function ...
        :return:
        """

        return self.definition.load_old_stars_component()

    # -----------------------------------------------------------------

    @property
    def old_disk_deprojection(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_component.deprojection

    # -----------------------------------------------------------------

    @property
    def old_disk_scaleheight(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_deprojection.scale_height

    # -----------------------------------------------------------------

    @lazyproperty
    def young_component(self):

        """
        This function ...
        :return:
        """

        return self.definition.load_young_stars_component()

    # -----------------------------------------------------------------

    @property
    def young_deprojection(self):

        """
        This function ...
        :return:
        """

        return self.young_component.deprojection

    # -----------------------------------------------------------------

    @property
    def young_scaleheight(self):

        """
        This function ...
        :return:
        """

        return self.young_deprojection.scale_height

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_component(self):

        """
        This function ...
        :return:
        """

        return self.definition.load_ionizing_stars_component()

    # -----------------------------------------------------------------

    @property
    def sfr_deprojection(self):

        """
        This function ...
        :return:
        """

        return self.sfr_component.deprojection

    # -----------------------------------------------------------------

    @property
    def sfr_scaleheight(self):

        """
        This function ...
        :return:
        """

        return self.sfr_deprojection.scale_height

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
    def dust_deprojection(self):

        """
        This function ...
        :return:
        """

        return self.dust_component.deprojection

    # -----------------------------------------------------------------

    @property
    def dust_scaleheight(self):

        """
        This function ...
        :return:
        """

        return self.dust_deprojection.scale_height

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
    # BULGE MAPS
    # -----------------------------------------------------------------

    @property
    def old_bulge_i1_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_i1_luminosity_map_earth(self):

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
    def old_bulge_i1_luminosity_map_faceon(self):

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
    def old_bulge_i1_luminosity_map_edgeon(self):

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
    def has_old_bulge_i1_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.has_old_bulge_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_i1_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_i1_luminosity_old_bulge and self.has_old_bulge_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_i1_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_i1_luminosity_old_bulge and self.has_old_bulge_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_i1_luminosity_map_edgeon(self):

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

        """
        This function ...
        :return:
        """

        return self.has_old_bulge_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_bolometric_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_old_bulge and self.has_old_bulge_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_bolometric_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_old_bulge and self.has_old_bulge_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_bulge_bolometric_luminosity_map_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_old_bulge and self.has_old_bulge_map_edgeon

    # -----------------------------------------------------------------
    # DISK MAPS
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
    def has_old_disk_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_disk_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.old_disk_map_path)

    # -----------------------------------------------------------------

    @property
    def old_disk_map_psf_filter(self):

        """
        Thisf unction ...
        :return:
        """

        return self.old_disk_map.psf_filter

    # -----------------------------------------------------------------

    @property
    def old_disk_map_fwhm(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_map.fwhm

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_map_wcs(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.old_disk_map_wcs.shape

    # -----------------------------------------------------------------

    @property
    def old_disk_map_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_map_wcs.pixelscale

    # -----------------------------------------------------------------

    @property
    def old_disk_i1_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def old_disk_i1_luminosity_map_earth(self):

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
    def old_disk_i1_luminosity_map_faceon(self):

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
    def old_disk_i1_luminosity_map_edgeon(self):

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
    def has_old_disk_i1_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.has_old_disk_i1_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_disk_i1_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_i1_luminosity_old_disk and self.has_old_disk_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_disk_i1_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_i1_luminosity_old_disk and self.has_old_disk_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_disk_i1_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_i1_luminosity_old_disk and self.has_old_disk_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_disk_bolometric_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.old_disk_bolometric_luminosity_map_earth

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
    def has_old_disk_bolometric_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.has_old_disk_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_disk_bolometric_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_old_disk and self.has_old_disk_map_earth

    # -----------------------------------------------------------------

    @property
    def has_old_disk_bolometric_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_old_disk and self.has_old_disk_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_disk_bolometric_luminosity_map_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_old_disk and self.has_old_disk_map_edgeon

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
    def has_young_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.young_map_path)

    # -----------------------------------------------------------------

    @property
    def young_map_psf_filter(self):

        """
        This function ...
        :return:
        """

        return self.young_map.psf_filter

    # -----------------------------------------------------------------

    @property
    def young_map_fwhm(self):

        """
        This function ...
        :return:
        """

        return self.young_map.fwhm

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map_wcs(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.young_map_wcs.shape

    # -----------------------------------------------------------------

    @property
    def young_map_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.young_map_wcs.pixelscale

    # -----------------------------------------------------------------

    @property
    def young_fuv_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.young_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def young_fuv_luminosity_map_earth(self):

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
    def young_fuv_luminosity_map_faceon(self):

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
    def young_fuv_luminosity_map_edgeon(self):

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
    def has_young_fuv_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.has_young_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_young_fuv_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_fuv_luminosity_young and self.has_young_map_earth

    # -----------------------------------------------------------------

    @property
    def has_young_fuv_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_fuv_luminosity_young and self.has_young_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_fuv_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_fuv_luminosity_young and self.has_young_map_edgeon

    # -----------------------------------------------------------------

    @property
    def young_bolometric_luminosity_map(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.has_young_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_young_bolometric_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_young and self.has_young_map_earth

    # -----------------------------------------------------------------

    @property
    def has_young_bolometric_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_young and self.has_young_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_bolometric_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_young and self.has_young_map_edgeon

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
    def has_sfr_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_map(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.sfr_map_path)

    # -----------------------------------------------------------------

    @property
    def sfr_map_psf_filter(self):

        """
        This function ...
        :return:
        """

        return self.sfr_map.psf_filter

    # -----------------------------------------------------------------

    @property
    def sfr_map_fwhm(self):

        """
        This function ...
        :return:
        """

        return self.sfr_map.fwhm

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_map_wcs(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.sfr_map_wcs.shape

    # -----------------------------------------------------------------

    @property
    def sfr_map_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.sfr_map_wcs.pixelscale

    # -----------------------------------------------------------------

    @property
    def sfr_fuv_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.sfr_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_fuv_luminosity_map_earth(self):

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
    def sfr_fuv_luminosity_map_faceon(self):

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
    def sfr_fuv_luminosity_map_edgeon(self):

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
    def has_sfr_fuv_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_fuv_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_fuv_luminosity_sfr and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_fuv_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_fuv_luminosity_sfr and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_sfr_fuv_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_fuv_luminosity_sfr and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_bolometric_luminosity_map(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.has_sfr_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_bolometric_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_sfr and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_bolometric_luminosity_map_faceon(self):

        """
        Thsi function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_sfr and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_sfr_bolometric_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_bolometric_luminosity_sfr and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def star_formation_rate_map(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.has_star_formation_rate_map_earth

    # -----------------------------------------------------------------

    @property
    def has_star_formation_rate_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_star_formation_rate_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_star_formation_rate_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_dust_mass_map(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.has_sfr_dust_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr_dust_mass and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr_dust_mass and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr_dust_mass and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_stellar_luminosity_map(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.has_sfr_stellar_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_stellar_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_stellar_luminosity_sfr and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_stellar_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_stellar_luminosity_sfr and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_sfr_stellar_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_stellar_luminosity_sfr and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_dust_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.sfr_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_luminosity_map_earth(self):

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
    def sfr_dust_luminosity_map_faceon(self):

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
    def sfr_dust_luminosity_map_edgeon(self):

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
    def has_sfr_dust_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.has_sfr_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_dust_luminosity_sfr and self.has_sfr_map_earth

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_dust_luminosity_sfr and self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_dust_luminosity_sfr and self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_fuv_luminosity_map(self):

        """
        Thisfunction ...
        :return:
        """

        return self.unevolved_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_fuv_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        # Get contributions
        young = self.young_fuv_luminosity_map_earth
        sfr = self.sfr_fuv_luminosity_map_earth

        # Uniformize
        young, sfr = convolve_and_rebin(young, sfr)

        # Sum the contributions
        return young + sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_fuv_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get contributions
        young = self.young_fuv_luminosity_map_faceon
        sfr = self.sfr_fuv_luminosity_map_faceon

        # Combine (no WCS): regrid and recenter??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_fuv_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get contributions
        young = self.young_fuv_luminosity_map_edgeon
        sfr = self.young_fuv_luminosity_map_edgeon

        # Combine (no WCS): regrid and recenter??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_fuv_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.has_unevolved_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_fuv_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_young_fuv_luminosity_map_earth and self.has_sfr_fuv_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_fuv_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_young_fuv_luminosity_map_faceon and self.has_sfr_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_fuv_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_young_fuv_luminosity_map_edgeon and self.has_sfr_fuv_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_bolometric_luminosity_map(self):

        """
        This function ...
        :return:
        """

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
        young, sfr = convolve_and_rebin(young, sfr)

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

        """
        This function ...
        :return:
        """

        return self.has_unevolved_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_bolometric_luminosity_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_young_bolometric_luminosity_map_earth and self.has_sfr_bolometric_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_bolometric_luminosity_map_faceon(self):

        """
        Thisn function ...
        :return:
        """

        return self.has_young_bolometric_luminosity_map_faceon and self.has_sfr_bolometric_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_bolometric_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_young_bolometric_luminosity_map_edgeon and self.has_sfr_bolometric_luminosity_map_edgeon

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

        # Get the dust map
        frame = self.dust_map_earth.copy()

        # Normalize to the dust mass
        frame.normalize(to=self.dust_mass)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the dust map
        frame = self.dust_map_faceon.copy()

        # Normalize to the dust mass
        frame.normalize(to=self.dust_mass)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the dust map
        frame = self.dust_map_edgeon.copy()

        # Normalize to the dust mass
        frame.normalize(to=self.dust_mass)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @property
    def has_dust_mass(self):

        """
        This function ...
        :return:
        """

        return True

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

        return self.has_dust_map_earth and self.has_dust_mass

    # -----------------------------------------------------------------

    @property
    def has_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_map_edgeon and self.has_dust_mass

    # -----------------------------------------------------------------

    @property
    def has_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_map_faceon and self.has_dust_mass

    # -----------------------------------------------------------------

    @property
    def total_dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return self.total_dust_mass_map_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_mass_map_earth(self):

        """
        This function ...
        :return:
        """

        # Get the contributions
        dust_mass = self.dust_mass_map
        sfr_dust_mass = self.sfr_dust_mass_map

        # Uniformize
        dust_mass, sfr_dust_mass = convolve_and_rebin(dust_mass, sfr_dust_mass)

        # Sum the contributions
        return dust_mass + sfr_dust_mass

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get the contributions
        dust_mass = self.dust_mass_map_faceon
        sfr_dust_mass = self.sfr_dust_mass_map_faceon

        # Combine (no WCS): regrid and center??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get the contributions
        dust_mass = self.dust_mass_map_edgeon
        sfr_dust_mass = self.sfr_dust_mass_map_edgeon

        # Combine (no WCS): regrid and center??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def has_total_dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return self.has_total_dust_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_dust_mass_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_mass_map_earth and self.has_sfr_dust_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def has_total_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_mass_map_faceon and self.has_sfr_dust_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_mass_map_edgeon and self.has_sfr_dust_mass_map_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity_map(self):

        """
        This function
        :return:
        """

        # WRONG!

        # # Get the dust map
        # frame = self.dust_map.copy()
        #
        # # Normalize to the (diffuse) dust luminosity
        # frame.normalize(to=self.diffuse_dust_luminosity)
        #
        # # Return
        # return frame

        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # WRONG!

        # # Get the dust map
        # frame = self.dust_map_faceon.copy()
        #
        # # Normalize to the (diffuse) dust luminosity
        # frame.normalize(to=self.diffuse_dust_luminosity)
        #
        # # Return
        # return frame

        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # WRONG!

        # # Get the dust map
        # frame = self.dust_map_edgeon.copy()
        #
        # # Normalize to the (diffuse) dust luminosity
        # frame.normalize(to=self.diffuse_dust_luminosity)
        #
        # # Return
        # return frame

        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def has_dust_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.has_diffuse_dust_luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_luminosity_map(self):

        """
        This function ...
        :return:
        """

        # WRONG!

        # # Get the contributions
        # dust_lum = self.dust_luminosity_map
        # sfr_dust_lum = self.sfr_dust_luminosity_map
        #
        # # Uniformize
        # dust_lum, sfr_dust_lum = convolve_and_rebin(dust_lum, sfr_dust_lum)
        #
        # # Sum the contributions
        # return dust_lum + sfr_dust_lum

        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_luminosity_map_faceon(self):

        """
        This function ...
        :return:
        """

        # WRONG!

        # Get the contributions
        dust_lum = self.dust_luminosity_map_faceon
        sfr_dust_lum = self.sfr_dust_luminosity_map_faceon

        # Combine (no WCS): regrid and recenter??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_luminosity_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # WRONG!

        # Get the contributions
        dust_lum = self.dust_luminosity_map_edgeon
        sfr_dust_lum = self.sfr_dust_luminosity_map_edgeon

        # Combine (no WCS): regrid and recenter??
        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def has_total_dust_luminosity_map(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_luminosity

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

        return ComponentProjections(bulge_simulation_name, self.old_bulge_component, path=self.old_bulge_projections_path,
                                    description="old bulge stellar component",
                                    projection=self.old_disk_projections.projection_earth,
                                    projection_faceon=self.old_disk_projections.projection_faceon,
                                    projection_edgeon=self.old_disk_projections.projection_edgeon, center=self.center)

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

        return ComponentProjections(disk_simulation_name, self.old_disk_component, path=self.old_disk_projections_path,
                                    earth=False, description="old disk stellar component",
                                    input_filepaths=[self.old_disk_map_path], center=self.center)

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

        return ComponentProjections(young_simulation_name, self.young_component, path=self.young_projections_path,
                                    earth=False, description="young stellar component",
                                    input_filepaths=[self.young_map_path], center=self.center)

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

        return ComponentProjections(sfr_simulation_name, self.sfr_component, path=self.sfr_projections_path,
                                    earth=False, description="SFR component", input_filepaths=[self.sfr_map_path],
                                    center=self.center)

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

        return ComponentProjections(dust_simulation_name, self.dust_component, path=self.dust_projections_path,
                                    earth=False, description="dust component", input_filepaths=[self.dust_map_path])

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

        return self.old_bulge_component.sed_filepath

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

    @lazyproperty
    def derived_parameter_values_total(self):

        """
        Total
        :return:
        """

        # Initialize
        values = OrderedDict()

        # Bolometric luminosity
        if self.has_observed_total_bolometric_luminosity: values[obs_total_bol_lum_name] = self.observed_total_bolometric_luminosity
        if self.has_intrinsic_total_bolometric_luminosity: values[intr_total_bol_lum_name] = self.intrinsic_total_bolometric_luminosity

        # Stellar luminosity
        if self.has_observed_stellar_luminosity: values[obs_total_stellar_bol_lum_name] = self.observed_stellar_luminosity
        if self.has_intrinsic_stellar_luminosity: values[intr_total_stellar_bol_lum_name] = self.intrinsic_stellar_luminosity

        # Total luminosity
        if self.has_bolometric_attenuation: values[bol_attenuation_name] = self.bolometric_attenuation

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_bulge(self):

        """
        Old bulge
        :return:
        """

        # Initialize
        values = OrderedDict()

        # I1 specific luminosity
        if self.has_observed_i1_luminosity_old_bulge: values[obs_bulge_spec_lum_name] = self.observed_i1_luminosity_old_bulge
        if self.has_intrinsic_i1_luminosity_old_bulge: values[intr_bulge_spec_lum_name] = self.intrinsic_i1_luminosity_old_bulge # part of parameter set

        # Bolometric luminosity
        if self.has_observed_bolometric_luminosity_old_bulge: values[obs_bulge_bol_lum_name] = self.observed_bolometric_luminosity_old_bulge
        if self.has_intrinsic_bolometric_luminosity_old_bulge: values[intr_bulge_bol_lum_name] = self.intrinsic_bolometric_luminosity_old_bulge

        # Attenuation
        if self.has_i1_attenuation_old_bulge: values[bulge_spec_attenuation_name] = self.i1_attenuation_old_bulge
        if self.has_bolometric_attenuation_old_bulge: values[bulge_bol_attenuation_name] = self.bolometric_attenuation_old_bulge

        # Dust
        if self.has_observed_dust_luminosity_old_bulge: values[obs_bulge_dust_lum_name] = self.observed_dust_luminosity_old_bulge

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_disk(self):

        """
        Old disk
        :return:
        """

        # Initialize
        values = OrderedDict()

        # I1 specific luminosity
        if self.has_observed_i1_luminosity_old_disk: values[obs_disk_spec_lum_name] = self.observed_i1_luminosity_old_disk
        if self.has_intrinsic_i1_luminosity_old_disk: values[intr_disk_spec_lum_name] = self.intrinsic_i1_luminosity_old_disk # part of parameter set

        # Bolometric luminosity
        if self.has_observed_bolometric_luminosity_old_disk: values[obs_disk_bol_lum_name] = self.observed_bolometric_luminosity_old_disk
        if self.has_intrinsic_bolometric_luminosity_old_disk: values[intr_disk_bol_lum_name] = self.intrinsic_bolometric_luminosity_old_disk

        # Attenuation
        if self.has_i1_attenuation_old_disk: values[disk_spec_attenuation_name] = self.i1_attenuation_old_disk
        if self.has_bolometric_attenuation_old_disk: values[disk_bol_attenuation_name] = self.bolometric_attenuation_old_disk

        # Dust
        if self.has_observed_dust_luminosity_old_disk: values[obs_disk_dust_lum_name] = self.observed_dust_luminosity_old_disk

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_old(self):

        """
        Old (evolved)
        :return:
        """

        # Initialize
        values = OrderedDict()

        # I1 specific luminosity
        if self.has_observed_i1_luminosity_old: values[obs_old_spec_lum_name] = self.observed_i1_luminosity_old
        if self.has_intrinsic_i1_luminosity_old: values[intr_old_spec_lum_name] = self.intrinsic_i1_luminosity_old

        # Bolometric luminosity
        if self.has_observed_bolometric_luminosity_old: values[obs_old_bol_lum_name] = self.observed_bolometric_luminosity_old
        if self.has_intrinsic_bolometric_luminosity_old: values[intr_old_bol_lum_name] = self.intrinsic_bolometric_luminosity_old

        # Attenuation
        if self.has_i1_attenuation_old: values[old_spec_attenuation_name] = self.i1_attenuation_old
        if self.has_bolometric_attenuation_old: values[old_bol_attenuation_name] = self.bolometric_attenuation_old

        # Dust
        if self.has_observed_dust_luminosity_old: values[obs_old_dust_lum_name] = self.observed_dust_luminosity_old

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_young(self):

        """
        Young stars
        :return:
        """

        # Initialize
        values = OrderedDict()

        # FUV specific luminosity
        if self.has_observed_fuv_luminosity_young: values[obs_young_spec_lum_name] = self.observed_fuv_luminosity_young
        if self.has_intrinsic_fuv_luminosity_young: values[intr_young_spec_lum_name] = self.intrinsic_fuv_luminosity_young # part of (free) parameter set

        # Bolometric luminosity
        if self.has_observed_bolometric_luminosity_young: values[obs_young_bol_lum_name] = self.observed_bolometric_luminosity_young
        if self.has_intrinsic_bolometric_luminosity_young: values[intr_young_bol_lum_name] = self.intrinsic_bolometric_luminosity_young

        # Attenuation
        if self.has_fuv_attenuation_young: values[young_spec_attenuation_name] = self.fuv_attenuation_young
        if self.has_bolometric_attenuation_young: values[young_bol_attenuation_name] = self.bolometric_attenuation_young

        # Dust
        if self.has_observed_dust_luminosity_young: values[obs_young_dust_lum_name] = self.observed_dust_luminosity_young

        # Return
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def derived_parameter_values_sfr(self):

        """
        Ionizing stars (SFR)
        :return:
        """

        # Initialize
        values = OrderedDict()

        # SFR
        if self.has_sfr: values[sfr_name] = self.sfr

        # FUV specific luminosity
        if self.has_observed_fuv_luminosity_sfr: values[obs_sfr_spec_lum_name] = self.observed_fuv_luminosity_sfr
        if self.has_intrinsic_fuv_luminosity_sfr: values[intr_sfr_spec_lum_name] = self.intrinsic_fuv_luminosity_sfr # part of the (free) parameter set

        # Bolometric luminosity
        if self.has_observed_bolometric_luminosity_sfr: values[obs_sfr_bol_lum_name] = self.observed_bolometric_luminosity_sfr
        if self.has_intrinsic_bolometric_luminosity_sfr: values[intr_sfr_bol_lum_name] = self.intrinsic_bolometric_luminosity_sfr

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
    def derived_parameter_values_unevolved(self):

        """
        Young + ionizing (unevolved)
        :return:
        """

        # Initialize
        values = OrderedDict()

        # FUV specific luminosity
        if self.has_observed_fuv_luminosity_unevolved: values[obs_unevolved_spec_lum_name] = self.observed_fuv_luminosity_unevolved
        if self.has_intrinsic_fuv_luminosity_unevolved: values[intr_unevolved_spec_lum_name] = self.intrinsic_fuv_luminosity_unevolved

        # Bolometric luminosity
        if self.has_observed_bolometric_luminosity_unevolved: values[obs_unevolved_bol_lum_name] = self.observed_bolometric_luminosity_unevolved
        if self.has_intrinsic_bolometric_luminosity_unevolved: values[intr_unevolved_bol_lum_name] = self.intrinsic_bolometric_luminosity_unevolved

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

        # Luminosity
        if self.has_dust_luminosity: values[dust_lum_name] = self.dust_luminosity
        if self.has_diffuse_dust_luminosity: values[diffuse_dust_lum_name] = self.diffuse_dust_luminosity

        # Mass
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

    @property
    def has_old_bulge_cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.has_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def old_bulge_cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.bulge_simulations.cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def has_old_disk_cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.has_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def old_disk_cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.disk_simulations.cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def has_young_cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def young_cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def has_sfr_cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.has_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def sfr_cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.sfr_simulations.cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def has_cell_properties(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_cell_properties

    # -----------------------------------------------------------------

    @property
    def cell_volumes(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.cell_volumes

    # -----------------------------------------------------------------

    @property
    def cell_dust_densities(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.cell_dust_densities

    # -----------------------------------------------------------------

    @property
    def cell_mass_fractions(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.cell_mass_fractions

    # -----------------------------------------------------------------

    @property
    def cell_optical_depths(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.cell_optical_depths

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_x_coordinates(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.cell_x_coordinates

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_y_coordinates(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.cell_y_coordinates

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_z_coordinates(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.cell_z_coordinates

    # -----------------------------------------------------------------

    @property
    def has_grid_files(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_grid_files

    # -----------------------------------------------------------------

    @property
    def grid_filepaths(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.grid_filepaths

    # -----------------------------------------------------------------

    @property
    def grid_xy_filepath(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.grid_xy_filepath

    # -----------------------------------------------------------------

    @property
    def grid_xz_filepath(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.grid_xz_filepath

    # -----------------------------------------------------------------

    @property
    def grid_yz_filepath(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.grid_yz_filepath

    # -----------------------------------------------------------------

    @property
    def grid_xyz_filepath(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.grid_xyz_filepath

# -----------------------------------------------------------------
