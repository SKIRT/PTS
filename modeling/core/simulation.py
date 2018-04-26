#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.model Contains the IntrinsicComponentSimulation, ObservedComponentSimulation,
#  and ComponentSimulation class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractproperty
from collections import OrderedDict
import numpy as np

# Import the relevant PTS classes and modules
from ...core.simulation.simulation import SkirtSimulation
from ...core.units.unit import parse_unit as u
from ...core.simulation.simulation import createsimulations
from ...core.tools.utils import lazyproperty
from ...core.simulation.output import SimulationOutput
from ...core.simulation.data import SimulationData
from ...core.tools import filesystem as fs
from ...core.tools import sequences
from ...core.data.sed import load_sed

# -----------------------------------------------------------------

stellar_dust_sed_split_wavelength = 5. * u("micron")

# -----------------------------------------------------------------

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

class ComponentSimulation(SkirtSimulation):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    @classmethod
    def from_output_path(cls, path, name=None):

        """
        This function ...
        :param path:
        :param name:
        :return:
        """

        return createsimulations(path, single=True, name=name, cls=cls)

    # -----------------------------------------------------------------

    @property
    def has_output(self):

        """
        This function ...
        :return:
        """

        return fs.is_empty(self.output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def output(self):

        """
        This function ...
        :return:
        """

        return SimulationOutput.from_directory(self.output_path, prefix=self.prefix())

    # -----------------------------------------------------------------

    @lazyproperty
    def data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.output)

    # -----------------------------------------------------------------

    @property
    def has_cell_properties(self):

        """
        This function ...
        :return:
        """

        return self.data.has_cell_properties

    # -----------------------------------------------------------------

    @property
    def cell_properties(self):

        """
        This function ...
        :return:
        """

        return self.data.cell_properties

    # -----------------------------------------------------------------

    @property
    def cell_properties_columns(self):

        """
        This function ...
        :return:
        """

        return self.cell_properties.colnames

    # -----------------------------------------------------------------

    @property
    def cell_volumes(self):

        """
        This function ...
        :return:
        """

        if "Volume" in self.cell_properties_columns: return np.asarray(self.cell_properties["Volume"])  # SKIRT 7
        elif "Cell volume" in self.cell_properties_columns: return np.asarray(self.cell_properties["Cell volume"])  # SKIRT 8
        else: raise IOError("")

    # -----------------------------------------------------------------

    @property
    def cell_dust_densities(self):

        """
        This function ...
        :return:
        """

        if "Density" in self.cell_properties_columns: return np.asarray(self.cell_properties["Density"])  # SKIRT 7
        elif "Average dust density in cell" in self.cell_properties_columns: return np.asarray(self.cell_properties["Average dust density in cell"])  # SKIRT 8
        else: raise IOError("")

    # -----------------------------------------------------------------

    @property
    def cell_mass_fractions(self):

        """
        This function ...
        :return:
        """

        return np.asarray(self.cell_properties["Mass fraction"])

    # -----------------------------------------------------------------

    @property
    def cell_optical_depths(self):

        """
        This function ...
        :return:
        """

        return np.asarray(self.cell_properties["Optical depth"])

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_x_coordinates(self):

        """
        This function ...
        :return:
        """

        return np.asarray(self.cell_properties["X coordinate of cell center"])

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_y_coordinates(self):

        """
        This function ...
        :return:
        """

        return np.asarray(self.cell_properties["Y coordinate of cell center"])

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_z_coordinates(self):

        """
        This function ...
        :return:
        """

        return np.asarray(self.cell_properties["Z coordinate of cell center"])

    # -----------------------------------------------------------------

    @property
    def has_grid_files(self):

        """
        This function ...
        :return:
        """

        return self.output.has_grids

    # -----------------------------------------------------------------

    @property
    def grid_filepaths(self):

        """
        This function ...
        :return:
        """

        return self.output.grids

    # -----------------------------------------------------------------

    @property
    def grid_xy_filepath(self):

        """
        This function ...
        :return:
        """

        return sequences.pick_contains(self.grid_filepaths, "gridxy.")

    # -----------------------------------------------------------------

    @property
    def grid_xz_filepath(self):

        """
        This function ...
        :return:
        """

        return sequences.pick_contains(self.grid_filepaths, "gridxz.")

    # -----------------------------------------------------------------

    @property
    def grid_yz_filepath(self):

        """
        This function ...
        :return:
        """

        return sequences.pick_contains(self.grid_filepaths, "gridyz.")

    # -----------------------------------------------------------------

    @property
    def grid_xyz_filepath(self):

        """
        This function ...
        :return:
        """

        return sequences.pick_contains(self.grid_filepaths, "gridxyz.")

    # -----------------------------------------------------------------

    @property
    def has_cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.data.has_stellar_density

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return np.asarray(self.data.stellar_density["Stellar density"])

# -----------------------------------------------------------------

class IntrinsicComponentSimulation(ComponentSimulation):

    """
    This class ...
    """

# -----------------------------------------------------------------

class ObservedComponentSimulation(ComponentSimulation):

    """
    This class ...
    """

# -----------------------------------------------------------------

class ComponentSimulations(object):

    """
    This function ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, name, observed):

        """
        The constructor ...
        :param name:
        :param observed:
        """

        # Set the name
        self.name = name

        # Set the observed simulation
        self.observed = observed

    # -----------------------------------------------------------------

    @property
    def observed_output_path(self):

        """
        This function ...
        :return:
        """

        return self.observed.output_path

    # -----------------------------------------------------------------

    @property
    def observed_output(self):

        """
        This function ...
        :return:
        """

        return self.observed.output

    # -----------------------------------------------------------------

    @property
    def observed_data(self):

        """
        This function ...
        :return:
        """

        return self.observed.data

    # -----------------------------------------------------------------

    @property
    def has_observed_output(self):

        """
        This function ...
        :return:
        """

        return self.observed.has_output

    # -----------------------------------------------------------------

    @property
    def has_observed_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_output and self.observed_data.has_seds

    # -----------------------------------------------------------------

    @property
    def observed_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.sed_paths_instruments[earth_name]

    # -----------------------------------------------------------------

    @property
    def has_other_observed_sed_orientations(self):

        """
        This function ...
        :return:
        """

        return len(self.observed_data.sed_paths_instruments) > 1

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_sed_orientation(self):

        """
        This function ...
        :return:
        """

        return faceon_name in self.observed_data.sed_paths_instruments

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_sed_orientation(self):

        """
        This function ...
        :return:
        """

        return edgeon_name in self.observed_data.sed_paths_instruments

    # -----------------------------------------------------------------

    @property
    def faceon_observed_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.sed_paths_instruments[faceon_name]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.sed_paths_instruments[edgeon_name]

    # -----------------------------------------------------------------

    @property
    def observed_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[faceon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[edgeon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def has_other_observed_sed_contributions(self):

        """
        This function ...
        :return:
        """

        return len(self.observed_data.seds[earth_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_transparent_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_other_observed_sed_contributions

    # -----------------------------------------------------------------

    @property
    def observed_sed_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][direct_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_sed_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_sed_dust(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][dust_contribution]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_dust_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_dust - self.observed_sed_dust_scattered

    # -----------------------------------------------------------------

    @property
    def observed_sed_dust_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][dust_scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_sed_transparent(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][transparent_contribution]

    # -----------------------------------------------------------------

    @abstractproperty
    def intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def has_cell_properties(self):

        """
        This function ...
        :return:
        """

        return self.observed.has_cell_properties

    # -----------------------------------------------------------------

    @property
    def cell_volumes(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_volumes

    # -----------------------------------------------------------------

    @property
    def cell_dust_densities(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_dust_densities

    # -----------------------------------------------------------------

    @property
    def cell_mass_fractions(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_mass_fractions

    # -----------------------------------------------------------------

    @property
    def cell_optical_depths(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_optical_depths

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_x_coordinates(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_x_coordinates

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_y_coordinates(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_y_coordinates

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_z_coordinates(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_z_coordinates

    # -----------------------------------------------------------------

    @property
    def has_grid_files(self):

        """
        This function ...
        :return:
        """

        return self.observed.has_grid_files

    # -----------------------------------------------------------------

    @property
    def grid_xy_filepath(self):

        """
        This function ...
        :return:
        """

        return self.observed.grid_xy_filepath

    # -----------------------------------------------------------------

    @property
    def grid_xz_filepath(self):

        """
        This function ...
        :return:
        """

        return self.observed.grid_xz_filepath

    # -----------------------------------------------------------------

    @property
    def grid_yz_filepath(self):

        """
        This function ...
        :return:
        """

        return self.observed.grid_yz_filepath

    # -----------------------------------------------------------------

    @property
    def grid_xyz_filepath(self):

        """
        This function ...
        :return:
        """

        return self.observed.grid_xyz_filepath

    # -----------------------------------------------------------------

    @property
    def has_cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.observed.has_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_stellar_density

# -----------------------------------------------------------------

class SingleComponentSimulations(ComponentSimulations):
    
    """
    Objects of this class describe the simulation(s) of radiative transfer model of a certain stellar component.
    """

    def __init__(self, name, observed, intrinsic=None):

        """
        This function ...
        :param name:
        :param intrinsic:
        :param observed:
        """

        # Call the constructor of the base class
        super(SingleComponentSimulations, self).__init__(name, observed)

        # Set the intrinsic simulation
        self.intrinsic = intrinsic

    # -----------------------------------------------------------------

    @classmethod
    def from_output_paths(cls, name, observed, intrinsic=None):

        """
        This function ...
        :param name:
        :param observed:
        :param intrinsic:
        :return:
        """

        # Load simulations
        observed = ObservedComponentSimulation.from_output_path(observed)
        intrinsic = IntrinsicComponentSimulation.from_output_path(intrinsic) if intrinsic is not None else None

        # Create and return
        return cls(name, observed, intrinsic=intrinsic)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic is not None

    # -----------------------------------------------------------------

    @property
    def intrinsic_output_path(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic.output_path if self.has_intrinsic else None

    # -----------------------------------------------------------------

    @property
    def intrinsic_output(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic.output

    # -----------------------------------------------------------------

    @property
    def intrinsic_data(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic.data

    # -----------------------------------------------------------------

    @property
    def intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # Check whether transparent SED is created
        if not self.has_transparent_sed: raise ValueError("Intrinsic SED cannot be calculated")

        # Return
        return self.observed_sed_transparent

    # -----------------------------------------------------------------

    @property
    def faceon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # ISOTROPIC RADIATION
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def edgeon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # ISOTROPIC RADIATION
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_other_observed_sed_contributions

# -----------------------------------------------------------------

class MultiComponentSimulations(ComponentSimulations):

    """
    Objects of this class describe the simulation(s) of a radiative transfer model containing of multiple stellar components.
    """

    def __init__(self, name, observed, intrinsic_seds=None):

        """
        This function ...
        :param name:
        :param observed:
        :param intrinsic_seds:
        """

        # Call the constructor of the base class
        super(MultiComponentSimulations, self).__init__(name, observed)

        # Set the SEDs of the components
        self.intrinsic_seds = intrinsic_seds

    # -----------------------------------------------------------------

    @classmethod
    def from_output_path(cls, name, observed, intrinsic_sed_paths=None):

        """
        This function ...
        :param name:
        :param observed:
        :param intrinsic_sed_paths:
        :return:
        """

        # Load observed simulation
        observed = ObservedComponentSimulation.from_output_path(observed)

        # Load intrinsic SEDs
        if intrinsic_sed_paths is not None:
            intrinsic_seds = OrderedDict()
            for component_name in intrinsic_sed_paths: intrinsic_seds[component_name] = load_sed(intrinsic_sed_paths[component_name])
        else: intrinsic_seds = None

        # Create and return
        return cls(name, observed, intrinsic_seds=intrinsic_seds)

    # -----------------------------------------------------------------

    @property
    def nintrinsic_seds(self):

        """
        This function ...
        :return:
        """

        return len(self.intrinsic_seds)

    # -----------------------------------------------------------------

    @property
    def component_names(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_seds.keys()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_seds(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_seds is not None and self.nintrinsic_seds > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # Transparent SED is written out
        if self.has_transparent_sed: return self.observed_sed_transparent

        # Has intrinsic SEDs, add them
        elif self.has_intrinsic_seds: return sum(self.intrinsic_seds.values())

        # Cannot be calculated
        else: raise ValueError("Intrinsic SED cannot be calculated")

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # ISOTROPIC RADIATION
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # ISOTROPIC RADIATION
        return self.intrinsic_sed

# -----------------------------------------------------------------
