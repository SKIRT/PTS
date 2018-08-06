#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.simulation.simulation Contains the ComponentSimulation, IntrinsicComponentSimulation, and
#  ObservedComponentSimulation class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractproperty
from collections import OrderedDict
import numpy as np

# Import the relevant PTS classes and modules
from ...core.simulation.simulation import SkirtSimulation
from ...core.simulation.simulation import createsimulations
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.simulation.output import SimulationOutput
from ...core.simulation.data import SimulationData
from ...core.tools import filesystem as fs
from ...core.tools import sequences

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

    def __init__(self, *args, **kwargs):

        """
        Thisf unction ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ComponentSimulation, self).__init__(*args, **kwargs)

        # Set earth coordinate system attribute
        self.earth_wcs = None

    # -----------------------------------------------------------------

    @classmethod
    def from_output_path(cls, path, name=None, earth_wcs=None):

        """
        This function ...
        :param path:
        :param name:
        :param earth_wcs:
        :return:
        """

        sim = createsimulations(path, single=True, name=name, cls=cls)
        if earth_wcs is not None: sim.earth_wcs = earth_wcs
        return sim

    # -----------------------------------------------------------------

    @lazyproperty
    def ski_file(self):

        """
        This function ...
        :return:
        """

        return self.parameters()

    # -----------------------------------------------------------------

    @lazyproperty
    def instrument_names(self):

        """
        This function ...
        :return:
        """

        return self.ski_file.instrumentnames()

    # -----------------------------------------------------------------

    @lazyproperty
    def instruments(self):

        """
        This function ...
        :return:
        """

        instruments = OrderedDict()
        for name in self.instrument_names: instruments[name] = self.ski_file.get_instrument_object(name)
        return instruments

    # -----------------------------------------------------------------

    @lazyproperty
    def instrument_distances(self):

        """
        This function ...
        :return:
        """

        distances = OrderedDict()
        for name in self.instrument_names: distances[name] = self.instruments[name].distance
        return distances

    # -----------------------------------------------------------------

    @lazyproperty
    def distance(self):

        """
        This function ...
        :return:
        """

        return sequences.get_all_close_value(self.instrument_distances.values(), ignore_none=True, return_none=True)

    # -----------------------------------------------------------------

    @property
    def has_output(self):

        """
        This function ...
        :return:
        """

        return not fs.is_empty(self.output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def output(self):

        """
        This function ...
        :return:
        """

        return SimulationOutput.from_directory(self.output_path, prefix=self.prefix())

    # -----------------------------------------------------------------

    @property
    def coordinate_systems(self):

        """
        This function ...
        :return:
        """

        systems = dict()
        if self.earth_wcs is not None: systems[earth_name] = self.earth_wcs
        return systems

    # -----------------------------------------------------------------

    @property
    def distances(self):

        """
        This function ...
        :return:
        """

        distances = dict()
        distances[earth_name] = self.distance
        distances[faceon_name] = self.distance
        distances[edgeon_name] = self.distance
        return distances

    # -----------------------------------------------------------------

    @lazyproperty
    def data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.output, coordinate_systems=self.coordinate_systems, distances=self.distances)

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

        if "X coordinate of cell center" in self.cell_properties_columns: return np.asarray(self.cell_properties["X coordinate of cell center"]) # SKIRT 8
        elif "X coordinate of cell center" in self.cell_absorptions_columns: return np.asarray(self.cell_absorptions["X coordinate of cell center"]) # SKIRT 7
        else: raise IOError("")

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_y_coordinates(self):

        """
        This function ...
        :return:
        """

        if "Y coordinate of cell center" in self.cell_properties_columns: return np.asarray(self.cell_properties["Y coordinate of cell center"]) # SKIRT 8
        elif "Y coordinate of cell center" in self.cell_absorptions_columns: return np.asarray(self.cell_absorptions["Y coordinate of cell center"]) # SKIRT 7
        else: raise IOError("")

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_z_coordinates(self):

        """
        This function ...
        :return:
        """

        if "Z coordinate of cell center" in self.cell_properties_columns: return np.asarray(self.cell_properties["Z coordinate of cell center"]) # SKIRT 8
        elif "Z coordinate of cell center" in self.cell_absorptions_columns: return np.asarray(self.cell_absorptions["Z coordinate of cell center"]) # SKIRT 7
        else: raise IOError("")

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

    @property
    def has_cell_absorptions(self):

        """
        This function ...
        :return:
        """

        return self.data.has_absorption

    # -----------------------------------------------------------------

    @property
    def cell_absorptions(self):

        """
        This function ...
        :return:
        """

        return self.data.absorption

    # -----------------------------------------------------------------

    @property
    def cell_absorptions_columns(self):

        """
        This function ...
        :return:
        """

        return self.cell_absorptions.colnames

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
