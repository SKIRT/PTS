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
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.simulation.output import SimulationOutput
from ...core.simulation.data import SimulationData
from ...core.tools import filesystem as fs
from ...core.tools import sequences
from ...core.data.sed import load_sed
from ...core.data.attenuation import AttenuationCurve
from ...magic.tools import extinction
from ...magic.core.frame import Frame
from ...magic.core.datacube import DataCube
from ...core.basics.log import log

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

    def __init__(self, name, observed, distance=None):

        """
        The constructor ...
        :param name:
        :param observed:
        :param distance:
        """

        # Set the name
        self.name = name

        # Set the observed simulation
        self.observed = observed

        # Set the distance
        if distance is not None: self.distance = distance

    # -----------------------------------------------------------------

    @property
    def has_observed(self):

        """
        This function ...
        :return:
        """

        return self.observed is not None

    # -----------------------------------------------------------------

    @property
    def distance(self):

        """
        This function ...
        :return:
        """

        return self.observed.distance

    # -----------------------------------------------------------------

    @distance.setter
    def distance(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        if not self.has_observed: log.warning("The simulation of the observed component does not exist: cannot set the distance")
        else: self.observed.distance = value

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

        return self.has_observed and self.observed.has_output

    # -----------------------------------------------------------------

    @property
    def has_observed_data(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_output and self.observed_data.has_any

    # -----------------------------------------------------------------

    @property
    def has_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed and self.has_observed_output and self.observed_data.has_seds

    # -----------------------------------------------------------------

    @property
    def has_observed_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_output and self.observed_data.has_images

    # -----------------------------------------------------------------

    @abstractproperty
    def has_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        pass

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

    @lazyproperty
    def wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed.wavelength_grid()

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

        return self.has_observed_data and len(self.observed_data.seds[earth_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_other_observed_sed_contributions_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and len(self.observed_data.seds[faceon_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_other_observed_sed_contributions_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and len(self.observed_data.seds[edgeon_name]) > 1

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

    @property
    def observed_cube_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.image_paths_instruments[earth_name]

    # -----------------------------------------------------------------

    @property
    def has_other_observed_cube_orientations(self):

        """
        This function ...
        :return:
        """

        return len(self.observed_data.image_paths_instruments) > 1

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_cube_orientation(self):

        """
        This function ...
        :return:
        """

        return faceon_name in self.observed_data.image_paths_instruments

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_cube_orientation(self):

        """
        This function ...
        :return:
        """

        return edgeon_name in self.observed_data.image_paths_instruments

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.image_paths_instruments[faceon_name]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.image_paths_instruments[edgeon_name]

    # -----------------------------------------------------------------

    @property
    def observed_cube(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def has_other_observed_cube_contributions(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and len(self.observed_data.images[earth_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_other_observed_cube_contributions_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and len(self.observed_data.images[faceon_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_other_observed_cube_contributions_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and len(self.observed_data.images[edgeon_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_transparent_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_other_observed_cube_contributions

    # -----------------------------------------------------------------

    @property
    def observed_cube_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][direct_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][direct_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][direct_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_cube_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_cube_dust(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][dust_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_dust(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][dust_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_dust(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][dust_contribution]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_cube_dust_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube_dust - self.observed_cube_dust_scattered

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_cube_dust_direct(self):

        """
        This function ...
        :return:
        """

        return self.faceon_observed_cube_dust - self.faceon_observed_cube_dust_scattered

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_observed_cube_dust_direct(self):

        """
        This function ...
        :return:
        """

        return self.edgeon_observed_cube_dust - self.edgeon_observed_cube_dust_scattered

    # -----------------------------------------------------------------

    @property
    def observed_cube_dust_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][dust_scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_dust_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][dust_scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_dust_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][dust_scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_cube_transparent(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][transparent_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_transparent(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][transparent_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_transparent(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][transparent_contribution]

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

    @lazyproperty
    def observed_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed.splice(max_wavelength=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_cube(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube.splice(max_wavelength=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_cube

    # -----------------------------------------------------------------

    @abstractproperty
    def intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def faceon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def edgeon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def faceon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def edgeon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_sed(self):

        """
        Thisf unction ...
        :return:
        """

        return self.intrinsic_sed.splice(max_wavelength=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_cube(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_cube.splice(max_wavelength=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bolometric_frame(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_bolometric_frame(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_frame(self):

        """
        This function ...
        :return:
        """

        return self.observed_stellar_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_frame(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_sed(self):

        """
        This function ...
        :return:
        """

        # TODO: use self.observed_sed_dust IF PRESENT
        return self.observed_sed.splice(min_wavelength=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_cube(self):

        """
        This function ...
        :return:
        """

        # TODO: use self.observed_cube_dust IF PRESENT
        return self.observed_cube.splice(min_wavelength=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_dust_sed(self):

        """
        This function ...
        :return:
        """

        # TODO: use self.faceon_observed_sed_dust IF PRESENT
        return self.faceon_observed_sed.splice(min_wavelength=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_dust_cube(self):

        """
        This function ...
        :return:
        """

        # TODO: use self.faceon_observed_cube_dust IF PRESENT
        return self.faceon_observed_cube.splice(min_wavelength=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed.splice(min_wavelength=stellar_dust_sed_split_wavelength)

    # ----------------------------------------------------------------

    @lazyproperty
    def intrinsic_dust_cube(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_cube.splice(min_wavelength=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_dust_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_dust_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_frame(self):

        """
        This function ...
        :return:
        """

        return self.observed_dust_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_dust_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_dust_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_dust_frame(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_dust_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_dust_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_dust_cube

    # -----------------------------------------------------------------

    @memoize_method
    def observed_photometry_at(self, wavelength, interpolate=True):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        return self.observed_sed.photometry_at(wavelength, interpolate=interpolate)

    # -----------------------------------------------------------------

    @memoize_method
    def observed_frame_at(self, wavelength, interpolate=False):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        if interpolate: raise NotImplementedError("Interpolating is not supported")
        return self.observed_cube.get_frame_for_wavelength(wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_photometry(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_cube

    # -----------------------------------------------------------------

    @memoize_method
    def intrinsic_photometry_at(self, wavelength, interpolate=True):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        return self.intrinsic_sed.photometry_at(wavelength, interpolate=interpolate)

    # -----------------------------------------------------------------

    @memoize_method
    def intrinsic_frame_at(self, wavelength, interpolate=False):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        if interpolate: raise NotImplementedError("Interpolating is not supported")
        return self.intrinsic_cube.get_frame_for_wavelength(wavelength)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_photometry(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def attenuation_curve(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_stellar_sed
        intrinsic = self.intrinsic_stellar_sed
        return AttenuationCurve.from_seds(observed, intrinsic)

    # -----------------------------------------------------------------

    @lazyproperty
    def attenuation_cube(self):

        """
        This function ...
        :return:
        """

        # Get the datacubes as arrays
        observed = self.observed_stellar_cube.asarray() # wavelength axis first
        intrinsic = self.intrinsic_stellar_cube.asarray()  # wavelength axis first

        # Calculate the attenuations
        attenuation = extinction.attenuation(observed, intrinsic)

        # Create and return the datacube of attenuations
        return DataCube.from_array(attenuation, wavelength_grid=self.wavelength_grid)

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_attenuation(self):

        """
        This function ...
        :return:
        """

        # Get the observed stellar luminosity
        observed = self.observed_stellar_luminosity.to("W", distance=self.distance).value

        # Get the intrinsic stellar luminosity
        #intrinsic = self.intrinsic_stellar_luminosity.to("W", distance=self.distance).value
        # includes intrinsic dust emission in template, so also ACTUAL intrinsic bolometric lum
        intrinsic = self.intrinsic_bolometric_luminosity.to("W", distance=self.distance).value

        # Calculate and return the extinction
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_attenuation_frame(self):

        """
        This function ...
        :return:
        """

        # Get the observed stellar frame
        observed = self.observed_stellar_frame
        observed.convert_to("W", distance=self.distance)

        # Get the intrinsic frame
        intrinsic = self.intrinsic_bolometric_frame
        intrinsic.convert_to("W", distance=self.distance)

        # Calculate and return the extinction
        attenuation = extinction.attenuation(observed.data, intrinsic.data)
        return Frame(attenuation, wcs=observed.wcs)

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_stellar_luminosity and self.has_intrinsic_bolometric_luminosity #self.has_intrinsic_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_stellar_frame

    # -----------------------------------------------------------------

    @memoize_method
    def attenuation_at(self, wavelength, interpolate=True):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        # Get the observed luminosity
        observed = self.observed_photometry_at(wavelength, interpolate=interpolate).to("W/micron", wavelength=wavelength, distance=self.distance).value

        # Get the intrinsic luminosity
        intrinsic = self.intrinsic_photometry_at(wavelength, interpolate=interpolate).to("W/micron", wavelength=wavelength, distance=self.distance).value

        # Calculate and return the extinction
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @memoize_method
    def attenuation_at_frame(self, wavelength, interpolate=True):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        # Get the observed frame
        observed = self.observed_frame_at(wavelength, interpolate=interpolate)
        observed.convert_to("W/micron", wavelength=wavelength, distance=self.distance)

        # Get the intrinsic frame
        intrinsic = self.intrinsic_frame_at(wavelength, interpolate=interpolate)
        intrinsic.convert_to("W/micron", wavelength=wavelength, distance=self.distance)

        # Calculate and return the extinction
        attenuation = extinction.attenuation(observed.data, intrinsic.data)
        return Frame(attenuation, wcs=observed.wcs, unit="W/micron", wavelength=wavelength)

    # -----------------------------------------------------------------

    @property
    def has_attenuation(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed and self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_attenuation_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_cube and self.has_intrinsic_cube

# -----------------------------------------------------------------

class SingleComponentSimulations(ComponentSimulations):
    
    """
    Objects of this class describe the simulation(s) of radiative transfer model of a certain stellar component.
    """

    def __init__(self, name, observed, intrinsic=None, distance=None):

        """
        This function ...
        :param name:
        :param intrinsic:
        :param observed:
        :param distance:
        """

        # Call the constructor of the base class
        super(SingleComponentSimulations, self).__init__(name, observed, distance=distance)

        # Set the intrinsic simulation
        self.intrinsic = intrinsic

    # -----------------------------------------------------------------

    @classmethod
    def from_output_paths(cls, name, observed, intrinsic=None, distance=None):

        """
        This function ...
        :param name:
        :param observed:
        :param intrinsic:
        :param distance:
        :return:
        """

        # Load observed simulation
        if not fs.is_directory(observed): raise ValueError("Observed simulation directory does not exist")
        if fs.has_files_in_path(observed): observed = ObservedComponentSimulation.from_output_path(observed)
        else:
            log.warning("Observed simulation has not been performed: no output files")
            observed = None

        # Load intrinsic simulation
        if intrinsic is not None:
            if not fs.is_directory(intrinsic): raise ValueError("Intrinsic simulation directory does not exist")
            if fs.has_files_in_path(intrinsic): intrinsic = IntrinsicComponentSimulation.from_output_path(intrinsic)
            else:
                log.warning("Intrinsic simulation has not been performed: no output files")
                intrinsic = None

        # Create and return
        return cls(name, observed, intrinsic=intrinsic, distance=distance)

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
    def has_intrinsic_sed(self):

        """
        This function ..
        :return:
        """

        return self.has_transparent_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_transparent_cube

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
    def intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Check whether transparant SED is created
        if not self.has_transparent_cube: raise ValueError("Intrinsic cube cannot be calculated")

        # Return
        return self.observed_cube_transparent

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
    def faceon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Check
        #if not self.has_faceon_transparent_cube: raise ValueError("Intrinsic cube from face-on orientation cannot be calculated")

        # Return
        #return self.faceon_observed_cube_transparent

        raise NotImplementedError("Not implemented")

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

    @property
    def edgeon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Check
        #if not self.has_edgeon_transparent_cube: raise ValueError("Intrinsic cube from edge-on orientation cannot be calculated")

        # Return
        #return self.edgeon_observed_cube_transparent

        raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

class MultiComponentSimulations(ComponentSimulations):

    """
    Objects of this class describe the simulation(s) of a radiative transfer model containing of multiple stellar components.
    """

    def __init__(self, name, observed, intrinsic_seds=None, intrinsic_cubes=None, distance=None):

        """
        This function ...
        :param name:
        :param observed:
        :param intrinsic_seds:
        :param intrinsic_cubes:
        :param distance:
        """

        # Call the constructor of the base class
        super(MultiComponentSimulations, self).__init__(name, observed, distance=distance)

        # Set the SEDs of the components
        self.intrinsic_seds = intrinsic_seds

        # Set the datacubes of the components
        self.intrinsic_cubes = intrinsic_cubes

    # -----------------------------------------------------------------

    @classmethod
    def from_output_path(cls, name, observed, intrinsic_sed_paths=None, intrinsic_cube_paths=None, distance=None):

        """
        This function ...
        :param name:
        :param observed:
        :param intrinsic_sed_paths:
        :param intrinsic_cube_paths:
        :param distance:
        :return:
        """

        # Load observed simulation
        observed = ObservedComponentSimulation.from_output_path(observed)

        # Load intrinsic SEDs
        if intrinsic_sed_paths is not None:
            intrinsic_seds = OrderedDict()
            for component_name in intrinsic_sed_paths: intrinsic_seds[component_name] = load_sed(intrinsic_sed_paths[component_name])
        else: intrinsic_seds = None

        # Load intrinsic cubes
        if intrinsic_cube_paths is not None:
            if intrinsic_seds is None: raise ValueError("When passing the filepaths of datacubes, the filepaths of corresponding simulated SEDs also have to be specified (for the wavelength grid)")
            intrinsic_cubes = OrderedDict()
            for component_name in intrinsic_cube_paths:
                if component_name not in intrinsic_seds: raise ValueError("SED of component '" + component_name + "' is not loaded")
                wavelength_grid = intrinsic_seds[component_name].wavelength_grid() # create wavelength grid from SED
                intrinsic_cubes[component_name] = DataCube.from_file(intrinsic_cube_paths[component_name], wavelength_grid=wavelength_grid)
        else: intrinsic_cubes = None

        # Create and return
        return cls(name, observed, intrinsic_seds=intrinsic_seds, intrinsic_cubes=intrinsic_cubes, distance=distance)

    # -----------------------------------------------------------------

    @property
    def component_names(self):

        """
        This function ...
        :return:
        """

        if self.has_intrinsic_seds: return self.intrinsic_seds.keys()
        elif self.has_intrinsic_cubes: return self.intrinsic_cubes.keys()
        else: raise ValueError("Component names are not defined")

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
    def nintrinsic_cubes(self):

        """
        This function ...
        :return:
        """

        return len(self.intrinsic_cubes)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_seds(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_seds is not None and self.nintrinsic_seds > 0

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cubes(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_cubes is not None and self.nintrinsic_cubes > 0

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        if self.has_transparent_sed: return True
        elif self.has_intrinsic_seds: return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        if self.has_transparent_cube: return True
        elif self.has_intrinsic_cubes: return True
        else: return False

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
        elif self.has_intrinsic_seds: return sequences.sum(self.intrinsic_seds.values())

        # Cannot be calculated
        else: raise ValueError("Intrinsic SED cannot be calculated")

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Transparent cube is written out
        if self.has_transparent_cube: return self.observed_cube_transparent

        # Has intrinsic cubes, add them
        elif self.has_intrinsic_cubes: return sequences.sum(self.intrinsic_cubes.values())

        # Cannot be calculated
        else: raise ValueError("Intrinsic datacube cannot be calculated")

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
    def faceon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Transparent cube is written out
        #if self.has_faceon_transparent_cube: return self.faceon_observed_cube_transparent

        # Cannot be calculated
        #else: raise ValueError("Intrinsic datacube from face-on orientation cannot be calculated")

        raise NotImplementedError("Not implemented")

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

    @lazyproperty
    def edgeon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Transparent cube is written out
        #if self.has_edgeon_transparent_cube: return self.edgeon_observed_cube_transparent

        # Cannot be calculated
        #else: raise ValueError("Intrinsic datacube from edge-on orientation cannot be calculated")

        raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------
