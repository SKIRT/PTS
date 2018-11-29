#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.component Contains the AnalysisComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from .context import AnalysisContext
from ...core.tools.utils import lazyproperty
from ...core.tools import filesystem as fs
from ...magic.core.frame import Frame
from ...core.tools import sequences
from ...core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

bol_map_name = "bol"
intr_stellar_map_name = "intr_stellar" # intrinsic stellar (bol) luminosity (transparent)
obs_stellar_map_name = "obs_stellar" # observed stellar (bol) luminosity
diffuse_dust_map_name = "diffuse_dust"
dust_map_name = "dust" # dust (bol) luminosity
#dust_with_internal_map_name = "dust_with_internal" # dust (bol) luminosity + internal dust (MAPPINGS)
scattered_map_name = "scattered" # scattered stellar luminosity
absorbed_diffuse_map_name = "absorbed_diffuse"
#absorbed_map_name = "absorbed" # absorbed stellar luminosity
#absorbed_with_internal_map_name = "absorbed_with_internal" # absorbed stellar luminosity + internal absorption (MAPPINGS)
fabs_diffuse_map_name = "fabs_diffuse"
fabs_map_name = "fabs"
attenuated_map_name = "attenuated" # attenuated stellar luminosity
direct_map_name = "direct" # direct stellar luminosity
stellar_mass_map_name = "stellar_mass"  # stellar mass
ssfr_map_name = "ssfr" # specific star formation rate

i1_map_name = "i1"
intr_i1_map_name = "intr_i1"

fuv_map_name = "fuv"
intr_fuv_map_name = "intr_fuv"

sfr_map_name = "sfr"
dust_mass_map_name = "dust_mass"
stellar_lum_map_name = "stellar_lum"
#dust_lum_map_name = "dust_lum"
intr_dust_map_name = "intr_dust"
diffuse_mass_map_name = "diffuse_mass"
mass_map_name = "mass"
#total_mass_map_name = "total_mass"
lum_map_name = "lum"
total_lum_map_name = "total_lum"

# -----------------------------------------------------------------

earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

total = "total"
old = "old"
young = "young"
ionizing = "ionizing"
unevolved = "unevolved"
extra = "extra"
contributions = [total, old, young, ionizing, unevolved, extra]

# -----------------------------------------------------------------

class AnalysisComponent(GalaxyModelingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(AnalysisComponent, self).__init__(*args, **kwargs)

        # The analysis context
        self.context = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisComponent, self).setup()

        # Create the analysis context
        self.context = AnalysisContext(self.analysis_path)

    # -----------------------------------------------------------------

    @property
    def timing_table_path(self):
        return self.context.timing_table_path

    # -----------------------------------------------------------------

    @property
    def memory_table_path(self):
        return self.context.memory_table_path

    # -----------------------------------------------------------------

    @property
    def timing_table(self):
        return self.context.timing_table

    # -----------------------------------------------------------------

    @property
    def memory_table(self):
        return self.context.memory_table

    # -----------------------------------------------------------------

    @property
    def cached_table(self):
        return self.context.cached_table

    # -----------------------------------------------------------------

    @property
    def analysis_run_names(self):
        return self.context.analysis_run_names

    # -----------------------------------------------------------------

    def get_run_path(self, run_name):
        return self.context.get_run_path(run_name)

    # -----------------------------------------------------------------

    def get_run_info(self, run_name):
        return self.context.get_run_info(run_name)

    # -----------------------------------------------------------------

    def get_run(self, run_name):
        return self.context.get_run(run_name)

# -----------------------------------------------------------------

class AnalysisRunComponent(AnalysisComponent):

    """
    This class ...
    """

    @lazyproperty
    def analysis_run(self):
        return self.get_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def analysis_run_path(self):
        return self.analysis_run.path

    # -----------------------------------------------------------------

    @property
    def model_name(self):
        return self.analysis_run.model_name

    # -----------------------------------------------------------------

    @property
    def model(self):
        return self.analysis_run.model

    # -----------------------------------------------------------------

    @property
    def cell_x_coordinates_colname(self):
        return "X coordinate of cell center"

    # -----------------------------------------------------------------

    @property
    def cell_y_coordinates_colname(self):
        return "Y coordinate of cell center"

    # -----------------------------------------------------------------

    @property
    def cell_z_coordinates_colname(self):
        return "Z coordinate of cell center"

    # -----------------------------------------------------------------

    @property
    def has_extra(self):
        return self.model.has_extra_component

    # -----------------------------------------------------------------
    # TOTAL SIMULATION
    # -----------------------------------------------------------------

    @property
    def total_contribution_simulation_path(self):
        return self.analysis_run.simulation_path_for_contribution(total)

    # -----------------------------------------------------------------

    @property
    def total_contribution_ski_path(self):
        return self.analysis_run.ski_path_for_contribution(total)

    # -----------------------------------------------------------------

    @property
    def total_contribution_output_path(self):
        return self.analysis_run.output_path_for_contribution(total)

    # -----------------------------------------------------------------

    @property
    def total_contribution_output(self):
        return self.model.total_simulation_output

    # -----------------------------------------------------------------

    @property
    def total_contribution_data(self):
        return self.model.total_simulation_data

    # -----------------------------------------------------------------
    # BULGE SIMULATION
    # -----------------------------------------------------------------

    @property
    def bulge_contribution_data(self):
        return self.model.bulge_simulation_data

    # -----------------------------------------------------------------
    # DISK SIMULATION
    # -----------------------------------------------------------------

    @property
    def disk_contribution_data(self):
        return self.model.disk_simulation_data

    # -----------------------------------------------------------------
    # OLD SIMULATION
    # -----------------------------------------------------------------

    @property
    def old_contribution_simulation_path(self):
        return self.analysis_run.simulation_path_for_contribution(old)

    # -----------------------------------------------------------------

    @property
    def old_contribution_ski_path(self):
        return self.analysis_run.ski_path_for_contribution(old)

    # -----------------------------------------------------------------

    @property
    def old_contribution_output_path(self):
        return self.analysis_run.output_path_for_contribution(old)

    # -----------------------------------------------------------------

    @property
    def old_contribution_output(self):
        return self.model.old_simulation_output

    # -----------------------------------------------------------------

    @property
    def old_contribution_data(self):
        return self.model.old_simulation_data

    # -----------------------------------------------------------------
    # YOUNG SIMULATION
    # -----------------------------------------------------------------

    @property
    def young_contribution_simulation_path(self):
        return self.analysis_run.simulation_path_for_contribution(young)

    # -----------------------------------------------------------------

    @property
    def young_contribution_ski_path(self):
        return self.analysis_run.ski_path_for_contribution(young)

    # -----------------------------------------------------------------

    @property
    def young_contribution_output_path(self):
        return self.analysis_run.output_path_for_contribution(young)

    # -----------------------------------------------------------------

    @property
    def young_contribution_output(self):
        return self.model.young_simulation_output

    # -----------------------------------------------------------------

    @property
    def young_contribution_data(self):
        return self.model.young_simulation_data

    # -----------------------------------------------------------------
    # IONIZING (SFR) SIMULATION
    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_simulation_path(self):
        return self.analysis_run.simulation_path_for_contribution(young)

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_ski_path(self):
        return self.analysis_run.ski_path_for_contribution(ionizing)

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_output_path(self):
        return self.analysis_run.output_path_for_contribution(ionizing)

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_output(self):
        return self.model.sfr_simulation_output

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_data(self):
        return self.model.sfr_simulation_data

    # -----------------------------------------------------------------
    # UNEVOLVED SIMULATION
    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_simulation_path(self):
        return self.analysis_run.simulation_path_for_contribution(unevolved)

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_ski_path(self):
        return self.analysis_run.ski_path_for_contribution(unevolved)

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_output_path(self):
        return self.analysis_run.output_path_for_contribution(unevolved)

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_output(self):
        return self.model.unevolved_simulation_output

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_data(self):
        return self.model.unevolved_simulation_data

    # -----------------------------------------------------------------
    # EXTRA SIMULATION
    # -----------------------------------------------------------------

    @property
    def extra_contribution_simulation_path(self):
        return self.analysis_run.simulation_path_for_contribution(extra)

    # -----------------------------------------------------------------

    @property
    def extra_contribution_ski_path(self):
        return self.analysis_run.ski_path_for_contribution(extra)

    # -----------------------------------------------------------------

    @property
    def extra_contribution_output_path(self):
        return self.analysis_run.output_path_for_contribution(extra)

    # -----------------------------------------------------------------

    @property
    def extra_contribution_output(self):
        return self.model.extra_simulation_output

    # -----------------------------------------------------------------

    @property
    def extra_contribution_data(self):
        return self.model.extra_simulation_data


    # -----------------------------------------------------------------
    # CELL PROPERTIES
    # -----------------------------------------------------------------

    @property
    def total_contribution_cell_properties_filepath(self):
        return self.total_contribution_data.cell_properties_path

    # -----------------------------------------------------------------

    @property
    def cell_properties_path(self):
        return self.total_contribution_cell_properties_filepath

    # -----------------------------------------------------------------

    @property
    def cell_properties(self):
        return self.total_contribution_data.cell_properties

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_mass_fractions(self):
        return np.asarray(self.cell_properties["Mass fraction"])

    # -----------------------------------------------------------------

    @property
    def cell_volumes(self):
        return self.model.cell_volumes  # is array

    # -----------------------------------------------------------------

    @property
    def cell_volume_unit(self):
        return self.model.cell_volume_unit

    # -----------------------------------------------------------------

    @property
    def total_contribution_absorption_filepath(self):
        #return self.total_contribution_data.absorption_path
        if self.total_contribution_data.has_absorption: return self.total_contribution_data.absorption_path
        elif self.total_contribution_data.has_isrf: return self.total_contribution_data.isrf_path
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_data(self):
        if self.total_contribution_data.has_absorption: return self.total_contribution_data.absorption
        elif self.total_contribution_data.has_isrf: return self.total_contribution_data.isrf
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_column_names(self):
        return fs.get_column_names(self.total_contribution_absorption_filepath, capitalize=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_column_units(self):
        return fs.get_column_units(self.total_contribution_absorption_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_column_name(self):
        abs_colnames = ["Absorbed bolometric luminosity", "Bolometric luminosity absorbed in cell"]
        #return sequences.find_single_in_both(abs_colnames, self.total_contribution_absorption_data.colnames)
        return sequences.find_single_in_both(abs_colnames, self.total_contribution_absorption_column_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_column_index(self):
        return self.total_contribution_absorption_column_names.index(self.total_contribution_absorption_column_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_unit(self):
        #return self.total_contribution_absorption_data.column_unit(self.total_contribution_absorption_column_name)
        return u(self.total_contribution_absorption_column_units[self.total_contribution_absorption_column_index])

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_luminosities(self):
        #return np.asarray(self.total_contribution_absorption_data[self.total_contribution_absorption_column_name])
        return fs.get_column(self.total_contribution_absorption_filepath, self.total_contribution_absorption_column_index, float, method="pandas")

    # -----------------------------------------------------------------

    @property
    def cell_coordinates_filepath(self):
        if self.total_contribution_data.has_absorption: return self.total_contribution_absorption_filepath # SKIRT7
        else: return self.cell_properties_path # SKIRT8

    # -----------------------------------------------------------------

    @property
    def cell_x_coordinates(self):
        if self.total_contribution_data.has_absorption: return np.asarray(self.total_contribution_absorption_data[self.cell_x_coordinates_colname]) # SKIRT7
        else: return self.model.cell_x_coordinates # SKIRT8

    # -----------------------------------------------------------------

    @property
    def cell_y_coordinates(self):
        if self.total_contribution_data.has_absorption: return np.asarray(self.total_contribution_absorption_data[self.cell_y_coordinates_colname]) # SKIRT7
        else: return self.model.cell_y_coordinates # SKIRT8

    # -----------------------------------------------------------------

    @property
    def cell_z_coordinates(self):
        if self.total_contribution_data.has_absorption: return np.asarray(self.total_contribution_absorption_data[self.cell_z_coordinates_colname]) # SKIRT7
        else: return self.model.cell_z_coordinates # SKIRT8

    # -----------------------------------------------------------------
    # DENSITY AND NORMALIZED MASS
    #   YOUNG
    # -----------------------------------------------------------------

    @property
    def young_cell_stellar_density(self):
        return self.model.young_cell_stellar_density  # is array

    # -----------------------------------------------------------------

    @lazyproperty
    def young_cell_normalized_mass(self):
        values = self.young_cell_stellar_density * self.cell_volumes
        values /= np.sum(values)
        return values

    # -----------------------------------------------------------------
    #   SFR
    # -----------------------------------------------------------------

    @property
    def sfr_cell_stellar_density(self):
        return self.model.sfr_cell_stellar_density  # is array

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_cell_normalized_mass(self):
        values = self.sfr_cell_stellar_density * self.cell_volumes
        values /= np.sum(values)
        return values

    # -----------------------------------------------------------------
    #   BULGE
    # -----------------------------------------------------------------

    @property
    def bulge_cell_stellar_density(self):
        return self.model.old_bulge_cell_stellar_density  # is array

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_cell_normalized_mass(self):
        values = self.bulge_cell_stellar_density * self.cell_volumes
        values /= np.sum(values)
        return values

    # -----------------------------------------------------------------
    #   DISK
    # -----------------------------------------------------------------

    @property
    def disk_cell_stellar_density(self):
        return self.model.old_disk_cell_stellar_density  # is array

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_cell_normalized_mass(self):
        values = self.disk_cell_stellar_density * self.cell_volumes
        values /= np.sum(values)
        return values

    # -----------------------------------------------------------------
    #   EXTRA
    # -----------------------------------------------------------------

    @property
    def extra_cell_stellar_density(self):
        return self.model.extra_cell_stellar_density  # is array

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_cell_normalized_mass(self):
        values = self.extra_cell_stellar_density * self.cell_volumes
        values /= np.sum(values)
        return values


    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def properties_path(self):
        return self.analysis_run.properties_path

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_path(self):
        return fs.create_directory_in(self.properties_path, "maps")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_total_path(self):
        return fs.create_directory_in(self.properties_maps_path, "total")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_total_earth_path(self):
        return fs.create_directory_in(self.properties_maps_total_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_total_faceon_path(self):
        return fs.create_directory_in(self.properties_maps_total_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_total_edgeon_path(self):
        return fs.create_directory_in(self.properties_maps_total_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_bulge_path(self):
        return fs.create_directory_in(self.properties_maps_path, "bulge")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_bulge_earth_path(self):
        return fs.create_directory_in(self.properties_maps_bulge_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_bulge_faceon_path(self):
        return fs.create_directory_in(self.properties_maps_bulge_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_bulge_edgeon_path(self):
        return fs.create_directory_in(self.properties_maps_bulge_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_disk_path(self):
        return fs.create_directory_in(self.properties_maps_path, "disk")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_disk_earth_path(self):
        return fs.create_directory_in(self.properties_maps_disk_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_disk_faceon_path(self):
        return fs.create_directory_in(self.properties_maps_disk_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_disk_edgeon_path(self):
        return fs.create_directory_in(self.properties_maps_disk_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_old_path(self):
        return fs.create_directory_in(self.properties_maps_path, "old")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_old_earth_path(self):
        return fs.create_directory_in(self.properties_maps_old_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_old_faceon_path(self):
        return fs.create_directory_in(self.properties_maps_old_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_old_edgeon_path(self):
        return fs.create_directory_in(self.properties_maps_old_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_young_path(self):
        return fs.create_directory_in(self.properties_maps_path, "young")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_young_earth_path(self):
        return fs.create_directory_in(self.properties_maps_young_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_young_faceon_path(self):
        return fs.create_directory_in(self.properties_maps_young_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_young_edgeon_path(self):
        return fs.create_directory_in(self.properties_maps_young_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_sfr_path(self):
        return fs.create_directory_in(self.properties_maps_path, "sfr")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_sfr_earth_path(self):
        return fs.create_directory_in(self.properties_maps_sfr_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_sfr_faceon_path(self):
        return fs.create_directory_in(self.properties_maps_sfr_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_sfr_edgeon_path(self):
        return fs.create_directory_in(self.properties_maps_sfr_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_unevolved_path(self):
        return fs.create_directory_in(self.properties_maps_path, "unevolved")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_unevolved_earth_path(self):
        return fs.create_directory_in(self.properties_maps_unevolved_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_unevolved_faceon_path(self):
        return fs.create_directory_in(self.properties_maps_unevolved_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_unevolved_edgeon_path(self):
        return fs.create_directory_in(self.properties_maps_unevolved_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_extra_path(self):
        return fs.create_directory_in(self.properties_maps_path, "extra")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_extra_earth_path(self):
        return fs.create_directory_in(self.properties_maps_extra_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_extra_faceon_path(self):
        return fs.create_directory_in(self.properties_maps_extra_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_extra_edgeon_path(self):
        return fs.create_directory_in(self.properties_maps_extra_path, edgeon_name)


    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_dust_path(self):
        return fs.create_directory_in(self.properties_maps_path, "dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_dust_earth_path(self):
        return fs.create_directory_in(self.properties_maps_dust_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_dust_faceon_path(self):
        return fs.create_directory_in(self.properties_maps_dust_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_dust_edgeon_path(self):
        return fs.create_directory_in(self.properties_maps_dust_path, edgeon_name)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------
    # TOTAL MAPS
    # -----------------------------------------------------------------

    def has_total_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return fs.has_file(self.properties_maps_total_earth_path, map_name, "fits")
        elif orientation == faceon_name: return fs.has_file(self.properties_maps_total_faceon_path, map_name, "fits")
        elif orientation == edgeon_name: return fs.has_file(self.properties_maps_total_edgeon_path, map_name, "fits")
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def load_total_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation: 
        :return: 
        """

        if orientation == earth_name: return Frame.from_file(fs.join(self.properties_maps_total_earth_path, map_name + ".fits"))
        elif orientation == faceon_name: return Frame.from_file(fs.join(self.properties_maps_total_faceon_path, map_name + ".fits"))
        elif orientation == edgeon_name: return Frame.from_file(fs.join(self.properties_maps_total_edgeon_path, map_name + ".fits"))
        else: raise ValueError("Invalid orientation: '" + orientation  + "'")

    # -----------------------------------------------------------------
    # TOTAL BOL
    # -----------------------------------------------------------------

    def has_total_bol_map(self, orientation=earth_name):
        return self.has_total_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_bol_map(self, orientation=earth_name):
        return self.load_total_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_bol_map(self, orientation=earth_name):
        if self.has_total_bol_map(orientation): return self.load_total_bol_map(orientation)
        if orientation == earth_name: return self.model.total_bolometric_luminosity_map_earth
        elif orientation == faceon_name: return self.model.total_bolometric_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.total_bolometric_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL INTR STELLAR
    # -----------------------------------------------------------------

    def has_total_intr_stellar_map(self, orientation=earth_name):
        return self.has_total_map(intr_stellar_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_intr_stellar_map(self, orientation=earth_name):
        return self.load_total_map(intr_stellar_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_intr_stellar_map(self, orientation=earth_name):
        if self.has_total_intr_stellar_map(orientation): return self.load_total_intr_stellar_map(orientation)
        if orientation == earth_name: return self.model.total_intrinsic_stellar_luminosity_map_earth
        elif orientation == faceon_name: return self.model.total_intrinsic_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.total_intrinsic_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL OBS STELLAR
    # -----------------------------------------------------------------

    def has_total_obs_stellar_map(self, orientation=earth_name):
        return self.has_total_map(obs_stellar_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_obs_stellar_map(self, orientation=earth_name):
        return self.load_total_map(obs_stellar_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_obs_stellar_map(self, orientation=earth_name):
        if self.has_total_obs_stellar_map(orientation): return self.load_total_obs_stellar_map(orientation)
        if orientation == earth_name: return self.model.total_observed_stellar_luminosity_map_earth
        elif orientation == faceon_name: return self.model.total_observed_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.total_observed_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL DIFFUSE DUST
    # -----------------------------------------------------------------

    def has_total_diffuse_dust_map(self, orientation=earth_name):
        return self.has_total_map(diffuse_dust_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_diffuse_dust_map(self, orientation=earth_name):
        return self.load_total_map(diffuse_dust_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_diffuse_dust_map(self, orientation=earth_name):
        if self.has_total_diffuse_dust_map(orientation): return self.load_total_diffuse_dust_map(orientation)
        if orientation == earth_name: return self.model.total_diffuse_dust_luminosity_map_earth
        elif orientation == faceon_name: return self.model.total_diffuse_dust_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.total_diffuse_dust_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL DUST
    # -----------------------------------------------------------------

    def has_total_dust_map(self, orientation=earth_name):
        return self.has_total_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_dust_map(self, orientation=earth_name):
        return self.load_total_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_dust_map(self, orientation=earth_name):
        if self.has_total_dust_map(orientation): return self.load_total_dust_map(orientation)
        if orientation == earth_name: return self.model.total_dust_luminosity_map_earth
        elif orientation == faceon_name: return self.model.total_dust_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.total_dust_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL SCATTERED
    # -----------------------------------------------------------------

    def has_total_scattered_map(self, orientation=earth_name):
        return self.has_total_map(scattered_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_scattered_map(self, orientation=earth_name):
        return self.load_total_map(scattered_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_scattered_map(self, orientation=earth_name):
        if self.has_total_scattered_map(orientation): return self.load_total_scattered_map(orientation)
        if orientation == earth_name: return self.model.total_scattered_stellar_luminosity_map_earth
        elif orientation == faceon_name: return self.model.total_scattered_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.total_scattered_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL ABSORBED DIFFUSE
    # -----------------------------------------------------------------

    def has_total_absorbed_diffuse_map(self, orientation=earth_name):
        return self.has_total_map(absorbed_diffuse_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_absorbed_diffuse_map(self, orientation=earth_name):
        return self.load_total_map(absorbed_diffuse_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_absorbed_diffuse_map(self, orientation=earth_name):
        if self.has_total_absorbed_diffuse_map(orientation): return self.load_total_absorbed_diffuse_map(orientation)
        if orientation == earth_name: return self.model.total_absorbed_diffuse_stellar_luminosity_map_earth
        elif orientation == faceon_name: return self.model.total_absorbed_diffuse_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.total_absorbed_diffuse_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL FABS DIFFUSE
    # -----------------------------------------------------------------

    def has_total_fabs_diffuse_map(self, orientation=earth_name):
        return self.has_total_map(fabs_diffuse_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_fabs_diffuse_map(self, orientation=earth_name):
        return self.load_total_map(fabs_diffuse_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_fabs_diffuse_map(self, orientation=earth_name):
        if self.has_total_fabs_diffuse_map(orientation): return self.load_total_fabs_diffuse_map(orientation)
        if orientation == earth_name: return self.model.total_fabs_diffuse_map_earth
        elif orientation == faceon_name: return self.model.total_fabs_diffuse_map_faceon
        elif orientation == edgeon_name: return self.model.total_fabs_diffuse_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL FABS
    # -----------------------------------------------------------------

    def has_total_fabs_map(self, orientation=earth_name):
        return self.has_total_map(fabs_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_fabs_map(self, orientation=earth_name):
        return self.load_total_map(fabs_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_fabs_map(self, orientation=earth_name):
        if self.has_total_fabs_map(orientation): return self.load_total_fabs_map(orientation)
        if orientation == earth_name: return self.model.total_fabs_map_earth
        elif orientation == faceon_name: return self.model.total_fabs_map_faceon
        elif orientation == edgeon_name: return self.model.total_fabs_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL ATTENUATION
    # -----------------------------------------------------------------

    def has_total_attenuated_map(self, orientation=earth_name):
        return self.has_total_map(attenuated_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_attenuated_map(self, orientation=earth_name):
        return self.load_total_map(attenuated_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_attenuated_map(self, orientation=earth_name):
        if self.has_total_attenuated_map(orientation): return self.load_total_attenuated_map(orientation)
        if orientation == earth_name: return self.model.total_attenuated_stellar_luminosity_map_earth
        elif orientation == faceon_name: return self.model.total_attenuated_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.total_attenuated_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL DIRECT
    # -----------------------------------------------------------------

    def has_total_direct_map(self, orientation=earth_name):
        return self.has_total_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_direct_map(self, orientation=earth_name):
        return self.load_total_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_direct_map(self, orientation=earth_name):
        if self.has_total_direct_map(orientation): return self.load_total_direct_map(orientation)
        if orientation == earth_name: return self.model.total_direct_stellar_luminosity_map_earth
        elif orientation == faceon_name: return self.model.total_direct_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.total_direct_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL SFR
    # -----------------------------------------------------------------

    def has_total_sfr_map(self, orientation=earth_name):
        return self.has_total_map(sfr_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_sfr_map(self, orientation=earth_name):
        return self.load_total_map(sfr_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_sfr_map(self, orientation=earth_name):
        if self.has_total_sfr_map(orientation): return self.load_total_sfr_map(orientation)
        if orientation == earth_name: return self.model.total_star_formation_rate_map_earth
        elif orientation == faceon_name: return self.model.total_star_formation_rate_map_faceon
        elif orientation == edgeon_name: return self.model.total_star_formation_rate_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL STELLAR MASS
    # -----------------------------------------------------------------

    def has_total_stellar_mass_map(self, orientation=earth_name):
        return self.has_total_map(stellar_mass_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_stellar_mass_map(self, orientation=earth_name):
        return self.load_total_map(stellar_mass_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_stellar_mass_map(self, orientation=earth_name):
        if self.has_total_stellar_mass_map(orientation): return self.load_total_stellar_mass_map(orientation)
        if orientation == earth_name: return self.model.total_stellar_mass_map_earth
        elif orientation == faceon_name: return self.model.total_stellar_mass_map_faceon
        elif orientation == edgeon_name: return self.model.total_stellar_mass_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # TOTAL SSFR
    # -----------------------------------------------------------------

    def has_total_ssfr_map(self, orientation=earth_name):
        return self.has_total_map(ssfr_map_name, orientation)

    # -----------------------------------------------------------------

    def load_total_ssfr_map(self, orientation=earth_name):
        return self.load_total_map(ssfr_map_name, orientation)

    # -----------------------------------------------------------------

    def get_total_ssfr_map(self, orientation=earth_name):
        if self.has_total_ssfr_map(orientation): return self.load_total_ssfr_map(orientation)
        if orientation == earth_name: return self.model.total_ssfr_map_earth
        elif orientation == faceon_name: return self.model.total_ssfr_map_faceon
        elif orientation == edgeon_name: return self.model.total_ssfr_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def get_total_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric luminosity
        if which == bol_map_name: return self.get_total_bol_map(orientation)

        # Intrinsic stellar luminosity (transparent luminosity)
        if which == intr_stellar_map_name: return self.get_total_intr_stellar_map(orientation)

        # Observed stellar luminosity
        elif which == obs_stellar_map_name: return self.get_total_obs_stellar_map(orientation)

        # Diffuse dust emission luminosity
        elif which == diffuse_dust_map_name: return self.get_total_diffuse_dust_map(orientation)

        # Dust emission luminosity
        elif which == dust_map_name: return self.get_total_dust_map(orientation)

        # Scattered stellar luminosity
        elif which == scattered_map_name: return self.get_total_scattered_map(orientation)

        # Absorbed stellar luminosity (by diffuse dust) (extinction)
        # absorbed = transparent - observed stellar (= observed - dust = direct + scattered)
        elif which == absorbed_diffuse_map_name: return self.get_total_absorbed_diffuse_map(orientation)

        # Absorbed stellar luminosity (extinction)
        # CUBE INFORMATION IS NOT AVAILABLE, SO MAP IS NOT USEFUL (IS JUST THE SAME AS DUST EMISSION MAP)
        #elif which == absorbed_map_name:

        # Fraction of energy absorbed by DIFFUSE dust
        elif which == fabs_diffuse_map_name: return self.get_total_fabs_diffuse_map(orientation)

        # Fraction of energy absorbed by dust
        elif which == fabs_map_name: return self.get_total_fabs_map(orientation)

        # Attenuated stellar luminosity (attenuation)
        elif which == attenuated_map_name: return self.get_total_attenuated_map(orientation) # attenuated = transparent - direct stellar

        # Direct luminosity
        elif which == direct_map_name: return self.get_total_direct_map(orientation)

        # Star formation rate
        elif which == sfr_map_name: return self.get_total_sfr_map(orientation)

        # Stellar mass
        elif which == stellar_mass_map_name: return self.get_total_stellar_mass_map(orientation)

        # Specific star formation rate
        elif which == ssfr_map_name: return self.get_total_ssfr_map(orientation)

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------
    # BULGE MAPS
    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def has_bulge_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return fs.has_file(self.properties_maps_bulge_earth_path, map_name, "fits")
        elif orientation == faceon_name: return fs.has_file(self.properties_maps_bulge_faceon_path, map_name, "fits")
        elif orientation == edgeon_name: return fs.has_file(self.properties_maps_bulge_edgeon_path, map_name, "fits")
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def load_bulge_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return Frame.from_file(fs.join(self.properties_maps_bulge_earth_path, map_name + ".fits"))
        elif orientation == faceon_name: return Frame.from_file(fs.join(self.properties_maps_bulge_faceon_path, map_name + ".fits"))
        elif orientation == edgeon_name: return Frame.from_file(fs.join(self.properties_maps_bulge_edgeon_path, map_name + ".fits"))
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # BULGE BOL
    # -----------------------------------------------------------------

    def has_bulge_bol_map(self, orientation=earth_name):
        return self.has_bulge_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def load_bulge_bol_map(self, orientation=earth_name):
        return self.load_bulge_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def get_bulge_bol_map(self, orientation=earth_name):
        if self.has_bulge_bol_map(orientation): return self.load_bulge_bol_map(orientation)
        if orientation == earth_name: return self.model.old_bulge_bolometric_luminosity_map
        elif orientation == faceon_name: return self.model.old_bulge_bolometric_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_bulge_bolometric_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # BULGE DIRECT
    # -----------------------------------------------------------------

    def has_bulge_direct_map(self, orientation=earth_name):
        return self.has_bulge_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def load_bulge_direct_map(self, orientation=earth_name):
        return self.load_bulge_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def get_bulge_direct_map(self, orientation=earth_name):
        if self.has_bulge_direct_map(orientation): return self.load_bulge_direct_map(orientation)
        if orientation == earth_name: return self.model.old_bulge_direct_stellar_luminosity_map
        elif orientation == faceon_name: return self.model.old_bulge_direct_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_bulge_direct_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation +"'")

    # -----------------------------------------------------------------
    # BULGE I1
    # -----------------------------------------------------------------

    def has_bulge_i1_map(self, orientation=earth_name):
        return self.has_bulge_map(i1_map_name, orientation)

    # -----------------------------------------------------------------

    def load_bulge_i1_map(self, orientation=earth_name):
        return self.load_bulge_map(i1_map_name, orientation)

    # -----------------------------------------------------------------

    def get_bulge_i1_map(self, orientation=earth_name):
        if self.has_bulge_i1_map(orientation): return self.load_bulge_i1_map(orientation)
        if orientation == earth_name: return self.model.old_bulge_i1_luminosity_map
        elif orientation == faceon_name: return self.model.old_bulge_i1_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_bulge_i1_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # BULGE INTR I1
    # -----------------------------------------------------------------

    def has_bulge_intr_i1_map(self, orientation=earth_name):
        return self.has_bulge_map(intr_i1_map_name, orientation)

    # -----------------------------------------------------------------

    def load_bulge_intr_i1_map(self, orientation=earth_name):
        return self.load_bulge_map(intr_i1_map_name, orientation)

    # -----------------------------------------------------------------

    def get_bulge_intr_i1_map(self, orientation=earth_name):
        if self.has_bulge_intr_i1_map(orientation): return self.load_bulge_intr_i1_map(orientation)
        if orientation == earth_name: return self.model.old_bulge_intrinsic_i1_luminosity_map
        elif orientation == faceon_name: return self.model.old_bulge_intrinsic_i1_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_bulge_intrinsic_i1_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # BULGE DUST
    # -----------------------------------------------------------------

    def has_bulge_dust_map(self, orientation=earth_name):
        return self.has_bulge_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def load_bulge_dust_map(self, orientation=earth_name):
        return self.load_bulge_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def get_bulge_dust_map(self, orientation=earth_name):
        if self.has_bulge_dust_map(orientation): return self.load_bulge_dust_map(orientation)
        if orientation == earth_name: return self.model.old_bulge_dust_luminosity_map
        elif orientation == faceon_name: return self.model.old_bulge_dust_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_bulge_dust_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def get_bulge_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric luminosity
        if which == bol_map_name: return self.get_bulge_bol_map(orientation)

        # Direct
        elif which == direct_map_name: return self.get_bulge_direct_map(orientation)

        # (observed) I1 lum
        elif which == i1_map_name: return self.get_bulge_i1_map(orientation)

        # Intrinsic I1
        elif which == intr_i1_map_name: return self.get_bulge_intr_i1_map(orientation)

        # Dust luminosity
        elif which == dust_map_name: return self.get_bulge_dust_map(orientation)

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------
    # DISK MAPS
    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def has_disk_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return fs.has_file(self.properties_maps_disk_earth_path, map_name, "fits")
        elif orientation == faceon_name: return fs.has_file(self.properties_maps_disk_faceon_path, map_name, "fits")
        elif orientation == edgeon_name: return fs.has_file(self.properties_maps_disk_edgeon_path, map_name, "fits")
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def load_disk_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return Frame.from_file(fs.join(self.properties_maps_disk_earth_path, map_name + ".fits"))
        elif orientation == faceon_name: return Frame.from_file(fs.join(self.properties_maps_disk_faceon_path, map_name + ".fits"))
        elif orientation == edgeon_name: return Frame.from_file(fs.join(self.properties_maps_disk_edgeon_path, map_name + ".fits"))
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # DISK BOL
    # -----------------------------------------------------------------

    def has_disk_bol_map(self, orientation=earth_name):
        return self.has_disk_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def load_disk_bol_map(self, orientation=earth_name):
        return self.load_disk_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def get_disk_bol_map(self, orientation=earth_name):
        if self.has_disk_bol_map(orientation): return self.load_disk_bol_map(orientation)
        if orientation == earth_name: return self.model.old_disk_bolometric_luminosity_map
        elif orientation == faceon_name: return self.model.old_disk_bolometric_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_disk_bolometric_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # DISK DIRECT
    # -----------------------------------------------------------------

    def has_disk_direct_map(self, orientation=earth_name):
        return self.has_disk_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def load_disk_direct_map(self, orientation=earth_name):
        return self.load_disk_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def get_disk_direct_map(self, orientation=earth_name):
        if self.has_disk_direct_map(orientation): return self.load_disk_direct_map(orientation)
        if orientation == earth_name: return self.model.old_disk_direct_stellar_luminosity_map
        elif orientation == faceon_name: return self.model.old_disk_direct_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_disk_direct_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # DISK I1
    # -----------------------------------------------------------------

    def has_disk_i1_map(self, orientation=earth_name):
        return self.has_disk_map(i1_map_name, orientation)

    # -----------------------------------------------------------------

    def load_disk_i1_map(self, orientation=earth_name):
        return self.load_disk_map(i1_map_name, orientation)

    # -----------------------------------------------------------------

    def get_disk_i1_map(self, orientation=earth_name):
        if self.has_disk_i1_map(orientation): return self.load_disk_i1_map(orientation)
        if orientation == earth_name: return self.model.old_disk_i1_luminosity_map
        elif orientation == faceon_name: return self.model.old_disk_i1_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_disk_i1_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # DISK INTRINSIC I1
    # -----------------------------------------------------------------

    def has_disk_intr_i1_map(self, orientation=earth_name):
        return self.has_disk_map(intr_i1_map_name, orientation)

    # -----------------------------------------------------------------

    def load_disk_intr_i1_map(self, orientation=earth_name):
        return self.load_disk_map(intr_i1_map_name, orientation)

    # -----------------------------------------------------------------

    def get_disk_intr_i1_map(self, orientation=earth_name):
        if self.has_disk_intr_i1_map(orientation): return self.load_disk_intr_i1_map(orientation)
        if orientation == earth_name: return self.model.old_disk_intrinsic_i1_luminosity_map
        elif orientation == faceon_name: return self.model.old_disk_intrinsic_i1_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_disk_intrinsic_i1_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # DISK DUST
    # -----------------------------------------------------------------

    def has_disk_dust_map(self, orientation=earth_name):
        return self.has_disk_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def load_disk_dust_map(self, orientation=earth_name):
        return self.load_disk_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def get_disk_dust_map(self, orientation=earth_name):
        if self.has_disk_dust_map(orientation): return self.load_disk_dust_map(orientation)
        if orientation == earth_name: return self.model.old_disk_dust_luminosity_map
        elif orientation == faceon_name: return self.model.old_disk_dust_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_disk_dust_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def get_disk_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric
        if which == bol_map_name: return self.get_disk_bol_map(orientation)

        # Direct
        if which == direct_map_name: return self.get_disk_direct_map(orientation)

        # (observed) I1
        elif which == i1_map_name: return self.get_disk_i1_map(orientation)

        # Intrinsic I1
        elif which == intr_i1_map_name: return self.get_disk_intr_i1_map(orientation)

        # Dust luminosity
        elif which == dust_map_name: return self.get_disk_dust_map(orientation)

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------
    # OLD MAPS
    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def has_old_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return fs.has_file(self.properties_maps_old_earth_path, map_name, "fits")
        elif orientation == faceon_name: return fs.has_file(self.properties_maps_old_faceon_path, map_name, "fits")
        elif orientation == edgeon_name: return fs.has_file(self.properties_maps_old_edgeon_path, map_name, "fits")
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def load_old_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return Frame.from_file(fs.join(self.properties_maps_old_earth_path, map_name + ".fits"))
        elif orientation == faceon_name: return Frame.from_file(fs.join(self.properties_maps_old_faceon_path, map_name + ".fits"))
        elif orientation == edgeon_name: return Frame.from_file(fs.join(self.properties_maps_old_edgeon_path, map_name + ".fits"))
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # OLD BOL
    # -----------------------------------------------------------------

    def has_old_bol_map(self, orientation=earth_name):
        return self.has_old_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def load_old_bol_map(self, orientation=earth_name):
        return self.load_old_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def get_old_bol_map(self, orientation=earth_name):
        if self.has_old_bol_map(orientation): return self.load_old_bol_map(orientation)
        if orientation == earth_name: return self.model.old_bolometric_luminosity_map
        elif orientation == faceon_name: return self.model.old_bolometric_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_bolometric_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # OLD DIRECT
    # -----------------------------------------------------------------

    def has_old_direct_map(self, orientation=earth_name):
        return self.has_old_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def load_old_direct_map(self, orientation=earth_name):
        return self.load_old_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def get_old_direct_map(self, orientation=earth_name):
        if self.has_old_direct_map(orientation): return self.load_old_direct_map(orientation)
        if orientation == earth_name: return self.model.old_direct_stellar_luminosity_map
        elif orientation == faceon_name: return self.model.old_direct_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_direct_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # OLD I1
    # -----------------------------------------------------------------

    def has_old_i1_map(self, orientation=earth_name):
        return self.has_old_map(i1_map_name, orientation)

    # -----------------------------------------------------------------

    def load_old_i1_map(self, orientation=earth_name):
        return self.load_old_map(i1_map_name, orientation)

    # -----------------------------------------------------------------

    def get_old_i1_map(self, orientation=earth_name):
        if self.has_old_i1_map(orientation): return self.load_old_i1_map(orientation)
        if orientation == earth_name: return self.model.old_i1_luminosity_map
        elif orientation == faceon_name: return self.model.old_i1_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_i1_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # OLD INTRINSIC I1
    # -----------------------------------------------------------------

    def has_old_intr_i1_map(self, orientation=earth_name):
        return self.has_old_map(intr_i1_map_name, orientation)

    # -----------------------------------------------------------------

    def load_old_intr_i1_map(self, orientation=earth_name):
        return self.load_old_map(intr_i1_map_name, orientation)

    # -----------------------------------------------------------------

    def get_old_intr_i1_map(self, orientation=earth_name):
        if self.has_old_intr_i1_map(orientation): return self.load_old_intr_i1_map(orientation)
        if orientation == earth_name: return self.model.old_intrinsic_i1_luminosity_map
        elif orientation == faceon_name: return self.model.old_intrinsic_i1_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_intrinsic_i1_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # OLD DUST
    # -----------------------------------------------------------------

    def has_old_dust_map(self, orientation=earth_name):
        return self.has_old_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def load_old_dust_map(self, orientation=earth_name):
        return self.load_old_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def get_old_dust_map(self, orientation=earth_name):
        if self.has_old_dust_map(orientation): return self.load_old_dust_map(orientation)
        if orientation == earth_name: return self.model.old_dust_luminosity_map
        elif orientation == faceon_name: return self.model.old_bulge_dust_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.old_bulge_dust_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def get_old_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric
        if which == bol_map_name: return self.get_old_bol_map(orientation)

        # Direct
        if which == direct_map_name: return self.get_old_direct_map(orientation)

        # (observed) I1
        elif which == i1_map_name: return self.get_old_i1_map(orientation)

        # Intrinsic I1
        elif which == intr_i1_map_name: return self.get_old_intr_i1_map(orientation)

        # Dust luminosity
        elif which == dust_map_name: return self.get_old_dust_map(orientation)

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------
    # YOUNG MAPS
    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def has_young_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :parma map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return fs.has_file(self.properties_maps_young_earth_path, map_name, "fits")
        elif orientation == faceon_name: return fs.has_file(self.properties_maps_young_faceon_path, map_name, "fits")
        elif orientation == edgeon_name: return fs.has_file(self.properties_maps_young_edgeon_path, map_name, "fits")
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def load_young_map(self, map_name, orientation=earth_name):

        """
        Thisf unction ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return Frame.from_file(fs.join(self.properties_maps_young_earth_path, map_name + ".fits"))
        elif orientation == faceon_name: return Frame.from_file(fs.join(self.properties_maps_young_faceon_path, map_name + ".fits"))
        elif orientation == edgeon_name: return Frame.from_file(fs.join(self.properties_maps_young_edgeon_path, map_name + ".fits"))
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # YOUNG BOL
    # -----------------------------------------------------------------

    def has_young_bol_map(self, orientation=earth_name):
        return self.has_young_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def load_young_bol_map(self, orientation=earth_name):
        return self.load_young_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def get_young_bol_map(self, orientation=earth_name):
        if self.has_young_bol_map(orientation): return self.load_young_bol_map(orientation)
        if orientation == earth_name: return self.model.young_bolometric_luminosity_map
        elif orientation == faceon_name: return self.model.young_bolometric_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.young_bolometric_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # YOUNG DIRECT
    # -----------------------------------------------------------------

    def has_young_direct_map(self, orientation=earth_name):
        return self.has_young_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def load_young_direct_map(self, orientation=earth_name):
        return self.load_young_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def get_young_direct_map(self, orientation=earth_name):
        if self.has_young_direct_map(orientation): return self.load_young_direct_map(orientation)
        if orientation == earth_name: return self.model.young_direct_stellar_luminosity_map
        elif orientation == faceon_name: return self.model.young_direct_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.young_direct_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # YOUNG FUV
    # -----------------------------------------------------------------

    def has_young_fuv_map(self, orientation=earth_name):
        return self.has_young_map(fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def load_young_fuv_map(self, orientation=earth_name):
        return self.load_young_map(fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def get_young_fuv_map(self, orientation=earth_name):
        if self.has_young_fuv_map(orientation): return self.load_young_fuv_map(orientation)
        if orientation == earth_name: return self.model.young_fuv_luminosity_map
        elif orientation == faceon_name: return self.model.young_fuv_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.young_fuv_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # YOUNG INTRINSIC FUV
    # -----------------------------------------------------------------

    def has_young_intr_fuv_map(self, orientation=earth_name):
        return self.has_young_map(intr_fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def load_young_intr_fuv_map(self, orientation=earth_name):
        return self.load_young_map(intr_fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def get_young_intr_fuv_map(self, orientation=earth_name):
        if self.has_young_intr_fuv_map(orientation): return self.load_young_intr_fuv_map(orientation)
        if orientation == earth_name: return self.model.young_intrinsic_fuv_luminosity_map
        elif orientation == faceon_name: return self.model.young_intrinsic_fuv_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.young_intrinsic_fuv_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # YOUNG DUST
    # -----------------------------------------------------------------

    def has_young_dust_map(self, orientation=earth_name):
        return self.has_young_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def load_young_dust_map(self, orientation=earth_name):
        return self.load_young_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def get_young_dust_map(self, orientation=earth_name):
        if self.has_young_dust_map(orientation): return self.load_young_dust_map(orientation)
        if orientation == earth_name: return self.model.young_dust_luminosity_map
        elif orientation == faceon_name: return self.model.young_dust_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.young_dust_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def get_young_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric
        if which == bol_map_name: return self.get_young_bol_map(orientation)

        # Direct
        elif which == direct_map_name: return self.get_young_direct_map(orientation)

        # (observed) FUV
        elif which == fuv_map_name: return self.get_young_fuv_map(orientation)

        # Intrinsic FUV
        elif which == intr_fuv_map_name: return self.get_young_intr_fuv_map(orientation)

        # Dust luminosity
        elif which == dust_map_name: return self.get_young_dust_map(orientation)

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------
    # SFR MAPS
    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def has_sfr_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return fs.has_file(self.properties_maps_sfr_earth_path, map_name, "fits")
        elif orientation == faceon_name: return fs.has_file(self.properties_maps_sfr_faceon_path, map_name, "fits")
        elif orientation == edgeon_name: return fs.has_file(self.properties_maps_sfr_edgeon_path, map_name, "fits")
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def load_sfr_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return Frame.from_file(fs.join(self.properties_maps_sfr_earth_path, map_name + ".fits"))
        elif orientation == faceon_name: return Frame.from_file(fs.join(self.properties_maps_sfr_faceon_path, map_name + ".fits"))
        elif orientation == edgeon_name: return Frame.from_file(fs.join(self.properties_maps_sfr_edgeon_path, map_name + ".fits"))
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # SFR BOL
    # -----------------------------------------------------------------

    def has_sfr_bol_map(self, orientation=earth_name):
        return self.has_sfr_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def load_sfr_bol_map(self, orientation=earth_name):
        return self.load_sfr_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def get_sfr_bol_map(self, orientation=earth_name):
        if self.has_sfr_bol_map(orientation): return self.load_sfr_bol_map(orientation)
        if orientation == earth_name: return self.model.sfr_bolometric_luminosity_map
        elif orientation == faceon_name: return self.model.sfr_bolometric_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.sfr_bolometric_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # SFR DIRECT
    # -----------------------------------------------------------------

    def has_sfr_direct_map(self, orientation=earth_name):
        return self.has_sfr_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def load_sfr_direct_map(self, orientation=earth_name):
        return self.load_sfr_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def get_sfr_direct_map(self, orientation=earth_name):
        if self.has_sfr_direct_map(orientation): return self.load_sfr_direct_map(orientation)
        if orientation == earth_name: return self.model.sfr_direct_stellar_luminosity_map
        elif orientation == faceon_name: return self.model.sfr_direct_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.sfr_direct_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # SFR FUV
    # -----------------------------------------------------------------

    def has_sfr_fuv_map(self, orientation=earth_name):
        return self.has_sfr_map(fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def load_sfr_fuv_map(self, orientation=earth_name):
        return self.load_sfr_map(fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def get_sfr_fuv_map(self, orientation=earth_name):
        if self.has_sfr_fuv_map(orientation): return self.load_sfr_fuv_map(orientation)
        if orientation == earth_name: return self.model.sfr_fuv_luminosity_map
        elif orientation == faceon_name: return self.model.sfr_fuv_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.sfr_fuv_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # SFR INTRINSIC FUV
    # -----------------------------------------------------------------

    def has_sfr_intr_fuv_map(self, orientation=earth_name):
        return self.has_sfr_map(intr_fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def load_sfr_intr_fuv_map(self, orientation=earth_name):
        return self.load_sfr_map(intr_fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def get_sfr_intr_fuv_map(self, orientation=earth_name):
        if self.has_sfr_intr_fuv_map(orientation): return self.load_sfr_intr_fuv_map(orientation)
        if orientation == earth_name: return self.model.sfr_intrinsic_fuv_luminosity_map
        elif orientation == faceon_name: return self.model.sfr_intrinsic_fuv_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.sfr_intrinsic_fuv_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # SFR SFR
    # -----------------------------------------------------------------

    def has_sfr_sfr_map(self, orientation=earth_name):
        return self.has_sfr_map(sfr_map_name, orientation)

    # -----------------------------------------------------------------

    def load_sfr_sfr_map(self, orientation=earth_name):
        return self.load_sfr_map(sfr_map_name, orientation)

    # -----------------------------------------------------------------

    def get_sfr_sfr_map(self, orientation=earth_name):
        if self.has_sfr_sfr_map(orientation): return self.load_sfr_sfr_map(orientation)
        if orientation == earth_name: return self.model.star_formation_rate_map
        elif orientation == faceon_name: return self.model.star_formation_rate_map_faceon
        elif orientation == edgeon_name: return self.model.star_formation_rate_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # SFR DUST MASS
    # -----------------------------------------------------------------

    def has_sfr_dust_mass_map(self, orientation=earth_name):
        return self.has_sfr_map(dust_mass_map_name, orientation)

    # -----------------------------------------------------------------

    def load_sfr_dust_mass_map(self, orientation=earth_name):
        return self.load_sfr_map(dust_mass_map_name, orientation)

    # -----------------------------------------------------------------

    def get_sfr_dust_mass_map(self, orientation=earth_name):
        if self.has_sfr_dust_mass_map(orientation): return self.load_sfr_dust_mass_map(orientation)
        if orientation == earth_name: return self.model.sfr_dust_mass_map
        elif orientation == faceon_name: return self.model.sfr_dust_mass_map_faceon
        elif orientation == edgeon_name: return self.model.sfr_dust_mass_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # SFR STELLAR LUM
    # -----------------------------------------------------------------

    def has_sfr_stellar_lum_map(self, orientation=earth_name):
        return self.has_sfr_map(stellar_lum_map_name, orientation)

    # -----------------------------------------------------------------

    def load_sfr_stellar_lum_map(self, orientation=earth_name):
        return self.load_sfr_map(stellar_lum_map_name, orientation)

    # -----------------------------------------------------------------

    def get_sfr_stellar_lum_map(self, orientation=earth_name):
        if self.has_sfr_stellar_lum_map(orientation): return self.load_sfr_stellar_lum_map(orientation)
        if orientation == earth_name: return self.model.sfr_stellar_luminosity_map
        elif orientation == faceon_name: return self.model.sfr_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.sfr_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # SFR INTRINSIC DUST
    # -----------------------------------------------------------------

    def has_sfr_intr_dust_map(self, orientation=earth_name):
        return self.has_sfr_map(intr_dust_map_name, orientation)

    # -----------------------------------------------------------------

    def load_sfr_intr_dust_map(self, orientation=earth_name):
        return self.load_sfr_map(intr_dust_map_name, orientation)

    # -----------------------------------------------------------------

    def get_sfr_intr_dust_map(self, orientation=earth_name):
        if self.has_sfr_intr_dust_map(orientation): return self.load_sfr_intr_dust_map(orientation)
        if orientation == earth_name: return self.model.sfr_intrinsic_dust_luminosity_map
        elif orientation == faceon_name: return self.model.sfr_intrinsic_dust_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.sfr_intrinsic_dust_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # SFR DUST
    # -----------------------------------------------------------------

    def has_sfr_dust_map(self, orientation=earth_name):
        return self.has_sfr_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def load_sfr_dust_map(self, orientation=earth_name):
        return self.load_sfr_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def get_sfr_dust_map(self, orientation=earth_name):
        if self.has_sfr_dust_map(orientation): return self.load_sfr_dust_map(orientation)
        if orientation == earth_name: return self.model.sfr_dust_luminosity_map
        elif orientation == faceon_name: return self.model.sfr_dust_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.sfr_dust_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def get_sfr_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric
        if which == bol_map_name: return self.get_sfr_bol_map(orientation)

        # Direct
        elif which == direct_map_name: return self.get_sfr_direct_map(orientation)

        # (observed) FUV
        elif which == fuv_map_name: return self.get_sfr_fuv_map(orientation)

        # Intrinsic FUV
        elif which == intr_fuv_map_name: return self.get_sfr_intr_fuv_map(orientation)

        # SFR
        elif which == sfr_map_name: return self.get_sfr_sfr_map(orientation)

        # Dust mass
        elif which == dust_mass_map_name: return self.get_sfr_dust_mass_map(orientation)

        # Stellar bolometric luminosity
        elif which == stellar_lum_map_name: return self.get_sfr_stellar_lum_map(orientation)

        # Intrinsic dust luminosity
        elif which == intr_dust_map_name: return self.get_sfr_intr_dust_map(orientation)

        # Dust bolometric luminosity
        elif which == dust_map_name: return self.get_sfr_dust_map(orientation)

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------
    # UNEVOLVED MAPS
    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def has_unevolved_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return fs.has_file(self.properties_maps_unevolved_earth_path, map_name, "fits")
        elif orientation == faceon_name: return fs.has_file(self.properties_maps_unevolved_faceon_path, map_name, "fits")
        elif orientation == edgeon_name: return fs.has_file(self.properties_maps_unevolved_edgeon_path, map_name, "fits")
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def load_unevolved_map(self, map_name, orientation=earth_name):

        """
        Thisn function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return Frame.from_file(fs.join(self.properties_maps_unevolved_earth_path, map_name + ".fits"))
        elif orientation == faceon_name: return Frame.from_file(fs.join(self.properties_maps_unevolved_faceon_path, map_name + ".fits"))
        elif orientation == edgeon_name: return Frame.from_file(fs.join(self.properties_maps_unevolved_edgeon_path, map_name + ".fits"))
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # UNEVOLVED BOL
    # -----------------------------------------------------------------

    def has_unevolved_bol_map(self, orientation=earth_name):
        return self.has_unevolved_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def load_unevolved_bol_map(self, orientation=earth_name):
        return self.load_unevolved_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def get_unevolved_bol_map(self, orientation=earth_name):
        if self.has_unevolved_bol_map(orientation): return self.load_unevolved_bol_map(orientation)
        if orientation == earth_name: return self.model.unevolved_bolometric_luminosity_map
        elif orientation == faceon_name: return self.model.unevolved_bolometric_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.unevolved_bolometric_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # UNEVOLVED DIRECT
    # -----------------------------------------------------------------

    def has_unevolved_direct_map(self, orientation=earth_name):
        return self.has_unevolved_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def load_unevolved_direct_map(self, orientation=earth_name):
        return self.load_unevolved_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def get_unevolved_direct_map(self, orientation=earth_name):
        if self.has_unevolved_direct_map(orientation): return self.load_unevolved_direct_map(orientation)
        if orientation == earth_name: return self.model.unevolved_direct_stellar_luminosity_map
        elif orientation == faceon_name: return self.model.unevolved_direct_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.unevolved_direct_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # UNEVOLVED FUV
    # -----------------------------------------------------------------

    def has_unevolved_fuv_map(self, orientation=earth_name):
        return self.has_unevolved_map(fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def load_unevolved_fuv_map(self, orientation=earth_name):
        return self.load_unevolved_map(fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def get_unevolved_fuv_map(self, orientation=earth_name):
        if self.has_unevolved_fuv_map(orientation): return self.load_unevolved_fuv_map(orientation)
        if orientation == earth_name: return self.model.unevolved_fuv_luminosity_map
        elif orientation == faceon_name: return self.model.unevolved_fuv_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.unevolved_fuv_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # UNEVOLVED INTRINSIC FUV
    # -----------------------------------------------------------------

    def has_unevolved_intr_fuv_map(self, orientation=earth_name):
        return self.has_unevolved_map(intr_fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def load_unevolved_intr_fuv_map(self, orientation=earth_name):
        return self.load_unevolved_map(intr_fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def get_unevolved_intr_fuv_map(self, orientation=earth_name):
        if self.has_unevolved_intr_fuv_map(orientation): return self.load_unevolved_intr_fuv_map(orientation)
        if orientation == earth_name: return self.model.unevolved_intrinsic_fuv_luminosity_map
        elif orientation == faceon_name: return self.model.unevolved_intrinsic_fuv_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.unevolved_intrinsic_fuv_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # UNEVOLVED SFR
    # -----------------------------------------------------------------

    def has_unevolved_sfr_map(self, orientation=earth_name):
        return self.has_unevolved_map(sfr_map_name, orientation)

    # -----------------------------------------------------------------

    def load_unevolved_sfr_map(self, orientation=earth_name):
        return self.load_unevolved_map(sfr_map_name, orientation)

    # -----------------------------------------------------------------

    def get_unevolved_sfr_map(self, orientation=earth_name):
        if self.has_unevolved_sfr_map(orientation): return self.load_unevolved_sfr_map(orientation)
        if orientation == earth_name: return self.model.unevolved_star_formation_rate_map
        elif orientation == faceon_name: return self.model.unevolved_star_formation_rate_map_faceon
        elif orientation == edgeon_name: return self.model.unevolved_star_formation_rate_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # UNEVOLVED DUST
    # -----------------------------------------------------------------

    def has_unevolved_dust_map(self, orientation=earth_name):
        return self.has_unevolved_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def load_unevolved_dust_map(self, orientation=earth_name):
        return self.load_unevolved_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def get_unevolved_dust_map(self, orientation=earth_name):
        if self.has_unevolved_dust_map(orientation): return self.load_unevolved_dust_map(orientation)
        if orientation == earth_name: return self.model.unevolved_dust_luminosity_map
        elif orientation == faceon_name: return self.model.unevolved_dust_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.unevolved_dust_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def get_unevolved_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric
        if which == bol_map_name: return self.get_unevolved_bol_map(orientation)

        # Direct
        elif which == direct_map_name: return self.get_unevolved_direct_map(orientation)

        # FUV
        elif which == fuv_map_name: return self.get_unevolved_fuv_map(orientation)

        # Intrinsic FUV
        elif which == intr_fuv_map_name: return self.get_unevolved_intr_fuv_map(orientation)

        # SFR
        elif which == sfr_map_name: return self.get_unevolved_sfr_map(orientation)

        # Dust luminosity
        elif which == dust_map_name: return self.get_unevolved_dust_map(orientation)

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------
    # EXTRA COMPONENT MAPS
    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def has_extra_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :parma map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return fs.has_file(self.properties_maps_extra_earth_path, map_name, "fits")
        elif orientation == faceon_name: return fs.has_file(self.properties_maps_extra_faceon_path, map_name, "fits")
        elif orientation == edgeon_name: return fs.has_file(self.properties_maps_extra_edgeon_path, map_name, "fits")
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def load_extra_map(self, map_name, orientation=earth_name):

        """
        Thisf unction ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return Frame.from_file(fs.join(self.properties_maps_extra_earth_path, map_name + ".fits"))
        elif orientation == faceon_name: return Frame.from_file(fs.join(self.properties_maps_extra_faceon_path, map_name + ".fits"))
        elif orientation == edgeon_name: return Frame.from_file(fs.join(self.properties_maps_extra_edgeon_path, map_name + ".fits"))
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # EXTRA BOL
    # -----------------------------------------------------------------

    def has_extra_bol_map(self, orientation=earth_name):
        return self.has_extra_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def load_extra_bol_map(self, orientation=earth_name):
        return self.load_extra_map(bol_map_name, orientation)

    # -----------------------------------------------------------------

    def get_extra_bol_map(self, orientation=earth_name):
        if self.has_extra_bol_map(orientation): return self.load_extra_bol_map(orientation)
        if orientation == earth_name: return self.model.extra_bolometric_luminosity_map
        elif orientation == faceon_name: return self.model.extra_bolometric_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.extra_bolometric_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # EXTRA DIRECT
    # -----------------------------------------------------------------

    def has_extra_direct_map(self, orientation=earth_name):
        return self.has_extra_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def load_extra_direct_map(self, orientation=earth_name):
        return self.load_extra_map(direct_map_name, orientation)

    # -----------------------------------------------------------------

    def get_extra_direct_map(self, orientation=earth_name):
        if self.has_extra_direct_map(orientation): return self.load_extra_direct_map(orientation)
        if orientation == earth_name: return self.model.extra_direct_stellar_luminosity_map
        elif orientation == faceon_name: return self.model.extra_direct_stellar_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.extra_direct_stellar_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # EXTRA FUV
    # -----------------------------------------------------------------

    def has_extra_fuv_map(self, orientation=earth_name):
        return self.has_extra_map(fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def load_extra_fuv_map(self, orientation=earth_name):
        return self.load_extra_map(fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def get_extra_fuv_map(self, orientation=earth_name):
        if self.has_extra_fuv_map(orientation): return self.load_extra_fuv_map(orientation)
        if orientation == earth_name: return self.model.extra_fuv_luminosity_map
        elif orientation == faceon_name: return self.model.extra_fuv_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.extra_fuv_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # EXTRA INTRINSIC FUV
    # -----------------------------------------------------------------

    def has_extra_intr_fuv_map(self, orientation=earth_name):
        return self.has_extra_map(intr_fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def load_extra_intr_fuv_map(self, orientation=earth_name):
        return self.load_extra_map(intr_fuv_map_name, orientation)

    # -----------------------------------------------------------------

    def get_extra_intr_fuv_map(self, orientation=earth_name):
        if self.has_extra_intr_fuv_map(orientation): return self.load_extra_intr_fuv_map(orientation)
        if orientation == earth_name: return self.model.extra_intrinsic_fuv_luminosity_map
        elif orientation == faceon_name: return self.model.extra_intrinsic_fuv_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.extra_intrinsic_fuv_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # EXTRA DUST
    # -----------------------------------------------------------------

    def has_extra_dust_map(self, orientation=earth_name):
        return self.has_extra_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def load_extra_dust_map(self, orientation=earth_name):
        return self.load_extra_map(dust_map_name, orientation)

    # -----------------------------------------------------------------

    def get_extra_dust_map(self, orientation=earth_name):
        if self.has_extra_dust_map(orientation): return self.load_extra_dust_map(orientation)
        if orientation == earth_name: return self.model.extra_dust_luminosity_map
        elif orientation == faceon_name: return self.model.extra_dust_luminosity_map_faceon
        elif orientation == edgeon_name: return self.model.extra_dust_luminosity_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def get_extra_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric
        if which == bol_map_name: return self.get_extra_bol_map(orientation)

        # Direct
        elif which == direct_map_name: return self.get_extra_direct_map(orientation)

        # (observed) FUV
        elif which == fuv_map_name: return self.get_extra_fuv_map(orientation)

        # Intrinsic FUV
        elif which == intr_fuv_map_name: return self.get_extra_intr_fuv_map(orientation)

        # Dust luminosity
        elif which == dust_map_name: return self.get_extra_dust_map(orientation)

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")


    # -----------------------------------------------------------------
    # DUST MAPS
    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def has_dust_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return fs.has_file(self.properties_maps_dust_earth_path, map_name, "fits")
        elif orientation == faceon_name: return fs.has_file(self.properties_maps_dust_faceon_path, map_name, "fits")
        elif orientation == edgeon_name: return fs.has_file(self.properties_maps_dust_edgeon_path, map_name, "fits")
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def load_dust_map(self, map_name, orientation=earth_name):

        """
        This function ...
        :param map_name:
        :param orientation:
        :return:
        """

        if orientation == earth_name: return Frame.from_file(fs.join(self.properties_maps_dust_earth_path, map_name + ".fits"))
        elif orientation == faceon_name: return Frame.from_file(fs.join(self.properties_maps_dust_faceon_path, map_name + ".fits"))
        elif orientation == edgeon_name: return Frame.from_file(fs.join(self.properties_maps_dust_edgeon_path, map_name + ".fits"))
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # DUST DIFFUSE MASS
    # -----------------------------------------------------------------

    def has_dust_diffuse_mass_map(self, orientation=earth_name):
        return self.has_dust_map(diffuse_mass_map_name, orientation)

    # -----------------------------------------------------------------

    def load_dust_diffuse_mass_map(self, orientation=earth_name):
        return self.load_dust_map(diffuse_mass_map_name, orientation)

    # -----------------------------------------------------------------

    def get_dust_diffuse_mass_map(self, orientation=earth_name):
        if self.has_dust_diffuse_mass_map(orientation): return self.load_dust_diffuse_mass_map(orientation)
        if orientation == earth_name: return self.model.diffuse_dust_mass_map
        elif orientation == faceon_name: return self.model.diffuse_dust_mass_map_faceon
        elif orientation == edgeon_name: return self.model.diffuse_dust_mass_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------
    # DUST MASS
    # -----------------------------------------------------------------

    def has_dust_mass_map(self, orientation=earth_name):
        return self.has_dust_map(mass_map_name, orientation)

    # -----------------------------------------------------------------

    def load_dust_mass_map(self, orientation=earth_name):
        return self.load_dust_map(mass_map_name, orientation)

    # -----------------------------------------------------------------

    def get_dust_mass_map(self, orientation=earth_name):
        if self.has_dust_mass_map(orientation): return self.load_dust_mass_map(orientation)
        if orientation == earth_name: return self.model.dust_mass_map
        elif orientation == faceon_name: return self.model.dust_mass_map_faceon
        elif orientation == edgeon_name: return self.model.dust_mass_map_edgeon
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # -----------------------------------------------------------------

    def get_dust_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Dust mass
        if which == diffuse_mass_map_name: return self.get_dust_diffuse_mass_map(orientation)

        # Total dust mass
        elif which == mass_map_name: return self.get_dust_mass_map(orientation)

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

# -----------------------------------------------------------------
