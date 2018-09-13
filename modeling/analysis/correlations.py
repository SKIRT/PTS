#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.correlations Contains the CorrelationsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import AnalysisRunComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.tools import introspection
from ...core.tools.utils import lazyproperty, lazyfileproperty
from ...core.basics.scatter import Scatter2D
from ..core.data import Data3D
from ...magic.core.frame import Frame
from ...magic.tools.plotting import plot_scatters, plot_stilts
from ...core.units.parsing import parse_unit as u
from ...core.tools import sequences

# -----------------------------------------------------------------

# Define paths to reference sSFR to Funev data files
correlations_data_path = fs.join(introspection.pts_dat_dir("modeling"), "HeatingCorrelations")
m51_data_path = fs.join(correlations_data_path, "DeLooze2014sSFRHeating.dat")
m31_data_path = fs.join(correlations_data_path, "M31_sSFR_Fyoung.dat")
m51_name = "M51"
m31_name = "M31"

# -----------------------------------------------------------------

# Names for correlations
ssfr_funev_name = "sSFR-Funev"
temperature_funev_name = "Temperature-Funev"

# -----------------------------------------------------------------

# Auxilary column names
#sfr_name = "SFR"
sfr_density_name = "vSFR"
#dust_mass_name = "Dust mass"
dust_density_name = "Dust density"
distance_center_name = "Distance from center"
bulge_disk_ratio_name = "Bulge disk ratio"
#fuv_ratio_name = "FUV ratio"
#fuv_h_name =
#fuv_i1_name =
temperature_name = "Dust temperature"
mean_age_name = "Mean stellar age"

# Auxilary column names for sSFR-Funev scatter data
#aux_colnames = [sfr_name, dust_mass_name, distance_center_name, bulge_disk_ratio_name]
aux_colnames = [sfr_density_name, dust_density_name, distance_center_name, bulge_disk_ratio_name]

# -----------------------------------------------------------------

ssfr_funev_base_colnames = ["sSFR", "Funev"]

# -----------------------------------------------------------------

class CorrelationsAnalyser(AnalysisRunComponent):
    
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
        super(CorrelationsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Writing
        self.write()

        # Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(CorrelationsAnalyser, self).setup()

        # Replot?
        if self.config.replot:
            self.config.replot_ssfr_funev = True

        # Replot sSFR-Funev
        if self.config.replot_ssfr_funev:
            self.config.replot_ssfr_salim_funev = True
            self.config.replot_ssfr_ke_funev = True
            self.config.replot_ssfr_mappings_ke_funev = True

        # Salim
        if self.config.replot_ssfr_salim_funev:
            fs.remove_file(self.ssfr_salim_funev_cells_plot_path)
            fs.remove_file(self.ssfr_salim_funev_pixels_plot_path)

        # KE
        if self.config.replot_ssfr_ke_funev:
            fs.remove_file(self.ssfr_ke_funev_cells_plot_path)
            fs.remove_file(self.ssfr_ke_funev_pixels_plot_path)

        # MAPPINGS + KE
        if self.config.replot_ssfr_mappings_ke_funev:
            fs.remove_file(self.ssfr_mappings_ke_funev_cells_plot_path)
            fs.remove_file(self.ssfr_mappings_ke_funev_pixels_plot_path)

    # -----------------------------------------------------------------

    @property
    def correlations_path(self):
        return self.analysis_run.correlations_path

    # -----------------------------------------------------------------

    @property
    def heating_path(self):
        return self.analysis_run.heating_path

    # -----------------------------------------------------------------

    @property
    def projected_heating_path(self):
        return fs.join(self.heating_path, "projected")

    # -----------------------------------------------------------------

    @property
    def cell_heating_path(self):
        return fs.join(self.heating_path, "cell")

    # -----------------------------------------------------------------

    @property
    def sfr_path(self):
        return self.analysis_run.sfr_path

    # -----------------------------------------------------------------

    @property
    def projected_sfr_path(self):
        return fs.join(self.sfr_path, "projected")

    # -----------------------------------------------------------------

    @property
    def cell_sfr_path(self):
        return fs.join(self.sfr_path, "cell")

    # -----------------------------------------------------------------
    # CELL PROPERTIES
    # -----------------------------------------------------------------

    @lazyproperty
    def cell_radii(self):
        return np.sqrt(self.cell_x_coordinates**2 + self.cell_y_coordinates**2 + self.cell_z_coordinates**2)

    # -----------------------------------------------------------------

    @lazyproperty
    def length_unit(self):
        return u("pc")

    # -----------------------------------------------------------------

    @lazyproperty
    def volume_unit(self):
        return self.length_unit**3

    # -----------------------------------------------------------------

    @property
    def cell_temperatures(self):
        return self.model.cell_temperatures

    # -----------------------------------------------------------------

    @property
    def temperature_unit(self):
        return self.model.cell_temperature_unit

    # -----------------------------------------------------------------
    # I1 LUMINOSITY / BULGE DISK RATIO
    # -----------------------------------------------------------------

    @lazyproperty
    def i1_luminosity_unit(self):
        return u("W/micron")

    # -----------------------------------------------------------------

    @property
    def bulge_intrinsic_i1_luminosity(self):
        return self.model.intrinsic_i1_luminosity_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_intrinsic_i1_luminosity_scalar(self):
        return self.bulge_intrinsic_i1_luminosity.to(self.i1_luminosity_unit).value

    # -----------------------------------------------------------------

    @property
    def disk_intrinsic_i1_luminosity(self):
        return self.model.intrinsic_i1_luminosity_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_intrinsic_i1_luminosity_scalar(self):
        return self.disk_intrinsic_i1_luminosity.to(self.i1_luminosity_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_cell_i1_luminosities(self):
        return self.bulge_cell_normalized_mass * self.bulge_intrinsic_i1_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_cell_i1_luminosities(self):
        return self.disk_cell_normalized_mass * self.disk_intrinsic_i1_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def old_cell_i1_luminosities(self):
        return self.bulge_cell_i1_luminosities + self.disk_cell_i1_luminosities

    # -----------------------------------------------------------------

    @property
    def cell_bd_ratio_path(self):
        return fs.join(self.correlations_path, "bd_ratio.dat")

    # -----------------------------------------------------------------

    @property
    def bd_ratio_name(self):
        return "BD_I1"

    # -----------------------------------------------------------------

    @property
    def bd_ratio_description(self):
        return "Ratio of bulge and disk I1 luminosity"

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_bd_ratio_path", True, write=True)
    def cell_bd_ratio_data(self):

        """
        This function ...
        :return:
        """

        # Calculate ratios in cells
        ratios = self.bulge_cell_i1_luminosities / self.disk_cell_i1_luminosities

        # Create the data with external xyz
        return Data3D.from_values(self.bd_ratio_name, ratios, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.bd_ratio_description,
                                  xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------

    @property
    def cell_bd_ratios(self):
        return self.cell_bd_ratio_data.values

    # -----------------------------------------------------------------
    # DUST MASS
    # -----------------------------------------------------------------

    @property
    def dust_mass_unit(self):
        return u("Msun")

    # -----------------------------------------------------------------

    @property
    def cell_diffuse_dust_mass_fractions(self):
        return self.model.cell_mass_fractions

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_dust_mass_scalar(self):
        return self.model.diffuse_dust_mass.to(self.dust_mass_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_diffuse_dust_masses(self):
        return self.cell_diffuse_dust_mass_fractions * self.diffuse_dust_mass_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_mass_scalar(self):
        return self.model.sfr_dust_mass.to(self.dust_mass_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_dust_masses(self):
        return self.sfr_cell_normalized_mass * self.sfr_dust_mass_scalar

    # -----------------------------------------------------------------

    @property
    def cell_dust_mass_path(self):
        return fs.join(self.correlations_path, "dust_mass.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_dust_mass(self):
        return fs.is_file(self.cell_dust_mass_path)

    # -----------------------------------------------------------------

    @property
    def dust_mass_name(self):
        return "Mdust"

    # -----------------------------------------------------------------

    @property
    def dust_mass_description(self):
        return "Total dust mass"

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_dust_mass_path", True, write=True)
    def cell_dust_mass_data(self):

        """
        This function ...
        :return:
        """

        # Calculate total dust masses in cells
        dust_masses = self.cell_diffuse_dust_masses + self.cell_sfr_dust_masses

        # Create the data with external xyz
        return Data3D.from_values(self.dust_mass_name, dust_masses, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, unit=self.dust_mass_unit,
                                  description=self.dust_mass_description, xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------

    @property
    def cell_dust_mass_values(self):
        return self.cell_dust_mass_data.values

    # -----------------------------------------------------------------

    @property
    def cell_dust_mass_unit(self):
        return self.cell_dust_mass_data.unit

    # -----------------------------------------------------------------

    @property
    def cell_dust_density_path(self):
        return fs.join(self.correlations_path, "dust_density.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_dust_density(self):
        return fs.is_file(self.cell_dust_density_path)

    # -----------------------------------------------------------------

    @property
    def dust_density_name(self):
        return "Rho_dust"

    # -----------------------------------------------------------------

    @property
    def dust_density_description(self):
        return "Total dust density"

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_dust_density_path", True, write=True)
    def cell_dust_density_data(self):

        """
        This function ...
        :return:
        """

        # Get values and unit
        densities = self.cell_dust_mass_values / self.cell_volumes
        unit = self.cell_dust_mass_unit / self.volume_unit

        # Create the data with external xyz
        return Data3D.from_values(self.dust_density_name, densities, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, unit=unit, description=self.dust_density_description,
                                  xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------

    @property
    def cell_dust_densities(self):
        return self.cell_dust_density_data.values

    # -----------------------------------------------------------------

    @property
    def cell_dust_density_unit(self):
        return self.cell_dust_density_data.unit

    # -----------------------------------------------------------------
    # M51 DATA
    # -----------------------------------------------------------------

    @lazyproperty
    def m51_column_names(self):
        return fs.get_column_names(m51_data_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def m51_data(self):
        return np.loadtxt(m51_data_path, unpack=True)

    # -----------------------------------------------------------------

    @property
    def m51_log10_ssfr(self):
        return self.m51_data[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def m51_ssfr(self):
        return 10**self.m51_log10_ssfr

    # -----------------------------------------------------------------

    @lazyproperty
    def m51_fractions(self):
        return 10**self.m51_data[1]

    # -----------------------------------------------------------------

    @property
    def m51_ssfr_funev_path(self):
        return fs.join(self.ssfr_funev_path, "m51.dat")

    # -----------------------------------------------------------------

    @property
    def has_m51_ssfr_funev(self):
        return fs.is_file(self.m51_ssfr_funev_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_unit(self):
        return u("Msun/yr")

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "m51_ssfr_funev_path", True, write=False)
    def m51_scatter(self):

        """
        This function ...
        :return:
        """

        return Scatter2D.from_xy(self.m51_ssfr, self.m51_fractions, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_unit,
                                 x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    # M31 data
    # -----------------------------------------------------------------

    @lazyproperty
    def m31_column_names(self):
        return fs.get_column_names(m31_data_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def m31_data(self):
        return np.loadtxt(m31_data_path, unpack=True)

    # -----------------------------------------------------------------

    @property
    def m31_log10_ssfr(self):
        return self.m31_data[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def m31_ssfr(self):
        return 10**self.m31_log10_ssfr

    # -----------------------------------------------------------------

    @lazyproperty
    def m31_fractions(self):
        return self.m31_data[1] / 100.

    # -----------------------------------------------------------------

    @property
    def m31_ssfr_funev_path(self):
        return fs.join(self.ssfr_funev_path, "m31.dat")

    # -----------------------------------------------------------------

    @property
    def has_m31_ssfr_funev(self):
        return fs.is_file(self.m31_ssfr_funev_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "m31_ssfr_funev_path", True, write=False)
    def m31_scatter(self):

        """
        This function ...
        :return:
        """

        return Scatter2D.from_xy(self.m31_ssfr, self.m31_fractions, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_unit,
                                 x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    # NAMES AND DESCRIPTIONS
    # -----------------------------------------------------------------

    @property
    def ssfr_name(self):
        return "sSFR"

    # -----------------------------------------------------------------

    @property
    def ssfr_description(self):
        return "specific star formation rate"

    # -----------------------------------------------------------------

    @property
    def temperature_name(self):
        return "Tdust"

    # -----------------------------------------------------------------

    @property
    def temperature_description(self):
        return "dust temperature"

    # -----------------------------------------------------------------

    @property
    def funev_name(self):
        return "Funev"

    # -----------------------------------------------------------------

    @property
    def funev_description(self):
        return "fraction of dust heating by unevolved stars"

    # -----------------------------------------------------------------
    # sSFR-Funev
    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_funev_path(self):
        return fs.create_directory_in(self.correlations_path, ssfr_funev_name)

    # -----------------------------------------------------------------
    # sSFR-Funev cell scatter data
    #   SALIM
    # -----------------------------------------------------------------

    @property
    def cell_funev_path(self):
        return fs.join(self.cell_heating_path, "total_fractions.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_funev(self):
        return fs.is_file(self.cell_funev_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_salim_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_salim.dat")

    # -----------------------------------------------------------------

    @property
    def cell_sfr_salim_path(self):
        return fs.join(self.cell_sfr_path, "sfr_salim.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_salim(self):
        return fs.is_file(self.cell_ssfr_salim_path)

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_salim(self):
        return fs.is_file(self.cell_sfr_salim_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_funev(self):
        return Data3D.from_file(self.cell_funev_path)

    # -----------------------------------------------------------------

    @property
    def cell_funev_values(self):
        return self.cell_funev.values

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_salim(self):
        return Data3D.from_file(self.cell_ssfr_salim_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_salim(self):
        return Data3D.from_file(self.cell_sfr_salim_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_salim_values(self):
        return self.cell_ssfr_salim.values

    # -----------------------------------------------------------------

    @property
    def cell_sfr_salim_values(self):
        return self.cell_sfr_salim.values

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_densities_salim(self):
        return self.cell_sfr_salim_values / self.cell_volumes

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_unit(self):
        return self.cell_ssfr_salim.unit

    # -----------------------------------------------------------------

    @property
    def sfr_salim_unit(self):
        return self.cell_sfr_salim.unit

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_density_salim_unit(self):
        return self.sfr_salim_unit / self.volume_unit

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_salim_mask(self):
        return np.isfinite(self.cell_ssfr_salim_values) * (self.cell_ssfr_salim_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_funev_mask(self):
        return np.isfinite(self.cell_funev_values) * (self.cell_funev_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mask_salim(self):
        return self.valid_cell_ssfr_salim_mask * self.valid_cell_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_values_salim(self):
        return self.cell_ssfr_salim_values[self.valid_cell_mask_salim]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_funev_values_salim(self):
        return self.cell_funev_values[self.valid_cell_mask_salim]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_sfr_values_salim(self):
        return self.cell_sfr_salim_values[self.valid_cell_mask_salim]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_sfr_densities_salim(self):
        return self.cell_sfr_densities_salim[self.valid_cell_mask_salim]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_dust_mass_values_salim(self):
        return self.cell_dust_mass_values[self.valid_cell_mask_salim]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_dust_densities_salim(self):
        return self.cell_dust_densities[self.valid_cell_mask_salim]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_radii_salim(self):
        return self.cell_radii[self.valid_cell_mask_salim]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_bd_ratios_salim(self):
        return self.cell_bd_ratios[self.valid_cell_mask_salim]

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells_salim.dat")

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_colnames(self):
        return fs.get_column_names(self.ssfr_salim_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.ssfr_salim_funev_cells_colnames, ssfr_funev_base_colnames)

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.ssfr_salim_funev_cells_aux_colnames, aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_ssfr_salim_funev_cells(self):
        if fs.is_file(self.ssfr_salim_funev_cells_path):
            if not self.ssfr_salim_funev_cells_has_all_aux_columns:
                colnames = self.ssfr_salim_funev_cells_aux_colnames
                #if sfr_name not in colnames: self.ssfr_salim_funev_cells.add_aux(sfr_name, self.valid_cell_sfr_values_salim, self.sfr_salim_unit, as_column=True)
                #if dust_mass_name not in colnames: self.ssfr_salim_funev_cells.add_aux(dust_mass_name, self.valid_cell_dust_mass_values_salim, self.cell_dust_mass_unit, as_column=True)
                if sfr_density_name not in colnames: self.ssfr_salim_funev_cells.add_aux(sfr_density_name, self.valid_cell_sfr_densities_salim, self.sfr_density_salim_unit, as_column=True)
                if dust_density_name not in colnames: self.ssfr_salim_funev_cells.add_aux(dust_density_name, self.valid_cell_dust_densities_salim, self.cell_dust_density_unit, as_column=True)
                if distance_center_name not in colnames: self.ssfr_salim_funev_cells.add_aux(distance_center_name, self.valid_cell_radii_salim, self.length_unit, as_column=True)
                if bulge_disk_ratio_name not in colnames: self.ssfr_salim_funev_cells.add_aux(bulge_disk_ratio_name, self.valid_cell_bd_ratios_salim, as_column=True)
                self.ssfr_salim_funev_cells.save() # save
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_salim_cells_aux(self):
        return {#sfr_name: self.valid_cell_sfr_values_salim,
                #dust_mass_name: self.valid_cell_dust_mass_values_salim,
                sfr_density_name: self.valid_cell_sfr_densities_salim,
                dust_density_name: self.valid_cell_dust_densities_salim,
                distance_center_name: self.valid_cell_radii_salim,
                bulge_disk_ratio_name: self.valid_cell_bd_ratios_salim}

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_salim_cells_aux_units(self):
        return {#sfr_name: self.sfr_salim_unit,
                #dust_mass_name: self.cell_dust_mass_unit,
                sfr_density_name: self.sfr_density_salim_unit,
                dust_density_name: self.cell_dust_density_unit,
                distance_center_name: self.length_unit}

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_salim_funev_cells_path", True, write=False)
    def ssfr_salim_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_ssfr_salim: raise IOError("The cell sSFR (Salim) data is not present: run the SFR analysis first")
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Create and return
        return Scatter2D.from_xy(self.valid_cell_ssfr_values_salim, self.valid_cell_funev_values_salim,
                                 x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_salim_unit,
                                 x_description=self.ssfr_description, y_description=self.funev_description,
                                 aux=self.ssfr_salim_cells_aux, aux_units=self.ssfr_salim_cells_aux_units) # auxilary axes

    # -----------------------------------------------------------------
    #   K&E
    # -----------------------------------------------------------------

    @property
    def cell_ssfr_ke_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_ke.dat")

    # -----------------------------------------------------------------

    @property
    def cell_sfr_ke_path(self):
        return fs.join(self.cell_sfr_path, "sfr_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_ke(self):
        return fs.is_file(self.cell_ssfr_ke_path)

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_ke(self):
        return fs.is_file(self.cell_sfr_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_ke(self):
        return Data3D.from_file(self.cell_ssfr_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_ke(self):
        return Data3D.from_file(self.cell_sfr_ke_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_unit(self):
        return self.cell_ssfr_ke.unit

    # -----------------------------------------------------------------

    @property
    def sfr_ke_unit(self):
        return self.cell_sfr_ke.unit

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_density_ke_unit(self):
        return self.sfr_ke_unit / self.volume_unit

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_ke_values(self):
        return self.cell_ssfr_ke.values

    # -----------------------------------------------------------------

    @property
    def cell_sfr_ke_values(self):
        return self.cell_sfr_ke.values

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_densities_ke(self):
        return self.cell_sfr_ke_values / self.cell_volumes

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_ke_mask(self):
        return np.isfinite(self.cell_ssfr_ke_values) * (self.cell_ssfr_ke_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_funev_mask(self):
        return np.isfinite(self.cell_funev_values) * (self.cell_funev_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mask_ke(self):
        return self.valid_cell_ssfr_ke_mask * self.valid_cell_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_values_ke(self):
        return self.cell_ssfr_ke_values[self.valid_cell_mask_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_funev_values_ke(self):
        return self.cell_funev_values[self.valid_cell_mask_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_sfr_values_ke(self):
        return self.cell_sfr_ke_values[self.valid_cell_mask_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_sfr_densities_ke(self):
        return self.cell_sfr_densities_ke[self.valid_cell_mask_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_dust_mass_values_ke(self):
        return self.cell_dust_mass_values[self.valid_cell_mask_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_dust_densities_ke(self):
        return self.cell_dust_densities[self.valid_cell_mask_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_radii_ke(self):
        return self.cell_radii[self.valid_cell_mask_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_bd_ratios_ke(self):
        return self.cell_bd_ratios[self.valid_cell_mask_ke]

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells_ke.dat")

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_colnames(self):
        return fs.get_column_names(self.ssfr_ke_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.ssfr_ke_funev_cells_colnames, ssfr_funev_base_colnames)

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.ssfr_ke_funev_cells_aux_colnames, aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_ssfr_ke_funev_cells(self):
        if fs.is_file(self.ssfr_ke_funev_cells_path):
            if not self.ssfr_ke_funev_cells_has_all_aux_columns:
                colnames = self.ssfr_ke_funev_cells_aux_colnames
                #if sfr_name not in colnames: self.ssfr_ke_funev_cells.add_aux(sfr_name, self.valid_cell_sfr_values_ke, self.sfr_ke_unit, as_column=True)
                #if dust_mass_name not in colnames: self.ssfr_ke_funev_cells.add_aux(dust_mass_name, self.valid_cell_dust_mass_values_ke, self.cell_dust_mass_unit, as_column=True)
                if sfr_density_name not in colnames: self.ssfr_ke_funev_cells.add_aux(sfr_density_name, self.valid_cell_sfr_densities_ke, self.sfr_density_ke_unit, as_column=True)
                if dust_density_name not in colnames: self.ssfr_ke_funev_cells.add_aux(dust_density_name, self.valid_cell_dust_densities_ke, self.cell_dust_density_unit, as_column=True)
                if distance_center_name not in colnames: self.ssfr_ke_funev_cells.add_aux(distance_center_name, self.valid_cell_radii_ke, self.length_unit, as_column=True)
                if bulge_disk_ratio_name not in colnames: self.ssfr_ke_funev_cells.add_aux(bulge_disk_ratio_name, self.valid_cell_bd_ratios_ke, as_column=True)
                self.ssfr_ke_funev_cells.save() # save
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_ke_cells_aux(self):
        return {#sfr_name: self.valid_cell_sfr_values_ke,
                #dust_mass_name: self.valid_cell_dust_mass_values_ke,
                sfr_density_name: self.valid_cell_sfr_densities_ke,
                dust_density_name: self.valid_cell_dust_densities_ke,
                distance_center_name: self.valid_cell_radii_ke,
                bulge_disk_ratio_name: self.valid_cell_bd_ratios_ke}

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_ke_cells_aux_units(self):
        return {#sfr_name: self.sfr_ke_unit,
                #dust_mass_name: self.cell_dust_mass_unit,
                sfr_density_name: self.sfr_density_ke_unit,
                dust_density_name: self.cell_dust_density_unit,
                distance_center_name: self.length_unit}

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_ke_funev_cells_path", True, write=False)
    def ssfr_ke_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_ssfr_ke: raise IOError("The cell sSFR (K&E) data is not present: run the SFR analysis first")
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Create and return
        return Scatter2D.from_xy(self.valid_cell_ssfr_values_ke, self.valid_cell_funev_values_ke,
                                 x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_ke_unit,
                                 x_description=self.ssfr_description, y_description=self.funev_description,
                                 aux=self.ssfr_ke_cells_aux, aux_units=self.ssfr_ke_cells_aux_units)  # auxilary axes

    # -----------------------------------------------------------------
    #   MAPPINGS
    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_mappings.dat")

    # -----------------------------------------------------------------

    @property
    def cell_sfr_mappings_path(self):
        return fs.join(self.cell_sfr_path, "sfr_mappings.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_mappings(self):
        return fs.is_file(self.cell_ssfr_mappings_path)

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_mappings(self):
        return fs.is_file(self.cell_sfr_mappings_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_mappings(self):
        return Data3D.from_file(self.cell_ssfr_mappings_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_mappings(self):
        return Data3D.from_file(self.cell_sfr_mappings_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_unit(self):
        return self.cell_ssfr_mappings.unit

    # -----------------------------------------------------------------

    @property
    def sfr_mappings_unit(self):
        return self.cell_sfr_mappings.unit

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_density_mappings_unit(self):
        return self.sfr_mappings_unit / self.volume_unit

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_values(self):
        return self.cell_ssfr_mappings.values

    # -----------------------------------------------------------------

    @property
    def cell_sfr_mappings_values(self):
        return self.cell_sfr_mappings.values

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_densities_mappings(self):
        return self.cell_sfr_mappings_values / self.cell_volumes

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_mappings_mask(self):
        return np.isfinite(self.cell_ssfr_mappings_values) * (self.cell_ssfr_mappings_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mask_mappings(self):
        return self.valid_cell_ssfr_mappings_mask * self.valid_cell_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_values_mappings(self):
        return self.cell_ssfr_mappings_values[self.valid_cell_mask_mappings]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_funev_values_mappings(self):
        return self.cell_funev_values[self.valid_cell_mask_mappings]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_sfr_values_mappings(self):
        return self.cell_sfr_mappings_values[self.valid_cell_mask_mappings]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_sfr_densities_mappings(self):
        return self.cell_sfr_densities_mappings[self.valid_cell_mask_mappings]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_dust_mass_values_mappings(self):
        return self.cell_dust_mass_values[self.valid_cell_mask_mappings]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_dust_densities_mappings(self):
        return self.cell_dust_densities[self.valid_cell_mask_mappings]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_radii_mappings(self):
        return self.cell_radii[self.valid_cell_mask_mappings]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_bd_ratios_mappings(self):
        return self.cell_bd_ratios[self.valid_cell_mask_mappings]

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells_mappings.dat")

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_cells_colnames(self):
        return fs.get_column_names(self.ssfr_mappings_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.ssfr_mappings_funev_cells_colnames, ssfr_funev_base_colnames)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.ssfr_mappings_funev_cells_aux_colnames, aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_funev_cells(self):
        if fs.is_file(self.ssfr_mappings_funev_cells_path):
            if not self.ssfr_mappings_funev_cells_has_all_aux_columns:
                colnames = self.ssfr_mappings_funev_cells_aux_colnames
                #if sfr_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(sfr_name, self.valid_cell_sfr_values_mappings, self.sfr_mappings_unit, as_column=True)
                #if dust_mass_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(dust_mass_name, self.valid_cell_dust_mass_values_mappings, self.cell_dust_mass_unit, as_column=True)
                if sfr_density_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(sfr_density_name, self.valid_cell_sfr_densities_mappings, self.sfr_density_mappings_unit, as_column=True)
                if dust_density_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(dust_density_name, self.valid_cell_dust_densities_mappings, self.cell_dust_density_unit, as_column=True)
                if distance_center_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(distance_center_name, self.valid_cell_radii_mappings, self.length_unit, as_column=True)
                if bulge_disk_ratio_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(bulge_disk_ratio_name, self.valid_cell_bd_ratios_mappings, as_column=True)
                self.ssfr_mappings_funev_cells.save() # save
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_mappings_cells_aux(self):
        return {#sfr_name: self.valid_cell_sfr_values_mappings,
                #dust_mass_name: self.valid_cell_dust_mass_values_mappings,
                sfr_density_name: self.valid_cell_sfr_densities_mappings,
                dust_density_name: self.valid_cell_dust_densities_mappings,
                distance_center_name: self.valid_cell_radii_mappings,
                bulge_disk_ratio_name: self.valid_cell_bd_ratios_mappings}

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_mappings_cells_aux_units(self):
        return {#sfr_name: self.sfr_mappings_unit,
                #dust_mass_name: self.cell_dust_mass_unit,
                sfr_density_name: self.sfr_density_mappings_unit,
                dust_density_name: self.cell_dust_density_unit,
                distance_center_name: self.length_unit}

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_mappings_funev_cells_path", True, write=False)
    def ssfr_mappings_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_ssfr_mappings: raise IOError("The cell sSFR (MAPPINGS) data is not present: run the SFR analysis first")
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Create and return
        return Scatter2D.from_xy(self.valid_cell_ssfr_values_mappings, self.valid_cell_funev_values_mappings,
                                 x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_mappings_unit,
                                 x_description=self.ssfr_description, y_description=self.funev_description,
                                 aux=self.ssfr_mappings_cells_aux, aux_units=self.ssfr_mappings_cells_aux_units)

    # -----------------------------------------------------------------
    #   MAPPINGS + K&E
    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_ke_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_mappings_ke.dat")

    # -----------------------------------------------------------------

    @property
    def cell_sfr_mappings_ke_path(self):
        return fs.join(self.cell_sfr_path, "sfr_mappings_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_mappings_ke(self):
        return fs.is_file(self.cell_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_mappings_ke(self):
        return fs.is_file(self.cell_sfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_mappings_ke(self):
        return Data3D.from_file(self.cell_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_mappings_ke(self):
        return Data3D.from_file(self.cell_sfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_unit(self):
        return self.cell_ssfr_mappings_ke.unit

    # -----------------------------------------------------------------

    @property
    def sfr_mappings_ke_unit(self):
        return self.cell_sfr_mappings_ke.unit

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_density_mappings_ke_unit(self):
        return self.sfr_mappings_ke_unit / self.volume_unit

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_ke_values(self):
        return self.cell_ssfr_mappings_ke.values

    # -----------------------------------------------------------------

    @property
    def cell_sfr_mappings_ke_values(self):
        return self.cell_sfr_mappings_ke.values

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_densities_mappings_ke(self):
        return self.cell_sfr_mappings_ke_values / self.cell_volumes

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_mappings_ke_mask(self):
        return np.isfinite(self.cell_ssfr_mappings_ke_values) * (self.cell_ssfr_mappings_ke_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mask_mappings_ke(self):
        return self.valid_cell_ssfr_mappings_ke_mask * self.valid_cell_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_values_mappings_ke(self):
        return self.cell_ssfr_mappings_ke_values[self.valid_cell_mask_mappings_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_funev_values_mappings_ke(self):
        return self.cell_funev_values[self.valid_cell_mask_mappings_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_sfr_values_mappings_ke(self):
        return self.cell_sfr_mappings_ke_values[self.valid_cell_mask_mappings_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_sfr_densities_mappings_ke(self):
        return self.cell_sfr_densities_mappings_ke[self.valid_cell_mask_mappings_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_dust_mass_values_mappings_ke(self):
        return self.cell_dust_mass_values[self.valid_cell_mask_mappings_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_dust_densities_mappings_ke(self):
        return self.cell_dust_densities[self.valid_cell_mask_mappings_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_radii_mappings_ke(self):
        return self.cell_radii[self.valid_cell_mask_mappings_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_bd_ratios_mappings_ke(self):
        return self.cell_bd_ratios[self.valid_cell_mask_mappings_ke]

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells_mappings_ke.dat")

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_colnames(self):
        return fs.get_column_names(self.ssfr_mappings_ke_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.ssfr_mappings_ke_funev_cells_colnames, ssfr_funev_base_colnames)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.ssfr_mappings_ke_funev_cells_aux_colnames, aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_ke_funev_cells(self):
        if fs.is_file(self.ssfr_mappings_ke_funev_cells_path):
            if not self.ssfr_mappings_ke_funev_cells_has_all_aux_columns:
                colnames = self.ssfr_mappings_ke_funev_cells_aux_colnames
                #if sfr_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(sfr_name, self.valid_cell_sfr_values_mappings_ke, self.sfr_mappings_ke_unit, as_column=True)
                #if dust_mass_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(dust_mass_name, self.valid_cell_dust_mass_values_mappings_ke, self.cell_dust_mass_unit, as_column=True)
                if sfr_density_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(sfr_density_name, self.valid_cell_sfr_densities_mappings_ke, self.sfr_density_mappings_ke_unit, as_column=True)
                if dust_density_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(dust_density_name, self.valid_cell_dust_densities_mappings_ke, self.cell_dust_density_unit, as_column=True)
                if distance_center_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(distance_center_name, self.valid_cell_radii_mappings_ke, self.length_unit, as_column=True)
                if bulge_disk_ratio_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(bulge_disk_ratio_name, self.valid_cell_bd_ratios_mappings_ke, as_column=True)
                self.ssfr_mappings_ke_funev_cells.save() # save
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_mappings_ke_cells_aux(self):
        return {#sfr_name: self.valid_cell_sfr_values_mappings_ke,
                #dust_mass_name: self.valid_cell_dust_mass_values_mappings_ke,
                sfr_density_name: self.valid_cell_sfr_densities_mappings_ke,
                dust_density_name: self.valid_cell_dust_densities_mappings_ke,
                distance_center_name: self.valid_cell_radii_mappings_ke,
                bulge_disk_ratio_name: self.valid_cell_bd_ratios_mappings_ke}

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_mappings_ke_cells_aux_units(self):
        return {#sfr_name: self.sfr_mappings_ke_unit,
                #dust_mass_name: self.cell_dust_mass_unit,
                sfr_density_name: self.sfr_density_mappings_ke_unit,
                dust_density_name: self.cell_dust_density_unit,
                distance_center_name: self.length_unit}

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_mappings_ke_funev_cells_path", True, write=False)
    def ssfr_mappings_ke_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_ssfr_mappings_ke: raise IOError("The cell sSFR (MAPPINGS + K&E) data is not present: run the SFR analysis first")
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Create and return
        return Scatter2D.from_xy(self.valid_cell_ssfr_values_mappings_ke, self.valid_cell_funev_values_mappings_ke,
                                 x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_mappings_ke_unit,
                                 x_description=self.ssfr_description, y_description=self.funev_description,
                                 aux=self.ssfr_mappings_ke_cells_aux, aux_units=self.ssfr_mappings_ke_cells_aux_units)

    # -----------------------------------------------------------------
    # sSFR-Funev pixel scatter data
    # -----------------------------------------------------------------

    @property
    def projected_heating_maps_path(self):
        return fs.join(self.projected_heating_path, "maps")

    # -----------------------------------------------------------------

    @property
    def pixel_funev_path(self):
        return fs.join(self.projected_heating_maps_path, "earth_diffuse.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_funev(self):
        return fs.is_file(self.pixel_funev_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_funev(self):
        return Frame.from_file(self.pixel_funev_path)

    # -----------------------------------------------------------------

    @property
    def pixel_funev_values(self):
        return self.pixel_funev.values

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_salim_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_salim_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_salim(self):
        return fs.is_file(self.pixel_ssfr_salim_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_salim(self):
        return Frame.from_file(self.pixel_ssfr_salim_path)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_salim_values(self):
        return self.pixel_ssfr_salim.values

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_salim_mask(self):
        return np.isfinite(self.pixel_ssfr_salim_values) * (self.pixel_ssfr_salim_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_funev_mask(self):
        return np.isfinite(self.pixel_funev_values) * (self.pixel_funev_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_mask_salim(self):
        return self.valid_pixel_ssfr_salim_mask * self.valid_pixel_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_values_salim(self):
        return self.pixel_ssfr_salim_values[self.valid_pixel_mask_salim]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_funev_values_salim(self):
        return self.pixel_funev_values[self.valid_pixel_mask_salim]

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_salim.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_salim_funev_pixels(self):
        return fs.is_file(self.ssfr_salim_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_salim_funev_pixels_path", True, write=False)
    def ssfr_salim_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_pixel_ssfr_salim: raise IOError("The sSFR (Salim) frame is not present: run the SFR analysis first")
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Create and return
        #return Scatter2D.from_xy(self.pixel_ssfr_values, self.pixel_funev_values, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_unit, x_description=self.ssfr_description, y_description=self.funev_description)
        return Scatter2D.from_xy(self.valid_pixel_ssfr_values_salim, self.valid_pixel_funev_values_salim, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_salim_unit, x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    #   K&E
    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_ke_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_ke_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_ke(self):
        return fs.is_file(self.pixel_ssfr_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_ke(self):
        return Frame.from_file(self.pixel_ssfr_ke_path)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_ke_values(self):
        return self.pixel_ssfr_ke.values

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_ke_mask(self):
        return np.isfinite(self.pixel_ssfr_ke_values) * (self.pixel_ssfr_ke_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_mask_ke(self):
        return self.valid_pixel_ssfr_ke_mask * self.valid_pixel_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_values_ke(self):
        return self.pixel_ssfr_ke_values[self.valid_pixel_mask_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_funev_values_ke(self):
        return self.pixel_funev_values[self.valid_pixel_mask_ke]

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_ke_funev_pixels(self):
        return fs.is_file(self.ssfr_ke_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_ke_funev_pixels_path", True, write=False)
    def ssfr_ke_funev_pixels(self):

        """
        This functino ...
        :return:
        """

        # Checks
        if not self.has_pixel_ssfr_ke: raise IOError("The sSFR (K&E) frame is not present: run the SFR analysis first")
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Create and return
        return Scatter2D.from_xy(self.valid_pixel_ssfr_values_ke, self.valid_pixel_funev_values_ke,
                                 x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_ke_unit,
                                 x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    #   MAPPINGS
    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_mappings_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_mappings_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_mappings(self):
        return fs.is_file(self.pixel_ssfr_mappings_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_mappings(self):
        return Frame.from_file(self.pixel_ssfr_mappings_path)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_mappings_values(self):
        return self.pixel_ssfr_mappings.values

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_mappings_mask(self):
        return np.isfinite(self.pixel_ssfr_mappings_values) * (self.pixel_ssfr_mappings_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_mask_mappings(self):
        return self.valid_pixel_ssfr_mappings_mask * self.valid_pixel_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_values_mappings(self):
        return self.pixel_ssfr_mappings_values[self.valid_pixel_mask_mappings]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_funev_values_mappings(self):
        return self.pixel_funev_values[self.valid_pixel_mask_mappings]

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_mappings.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_funev_pixels(self):
        return fs.is_file(self.ssfr_mappings_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_mappings_funev_pixels_path", True, write=False)
    def ssfr_mappings_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_pixel_ssfr_mappings: raise IOError("The sSFR (MAPPINGS) frame is not present: run the SFR analysis first")
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Create and return
        return Scatter2D.from_xy(self.valid_pixel_ssfr_values_mappings, self.valid_pixel_funev_values_mappings,
                                 x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_mappings_unit,
                                 x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    #   MAPPINGS + K&E
    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_mappings_ke_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_mappings_ke_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_mappings_ke(self):
        return fs.is_file(self.pixel_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_mappings_ke(self):
        return Frame.from_file(self.pixel_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_mappings_ke_values(self):
        return self.pixel_ssfr_mappings_ke.values

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_mappings_ke_mask(self):
        return np.isfinite(self.pixel_ssfr_mappings_ke_values) * (self.pixel_ssfr_mappings_ke_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_mask_mappings_ke(self):
        return self.valid_pixel_ssfr_mappings_ke_mask * self.valid_pixel_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_values_mappings_ke(self):
        return self.pixel_ssfr_mappings_ke_values[self.valid_pixel_mask_mappings_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_funev_values_mappings_ke(self):
        return self.pixel_funev_values[self.valid_pixel_mask_mappings_ke]

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_mappings_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_ke_funev_pixels(self):
        return fs.is_file(self.ssfr_mappings_ke_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_mappings_ke_funev_pixels_path", True, write=False)
    def ssfr_mappings_ke_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_pixel_ssfr_mappings_ke: raise IOError("The sSFR (MAPPINGS + K&E) frame is not present: run the SFR analysis first")
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Create and return
        return Scatter2D.from_xy(self.valid_pixel_ssfr_values_mappings_ke, self.valid_pixel_funev_values_mappings_ke, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_mappings_ke_unit, x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    # Temperature-Funev
    # -----------------------------------------------------------------

    @property
    def has_temperature_cells(self):
        return True

    # -----------------------------------------------------------------

    @property
    def has_temperature_pixels(self):
        return False # currently, I know no direct way of getting a temperature map for an instrument

    # -----------------------------------------------------------------

    @lazyproperty
    def temperature_funev_path(self):
        return fs.create_directory_in(self.correlations_path, temperature_funev_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_temperature_mask(self):
        return self.cell_temperatures != 0

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mask_temperature_funev(self):
        return self.valid_cell_temperature_mask * self.valid_cell_funev_mask

    # -----------------------------------------------------------------

    @property
    def temperature_funev_cells_path(self):
        return fs.join(self.temperature_funev_path, "cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_temperature_funev_cells(self):
        return fs.is_file(self.temperature_funev_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "temperature_funev_cells_path", True, write=False)
    def temperature_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_temperature_cells: raise IOError("The cell temperature data is not present")
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Get values
        valid_cell_temperatures = self.cell_temperatures[self.valid_cell_mask_temperature_funev]
        valid_cell_funev_values = self.cell_funev_values[self.valid_cell_mask_temperature_funev]

        # Create and return
        return Scatter2D.from_xy(valid_cell_temperatures, valid_cell_funev_values,
                                 x_name=self.temperature_name, y_name=self.funev_name, x_unit=self.temperature_unit,
                                 x_description=self.temperature_description, y_description=self.funev_description)

    # -----------------------------------------------------------------

    @property
    def temperature_funev_pixels_path(self):
        return fs.join(self.temperature_funev_path, "pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_temperature_funev_pixels(self):
        return fs.is_file(self.temperature_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "temperature_funev_pixels_path", True, write=False)
    def temperature_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_temperature_pixels: raise IOError("The temperature frame cannot be obtained")
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        valid_pixel_temperatures = valid_pixel_funev_values = None

        # Create and return
        return Scatter2D.from_xy(valid_pixel_temperatures, valid_pixel_funev_values,
                                 x_name=self.temperature_name, y_name=self.funev_name, x_unit=self.temperature_unit,
                                 x_description=self.temperature_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def do_ssfr_salim_funev_cells(self):
        return self.has_cell_ssfr_salim and self.has_cell_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_salim_funev_pixels(self):
        return self.has_pixel_ssfr_salim and self.has_pixel_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_ke_funev_cells(self):
        return self.has_cell_ssfr_ke and self.has_cell_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_ke_funev_pixels(self):
        return self.has_pixel_ssfr_ke and self.has_pixel_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_mappings_funev_cells(self):
        return self.has_cell_ssfr_mappings and self.has_cell_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_mappings_funev_pixels(self):
        return self.has_pixel_ssfr_mappings and self.has_pixel_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_mappings_ke_funev_cells(self):
        return self.has_cell_ssfr_mappings_ke and self.has_cell_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_mappings_ke_funev_pixels(self):
        return self.has_pixel_ssfr_mappings_ke and self.has_pixel_funev

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # sSFR Funev scatter data
        self.write_ssfr_funev()

        # Temperature Funev scatter data
        self.write_temperature_funev()

    # -----------------------------------------------------------------

    def write_ssfr_funev(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR to Funev scatter data ...")

        # Salim
        self.write_ssfr_salim_funev()

        # K&E
        self.write_ssfr_ke_funev()

        # MAPPINGS
        self.write_ssfr_mappings_funev()

        # MAPPINGS + K&E
        self.write_ssfr_mappings_ke_funev()

        # Write the M51 scatter data
        if self.do_write_m51_ssfr_funev: self.write_m51_ssfr_funev()

        # Write the M31 scatter data
        if self.do_write_m31_ssfr_funev: self.write_m31_ssfr_funev()

    # -----------------------------------------------------------------
    # SALIM
    # -----------------------------------------------------------------

    def write_ssfr_salim_funev(self):

        """
        This function ...
        :return:
        """

        # Write the sSFR to Funev cell scatter data
        if self.do_write_ssfr_salim_funev_cells: self.write_ssfr_salim_funev_cells()

        # Write the sSFR to Funev pixel scatter data
        if self.do_write_ssfr_salim_funev_pixels: self.write_ssfr_salim_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_salim_funev_cells(self):
        return self.do_ssfr_salim_funev_cells and not self.has_ssfr_salim_funev_cells

    # -----------------------------------------------------------------

    def write_ssfr_salim_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR to Funev dust cell scatter data ...")

        # Write
        self.ssfr_salim_funev_cells.saveto(self.ssfr_salim_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_salim_funev_pixels(self):
        return self.do_ssfr_salim_funev_pixels and not self.has_ssfr_salim_funev_pixels

    # -----------------------------------------------------------------

    def write_ssfr_salim_funev_pixels(self):

        """
        This fnction ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (Salim) to Funev pixel scatter data ...")

        # Write
        self.ssfr_salim_funev_pixels.saveto(self.ssfr_salim_funev_pixels_path)

    # -----------------------------------------------------------------
    # KENNICUTT & EVANS
    # -----------------------------------------------------------------

    def write_ssfr_ke_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        self.write_ssfr_ke_funev_cells()

        # Pixels
        self.write_ssfr_ke_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_ke_funev_cells(self):
        return self.do_ssfr_ke_funev_cells and not self.has_ssfr_ke_funev_cells

    # -----------------------------------------------------------------

    def write_ssfr_ke_funev_cells(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (K&E) to Funev cell scatter data ...")

        # Write
        self.ssfr_ke_funev_cells.saveto(self.ssfr_ke_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_ke_funev_pixels(self):
        return self.do_ssfr_ke_funev_pixels and not self.has_ssfr_ke_funev_pixels

    # -----------------------------------------------------------------

    def write_ssfr_ke_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (K&E) to Funev pixel scatter data ...")

        # Write
        self.ssfr_ke_funev_pixels.saveto(self.ssfr_ke_funev_pixels_path)

    # -----------------------------------------------------------------
    # MAPPINGS
    # -----------------------------------------------------------------

    def write_ssfr_mappings_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_ssfr_mappings_funev_cells: self.write_ssfr_mappings_funev_cells()

        # Pixels
        if self.do_write_ssfr_mappings_funev_pixels: self.write_ssfr_mappings_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_mappings_funev_cells(self):
        return self.do_ssfr_mappings_funev_cells and not self.has_ssfr_mappings_funev_cells

    # -----------------------------------------------------------------

    def write_ssfr_mappings_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (MAPPINGS) to Funev cell scatter data ...")

        # Write
        self.ssfr_mappings_funev_cells.saveto(self.ssfr_mappings_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_mappings_funev_pixels(self):
        return self.do_ssfr_mappings_funev_pixels and not self.has_ssfr_mappings_funev_pixels

    # -----------------------------------------------------------------

    def write_ssfr_mappings_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (MAPPINGS) to Funev pixel scatter data ...")

        # Write
        self.ssfr_mappings_funev_pixels.saveto(self.ssfr_mappings_funev_pixels_path)

    # -----------------------------------------------------------------
    # MAPPINGS + KENNICUTT & EVANS
    # -----------------------------------------------------------------

    def write_ssfr_mappings_ke_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_ssfr_mappings_ke_funev_cells: self.write_ssfr_mappings_ke_funev_cells()

        # Pixels
        if self.do_write_ssfr_mappings_ke_funev_pixels: self.write_ssfr_mappings_ke_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_mappings_ke_funev_cells(self):
        return self.do_ssfr_mappings_ke_funev_cells and not self.has_ssfr_mappings_ke_funev_cells

    # -----------------------------------------------------------------

    def write_ssfr_mappings_ke_funev_cells(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (MAPPINGS + K&E) to Funev cell scatter data ...")

        # Write
        self.ssfr_mappings_ke_funev_cells.saveto(self.ssfr_mappings_ke_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_mappings_ke_funev_pixels(self):
        return self.do_ssfr_mappings_ke_funev_pixels and not self.has_ssfr_mappings_ke_funev_pixels

    # -----------------------------------------------------------------

    def write_ssfr_mappings_ke_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (MAPPINGS + K&E) to Funev pixel scatter data ...")

        # Write
        self.ssfr_mappings_ke_funev_pixels.saveto(self.ssfr_mappings_ke_funev_pixels_path)

    # -----------------------------------------------------------------
    # M51
    # -----------------------------------------------------------------

    @property
    def do_write_m51_ssfr_funev(self):
        return not self.has_m51_ssfr_funev

    # -----------------------------------------------------------------

    def write_m51_ssfr_funev(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the M51 sSFR to Funev scatter data ...")

        # Write
        self.m51_scatter.saveto(self.m51_ssfr_funev_path)

    # -----------------------------------------------------------------
    # M31
    # -----------------------------------------------------------------

    @property
    def do_write_m31_ssfr_funev(self):
        return not self.has_m31_ssfr_funev

    # -----------------------------------------------------------------

    def write_m31_ssfr_funev(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the M31 sSFR to Funev scatter data ...")

        # Write
        self.m31_scatter.saveto(self.m31_ssfr_funev_path)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def do_temperature_cells(self):
        return self.has_temperature_cells

    # -----------------------------------------------------------------

    @property
    def do_temperature_pixels(self):
        return self.has_temperature_pixels

    # -----------------------------------------------------------------

    @property
    def do_write_temperature_funev_cells(self):
        return self.do_temperature_cells and not self.has_temperature_funev_cells

    # -----------------------------------------------------------------

    @property
    def do_write_temperature_funev_pixels(self):
        return self.do_temperature_pixels and not self.has_temperature_funev_pixels

    # -----------------------------------------------------------------

    def write_temperature_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_temperature_funev_cells: self.write_temperature_funev_cells()

        # Pixels
        if self.do_write_temperature_funev_pixels: self.write_temperature_funev_pixels()

    # -----------------------------------------------------------------

    def write_temperature_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the temperature to Funev dust cell scatter data ...")

        # Write
        self.temperature_funev_cells.saveto(self.temperature_funev_cells_path)

    # -----------------------------------------------------------------

    def write_temperature_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the temperature to Funev pixel scatter data ...")

        # Write
        self.temperature_funev_pixels.saveto(self.temperature_funev_pixels_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot
        self.plot_ssfr_funev()

        # Plot
        self.plot_temperature_funev()

    # -----------------------------------------------------------------

    def plot_ssfr_funev(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR to Funev scatter data ...")

        # Salim
        self.plot_ssfr_salim_funev()

        # K&E
        self.plot_ssfr_ke_funev()

        # MAPPINGS
        self.plot_ssfr_mappings_funev()

        # MAPPINGS + K&E
        self.plot_ssfr_mappings_ke_funev()

    # -----------------------------------------------------------------

    def plot_ssfr_salim_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_plot_ssfr_salim_funev_cells: self.plot_ssfr_salim_funev_cells()

        # Pixels
        if self.do_plot_ssfr_salim_funev_pixels: self.plot_ssfr_salim_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_plot_path(self):
        return fs.join(self.ssfr_funev_path, "cells_salim.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_salim_funev_cells_plot(self):
        return fs.is_file(self.ssfr_salim_funev_cells_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_salim_funev_cells(self):
        return self.do_ssfr_salim_funev_cells and not self.has_ssfr_salim_funev_cells_plot

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_salim_funev_cells
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_scatter_paths(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_salim_funev_cells_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_title(self):
        return "Correlation between sSFR (Salim) and Funev of dust cells"

    # -----------------------------------------------------------------

    @property
    def ssfr_limits(self):
        return (1e-13,1e-9,)

    # -----------------------------------------------------------------

    @property
    def funev_limits(self):
        return (0,1,)

    # -----------------------------------------------------------------

    def plot_ssfr_salim_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR (Salim) to Funev dust cell scatter data ...")

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_salim_funev_cells_scatter_paths, self.ssfr_name, self.funev_name,
                                           self.ssfr_description, self.funev_description,
                                           title=self.ssfr_salim_funev_cells_title, path=self.ssfr_salim_funev_cells_plot_path,
                                           ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_salim_funev_cells_scatters, title=self.ssfr_salim_funev_cells_title, xlog=True,
                      path=self.ssfr_salim_funev_cells_plot_path, xlimits=self.ssfr_limits, ylimits=self.funev_limits,
                      density=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_pixels_plot_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_salim.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_salim_funev_pixels_plot(self):
        return fs.is_file(self.ssfr_salim_funev_pixels_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_salim_funev_pixels(self):
        return self.do_ssfr_salim_funev_pixels and not self.has_ssfr_salim_funev_pixels_plot

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_pixels_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_salim_funev_pixels
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_pixels_scatter_paths(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_salim_funev_pixels_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_pixels_title(self):
        return "Correlation between sSFR (Salim) and Funev of model pixels"

    # -----------------------------------------------------------------

    def plot_ssfr_salim_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR (Salim) to Funev pixel scatter data ...")

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_salim_funev_pixels_scatter_paths, self.ssfr_name, self.funev_name,
                                           self.ssfr_description, self.funev_description,
                                           title=self.ssfr_salim_funev_pixels_title, path=self.ssfr_salim_funev_pixels_plot_path,
                                           ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_salim_funev_pixels_scatters, title=self.ssfr_salim_funev_pixels_title, xlog=True,
                      path=self.ssfr_salim_funev_pixels_plot_path, xlimits=self.ssfr_limits, ylimits=self.funev_limits,
                      density=True)

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_ke_funev_cells(self):
        return self.do_ssfr_ke_funev_cells and not self.has_ssfr_ke_funev_cells_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_ke_funev_pixels(self):
        return self.do_ssfr_ke_funev_pixels and not self.has_ssfr_ke_funev_pixels_plot

    # -----------------------------------------------------------------

    def plot_ssfr_ke_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_plot_ssfr_ke_funev_cells: self.plot_ssfr_ke_funev_cells()

        # Pixels
        if self.do_plot_ssfr_ke_funev_pixels: self.plot_ssfr_ke_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_plot_path(self):
        return fs.join(self.ssfr_funev_path, "cells_ke.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_ke_funev_cells_plot(self):
        return fs.is_file(self.ssfr_ke_funev_cells_plot_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_ke_funev_cells
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_scatter_paths(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_ke_funev_cells_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_title(self):
        return "Correlation between sSFR (K&E) and Funev of dust cells"

    # -----------------------------------------------------------------

    def plot_ssfr_ke_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR (K&E) to Funev cell scatter data ...")

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_ke_funev_cells_scatter_paths, self.ssfr_name, self.funev_name,
                                            self.ssfr_description, self.funev_description,
                                            title=self.ssfr_ke_funev_cells_title, path=self.ssfr_ke_funev_cells_plot_path,
                                            ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_ke_funev_cells_scatters, title=self.ssfr_ke_funev_cells_title, xlog=True,
                          path=self.ssfr_ke_funev_cells_plot_path, xlimits=self.ssfr_limits,
                          ylimits=self.funev_limits,
                          density=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_pixels_plot_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_ke.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_ke_funev_pixels_plot(self):
        return fs.is_file(self.ssfr_ke_funev_pixels_plot_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_pixels_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_ke_funev_pixels
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_pixels_scatter_paths(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_ke_funev_pixels_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_pixels_title(self):
        return "Correlation between sSFR (K&E) and Funev of model pixels"

    # -----------------------------------------------------------------

    def plot_ssfr_ke_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_ke_funev_pixels_scatter_paths, self.ssfr_name, self.funev_name,
                                           self.ssfr_description, self.funev_description,
                                           title=self.ssfr_ke_funev_pixels_title, path=self.ssfr_ke_funev_pixels_plot_path,
                                           ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_ke_funev_pixels_scatters, title=self.ssfr_ke_funev_pixels_title, xlog=True,
                          path=self.ssfr_ke_funev_pixels_plot_path, xlimits=self.ssfr_limits,
                          ylimits=self.funev_limits,
                          density=True)

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_mappings_funev_cells(self):
        return self.do_ssfr_mappings_funev_cells and not self.has_ssfr_mappings_funev_cells_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_mappings_funev_pixels(self):
        return self.do_ssfr_mappings_funev_pixels and not self.has_ssfr_mappings_funev_pixels_plot

    # -----------------------------------------------------------------

    def plot_ssfr_mappings_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_plot_ssfr_mappings_funev_cells: self.plot_ssfr_mappings_funev_cells()

        # Pixels
        if self.do_plot_ssfr_mappings_funev_pixels: self.plot_ssfr_mappings_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_cells_plot_path(self):
        return fs.join(self.ssfr_funev_path, "cells_mappings.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_funev_cells_plot(self):
        return fs.is_file(self.ssfr_mappings_funev_cells_plot_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_cells_scatters(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_mappings_funev_cells
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_cells_scatter_paths(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_mappings_funev_cells_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_cells_title(self):
        return "Correlation between sSFR (MAPPINGS) and Funev of dust cells"

    # -----------------------------------------------------------------

    def plot_ssfr_mappings_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat:

            plot_stilts(self.ssfr_mappings_funev_cells_scatter_paths, self.ssfr_name, self.funev_name,
                        self.ssfr_description, self.funev_description,
                        title=self.ssfr_mappings_funev_cells_title, path=self.ssfr_mappings_funev_cells_plot_path,
                        ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else:
            plot_scatters(self.ssfr_mappings_funev_cells_scatters, title=self.ssfr_mappings_funev_cells_title,
                          xlog=True, path=self.ssfr_mappings_funev_cells_plot_path, xlimits=self.ssfr_limits,
                          ylimits=self.funev_limits, density=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_pixels_plot_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_mappings.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_funev_pixels_plot(self):
        return fs.is_file(self.ssfr_mappings_funev_pixels_plot_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_pixels_scatters(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_mappings_funev_pixels
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_pixels_scatter_paths(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_mappings_funev_pixels_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_funev_pixels_title(self):
        return "Correlation between sSFR (MAPPINGS) and Funev of model pixels"

    # -----------------------------------------------------------------

    def plot_ssfr_mappings_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat:

            plot_stilts(self.ssfr_mappings_funev_pixels_scatter_paths, self.ssfr_name, self.funev_name,
                        self.ssfr_description, self.funev_description,
                        title=self.ssfr_mappings_funev_pixels_title,
                        path=self.ssfr_mappings_funev_pixels_plot_path,
                        ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else:
            plot_scatters(self.ssfr_mappings_funev_pixels_scatters, title=self.ssfr_mappings_funev_pixels_title,
                          xlog=True, path=self.ssfr_mappings_funev_pixels_plot_path, xlimits=self.ssfr_limits,
                          ylimits=self.funev_limits, density=True)

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_mappings_ke_funev_cells(self):
        return self.do_ssfr_mappings_ke_funev_cells and not self.has_ssfr_mappings_ke_funev_cells_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_mappings_ke_funev_pixels(self):
        return self.do_ssfr_mappings_ke_funev_pixels and not self.has_ssfr_mappings_ke_funev_pixels_plot

    # -----------------------------------------------------------------

    def plot_ssfr_mappings_ke_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_plot_ssfr_mappings_ke_funev_cells: self.plot_ssfr_mappings_ke_funev_cells()

        # Pixels
        if self.do_plot_ssfr_mappings_ke_funev_pixels: self.plot_ssfr_mappings_ke_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_plot_path(self):
        return fs.join(self.ssfr_funev_path, "cells_mappings_ke.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_ke_funev_cells_plot(self):
        return fs.is_file(self.ssfr_mappings_ke_funev_cells_plot_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_scatters(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_mappings_ke_funev_cells
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_scatter_paths(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_mappings_ke_funev_cells_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_title(self):
        return "Correlation between sSFR (MAPPINGS + K&E) and Funev of dust cells"

    # -----------------------------------------------------------------

    def plot_ssfr_mappings_ke_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_mappings_ke_funev_cells_scatter_paths, self.ssfr_name, self.funev_name,
                                            self.ssfr_description, self.funev_description,
                                            title=self.ssfr_mappings_ke_funev_cells_title, path=self.ssfr_mappings_ke_funev_cells_plot_path,
                                            ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_mappings_ke_funev_cells_scatters, title=self.ssfr_mappings_ke_funev_cells_title, xlog=True,
                              path=self.ssfr_mappings_ke_funev_cells_plot_path, xlimits=self.ssfr_limits,
                              ylimits=self.funev_limits,
                              density=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_pixels_plot_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_mappings_ke.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_ke_funev_pixels_plot(self):
        return fs.is_file(self.ssfr_mappings_ke_funev_pixels_plot_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_pixels_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_mappings_ke_funev_pixels
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_pixels_scatter_paths(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_mappings_ke_funev_pixels_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_pixels_title(self):
        return "Correlation between sSFR (MAPPINGS + K&E) and Funev of model pixels"

    # -----------------------------------------------------------------

    def plot_ssfr_mappings_ke_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_mappings_ke_funev_pixels_scatter_paths, self.ssfr_name, self.funev_name,
                                        self.ssfr_description, self.funev_description,
                                        title=self.ssfr_mappings_ke_funev_pixels_title, path=self.ssfr_mappings_ke_funev_pixels_plot_path,
                                        ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_mappings_ke_funev_pixels_scatters, title=self.ssfr_mappings_ke_funev_pixels_title, xlog=True,
                              path=self.ssfr_mappings_ke_funev_pixels_plot_path, xlimits=self.ssfr_limits,
                              ylimits=self.funev_limits,
                              density=True)

    # -----------------------------------------------------------------

    @property
    def do_plot_temperature_funev_cells(self):
        return self.do_temperature_cells and not self.has_temperature_funev_cells_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_temperature_funev_pixels(self):
        return self.do_temperature_pixels and not self.has_temperature_funev_pixels_plot

    # -----------------------------------------------------------------

    @property
    def temperature_limits(self):
        return None

    # -----------------------------------------------------------------

    def plot_temperature_funev(self):

        """
        This function ...
        :return:
        """
        
        # Cells
        if self.do_plot_temperature_funev_cells: self.plot_temperature_funev_cells()
        
        # Pixels
        if self.do_plot_temperature_funev_pixels: self.plot_temperature_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def temperature_funev_cells_plot_path(self):
        return fs.join(self.temperature_funev_path, "cells.pdf")

    # -----------------------------------------------------------------

    @property
    def has_temperature_funev_cells_plot(self):
        return fs.is_file(self.temperature_funev_cells_plot_path)

    # -----------------------------------------------------------------

    @property
    def temperature_funev_cells_title(self):
        return "Correlation between dust temperature and Funev of dust cells"

    # -----------------------------------------------------------------

    @property
    def temperature_funev_cells_scatters(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.temperature_funev_cells
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    @property
    def temperature_funev_cells_scatter_paths(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.temperature_funev_cells_path
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    def plot_temperature_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat:

            plot_stilts(self.temperature_funev_cells_scatter_paths, self.temperature_name, self.funev_name,
                        self.temperature_description, self.funev_description,
                        title=self.temperature_funev_cells_title, path=self.temperature_funev_cells_plot_path,
                        ylimits=self.funev_limits, xlog=True, xlimits=self.temperature_limits)

        # Plot using Matplotlib
        else: raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def temperature_funev_pixels_plot_path(self):
        return fs.join(self.temperature_funev_path, "pixels.pdf")

    # -----------------------------------------------------------------

    @property
    def has_temperature_funev_pixels_plot(self):
        return fs.is_file(self.temperature_funev_pixels_plot_path)

    # -----------------------------------------------------------------

    @property
    def temperature_funev_pixels_title(self):
        return "Correlation between dust temperature and Funev of pixels"

    # -----------------------------------------------------------------

    @property
    def temperature_funev_pixels_scatters(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.temperature_funev_pixels
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    @property
    def temperature_funev_pixels_scatter_paths(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.temperature_funev_pixels_path
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    def plot_temperature_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat:

            plot_stilts(self.temperature_funev_pixels_scatter_paths, self.temperature_name, self.funev_name,
                        self.temperature_description, self.funev_description,
                        title=self.temperature_funev_pixels_title,
                        path=self.temperature_funev_pixels_plot_path,
                        ylimits=self.funev_limits, xlog=True, xlimits=self.temperature_limits)

        # Plot using Matplotlib
        else: raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------
