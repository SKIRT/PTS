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

# Import astronomical modules
from astropy.units import DexUnit

# Import the relevant PTS classes and modules
from .component import AnalysisRunComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.tools import introspection
from ...core.tools.utils import lazyproperty, lazyfileproperty
from ...core.basics.scatter import Scatter2D
from ..core.data import Data3D
from ...magic.core.frame import Frame
from ...magic.tools.plotting import plot_stilts, plot_scatters_density
from ...core.units.parsing import parse_unit as u
from ...core.tools import sequences
from ...magic.core.list import uniformize

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
sfr_sfr_name = "SFR-SFR"
mean_age_funev_name = "Mean age-Funev"
mean_age_ssfr_name = "Mean age-sSFR"

# -----------------------------------------------------------------

# Auxilary column names
sfr_density_name = "vSFR"
dust_density_name = "Dust density"
distance_center_name = "Distance from center"
radius_name = "Radius"
dust_heights_name = "Dust scale heights"
bulge_disk_ratio_name = "Bulge disk ratio"
temperature_name = "Dust temperature"
luminosity_density_name = "Dust luminosity density"
mean_age_name = "Mean stellar age"
funev_name = "Fraction of heating by unevolved stars"
sfr_density_salim_name = "vSFR (Salim)"
sfr_density_ke_name = "vSFR (K&E)"
sfr_density_mappings_name = "vSFR (MAPPINGS)"
sfr_density_mappings_ke_name = "vSFR (MAPPINGS + K&E)"

# -----------------------------------------------------------------

# Auxilary column names for sSFR-Funev scatter data
ssfr_funev_aux_colnames = [sfr_density_name, dust_density_name, distance_center_name, radius_name, dust_heights_name, bulge_disk_ratio_name, temperature_name, luminosity_density_name, mean_age_name]

# Auxilary column names for SFR-SFR scatter
sfr_sfr_cells_aux_colnames = [temperature_name, mean_age_name, funev_name]

# Auxilary column names for colour (sSFR) - Funev scatter data
colour_funev_aux_colnames = [sfr_density_salim_name, sfr_density_ke_name, sfr_density_mappings_name, sfr_density_mappings_ke_name, dust_density_name, distance_center_name, bulge_disk_ratio_name, temperature_name, mean_age_name]

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

    @property
    def energy_path(self):
        return self.analysis_run.energy_path

    # -----------------------------------------------------------------

    @property
    def projected_energy_path(self):
        return fs.join(self.energy_path, "projected")

    # -----------------------------------------------------------------

    @property
    def cell_energy_path(self):
        return fs.join(self.energy_path, "cell")

    # -----------------------------------------------------------------

    @property
    def cell_dust_luminosity_density_data_path(self):
        return fs.join(self.cell_energy_path, "total_density.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_dust_luminosity_density_data(self):
        return Data3D.from_file(self.cell_dust_luminosity_density_data_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_dust_luminosity_densities(self):
        return self.cell_dust_luminosity_density_data.values

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_dust_luminosity_density_unit(self):
        return self.cell_dust_luminosity_density_data.unit

    # -----------------------------------------------------------------
    # CELL PROPERTIES
    # -----------------------------------------------------------------

    @lazyproperty
    def cell_radii(self):
        return np.sqrt(self.cell_x_coordinates**2 + self.cell_y_coordinates**2 + self.cell_z_coordinates**2)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_polar_radii(self):
        return np.sqrt(self.cell_x_coordinates**2 + self.cell_y_coordinates**2)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_dust_heights(self):
        return np.abs(self.cell_z_coordinates) / self.model.dust_scaleheight.to(self.length_unit).value

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
    # BOLOMETRIC LUMINOSITY / STELLAR CONTRIBUTION
    # -----------------------------------------------------------------

    @lazyproperty
    def bol_luminosity_unit(self):
        return u("Lsun")

    # -----------------------------------------------------------------

    @property
    def bulge_bolometric_luminosity(self):
        return self.model.intrinsic_bolometric_luminosity_old_bulge
        # SHOULD BE THE SAME (not tested):
        #return self.model.observed_bolometric_luminosity_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_bolometric_luminosity_scalar(self):
        return self.bulge_bolometric_luminosity.to(self.bol_luminosity_unit, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @property
    def disk_bolometric_luminosity(self):
        return self.model.intrinsic_bolometric_luminosity_old_disk
        # SHOULD BE THE SAME (not tested):
        #return self.model.observed_bolometric_luminosity_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_bolometric_luminosity_scalar(self):
        return self.disk_bolometric_luminosity.to(self.bol_luminosity_unit, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @property
    def young_bolometric_luminosity(self):
        return self.model.intrinsic_bolometric_luminosity_young
        # SHOULD BE THE SAME (not tested):
        #return self.model.observed_bolometric_luminosity_young

    # -----------------------------------------------------------------

    @lazyproperty
    def young_bolometric_luminosity_scalar(self):
        return self.young_bolometric_luminosity.to(self.bol_luminosity_unit, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @property
    def sfr_bolometric_luminosity(self):
        return self.model.intrinsic_bolometric_luminosity_sfr
        # SHOULD BE THE SAME (not tested):
        #return self.model.observed_bolometric_luminosity_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_bolometric_luminosity_scalar(self):
        return self.sfr_bolometric_luminosity.to(self.bol_luminosity_unit, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_cell_bol_luminosities(self):
        return self.bulge_cell_normalized_mass * self.bulge_bolometric_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_cell_bol_luminosities(self):
        return self.disk_cell_normalized_mass * self.disk_bolometric_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def old_cell_bol_luminosities(self):
        return self.bulge_cell_bol_luminosities + self.disk_cell_bol_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def young_cell_bol_luminosities(self):
        return self.young_cell_normalized_mass * self.young_bolometric_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_cell_bol_luminosities(self):
        return self.sfr_cell_normalized_mass * self.sfr_bolometric_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_cell_bol_luminosities(self):
        return self.young_cell_bol_luminosities + self.sfr_cell_bol_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def total_cell_bol_luminosities(self):
        return self.old_cell_bol_luminosities + self.unevolved_cell_bol_luminosities

    # -----------------------------------------------------------------
    # MEAN STELLAR AGE
    # -----------------------------------------------------------------

    @property
    def age_unit(self):
        return u("Myr")

    # -----------------------------------------------------------------

    @lazyproperty
    def log_age_unit(self):
        return DexUnit("Myr")

    # -----------------------------------------------------------------

    @property
    def cell_mean_age_path(self):
        return fs.join(self.correlations_path, "mean_age.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_mean_age(self):
        return fs.is_file(self.cell_mean_age_path)

    # -----------------------------------------------------------------

    @property
    def mean_age_name(self):
        return "Age"

    # -----------------------------------------------------------------

    @property
    def mean_age_description(self):
        return "Mean stellar age"

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stellar_age_scalar(self):
        return self.model.old_mean_stellar_age.to(self.age_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def log_old_stellar_age(self):
        return np.log10(self.old_stellar_age_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stellar_age_scalar(self):
        return self.model.young_mean_stellar_age.to(self.age_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def log_young_stellar_age(self):
        return np.log10(self.young_stellar_age_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_stellar_age_scalar(self):
        return self.model.sfr_mean_stellar_age.to(self.age_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def log_sfr_stellar_age(self):
        return np.log10(self.sfr_stellar_age_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_ages(self):
        return np.array([self.old_stellar_age_scalar, self.young_stellar_age_scalar, self.sfr_stellar_age_scalar])

    # -----------------------------------------------------------------

    @lazyproperty
    def log_stellar_ages(self):
        return np.array([self.log_old_stellar_age, self.log_young_stellar_age, self.log_sfr_stellar_age])

    # -----------------------------------------------------------------

    @property
    def cell_old_bol_fraction_path(self):
        return fs.join(self.correlations_path, "old_bol_fraction.dat")

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_old_bol_fraction_path", True, write=True)
    def cell_old_bol_fraction_data(self):

        """
        This function ...
        :return:
        """

        # Calculate fractions
        fractions = self.old_cell_bol_luminosities / self.total_cell_bol_luminosities

        # Create the data with external xyz
        return Data3D.from_values("fbol_old", fractions, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_old_luminosity_contributions(self):
        return self.cell_old_bol_fraction_data.values

    # -----------------------------------------------------------------

    @property
    def cell_young_bol_fraction_path(self):
        return fs.join(self.correlations_path, "young_bol_fraction.dat")

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_young_bol_fraction_path", True, write=True)
    def cell_young_bol_fraction_data(self):

        """
        This function ...
        :return:
        """

        # Calculate fractions
        fractions = self.young_cell_bol_luminosities / self.total_cell_bol_luminosities

        # Create the data with external xyz
        return Data3D.from_values("fbol_young", fractions, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_young_luminosity_contributions(self):
        return self.cell_young_bol_fraction_data.values

    # -----------------------------------------------------------------

    @property
    def cell_sfr_bol_fraction_path(self):
        return fs.join(self.correlations_path, "sfr_bol_fraction.dat")

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_sfr_bol_fraction_path", True, write=True)
    def cell_sfr_bol_fraction_data(self):

        """
        This function ...
        :return:
        """

        # Calculate fractions
        fractions = self.sfr_cell_bol_luminosities / self.total_cell_bol_luminosities

        # Create the data with external xyz
        return Data3D.from_values("fbol_sfr", fractions, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_luminosity_contributions(self):
        return self.cell_sfr_bol_fraction_data.values

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_stellar_luminosity_contributions(self):
        return np.stack([self.cell_old_luminosity_contributions, self.cell_young_luminosity_contributions, self.cell_sfr_luminosity_contributions])

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_mean_age_path", True, write=True)
    def cell_mean_age_data(self):

        """
        This function ...
        :return:
        """

        # Check normalization of the weights
        total = self.cell_old_luminosity_contributions + self.cell_young_luminosity_contributions + self.cell_sfr_luminosity_contributions
        valid_total = total[np.isfinite(total)]
        close = np.isclose(valid_total, 1.)

        if not np.all(close): raise RuntimeError("Something went wrong")

        # Calculate the mean age in dust cell, by weighing each stellar component age (old, young, ionizing) by the respective contribution to the global luminosity in that cell
        # DOESN'T WORK?? GIVES WAY TOO LOW VALUES
        #log_mean_ages = numbers.weighed_arithmetic_mean_numpy(self.log_stellar_ages, weights=self.cell_stellar_luminosity_contributions)
        log_mean_ages = self.log_old_stellar_age * self.cell_old_luminosity_contributions + self.log_young_stellar_age * self.cell_young_luminosity_contributions + self.log_sfr_stellar_age * self.cell_sfr_luminosity_contributions

        # Create the data with external xyz
        return Data3D.from_values(self.mean_age_name, log_mean_ages, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, unit=self.log_age_unit,
                                  description=self.mean_age_description, xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------

    @property
    def cell_mean_ages(self):
        return self.cell_mean_age_data.values

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
        return u("1/yr")

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "m51_ssfr_funev_path", True, write=False)
    def m51_scatter(self):

        """
        This function ...
        :return:
        """

        return Scatter2D.from_xy(self.m51_ssfr, self.m51_fractions, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_unit, x_description=self.ssfr_description, y_description=self.funev_description)

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

        return Scatter2D.from_xy(self.m31_ssfr, self.m31_fractions, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_unit, x_description=self.ssfr_description, y_description=self.funev_description)

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

    def create_ssfr_funev_pixel_scatter(self, name, ssfr_frame, funev_frame):

        """
        This function ...
        :param name:
        :param ssfr_frame:
        :param funev_frame:
        :return:
        """

        # Inform the user
        log.debug("Creating " + name + " sSFR - Funev pixel scatter ...")

        # Create and return
        return create_pixel_scatter(self.ssfr_name, self.funev_name, ssfr_frame, funev_frame, self.ssfr_description, self.funev_description, x_unit="1/yr")

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_cells_aux(self):
        return {sfr_density_name: self.valid_cell_sfr_densities_mappings_ke,
                dust_density_name: self.valid_cell_dust_densities_mappings_ke,
                distance_center_name: self.valid_cell_radii_mappings_ke,
                bulge_disk_ratio_name: self.valid_cell_bd_ratios_mappings_ke,
                temperature_name: self.valid_cell_temperatures_mappings_ke,
                mean_age_name: self.valid_cell_mean_ages_mappings_ke}

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_cells_aux_units(self):
        return {sfr_density_name: self.sfr_density_mappings_ke_unit,
                dust_density_name: self.cell_dust_density_unit,
                distance_center_name: self.length_unit,
                temperature_name: self.temperature_unit,
                mean_age_name: self.log_age_unit}

    # -----------------------------------------------------------------

    def create_ssfr_funev_cell_scatter(self, name, ssfr_data, funev_data, sfr_densities=None, sfr_density_unit=None):

        """
        This function ...
        :param name:
        :param ssfr_data:
        :param funev_data:
        :param sfr_densities:
        :param sfr_density_unit:
        :return:
        """

        # Debugging
        log.debug("Creating " + name + " sSFR - Funev cell scatter ...")

        # Set auxilary columns
        aux = OrderedDict()
        if sfr_densities is not None: aux[sfr_density_name] = sfr_densities
        aux[dust_density_name] = self.cell_dust_densities
        aux[distance_center_name] = self.cell_radii
        aux[radius_name] = self.cell_polar_radii
        aux[dust_heights_name] = self.cell_dust_heights
        aux[bulge_disk_ratio_name] = self.cell_bd_ratios
        aux[temperature_name] = self.cell_temperatures
        aux[luminosity_density_name] = self.cell_dust_luminosity_densities
        aux[mean_age_name] = self.cell_mean_ages

        # Set auxilary column units
        aux_units = OrderedDict()
        if sfr_densities is not None and sfr_density_unit is not None: aux_units[sfr_density_name] = sfr_density_unit
        aux_units[dust_density_name] = self.cell_dust_density_unit
        aux_units[distance_center_name] = self.length_unit
        aux_units[radius_name] = self.length_unit
        aux_units[temperature_name] = self.temperature_unit
        aux_units[luminosity_density_name] = self.cell_dust_luminosity_density_unit
        aux_units[mean_age_name] = self.log_age_unit

        # Create
        return create_cell_scatter(self.ssfr_name, self.funev_name, ssfr_data, funev_data, self.ssfr_description, self.funev_description, aux=aux, aux_units=aux_units, aux_is_arrays=True)

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
        #return fs.join(self.cell_sfr_path, "ssfr_salim.dat")
        return fs.join(self.cell_sfr_path, "ssfr_salim_corrected.dat")

    # -----------------------------------------------------------------

    @property
    def cell_sfr_salim_path(self):
        #return fs.join(self.cell_sfr_path, "sfr_salim.dat")
        return fs.join(self.cell_sfr_path, "sfr_salim_corrected.dat")

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

    @lazyproperty
    def valid_cell_temperatures_salim(self):
        return self.cell_temperatures[self.valid_cell_mask_salim]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mean_ages_salim(self):
        return self.cell_mean_ages[self.valid_cell_mask_salim]

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
        return sequences.contains_all(self.ssfr_salim_funev_cells_aux_colnames, ssfr_funev_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_ssfr_salim_funev_cells(self):
        if fs.is_file(self.ssfr_salim_funev_cells_path):
            if not self.ssfr_salim_funev_cells_has_all_aux_columns:
                #colnames = self.ssfr_salim_funev_cells_aux_colnames
                #if sfr_density_name not in colnames: self.ssfr_salim_funev_cells.add_aux(sfr_density_name, self.valid_cell_sfr_densities_salim, self.sfr_density_salim_unit, as_column=True)
                #if dust_density_name not in colnames: self.ssfr_salim_funev_cells.add_aux(dust_density_name, self.valid_cell_dust_densities_salim, self.cell_dust_density_unit, as_column=True)
                #if distance_center_name not in colnames: self.ssfr_salim_funev_cells.add_aux(distance_center_name, self.valid_cell_radii_salim, self.length_unit, as_column=True)
                #if bulge_disk_ratio_name not in colnames: self.ssfr_salim_funev_cells.add_aux(bulge_disk_ratio_name, self.valid_cell_bd_ratios_salim, as_column=True)
                #if temperature_name not in colnames: self.ssfr_salim_funev_cells.add_aux(temperature_name, self.valid_cell_temperatures_salim, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.ssfr_salim_funev_cells.add_aux(mean_age_name, self.valid_cell_mean_ages_salim, self.log_age_unit, as_column=True)
                #self.ssfr_salim_funev_cells.save() # save
                log.warning("The Salim sSFR - Funev cell scatter file does not contain all auxilary columns: removing and recalculating ...")
                fs.remove_file(self.ssfr_salim_funev_cells_path)
                return False
            return True
        else: return False

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

        # Create
        return self.create_ssfr_funev_cell_scatter("Salim", self.cell_ssfr_salim, self.cell_funev, sfr_densities=self.cell_sfr_densities_salim, sfr_density_unit=self.sfr_density_salim_unit)

    # -----------------------------------------------------------------
    #   K&E
    # -----------------------------------------------------------------

    @property
    def cell_ssfr_ke_path(self):
        #return fs.join(self.cell_sfr_path, "ssfr_ke.dat")
        return fs.join(self.cell_sfr_path, "ssfr_ke_corrected.dat")

    # -----------------------------------------------------------------

    @property
    def cell_sfr_ke_path(self):
        #return fs.join(self.cell_sfr_path, "sfr_ke.dat")
        return fs.join(self.cell_sfr_path, "sfr_ke_corrected.dat")

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

    @lazyproperty
    def valid_cell_temperatures_ke(self):
        return self.cell_temperatures[self.valid_cell_mask_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mean_ages_ke(self):
        return self.cell_mean_ages[self.valid_cell_mask_ke]

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
        return sequences.contains_all(self.ssfr_ke_funev_cells_aux_colnames, ssfr_funev_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_ssfr_ke_funev_cells(self):
        if fs.is_file(self.ssfr_ke_funev_cells_path):
            if not self.ssfr_ke_funev_cells_has_all_aux_columns:
                #colnames = self.ssfr_ke_funev_cells_aux_colnames
                #if sfr_density_name not in colnames: self.ssfr_ke_funev_cells.add_aux(sfr_density_name, self.valid_cell_sfr_densities_ke, self.sfr_density_ke_unit, as_column=True)
                #if dust_density_name not in colnames: self.ssfr_ke_funev_cells.add_aux(dust_density_name, self.valid_cell_dust_densities_ke, self.cell_dust_density_unit, as_column=True)
                #if distance_center_name not in colnames: self.ssfr_ke_funev_cells.add_aux(distance_center_name, self.valid_cell_radii_ke, self.length_unit, as_column=True)
                #if bulge_disk_ratio_name not in colnames: self.ssfr_ke_funev_cells.add_aux(bulge_disk_ratio_name, self.valid_cell_bd_ratios_ke, as_column=True)
                #if temperature_name not in colnames: self.ssfr_ke_funev_cells.add_aux(temperature_name, self.valid_cell_temperatures_ke, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.ssfr_ke_funev_cells.add_aux(mean_age_name, self.valid_cell_mean_ages_ke, self.log_age_unit, as_column=True)
                #self.ssfr_ke_funev_cells.save() # save
                log.warning("The K&E sSFR - Funev cell scatter file does not contain all axuilary columns: removing and recalculating ...")
                fs.remove_file(self.ssfr_ke_funev_cells_path)
                return False
            return True
        else: return False

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

        # Create
        return self.create_ssfr_funev_cell_scatter("K&E", self.cell_ssfr_ke, self.cell_funev, sfr_densities=self.cell_sfr_densities_ke, sfr_density_unit=self.sfr_density_ke_unit)

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

    @lazyproperty
    def valid_cell_temperatures_mappings(self):
        return self.cell_temperatures[self.valid_cell_mask_mappings]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mean_ages_mappings(self):
        return self.cell_mean_ages[self.valid_cell_mask_mappings]

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
        return sequences.contains_all(self.ssfr_mappings_funev_cells_aux_colnames, ssfr_funev_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_funev_cells(self):
        if fs.is_file(self.ssfr_mappings_funev_cells_path):
            if not self.ssfr_mappings_funev_cells_has_all_aux_columns:
                #colnames = self.ssfr_mappings_funev_cells_aux_colnames
                #if sfr_density_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(sfr_density_name, self.valid_cell_sfr_densities_mappings, self.sfr_density_mappings_unit, as_column=True)
                #if dust_density_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(dust_density_name, self.valid_cell_dust_densities_mappings, self.cell_dust_density_unit, as_column=True)
                #if distance_center_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(distance_center_name, self.valid_cell_radii_mappings, self.length_unit, as_column=True)
                #if bulge_disk_ratio_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(bulge_disk_ratio_name, self.valid_cell_bd_ratios_mappings, as_column=True)
                #if temperature_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(temperature_name, self.valid_cell_temperatures_mappings, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.ssfr_mappings_funev_cells.add_aux(mean_age_name, self.valid_cell_mean_ages_mappings, self.log_age_unit, as_column=True)
                #self.ssfr_mappings_funev_cells.save() # save
                log.warning("The MAPPINGS sSFR - Funev cell scatter file does not contain all axuilary columns: removing and recalculating ...")
                fs.remove_file(self.ssfr_mappings_funev_cells_path)
                return False
            return True
        else: return False

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

        # Create
        return self.create_ssfr_funev_cell_scatter("MAPPINGS", self.cell_ssfr_mappings, self.cell_funev, sfr_densities=self.cell_sfr_densities_mappings, sfr_density_unit=self.sfr_density_mappings_unit)

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

    @lazyproperty
    def valid_cell_temperatures_mappings_ke(self):
        return self.cell_temperatures[self.valid_cell_mask_mappings_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mean_ages_mappings_ke(self):
        return self.cell_mean_ages[self.valid_cell_mask_mappings_ke]

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
        return sequences.contains_all(self.ssfr_mappings_ke_funev_cells_aux_colnames, ssfr_funev_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_ke_funev_cells(self):
        if fs.is_file(self.ssfr_mappings_ke_funev_cells_path):
            if not self.ssfr_mappings_ke_funev_cells_has_all_aux_columns:
                #colnames = self.ssfr_mappings_ke_funev_cells_aux_colnames
                #if sfr_density_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(sfr_density_name, self.valid_cell_sfr_densities_mappings_ke, self.sfr_density_mappings_ke_unit, as_column=True)
                #if dust_density_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(dust_density_name, self.valid_cell_dust_densities_mappings_ke, self.cell_dust_density_unit, as_column=True)
                #if distance_center_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(distance_center_name, self.valid_cell_radii_mappings_ke, self.length_unit, as_column=True)
                #if bulge_disk_ratio_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(bulge_disk_ratio_name, self.valid_cell_bd_ratios_mappings_ke, as_column=True)
                #if temperature_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(temperature_name, self.valid_cell_temperatures_mappings_ke, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.ssfr_mappings_ke_funev_cells.add_aux(mean_age_name, self.valid_cell_mean_ages_mappings_ke, self.log_age_unit, as_column=True)
                #self.ssfr_mappings_ke_funev_cells.save() # save
                log.warning("The MAPPINGS + K&E sSFR - Funev cell scatter file does not contain all axuilary columns: removing and recalculating ...")
                fs.remove_file(self.ssfr_mappings_ke_funev_cells_path)
                return False
            return True
        else: return False

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

        # Create
        return self.create_ssfr_funev_cell_scatter("MAPPINGS + K&E", self.cell_ssfr_mappings_ke, self.cell_funev, sfr_densities=self.cell_sfr_densities_mappings_ke, sfr_density_unit=self.sfr_density_mappings_ke_unit)

    # -----------------------------------------------------------------
    # sSFR-Funev pixel scatter data
    # -----------------------------------------------------------------

    @property
    def projected_heating_maps_path(self):
        return fs.join(self.projected_heating_path, "maps")

    # -----------------------------------------------------------------

    @property
    def pixel_funev_path(self):
        # is there a reason this was earth_diffuse first?
        #return fs.join(self.projected_heating_maps_path, "earth_diffuse.fits")
        return fs.join(self.projected_heating_maps_path, "earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_funev(self):
        return fs.is_file(self.pixel_funev_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_funev(self):
        return Frame.from_file(self.pixel_funev_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_funev_values(self):
        return self.pixel_funev.values

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_salim_path(self):
        #return fs.join(self.projected_sfr_path, "ssfr_salim_earth.fits")
        return fs.join(self.projected_sfr_path, "ssfr_salim_earth_corrected.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_salim(self):
        return fs.is_file(self.pixel_ssfr_salim_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_salim(self):
        return Frame.from_file(self.pixel_ssfr_salim_path)

    # -----------------------------------------------------------------

    @lazyproperty
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

        # Create scatter
        return self.create_ssfr_funev_pixel_scatter("Salim", self.pixel_ssfr_salim, self.pixel_funev)

    # -----------------------------------------------------------------
    #   K&E
    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_ke_path(self):
        #return fs.join(self.projected_sfr_path, "ssfr_ke_earth.fits")
        return fs.join(self.projected_sfr_path, "ssfr_ke_earth_corrected.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_ke(self):
        return fs.is_file(self.pixel_ssfr_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_ke(self):
        return Frame.from_file(self.pixel_ssfr_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
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
        return self.create_ssfr_funev_pixel_scatter("K&E", self.pixel_ssfr_ke, self.pixel_funev)

    # -----------------------------------------------------------------
    #   MAPPINGS
    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_mappings_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_mappings_earth.fits")

    # -----------------------------------------------------------------

    @property
    def pixel_sfr_mappings_path(self):
        return fs.join(self.projected_sfr_path, "sfr_mappings_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_mappings(self):
        return fs.is_file(self.pixel_ssfr_mappings_path)

    # -----------------------------------------------------------------

    @property
    def has_pixel_sfr_mappings(self):
        return fs.is_file(self.pixel_sfr_mappings_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_mappings(self):
        return Frame.from_file(self.pixel_ssfr_mappings_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_sfr_mappings(self):
        return Frame.from_file(self.pixel_sfr_mappings_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_mappings_values(self):
        return self.pixel_ssfr_mappings.values

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_sfr_mappings_values(self):
        return self.pixel_sfr_mappings.values

    # -----------------------------------------------------------------

    @property
    def pixel_sfr_mappings_pixelarea(self):
        return self.pixel_sfr_mappings.pixelarea

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
        return self.create_ssfr_funev_pixel_scatter("MAPPINGS", self.pixel_ssfr_mappings, self.pixel_funev)

    # -----------------------------------------------------------------
    #   MAPPINGS + K&E
    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_mappings_ke_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_mappings_ke_earth.fits")

    # -----------------------------------------------------------------

    @property
    def pixel_sfr_mappings_ke_path(self):
        return fs.join(self.projected_sfr_path, "sfr_mappings_ke_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_mappings_ke(self):
        return fs.is_file(self.pixel_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @property
    def has_pixel_sfr_mappings_ke(self):
        return fs.is_file(self.pixel_sfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_mappings_ke(self):
        return Frame.from_file(self.pixel_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_sfr_mappings_ke(self):
        return Frame.from_file(self.pixel_sfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_mappings_ke_values(self):
        return self.pixel_ssfr_mappings_ke.values

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_sfr_mappings_ke_values(self):
        return self.pixel_sfr_mappings_ke.values

    # -----------------------------------------------------------------

    @property
    def pixel_sfr_mappings_ke_pixelarea(self):
        return self.pixel_sfr_mappings_ke.pixelarea

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
        return self.create_ssfr_funev_pixel_scatter("MAPPINGS + K&E", self.pixel_ssfr_mappings_ke, self.pixel_funev)

    # -----------------------------------------------------------------
    # FUV-H / Funev
    # -----------------------------------------------------------------

    @property
    def fuv_h_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells_fuv_h.dat")

    # -----------------------------------------------------------------

    @property
    def fuv_h_funev_cells_colnames(self):
        return fs.get_column_names(self.fuv_h_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def fuv_h_funev_base_colnames(self):
        return [self.fuv_h_ssfr_name, self.funev_name]

    # -----------------------------------------------------------------

    @property
    def fuv_h_funev_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.fuv_h_funev_cells_colnames, self.fuv_h_funev_base_colnames)

    # -----------------------------------------------------------------

    @property
    def fuv_h_funev_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.fuv_h_funev_cells_aux_colnames, colour_funev_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_fuv_h_funev_cells(self):
        if fs.is_file(self.fuv_h_funev_cells_path):
            if not self.fuv_h_funev_cells_has_all_aux_columns:
                #colnames = self.fuv_h_funev_cells_aux_colnames
                #if sfr_density_salim_name not in colnames: self.fuv_h_funev_cells.add_aux(sfr_density_salim_name, self.cell_sfr_densities_salim, self.sfr_density_salim_unit, as_column=True)
                #if sfr_density_ke_name not in colnames: self.fuv_h_funev_cells.add_aux(sfr_density_ke_name, self.cell_sfr_densities_ke, self.sfr_density_ke_unit, as_column=True)
                #if sfr_density_mappings_name not in colnames: self.fuv_h_funev_cells.add_aux(sfr_density_mappings_name, self.cell_sfr_densities_mappings, self.sfr_density_mappings_unit, as_column=True)
                #if sfr_density_mappings_ke_name not in colnames: self.fuv_h_funev_cells.add_aux(sfr_density_mappings_ke_name, self.cell_sfr_densities_mappings_ke, self.sfr_density_mappings_ke_unit, as_column=True)
                #if dust_density_name not in colnames: self.fuv_h_funev_cells.add_aux(dust_density_name, self.cell_dust_densities, self.cell_dust_density_unit, as_column=True)
                #if distance_center_name not in colnames: self.fuv_h_funev_cells.add_aux(distance_center_name, self.cell_radii, self.length_unit, as_column=True)
                #if bulge_disk_ratio_name not in colnames: self.fuv_h_funev_cells.add_aux(bulge_disk_ratio_name, self.cell_bd_ratios, as_column=True)
                #if temperature_name not in colnames: self.fuv_h_funev_cells.add_aux(temperature_name, self.cell_temperatures, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.fuv_h_funev_cells.add_aux(mean_age_name, self.cell_mean_ages, self.log_age_unit, as_column=True)
                #self.fuv_h_funev_cells.save() # save
                fs.remove_file(self.fuv_h_funev_cells_path)
                return False
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyproperty
    def colour_funev_cells_aux(self):
        return {sfr_density_salim_name: self.cell_sfr_densities_salim,
                sfr_density_ke_name: self.cell_sfr_densities_ke,
                sfr_density_mappings_name: self.cell_sfr_densities_mappings,
                sfr_density_mappings_ke_name: self.cell_sfr_densities_mappings_ke,
                dust_density_name: self.cell_dust_densities,
                distance_center_name: self.cell_radii,
                bulge_disk_ratio_name: self.cell_bd_ratios,
                temperature_name: self.cell_temperatures,
                mean_age_name: self.cell_mean_ages}

    # -----------------------------------------------------------------

    @lazyproperty
    def colour_funev_cells_aux_units(self):
        return {sfr_density_salim_name: self.sfr_density_salim_unit,
                sfr_density_ke_name: self.sfr_density_ke_unit,
                sfr_density_mappings_name: self.sfr_density_mappings_unit,
                sfr_density_mappings_ke_name: self.sfr_density_mappings_ke_unit,
                dust_density_name: self.cell_dust_density_unit,
                distance_center_name: self.length_unit,
                temperature_name: self.temperature_unit,
                mean_age_name: self.log_age_unit}

    # -----------------------------------------------------------------

    def create_colour_funev_cell_scatter(self, colour_name, colour_data, colour_description):

        """
        This function ...
        :param colour_name:
        :param colour_data:
        :param colour_description:
        :return:
        """

        # Debugging
        log.debug("Creating " + colour_name + " - Funev cell scatter ...")

        # Create
        return create_cell_scatter(colour_name, self.funev_name, colour_data, self.cell_funev, colour_description, self.funev_description,
                                   aux=self.colour_funev_cells_aux, aux_units=self.colour_funev_cells_aux_units, aux_is_arrays=True,
                                   x_unit=self.magnitude_unit)

    # -----------------------------------------------------------------

    def create_colour_funev_pixel_scatter(self, colour_name, ssfr_frame, colour_description):

        """
        This function ...
        :param colour_name:
        :param ssfr_frame:
        :param colour_description:
        :return:
        """

        # Debugging
        log.debug("Creating " + colour_name + " - Funev pixel scatter ...")

        # Create
        return create_pixel_scatter(colour_name, self.funev_name, ssfr_frame, self.pixel_funev, colour_description, self.funev_description,
                                    x_unit=self.magnitude_unit)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "fuv_h_funev_cells_path", True, write=False)
    def fuv_h_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Create scatter
        return self.create_colour_funev_cell_scatter(self.fuv_h_ssfr_name, self.cell_ssfr_fuv_h, self.fuv_h_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def fuv_h_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_fuv_h.dat")

    # -----------------------------------------------------------------

    @property
    def has_fuv_h_funev_pixels(self):
        return fs.is_file(self.fuv_h_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "fuv_h_funev_pixels_path", True, write=False)
    def fuv_h_funev_pixels(self):

        """
        Thisn function ...
        :return:
        """

        # Checks
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Create scatter
        return self.create_colour_funev_pixel_scatter(self.fuv_h_ssfr_name, self.pixel_ssfr_fuv_h, self.fuv_h_ssfr_description)

    # -----------------------------------------------------------------
    # FUV-R / Funev
    # -----------------------------------------------------------------

    @property
    def fuv_r_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells_fuv_r.dat")

    # -----------------------------------------------------------------

    @property
    def fuv_r_funev_cells_colnames(self):
        return fs.get_column_names(self.fuv_r_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def fuv_r_funev_base_colnames(self):
        return [self.fuv_r_ssfr_name, self.funev_name]

    # -----------------------------------------------------------------

    @property
    def fuv_r_funev_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.fuv_r_funev_cells_colnames, self.fuv_r_funev_base_colnames)

    # -----------------------------------------------------------------

    @property
    def fuv_r_funev_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.fuv_r_funev_cells_aux_colnames, colour_funev_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_fuv_r_funev_cells(self):
        if fs.is_file(self.fuv_r_funev_cells_path):
            if not self.fuv_r_funev_cells_has_all_aux_columns:
                #colnames = self.fuv_r_funev_cells_aux_colnames
                #if sfr_density_salim_name not in colnames: self.fuv_r_funev_cells.add_aux(sfr_density_salim_name, self.cell_sfr_densities_salim, self.sfr_density_salim_unit, as_column=True)
                #if sfr_density_ke_name not in colnames: self.fuv_r_funev_cells.add_aux(sfr_density_ke_name, self.cell_sfr_densities_ke, self.sfr_density_ke_unit, as_column=True)
                #if sfr_density_mappings_name not in colnames: self.fuv_r_funev_cells.add_aux(sfr_density_mappings_name, self.cell_sfr_densities_mappings, self.sfr_density_mappings_unit, as_column=True)
                #if sfr_density_mappings_ke_name not in colnames: self.fuv_r_funev_cells.add_aux(sfr_density_mappings_ke_name, self.cell_sfr_densities_mappings_ke, self.sfr_density_mappings_ke_unit, as_column=True)
                #if dust_density_name not in colnames: self.fuv_r_funev_cells.add_aux(dust_density_name, self.cell_dust_densities, self.cell_dust_density_unit, as_column=True)
                #if distance_center_name not in colnames: self.fuv_r_funev_cells.add_aux(distance_center_name, self.cell_radii, self.length_unit, as_column=True)
                #if bulge_disk_ratio_name not in colnames: self.fuv_r_funev_cells.add_aux(bulge_disk_ratio_name, self.cell_bd_ratios, as_column=True)
                #if temperature_name not in colnames: self.fuv_r_funev_cells.add_aux(temperature_name, self.cell_temperatures, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.fuv_r_funev_cells.add_aux(mean_age_name, self.cell_mean_ages, self.log_age_unit, as_column=True)
                #self.fuv_r_funev_cells.save() # save
                fs.remove_file(self.fuv_r_funev_cells_path)
                return False
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "fuv_r_funev_cells_path", True, write=False)
    def fuv_r_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Create scatter
        return self.create_colour_funev_cell_scatter(self.fuv_r_ssfr_name, self.cell_ssfr_fuv_r, self.fuv_r_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def fuv_r_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_fuv_r.dat")

    # -----------------------------------------------------------------

    @property
    def has_fuv_r_funev_pixels(self):
        return fs.is_file(self.fuv_r_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "fuv_r_funev_pixels_path", True, write=False)
    def fuv_r_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Create scatter
        return self.create_colour_funev_pixel_scatter(self.fuv_r_ssfr_name, self.pixel_ssfr_fuv_r, self.fuv_r_ssfr_description)

    # -----------------------------------------------------------------
    # NUV-H / Funev
    # -----------------------------------------------------------------

    @property
    def nuv_h_ssfr_name(self):
        return "FUV-H"

    # -----------------------------------------------------------------

    @property
    def nuv_h_ssfr_description(self):
        return "NUV-H colour as specific star formation rate proxy"

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_nuv_h_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_nuv_h.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_nuv_h(self):
        return fs.is_file(self.cell_ssfr_nuv_h_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_nuv_h(self):
        return Data3D.from_file(self.cell_ssfr_nuv_h_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_nuv_h_values(self):
        return self.cell_ssfr_nuv_h.values

    # -----------------------------------------------------------------

    @property
    def nuv_h_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells_nuv_h.dat")

    # -----------------------------------------------------------------

    @property
    def nuv_h_funev_cells_colnames(self):
        return fs.get_column_names(self.nuv_h_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def nuv_h_funev_base_colnames(self):
        return [self.nuv_h_ssfr_name, self.funev_name]

    # -----------------------------------------------------------------

    @property
    def nuv_h_funev_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.nuv_h_funev_cells_colnames, self.nuv_h_funev_base_colnames)

    # -----------------------------------------------------------------

    @property
    def nuv_h_funev_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.nuv_h_funev_cells_aux_colnames, colour_funev_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_nuv_h_funev_cells(self):
        if fs.is_file(self.nuv_h_funev_cells_path):
            if not self.nuv_h_funev_cells_has_all_aux_columns:
                #colnames = self.nuv_h_funev_cells_aux_colnames
                #if sfr_density_salim_name not in colnames: self.nuv_h_funev_cells.add_aux(sfr_density_salim_name, self.cell_sfr_densities_salim, self.sfr_density_salim_unit, as_column=True)
                #if sfr_density_ke_name not in colnames: self.nuv_h_funev_cells.add_aux(sfr_density_ke_name, self.cell_sfr_densities_ke, self.sfr_density_ke_unit, as_column=True)
                #if sfr_density_mappings_name not in colnames: self.nuv_h_funev_cells.add_aux(sfr_density_mappings_name, self.cell_sfr_densities_mappings, self.sfr_density_mappings_unit, as_column=True)
                #if sfr_density_mappings_ke_name not in colnames: self.nuv_h_funev_cells.add_aux(sfr_density_mappings_ke_name, self.cell_sfr_densities_mappings_ke, self.sfr_density_mappings_ke_unit, as_column=True)
                #if dust_density_name not in colnames: self.nuv_h_funev_cells.add_aux(dust_density_name, self.cell_dust_densities, self.cell_dust_density_unit, as_column=True)
                #if distance_center_name not in colnames: self.nuv_h_funev_cells.add_aux(distance_center_name, self.cell_radii, self.length_unit, as_column=True)
                #if bulge_disk_ratio_name not in colnames: self.nuv_h_funev_cells.add_aux(bulge_disk_ratio_name, self.cell_bd_ratios, as_column=True)
                #if temperature_name not in colnames: self.nuv_h_funev_cells.add_aux(temperature_name, self.cell_temperatures, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.nuv_h_funev_cells.add_aux(mean_age_name, self.cell_mean_ages, self.log_age_unit, as_column=True)
                #self.nuv_h_funev_cells.save() # save
                fs.remove_file(self.nuv_h_funev_cells_path)
                return False
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "nuv_h_funev_cells_path", True, write=False)
    def nuv_h_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Create scatter
        return self.create_colour_funev_cell_scatter(self.nuv_h_ssfr_name, self.cell_ssfr_nuv_h, self.nuv_h_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_nuv_h_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_nuv_h_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_nuv_h(self):
        return fs.is_file(self.pixel_ssfr_nuv_h_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_nuv_h(self):
        return Frame.from_file(self.pixel_ssfr_nuv_h_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_nuv_h_values(self):
        return self.pixel_ssfr_nuv_h.values

    # -----------------------------------------------------------------

    @property
    def nuv_h_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_nuv_h.dat")

    # -----------------------------------------------------------------

    @property
    def has_nuv_h_funev_pixels(self):
        return fs.is_file(self.nuv_h_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "nuv_h_funev_pixels_path", True, write=False)
    def nuv_h_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Create scatter
        return self.create_colour_funev_pixel_scatter(self.nuv_h_ssfr_name, self.pixel_ssfr_nuv_h, self.nuv_h_ssfr_description)

    # -----------------------------------------------------------------
    # NUV-R / Funev
    # -----------------------------------------------------------------

    @property
    def nuv_r_ssfr_name(self):
        return "NUV-r"

    # -----------------------------------------------------------------

    @property
    def nuv_r_ssfr_description(self):
        return "NUV-r colour as specific star formation rate proxy"

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_nuv_r_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_nuv_r.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_nuv_r(self):
        return fs.is_file(self.cell_ssfr_nuv_r_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_nuv_r(self):
        return Data3D.from_file(self.cell_ssfr_nuv_r_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_nuv_r_values(self):
        return self.cell_ssfr_nuv_r.values

    # -----------------------------------------------------------------

    @property
    def nuv_r_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells_nuv_r.dat")

    # -----------------------------------------------------------------

    @property
    def nuv_r_funev_cells_colnames(self):
        return fs.get_column_names(self.nuv_r_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def nuv_r_funev_base_colnames(self):
        return [self.nuv_r_ssfr_name, self.funev_name]

    # -----------------------------------------------------------------

    @property
    def nuv_r_funev_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.nuv_r_funev_cells_colnames, self.nuv_r_funev_base_colnames)

    # -----------------------------------------------------------------

    @property
    def nuv_r_funev_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.nuv_r_funev_cells_aux_colnames, colour_funev_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_nuv_r_funev_cells(self):
        if fs.is_file(self.nuv_r_funev_cells_path):
            if not self.nuv_r_funev_cells_has_all_aux_columns:
                #colnames = self.nuv_r_funev_cells_aux_colnames
                #if sfr_density_salim_name not in colnames: self.nuv_r_funev_cells.add_aux(sfr_density_salim_name, self.cell_sfr_densities_salim, self.sfr_density_salim_unit, as_column=True)
                #if sfr_density_ke_name not in colnames: self.nuv_r_funev_cells.add_aux(sfr_density_ke_name, self.cell_sfr_densities_ke, self.sfr_density_ke_unit, as_column=True)
                #if sfr_density_mappings_name not in colnames: self.nuv_r_funev_cells.add_aux(sfr_density_mappings_name, self.cell_sfr_densities_mappings, self.sfr_density_mappings_unit, as_column=True)
                #if sfr_density_mappings_ke_name not in colnames: self.nuv_r_funev_cells.add_aux(sfr_density_mappings_ke_name, self.cell_sfr_densities_mappings_ke, self.sfr_density_mappings_ke_unit, as_column=True)
                #if dust_density_name not in colnames: self.nuv_r_funev_cells.add_aux(dust_density_name, self.cell_dust_densities, self.cell_dust_density_unit, as_column=True)
                #if distance_center_name not in colnames: self.nuv_r_funev_cells.add_aux(distance_center_name, self.cell_radii, self.length_unit, as_column=True)
                #if bulge_disk_ratio_name not in colnames: self.nuv_r_funev_cells.add_aux(bulge_disk_ratio_name, self.cell_bd_ratios, as_column=True)
                #if temperature_name not in colnames: self.nuv_r_funev_cells.add_aux(temperature_name, self.cell_temperatures, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.nuv_r_funev_cells.add_aux(mean_age_name, self.cell_mean_ages, self.log_age_unit, as_column=True)
                #self.nuv_r_funev_cells.save() # save
                fs.remove_file(self.nuv_r_funev_cells_path)
                return False
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "nuv_r_funev_cells_path", True, write=False)
    def nuv_r_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Create scatter
        return self.create_colour_funev_cell_scatter(self.nuv_r_ssfr_name, self.cell_ssfr_nuv_r, self.nuv_r_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_nuv_r_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_nuv_r_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_nuv_r(self):
        return fs.is_file(self.pixel_ssfr_nuv_r_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_nuv_r(self):
        return Frame.from_file(self.pixel_ssfr_nuv_r_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_nuv_r_values(self):
        return self.pixel_ssfr_nuv_r.values

    # -----------------------------------------------------------------

    @property
    def nuv_r_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_nuv_r.dat")

    # -----------------------------------------------------------------

    @property
    def has_nuv_r_funev_pixels(self):
        return fs.is_file(self.nuv_r_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "nuv_r_funev_pixels_path", True, write=False)
    def nuv_r_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Create scatter
        return self.create_colour_funev_pixel_scatter(self.nuv_r_ssfr_name, self.pixel_ssfr_nuv_r, self.nuv_r_ssfr_description)

    # -----------------------------------------------------------------
    # TIR SFR
    # -----------------------------------------------------------------

    @property
    def cell_ssfr_tir_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_tir.dat")

    # -----------------------------------------------------------------

    @property
    def cell_sfr_tir_path(self):
        return fs.join(self.cell_sfr_path, "sfr_tir.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_tir(self):
        return fs.is_file(self.cell_ssfr_tir_path)

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_tir(self):
        return fs.is_file(self.cell_sfr_tir_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_tir(self):
        return Data3D.from_file(self.cell_ssfr_tir_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_tir(self):
        return Data3D.from_file(self.cell_sfr_tir_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_tir_values(self):
        return self.cell_ssfr_tir.values

    # -----------------------------------------------------------------

    @property
    def ssfr_tir_units(self):
        return self.cell_ssfr_tir.unit

    # -----------------------------------------------------------------

    @property
    def cell_sfr_tir_values(self):
        return self.cell_sfr_tir.values

    # -----------------------------------------------------------------

    @property
    def sfr_tir_unit(self):
        return self.cell_sfr_tir.unit

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_tir_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_tir_earth.fits")

    # -----------------------------------------------------------------

    @property
    def pixel_sfr_tir_path(self):
        return fs.join(self.projected_sfr_path, "sfr_tir_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_tir(self):
        return fs.is_file(self.pixel_ssfr_tir_path)

    # -----------------------------------------------------------------

    @property
    def has_pixel_sfr_tir(self):
        return fs.is_file(self.pixel_sfr_tir_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_tir(self):
        return Frame.from_file(self.pixel_ssfr_tir_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_sfr_tir(self):
        return Frame.from_file(self.pixel_sfr_tir_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_tir_values(self):
        return self.pixel_ssfr_tir.values

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_sfr_tir_values(self):
        return self.pixel_sfr_tir.values

    # -----------------------------------------------------------------

    @property
    def pixel_sfr_tir_pixelarea(self):
        return self.pixel_sfr_tir.pixelarea

    # -----------------------------------------------------------------
    # 24 micron SFR
    # -----------------------------------------------------------------

    @property
    def cell_ssfr_24um_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_24um.dat")

    # -----------------------------------------------------------------

    @property
    def cell_sfr_24um_path(self):
        return fs.join(self.cell_sfr_path, "sfr_24um.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_24um(self):
        return fs.is_file(self.cell_ssfr_24um_path)

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_24um(self):
        return fs.is_file(self.cell_sfr_24um_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_24um(self):
        return Data3D.from_file(self.cell_ssfr_24um_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_24um(self):
        return Data3D.from_file(self.cell_sfr_24um_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_24um_values(self):
        return self.cell_ssfr_24um.values

    # -----------------------------------------------------------------

    @property
    def cell_sfr_24um_values(self):
        return self.cell_sfr_24um.values

    # -----------------------------------------------------------------

    @property
    def ssfr_24um_unit(self):
        return self.cell_ssfr_24um.unit

    # -----------------------------------------------------------------

    @property
    def sfr_24um_unit(self):
        return self.cell_sfr_24um.unit

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_24um_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_24um_earth.fits")

    # -----------------------------------------------------------------

    @property
    def pixel_sfr_24um_path(self):
        return fs.join(self.projected_sfr_path, "sfr_24um_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_24um(self):
        return fs.is_file(self.pixel_ssfr_24um_path)

    # -----------------------------------------------------------------

    @property
    def has_pixel_sfr_24um(self):
        return fs.is_file(self.pixel_sfr_24um_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_24um(self):
        return Frame.from_file(self.pixel_ssfr_24um_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_sfr_24um(self):
        return Frame.from_file(self.pixel_sfr_24um_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_24um_values(self):
        return self.pixel_ssfr_24um.values

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_sfr_24um_values(self):
        return self.pixel_sfr_24um.values

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_sfr_24um_pixelarea(self):
        return self.pixel_sfr_24um.pixelarea

    # -----------------------------------------------------------------
    # Temperature - Funev
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
    # MEAN AGE - FUNEV
    # -----------------------------------------------------------------

    @property
    def has_mean_age_cells(self):
        return True

    # -----------------------------------------------------------------

    @property
    def has_mean_age_pixels(self):
        return False # not implemented yet

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_age_funev_path(self):
        return fs.create_directory_in(self.correlations_path, mean_age_funev_name)

    # -----------------------------------------------------------------

    @property
    def mean_age_funev_cells_path(self):
        return fs.join(self.mean_age_funev_path, "cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_mean_age_funev_cells(self):
        return fs.is_file(self.mean_age_funev_cells_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mean_age_mask(self):
        return self.cell_mean_ages != 0

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mask_mean_age_funev(self):
        return self.valid_cell_mean_age_mask * self.valid_cell_funev_mask

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mean_age_funev_cells_path", True, write=False)
    def mean_age_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Get values
        valid_cell_mean_ages = self.cell_mean_ages[self.valid_cell_mask_mean_age_funev]
        valid_cell_funev_values = self.cell_funev_values[self.valid_cell_mask_mean_age_funev]

        # Create and return
        return Scatter2D.from_xy(valid_cell_mean_ages, valid_cell_funev_values,
                                 x_name=self.mean_age_name, y_name=self.funev_name, x_unit=self.log_age_unit,
                                 x_description=self.mean_age_description, y_description=self.funev_description)

    # -----------------------------------------------------------------

    @property
    def mean_age_funev_pixels_path(self):
        return fs.join(self.mean_age_funev_path, "pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mean_age_funev_pixels(self):
        return fs.is_file(self.mean_age_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mean_age_funev_pixels_path", True, write=False)
    def mean_age_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Get values
        valid_pixel_mean_ages = None
        valid_pixel_funev_values = None

        # Return
        return Scatter2D.from_xy(valid_pixel_mean_ages, valid_pixel_funev_values,
                                 x_name=self.mean_age_name, y_name=self.funev_name, x_unit=self.log_age_unit,
                                 x_description=self.mean_age_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    # MEAN AGE - sSFR
    # -----------------------------------------------------------------

    @lazyproperty
    def mean_age_ssfr_path(self):
        return fs.create_directory_in(self.correlations_path, mean_age_ssfr_name)

    # -----------------------------------------------------------------

    @property
    def mean_age_ssfr_mappings_cells_path(self):
        return fs.join(self.mean_age_ssfr_path, "mappings_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_mean_age_ssfr_mappings_cells(self):
        return fs.is_file(self.mean_age_ssfr_mappings_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mean_age_ssfr_mappings_cells_path", True, write=False)
    def mean_age_ssfr_mappings_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_ssfr_mappings: raise IOError("The cell sSFR MAPPINGS data is not present: run the sfr analysis first")

        # Get values
        mean_ages = self.cell_mean_ages
        ssfr = self.cell_ssfr_mappings_values

        # Create and return
        return Scatter2D.from_xy(mean_ages, ssfr, x_name=self.mean_age_name, y_name=self.mappings_ssfr_name,
                                 x_unit=self.log_age_unit, y_unit=self.ssfr_mappings_unit,
                                 x_description=self.mean_age_description, y_description=self.mappings_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mean_age_ssfr_mappings_pixels_path(self):
        return fs.join(self.mean_age_ssfr_path, "mappings_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mean_age_ssfr_mappings_pixels(self):
        return fs.is_file(self.mean_age_ssfr_mappings_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mean_age_ssfr_mappings_pixels_path", True, write=False)
    def mean_age_ssfr_mappings_pixels(self):

        """
        This function ...
        :return:
        """

        # Get values
        mean_ages = None
        ssfr = None

        return None

    # -----------------------------------------------------------------

    @property
    def mean_age_ssfr_mappings_ke_cells_path(self):
        return fs.join(self.mean_age_ssfr_path, "mappings_ke_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_mean_age_ssfr_mappings_ke_cells(self):
        return fs.is_file(self.mean_age_ssfr_mappings_ke_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mean_age_ssfr_mappings_ke_cells_path", True, write=False)
    def mean_age_ssfr_mappings_ke_cells(self):

        """
        Thisf unction ...
        :return:
        """

        # Checks
        if not self.has_cell_ssfr_mappings_ke: raise IOError("The cell sSFR MAPPINGS K&E data is not present: run the sfr analysis first")

        # Get values
        mean_ages = self.cell_mean_ages
        ssfr = self.cell_ssfr_mappings_ke_values

        # Create and return
        return Scatter2D.from_xy(mean_ages, ssfr, x_name=self.mean_age_name, y_name=self.mappings_ke_ssfr_name,
                                 x_unit=self.log_age_unit, y_unit=self.ssfr_mappings_ke_unit,
                                 x_description=self.mean_age_description, y_description=self.mappings_ke_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mean_age_ssfr_mappings_ke_pixels_path(self):
        return fs.join(self.mean_age_ssfr_path, "mappings_ke_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mean_age_ssfr_mappings_ke_pixels(self):
        return fs.is_file(self.mean_age_ssfr_mappings_ke_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mean_age_ssfr_mappings_ke_pixels_path", True, write=False)
    def mean_age_ssfr_mappings_ke_pixels(self):

        """
        This function ...
        :return:
        """

        # Get values
        mean_ages = None
        ssfr = None

        return None

    # -----------------------------------------------------------------
    # MAPPINGS SFR - TIR SFR
    # -----------------------------------------------------------------

    @property
    def mappings_sfr_name(self):
        return "SFR_MAPPINGS"

    # -----------------------------------------------------------------

    @property
    def mappings_sfr_description(self):
        return "Star formation rate (MAPPINGS)"

    # -----------------------------------------------------------------

    @property
    def tir_sfr_name(self):
        return "SFR_TIR"

    # -----------------------------------------------------------------

    @property
    def tir_sfr_description(self):
        return "Star formation rate (TIR)"

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_sfr_path(self):
        return fs.create_directory_in(self.correlations_path, "SFR-SFR")

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_densities_tir(self):
        return self.cell_sfr_tir_values / self.cell_volumes

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_density_tir_unit(self):
        return self.sfr_tir_unit / self.volume_unit

    # -----------------------------------------------------------------

    @property
    def mappings_sfr_density_name(self):
        return "vSFR_MAPPINGS"

    # -----------------------------------------------------------------

    @property
    def mappings_sfr_density_description(self):
        return "Star formation rate density (MAPPINGS)"

    # -----------------------------------------------------------------

    @property
    def tir_sfr_density_name(self):
        return "vSFR_TIR"

    # -----------------------------------------------------------------

    @property
    def tir_sfr_density_description(self):
        return "Star formation rate density (TIR)"

    # -----------------------------------------------------------------

    @property
    def mappings_tir_sfr_cells_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_tir_cells.dat")

    # -----------------------------------------------------------------

    @property
    def mappings_tir_sfr_cells_colnames(self):
        return fs.get_column_names(self.mappings_tir_sfr_cells_path)

    # -----------------------------------------------------------------

    @property
    def mappings_tir_sfr_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.mappings_tir_sfr_cells_colnames, [self.mappings_sfr_density_name, self.tir_sfr_density_name])

    # -----------------------------------------------------------------

    @property
    def mappings_tir_sfr_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.mappings_tir_sfr_cells_aux_colnames, sfr_sfr_cells_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_mappings_tir_sfr_cells(self):
        if fs.is_file(self.mappings_tir_sfr_cells_path):
            if not self.mappings_tir_sfr_cells_has_all_aux_columns:
                #colnames = self.mappings_tir_sfr_cells_aux_colnames
                #if temperature_name not in colnames: self.mappings_tir_sfr_cells.add_aux(temperature_name, self.cell_temperatures, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.mappings_tir_sfr_cells.add_aux(mean_age_name, self.cell_mean_ages, self.log_age_unit, as_column=True)
                #if funev_name not in colnames: self.mappings_tir_sfr_cells.add_aux(funev_name, self.cell_funev_values, as_column=True)
                #self.mappings_tir_sfr_cells.save() # save
                fs.remove_file(self.mappings_tir_sfr_cells_path)
                return False
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_sfr_cells_aux(self):
        return {temperature_name: self.cell_temperatures,
                mean_age_name: self.cell_mean_ages,
                funev_name: self.cell_funev_values}

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_sfr_cells_aux_units(self):
        return {temperature_name: self.temperature_unit,
                mean_age_name: self.log_age_unit}

    # -----------------------------------------------------------------

    def create_sfr_sfr_cell_scatter(self, sfr0_data, sfr1_data, sfr0_name, sfr1_name, sfr0_description=None, sfr1_description=None):

        """
        This function ...
        :param sfr0_data:
        :param sfr1_data:
        :param sfr0_name:
        :param sfr1_name:
        :param sfr0_description:
        :param sfr1_description:
        :return:
        """

        # Inform the user
        log.debug("Creating the " + sfr0_name + " to " + sfr1_name + " cell scatter ...")

        # Get sfr densities 0
        sfr_density0 = sfr0_data.values / self.cell_volumes
        sfr_density0_unit = sfr0_data.unit / self.volume_unit

        # Get sfr densities 1
        sfr_density1 = sfr1_data.values / self.cell_volumes
        sfr_density1_unit = sfr1_data.unit / self.volume_unit

        # Create and return
        return create_cell_scatter(sfr0_name, sfr1_name, sfr_density0, sfr_density1, sfr0_description, sfr1_description,
                                   x_unit=sfr_density0_unit, y_unit=sfr_density1_unit, aux=self.sfr_sfr_cells_aux,
                                   aux_units=self.sfr_sfr_cells_aux_units, is_arrays=True, aux_is_arrays=True)

    # -----------------------------------------------------------------

    def create_sfr_sfr_pixel_scatter(self, sfr0_frame, sfr1_frame, sfr0_name, sfr1_name, sfr0_description=None, sfr1_description=None):

        """
        This function ...
        :param sfr0_frame:
        :param sfr1_frame:
        :param sfr0_name:
        :param sfr1_name:
        :param sfr0_description:
        :param sfr1_description:
        :return:
        """

        # Inform the user
        log.info("Creating the " + sfr0_name + " to " + sfr1_name + " pixel scatter ...")

        # Get frames per unit of pixelarea
        sfr0_frame = sfr0_frame / sfr0_frame.pixelarea
        sfr1_frame = sfr1_frame / sfr1_frame.pixelarea

        # Create and return
        return create_pixel_scatter(sfr0_name, sfr1_name, sfr0_frame, sfr1_frame, sfr0_description, sfr1_description, same_units=True, convolve=True)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_tir_sfr_cells_path", True, write=False)
    def mappings_tir_sfr_cells(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_sfr_sfr_cell_scatter(self.cell_sfr_mappings, self.cell_sfr_tir, self.mappings_sfr_density_name, self.tir_sfr_density_name, self.mappings_sfr_density_description, self.tir_sfr_density_description)

    # -----------------------------------------------------------------

    @property
    def mappings_tir_sfr_pixels_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_tir_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_tir_sfr_pixels(self):
        return fs.is_file(self.mappings_tir_sfr_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_tir_sfr_pixels_path", True, write=False)
    def mappings_tir_sfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_sfr_sfr_pixel_scatter(self.pixel_sfr_mappings, self.pixel_sfr_tir, self.mappings_sfr_density_name, self.tir_sfr_density_name, self.mappings_sfr_density_description, self.tir_sfr_density_description)

    # -----------------------------------------------------------------
    # MAPPINGS K&E SFR - TIR SFR
    # -----------------------------------------------------------------

    @property
    def mappings_ke_sfr_name(self):
        return "SFR_MAPPINGS_KE"

    # -----------------------------------------------------------------

    @property
    def mappings_ke_sfr_description(self):
        return "Star formation rate (MAPPINGS + K&E)"

    # -----------------------------------------------------------------

    @property
    def mappings_ke_sfr_density_name(self):
        return "vSFR_MAPPINGS_KE"

    # -----------------------------------------------------------------

    @property
    def mappings_ke_sfr_density_description(self):
        return "Star formation rate density (MAPPINGS + K&E)"

    # -----------------------------------------------------------------

    @property
    def mappings_ke_tir_sfr_cells_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_ke_tir_cells.dat")

    # -----------------------------------------------------------------

    @property
    def mappings_ke_tir_sfr_cells_colnames(self):
        return fs.get_column_names(self.mappings_ke_tir_sfr_cells_path)

    # -----------------------------------------------------------------

    @property
    def mappings_ke_tir_sfr_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.mappings_ke_tir_sfr_cells_colnames, [self.mappings_ke_sfr_density_name, self.tir_sfr_density_name])

    # -----------------------------------------------------------------

    @property
    def mappings_ke_tir_sfr_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.mappings_ke_tir_sfr_cells_aux_colnames, sfr_sfr_cells_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_tir_sfr_cells(self):
        if fs.is_file(self.mappings_ke_tir_sfr_cells_path):
            if not self.mappings_ke_tir_sfr_cells_has_all_aux_columns:
                #colnames = self.mappings_ke_tir_sfr_cells_aux_colnames
                #if temperature_name not in colnames: self.mappings_ke_tir_sfr_cells.add_aux(temperature_name, self.cell_temperatures, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.mappings_ke_tir_sfr_cells.add_aux(mean_age_name, self.cell_mean_ages, self.log_age_unit, as_column=True)
                #if funev_name not in colnames: self.mappings_ke_tir_sfr_cells.add_aux(funev_name, self.cell_funev_values, as_column=True)
                #self.mappings_tir_sfr_cells.save() # save
                fs.remove_file(self.mappings_ke_tir_sfr_cells_path)
                return False
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_tir_sfr_cells_path", True, write=False)
    def mappings_ke_tir_sfr_cells(self):

        """
        This function ...
        :return:
        """

        # Create the scatter
        return self.create_sfr_sfr_cell_scatter(self.cell_sfr_mappings_ke, self.cell_sfr_tir, self.mappings_ke_sfr_density_name, self.tir_sfr_density_name, self.mappings_ke_sfr_density_description, self.tir_sfr_density_description)

    # -----------------------------------------------------------------

    @property
    def mappings_ke_tir_sfr_pixels_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_ke_tir_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_tir_sfr_pixels(self):
        return fs.is_file(self.mappings_ke_tir_sfr_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_tir_sfr_pixels_path", True, write=False)
    def mappings_ke_tir_sfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_sfr_sfr_pixel_scatter(self.pixel_sfr_mappings_ke, self.pixel_sfr_tir, self.mappings_ke_sfr_density_name, self.tir_sfr_density_name, self.mappings_ke_sfr_density_description, self.tir_sfr_density_description)

    # -----------------------------------------------------------------
    # MAPPINGS SFR - 24 um SFR
    # -----------------------------------------------------------------

    @property
    def mips24_sfr_name(self):
        return "SFR_24um"

    # -----------------------------------------------------------------

    @property
    def mips24_sfr_description(self):
        return "Star formation rate (MIPS 24 micron)"

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_sfr_densities_24um(self):
        return self.cell_sfr_24um_values / self.cell_volumes

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_density_24um_unit(self):
        return self.sfr_24um_unit / self.volume_unit

    # -----------------------------------------------------------------

    @property
    def mips24_sfr_density_name(self):
        return "vSFR_24um"

    # -----------------------------------------------------------------

    @property
    def mips24_sfr_density_description(self):
        return "Star formation rate density (24 micron)"

    # -----------------------------------------------------------------

    @property
    def mappings_24um_sfr_cells_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_24um_cells.dat")

    # -----------------------------------------------------------------

    @property
    def mappings_24um_sfr_cells_colnames(self):
        return fs.get_column_names(self.mappings_24um_sfr_cells_path)

    # -----------------------------------------------------------------

    @property
    def mappings_24um_sfr_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.mappings_24um_sfr_cells_colnames, [self.mappings_sfr_density_name, self.mips24_sfr_density_name])

    # -----------------------------------------------------------------

    @property
    def mappings_24um_sfr_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.mappings_24um_sfr_cells_aux_colnames, sfr_sfr_cells_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_mappings_24um_sfr_cells(self):
        if fs.is_file(self.mappings_24um_sfr_cells_path):
            if not self.mappings_24um_sfr_cells_has_all_aux_columns:
                #colnames = self.mappings_24um_sfr_cells_aux_colnames
                #if temperature_name not in colnames: self.mappings_24um_sfr_cells.add_aux(temperature_name, self.cell_temperatures, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.mappings_24um_sfr_cells.add_aux(mean_age_name, self.cell_mean_ages, self.log_age_unit, as_column=True)
                #if funev_name not in colnames: self.mappings_24um_sfr_cells.add_aux(funev_name, self.cell_funev_values, as_column=True)
                #self.mappings_24um_sfr_cells.save() # save
                fs.remove_file(self.mappings_24um_sfr_cells_path)
                return False
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_24um_sfr_cells_path", True, write=False)
    def mappings_24um_sfr_cells(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_sfr_sfr_cell_scatter(self.cell_sfr_mappings, self.cell_sfr_24um, self.mappings_sfr_density_name, self.mips24_sfr_density_name, self.mappings_sfr_density_description, self.mips24_sfr_density_description)

    # -----------------------------------------------------------------

    @property
    def mappings_24um_sfr_pixels_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_24um_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_24um_sfr_pixels(self):
        return fs.is_file(self.mappings_24um_sfr_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_24um_sfr_pixels_path", True, write=False)
    def mappings_24um_sfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_sfr_sfr_pixel_scatter(self.pixel_sfr_mappings, self.pixel_sfr_24um, self.mappings_sfr_density_name, self.mips24_sfr_density_name, self.mappings_sfr_density_description, self.mips24_sfr_density_description)

    # -----------------------------------------------------------------
    # MAPPINGS K&E SFR - 24um SFR
    # -----------------------------------------------------------------

    @property
    def mappings_ke_24um_sfr_cells_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_ke_24um_cells.dat")

    # -----------------------------------------------------------------

    @property
    def mappings_ke_24um_sfr_cells_colnames(self):
        return fs.get_column_names(self.mappings_ke_24um_sfr_cells_path)

    # -----------------------------------------------------------------

    @property
    def mappings_ke_24um_sfr_cells_aux_colnames(self):
        return sequences.elements_not_in_other(self.mappings_ke_24um_sfr_cells_colnames, [self.mappings_ke_sfr_density_name, self.mips24_sfr_density_name])

    # -----------------------------------------------------------------

    @property
    def mappings_ke_24um_sfr_cells_has_all_aux_columns(self):
        return sequences.contains_all(self.mappings_ke_24um_sfr_cells_aux_colnames, sfr_sfr_cells_aux_colnames)

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_24um_sfr_cells(self):
        if fs.is_file(self.mappings_ke_24um_sfr_cells_path):
            if not self.mappings_ke_24um_sfr_cells_has_all_aux_columns:
                #colnames = self.mappings_ke_24um_sfr_cells_aux_colnames
                #if temperature_name not in colnames: self.mappings_ke_24um_sfr_cells.add_aux(temperature_name, self.cell_temperatures, self.temperature_unit, as_column=True)
                #if mean_age_name not in colnames: self.mappings_ke_24um_sfr_cells.add_aux(mean_age_name, self.cell_mean_ages, self.log_age_unit, as_column=True)
                #if funev_name not in colnames: self.mappings_ke_24um_sfr_cells.add_aux(funev_name, self.cell_funev_values, as_column=True)
                #self.mappings_ke_24um_sfr_cells.save() # save
                fs.remove_file(self.mappings_ke_24um_sfr_cells_path)
                return False
            return True
        else: return False

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_24um_sfr_cells_path", True, write=False)
    def mappings_ke_24um_sfr_cells(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_sfr_sfr_cell_scatter(self.cell_sfr_mappings_ke, self.cell_sfr_24um, self.mappings_ke_sfr_density_name, self.mips24_sfr_density_name, self.mappings_ke_sfr_density_description, self.mips24_sfr_density_description)

    # -----------------------------------------------------------------

    @property
    def mappings_ke_24um_sfr_pixels_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_ke_24um_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_24um_sfr_pixels(self):
        return fs.is_file(self.mappings_ke_24um_sfr_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_24um_sfr_pixels_path", True, write=False)
    def mappings_ke_24um_sfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_sfr_sfr_pixel_scatter(self.pixel_sfr_mappings_ke, self.pixel_sfr_24um, self.mappings_ke_sfr_density_name, self.mips24_sfr_density_name, self.mappings_ke_sfr_density_description, self.mips24_sfr_density_description)

    # -----------------------------------------------------------------
    # MAPPINGS sSFR - FUV-H sSFR
    # -----------------------------------------------------------------

    def create_ssfr_colour_cell_scatter(self, ssfr_data, colour_data, ssfr_name, colour_name, ssfr_description, colour_description):

        """
        This function ...
        :param ssfr_data:
        :param colour_data:
        :param ssfr_name:
        :param colour_name:
        :param ssfr_description:
        :param colour_description:
        :return:
        """

        # Inform the user
        log.debug("Creating the " + ssfr_name + " to " + colour_name + " cell scatter ...")

        # Create and return
        return create_cell_scatter(ssfr_name, colour_name, ssfr_data, colour_data, ssfr_description, colour_description,
                                   x_unit=self.ssfr_unit, y_unit=self.magnitude_unit)

    # -----------------------------------------------------------------

    def create_ssfr_colour_pixel_scatter(self, ssfr_frame, colour_frame, ssfr_name, colour_name, ssfr_description, colour_description):

        """
        This function ...
        :param ssfr_frame:
        :param colour_frame:
        :param ssfr_name:
        :param colour_name:
        :param ssfr_description:
        :param colour_description:
        :return:
        """

        # Inform the user
        log.info("Creating the " + ssfr_name + " to " + colour_name + " pixel scatter ...")

        # Create and return
        return create_pixel_scatter(ssfr_name, colour_name, ssfr_frame, colour_frame, ssfr_description, colour_description, x_unit=self.ssfr_unit, y_unit=self.magnitude_unit, same_units=False, convolve=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_ssfr_path(self):
        return fs.create_directory_in(self.correlations_path, "sSFR-sSFR")

    # -----------------------------------------------------------------

    @property
    def mappings_fuv_h_ssfr_cells_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_fuv_h_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_fuv_h_ssfr_cells(self):
        return fs.is_file(self.mappings_fuv_h_ssfr_cells_path)

    # -----------------------------------------------------------------

    @property
    def mappings_ssfr_name(self):
        return "sSFR_MAPPINGS"

    # -----------------------------------------------------------------

    @property
    def mappings_ssfr_description(self):
        return "Specific star formation rate (MAPPINGS)"

    # -----------------------------------------------------------------

    @property
    def fuv_h_ssfr_name(self):
        return "FUV-H"

    # -----------------------------------------------------------------

    @property
    def fuv_h_ssfr_description(self):
        return "FUV-H colour as specific star formation rate proxy"

    # -----------------------------------------------------------------

    @lazyproperty
    def magnitude_unit(self):
        return u("mag")

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_fuv_h_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_fuv_h.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_fuv_h(self):
        return fs.is_file(self.cell_ssfr_fuv_h_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_fuv_h(self):
        return Data3D.from_file(self.cell_ssfr_fuv_h_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_fuv_h_values(self):
        return self.cell_ssfr_fuv_h.values

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_fuv_h_ssfr_cells_path", True, write=False)
    def mappings_fuv_h_ssfr_cells(self):

        """
        Thisfunction ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_cell_scatter(self.cell_ssfr_mappings, self.cell_ssfr_fuv_h, self.mappings_ssfr_name, self.fuv_h_ssfr_name, self.mappings_ssfr_description, self.fuv_h_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_fuv_h_ssfr_pixels_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_fuv_h_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_fuv_h_ssfr_pixels(self):
        return fs.is_file(self.mappings_fuv_h_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_fuv_h_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_fuv_h_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_fuv_h(self):
        return fs.is_file(self.pixel_ssfr_fuv_h_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_fuv_h(self):
        return Frame.from_file(self.pixel_ssfr_fuv_h_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_fuv_h_values(self):
        return self.pixel_ssfr_fuv_h.values

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_fuv_h_ssfr_pixels_path", True, write=False)
    def mappings_fuv_h_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_pixel_scatter(self.pixel_ssfr_mappings, self.pixel_ssfr_fuv_h, self.mappings_ssfr_name, self.fuv_h_ssfr_name, self.mappings_ssfr_description, self.fuv_h_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_fuv_r_ssfr_cells_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_fuv_r_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_fuv_r_ssfr_cells(self):
        return fs.is_file(self.mappings_fuv_r_ssfr_cells_path)

    # -----------------------------------------------------------------

    @property
    def fuv_r_ssfr_name(self):
        return "FUV-r"

    # -----------------------------------------------------------------

    @property
    def fuv_r_ssfr_description(self):
        return "FUV-r colour as specific star formation rate proxy"

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_fuv_r_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_fuv_r.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_fuv_r(self):
        return fs.is_file(self.cell_ssfr_fuv_r_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_fuv_r(self):
        return Data3D.from_file(self.cell_ssfr_fuv_r_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_fuv_r_values(self):
        return self.cell_ssfr_fuv_r.values

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_fuv_r_ssfr_cells_path", True, write=False)
    def mappings_fuv_r_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_cell_scatter(self.cell_ssfr_mappings, self.cell_ssfr_fuv_r, self.mappings_ssfr_name, self.fuv_r_ssfr_name, self.mappings_ssfr_description, self.fuv_r_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_fuv_r_ssfr_pixels_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_fuv_r_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_fuv_r_ssfr_pixels(self):
        return fs.is_file(self.mappings_fuv_r_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_fuv_r_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_fuv_r_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_fuv_r(self):
        return fs.is_file(self.pixel_ssfr_fuv_r_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_fuv_r(self):
        return Frame.from_file(self.pixel_ssfr_fuv_r_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_fuv_r_values(self):
        return self.pixel_ssfr_fuv_r.values

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_fuv_r_ssfr_pixels_path", True, write=False)
    def mappings_fuv_r_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_pixel_scatter(self.pixel_ssfr_mappings, self.pixel_ssfr_fuv_r, self.mappings_ssfr_name, self.fuv_r_ssfr_name, self.mappings_ssfr_description, self.fuv_r_ssfr_description)

    # -----------------------------------------------------------------
    # MAPPINGS K&E - FUV-H/R
    # -----------------------------------------------------------------

    @property
    def mappings_ke_ssfr_name(self):
        return "sSFR_MAPPINGS_KE"

    # -----------------------------------------------------------------

    @property
    def mappings_ke_ssfr_description(self):
        return "Specific star formation rate (MAPPINGS + K&E)"

    # -----------------------------------------------------------------

    @property
    def mappings_ke_fuv_h_ssfr_cells_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_ke_fuv_h_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_fuv_h_ssfr_cells(self):
        return fs.is_file(self.mappings_ke_fuv_h_ssfr_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_fuv_h_ssfr_cells_path", True, write=False)
    def mappings_ke_fuv_h_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_cell_scatter(self.cell_ssfr_mappings_ke, self.cell_ssfr_fuv_h, self.mappings_ke_ssfr_name, self.fuv_h_ssfr_name, self.mappings_ke_ssfr_description, self.fuv_h_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_ke_fuv_h_ssfr_pixels_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_ke_fuv_h_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_fuv_h_ssfr_pixels(self):
        return fs.is_file(self.mappings_ke_fuv_h_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_fuv_h_ssfr_pixels_path", True, write=False)
    def mappings_ke_fuv_h_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_pixel_scatter(self.pixel_ssfr_mappings_ke, self.pixel_ssfr_fuv_h, self.mappings_ke_ssfr_name, self.fuv_h_ssfr_name, self.mappings_ke_ssfr_description, self.fuv_h_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_ke_fuv_r_ssfr_cells_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_ke_fuv_r_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_fuv_r_ssfr_cells(self):
        return fs.is_file(self.mappings_ke_fuv_r_ssfr_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_fuv_r_ssfr_cells_path", True, write=False)
    def mappings_ke_fuv_r_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_cell_scatter(self.cell_ssfr_mappings_ke, self.cell_ssfr_fuv_r, self.mappings_ke_ssfr_name, self.fuv_r_ssfr_name, self.mappings_ke_ssfr_description, self.fuv_r_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_ke_fuv_r_ssfr_pixels_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_ke_fuv_r_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_fuv_r_ssfr_pixels(self):
        return fs.is_file(self.mappings_ke_fuv_r_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_fuv_r_ssfr_pixels_path", True, write=False)
    def mappings_ke_fuv_r_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_pixel_scatter(self.pixel_ssfr_mappings_ke, self.pixel_ssfr_fuv_r, self.mappings_ke_ssfr_name, self.fuv_r_ssfr_name, self.mappings_ke_ssfr_description, self.fuv_r_ssfr_description)

    # -----------------------------------------------------------------
    # MAPPINGS - NUV-H/R
    # -----------------------------------------------------------------

    @property
    def mappings_nuv_h_ssfr_cells_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_nuv_h_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_nuv_h_ssfr_cells(self):
        return fs.is_file(self.mappings_nuv_h_ssfr_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_nuv_h_ssfr_cells_path", True, write=False)
    def mappings_nuv_h_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_cell_scatter(self.cell_ssfr_mappings, self.cell_ssfr_nuv_h, self.mappings_ssfr_name, self.nuv_h_ssfr_name, self.mappings_ssfr_description, self.nuv_h_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_nuv_h_ssfr_pixels_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_nuv_h_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_nuv_h_ssfr_pixels(self):
        return fs.is_file(self.mappings_nuv_h_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_nuv_h_ssfr_pixels_path", True, write=False)
    def mappings_nuv_h_ssfr_pixels(self):

        """
        Thisf unction ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_pixel_scatter(self.pixel_ssfr_mappings, self.pixel_ssfr_nuv_h, self.mappings_ssfr_name, self.nuv_h_ssfr_name, self.mappings_ssfr_description, self.nuv_h_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_nuv_r_ssfr_cells_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_nuv_r_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_nuv_r_ssfr_cells(self):
        return fs.is_file(self.mappings_nuv_r_ssfr_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_nuv_r_ssfr_cells_path", True, write=False)
    def mappings_nuv_r_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_cell_scatter(self.cell_ssfr_mappings, self.cell_ssfr_nuv_r, self.mappings_ssfr_name, self.nuv_r_ssfr_name, self.mappings_ssfr_description, self.nuv_r_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_nuv_r_ssfr_pixels_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_nuv_r_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_nuv_r_ssfr_pixels(self):
        return fs.is_file(self.mappings_nuv_r_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_nuv_r_ssfr_pixels_path", True, write=False)
    def mappings_nuv_r_ssfr_pixels(self):

        """
        Thisf unction ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_pixel_scatter(self.pixel_ssfr_mappings, self.pixel_ssfr_nuv_r, self.mappings_ssfr_name, self.nuv_r_ssfr_name, self.mappings_ssfr_description, self.nuv_r_ssfr_description)

    # -----------------------------------------------------------------
    # MAPPINGS KE - NUV-H/R
    # -----------------------------------------------------------------

    @property
    def mappings_ke_nuv_h_ssfr_cells_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_ke_nuv_h_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_nuv_h_ssfr_cells(self):
        return fs.is_file(self.mappings_ke_nuv_h_ssfr_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_nuv_h_ssfr_cells_path", True, write=False)
    def mappings_ke_nuv_h_ssfr_cells(self):
        
        """
        This fnuction ...
        :return: 
        """

        # Create scatter
        return self.create_ssfr_colour_cell_scatter(self.cell_ssfr_mappings_ke, self.cell_ssfr_nuv_h, self.mappings_ke_ssfr_name, self.nuv_h_ssfr_name, self.mappings_ke_ssfr_description, self.nuv_h_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_ke_nuv_h_ssfr_pixels_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_ke_nuv_h_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_nuv_h_ssfr_pixels(self):
        return fs.is_file(self.mappings_ke_nuv_h_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_nuv_h_ssfr_pixels_path", True, write=False)
    def mappings_ke_nuv_h_ssfr_pixels(self):
        
        """
        This function ...
        :return: 
        """

        # Create scatter
        return self.create_ssfr_colour_pixel_scatter(self.pixel_ssfr_mappings_ke, self.pixel_ssfr_nuv_h, self.mappings_ke_ssfr_name, self.nuv_h_ssfr_name, self.mappings_ke_ssfr_description, self.nuv_h_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_ke_nuv_r_ssfr_cells_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_ke_nuv_r_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_nuv_r_ssfr_cells(self):
        return fs.is_file(self.mappings_ke_nuv_r_ssfr_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_nuv_r_ssfr_cells_path", True, write=False)
    def mappings_ke_nuv_r_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_cell_scatter(self.cell_ssfr_mappings_ke, self.cell_ssfr_nuv_r, self.mappings_ke_ssfr_name, self.nuv_r_ssfr_name, self.mappings_ke_ssfr_description, self.nuv_r_ssfr_description)

    # -----------------------------------------------------------------

    @property
    def mappings_ke_nuv_r_ssfr_pixels_path(self):
        return fs.join(self.ssfr_ssfr_path, "mappings_ke_nuv_r_pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_mappings_ke_nuv_r_ssfr_pixels(self):
        return fs.is_file(self.mappings_ke_nuv_r_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "mappings_ke_nuv_r_ssfr_pixels_path", True, write=False)
    def mappings_ke_nuv_r_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Create scatter
        return self.create_ssfr_colour_pixel_scatter(self.pixel_ssfr_mappings_ke, self.pixel_ssfr_nuv_r, self.mappings_ke_ssfr_name, self.nuv_r_ssfr_name, self.mappings_ke_ssfr_description, self.nuv_r_ssfr_description)

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

        # SFR SFR scatter data
        self.write_sfr_sfr()

        # sSFR sSFR scatter data
        self.write_ssfr_ssfr()

        # Mean age Funev scatter data
        self.write_mean_age_funev()
        
        # Mean age sSFR scatter data
        self.write_mean_age_ssfr()

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

        # FUV-H
        self.write_fuv_h_funev()

        # FUV-R
        self.write_fuv_r_funev()

        # NUV-H
        self.write_nuv_h_funev()

        # NUV-R
        self.write_nuv_r_funev()

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
        log.info("Writing the sSFR (Salim) to Funev dust cell scatter data ...")

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
        if self.do_write_ssfr_ke_funev_cells: self.write_ssfr_ke_funev_cells()

        # Pixels
        if self.do_write_ssfr_ke_funev_pixels: self.write_ssfr_ke_funev_pixels()

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

    @property
    def do_write_fuv_h_funev_cells(self):
        return not self.has_fuv_h_funev_cells

    # -----------------------------------------------------------------

    @property
    def do_write_fuv_h_funev_pixels(self):
        return not self.has_fuv_h_funev_pixels

    # -----------------------------------------------------------------

    def write_fuv_h_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_fuv_h_funev_cells: self.write_fuv_h_funev_cells()

        # Pixels
        if self.do_write_fuv_h_funev_pixels: self.write_fuv_h_funev_pixels()

    # -----------------------------------------------------------------

    def write_fuv_h_funev_cells(self):

        """
        This function ...
        :return:
        """

        self.fuv_h_funev_cells.saveto(self.fuv_h_funev_cells_path)

    # -----------------------------------------------------------------

    def write_fuv_h_funev_pixels(self):

        """
        This function ...
        :return:
        """

        self.fuv_h_funev_pixels.saveto(self.fuv_h_funev_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_fuv_r_funev_cells(self):
        return not self.has_fuv_r_funev_cells

    # -----------------------------------------------------------------

    @property
    def do_write_fuv_r_funev_pixels(self):
        return not self.has_fuv_r_funev_pixels

    # -----------------------------------------------------------------

    def write_fuv_r_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_fuv_r_funev_cells: self.write_fuv_r_funev_cells()

        # Pixels
        if self.do_write_fuv_r_funev_pixels: self.write_fuv_r_funev_pixels()

    # -----------------------------------------------------------------

    def write_fuv_r_funev_cells(self):

        """
        This function ...
        :return:
        """

        self.fuv_r_funev_cells.saveto(self.fuv_r_funev_cells_path)

    # -----------------------------------------------------------------

    def write_fuv_r_funev_pixels(self):

        """
        This function ...
        :return:
        """

        self.fuv_r_funev_pixels.saveto(self.fuv_r_funev_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_nuv_h_funev_cells(self):
        return not self.has_nuv_h_funev_cells

    # -----------------------------------------------------------------

    @property
    def do_write_nuv_h_funev_pixels(self):
        return not self.has_nuv_h_funev_pixels

    # -----------------------------------------------------------------

    def write_nuv_h_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_nuv_h_funev_cells: self.write_nuv_h_funev_cells()

        # Pixels
        if self.do_write_nuv_h_funev_pixels: self.write_nuv_h_funev_pixels()

    # -----------------------------------------------------------------

    def write_nuv_h_funev_cells(self):

        """
        This function ...
        :return:
        """

        self.nuv_h_funev_cells.saveto(self.nuv_h_funev_cells_path)

    # -----------------------------------------------------------------

    def write_nuv_h_funev_pixels(self):

        """
        This function ...
        :return:
        """

        self.nuv_h_funev_pixels.saveto(self.nuv_h_funev_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_nuv_r_funev_cells(self):
        return not self.has_nuv_r_funev_cells

    # -----------------------------------------------------------------

    @property
    def do_write_nuv_r_funev_pixels(self):
        return not self.has_nuv_r_funev_pixels

    # -----------------------------------------------------------------

    def write_nuv_r_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_nuv_r_funev_cells: self.write_nuv_r_funev_cells()

        # Pixels
        if self.do_write_nuv_r_funev_pixels: self.write_nuv_r_funev_pixels()

    # -----------------------------------------------------------------

    def write_nuv_r_funev_cells(self):

        """
        This function ...
        :return:
        """

        self.nuv_r_funev_cells.saveto(self.nuv_r_funev_cells_path)

    # -----------------------------------------------------------------

    def write_nuv_r_funev_pixels(self):

        """
        This function ...
        :return:
        """

        self.nuv_r_funev_pixels.saveto(self.nuv_r_funev_pixels_path)

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

        # Inform the user
        log.info("Writing the dust temperature to Funev scatter data ...")

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

    def write_sfr_sfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SFR-SFR correlation data ...")

        # MAPPINGS to TIR
        self.write_mappings_tir_sfr()

        # MAPPINGS to 24um
        self.write_mappings_24um_sfr()

        # MAPPINGS + K&E to TIR
        self.write_mappings_ke_tir_sfr()

        # MAPPINGS + K&E to 24um
        self.write_mappings_ke_24um_sfr()

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_tir_sfr_cells(self):
        return not self.has_mappings_tir_sfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_tir_sfr_pixels(self):
        return not self.has_mappings_tir_sfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_tir_sfr(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_mappings_tir_sfr_cells: self.write_mappings_tir_sfr_cells()

        # Pixels
        if self.do_write_mappings_tir_sfr_pixels: self.write_mappings_tir_sfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_tir_sfr_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the MAPPINGS SFR to TIR SFR dust cell scatter data ...")

        # Write
        self.mappings_tir_sfr_cells.saveto(self.mappings_tir_sfr_cells_path)

    # -----------------------------------------------------------------

    def write_mappings_tir_sfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the MAPPINGS SFR to TIR SFR pixel scatter data ...")

        # Write
        self.mappings_tir_sfr_pixels.saveto(self.mappings_tir_sfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_24um_sfr_cells(self):
        return not self.has_mappings_24um_sfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_24um_sfr_pixels(self):
        return not self.has_mappings_24um_sfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_24um_sfr(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_mappings_24um_sfr_cells: self.write_mappings_24um_sfr_cells()

        # Pixels
        if self.do_write_mappings_24um_sfr_pixels: self.write_mappings_24um_sfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_24um_sfr_cells(self):

        """
        THins function ...
        :return:
        """

        # Inform the user
        log.info("Writing the MAPPINGS SFR to 24 micron SFR dust cell scatter data ...")

        # Write
        self.mappings_24um_sfr_cells.saveto(self.mappings_24um_sfr_cells_path)

    # -----------------------------------------------------------------

    def write_mappings_24um_sfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the MAPPINGS SFR to 24 micron SFR dust pixel scatter data ...")

        # Write
        self.mappings_24um_sfr_pixels.saveto(self.mappings_24um_sfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_tir_sfr_cells(self):
        return not self.has_mappings_ke_tir_sfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_tir_sfr_pixels(self):
        return not self.has_mappings_ke_tir_sfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_ke_tir_sfr(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_mappings_ke_tir_sfr_cells: self.write_mappings_ke_tir_sfr_cells()

        # Pixels
        if self.do_write_mappings_ke_tir_sfr_pixels: self.write_mappings_ke_tir_sfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_ke_tir_sfr_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the MAPPINGS K&E SFR to TIR SFR dust cell scatter data ...")

        # Write
        self.mappings_ke_tir_sfr_cells.saveto(self.mappings_ke_tir_sfr_cells_path)

    # -----------------------------------------------------------------

    def write_mappings_ke_tir_sfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the MAPPINGS K&E SFR to TIR SFR pixel scatter data ...")

        # Write
        self.mappings_ke_tir_sfr_pixels.saveto(self.mappings_ke_tir_sfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_24um_sfr_cells(self):
        return not self.has_mappings_ke_24um_sfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_24um_sfr_pixels(self):
        return not self.has_mappings_ke_24um_sfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_ke_24um_sfr(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_mappings_ke_24um_sfr_cells: self.write_mappings_ke_24um_sfr_cells()

        # Pixels
        if self.do_write_mappings_ke_24um_sfr_pixels: self.write_mappings_ke_24um_sfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_ke_24um_sfr_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the MAPPINGS + K&E SFR to 24 micron SFR dust cell scatter data ...")

        # Write
        self.mappings_ke_24um_sfr_cells.saveto(self.mappings_ke_24um_sfr_cells_path)

    # -----------------------------------------------------------------

    def write_mappings_ke_24um_sfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the MAPPINGS + K&E SFR to 24 micron SFR dust pixel scatter data ...")

        # Write
        self.mappings_ke_24um_sfr_pixels.saveto(self.mappings_ke_24um_sfr_pixels_path)

    # -----------------------------------------------------------------

    def write_ssfr_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR-sSFR correlation data ...")

        # MAPPINGS to FUV-H
        self.write_mappings_fuv_h_ssfr()

        # MAPPINGS to FUV-R
        self.write_mappings_fuv_r_ssfr()

        # MAPPINGS to NUV-H
        self.write_mappings_nuv_h_ssfr()

        # MAPPINGS to NUV-R
        self.write_mappings_nuv_r_ssfr()

        # MAPPINGS KE to FUV-H
        self.write_mappings_ke_fuv_h_ssfr()

        # MAPPINGS KE to FUV-R
        self.write_mappings_ke_fuv_r_ssfr()

        # MAPPINGS KE to NUV-H
        self.write_mappings_ke_nuv_h_ssfr()

        # MAPPINGS KE to NUV-R
        self.write_mappings_ke_nuv_r_ssfr()

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_fuv_h_ssfr_cells(self):
        return not self.has_mappings_fuv_h_ssfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_fuv_h_ssfr_pixels(self):
        return not self.has_mappings_fuv_h_ssfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_fuv_h_ssfr(self):

        """
        Thisf unction ...
        :return:
        """

        # Cells
        if self.do_write_mappings_fuv_h_ssfr_cells: self.write_mappings_fuv_h_ssfr_cells()

        # Pixels
        if self.do_write_mappings_fuv_h_ssfr_pixels: self.write_mappings_fuv_h_ssfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_fuv_h_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_fuv_h_ssfr_cells.saveto(self.mappings_fuv_h_ssfr_cells_path)

    # -----------------------------------------------------------------

    def write_mappings_fuv_h_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_fuv_h_ssfr_pixels.saveto(self.mappings_fuv_h_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_fuv_r_ssfr_cells(self):
        return not self.has_mappings_fuv_r_ssfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_fuv_r_ssfr_pixels(self):
        return not self.has_mappings_fuv_r_ssfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_fuv_r_ssfr(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_mappings_fuv_r_ssfr_cells: self.write_mappings_fuv_r_ssfr_cells()

        # Pixels
        if self.do_write_mappings_fuv_r_ssfr_pixels: self.write_mappings_fuv_r_ssfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_fuv_r_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_fuv_r_ssfr_cells.saveto(self.mappings_fuv_r_ssfr_cells_path)

    # -----------------------------------------------------------------

    def write_mappings_fuv_r_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_fuv_r_ssfr_pixels.saveto(self.mappings_fuv_r_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_nuv_h_ssfr_cells(self):
        return not self.has_mappings_nuv_h_ssfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_nuv_h_ssfr_pixels(self):
        return not self.has_mappings_nuv_h_ssfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_nuv_h_ssfr(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_mappings_nuv_h_ssfr_cells: self.write_mappings_nuv_h_ssfr_cells()

        # Pixels
        if self.do_write_mappings_nuv_h_ssfr_pixels: self.write_mappings_nuv_h_ssfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_nuv_h_ssfr_cells(self):
        
        """
        This function ...
        :return: 
        """

        # Write
        self.mappings_nuv_h_ssfr_cells.saveto(self.mappings_nuv_h_ssfr_cells_path)
        
    # -----------------------------------------------------------------
        
    def write_mappings_nuv_h_ssfr_pixels(self):
        
        """
        Ths function ...
        :return: 
        """

        # Write
        self.mappings_nuv_h_ssfr_pixels.saveto(self.mappings_nuv_h_ssfr_pixels_path)
        
    # -----------------------------------------------------------------

    @property
    def do_write_mappings_nuv_r_ssfr_cells(self):
        return not self.has_mappings_nuv_r_ssfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_nuv_r_ssfr_pixels(self):
        return not self.has_mappings_nuv_r_ssfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_nuv_r_ssfr(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_mappings_nuv_r_ssfr_cells: self.write_mappings_nuv_r_ssfr_cells()

        # Pixels
        if self.do_write_mappings_nuv_r_ssfr_pixels: self.write_mappings_nuv_r_ssfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_nuv_r_ssfr_cells(self):

        """
        Thisf unction ...
        :return:
        """

        # Write
        self.mappings_nuv_r_ssfr_cells.saveto(self.mappings_nuv_r_ssfr_cells_path)

    # -----------------------------------------------------------------

    def write_mappings_nuv_r_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_nuv_r_ssfr_pixels.saveto(self.mappings_nuv_r_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_fuv_h_ssfr_cells(self):
        return not self.has_mappings_ke_fuv_h_ssfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_fuv_h_ssfr_pixels(self):
        return not self.has_mappings_ke_fuv_h_ssfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_ke_fuv_h_ssfr(self):

        """
        Thisf unction ...
        :return:
        """

        # Cells
        if self.do_write_mappings_ke_fuv_h_ssfr_cells: self.write_mappings_ke_fuv_h_ssfr_cells()

        # Pixels
        if self.do_write_mappings_ke_fuv_h_ssfr_pixels: self.write_mappings_ke_fuv_h_ssfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_ke_fuv_h_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_ke_fuv_h_ssfr_cells.saveto(self.mappings_ke_fuv_h_ssfr_cells_path)

    # -----------------------------------------------------------------

    def write_mappings_ke_fuv_h_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_ke_fuv_h_ssfr_pixels.saveto(self.mappings_ke_fuv_h_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_fuv_r_ssfr_cells(self):
        return not self.has_mappings_ke_fuv_r_ssfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_fuv_r_ssfr_pixels(self):
        return not self.has_mappings_ke_fuv_r_ssfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_ke_fuv_r_ssfr(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_mappings_ke_fuv_r_ssfr_cells: self.write_mappings_ke_fuv_r_ssfr_cells()

        # Pixels
        if self.do_write_mappings_ke_fuv_r_ssfr_pixels: self.write_mappings_ke_fuv_r_ssfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_ke_fuv_r_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_ke_fuv_r_ssfr_cells.saveto(self.mappings_ke_fuv_r_ssfr_cells_path)

    # -----------------------------------------------------------------

    def write_mappings_ke_fuv_r_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_ke_fuv_r_ssfr_pixels.saveto(self.mappings_ke_fuv_r_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_nuv_h_ssfr_cells(self):
        return not self.has_mappings_ke_nuv_h_ssfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_nuv_h_ssfr_pixels(self):
        return not self.has_mappings_ke_nuv_h_ssfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_ke_nuv_h_ssfr(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_mappings_ke_nuv_h_ssfr_cells: self.write_mappings_ke_nuv_h_ssfr_cells()

        # Pixels
        if self.do_write_mappings_ke_nuv_h_ssfr_pixels: self.write_mappings_ke_nuv_h_ssfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_ke_nuv_h_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_ke_nuv_h_ssfr_cells.saveto(self.mappings_ke_nuv_h_ssfr_cells_path)

    # -----------------------------------------------------------------

    def write_mappings_ke_nuv_h_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_ke_nuv_h_ssfr_pixels.saveto(self.mappings_ke_nuv_h_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_nuv_r_ssfr_cells(self):
        return not self.has_mappings_ke_nuv_r_ssfr_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mappings_ke_nuv_r_ssfr_pixels(self):
        return not self.has_mappings_ke_nuv_r_ssfr_pixels

    # -----------------------------------------------------------------

    def write_mappings_ke_nuv_r_ssfr(self):

        """
        This function ..
        :return:
        """

        # Cells
        if self.do_write_mappings_ke_nuv_r_ssfr_cells: self.write_mappings_ke_nuv_r_ssfr_cells()

        # Pixels
        if self.do_write_mappings_ke_nuv_r_ssfr_pixels: self.write_mappings_ke_nuv_r_ssfr_pixels()

    # -----------------------------------------------------------------

    def write_mappings_ke_nuv_r_ssfr_cells(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_ke_nuv_r_ssfr_cells.saveto(self.mappings_ke_nuv_r_ssfr_cells_path)

    # -----------------------------------------------------------------

    def write_mappings_ke_nuv_r_ssfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mappings_ke_nuv_r_ssfr_pixels.saveto(self.mappings_ke_nuv_r_ssfr_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_mean_age_cells(self):
        return self.has_mean_age_cells

    # -----------------------------------------------------------------

    @property
    def do_mean_age_pixels(self):
        return self.has_mean_age_pixels

    # -----------------------------------------------------------------

    @property
    def do_write_mean_age_funev_cells(self):
        return self.do_mean_age_cells and not self.has_mean_age_funev_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mean_age_funev_pixels(self):
        return self.do_mean_age_pixels and not self.has_mean_age_funev_pixels

    # -----------------------------------------------------------------

    def write_mean_age_funev(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the mean stellar age to Funev scatter data ...")

        # Cells
        if self.do_write_mean_age_funev_cells: self.write_mean_age_funev_cells()

        # Pixels
        if self.do_write_mean_age_funev_pixels: self.write_mean_age_funev_pixels()

    # -----------------------------------------------------------------

    def write_mean_age_funev_cells(self):

        """
        This functio n...
        :return:
        """

        # Write
        self.mean_age_funev_cells.saveto(self.mean_age_funev_cells_path)

    # -----------------------------------------------------------------

    def write_mean_age_funev_pixels(self):

        """
        Thisn function ...
        :return:
        """

        # Write
        self.mean_age_funev_pixels.saveto(self.mean_age_funev_pixels_path)

    # -----------------------------------------------------------------

    def write_mean_age_ssfr(self):
        
        """
        This function ...
        :return: 
        """
        
        # MAPPINGS
        self.write_mean_age_ssfr_mappings()

        # MAPPINGS + K&E
        self.write_mean_age_ssfr_mappings_ke()

    # -----------------------------------------------------------------

    @property
    def do_write_mean_age_ssfr_mappings_cells(self):
        return self.do_mean_age_cells and not self.has_mean_age_ssfr_mappings_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mean_age_ssfr_mappings_pixels(self):
        return self.do_mean_age_pixels and not self.has_mean_age_ssfr_mappings_pixels

    # -----------------------------------------------------------------

    def write_mean_age_ssfr_mappings(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_mean_age_ssfr_mappings_cells: self.write_mean_age_ssfr_mappings_cells()

        # Pixels
        if self.do_write_mean_age_ssfr_mappings_pixels: self.write_mean_age_ssfr_mappings_pixels()

    # -----------------------------------------------------------------

    def write_mean_age_ssfr_mappings_cells(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mean_age_ssfr_mappings_cells.saveto(self.mean_age_ssfr_mappings_cells_path)

    # -----------------------------------------------------------------

    def write_mean_age_ssfr_mappings_pixels(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mean_age_ssfr_mappings_pixels.saveto(self.mean_age_ssfr_mappings_pixels_path)

    # -----------------------------------------------------------------

    @property
    def do_write_mean_age_ssfr_mappings_ke_cells(self):
        return self.do_mean_age_cells and not self.has_mean_age_ssfr_mappings_ke_cells

    # -----------------------------------------------------------------

    @property
    def do_write_mean_age_ssfr_mappings_ke_pixels(self):
        return self.do_mean_age_pixels and not self.has_mean_age_ssfr_mappings_ke_pixels

    # -----------------------------------------------------------------

    def write_mean_age_ssfr_mappings_ke(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_write_mean_age_ssfr_mappings_ke_cells: self.write_mean_age_ssfr_mappings_ke_cells()

        # Pixels
        if self.do_write_mean_age_ssfr_mappings_ke_pixels: self.write_mean_age_ssfr_mappings_ke_pixels()

    # -----------------------------------------------------------------

    def write_mean_age_ssfr_mappings_ke_cells(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mean_age_ssfr_mappings_ke_cells.saveto(self.mean_age_ssfr_mappings_ke_cells_path)

    # -----------------------------------------------------------------

    def write_mean_age_ssfr_mappings_ke_pixels(self):

        """
        This function ...
        :return:
        """

        # Write
        self.mean_age_ssfr_mappings_ke_pixels.saveto(self.mean_age_ssfr_mappings_ke_pixels_path)

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

        # SFR SFR scatter data
        self.plot_sfr_sfr()

        # sSFR sSFR scatter data
        self.plot_ssfr_ssfr()

        # Plot mean age - Funev
        self.plot_mean_age_funev()

        # Mean age sSFR scatter data
        self.plot_mean_age_ssfr()

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
        else: plot_scatters_density(self.ssfr_salim_funev_cells_scatters, title=self.ssfr_salim_funev_cells_title, xlog=True,
                      path=self.ssfr_salim_funev_cells_plot_path, xlimits=self.ssfr_limits, ylimits=self.funev_limits)

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
        else: plot_scatters_density(self.ssfr_salim_funev_pixels_scatters, title=self.ssfr_salim_funev_pixels_title, xlog=True,
                      path=self.ssfr_salim_funev_pixels_plot_path, xlimits=self.ssfr_limits, ylimits=self.funev_limits)

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
        else: plot_scatters_density(self.ssfr_ke_funev_cells_scatters, title=self.ssfr_ke_funev_cells_title, xlog=True,
                          path=self.ssfr_ke_funev_cells_plot_path, xlimits=self.ssfr_limits,
                          ylimits=self.funev_limits)

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
        else: plot_scatters_density(self.ssfr_ke_funev_pixels_scatters, title=self.ssfr_ke_funev_pixels_title, xlog=True,
                          path=self.ssfr_ke_funev_pixels_plot_path, xlimits=self.ssfr_limits,
                          ylimits=self.funev_limits)

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
            plot_scatters_density(self.ssfr_mappings_funev_cells_scatters, title=self.ssfr_mappings_funev_cells_title,
                          xlog=True, path=self.ssfr_mappings_funev_cells_plot_path, xlimits=self.ssfr_limits,
                          ylimits=self.funev_limits)

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
            plot_scatters_density(self.ssfr_mappings_funev_pixels_scatters, title=self.ssfr_mappings_funev_pixels_title,
                          xlog=True, path=self.ssfr_mappings_funev_pixels_plot_path, xlimits=self.ssfr_limits,
                          ylimits=self.funev_limits)

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
        else: plot_scatters_density(self.ssfr_mappings_ke_funev_cells_scatters, title=self.ssfr_mappings_ke_funev_cells_title, xlog=True,
                              path=self.ssfr_mappings_ke_funev_cells_plot_path, xlimits=self.ssfr_limits,
                              ylimits=self.funev_limits)

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
        else: plot_scatters_density(self.ssfr_mappings_ke_funev_pixels_scatters, title=self.ssfr_mappings_ke_funev_pixels_title, xlog=True,
                              path=self.ssfr_mappings_ke_funev_pixels_plot_path, xlimits=self.ssfr_limits,
                              ylimits=self.funev_limits)

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

    @property
    def sfr_limits(self):
        return None

    # -----------------------------------------------------------------

    def plot_sfr_sfr(self):

        """
        This function ...
        :return:
        """

        # MAPPINGS to TIR
        self.plot_mappings_tir_sfr()

        # MAPPINGS to 24um
        self.plot_mappings_24um_sfr()

        # MAPPINGS KE to TIR
        self.plot_mappings_ke_tir_sfr()

        # MAPPINGS KE to 24um
        self.plot_mappings_ke_24um_sfr()

    # -----------------------------------------------------------------

    @property
    def do_plot_mappings_tir_sfr_cells(self):
        return not self.has_mappings_tir_sfr_cells_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_mappings_tir_sfr_pixels(self):
        return not self.has_mappings_tir_sfr_pixels_plot

    # -----------------------------------------------------------------

    def plot_mappings_tir_sfr(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_plot_mappings_tir_sfr_cells: self.plot_mappings_tir_sfr_cells()

        # Pixels
        if self.do_plot_mappings_tir_sfr_pixels: self.plot_mappings_tir_sfr_pixels()

    # -----------------------------------------------------------------

    @property
    def mappings_tir_sfr_cells_scatter_paths(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.mappings_tir_sfr_cells_path
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    @property
    def mappings_tir_sfr_cells_title(self):
        return "Correlation between MAPPINGS SFR and TIR SFR of dust cells"

    # -----------------------------------------------------------------

    @property
    def mappings_tir_sfr_cells_plot_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_tir_cells.pdf")

    # -----------------------------------------------------------------

    @property
    def has_mappings_tir_sfr_cells_plot(self):
        return fs.is_file(self.mappings_tir_sfr_cells_plot_path)

    # -----------------------------------------------------------------

    def plot_mappings_tir_sfr_cells(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat:

            plot_stilts(self.mappings_tir_sfr_cells_scatter_paths, self.mappings_sfr_name, self.tir_sfr_name,
                        self.mappings_sfr_description, self.tir_sfr_description,
                        title=self.mappings_tir_sfr_cells_title, path=self.mappings_tir_sfr_cells_plot_path,
                        ylimits=self.sfr_limits, xlog=True, xlimits=self.sfr_limits)

        # Plot using Matplotlib
        else: raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def mappings_tir_sfr_pixels_scatter_paths(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.mappings_tir_sfr_pixels_path
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    @property
    def mappings_tir_sfr_pixels_title(self):
        return "Correlation between MAPPINGS SFR and TIR SFR of model pixels"

    # -----------------------------------------------------------------

    @property
    def mappings_tir_sfr_pixels_plot_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_tir_pixels.pdf")

    # -----------------------------------------------------------------

    @property
    def has_mappings_tir_sfr_pixels_plot(self):
        return fs.is_file(self.mappings_tir_sfr_pixels_plot_path)

    # -----------------------------------------------------------------

    def plot_mappings_tir_sfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat:

            plot_stilts(self.mappings_tir_sfr_pixels_scatter_paths, self.mappings_sfr_name, self.tir_sfr_name,
                        self.mappings_sfr_description, self.tir_sfr_description,
                        title=self.mappings_tir_sfr_pixels_title, path=self.mappings_tir_sfr_pixels_plot_path,
                        ylimits=self.sfr_limits, xlog=True, xlimits=self.sfr_limits)

        # Plot using Matplotlib
        else: raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def do_plot_mappings_24um_sfr_cells(self):
        return not self.has_mappings_24um_sfr_cells_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_mappings_24um_sfr_pixels(self):
        return not self.has_mappings_24um_sfr_pixels_plot

    # -----------------------------------------------------------------

    def plot_mappings_24um_sfr(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_plot_mappings_24um_sfr_cells: self.plot_mappings_24um_sfr_cells()

        # Pixels
        if self.do_plot_mappings_24um_sfr_pixels: self.plot_mappings_24um_sfr_pixels()

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings_24um_sfr_cells_scatter_paths(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.mappings_24um_sfr_cells_path
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    @property
    def mappings_24um_sfr_cells_title(self):
        return "Correlation between MAPPINGS SFR and 24 micron SFR of dust cells"

    # -----------------------------------------------------------------

    @property
    def mappings_24um_sfr_cells_plot_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_24um_cells.pdf")

    # -----------------------------------------------------------------

    @property
    def has_mappings_24um_sfr_cells_plot(self):
        return fs.is_file(self.mappings_24um_sfr_cells_plot_path)

    # -----------------------------------------------------------------

    def plot_mappings_24um_sfr_cells(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat:

            plot_stilts(self.mappings_24um_sfr_cells_scatter_paths, self.mappings_sfr_name, self.mips24_sfr_name,
                        self.mappings_sfr_description, self.mips24_sfr_description,
                        title=self.mappings_24um_sfr_cells_title, path=self.mappings_24um_sfr_cells_plot_path,
                        ylimits=self.sfr_limits, xlog=True, xlimits=self.sfr_limits)

        # Plot using Matplotlib
        else: raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings_24um_sfr_pixels_scatter_paths(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.mappings_24um_sfr_pixels_path
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    @property
    def mappings_24um_sfr_pixels_title(self):
        return "Correlation between MAPPINGS SFR and 24 micron SFR of model pixels"

    # -----------------------------------------------------------------

    @property
    def mappings_24um_sfr_pixels_plot_path(self):
        return fs.join(self.sfr_sfr_path, "mappings_24um_pixels.pdf")

    # -----------------------------------------------------------------

    @property
    def has_mappings_24um_sfr_pixels_plot(self):
        return fs.is_file(self.mappings_24um_sfr_pixels_plot_path)

    # -----------------------------------------------------------------

    def plot_mappings_24um_sfr_pixels(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat:

            plot_stilts(self.mappings_24um_sfr_pixels_scatter_paths, self.mappings_sfr_name, self.mips24_sfr_name,
                        self.mappings_sfr_description, self.mips24_sfr_description,
                        title=self.mappings_24um_sfr_pixels_title, path=self.mappings_24um_sfr_pixels_plot_path,
                        ylimits=self.sfr_limits, xlog=True, xlimits=self.sfr_limits)

        # Plot using Matplotlib
        else: raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    def plot_mappings_ke_tir_sfr(self):

        """
        This function ..
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def plot_mappings_ke_24um_sfr(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def plot_ssfr_ssfr(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def do_plot_mean_age_funev_cells(self):
        return self.do_mean_age_cells and not self.has_mean_age_funev_cells_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_mean_age_funev_pixels(self):
        return self.do_mean_age_pixels and not self.has_mean_age_funev_pixels_plot

    # -----------------------------------------------------------------

    @property
    def mean_age_limits(self):
        return None

    # -----------------------------------------------------------------

    def plot_mean_age_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_plot_mean_age_funev_cells: self.plot_mean_age_funev_cells()

        # Pixels
        if self.do_plot_mean_age_funev_pixels: self.plot_mean_age_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def mean_age_funev_cells_scatters(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.mean_age_funev_cells
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    @property
    def mean_age_funev_cells_scatter_paths(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.mean_age_funev_cells_path
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    @property
    def mean_age_funev_cells_plot_path(self):
        return fs.join(self.mean_age_funev_path, "cells.pdf")

    # -----------------------------------------------------------------

    @property
    def has_mean_age_funev_cells_plot(self):
        return fs.is_file(self.mean_age_funev_cells_plot_path)

    # -----------------------------------------------------------------

    @property
    def mean_age_funev_cells_title(self):
        return "Correlation between mean stellar age and Funev of cells"

    # -----------------------------------------------------------------

    def plot_mean_age_funev_cells(self):

        """
        Thisn function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat:

            plot_stilts(self.mean_age_funev_cells_scatter_paths, self.mean_age_name, self.funev_name,
                        self.mean_age_description, self.funev_description,
                        title=self.mean_age_funev_cells_title, path=self.mean_age_funev_cells_plot_path,
                        ylimits=self.funev_limits, xlog=True, xlimits=self.mean_age_limits)

        # Plot using Matplotlib
        else: raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def mean_age_funev_pixels_scatters(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.mean_age_funev_pixels
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    @property
    def mean_age_funev_pixels_scatter_paths(self):
        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.mean_age_funev_pixels_path
        # no references (yet)
        return scatters

    # -----------------------------------------------------------------

    @property
    def mean_age_funev_pixels_plot_path(self):
        return fs.join(self.mean_age_funev_path, "cells.pdf")

    # -----------------------------------------------------------------

    @property
    def has_mean_age_funev_pixels_plot(self):
        return fs.is_file(self.mean_age_funev_pixels_plot_path)

    # -----------------------------------------------------------------

    @property
    def mean_age_funev_pixels_title(self):
        return "Correlation between mean stellar age and Funev of pixels"

    # -----------------------------------------------------------------

    def plot_mean_age_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat:

            plot_stilts(self.mean_age_funev_pixels_scatter_paths, self.mean_age_name, self.funev_name,
                        self.mean_age_description, self.funev_description,
                        title=self.mean_age_funev_pixels_title, path=self.mean_age_funev_pixels_plot_path,
                        ylimits=self.funev_limits, xlog=True, xlimits=self.mean_age_limits)

        # Plot using Matplotlib
        else: raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    def plot_mean_age_ssfr(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

def create_pixel_scatter(x_name, y_name, x_frame, y_frame, x_description, y_description, x_unit=None, y_unit=None,
                         aux=None, aux_units=None, aux_is_arrays=False, same_units=False, convolve=False):

    """
    This function ...
    :param x_name:
    :param y_name:
    :param x_frame:
    :param y_frame:
    :param x_description:
    :param y_description:
    :param x_unit:
    :param y_unit:
    :param aux:
    :param aux_units:
    :param aux_is_arrays:
    :param same_units:
    :param convolve:
    :return:
    """

    # Get units
    if x_unit is None: x_unit = x_frame.unit
    if y_unit is None: y_unit = y_frame.unit

    # Rebin to the same pixelscale
    x_frame, y_frame = uniformize(x_frame, y_frame, convolve=convolve, convert=same_units)

    # Get values
    x_values = x_frame.values
    y_values = y_frame.values

    # Set aux
    if aux is not None and not aux_is_arrays: aux = OrderedDict((key, data.values) for key, data in aux.items())

    # Create scatter
    return create_scatter(x_name, y_name, x_values, y_values, x_description, y_description, x_unit, y_unit, aux=aux, aux_units=aux_units)

# -----------------------------------------------------------------

def create_cell_scatter(x_name, y_name, x_data, y_data, x_description, y_description, x_unit=None, y_unit=None,
                        aux=None, aux_units=None, is_arrays=False, aux_is_arrays=False):

    """
    This function ....
    :param x_name:
    :param y_name:
    :param x_data:
    :param y_data:
    :param x_description:
    :param y_description:
    :param x_unit:
    :param y_unit:
    :param aux:
    :param aux_units:
    :param is_arrays:
    :param aux_is_arrays:
    :return:
    """

    # Get units
    if x_unit is None and not is_arrays: x_unit = x_data.unit
    if y_unit is None and not is_arrays: y_unit = y_data.unit

    # Check dimensions
    if len(x_data) != len(y_data): raise ValueError("x and y data must have the same size")

    # Get values
    if is_arrays:
        x_values = x_data
        y_values = y_data
    else:
        x_values = x_data.values
        y_values = y_data.values

    # Set aux
    if aux is not None and not aux_is_arrays: aux = OrderedDict((key, data.values) for key, data in aux.items())

    # Create
    return create_scatter(x_name, y_name, x_values, y_values, x_description, y_description, x_unit, y_unit, aux=aux, aux_units=aux_units)

# -----------------------------------------------------------------

def create_scatter(x_name, y_name, x_values, y_values, x_description=None, y_description=None, x_unit=None, y_unit=None, aux=None, aux_units=None):

    """
    Thisf unction ...
    :param x_name:
    :param y_name:
    :param x_values:
    :param y_values:
    :param x_description:
    :param y_description:
    :param x_unit:
    :param y_unit:
    :param aux:
    :param aux_units:
    :return:
    """

    # Create mask of valid values
    valid_ssfr_mask = np.isfinite(x_values) * (x_values != 0)
    valid_funev_mask = np.isfinite(y_values) * (y_values != 0)
    valid_mask = valid_ssfr_mask * valid_funev_mask

    # Get valid values
    valid_x = x_values[valid_mask]
    valid_y = y_values[valid_mask]

    # Get valid aux values
    if aux is not None: aux = OrderedDict((key, values[valid_mask]) for key, values in aux.items())

    # Create and return
    return Scatter2D.from_xy(valid_x, valid_y, x_name=x_name, y_name=y_name, x_unit=x_unit, y_unit=y_unit,
                             x_description=x_description, y_description=y_description, aux=aux, aux_units=aux_units)

# -----------------------------------------------------------------
