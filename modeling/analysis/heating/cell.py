#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.celldustheating Contains the CellDustHeatingAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.basics.distribution import Distribution, Distribution2D
from .tables import AbsorptionTable
from ....core.tools.utils import lazyproperty
from ....core.plot.distribution import plot_distribution, plot_distributions
from ....core.basics.containers import DefaultOrderedDict

# -----------------------------------------------------------------

class CellDustHeatingAnalyser(DustHeatingAnalysisComponent):
    
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
        super(CellDustHeatingAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The table with the absorbed luminosities
        self.absorptions = None

        # The heating fraction of the unevolved stellar population for each dust cell
        self.heating_fractions = None
        self.diffuse_heating_fractions = None

        # The distribution of heating fractions
        self.distribution = None
        self.distribution_diffuse = None

        # The 2D distribution of heating fractions
        self.radial_distribution = None

    # -----------------------------------------------------------------

    @property
    def needs_absorption_table(self):

        """
        This function ...
        :return:
        """

        #return (not self.has_heating_fractions) or (not self.has_heating_fractions_diffuse)
        return True # x and y coordinates are required to create radial distribution

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the absorption table
        if self.needs_absorption_table: self.get_absorption_table()

        # 3. Calculate the heating fraction of the unevolved stellar population
        self.get_fractions()

        # 4. Calculate the distribution of the heating fraction of the unevolved stellar population
        self.get_distributions()

        # 6. Writing
        self.write()

        # 7. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(CellDustHeatingAnalyser, self).setup(**kwargs)

        if self.config.recreate_table: self.config.recalculate_fractions = True
        if self.config.recalculate_fractions: self.config.recalculate_distributions = True

        # Set flags
        if self.config.recalculate_distributions:
            self.config.recalculate_distribution = True
            self.config.recalculate_distribution_diffuse = True
            self.config.recalculate_radial_distribution = True

        # Set replot flags
        if self.config.replot:
            self.config.replot_distribution = True
            self.config.replot_radial_distribution = True
            self.config.replot_map = True

        if self.config.recalculate_distribution or self.config.recalculate_distribution_diffuse: self.config.replot_distribution = True
        if self.config.recalculate_radial_distribution: self.config.replot_radial_distribution = True
        if self.config.recalculate_fractions: self.config.replot_map = True

        # Remove
        if self.config.recreate_table and self.has_absorptions: self.remove_absorptions()
        if self.config.recalculate_fractions and self.has_heating_fractions: self.remove_heating_fractions()
        if self.config.recalculate_fractions and self.has_heating_fractions_diffuse: self.remove_heating_fractions_diffuse()
        if self.config.recalculate_distribution and self.has_distribution: self.remove_distribution()
        if self.config.recalculate_distribution_diffuse and self.has_distribution_diffuse: self.remove_distribution_diffuse()
        if self.config.recalculate_radial_distribution and self.has_radial_distribution: self.remove_radial_distribution()
        if self.config.replot_distribution and self.has_distribution_plot: self.remove_distribution_plot()
        if self.config.replot_radial_distribution and self.has_radial_distribution_plot: self.remove_radial_distribution_plot()
        if self.config.replot_map and self.has_heating_map_plot: self.remove_heating_map_plot()

    # -----------------------------------------------------------------

    @lazyproperty
    def ncells(self):

        """
        This function ...
        :return:
        """

        return len(self.total_contribution_absorption_data)

    # -----------------------------------------------------------------

    def get_absorption_table(self):

        """
        This function ...
        :return:
        """

        # Create
        if not self.has_absorptions: self.create_absorption_table()

        # Load the table
        else: self.load_absorption_table()

    # -----------------------------------------------------------------

    def create_absorption_table(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Creating the absorption table ...")

        # Get the coordinates (should be in pc)
        x = self.cell_x_coordinates
        y = self.cell_y_coordinates
        z = self.cell_z_coordinates

        # Get luminosity for total stellar population
        #total_absorptions = self.total_contribution_absorption_data["Absorbed bolometric luminosity"]
        total_absorptions = self.total_contribution_absorption_luminosities
        total_unit = self.total_contribution_absorption_unit
        total_conversion = total_unit.conversion_factor("W")
        total_absorptions_watt = total_absorptions * total_conversion

        # Get luminosity for old stellar population
        #old_absorptions = self.old_contribution_absorption_data["Absorbed bolometric luminosity"]
        old_absorptions = self.old_contribution_absorption_luminosities
        old_unit = self.old_contribution_absorption_unit
        old_conversion = old_unit.conversion_factor("W")
        old_absorptions_watt = old_absorptions * old_conversion

        # Get luminosity for young stellar population
        #young_absorptions = self.young_contribution_absorption_data["Absorbed bolometric luminosity"]
        young_absorptions = self.young_contribution_absorption_luminosities
        young_unit = self.young_contribution_absorption_unit
        young_conversion = young_unit.conversion_factor("W")
        young_absorptions_watt = young_absorptions * young_conversion

        # Get luminosity for ionizing stellar population
        #ionizing_absorptions = self.ionizing_contribution_absorption_data["Absorbed bolometric luminosity"]
        ionizing_absorptions = self.ionizing_contribution_absorption_luminosities
        ionizing_unit = self.ionizing_contribution_absorption_unit
        ionizing_conversion = ionizing_unit.conversion_factor("W")
        ionizing_absorptions_watt = ionizing_absorptions * ionizing_conversion

        # Create the table
        self.absorptions = AbsorptionTable.from_columns(x, y, z, total_absorptions_watt, old_absorptions_watt, young_absorptions_watt, ionizing_absorptions_watt)

        # TEMP? SAVE
        #self.absorptions.saveto(self.absorption_table_path)

    # -----------------------------------------------------------------

    def load_absorption_table(self):

        """
        This function ...
        :return:
        """

        # Success
        log.success("Absorption table has already been created: loading from file ...")

        # Load
        self.absorptions = AbsorptionTable.from_file(self.absorption_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def zero_absorption_mask(self):

        """
        This function ...
        :return:
        """

        return self.absorptions.zero_absorption_mask

    # -----------------------------------------------------------------

    def get_fractions(self):

        self.get_heating_fractions()

        self.get_heating_fractions_diffuse()

    # -----------------------------------------------------------------

    def get_heating_fractions(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_heating_fractions: self.load_heating_fractions()

        # Calculate
        else: self.calculate_heating_fractions()

    # -----------------------------------------------------------------

    def load_heating_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the heating fractions ...")

        # Load
        self.heating_fractions = np.loadtxt(self.heating_fractions_path)

    # -----------------------------------------------------------------

    def calculate_heating_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the heating fraction of the unevolved stellar population ...")

        # SEBA:

        ## Total energy absorbed in the new component
        ## Derived from ../../SKIRTrun/models/testHeating/MappingsHeating/plotSEDs.py
        #Lnew = 176495776.676  # in Lsun

        #energy_new = volume * density_new * Lnew
        #F_abs_yng = (yng + new + energy_new) / (old + yng + new + energy_new)

        #cell_properties = self.model.cell_properties
        #volumes = cell_properties["Volume"]

        # For intrinsic SFR heating
        volumes = self.model.cell_volumes
        density = self.model.sfr_cell_stellar_density
        absorbed_energy = self.model.intrinsic_sfr_dust_luminosity
        print("absorbed energy", absorbed_energy)
        absorbed_energy_watt = absorbed_energy.to("W", distance=self.galaxy_distance).value
        print("absorbed energy (watt)", absorbed_energy_watt)

        #absorptions_unevolved_diffuse = self.absorptions["Absorbed bolometric luminosity of the young stellar population"] + self.absorptions["Absorbed bolometric luminosity of the ionizing stellar population"]
        #absorptions_unevolved_diffuse = self.absorptions["young"] + self.absorptions["ionizing"]

        #absorptions_unevolved_diffuse = self.absorptions.unevolved(unit="W", asarray=True)

        young_unit = self.absorptions.column_unit("young")
        ionizing_unit = self.absorptions.column_unit("ionizing")
        young_conversion = young_unit.conversion_factor("W")
        ionizing_conversion = ionizing_unit.conversion_factor("W")
        print("young conversion", young_conversion)
        print("ionizing conversion", ionizing_conversion)

        # Unevolved,  diffuse dust
        absorptions_young = np.asarray(self.absorptions["young"])
        absorptions_ionizing = np.asarray(self.absorptions["ionizing"])
        absorptions_young_watt = absorptions_young * young_conversion
        absorptions_ionizing_watt = absorptions_ionizing * ionizing_conversion
        #absorptions_unevolved_diffuse = self.absorptions[""]
        absorptions_unevolved_diffuse_watt = absorptions_young_watt + absorptions_ionizing_watt

        # Mass of ionizing stars for each cell
        relative_mass_ionizing = volumes * density
        normalization = sum(relative_mass_ionizing)
        print("norm", normalization)
        relative_mass_ionizing /= normalization

        #
        #absorptions_ionizing_internal_watt = volumes * density * absorbed_energy
        absorptions_ionizing_internal_watt = relative_mass_ionizing * absorbed_energy_watt

        #absorptions_unevolved = absorptions_unevolved_diffuse + absorptions_ionizing_internal
        absorptions_unevolved_watt = absorptions_unevolved_diffuse_watt + absorptions_ionizing_internal_watt

        #absorptions_total = self.absorptions["Absorbed bolometric luminosity of the total stellar population"]
        #absorptions_total = self.absorptions["total"]
        #absorptions_total = self.absorptions.total(unit="W", asarray=True)
        #absorptions_total = absorptions_unevolved_diffuse + absorptions_ionizing_internal + absorptions_evolved # TODO !!

        # Evolved: old
        #absorptions_evolved = self.absorptions.old(unit="W", asarray=True)
        old_unit = self.absorptions.column_unit("old")
        old_conversion = old_unit.conversion_factor("W")
        print("old conversion", old_conversion)
        absorptions_old = np.asarray(self.absorptions["old"])
        absorptions_old_watt = absorptions_old * old_conversion

        # Total absorptions
        #absorptions_total = absorptions_unevolved + absorptions_evolved

        # TOTAL
        absorptions_total_watt = absorptions_unevolved_watt + absorptions_old_watt

        # Calculate the heating fraction of the unevolved stellar population in each dust cell
        #self.heating_fractions = absorptions_unevolved_diffuse / absorptions_total
        #self.heating_fractions = absorptions_unevolved / absorptions_total
        self.heating_fractions = absorptions_unevolved_watt / absorptions_total_watt

    # -----------------------------------------------------------------

    def get_heating_fractions_diffuse(self):

        """
        This function ...
        :return:
        """

        if self.has_heating_fractions_diffuse: self.load_heating_fractions_diffuse()
        else: self.calculate_heating_fractions_diffuse()

    # -----------------------------------------------------------------

    def load_heating_fractions_diffuse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the diffuse heating fractions ...")

        # Load
        self.diffuse_heating_fractions = np.loadtxt(self.heating_fractions_diffuse_path)

    # -----------------------------------------------------------------

    def calculate_heating_fractions_diffuse(self):

        """
        This function ...
        :return:
        """

        young_unit = self.absorptions.column_unit("young")
        ionizing_unit = self.absorptions.column_unit("ionizing")
        young_conversion = young_unit.conversion_factor("W")
        ionizing_conversion = ionizing_unit.conversion_factor("W")
        print("young conversion", young_conversion)
        print("ionizing conversion", ionizing_conversion)

        # Unevolved,  diffuse dust
        absorptions_young = np.asarray(self.absorptions["young"])
        absorptions_ionizing = np.asarray(self.absorptions["ionizing"])
        absorptions_young_watt = absorptions_young * young_conversion
        absorptions_ionizing_watt = absorptions_ionizing * ionizing_conversion
        # absorptions_unevolved_diffuse = self.absorptions[""]
        absorptions_unevolved_diffuse_watt = absorptions_young_watt + absorptions_ionizing_watt

        total_unit = self.absorptions.column_unit("total")
        total_conversion = total_unit.conversion_factor("W")
        absorptions_total_diffuse = np.asarray(self.absorptions["total"])
        absorptions_total_diffuse_watt = absorptions_total_diffuse * total_conversion

        # Diffuse
        self.diffuse_heating_fractions = absorptions_unevolved_diffuse_watt / absorptions_total_diffuse_watt

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_fraction_nans(self):

        """
        Thisn function ...
        :return:
        """

        return np.isnan(self.heating_fractions)

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_heating_fraction_nans(self):

        """
        This function ...
        :return:
        """

        return np.isnan(self.diffuse_heating_fractions)

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_fraction_infs(self):

        """
        This function ...
        :return:
        """

        return np.isinf(self.heating_fractions)

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_heating_fraction_infs(self):

        """
        This function ...
        :return:
        """

        return np.isinf(self.diffuse_heating_fractions)

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_fraction_unphysical(self):

        """
        Thisn function ...
        :return:
        """

        greater_than_one_mask = self.heating_fractions > 1.0
        ngreater_than_one = np.sum(greater_than_one_mask)
        relative_ngreater_than_one = float(ngreater_than_one) / len(self.heating_fractions)

        # Debugging
        log.debug(str(ngreater_than_one) + " pixels have a heating fraction greater than unity (" + str(relative_ngreater_than_one * 100) + "%)")

        # Return
        return greater_than_one_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_heating_fraction_unphysical(self):

        """
        This function ...
        :return:
        """

        greater_than_one_mask = self.diffuse_heating_fractions > 1.0
        #ngreater_than_one = np.sum(greater_than_one_mask)
        #relative_ngreater_than_one = float(ngreater_than_one) / len(self.diffuse_heating_fractions)

        # Debugging
        #log.debug(str(ngreater_than_one) + " pixels have a heating fraction greater than unity (" + str(
        #    relative_ngreater_than_one * 100) + "%)")

        # Return
        return greater_than_one_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_fractions_mask(self):

        """
        This function ...
        :return:
        """

        return self.heating_fraction_nans + self.heating_fraction_infs + self.heating_fraction_unphysical

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_heating_fractions_mask(self):

        """
        This function ...
        :return:
        """

        return self.diffuse_heating_fraction_nans + self.diffuse_heating_fraction_infs + self.diffuse_heating_fraction_unphysical

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_heating_fractions(self):

        """
        This function ...
        :return:
        """

        return np.ma.MaskedArray(self.heating_fractions, self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_heating_fractions_diffuse(self):

        """
        This function ...
        :return:
        """

        return np.ma.MaskedArray(self.diffuse_heating_fractions, self.diffuse_heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_weights(self):

        """
        This function ...
        :return:
        """

        return np.asarray(self.cell_properties["Mass fraction"])

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_weights(self):

        """
        This function ...
        :return:
        """

        return np.ma.MaskedArray(self.cell_weights, self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_weights_diffuse(self):

        """
        This function ...
        :return:
        """

        return np.ma.MaskedArray(self.cell_weights, self.diffuse_heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def radii(self):

        """
        This function ...
        :return:
        """

        return np.sqrt(self.x_coordinates**2 + self.y_coordinates**2)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_radii(self):

        """
        This fuction ...
        :return:
        """

        return np.ma.MaskedArray(self.radii, self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def x_coordinates(self):

        """
        Thisn function ...
        :return:
        """

        return np.asarray(self.absorptions["x"])

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_x_coordinates(self):

        """
        This function ...
        :return:
        """

        return np.ma.MaskedArray(self.x_coordinates, mask=self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def y_coordinates(self):

        """
        Thisn function ...
        :return:
        """

        return np.asarray(self.absorptions["y"])

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_y_coordinates(self):

        """
        Thisn function ...
        :return:
        """

        return np.ma.MaskedArray(self.y_coordinates, mask=self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def z_coordinates(self):

        """
        Thisnf unction ...
        :return:
        """

        return np.asarray(self.absorptions["z"])

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_z_coordinates(self):

        """
        This function ...
        :return:
        """

        return np.ma.MaskedArray(self.z_coordinates, mask=self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @property
    def needs_distribution(self):

        """
        This function ...
        :return:
        """

        return self.do_plot_distribution

    # -----------------------------------------------------------------

    @property
    def needs_distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        return self.do_plot_distribution

    # -----------------------------------------------------------------

    @property
    def needs_radial_distribution(self):

        """
        This function ...
        :return:
        """

        return self.do_plot_radial_distribution

    # -----------------------------------------------------------------

    def get_distributions(self):

        """
        This function ...
        :return:
        """

        # Distribution
        if self.needs_distribution: self.get_distribution()

        # Distribution diffuse
        if self.needs_distribution_diffuse: self.get_distribution_diffuse()

        # Radial distribution
        if self.needs_radial_distribution: self.get_radial_distribution()

    # -----------------------------------------------------------------

    def get_distribution(self):

        """
        This function ...
        :return:
        """

        if self.has_distribution: self.load_distribution()
        else: self.calculate_distribution()

    # -----------------------------------------------------------------

    def load_distribution(self):

        """
        This function ...
        :return:
        """

        self.distribution = Distribution.from_file(self.distribution_path)

    # -----------------------------------------------------------------

    def get_distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        if self.has_distribution_diffuse: self.load_distribution_diffuse()
        else: self.calculate_distribution_diffuse()

    # -----------------------------------------------------------------

    def load_distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        self.distribution_diffuse = Distribution.from_file(self.distribution_diffuse_path)

    # -----------------------------------------------------------------

    def get_radial_distribution(self):

        """
        This function ...
        :return:
        """

        #if self.has_radial_distribution: self.load_radial_distribution()
        #else: self.calculate_radial_distribution()

        # Cannot save the radial distribution for the moment
        self.calculate_radial_distribution()

    # -----------------------------------------------------------------

    def calculate_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the distribution of heating fractions of the unevolved stellar population ...")

        # Generate the distribution
        # Weights are dust mass fraction
        self.distribution = Distribution.from_values("Heating fraction", self.valid_heating_fractions,
                                                     nbins=self.config.nbins, weights=self.valid_cell_weights)

    # -----------------------------------------------------------------

    def calculate_distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the distribution of heating fractions of the unevolved stellar population (diffuse) ...")

        # Generate the distribution
        # Weights are dust mass fraction
        self.distribution_diffuse = Distribution.from_values("Heating fraction", self.valid_heating_fractions_diffuse,
                                                             nbins=self.config.nbins, weights=self.valid_cell_weights_diffuse)

    # -----------------------------------------------------------------

    def calculate_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the radial distribution of heating fractions of the unevolved stellar population ...")

        # Generate the radial distribution
        self.radial_distribution = Distribution2D.from_values(self.valid_radii, self.valid_heating_fractions,
                                                              weights=self.valid_cell_weights, x_name="radius (pc)",
                                                              y_name="Heating fraction of unevolved stars", nbins=self.config.nradial_bins)

    # -----------------------------------------------------------------

    @property
    def do_write_absorptions(self):

        """
        This function ...
        :return:
        """

        return not self.has_absorptions

    # -----------------------------------------------------------------

    @property
    def do_write_heating_fractions(self):

        """
        This function ...
        :return:
        """

        return not self.has_heating_fractions

    # -----------------------------------------------------------------

    @property
    def do_write_heating_fractions_diffuse(self):

        """
        This function ...
        :return:
        """

        return not self.has_heating_fractions_diffuse

    # -----------------------------------------------------------------

    @property
    def do_write_distribution(self):

        """
        This function ...
        :return:
        """

        return not self.has_distribution

    # -----------------------------------------------------------------

    @property
    def do_write_distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        return not self.has_distribution_diffuse

    # -----------------------------------------------------------------

    @property
    def do_write_radial_distribution(self):

        """
        This function ...
        :return:
        """

        return not self.has_radial_distribution

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the absorption table
        if self.do_write_absorptions: self.write_absorptions()

        # Write the heating fractions
        if self.do_write_heating_fractions: self.write_heating_fractions()

        # Write the diffuse heating fractions
        if self.do_write_heating_fractions_diffuse: self.write_heating_fractions_diffuse()

        # Write the distribution of heating fractions
        if self.do_write_distribution: self.write_distribution()

        # Write diffuse distribution
        if self.do_write_distribution_diffuse: self.write_distribution_diffuse()

        # Write the radial distribution of heating fractions
        if self.do_write_radial_distribution: self.write_radial_distribution()

    # -----------------------------------------------------------------

    @property
    def absorption_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "absorptions.dat")

    # -----------------------------------------------------------------

    @property
    def has_absorptions(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.absorption_table_path)

    # -----------------------------------------------------------------

    def remove_absorptions(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.absorption_table_path)

    # -----------------------------------------------------------------

    def write_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the absorption table ...")

        # Save the table
        self.absorptions.saveto(self.absorption_table_path)

    # -----------------------------------------------------------------

    @property
    def heating_fractions_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "heating_fractions.dat")

    # -----------------------------------------------------------------

    @property
    def has_heating_fractions(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.heating_fractions_path)

    # -----------------------------------------------------------------

    def remove_heating_fractions(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.heating_fractions_path)

    # -----------------------------------------------------------------

    def write_heating_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the heating fractions ...")

        # Write
        np.savetxt(self.heating_fractions_path, self.heating_fractions)

    # -----------------------------------------------------------------

    @property
    def heating_fractions_diffuse_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "heating_fractions_diffuse.dat")

    # -----------------------------------------------------------------

    @property
    def has_heating_fractions_diffuse(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.heating_fractions_diffuse_path)

    # -----------------------------------------------------------------

    def remove_heating_fractions_diffuse(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.heating_fractions_diffuse_path)

    # -----------------------------------------------------------------

    def write_heating_fractions_diffuse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the diffuse heating fractions ...")

        # Write
        np.savetxt(self.heating_fractions_diffuse_path, self.diffuse_heating_fractions)

    # -----------------------------------------------------------------

    @property
    def distribution_path(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the distribution file
        return fs.join(self.cell_heating_path, "distribution.dat")

    # -----------------------------------------------------------------

    @property
    def has_distribution(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.distribution_path)

    # -----------------------------------------------------------------

    def remove_distribution(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.distribution_path)

    # -----------------------------------------------------------------

    def write_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the distribution of heating fractions of the unevolved stellar population ...")

        # Save the distribution
        self.distribution.saveto(self.distribution_path)

    # -----------------------------------------------------------------

    @property
    def distribution_diffuse_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "distribution_diffuse.dat")

    # -----------------------------------------------------------------

    @property
    def has_distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.distribution_diffuse_path)

    # -----------------------------------------------------------------

    def remove_distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.distribution_diffuse_path)

    # -----------------------------------------------------------------

    def write_distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the distribution of heating fractions (diffuse) ...")

        # Save
        self.distribution_diffuse.saveto(self.distribution_diffuse_path)

    # -----------------------------------------------------------------

    @property
    def radial_distribution_path(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the radial distribution file
        return fs.join(self.cell_heating_path, "radial_distribution.dat")

    # -----------------------------------------------------------------

    @property
    def has_radial_distribution(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.radial_distribution_path)

    # -----------------------------------------------------------------

    def remove_radial_distribution(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.radial_distribution_path)

    # -----------------------------------------------------------------

    def write_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the radial distribution of heating fractions of the unevolved stellar population ...")

        # Save the radial distribution
        self.radial_distribution.saveto(self.radial_distribution_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_distribution(self):

        """
        This function ...
        :return:
        """

        return self.config.plot_distribution and not self.has_distribution_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_radial_distribution(self):

        """
        This function ...
        :return:
        """

        return self.config.plot_radial_distribution and not self.has_radial_distribution_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_map(self):

        """
        This function ...
        :return:
        """

        return self.config.plot_map and not self.has_heating_map_plot

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the distribution of heating fractions
        if self.do_plot_distribution: self.plot_distribution()

        # Plot the radial distribution of heating fractions
        if self.do_plot_radial_distribution: self.plot_radial_distribution()

        # Plot a map of the heating fraction for a face-on view of the galaxy
        if self.do_plot_map: self.plot_map()

    # -----------------------------------------------------------------

    @property
    def distribution_plot_path(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the plot file
        return fs.join(self.cell_heating_path, "distribution.pdf")

    # -----------------------------------------------------------------

    @property
    def has_distribution_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.distribution_plot_path)

    # -----------------------------------------------------------------

    def remove_distribution_plot(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.distribution_plot_path)

    # -----------------------------------------------------------------

    def plot_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a histogram of the distribution of heating fractions of the unevolved stellar population ...")

        # Create the plot file
        #self.distribution.plot(title="Distribution of the heating fraction of the unevolved stellar population", path=path)
        #plot_distribution(self.distribution, title="Distribution of the heating fraction of the unevolved stellar population", path=self.distribution_plot_path)

        distributions = OrderedDict()
        distributions["unevolved"] = self.distribution
        distributions["unevolved (diffuse)"] = self.distribution_diffuse

        # Plot
        plot_distributions(distributions, path=self.distribution_plot_path, alpha=0.5)

    # -----------------------------------------------------------------

    @property
    def radial_distribution_plot_path(self):

        """
        Thisn function ...
        :return:
        """

        # Determine the path to the plot file
        return fs.join(self.cell_heating_path, "radial_distribution.pdf")

    # -----------------------------------------------------------------

    @property
    def has_radial_distribution_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.radial_distribution_plot_path)

    # -----------------------------------------------------------------

    def remove_radial_distribution_plot(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.radial_distribution_plot_path)

    # -----------------------------------------------------------------

    def plot_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a 2D histogram of the radial distribution of the heating fractions of the unevolved stellar population ...")

        # Radii
        # FOR M81
        #radii = [2.22953938405 * u("kpc"), 3.34430907608 * u("kpc"), 5.90827936773 * u("kpc"), 8.9181575362 * u("kpc")]
        radii = [2.22953938405, 3.34430907608, 5.90827936773, 8.9181575362]  # kpc
        radii_pc = [radius*1000 for radius in radii]

        # Set title
        title = "Radial distribution of the heating fraction of the unevolved stellar population"

        # Create the plot file
        self.radial_distribution.plot(radii=radii_pc, title=title, path=self.radial_distribution_plot_path)

    # -----------------------------------------------------------------

    @property
    def heating_map_plot_path(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the plot file
        return fs.join(self.cell_heating_path, "map.pdf")

    # -----------------------------------------------------------------

    @property
    def has_heating_map_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.heating_map_plot_path)

    # -----------------------------------------------------------------

    def remove_heating_map_plot(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.heating_map_plot_path)

    # -----------------------------------------------------------------

    def plot_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a map of the heating fraction of the unevolved stellar population for a face-on view of the galaxy ...")

        # Create figure
        plt.figure()

        x_coordinates = self.valid_x_coordinates
        y_coordinates = self.valid_y_coordinates
        fractions = self.valid_heating_fractions
        weights = self.valid_cell_weights

        z_indices = DefaultOrderedDict(list)

        ncoordinates = len(x_coordinates)
        for i in range(ncoordinates):
            x = x_coordinates[i]
            y = y_coordinates[i]
            z_indices[(x,y)].append(i)

        #print(x, x.shape)
        #print(y, y.shape)
        #print(z, z.shape)

        projected_x = []
        projected_y = []
        projected_fractions = []

        for x,y in z_indices:

            indices = z_indices[(x,y)]

            z_values = self.valid_z_coordinates[indices]

            fractions_column = fractions[indices]
            weights_column = weights[indices]

            normalization = np.sum(fractions_column)
            if normalization == 0: continue

            fraction = np.sum(fractions_column * weights_column) / normalization

            projected_x.append(x)
            projected_y.append(y)
            projected_fractions.append(fraction)

        #plt.pcolormesh(x, y, z, cmap='RdBu', vmin=0.0, vmax=1.0)

        ax = plt.gca()
        #ax.pcolormesh(x, y, fractions)

        ax.pcolormesh(projected_x, projected_y, projected_fractions)

        # Plot
        plt.savefig(self.heating_map_plot_path)
        plt.close()

    # -----------------------------------------------------------------

    def plot_map2(self):

        """
        This function ...
        :return:
        """

        from matplotlib import mlab

        x_ticks = x
        y_ticks = y

        z_grid = mlab.griddata(x, y, z, x, y)

        from matplotlib.backends import backend_agg as agg
        from matplotlib import cm

        # plot
        # fig = Figure()  # create the figure
        fig = plt.figure()
        agg.FigureCanvasAgg(fig)  # attach the rasterizer
        ax = fig.add_subplot(1, 1, 1)  # make axes to plot on
        ax.set_title("Interpolated Contour Plot of Experimental Data")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        cmap = cm.get_cmap("hot")  # get the "hot" color map
        contourset = ax.contourf(x_ticks, y_ticks, z_grid, 10, cmap=cmap)

        # Plot
        plt.savefig(self.heating_map_plot_path)
        plt.close()

# -----------------------------------------------------------------
