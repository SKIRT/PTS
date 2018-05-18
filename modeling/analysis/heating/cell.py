#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating.cell Contains the CellDustHeatingAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import matplotlib.mlab as mlab
from scipy.interpolate import griddata

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.basics.distribution import Distribution, Distribution2D
from .tables import AbsorptionTable
from ....core.tools.utils import lazyproperty
from ....core.plot.distribution import plot_distribution, plot_distributions
from ....core.basics.containers import DefaultOrderedDict
from ....magic.core.frame import Frame
from ....core.tools import nr
from ....magic.basics.vector import PixelShape
from ....core.basics.range import QuantityRange, RealRange
from ....core.tools.stringify import tostr
from ....core.tools import numbers, sequences
from ....core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

old = "old"
young = "young"
ionizing = "ionizing"
dust = "dust"
disk_components = [old, young, ionizing, dust]

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

        # The heating maps
        self.map = None
        self.map_stddev = None
        self.map_ncells = None
        self.map_midplane = None
        self.map_midplane_stddev = None
        self.map_midplane_ncells = None

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

        # 5. Create the maps of the heating fraction
        self.get_maps()

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
        if self.config.recreate_map and self.has_map: self.remove_map()
        if self.config.recreate_map_midplane and self.has_map_midplane: self.remove_map_midplane()
        if self.config.replot_distribution and self.has_distribution_plot: self.remove_distribution_plot()
        if self.config.replot_radial_distribution and self.has_radial_distribution_plot: self.remove_radial_distribution_plot()
        if self.config.replot_map and self.has_map_plot: self.remove_map_plot()
        if self.config.replot_map_midplane and self.has_map_midplane_plot: self.remove_map_midplane_plot()

        # Check for consistency
        if self.config.consistency:

            # Remove the map of heating fraction and/or stddev and/or ncells if either one of them is not present
            if sequences.some_true_but_not_all([self.has_map, self.has_map_stddev, self.has_map_ncells]): self.remove_map_stddev_and_ncells()

            # Remove the map of the heating fraction in the midplane and/or stddev and/or ncells if either one of them is not present
            if sequences.some_true_but_not_all([self.has_map_midplane, self.has_map_midplane_stddev, self.has_map_midplane_ncells]): self.remove_map_midplane_stddev_and_ncells()

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

    @property
    def needs_heating_fractions(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def needs_heating_fractions_diffuse(self):

        """
        This function ...
        :return:
        """

        return self.do_write_distribution_diffuse or self.do_plot_distribution

    # -----------------------------------------------------------------

    def get_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the heating fractions ...")

        # Total heating fractions
        if self.needs_heating_fractions: self.get_heating_fractions()

        # Diffuse heating fractions
        if self.needs_heating_fractions_diffuse: self.get_heating_fractions_diffuse()

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
        #print("absorbed energy", absorbed_energy)
        absorbed_energy_watt = absorbed_energy.to("W", distance=self.galaxy_distance).value
        #print("absorbed energy (watt)", absorbed_energy_watt)

        #absorptions_unevolved_diffuse = self.absorptions["Absorbed bolometric luminosity of the young stellar population"] + self.absorptions["Absorbed bolometric luminosity of the ionizing stellar population"]
        #absorptions_unevolved_diffuse = self.absorptions["young"] + self.absorptions["ionizing"]

        #absorptions_unevolved_diffuse = self.absorptions.unevolved(unit="W", asarray=True)

        young_unit = self.absorptions.column_unit("young")
        ionizing_unit = self.absorptions.column_unit("ionizing")
        young_conversion = young_unit.conversion_factor("W")
        ionizing_conversion = ionizing_unit.conversion_factor("W")
        #print("young conversion", young_conversion)
        #print("ionizing conversion", ionizing_conversion)

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
        #print("norm", normalization)
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
        #print("old conversion", old_conversion)
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
        #print("young conversion", young_conversion)
        #print("ionizing conversion", ionizing_conversion)

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
    def coordinates_unit(self):

        """
        This function ...
        :return:
        """

        return sequences.get_all_equal_value([self.x_coordinates_unit, self.y_coordinates_unit, self.z_coordinates_unit])

    # -----------------------------------------------------------------

    @lazyproperty
    def coordinates_unit_string(self):

        """
        Thins function ...
        :return:
        """

        return tostr(self.coordinates_unit)

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
    def x_coordinates_unit(self):

        """
        This function ...
        :return:
        """

        return self.absorptions.column_unit("x")

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
    def y_coordinates_unit(self):

        """
        This function ...
        :return:
        """

        return self.absorptions.column_unit("y")

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
    def z_coordinates_unit(self):

        """
        Thisn function ...
        :return:
        """

        return self.absorptions.column_unit("z")

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
    def needs_map(self):

        """
        This function ...
        :return:
        """

        return not self.has_map or self.do_plot_map

    # -----------------------------------------------------------------

    @property
    def needs_map_midplane(self):

        """
        This function ...
        :return:
        """

        return not self.has_map_midplane or self.do_plot_map_midplane

    # -----------------------------------------------------------------

    def get_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the maps of the heating fraction ...")

        # Get the map
        if self.needs_map: self.get_map()

        # Get the map in the midplane
        if self.needs_map_midplane: self.get_map_midplane()

    # -----------------------------------------------------------------

    def get_map(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_map: self.load_map()

        # Create
        else: self.create_map()

    # -----------------------------------------------------------------

    def load_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of the heating fraction ...")

        # Load
        self.map = Frame.from_file(self.map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_unique_x_coordinates(self):

        """
        This function ...
        :return:
        """

        return np.sort(np.unique(self.valid_x_coordinates))

    # -----------------------------------------------------------------

    @lazyproperty
    def min_x_coordinate(self):

        """
        This function ...
        :return:
        """

        return self.sorted_unique_x_coordinates[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def max_x_coordinate(self):

        """
        This function ...
        :return:
        """

        return self.sorted_unique_x_coordinates[-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_spacings(self):

        """
        This function ...
        :return:
        """

        return np.diff(self.sorted_unique_x_coordinates)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_average_spacing(self):

        """
        This function ...
        :return:
        """

        return np.mean(self.unique_x_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_median_spacing(self):

        """
        This function ...
        :return:
        """

        return np.median(self.unique_x_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_stddev_spacing(self):

        """
        This function ...
        :return:
        """

        return np.std(self.unique_x_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_min_spacing(self):

        """
        Thisf unction ...
        :return:
        """

        return np.min(self.unique_x_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_max_spacing(self):

        """
        This function ...
        :return:
        """

        return np.max(self.unique_x_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_unique_y_coordinates(self):

        """
        This function ...
        :return:
        """

        return np.sort(np.unique(self.valid_y_coordinates))

    # -----------------------------------------------------------------

    @lazyproperty
    def min_y_coordinate(self):

        """
        This function ...
        :return:
        """

        return self.sorted_unique_y_coordinates[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def max_y_coordinate(self):

        """
        This function ...
        :return:
        """

        return self.sorted_unique_y_coordinates[-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_spacings(self):

        """
        This function ...
        :return:
        """

        return np.diff(self.sorted_unique_y_coordinates)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_average_spacing(self):

        """
        This function ...
        :return:
        """

        return np.mean(self.unique_y_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_median_spacing(self):

        """
        This function ...
        :return:
        """

        return np.median(self.unique_y_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_stddev_spacing(self):

        """
        This function ...
        :return:
        """

        return np.std(self.unique_y_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_min_spacing(self):

        """
        Thisf unction ...
        :return:
        """

        return np.min(self.unique_y_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_max_spacing(self):

        """
        This function ...
        :return:
        """

        return np.max(self.unique_y_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def map_spacing(self):

        """
        This function ...
        :return:
        """

        if self.config.midplane_spacing_measure == "min": spacing = np.mean([self.unique_x_coordinates_min_spacing, self.unique_y_coordinates_min_spacing])
        elif self.config.midplane_spacing_measure == "max": spacing = np.mean([self.unique_x_coordinates_max_spacing, self.unique_y_coordinates_max_spacing])
        elif self.config.midplane_spacing_measure == "mean": spacing = np.mean([self.unique_x_coordinates_average_spacing, self.unique_y_coordinates_average_spacing])
        elif self.config.midplane_spacing_measure == "median": spacing = np.mean([self.unique_x_coordinates_median_spacing, self.unique_y_coordinates_median_spacing])
        else: raise ValueError("Invalid value for 'map_spacing_measure'")

        # Return
        return spacing * self.config.map_spacing_factor

    # -----------------------------------------------------------------

    @lazyproperty
    def map_radius(self):

        """
        This function ...
        :return:
        """

        return max(abs(self.min_x_coordinate), self.max_x_coordinate, abs(self.min_y_coordinate), self.max_y_coordinate)

    # -----------------------------------------------------------------

    @lazyproperty
    def map_coordinates(self):

        """
        This function ...
        :return:
        """

        return np.arange(-self.map_radius, self.map_radius, self.map_spacing)

    # -----------------------------------------------------------------

    @lazyproperty
    def half_map_spacing(self):

        """
        This function ...
        :return:
        """

        return 0.5 * self.map_spacing

    # -----------------------------------------------------------------

    def get_coordinate_range_around_pixel_map(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        x_range = RealRange(x-self.half_map_spacing,x+self.half_map_spacing)
        y_range = RealRange(y-self.half_map_spacing,y+self.half_map_spacing)
        return x_range, y_range

    # -----------------------------------------------------------------

    @lazyproperty
    def map_shape(self):

        """
        This function ...
        :return:
        """

        return PixelShape.square(self.nmap_coordinates)

    # -----------------------------------------------------------------

    @lazyproperty
    def nmap_coordinates(self):

        """
        This function ...
        :return:
        """

        return len(self.map_coordinates)

    # -----------------------------------------------------------------

    @property
    def map_npixels(self):

        """
        This function ...
        :return:
        """

        return self.map_shape.nxy

    # -----------------------------------------------------------------

    def get_coordinate_mask_x_for_map(self, x_range):

        """
        This function .....
        :param x_range:
        :return:
        """

        x_mask = (x_range.min < self.valid_x_coordinates) * (self.valid_x_coordinates <= x_range.max)
        return x_mask

    # -----------------------------------------------------------------

    def get_coordinate_mask_y_for_map(self, y_range):

        """
        This function ...
        :param y_range:
        :return:
        """

        y_mask = (y_range.min < self.valid_y_coordinates) * (self.valid_y_coordinates <= y_range.max)
        return y_mask

    # -----------------------------------------------------------------

    def get_coordinate_mask_for_map(self, x_range, y_range):

        """
        This function ...
        :param x_range:
        :param y_range:
        :return:
        """

        x_mask = self.get_coordinate_mask_x_for_map(x_range)
        y_mask = self.get_coordinate_mask_y_for_map(y_range)
        return x_mask * y_mask

    # -----------------------------------------------------------------

    def get_coordinate_indices_in_column_for_map(self, x_range, y_range):

        """
        This function ...
        :param x_range:
        :param y_range:
        :return:
        """

        mask = self.get_coordinate_mask_for_map(x_range, y_range)
        return np.where(mask)[0]

    # -----------------------------------------------------------------

    def create_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of the heating fraction ...")

        # Initialize map of the midplane heating fraction
        self.map = Frame.initialize_nans(self.map_shape)
        self.map_stddev = Frame.initialize_nans(self.map_shape)
        self.map_ncells = Frame.initialize_nans(self.map_shape)

        # Set the pixelscale and the coordinate info
        self.map.pixelscale = self.map_spacing * self.coordinates_unit
        self.map.distance = self.galaxy_distance
        self.map.set_meta("x_min", repr(-self.map_radius) + " " + self.coordinates_unit_string)
        self.map.set_meta("x_max", repr(self.map_radius) + " " + self.coordinates_unit_string)
        self.map.set_meta("y_min", repr(-self.map_radius) + " " + self.coordinates_unit_string)
        self.map.set_meta("y_max", repr(self.map_radius) + " " + self.coordinates_unit_string)

        # Loop over the x and y coordinates
        index = 1
        for i, x in enumerate(self.map_coordinates):
            for j, y in enumerate(self.map_coordinates):

                # Determine range of x and y
                x_range, y_range = self.get_coordinate_range_around_pixel_map(x, y)

                # Show
                if index % 100 == 0: log.debug("Calculating heating fraction in the pixel " + str(index) + " of " + str(self.map_npixels) + " (" + tostr(float(index) / self.map_npixels * 100, decimal_places=1, round=True) + "%) ...")

                # Get the indices
                indices = self.get_coordinate_indices_in_column_for_map(x_range, y_range)
                nindices = indices.shape[0]

                # Set number of cells
                self.map_ncells[j, i] = nindices

                # If any cells
                if nindices > 0:

                    # Calculate the heating fraction
                    fractions = self.valid_heating_fractions[indices]
                    weights = self.valid_cell_weights[indices]

                    # Calculate the mean heating fraction
                    fraction = numbers.weighed_arithmetic_mean_numpy(fractions, weights)
                    fraction_stddev = numbers.weighed_standard_deviation_numpy(fractions, weights, mean=fraction)

                    # Set fraction
                    self.map[j, i] = fraction
                    self.map_stddev[j, i] = fraction_stddev

                # Increment
                index += 1

    # -----------------------------------------------------------------

    @lazyproperty
    def absolute_z_coordinates(self):

        """
        This function ...
        :return:
        """

        return abs(self.valid_z_coordinates)

    # -----------------------------------------------------------------

    @lazyproperty
    def smallest_absolute_z(self):

        """
        This function ...
        :return:
        """

        return min(self.absolute_z_coordinates)

    # -----------------------------------------------------------------

    @lazyproperty
    def where_smallest_absolute_z(self):

        """
        This function ...
        :return:
        """

        return self.absolute_z_coordinates == self.smallest_absolute_z

    # -----------------------------------------------------------------

    @lazyproperty
    def ncells_with_smallest_absolute_z(self):

        """
        This function ...
        :return:
        """

        return int(np.sum(self.where_smallest_absolute_z))

    # -----------------------------------------------------------------

    @lazyproperty
    def midplane_component_scaleheight(self):

        """
        This function ...
        :return:
        """

        # Old stellar disk
        if self.config.midplane_component == old: return self.model.old_disk_scaleheight

        # Young stellar disk
        elif self.config.midplane_component == young: return self.model.young_scaleheight

        # Ionizing stellar disk
        elif self.config.midplane_component == ionizing: return self.model.sfr_scaleheight

        # Dust disk
        elif self.config.midplane_component == dust: return self.model.dust_scaleheight

        # Invalid
        else: raise ValueError("Invalid midplane component: '" + self.config.midplane_component + "'")

    # -----------------------------------------------------------------

    @lazyproperty
    def midplane_height(self):

        """
        This function ...
        :return:
        """

        return self.midplane_component_scaleheight * self.config.midplane_factor

    # -----------------------------------------------------------------

    @lazyproperty
    def midplane_height_pc(self):

        """
        This function ...
        :return:
        """

        return self.midplane_height.to("pc").value

    # -----------------------------------------------------------------

    @lazyproperty
    def where_midplane_z(self):

        """
        This function ...
        :return:
        """

        return self.absolute_z_coordinates <= self.midplane_height_pc

    # -----------------------------------------------------------------

    @lazyproperty
    def ncells_midplane(self):

        """
        Thisn function ...
        :return:
        """

        return int(np.sum(self.where_midplane_z))

    # -----------------------------------------------------------------

    @lazyproperty
    def x_coordinates_midplane(self):

        """
        This function ...
        :return:
        """

        return self.valid_x_coordinates[self.where_midplane_z]
        #return self.valid_x_coordinates[self.where_smallest_absolute_z]

    # -----------------------------------------------------------------

    @lazyproperty
    def y_coordinates_midplane(self):

        """
        This function ...
        :return:
        """

        return self.valid_y_coordinates[self.where_midplane_z]
        #return self.valid_y_coordinates[self.where_smallest_absolute_z]

    # -----------------------------------------------------------------

    @lazyproperty
    def z_coordinates_midplane(self):

        """
        This function ...
        :return:
        """

        return self.valid_z_coordinates[self.where_midplane_z]
        #return self.valid_z_coordinates[self.where_smallest_absolute_z]

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_fractions_midplane(self):

        """
        This function ...
        :return:
        """

        return self.valid_heating_fractions[self.where_midplane_z]
        #return self.valid_heating_fractions[self.where_smallest_absolute_z]

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_weights_midplane(self):

        """
        This function ...
        :return:
        """

        return self.valid_cell_weights[self.where_midplane_z]

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_unique_x_coordinates_midplane(self):

        """
        This function ...
        :return:
        """

        return np.sort(np.unique(self.x_coordinates_midplane))

    # -----------------------------------------------------------------

    @lazyproperty
    def min_x_coordinate_midplane(self):

        """
        This function ...
        :return:
        """

        return self.sorted_unique_x_coordinates_midplane[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def max_x_coordinate_midplane(self):

        """
        This function ...
        :return:
        """

        return self.sorted_unique_x_coordinates_midplane[-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_midplane_spacings(self):

        """
        This function ...
        :return:
        """

        return np.diff(self.sorted_unique_x_coordinates_midplane)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_midplane_average_spacing(self):

        """
        This function ...
        :return:
        """

        return np.mean(self.unique_x_coordinates_midplane_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_midplane_median_spacing(self):

        """
        Thisf unction ...
        :return:
        """

        return np.median(self.unique_x_coordinates_midplane_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_midplane_stddev_spacing(self):

        """
        This function ...
        :return:
        """

        return np.std(self.unique_x_coordinates_midplane_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_midplane_min_spacing(self):

        """
        This function ...
        :return:
        """

        return np.min(self.unique_x_coordinates_midplane_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_midplane_max_spacing(self):

        """
        This function ...
        :return:
        """

        return np.max(self.unique_x_coordinates_midplane_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_unique_y_coordinates_midplane(self):

        """
        This function ...
        :return:
        """

        return np.sort(np.unique(self.y_coordinates_midplane))

    # -----------------------------------------------------------------

    @lazyproperty
    def min_y_coordinate_midplane(self):

        """
        This function ...
        :return:
        """

        return self.sorted_unique_y_coordinates_midplane[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def max_y_coordinate_midplane(self):

        """
        This function ...
        :return:
        """

        return self.sorted_unique_y_coordinates_midplane[-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_midplane_spacings(self):

        """
        This function ...
        :return:
        """

        return np.diff(self.sorted_unique_y_coordinates_midplane)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_midplane_average_spacing(self):

        """
        This function ...
        :return:
        """

        return np.mean(self.unique_y_coordinates_midplane_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_midplane_median_spacing(self):

        """
        This function ...
        :return:
        """

        return np.median(self.unique_y_coordinates_midplane_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_midplane_stddev_spacing(self):

        """
        This function ...
        :return:
        """

        return np.std(self.unique_y_coordinates_midplane_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_midplane_min_spacing(self):

        """
        Thisf unction ...
        :return:
        """

        return np.min(self.unique_y_coordinates_midplane_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_midplane_max_spacing(self):

        """
        This function ...
        :return:
        """

        return np.max(self.unique_y_coordinates_midplane_spacings)

    # -----------------------------------------------------------------

    def get_map_midplane(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_map_midplane: self.load_map_midplane()

        # Create
        else: self.create_map_midplane()

    # -----------------------------------------------------------------

    def load_map_midplane(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of the heating fraction in the midplane ...")

        # Load
        self.map_midplane = Frame.from_file(self.map_midplane_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def midplane_map_spacing(self):

        """
        This function ...
        :return:
        """

        if self.config.midplane_spacing_measure == "min": spacing = np.mean([self.unique_x_coordinates_midplane_min_spacing, self.unique_y_coordinates_midplane_min_spacing])
        elif self.config.midplane_spacing_measure == "max": spacing = np.mean([self.unique_x_coordinates_midplane_max_spacing, self.unique_y_coordinates_midplane_max_spacing])
        elif self.config.midplane_spacing_measure == "mean": spacing = np.mean([self.unique_x_coordinates_midplane_average_spacing, self.unique_y_coordinates_midplane_average_spacing])
        elif self.config.midplane_spacing_measure == "median": spacing = np.mean([self.unique_x_coordinates_midplane_median_spacing, self.unique_y_coordinates_midplane_median_spacing])
        else: raise ValueError("Invalid value for 'midplane_spacing_measure'")

        # Return
        return spacing * self.config.midplane_spacing_factor

    # -----------------------------------------------------------------

    @lazyproperty
    def midplane_map_radius(self):

        """
        This function ...
        :return:
        """

        return max(abs(self.min_x_coordinate_midplane), self.max_x_coordinate_midplane, abs(self.min_y_coordinate_midplane), self.max_y_coordinate_midplane)

    # -----------------------------------------------------------------

    @lazyproperty
    def midplane_map_coordinates(self):

        """
        This function ...
        :return:
        """

        return np.arange(-self.midplane_map_radius, self.midplane_map_radius, self.midplane_map_spacing)

    # -----------------------------------------------------------------

    @lazyproperty
    def half_midplane_map_spacing(self):

        """
        This function ...
        :return:
        """

        return 0.5 * self.midplane_map_spacing

    # -----------------------------------------------------------------

    def get_coordinate_range_around_pixel_midplane_map(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        x_range = RealRange(x-self.half_midplane_map_spacing,x+self.half_midplane_map_spacing)
        y_range = RealRange(y-self.half_midplane_map_spacing,y+self.half_midplane_map_spacing)
        return x_range, y_range

    # -----------------------------------------------------------------

    @lazyproperty
    def nmidplane_map_coordinates(self):

        """
        This function ...
        :return:
        """

        return len(self.midplane_map_coordinates)

    # -----------------------------------------------------------------

    @lazyproperty
    def midplane_map_shape(self):

        """
        This function ...
        :return:
        """

        return PixelShape.square(self.nmidplane_map_coordinates)

    # -----------------------------------------------------------------

    @property
    def midplane_map_npixels(self):

        """
        This function ...
        :return:
        """

        return self.midplane_map_shape.nxy

    # -----------------------------------------------------------------

    def get_pixel_in_midplane_map(self, x, y):

        """
        This function returns i, j
        :param x:
        :param y:
        :return:
        """

        i = nr.find_nearest_index(self.midplane_map_coordinates, x)
        j = nr.find_nearest_index(self.midplane_map_coordinates, y)
        return i, j

    # -----------------------------------------------------------------

    def get_coordinate_mask_x_for_map_midplane(self, x_range):

        """
        This function .....
        :param x_range:
        :return:
        """

        # x_range_pc = x_range.to("pc").value
        # x_mask = x_range.min < self.x_coordinates_midplane <= x_range.max
        x_mask = (x_range.min < self.x_coordinates_midplane) * (self.x_coordinates_midplane <= x_range.max)
        return x_mask

    # -----------------------------------------------------------------

    def get_coordinate_mask_y_for_map_midplane(self, y_range):

        """
        This function ...
        :param y_range:
        :return:
        """

        # y_range_pc = y_range.to("pc").value
        # y_mask = y_range.min < self.y_coordinates_midplane <= y_range.max
        y_mask = (y_range.min < self.y_coordinates_midplane) * (self.y_coordinates_midplane <= y_range.max)
        return y_mask

    # -----------------------------------------------------------------

    def get_coordinate_mask_for_map_midplane(self, x_range, y_range):

        """
        This function ...
        :param x_range:
        :param y_range:
        :return:
        """

        x_mask = self.get_coordinate_mask_x_for_map_midplane(x_range)
        y_mask = self.get_coordinate_mask_y_for_map_midplane(y_range)
        #print(x_mask, x_mask.shape)
        #print(y_mask, y_mask.shape)
        return x_mask * y_mask

    # -----------------------------------------------------------------

    def get_coordinate_indices_in_column_for_map_midplane(self, x_range, y_range):

        """
        This function ...
        :param x_range:
        :param y_range:
        :return:
        """

        mask = self.get_coordinate_mask_for_map_midplane(x_range, y_range)
        #print(mask.shape, np.sum(mask))

        #return np.where(x_mask * y_mask)
        return np.where(mask)[0]

    # -----------------------------------------------------------------

    def create_map_midplane(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of the heating fraction in the midplane ...")

        # INitialize map of the midplane heating fraction
        self.map_midplane = Frame.initialize_nans(self.midplane_map_shape)
        self.map_midplane_stddev = Frame.initialize_nans(self.midplane_map_shape)
        self.map_midplane_ncells = Frame.initialize_nans(self.midplane_map_shape)

        # Set the pixelscale and the coordinate info
        self.map_midplane.pixelscale = self.midplane_map_spacing * self.coordinates_unit
        self.map_midplane.distance = self.galaxy_distance
        self.map_midplane.set_meta("x_min", repr(-self.midplane_map_radius) + " " + self.coordinates_unit_string)
        self.map_midplane.set_meta("x_max", repr(self.midplane_map_radius) + " " + self.coordinates_unit_string)
        self.map_midplane.set_meta("y_min", repr(-self.midplane_map_radius) + " " + self.coordinates_unit_string)
        self.map_midplane.set_meta("y_max", repr(self.midplane_map_radius) + " " + self.coordinates_unit_string)

        # Loop over the x and y coordinates
        index = 1
        for i, x in enumerate(self.midplane_map_coordinates):
            for j, y in enumerate(self.midplane_map_coordinates):

                # Determine range of x and y
                x_range, y_range = self.get_coordinate_range_around_pixel_midplane_map(x, y)

                # Show
                if index % 100 == 0: log.debug("Calculating heating fraction of the midplane in the pixel " + str(index) + " of " + str(self.midplane_map_npixels) + " (" + tostr(float(index) / self.midplane_map_npixels * 100, decimal_places=1, round=True) + "%) ...")

                # Get the indices
                indices = self.get_coordinate_indices_in_column_for_map_midplane(x_range, y_range)
                nindices = indices.shape[0]

                # Set number of cells
                self.map_midplane_ncells[j,i] = nindices

                # If any cells
                if nindices > 0:

                    # Calculate the heating fraction
                    fractions = self.heating_fractions_midplane[indices]
                    weights = self.cell_weights_midplane[indices]

                    # Calculate the mean heatin fraction
                    fraction = numbers.weighed_arithmetic_mean_numpy(fractions, weights)
                    fraction_stddev = numbers.weighed_standard_deviation_numpy(fractions, weights, mean=fraction)

                    # Set fraction
                    self.map_midplane[j, i] = fraction
                    self.map_midplane_stddev[j, i] = fraction_stddev

                # Increment
                index += 1

    # -----------------------------------------------------------------

    @property
    def do_write_absorptions(self):

        """
        This function ...
        :return:
        """

        return not self.has_absorptions and self.absorptions is not None

    # -----------------------------------------------------------------

    @property
    def do_write_heating_fractions(self):

        """
        This function ...
        :return:
        """

        return not self.has_heating_fractions and self.heating_fractions is not None

    # -----------------------------------------------------------------

    @property
    def do_write_heating_fractions_diffuse(self):

        """
        This function ...
        :return:
        """

        return not self.has_heating_fractions_diffuse and self.diffuse_heating_fractions is not None

    # -----------------------------------------------------------------

    @property
    def do_write_distribution(self):

        """
        This function ...
        :return:
        """

        return not self.has_distribution and self.distribution is not None

    # -----------------------------------------------------------------

    @property
    def do_write_distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        return not self.has_distribution_diffuse and self.distribution_diffuse is not None

    # -----------------------------------------------------------------

    @property
    def do_write_radial_distribution(self):

        """
        This function ...
        :return:
        """

        return not self.has_radial_distribution and self.radial_distribution is not None

    # -----------------------------------------------------------------

    @property
    def do_write_map(self):

        """
        This function ...
        :return:
        """

        return not self.has_map and self.map is not None

    # -----------------------------------------------------------------

    @property
    def do_write_map_stddev(self):

        """
        This function ...
        :return:
        """

        return not self.has_map_stddev and self.map_stddev is not None

    # -----------------------------------------------------------------

    @property
    def do_write_map_ncells(self):

        """
        This function ...
        :return:
        """

        return not self.has_map_ncells and self.map_ncells is not None

    # -----------------------------------------------------------------

    @property
    def do_write_map_midplane(self):

        """
        This function ...
        :return:
        """

        return not self.has_map_midplane and self.map_midplane is not None

    # -----------------------------------------------------------------

    @property
    def do_write_map_midplane_stddev(self):

        """
        This function ...
        :return:
        """

        return not self.has_map_midplane_stddev and self.map_midplane_stddev is not None

    # -----------------------------------------------------------------

    @property
    def do_write_map_midplane_ncells(self):

        """
        This function ...
        :return:
        """

        return not self.has_map_midplane_ncells and self.map_midplane_ncells is not None

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

        # Write the map
        if self.do_write_map: self.write_map()

        # Write the map of the stddev
        if self.do_write_map_stddev: self.write_map_stddev()

        # Write the map of the number of cells of the map
        if self.do_write_map_ncells: self.write_map_ncells()

        # Write the map of the midplane
        if self.do_write_map_midplane: self.write_map_midplane()

        # Write the map of the stddev of the heating fraction in the midplane
        if self.do_write_map_midplane_stddev: self.write_map_midplane_stddev()

        # Write the map of the number of cells of the midplane map
        if self.do_write_map_midplane_ncells: self.write_map_midplane_ncells()

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
    def map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "map.fits")

    # -----------------------------------------------------------------

    @property
    def has_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.map_path)

    # -----------------------------------------------------------------

    def remove_map(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_path)

    # -----------------------------------------------------------------

    def write_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the heating fraction ...")

        # Save the map
        self.map.saveto(self.map_path)

    # -----------------------------------------------------------------

    @property
    def map_stddev_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "map_stddev.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_stddev(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.map_stddev_path)

    # -----------------------------------------------------------------

    def remove_map_stddev(self):

        """
        Thisn function ...
        :return:
        """

        fs.remove_file(self.map_stddev_path)

    # -----------------------------------------------------------------

    def write_map_stddev(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the standard deviation on the heating fraction ...")

        # Save the map
        self.map_stddev.saveto(self.map_stddev_path)

    # -----------------------------------------------------------------

    @property
    def map_ncells_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "map_ncells.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_ncells(self):

        """
        Thisn function ...
        :return:
        """

        return fs.is_file(self.map_ncells_path)

    # -----------------------------------------------------------------

    def remove_map_ncells(self):

        """
        Thisn function ...
        :return:
        """

        fs.remove_file(self.map_ncells_path)

    # -----------------------------------------------------------------

    def write_map_ncells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the number of cells used for the heating fraction map ...")

        # Save the map
        self.map_ncells.saveto(self.map_ncells_path)

    # -----------------------------------------------------------------

    @property
    def map_midplane_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "map_midplane.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_midplane(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.map_midplane_path)

    # -----------------------------------------------------------------

    def remove_map_midplane(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_midplane_path)

    # -----------------------------------------------------------------

    def write_map_midplane(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the heating fraction in the midplane ...")

        # Save the map
        self.map_midplane.saveto(self.map_midplane_path)

    # -----------------------------------------------------------------

    @property
    def map_midplane_stddev_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "map_midplane_stddev.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_midplane_stddev(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.map_midplane_stddev_path)

    # -----------------------------------------------------------------

    def remove_map_midplane_stddev(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_midplane_stddev_path)

    # -----------------------------------------------------------------

    def write_map_midplane_stddev(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the standard deviation on the heating fraction map in the midplane ...")

        # Save the map
        self.map_midplane_stddev.saveto(self.map_midplane_stddev_path)

    # -----------------------------------------------------------------

    @property
    def map_midplane_ncells_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "map_midplane_ncells.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_midplane_ncells(self):

        """
        THis function ...
        :return:
        """

        return fs.is_file(self.map_midplane_ncells_path)

    # -----------------------------------------------------------------

    def remove_map_midplane_ncells(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_midplane_ncells_path)

    # -----------------------------------------------------------------

    def write_map_midplane_ncells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the number of cells used for the heating fraction map in the midplane ...")

        # Save the map
        self.map_midplane_ncells.saveto(self.map_midplane_ncells_path)

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

        return self.config.plot_map and not self.has_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_map_midplane(self):

        """
        This function ...
        :return:
        """

        return self.config.plot_map_midplane and not self.has_map_midplane_plot

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
        #if self.do_plot_map: self.plot_map()

        # Plot
        #if self.do_plot_map_midplane: self.plot_map_midplane()

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
    def map_plot_path(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the plot file
        #return fs.join(self.cell_heating_path, "map.pdf")
        return fs.join(self.cell_heating_path, "map.png")

    # -----------------------------------------------------------------

    @property
    def has_map_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.map_plot_path)

    # -----------------------------------------------------------------

    def remove_map_plot(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_plot_path)

    # -----------------------------------------------------------------

    def plot_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a map of the heating fraction of the unevolved stellar population for a face-on view of the galaxy ...")

    # -----------------------------------------------------------------

    @property
    def map_midplane_plot_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "map_midplane.png")

    # -----------------------------------------------------------------

    @property
    def has_map_midplane_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.map_midplane_plot_path)

    # -----------------------------------------------------------------

    def remove_map_midplane_plot(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_midplane_plot_path)

    # -----------------------------------------------------------------

    def plot_map_midplane(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a map of the heating fraction in the midplane of the galaxy ...")

    # -----------------------------------------------------------------

    def remove_map_stddev_and_ncells(self):

        """
        This function ...
        :return:
        """

        if self.has_map: self.remove_map()
        if self.has_map_stddev: self.remove_map_stddev()
        if self.has_map_ncells: self.remove_map_ncells()

    # -----------------------------------------------------------------

    def remove_map_midplane_stddev_and_ncells(self):

        """
        This function ...
        :return:
        """

        if self.has_map_midplane: self.remove_map_midplane()
        if self.has_map_midplane_stddev: self.remove_map_midplane_stddev()
        if self.has_map_midplane_ncells: self.remove_map_midplane_ncells()

# -----------------------------------------------------------------
