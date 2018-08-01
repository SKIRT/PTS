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

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.tools.utils import lazyproperty, lazyfileproperty
from ...projection.data import Data3D, project_data
from ....core.basics.distribution import Distribution, Distribution2D

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

        # The heating map and related maps
        #self.map = None
        #self.map_stddev = None
        #self.map_ncells = None
        #self.map_interpolated = None

        # The heating map in the midplane and related maps
        #self.map_midplane = None
        #self.map_midplane_stddev = None
        #self.map_midplane_ncells = None
        #self.map_midplane_interpolated = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

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

    # -----------------------------------------------------------------

    @property
    def length_unit(self):
        return "pc"

    # -----------------------------------------------------------------

    @property
    def luminosity_unit(self):
        return "W"

    # -----------------------------------------------------------------
    # TOTAL
    # -----------------------------------------------------------------

    @property
    def total_absorptions_name(self):
        return "Labs_total"

    # -----------------------------------------------------------------

    @property
    def total_absorptions_description(self):
        return "Cell absorbed luminosities in the total simulation"

    # -----------------------------------------------------------------

    @property
    def total_absorptions_path(self):
        return fs.join(self.cell_heating_path, "absorptions_total.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_absorptions(self):
        return fs.is_file(self.total_absorptions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "total_absorptions_path", True, write=False)
    def total_absorptions(self):

        """
        This function ...
        :return:
        """

        # Calculate luminosities
        conversion_factor = self.total_contribution_absorption_unit.conversion_factor(self.luminosity_unit)
        luminosities_watt = self.total_contribution_absorption_luminosities * conversion_factor

        # Create the data
        return Data3D(self.total_absorptions_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                       luminosities_watt, length_unit=self.length_unit, unit=self.luminosity_unit,
                       description=self.total_absorptions_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # OLD
    # -----------------------------------------------------------------

    @property
    def old_absorptions_name(self):
        return "Labs_old"

    # -----------------------------------------------------------------

    @property
    def old_absorptions_description(self):
        return "Cell absorbed luminosities of the old stars"

    # -----------------------------------------------------------------

    @property
    def old_absorptions_path(self):
        return fs.join(self.cell_heating_path, "absorptions_old.dat")

    # -----------------------------------------------------------------

    @property
    def has_old_absorptions(self):
        return fs.is_file(self.old_absorptions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "old_absorptions_path", True, write=False)
    def old_absorptions(self):

        """
        This function ...
        :return:
        """

        # Calculate luminosities
        conversion_factor = self.old_contribution_absorption_unit.conversion_factor(self.luminosity_unit)
        luminosities_watt = self.old_contribution_absorption_luminosities * conversion_factor

        # Create the data
        return Data3D(self.old_absorptions_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                      luminosities_watt, length_unit=self.length_unit, unit=self.luminosity_unit,
                      description=self.old_absorptions_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # YOUNG
    # -----------------------------------------------------------------

    @property
    def young_absorptions_name(self):
        return "Labs_young"

    # -----------------------------------------------------------------

    @property
    def young_absorptions_description(self):
        return "Cell absorbed luminosities of the young stars"

    # -----------------------------------------------------------------

    @property
    def young_absorptions_path(self):
        return fs.join(self.cell_heating_path, "absorptions_young.dat")

    # -----------------------------------------------------------------

    @property
    def has_young_absorptions(self):
        return fs.is_file(self.young_absorptions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "young_absorptions_path", True, write=False)
    def young_absorptions(self):

        """
        This function ...
        :return:
        """

        # Calculate luminosities
        conversion_factor = self.young_contribution_absorption_unit.conversion_factor(self.luminosity_unit)
        luminosities_watt = self.young_contribution_absorption_luminosities * conversion_factor

        # Create the data
        return Data3D(self.young_absorptions_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                      luminosities_watt, length_unit=self.length_unit, unit=self.luminosity_unit,
                      description=self.young_absorptions_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # IONIZING
    # -----------------------------------------------------------------

    @property
    def ionizing_absorptions_name(self):
        return "Labs_ionizing"

    # -----------------------------------------------------------------

    @property
    def ionizing_absorptions_description(self):
        return "Cell absorbed luminosities of the ionizing stars"

    # -----------------------------------------------------------------

    @property
    def ionizing_absorptions_path(self):
        return fs.join(self.cell_heating_path, "absorptions_ionizing.dat")

    # -----------------------------------------------------------------

    @property
    def has_ionizing_absorptions(self):
        return fs.is_file(self.ionizing_absorptions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "ionizing_absorptions_path", True, write=False)
    def ionizing_absorptions(self):

        """
        This function ...
        :return:
        """

        # Calculate luminosities
        conversion_factor = self.ionizing_contribution_absorption_unit.conversion_factor(self.luminosity_unit)
        luminosities_watt = self.ionizing_contribution_absorption_luminosities * conversion_factor

        # Create the data
        return Data3D(self.ionizing_absorptions_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                      luminosities_watt, length_unit=self.length_unit, unit=self.luminosity_unit,
                      description=self.ionizing_absorptions_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # INTERNAL
    # -----------------------------------------------------------------

    @property
    def internal_absorptions_name(self):
        return "Labs_internal"

    # -----------------------------------------------------------------

    @property
    def internal_absorptions_description(self):
        return "Cell absorbed luminosities in star formation regions dust"

    # -----------------------------------------------------------------

    @property
    def internal_absorptions_path(self):
        return fs.join(self.cell_heating_path, "absorptions_internal.dat")

    # -----------------------------------------------------------------

    @property
    def has_internal_absorptions(self):
        return fs.is_file(self.internal_absorptions_path)

    # -----------------------------------------------------------------

    @property
    def volumes(self):
        return self.model.cell_volumes

    # -----------------------------------------------------------------

    @property
    def ionizing_densities(self):
        return self.model.sfr_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def internal_absorbed_luminosity(self):
        return self.model.intrinsic_dust_luminosity_sfr.to(self.luminosity_unit, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def relative_masses_ionizing(self):
        masses = self.volumes * self.ionizing_densities
        masses /= sum(masses) # normalized
        return masses

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "internal_absorptions_path", True, write=False)
    def internal_absorptions(self):

        """
        This function ...
        :return:
        """

        # absorptions_ionizing_internal_watt = volumes * density * absorbed_energy
        luminosities_watt = self.relative_masses_ionizing * self.internal_absorbed_luminosity

        # Create the data
        return Data3D(self.internal_absorptions_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                      luminosities_watt, length_unit=self.length_unit, unit=self.luminosity_unit,
                      description=self.internal_absorptions_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # DIFFUSE HEATING FRACTIONS
    # -----------------------------------------------------------------

    @property
    def diffuse_fractions_name(self):
        return "Funev_diffuse"

    # -----------------------------------------------------------------

    @property
    def diffuse_fractions_description(self):
        return "Heating fraction by the unevolved stellar populations (diffuse dust)"

    # -----------------------------------------------------------------

    @property
    def young_absorption_values(self):
        return self.young_absorptions.values

    # -----------------------------------------------------------------

    @property
    def ionizing_absorption_values(self):
        return self.ionizing_absorptions.values

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_values_diffuse(self):
        return self.young_absorption_values + self.ionizing_absorption_values

    # -----------------------------------------------------------------

    @property
    def total_absorption_values_diffuse(self):
        return self.total_absorptions.values

    # -----------------------------------------------------------------

    @lazyproperty
    def fraction_values_diffuse(self):
        return self.unevolved_absorption_values_diffuse / self.total_absorption_values_diffuse

    # -----------------------------------------------------------------

    @property
    def diffuse_fractions_path(self):
        return fs.join(self.cell_heating_path, "diffuse_fractions.dat")

    # -----------------------------------------------------------------

    @property
    def has_diffuse_fractions(self):
        return fs.is_file(self.diffuse_fractions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "diffuse_fractions_path", True, write=False)
    def diffuse_fractions(self):

        """
        This function ...
        :return:
        """

        # Create and return the data
        return Data3D(self.diffuse_fractions_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                      self.fraction_values_diffuse, length_unit=self.length_unit,
                      description=self.diffuse_fractions_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # TOTAL HEATING FRACTIONS
    # -----------------------------------------------------------------

    @property
    def total_fractions_name(self):
        return "Funev_total"

    # -----------------------------------------------------------------

    @property
    def total_fractions_description(self):
        return "Heating fraction by the unevolved stellar population (all dust)"

    # -----------------------------------------------------------------

    @property
    def internal_absorption_values(self):
        return self.internal_absorptions.values

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_values(self):
        return self.unevolved_absorption_values_diffuse + self.internal_absorption_values

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_values(self):
        return self.total_absorption_values_diffuse + self.internal_absorption_values

    # -----------------------------------------------------------------

    @lazyproperty
    def fraction_values(self):
        return self.unevolved_absorption_values / self.total_absorption_values

    # -----------------------------------------------------------------

    @property
    def total_fractions_path(self):
        return fs.join(self.cell_heating_path, "total_fractions.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_fractions(self):
        return fs.is_file(self.total_fractions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "total_fractions_path", True, write=False)
    def total_fractions(self):

        """
        This function ...
        :return:
        """

        # Create and return the data
        return Data3D(self.total_fractions_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                      self.fraction_values, length_unit=self.length_unit, description=self.total_fractions_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # DISTRIBUTION DIFFUSE
    # -----------------------------------------------------------------

    @property
    def distribution_diffuse_name(self):
        return "Heating fraction (diffuse dust)"

    # -----------------------------------------------------------------

    @property
    def distribution_diffuse_path(self):
        return fs.join(self.cell_heating_path, "distribution_diffuse.dat")

    # -----------------------------------------------------------------

    @property
    def has_distribution_diffuse(self):
        return fs.is_file(self.distribution_diffuse_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Distribution, "distribution_diffuse_path", True, write=False)
    def distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the distribution of heating fractions of the unevolved stellar population (diffuse) ...")

        # Generate the distribution
        # Weights are dust mass fraction
        return Distribution.from_values(self.distribution_diffuse_name, self.valid_heating_fractions_diffuse, nbins=self.config.nbins, weights=self.valid_cell_weights_diffuse)

    # -----------------------------------------------------------------
    # DISTRIBUTION TOTAL
    # -----------------------------------------------------------------

    @property
    def distribution_total_name(self):
        return "Heating fraction (all dust)"

    # -----------------------------------------------------------------

    @property
    def distribution_total_path(self):
        return fs.join(self.cell_heating_path, "distribution_total.dat")

    # -----------------------------------------------------------------

    @property
    def has_distribution_total(self):
        return fs.is_file(self.distribution_total_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Distribution, "distribution_total_path", True, write=False)
    def distribution_total(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the distribution of heating fractions of the unevolved stellar population ...")

        # Generate the distribution
        # Weights are dust mass fraction
        return Distribution.from_values(self.distribution_total_name, self.valid_heating_fractions, nbins=self.config.nbins, weights=self.valid_cell_weights)

    # -----------------------------------------------------------------
    # RADIAL DISTRIBUTION
    # -----------------------------------------------------------------

    @property
    def radial_distribution_name(self):
        return "Heating fraction"

    # -----------------------------------------------------------------

    @property
    def radial_distribution_path(self):
        return fs.join(self.cell_heating_path, "radial_distribution.dat")

    # -----------------------------------------------------------------

    @property
    def has_radial_distribution(self):
        return fs.is_file(self.radial_distribution_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Distribution2D, "radial_distribution_path", True, write=False)
    def radial_distribution(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Calculating the radial distribution of heating fractions of the unevolved stellar population ...")

        # Generate the radial distribution
        return Distribution2D.from_values(self.valid_radii, self.valid_heating_fractions,
                                          weights=self.valid_cell_weights, x_name="radius (pc)",
                                          y_name=self.radial_distribution_name,
                                          nbins=self.config.nradial_bins)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write absorptions
        self.write_absorptions()

        # Write the heating fractions
        self.write_fractions()

        # Write the distributions
        self.write_distributions()

    # -----------------------------------------------------------------

    def write_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the absorption data ...")

        # Total
        if self.do_write_total_absorptions: self.write_total_absorptions()

        # Old
        if self.do_write_old_absorptions: self.write_old_absorptions()

        # Young
        if self.do_write_young_absorptions: self.write_young_absorptions()

        # Ionizing
        if self.do_write_ionizing_absorptions: self.write_ionizing_absorptions()

        # Internal
        if self.do_write_internal_absorptions: self.write_internal_absorptions()

    # -----------------------------------------------------------------

    @property
    def do_write_total_absorptions(self):
        return not self.has_total_absorptions

    # -----------------------------------------------------------------

    def write_total_absorptions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.total_absorptions.saveto(self.total_absorptions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_old_absorptions(self):
        return not self.has_old_absorptions

    # -----------------------------------------------------------------

    def write_old_absorptions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.old_absorptions.saveto(self.old_absorptions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_young_absorptions(self):
        return not self.has_young_absorptions

    # -----------------------------------------------------------------

    def write_young_absorptions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.young_absorptions.saveto(self.young_absorptions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_ionizing_absorptions(self):
        return not self.has_ionizing_absorptions

    # -----------------------------------------------------------------

    def write_ionizing_absorptions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.ionizing_absorptions.saveto(self.ionizing_absorptions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_internal_absorptions(self):
        return not self.has_internal_absorptions

    # -----------------------------------------------------------------

    def write_internal_absorptions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.internal_absorptions.saveto(self.internal_absorptions_path)

    # -----------------------------------------------------------------

    def write_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the heating fraction data ...")

        # Diffuse
        if self.do_write_diffuse_fractions: self.write_diffuse_fractions()

        # Total
        if self.do_write_total_fractions: self.write_total_fractions()

    # -----------------------------------------------------------------

    @property
    def do_write_diffuse_fractions(self):
        return not self.has_diffuse_fractions

    # -----------------------------------------------------------------

    def write_diffuse_fractions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.diffuse_fractions.saveto(self.diffuse_fractions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_total_fractions(self):
        return not self.has_total_fractions

    # -----------------------------------------------------------------

    def write_total_fractions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.total_fractions.saveto(self.total_fractions_path)

    # -----------------------------------------------------------------

    def write_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the distributions of the heating fraction ...")

        # Diffuse
        if self.do_write_distribution_diffuse: self.write_distribution_diffuse()

        # Total
        if self.do_write_distribution_total: self.write_distribution_total()

        # Radial
        if self.do_write_radial_distribution: self.write_radial_distribution()

    # -----------------------------------------------------------------

    @property
    def do_write_distribution_diffuse(self):
        return not self.has_distribution_diffuse

    # -----------------------------------------------------------------

    def write_distribution_diffuse(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @property
    def do_write_distribution_total(self):
        return not self.has_distribution_total

    # -----------------------------------------------------------------

    def write_distribution_total(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @property
    def do_write_radial_distribution(self):
        return not self.has_radial_distribution

    # -----------------------------------------------------------------

    def write_radial_distribution(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
