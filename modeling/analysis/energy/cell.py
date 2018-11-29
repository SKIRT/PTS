#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.energy.cell Contains the CellEnergyAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ..component import AnalysisRunComponent
from ....core.basics.configuration import open_box
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.tools.utils import lazyproperty, lazyfileproperty
from ...core.data import Data3D
from ....core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

class CellEnergyAnalyser(AnalysisRunComponent):
    
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
        super(CellEnergyAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Show
        self.show()

        # Write
        self.write()

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_energy_path(self):
        return fs.create_directory_in(self.analysis_run.energy_path, "cell")

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_luminosity_unit(self):
        return u("W")

    # -----------------------------------------------------------------

    @property
    def length_unit(self):
        return "pc"

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def total_diffuse_luminosities(self):
        return self.total_contribution_absorption_luminosities

    # -----------------------------------------------------------------

    @property
    def total_diffuse_luminosities_unit(self):
        return self.total_contribution_absorption_unit

    # -----------------------------------------------------------------

    @property
    def total_diffuse_luminosity(self):
        return np.sum(self.total_diffuse_luminosities) * self.total_diffuse_luminosities_unit

    # -----------------------------------------------------------------

    @property
    def sfr_stellar_densities(self):
        return self.sfr_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def volumes(self):
        return self.model.cell_volumes

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_relative_masses(self):
        masses = self.volumes * self.sfr_stellar_densities
        masses /= sum(masses) # normalized
        return masses

    # -----------------------------------------------------------------

    @lazyproperty
    def simple_internal_absorbed_luminosity(self):
        return self.model.intrinsic_dust_luminosity_sfr.to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def absorption_path(self):
        return self.analysis_run.absorption_path

    # -----------------------------------------------------------------
    # TOTAL
    # -----------------------------------------------------------------

    @property
    def total_absorption_properties_path(self):
        return fs.join(self.absorption_path, "total.txt")

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_properties(self):
        return open_box(self.total_absorption_properties_path)

    # -----------------------------------------------------------------

    @property
    def total_absorbed_luminosity_diffuse(self):
        return self.total_absorption_properties.diffuse.absorbed

    # -----------------------------------------------------------------

    @property
    def total_dust_luminosity_diffuse(self):
        return self.total_absorption_properties.diffuse.dust

    # -----------------------------------------------------------------

    @property
    def total_absorbed_luminosity_all(self):
        return self.total_absorption_properties.all.absorbed

    # -----------------------------------------------------------------

    @property
    def total_dust_luminosity_all(self):
        return self.total_absorption_properties.all.dust

    # -----------------------------------------------------------------
    # SFR
    # -----------------------------------------------------------------

    @property
    def sfr_absorption_properties_path(self):
        return fs.join(self.absorption_path, "sfr.txt")

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_properties(self):
        return open_box(self.sfr_absorption_properties_path)

    # -----------------------------------------------------------------

    @property
    def internal_absorbed_luminosity(self):
        return self.sfr_absorption_properties.internal.absorbed

    # -----------------------------------------------------------------

    @property
    def internal_dust_luminosity(self):
        return self.sfr_absorption_properties.internal.dust

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def internal_luminosities(self):
        return self.sfr_relative_masses * self.internal_absorbed_luminosity.to(self.bolometric_luminosity_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def total_luminosities_path(self):
        return fs.join(self.cell_energy_path, "total.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_luminosities(self):
        return fs.is_file(self.total_luminosities_path)

    # -----------------------------------------------------------------

    @property
    def dust_lum_name(self):
        return "Ldust"

    # -----------------------------------------------------------------

    @property
    def dust_lum_description(self):
        return "Bolometric dust luminosity"

    # -----------------------------------------------------------------

    @lazyproperty
    def total_luminosity_values(self):
        return self.total_diffuse_luminosities + self.internal_luminosities

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "total_luminosities_path", True, write=False)
    def total_luminosities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the total dust cell luminosities ...")

        # Create the data with external xyz
        return Data3D.from_values(self.dust_lum_name, self.total_luminosity_values, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.dust_lum_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.bolometric_luminosity_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_luminosity_density_values(self):
        return self.total_luminosity_values / self.cell_volumes

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_luminosity_density_unit(self):
        return self.bolometric_luminosity_unit / self.cell_volume_unit

    # -----------------------------------------------------------------

    @property
    def total_luminosity_densities_path(self):
        return fs.join(self.cell_energy_path, "total_density.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_luminosity_densities(self):
        return fs.is_file(self.total_luminosity_densities_path)

    # -----------------------------------------------------------------

    @property
    def dust_lum_density_name(self):
        return "vLdust"

    # -----------------------------------------------------------------

    @property
    def dust_lum_density_description(self):
        return "Bolometric dust luminosity density"

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "total_luminosity_densities_path", True, write=False)
    def total_luminosity_densities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the total dust cell luminosity densities ...")

        # Create the data with external xyz
        return Data3D.from_values(self.dust_lum_density_name, self.total_luminosity_density_values, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.dust_lum_density_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.bolometric_luminosity_density_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        print("")

        print("TOTAL:")
        print(self.total_diffuse_luminosity) # USE THIS?
        print(self.total_absorbed_luminosity_diffuse)
        print(self.total_dust_luminosity_diffuse)
        print(self.total_absorbed_luminosity_all)
        print(self.total_dust_luminosity_all)
        print("")

        print("SFR:")
        print("")
        print(self.simple_internal_absorbed_luminosity) # NAIVE, BUT IN BETWEEN THE OTHER TWO?
        print(self.internal_absorbed_luminosity) # USE THIS?
        print(self.internal_dust_luminosity)
        print("")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the dust luminosities in each cell
        self.write_luminosities()

    # -----------------------------------------------------------------

    @property
    def do_write_total_luminosities(self):
        return not self.has_total_luminosities

    # -----------------------------------------------------------------

    @property
    def do_write_total_luminosity_densities(self):
        return not self.has_total_luminosity_densities

    # -----------------------------------------------------------------

    def write_luminosities(self):

        """
        This function ...
        :return:
        """

        # Total
        if self.do_write_total_luminosities: self.write_total_luminosities()

        # Total luminosity density
        if self.do_write_total_luminosity_densities: self.write_total_luminosity_densities()

    # -----------------------------------------------------------------

    def write_total_luminosities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the total dust cell luminosities ...")

        # Write
        self.total_luminosities.saveto(self.total_luminosities_path)

    # -----------------------------------------------------------------

    def write_total_luminosity_densities(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the total dust cell luminosity densities ...")

        # Write
        self.total_luminosity_densities.saveto(self.total_luminosity_densities_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
