#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.sfr Contains the SFRAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import AnalysisComponent, AnalysisRunComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.core.frame import Frame
from ...core.tools.utils import lazyproperty
from ..core.data import Data3D
from ..projection.data import project_data
from ..core.model import oliver_stellar_mass, salim_fuv_to_sfr
from ...core.units.parsing import parse_unit as u
from ...magic.tools.plotting import plot_map

# -----------------------------------------------------------------

projected_name = "projected"
cell_name = "cell"

# -----------------------------------------------------------------

class SFRAnalyser(AnalysisRunComponent):
    
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
        super(SFRAnalyser, self).__init__(*args, **kwargs)

        # The projected maps
        self.sfr_earth_map = None
        self.sfr_faceon_map = None
        self.stellar_mass_earth_map = None
        self.stellar_mass_faceon_map = None
        self.ssfr_earth_map = None
        self.ssfr_faceon_map = None

        # The 3D data
        self.fuv_data = None # intrinsic
        self.i1_data = None  # also intrinsic
        self.sfr_data = None
        self.stellar_mass_data = None
        self.ssfr_data = None

        # The cell maps
        self.sfr_data_faceon_map = None
        self.stellar_mass_data_faceon_map = None
        self.ssfr_data_faceon_map = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get projected star formation rates
        self.get_projected()

        # Get cell star formation rates
        self.get_cell()

        # Get the projected cell star formation rates
        self.get_cell_maps()

        # Writing
        self.write()

        # Plotting
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SFRAnalyser, self).setup()

    # -----------------------------------------------------------------

    @property
    def earth_projection(self):
        return self.analysis_run.earth_projection

    # -----------------------------------------------------------------

    @property
    def faceon_projection(self):
        return self.analysis_run.faceon_projection

    # -----------------------------------------------------------------

    @property
    def edgeon_projection(self):
        return self.analysis_run.edgeon_projection

    # -----------------------------------------------------------------

    def get_projected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the projected star formation rates ...")

        # Star formation rate
        self.get_projected_sfr()

        # Stellar mass
        self.get_projected_mass()

        # Specific star formation rate
        self.get_projected_ssfr()

    # -----------------------------------------------------------------

    @property
    def needs_projected_sfr_earth(self):
        return not self.has_projected_sfr_earth or self.do_plot_projected_sfr_earth_map

    # -----------------------------------------------------------------

    @property
    def needs_projected_sfr_faceon(self):
        return not self.has_projected_sfr_faceon or self.do_plot_projected_sfr_faceon_map

    # -----------------------------------------------------------------

    def get_projected_sfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the projected star formation rate ...")

        # Earth
        if self.needs_projected_sfr_earth: self.get_projected_sfr_earth()

        # Faceon
        if self.needs_projected_sfr_faceon: self.get_projected_sfr_faceon()

    # -----------------------------------------------------------------

    def get_projected_sfr_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the projected star formation rate from the earth projection ...")

        # Load
        if self.has_projected_sfr_earth: self.sfr_earth_map = Frame.from_file(self.projected_sfr_earth_path)

        # Calculate
        else: self.sfr_earth_map = self.model.total_star_formation_rate_map_earth

    # -----------------------------------------------------------------

    def get_projected_sfr_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the projected star formation rate from the faceon projection ...")

        # Load
        if self.has_projected_sfr_faceon: self.sfr_faceon_map = Frame.from_file(self.projected_sfr_faceon_path)

        # Calculate
        else: self.sfr_faceon_map = self.model.total_star_formation_rate_map_faceon

    # -----------------------------------------------------------------

    @property
    def needs_projected_mass_earth(self):
        return not self.has_projected_mass_earth or self.do_plot_projected_mass_earth

    # -----------------------------------------------------------------

    @property
    def needs_projected_mass_faceon(self):
        return not self.has_projected_mass_faceon or self.do_plot_projected_mass_faceon

    # -----------------------------------------------------------------

    def get_projected_mass(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the projected stellar mass ...")

        # Earth
        if self.needs_projected_mass_earth: self.get_projected_mass_earth()

        # Faceon
        if self.needs_projected_mass_faceon: self.get_projected_mass_faceon()

    # -----------------------------------------------------------------

    def get_projected_mass_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the projected stellar mass from the earth projection ...")

        # Load
        if self.has_projected_mass_earth: self.stellar_mass_earth_map = Frame.from_file(self.projected_mass_earth_path)

        # Calculate
        else: self.stellar_mass_earth_map = self.model.total_stellar_mass_map_earth

    # -----------------------------------------------------------------

    def get_projected_mass_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the projected stellar mass from the faceon projection ...")

        # Load
        if self.has_projected_mass_faceon: self.stellar_mass_faceon_map = Frame.from_file(self.projected_mass_faceon_path)

        # Calculate
        else: self.stellar_mass_faceon_map = self.model.total_stellar_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def needs_projected_ssfr_earth(self):
        return not self.has_projected_ssfr_earth or self.do_plot_projected_ssfr_earth

    # -----------------------------------------------------------------

    @property
    def needs_projected_ssfr_faceon(self):
        return not self.has_projected_ssfr_faceon or self.do_plot_projected_ssfr_faceon

    # -----------------------------------------------------------------

    def get_projected_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the projected specific star formation rate ...")

        # Earth
        if self.needs_projected_ssfr_earth: self.get_projected_ssfr_earth()

        # Faceon
        if self.needs_projected_ssfr_faceon: self.get_projected_ssfr_faceon()

    # -----------------------------------------------------------------

    def get_projected_ssfr_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the projected specific star formation rate from the earth projection ...")

        # Load
        if self.has_projected_ssfr_earth: self.ssfr_earth_map = Frame.from_file(self.projected_ssfr_earth_path)

        # Calculate
        else: self.ssfr_earth_map = self.model.total_ssfr_map_earth

    # -----------------------------------------------------------------

    def get_projected_ssfr_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_projected_ssfr_faceon: self.ssfr_faceon_map = Frame.from_file(self.projected_ssfr_faceon_path)

        # Calculate
        else: self.ssfr_faceon_map = self.model.total_ssfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def needs_cell_fuv(self):
        return not self.has_cell_fuv or not self.has_cell_sfr

    # -----------------------------------------------------------------

    @property
    def needs_cell_i1(self):
        return not self.has_cell_i1 or not self.has_cell_mass

    # -----------------------------------------------------------------

    @property
    def needs_cell_sfr(self):
        return not self.has_cell_sfr or (not self.has_cell_ssfr and self.needs_cell_ssfr) or (not self.has_cell_sfr_map and self.needs_cell_sfr_map)

    # -----------------------------------------------------------------

    @property
    def needs_cell_mass(self):
        return not self.has_cell_mass or (not self.has_cell_ssfr and self.needs_cell_ssfr) or (not self.has_cell_mass_map and self.needs_cell_mass_map)

    # -----------------------------------------------------------------

    @property
    def needs_cell_ssfr(self):
        return not self.has_cell_ssfr or (not self.has_cell_ssfr_map and self.needs_cell_ssfr_map)

    # -----------------------------------------------------------------

    def get_cell(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the cell star formation rates ...")

        # FUV luminosity
        if self.needs_cell_fuv: self.get_cell_fuv()

        # I1 luminosity
        if self.needs_cell_i1: self.get_cell_i1()

        # Star formation rate
        if self.needs_cell_sfr: self.get_cell_sfr()

        # Stellar mass
        if self.needs_cell_mass: self.get_cell_mass()

        # Specific star formation rate
        if self.needs_cell_ssfr: self.get_cell_ssfr()

    # -----------------------------------------------------------------

    @property
    def cell_x_coordinates(self):
        return self.model.cell_x_coordinates

    # -----------------------------------------------------------------

    @property
    def cell_y_coordinates(self):
        return self.model.cell_y_coordinates

    # -----------------------------------------------------------------

    @property
    def cell_z_coordinates(self):
        return self.model.cell_z_coordinates

    # -----------------------------------------------------------------

    @property
    def cell_volumes(self):
        return self.model.cell_volumes # is array

    # -----------------------------------------------------------------

    @property
    def young_cell_stellar_density(self):
        return self.model.young_cell_stellar_density # is array

    # -----------------------------------------------------------------

    @lazyproperty
    def young_cell_normalized_mass(self):
        values = self.young_cell_stellar_density * self.cell_volumes
        values /= np.sum(values)
        return values

    # -----------------------------------------------------------------

    @property
    def sfr_cell_stellar_density(self):
        return self.model.sfr_cell_stellar_density # is array

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_cell_normalized_mass(self):
        values = self.sfr_cell_stellar_density * self.cell_volumes
        values /= np.sum(values)
        return values

    # -----------------------------------------------------------------

    @property
    def bulge_cell_stellar_density(self):
        return self.model.old_bulge_cell_stellar_density # is array

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_cell_normalized_mass(self):
        values = self.bulge_cell_stellar_density * self.cell_volumes
        values /= np.sum(values)
        return values

    # -----------------------------------------------------------------

    @property
    def disk_cell_stellar_density(self):
        return self.model.old_disk_cell_stellar_density # is array

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_cell_normalized_mass(self):
        values = self.disk_cell_stellar_density * self.cell_volumes
        values /= np.sum(values)
        return values

    # -----------------------------------------------------------------

    @property
    def young_intrinsic_fuv_luminosity(self):
        return self.model.intrinsic_fuv_luminosity_young

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_fuv_luminosity_scalar(self):
        return self.young_intrinsic_fuv_luminosity.to(self.fuv_luminosity_unit).value

    # -----------------------------------------------------------------

    @property
    def sfr_intrinsic_fuv_luminosity(self):
        return self.model.intrinsic_fuv_luminosity_sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_fuv_luminosity_scalar(self):
        return self.sfr_intrinsic_fuv_luminosity.to(self.fuv_luminosity_unit).value

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
    def fuv_luminosity_unit(self):
        return u("W/micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def young_cell_fuv_luminosities(self):
        return self.young_cell_normalized_mass * self.young_intrinsic_fuv_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_cell_fuv_luminosities(self):
        return self.sfr_cell_normalized_mass * self.sfr_intrinsic_fuv_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_cell_fuv_luminosities(self):
        return self.young_cell_fuv_luminosities + self.sfr_cell_fuv_luminosities

    # -----------------------------------------------------------------

    @property
    def fuv_name(self):
        return "FUV"

    # -----------------------------------------------------------------

    @property
    def fuv_description(self):
        return "Intrinsic FUV luminosity of unevolved stars"

    # -----------------------------------------------------------------

    def get_cell_fuv(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the cell FUV luminosity ...")

        # Create the data
        self.fuv_data = Data3D(self.fuv_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates, self.unevolved_cell_fuv_luminosities, length_unit="pc", unit=self.fuv_luminosity_unit, description=self.fuv_description, distance=self.galaxy_distance, wavelength=self.fuv_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_luminosity_unit(self):
        return u("W/micron")

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

    def get_cell_i1(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the cell I1 luminosity ...")

        # Create the data
        self.i1_data = Data3D("I1", self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates, self.old_cell_i1_luminosities, length_unit="pc", unit=self.i1_luminosity_unit, description="Intrinsic I1 luminosity of evolved stars", distance=self.galaxy_distance, wavelength=self.i1_wavelength)

    # -----------------------------------------------------------------

    def get_cell_sfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the cell star formation rate ...")

        # Load
        if self.has_cell_sfr: self.sfr_data = Data3D.from_file(self.cell_sfr_path)

        # Calculate
        else: self.calculate_cell_sfr()

    # -----------------------------------------------------------------

    def calculate_cell_sfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell star formation rate ...")

        # Calculate the SFR data from the intrinsic FUV data
        self.sfr_data = salim_fuv_to_sfr(self.fuv_data)

    # -----------------------------------------------------------------

    def get_cell_mass(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the cell stellar mass ...")

        # Load
        if self.has_cell_mass: self.stellar_mass_data = Data3D.from_file(self.cell_mass_path)

        # Calculate
        else: self.calculate_cell_mass()

    # -----------------------------------------------------------------

    def calculate_cell_mass(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell stellar mass ...")

        # Calculate the stellar mass data from the I1 data
        self.stellar_mass_data = oliver_stellar_mass(self.i1_data, hubble_type=self.hubble_stage_type, hubble_subtype=self.hubble_stage_subtype)

    # -----------------------------------------------------------------

    def get_cell_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the cell specific star formation rate ...")

        # Load
        if self.has_cell_ssfr: self.ssfr_data = Data3D.from_file(self.cell_ssfr_path)

        # Calculate
        else: self.calculate_cell_ssfr()

    # -----------------------------------------------------------------

    @property
    def cell_sfrs(self):
        return self.sfr_data.values

    # -----------------------------------------------------------------

    @property
    def sfr_unit(self):
        return self.sfr_data.unit

    # -----------------------------------------------------------------

    @property
    def cell_masses(self):
        return self.stellar_mass_data.values

    # -----------------------------------------------------------------

    @property
    def stellar_mass_unit(self):
        return self.stellar_mass_data.unit

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfrs(self):
        return self.cell_sfrs / self.cell_masses

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_unit(self):
        return self.sfr_unit / self.stellar_mass_unit

    # -----------------------------------------------------------------

    def calculate_cell_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell specific star formation rate ...")

        # Create the data
        self.ssfr_data = Data3D("sSFR", self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                              self.cell_ssfrs, length_unit="pc", unit=self.ssfr_unit, description="specific star formation rate", distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def needs_cell_sfr_map(self):
        return not self.has_cell_sfr_map or self.do_plot_cell_sfr_map

    # -----------------------------------------------------------------

    @property
    def needs_cell_mass_map(self):
        return not self.has_cell_sfr_map or self.do_plot_cell_mass_map

    # -----------------------------------------------------------------

    @property
    def needs_cell_ssfr_map(self):
        return not self.has_cell_sfr_map or self.do_plot_cell_ssfr_map

    # -----------------------------------------------------------------

    def get_cell_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the maps of the cell star formation rates ...")

        # Star formation rate
        if self.needs_cell_sfr_map: self.get_cell_sfr_map()

        # Stellar mass
        if self.needs_cell_mass_map: self.get_cell_mass_map()

        # Specific star formation rate
        if self.needs_cell_ssfr_map: self.get_cell_ssfr_map()

    # -----------------------------------------------------------------

    def get_cell_sfr_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the map of the cell star formation rate ...")

        # Load?
        if self.has_cell_sfr_map: self.sfr_data_faceon_map = Frame.from_file(self.cell_sfr_map_path)

        # Create
        else: self.create_cell_sfr_map()

    # -----------------------------------------------------------------

    @property
    def sfr_name(self):
        return "SFR"

    # -----------------------------------------------------------------

    @property
    def sfr_description(self):
        return "Star formation rate"

    # -----------------------------------------------------------------

    def create_cell_sfr_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of the cell star formation rate ...")

        # Create
        self.sfr_data_faceon_map = project_data(self.sfr_name, self.sfr_data, self.faceon_projection, description=self.sfr_description)

    # -----------------------------------------------------------------

    def get_cell_mass_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the map of the cell stellar mass ...")

        # Load?
        if self.has_cell_mass_map: self.stellar_mass_data_faceon_map = Frame.from_file(self.cell_mass_map_path)

        # Create
        else: self.create_cell_mass_map()

    # -----------------------------------------------------------------

    @property
    def stellar_mass_name(self):
        return "Mstellar"

    # -----------------------------------------------------------------

    @property
    def stellar_mass_description(self):
        return "Stellar mass"

    # -----------------------------------------------------------------

    def create_cell_mass_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of the cell stellar mass ...")

    # -----------------------------------------------------------------

    def get_cell_ssfr_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the map of the cell specific star formation rate ...")

        # Load?
        if self.has_cell_ssfr_map: self.ssfr_data_faceon_map = Frame.from_file(self.cell_ssfr_map_path)

        # Create
        else: self.create_cell_ssfr_map()

    # -----------------------------------------------------------------

    @property
    def ssfr_name(self):
        return "sSFR"

    # -----------------------------------------------------------------

    @property
    def ssfr_description(self):
        return "Specific star formation rate"

    # -----------------------------------------------------------------

    def create_cell_ssfr_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of the cell specific star formation rate ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Projected
        self.write_projected()

        # Cell
        self.write_cell()

        # Cell maps
        self.write_cell_maps()

    # -----------------------------------------------------------------

    @property
    def sfr_path(self):
        return self.analysis_run.sfr_path

    # -----------------------------------------------------------------

    @lazyproperty
    def projected_path(self):
        return fs.create_directory_in(self.sfr_path, projected_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_path(self):
        return fs.create_directory_in(self.sfr_path, cell_name)

    # -----------------------------------------------------------------

    def write_projected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the projected star formation rates ...")

        # Star formation rate
        self.write_projected_sfr()

        # Stellar mass
        self.write_projected_mass()

        # Specific star formation rate
        self.write_projected_ssfr()

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_earth(self):
        return not self.has_projected_sfr_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_faceon(self):
        return not self.has_projected_sfr_faceon

    # -----------------------------------------------------------------

    def write_projected_sfr(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_projected_sfr_earth: self.write_projected_sfr_earth()

        # Faceon
        if self.do_write_projected_sfr_faceon: self.write_projected_sfr_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_sfr_earth_path(self):
        return fs.join(self.projected_path, "sfr_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_earth(self):
        return fs.is_file(self.projected_sfr_earth_path)

    # -----------------------------------------------------------------

    def write_projected_sfr_earth(self):

        """
        Thisn function ...
        :return:
        """

        # Write
        self.sfr_earth_map.saveto(self.projected_sfr_earth_path)

    # -----------------------------------------------------------------

    @property
    def projected_sfr_faceon_path(self):
        return fs.join(self.projected_path, "sfr_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_faceon(self):
        return fs.is_file(self.projected_sfr_faceon_path)

    # -----------------------------------------------------------------

    def write_projected_sfr_faceon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_faceon_map.saveto(self.projected_sfr_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_mass_earth(self):
        return not self.has_projected_mass_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_mass_faceon(self):
        return not self.has_projected_mass_faceon

    # -----------------------------------------------------------------

    def write_projected_mass(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_projected_mass_earth: self.write_projected_mass_earth()

        # Faceon
        if self.do_write_projected_mass_faceon: self.write_projected_mass_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_mass_earth_path(self):
        return fs.join(self.projected_path, "mass_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_mass_earth(self):
        return fs.is_file(self.projected_mass_earth_path)

    # -----------------------------------------------------------------

    def write_projected_mass_earth(self):

        """
        This function ...
        :return:
        """

        # Write
        self.stellar_mass_earth_map.saveto(self.projected_mass_earth_path)

    # -----------------------------------------------------------------

    @property
    def projected_mass_faceon_path(self):
        return fs.join(self.projected_path, "mass_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_mass_faceon(self):
        return fs.is_file(self.projected_mass_faceon_path)

    # -----------------------------------------------------------------

    def write_projected_mass_faceon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.stellar_mass_faceon_map.saveto(self.projected_mass_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_earth(self):
        return not self.has_projected_ssfr_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_faceon(self):
        return not self.has_projected_ssfr_faceon

    # -----------------------------------------------------------------

    def write_projected_ssfr(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_projected_ssfr_earth: self.write_projected_ssfr_earth()

        # Faceon
        if self.do_write_projected_ssfr_faceon: self.write_projected_ssfr_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_earth_path(self):
        return fs.join(self.projected_path, "ssfr_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_earth(self):
        return fs.is_file(self.projected_ssfr_earth_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_earth(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_earth_map.saveto(self.projected_ssfr_earth_path)

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_faceon_path(self):
        return fs.join(self.projected_path, "ssfr_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_faceon(self):
        return fs.is_file(self.projected_ssfr_faceon_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_faceon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_faceon_map.saveto(self.projected_ssfr_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_cell_fuv(self):
        return not self.has_cell_fuv

    # -----------------------------------------------------------------

    @property
    def do_write_cell_i1(self):
        return not self.has_cell_i1

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr(self):
        return not self.has_cell_sfr

    # -----------------------------------------------------------------

    @property
    def do_write_cell_mass(self):
        return not self.has_cell_mass

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr(self):
        return not self.has_cell_ssfr

    # -----------------------------------------------------------------

    def write_cell(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cell star formation rates ...")

        # Intrinsic FUV luminosity
        if self.do_write_cell_fuv: self.write_cell_fuv()

        # Observed I1 luminosity
        if self.do_write_cell_i1: self.write_cell_i1()

        # Star formation rate
        if self.do_write_cell_sfr: self.write_cell_sfr()

        # Stellar mass
        if self.do_write_cell_mass: self.write_cell_mass()

        # Specific star formation rate
        if self.do_write_cell_ssfr: self.write_cell_ssfr()

    # -----------------------------------------------------------------

    @property
    def cell_fuv_path(self):
        return fs.join(self.cell_path, "fuv.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_fuv(self):
        return fs.is_file(self.cell_fuv_path)

    # -----------------------------------------------------------------

    def write_cell_fuv(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cell intrinsic FUV luminosity ...")

        # Write
        self.fuv_data.saveto(self.cell_fuv_path)

    # -----------------------------------------------------------------

    @property
    def cell_i1_path(self):
        return fs.join(self.cell_path, "i1.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_i1(self):
        return fs.is_file(self.cell_i1_path)

    # -----------------------------------------------------------------

    def write_cell_i1(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cell observed I1 luminosity ...")

        # Write
        self.i1_data.saveto(self.cell_i1_path)

    # -----------------------------------------------------------------

    @property
    def cell_sfr_path(self):
        return fs.join(self.cell_path, "sfr.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr(self):
        return fs.is_file(self.cell_sfr_path)

    # -----------------------------------------------------------------

    def write_cell_sfr(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Writing the cell star formation rate ...")

        # Write
        self.sfr_data.saveto(self.cell_sfr_path)

    # -----------------------------------------------------------------

    @property
    def cell_mass_path(self):
        return fs.join(self.cell_path, "mass.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_mass(self):
        return fs.is_file(self.cell_mass_path)

    # -----------------------------------------------------------------

    def write_cell_mass(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cell stellar mass ...")

        # Write
        self.stellar_mass_data.saveto(self.cell_mass_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_path(self):
        return fs.join(self.cell_path, "ssfr.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr(self):
        return fs.is_file(self.cell_ssfr_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cell specific star formation rate ...")

        # Write
        self.ssfr_data.saveto(self.cell_ssfr_path)

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_map(self):
        return not self.has_cell_sfr_map

    # -----------------------------------------------------------------

    @property
    def do_write_cell_mass_map(self):
        return not self.has_cell_mass_map

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_map(self):
        return not self.has_cell_ssfr_map

    # -----------------------------------------------------------------

    def write_cell_maps(self):

        """
        Thsi function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps of the cell star formation rates ...")

        # Star formation rate
        if self.do_write_cell_sfr_map: self.write_cell_sfr_map()

        # Stellar mass
        if self.do_write_cell_mass_map: self.write_cell_mass_map()

        # Specific star formation rate
        if self.do_write_cell_ssfr_map: self.write_cell_ssfr_map()

    # -----------------------------------------------------------------

    @property
    def cell_sfr_map_path(self):
        return fs.join(self.cell_path, "sfr_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_map(self):
        return fs.is_file(self.cell_sfr_map_path)

    # -----------------------------------------------------------------

    def write_cell_sfr_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the cell star formation rate ...")

        # Write
        self.sfr_data_faceon_map.saveto(self.cell_sfr_map_path)

    # -----------------------------------------------------------------

    @property
    def cell_mass_map_path(self):
        return fs.join(self.cell_path, "mass_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cell_mass_map(self):
        return fs.is_file(self.cell_mass_map_path)

    # -----------------------------------------------------------------

    def write_cell_mass_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the cell stellar mass ...")

        # Write
        self.stellar_mass_data_faceon_map.saveto(self.cell_mass_map_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_map_path(self):
        return fs.join(self.cell_path, "ssfr_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_map(self):
        return fs.is_file(self.cell_ssfr_map_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the cell specific star formation rate ...")

        # Write
        self.ssfr_data_faceon_map.saveto(self.cell_ssfr_map_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Projected maps
        self.plot_projected()

        # Cell maps
        self.plot_cell_maps()

    # -----------------------------------------------------------------

    def plot_projected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the projected star formation rates ...")

        # Star formation rate
        self.plot_projected_sfr()

        # Stellar mass
        self.plot_projected_mass()

        # Specific star formation rate
        self.plot_projected_ssfr()

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_earth_map(self):
        return not self.has_projected_sfr_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_faceon_map(self):
        return not self.has_projected_sfr_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_sfr(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_sfr_earth_map: self.plot_projected_sfr_earth()

        # Faceon
        if self.do_plot_projected_sfr_faceon_map: self.plot_projected_sfr_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_sfr_earth_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_earth_map_plot(self):
        return fs.is_file(self.projected_sfr_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_earth(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_map(self.sfr_earth_map, path=self.projected_sfr_earth_map_plot_path, cmap="inferno", colorbar=True)

    # -----------------------------------------------------------------

    @property
    def projected_sfr_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_faceon_map_plot(self):
        return fs.is_file(self.projected_sfr_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_faceon(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_map(self.sfr_faceon_map, path=self.projected_sfr_faceon_map_plot_path, cmap="inferno", colorbar=True)

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_mass_earth(self):
        return not self.has_projected_mass_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_mass_faceon(self):
        return not self.has_projected_mass_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_mass(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_mass_earth: self.plot_projected_mass_earth()

        # Faceon
        if self.do_plot_projected_mass_faceon: self.plot_projected_mass_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_mass_earth_map_plot_path(self):
        return fs.join(self.projected_path, "mass_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_mass_earth_map_plot(self):
        return fs.is_file(self.projected_mass_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_mass_earth(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_map(self.stellar_mass_earth_map, path=self.projected_mass_earth_map_plot_path, cmap="inferno", colorbar=True)

    # -----------------------------------------------------------------

    @property
    def projected_mass_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "mass_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_mass_faceon_map_plot(self):
        return fs.is_file(self.projected_mass_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_mass_faceon(self):

        """
        Thisf unction ...
        :return:
        """

        # Plot
        plot_map(self.stellar_mass_faceon_map, path=self.projected_mass_faceon_map_plot_path, cmap="inferno", colorbar=True)

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_ssfr_earth(self):
        return not self.has_projected_ssfr_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_ssfr_faceon(self):
        return not self.has_projected_ssfr_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_ssfr(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_ssfr_earth: self.plot_projected_ssfr_earth()

        # Faceon
        if self.do_plot_projected_ssfr_faceon: self.plot_projected_ssfr_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_earth_map_plot_path(self):
        return fs.join(self.projected_path, "ssfr_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_earth_map_plot(self):
        return fs.is_file(self.projected_ssfr_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_ssfr_earth(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_map(self.ssfr_earth_map, path=self.projected_ssfr_earth_map_plot_path, cmap="inferno", colorbar=True)

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "ssfr_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_faceon_map_plot(self):
        return fs.is_file(self.projected_ssfr_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_ssfr_faceon(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_map(self.ssfr_faceon_map, path=self.projected_ssfr_faceon_map_plot_path, cmap="inferno", colorbar=True)

    # -----------------------------------------------------------------

    @property
    def do_plot_cell_sfr_map(self):
        return not self.has_cell_sfr_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_cell_mass_map(self):
        return not self.has_cell_mass_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_cell_ssfr_map(self):
        return not self.has_cell_ssfr_map_plot

    # -----------------------------------------------------------------

    def plot_cell_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the cell star formation rates ...")

        # Star formation rate
        if self.do_plot_cell_sfr_map: self.plot_cell_sfr_map()

        # Stellar mass
        if self.do_plot_cell_mass_map: self.plot_cell_mass_map()

        # Specific star formation rate
        if self.do_plot_cell_ssfr_map: self.plot_cell_ssfr_map()

    # -----------------------------------------------------------------

    @property
    def cell_sfr_map_plot_path(self):
        return fs.join(self.cell_path, "sfr_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_map_plot(self):
        return fs.is_file(self.cell_sfr_map_plot_path)

    # -----------------------------------------------------------------

    def plot_cell_sfr_map(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_map(self.sfr_data_faceon_map, path=self.cell_sfr_map_plot_path, cmap="inferno", colorbar=True)

    # -----------------------------------------------------------------

    @property
    def cell_mass_map_plot_path(self):
        return fs.join(self.cell_path, "mass_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_cell_mass_map_plot(self):
        return fs.is_file(self.cell_mass_map_plot_path)

    # -----------------------------------------------------------------

    def plot_cell_mass_map(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_map(self.stellar_mass_data_faceon_map, path=self.cell_mass_map_plot_path, cmap="inferno", colorbar=True)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_map_plot_path(self):
        return fs.join(self.cell_path, "ssfr_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_map_plot(self):
        return fs.is_file(self.cell_ssfr_map_plot_path)

    # -----------------------------------------------------------------

    def plot_cell_ssfr_map(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_map(self.ssfr_data_faceon_map, path=self.cell_ssfr_map_plot_path, cmap="inferno", colorbar=True)

# -----------------------------------------------------------------
