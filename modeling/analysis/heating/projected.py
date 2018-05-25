#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating.projected Contains the ProjectedDustHeatingAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.tools import tables
from ....magic.plot.imagegrid import StandardImageGridPlotter
from ....core.tools.utils import lazyproperty
from ....magic.core.frame import Frame

# -----------------------------------------------------------------

class ProjectedDustHeatingAnalyser(DustHeatingAnalysisComponent):
    
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
        super(ProjectedDustHeatingAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Total
        self.total_absorptions_earth = None
        self.total_absorptions_faceon = None
        self.total_absorptions_edgeon = None

        # Young
        self.young_absorptions_earth = None
        self.young_absorptions_faceon = None
        self.young_absorptions_edgeon = None

        # Ionizing
        self.ionizing_absorptions_earth = None
        self.ionizing_absorptions_faceon = None
        self.ionizing_absorptions_edgeon = None

        # Internal SFR absorption
        self.internal_absorptions_earth = None
        self.internal_absorptions_faceon = None
        self.internal_absorptions_edgeon = None

        # Maps of heating fraction
        self.map_earth = None
        self.map_faceon = None
        self.map_edgeon = None

        # Cubes of spectral heating
        self.cube_earth = None
        self.cube_faceon = None
        self.cube_edgeon = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get the absorption maps
        self.get_absorptions()

        # Get the maps of the heating fraction
        self.get_maps()

        # Get the cube of the heating fraction per wavelength
        self.get_cubes()

        # Show
        self.show()

        # 4. Writing
        self.write()

        # 5. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ProjectedDustHeatingAnalyser, self).setup()

        # Check
        if (not self.do_earth) and (not self.do_faceon) and (not self.do_edgeon): raise RuntimeError("Cannot create any map: not enough simulation data is present")

    # -----------------------------------------------------------------

    @property
    def do_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_earth_cube_contributions_all

    # -----------------------------------------------------------------

    @lazyproperty
    def do_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_faceon_cube_contributions_all

    # -----------------------------------------------------------------

    @lazyproperty
    def do_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_edgeon_cube_contributions_all

    # -----------------------------------------------------------------

    @lazyproperty
    def has_earth_cube_all(self):

        """
        This function ...
        :return:
        """

        return self.has_earth_cube_total and self.has_earth_cube_young and self.has_earth_cube_ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def has_earth_cube_contributions_all(self):

        """
        Thisn function ...
        :return:
        """

        return self.has_earth_cube_contributions_total and self.has_earth_cube_contributions_young and self.has_earth_cube_contributions_ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def has_faceon_cube_all(self):

        """
        This function ...
        :return:
        """

        return self.has_faceon_cube_total and self.has_faceon_cube_young and self.has_faceon_cube_ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def has_faceon_cube_contributions_all(self):

        """
        Thisfunction ...
        :return:
        """

        return self.has_faceon_cube_contributions_total and self.has_faceon_cube_contributions_young and self.has_faceon_cube_contributions_ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def has_edgeon_cube_all(self):

        """
        Thisn function ...
        :return:
        """

        return self.has_edgeon_cube_total and self.has_edgeon_cube_young and self.has_edgeon_cube_ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def has_edgeon_cube_contributions_all(self):

        """
        This function ...
        :return:
        """

        return self.has_edgeon_cube_contributions_total and self.has_edgeon_cube_contributions_young and self.has_edgeon_cube_contributions_ionizing

    # -----------------------------------------------------------------

    @property
    def total_simulations(self):

        """
        This function ...
        :return:
        """

        return self.model.total_simulations

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_total(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_contributions_total(self):

        """
        This function ...
        :return:
        """

        return self.has_earth_cube_total and self.total_simulations.has_other_observed_cube_contributions

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_total(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_faceon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_contributions_total(self):

        """
        This function ...
        :return:
        """

        return self.has_faceon_cube_total and self.total_simulations.has_other_observed_cube_contributions_faceon

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_total(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.has_edgeon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_contributions_total(self):

        """
        This function ...
        :return:
        """

        return self.has_edgeon_cube_total and self.total_simulations.has_other_observed_cube_contributions_edgeon

    # -----------------------------------------------------------------

    @property
    def young_simulations(self):

        """
        This function ...
        :return:
        """

        return self.model.young_simulations

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_young(self):

        """
        Thisfunction ...
        :return:
        """

        return self.young_simulations.has_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_contributions_young(self):

        """
        This function ...
        :return:
        """

        return self.has_earth_cube_young and self.young_simulations.has_other_observed_cube_contributions

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_faceon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_contributions_young(self):

        """
        This function ...
        :return:
        """

        return self.has_faceon_cube_young and self.young_simulations.has_other_observed_cube_contributions_faceon

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_young(self):

        """
        This function ...
        :return:
        """

        return self.young_simulations.has_edgeon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_contributions_young(self):

        """
        This function ...
        :return:
        """

        return self.has_edgeon_cube_young and self.young_simulations.has_other_observed_cube_contributions_edgeon

    # -----------------------------------------------------------------

    @property
    def ionizing_simulations(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_simulations

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_ionizing(self):

        """
        Thisf unction ...
        :return:
        """

        return self.ionizing_simulations.has_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_contributions_ionizing(self):

        """
        Thisfunction ...
        :return:
        """

        return self.has_earth_cube_ionizing and self.ionizing_simulations.has_other_observed_cube_contributions

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_simulations.has_faceon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_contributions_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.has_faceon_cube_ionizing and self.ionizing_simulations.has_other_observed_cube_contributions_faceon

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_simulations.has_edgeon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_contributions_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.has_edgeon_cube_ionizing and self.ionizing_simulations.has_other_observed_cube_contributions_edgeon

    # -----------------------------------------------------------------

    def get_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the absorption maps ...")

        # Total
        self.get_total_absorptions()

        # Young
        self.get_young_absorptions()

        # Ionizing
        self.get_ionizing_absorptions()

        # Internal
        self.get_internal_absorptions()

    # -----------------------------------------------------------------

    def get_total_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the total absorption maps ..")

        # Earth
        self.get_total_absorptions_earth()

        # Face-on
        self.get_total_absorptions_faceon()

        # Edge-on
        self.get_total_absorptions_edgeon()

    # -----------------------------------------------------------------

    def get_total_absorptions_earth(self):
        
        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Getting the total absorption map from the earth projection ...")

        # Get total absorption map
        if self.has_total_absorptions_earth: self.load_total_absorptions_earth()

        # Calculate
        else: self.calculate_total_absorptions_earth()

    # -----------------------------------------------------------------

    def load_total_absorptions_earth(self):

        """
        Thisf unction ...
        :return:
        """

        # Load
        self.total_absorptions_earth = Frame.from_file(self.total_absorptions_earth_path)

    # -----------------------------------------------------------------

    def calculate_total_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.total_absorptions_earth = self.total_simulations.observed_dust_frame

    # -----------------------------------------------------------------

    def get_total_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the total absorption map from the faceon projection ...")

        # Load
        if self.has_total_absorptions_faceon: self.load_total_absorptions_faceon()

        # Calculate
        else: self.calculate_total_absorptions_faceon()

    # -----------------------------------------------------------------

    def load_total_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        self.total_absorptions_faceon = Frame.from_file(self.total_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def calculate_total_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.total_absorptions_faceon = self.total_simulations.faceon_observed_dust_sed

    # -----------------------------------------------------------------

    def get_total_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the total absorption maps from the edgeon projection ...")

        # Load
        if self.has_total_absorptions_edgeon: self.load_total_absorptions_edgeon()

        # Calculate
        else: self.calculate_total_absorptions_edgeon()

    # -----------------------------------------------------------------

    def load_total_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def calculate_total_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def get_young_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the young stellar absorption maps ...")

        # Earth
        self.get_young_absorptions_earth()

        # Face-on
        self.get_young_absorptions_faceon()

        # Edge-on
        self.get_young_absorptions_edgeon()

    # -----------------------------------------------------------------

    def get_young_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_young_absorptions_earth: self.load_young_absorptions()

        # Calculate
        else: self.calculate_young_absorptions()

    # -----------------------------------------------------------------

    def load_young_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        # Load
        self.young_absorptions_earth = Frame.from_file(self.young_absorptions_earth_path)

    # -----------------------------------------------------------------

    def calculate_young_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.young_absorptions_earth = self.young_simulations.observed_dust_frame

    # -----------------------------------------------------------------

    def get_young_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_young_absorptions_faceon: self.load_young_absorptions_faceon()

        # Calculate
        else: self.calculate_young_absorptions_faceon()

    # -----------------------------------------------------------------

    def load_young_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        self.young_absorptions_faceon = Frame.from_file(self.young_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def calculate_young_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.young_absorptions_faceon = self.young_simulations.faceon_observed_dust_frame

    # -----------------------------------------------------------------

    def get_young_absorption_edgeon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_young_absorptions_edgeon: self.load_young_absorptions_edgeon()

        # Calculate
        else: self.calculate_young_absorptions_edgeon()

    # -----------------------------------------------------------------

    def load_young_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def calculate_young_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def get_ionizing_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the ionizing stellar absorption maps ...")

        # Earth
        self.get_ionizing_absorptions_earth()

        # Face-on
        self.get_ionizing_absorptions_faceon()

        # Edge-on
        self.get_ionizing_absorptions_edgeon()

    # -----------------------------------------------------------------

    def get_ionizing_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_ionizing_absorptions_earth: self.load_ionizing_absorptions_earth()

        # Calculate
        else: self.calculate_ionizing_absorptions_earth()

    # -----------------------------------------------------------------

    def load_ionizing_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        # Load
        self.ionizing_absorptions_earth = Frame.from_file(self.ionizing_absorptions_earth_path)

    # -----------------------------------------------------------------

    def calculate_ionizing_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.ionizing_absorptions_earth = self.ionizing_simulations.observed_dust_frame

    # -----------------------------------------------------------------

    def get_ionizing_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_ionizing_absorptions_faceon: self.load_ionizing_absorptions_faceon()

        # Calculate
        else: self.calculate_ionizing_absorptions_faceon()

    # -----------------------------------------------------------------

    def load_ionizing_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        self.ionizing_absorptions_faceon = Frame.from_file(self.ionizing_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def calculate_ionizing_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.ionizing_absorptions_faceon = self.ionizing_simulations.faceon_observed_dust_frame

    # -----------------------------------------------------------------

    def get_internal_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.get_internal_absorptions_earth()

        # Faceon
        self.get_internal_absorptions_faceon()

        # Edge-on
        self.get_internal_absorptions_edgeon()

    # -----------------------------------------------------------------

    def get_internal_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_internal_absorptions_earth: self.load_internal_absorptions_earth()

        # Calculate
        else: self.calculate_internal_absorptions_earth()

    # -----------------------------------------------------------------

    def load_internal_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        # Load
        self.internal_absorptions_earth = Frame.from_file(self.internal_absorptions_earth_path)

    # -----------------------------------------------------------------

    def calculate_internal_absorption(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get_internal_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_internal_absorptions_faceon: self.load_internal_absorptions_faceon()

        # Calculate
        else: self.calculate_internal_absorptions_faceon()

    # -----------------------------------------------------------------

    def load_internal_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        self.internal_absorptions_faceon = Frame.from_file(self.internal_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def calculate_internal_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def get_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.get_map_earth()

        # Face-on
        self.get_map_faceon()

        # Edge-on
        self.get_map_edgeon()

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        return self.young_absorptions_earth + self.ionizing_absorptions_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_earth_with_internal(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_absorptions_earth + self.internal_absorptions_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorptions_earth_with_internal(self):

        """
        This function ...
        :return:
        """

        return self.total_absorptions_earth + self.internal_absorptions_earth

    # -----------------------------------------------------------------

    def get_earth_map(self):

        """
        This function ...
        :return:
        """

        if self.has_map_earth: self.load_map_earth()

        else: self.calculate_map_earth()

    # -----------------------------------------------------------------

    def calculate_map_earth(self):

        """
        Thisn function ...
        :return:
        """

        # Calculate
        self.map_earth = self.unevolved_absorptions_earth_with_internal / self.total_absorptions_earth_with_internal

    # -----------------------------------------------------------------

    def get_faceon_map(self):

        """
        This function ...
        :return:
        """

        # Calculate
        self.map_faceon = self.unevolved_absorptions_faceon_with_internal /

    # -----------------------------------------------------------------

    def get_edgeon_map(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def get_cubes(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.get_earth_cube()

        # Face-on
        self.get_faceon_cube()

    # -----------------------------------------------------------------

    def get_earth_cube(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def get_faceon_cube(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write
        self.write_absorptions()

    # -----------------------------------------------------------------

    @lazyproperty
    def absorptions_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.projected_heating_path, "absorptions")

    # -----------------------------------------------------------------

    def write_absorptions(self):

        """
        This function ...
        :return:
        """

        # Total
        self.write_total_absorptions()

        # Young
        self.write_young_absorptions()

        # Ionizing
        self.write_ionizing_absorptions()

        # Internal
        self.write_internal_absorptions()

    # -----------------------------------------------------------------

    def write_total_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the total absorption maps ...")

        # Earth
        self.write_total_absorption_earth()

        # Face-on
        self.write_total_absorption_faceon()

        # Edge-on
        self.write_total_absorption_edgeon()

    # -----------------------------------------------------------------

    @property
    def total_absorption_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "total_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_absorption_earth(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.total_absorption_earth_path)

    # -----------------------------------------------------------------

    def remove_total_absorption(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.total_absorption_earth_path)

    # -----------------------------------------------------------------

    def write_total_absorption_earth(self):

        """
        This function ...
        :return:
        """

        self.total_absorption.saveto(self.total_absorption_earth_path)

    # -----------------------------------------------------------------

    @property
    def total_absorption_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "total_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.total_absorption_faceon_path)

    # -----------------------------------------------------------------

    def remove_total_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.total_absorption_faceon_path)

    # -----------------------------------------------------------------

    def write_total_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.total_absorption_faceon.saveto(self.total_absorption_faceon_path)

    # -----------------------------------------------------------------

    def write_total_absorption_edgeon(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_young_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.write_young_absorption()

        # Faceon
        self.write_young_absorption_faceon()

    # -----------------------------------------------------------------

    @property
    def young_absorption_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "young_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_absorption_earth(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_absorption_earth_path)

    # -----------------------------------------------------------------

    def remove_young_absorption(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.young_absorption_earth_path)

    # -----------------------------------------------------------------

    def write_young_absorption(self):

        """
        This function ...
        :return:
        """

        self.young_absorption.saveto(self.young_absorption_path)

    # -----------------------------------------------------------------

    @property
    def young_absorption_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "young_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_absorption_faceon_path)

    # -----------------------------------------------------------------

    def remove_young_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.young_absorption_faceon_path)

    # -----------------------------------------------------------------

    def write_young_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        self.young_absorption_faceon.saveto(self.young_absorption_faceon_path)

    # -----------------------------------------------------------------

    def write_ionizing_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.write_ionizing_absorption()

        # Face-on
        self.write_ionizing_absorption_faceon()

    # -----------------------------------------------------------------

    @property
    def ionizing_absorption_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "ionizing.fits")

    # -----------------------------------------------------------------

    @property
    def has_ionizing_absorption(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ionizing_absorption_path)

    # -----------------------------------------------------------------

    def remove_ionizing_absorption(self):

        """
        This function ...
        :return:
        """

        return fs.remove_file(self.ionizing_absorption_path)

    # -----------------------------------------------------------------

    def write_ionizing_absorption(self):

        """
        This function ...
        :return:
        """

        self.ionizing_absorption.saveto(self.ionizing_absorption_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_absorption_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "ionizing_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_ionizing_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ionizing_absorption_faceon_path)

    # -----------------------------------------------------------------

    def remove_ionizing_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.ionizing_absorption_faceon_path)

    # -----------------------------------------------------------------

    def write_ionizing_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        self.ionizing_absorption_faceon.saveto(self.ionizing_absorption_faceon_path)

    # -----------------------------------------------------------------

    def write_internal_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.write_internal_absorption()

        # Face-on
        self.write_internal_absorption_faceon()

    # -----------------------------------------------------------------

    @property
    def internal_absorption_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.absorptions_path, "internal.fits")

    # -----------------------------------------------------------------

    @property
    def has_internal_absorption(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.internal_absorption_path)

    # -----------------------------------------------------------------

    def write_internal_absorption(self):

        """
        This function ...
        :return:
        """

        self.internal_absorption.saveto(self.internal_absorption_path)

    # -----------------------------------------------------------------

    @property
    def internal_absorption_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "internal_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_internal_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.internal_absorption_faceon_path)

    # -----------------------------------------------------------------

    def remove_internal_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.remove_file(self.internal_absorption_faceon_path)

    # -----------------------------------------------------------------

    def write_internal_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.internal_absorption_faceon.saveto(self.internal_absorption_faceon_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
