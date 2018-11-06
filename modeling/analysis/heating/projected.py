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

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.tools.utils import lazyproperty, memoize_method
from ....magic.core.frame import Frame
from ....magic.core.datacube import DataCube
from ....magic.core.list import uniformize
from ....core.units.parsing import parse_quantity
from ....core.basics.curve import WavelengthCurve
from ....magic.tools import plotting

# -----------------------------------------------------------------

max_wavelength_absorption = parse_quantity("5 micron")
min_wavelength_emission = parse_quantity("10 micron")

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

        # Extra
        self.extra_absorptions_earth = None
        self.extra_absorptions_faceon = None
        self.extra_absorptions_edgeon = None

        # Maps of heating fraction
        self.map_earth = None
        self.map_faceon = None
        self.map_edgeon = None

        # Maps of heating fraction (diffuse)
        self.map_earth_diffuse = None
        self.map_faceon_diffuse = None
        self.map_edgeon_diffuse = None

        # Maps of heating fraction (extra)
        self.map_earth_extra = None
        self.map_faceon_extra = None
        self.map_edgeon_extra = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Get the absorption maps
        self.get_absorptions()

        # 2. Get the maps of the heating fraction
        self.get_maps()

        # 3. Show
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
        return self.has_earth_cube_all
        #return self.has_earth_cube_contributions_all

    # -----------------------------------------------------------------

    @lazyproperty
    def do_faceon(self):
        return self.has_faceon_cube_all
        #return self.has_faceon_cube_contributions_all

    # -----------------------------------------------------------------

    @lazyproperty
    def do_edgeon(self):
        return self.has_edgeon_cube_all
        #return self.has_edgeon_cube_contributions_all

    # -----------------------------------------------------------------
    # HAS EARTH CUBES
    # -----------------------------------------------------------------

    @lazyproperty
    def has_earth_cube_all(self):
        return self.has_earth_cube_total and self.has_earth_cube_young and self.has_earth_cube_ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def has_earth_cube_contributions_all(self):
        return self.has_earth_cube_contributions_total and self.has_earth_cube_contributions_young and self.has_earth_cube_contributions_ionizing

    # -----------------------------------------------------------------
    # HAS FACEON CUBES
    # -----------------------------------------------------------------

    @lazyproperty
    def has_faceon_cube_all(self):
        return self.has_faceon_cube_total and self.has_faceon_cube_young and self.has_faceon_cube_ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def has_faceon_cube_contributions_all(self):
        return self.has_faceon_cube_contributions_total and self.has_faceon_cube_contributions_young and self.has_faceon_cube_contributions_ionizing

    # -----------------------------------------------------------------
    # HAS EDGEON CUBES
    # -----------------------------------------------------------------

    @lazyproperty
    def has_edgeon_cube_all(self):
        return self.has_edgeon_cube_total and self.has_edgeon_cube_young and self.has_edgeon_cube_ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def has_edgeon_cube_contributions_all(self):
        return self.has_edgeon_cube_contributions_total and self.has_edgeon_cube_contributions_young and self.has_edgeon_cube_contributions_ionizing

    # -----------------------------------------------------------------
    # TOTAL SIMULATIONS
    # -----------------------------------------------------------------

    @property
    def total_simulations(self):
        return self.model.total_simulations

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_total(self):
        return self.total_simulations.has_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_contributions_total(self):
        return self.has_earth_cube_total and self.total_simulations.has_other_observed_cube_contributions

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_total(self):
        return self.total_simulations.has_faceon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_contributions_total(self):
        return self.has_faceon_cube_total and self.total_simulations.has_other_observed_cube_contributions_faceon

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_total(self):
        return self.total_simulations.has_edgeon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_contributions_total(self):
        return self.has_edgeon_cube_total and self.total_simulations.has_other_observed_cube_contributions_edgeon

    # -----------------------------------------------------------------
    # YOUNG SIMULATIONS
    # -----------------------------------------------------------------

    @property
    def young_simulations(self):
        return self.model.young_simulations

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_young(self):
        return self.young_simulations.has_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_contributions_young(self):
        return self.has_earth_cube_young and self.young_simulations.has_other_observed_cube_contributions

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_young(self):
        return self.young_simulations.has_faceon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_contributions_young(self):
        return self.has_faceon_cube_young and self.young_simulations.has_other_observed_cube_contributions_faceon

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_young(self):
        return self.young_simulations.has_edgeon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_contributions_young(self):
        return self.has_edgeon_cube_young and self.young_simulations.has_other_observed_cube_contributions_edgeon

    # -----------------------------------------------------------------
    # IONIZING SIMULATIONS
    # -----------------------------------------------------------------

    @property
    def ionizing_simulations(self):
        return self.model.sfr_simulations

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_ionizing(self):
        return self.ionizing_simulations.has_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_earth_cube_contributions_ionizing(self):
        return self.has_earth_cube_ionizing and self.ionizing_simulations.has_other_observed_cube_contributions

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_ionizing(self):
        return self.ionizing_simulations.has_faceon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_contributions_ionizing(self):
        return self.has_faceon_cube_ionizing and self.ionizing_simulations.has_other_observed_cube_contributions_faceon

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_ionizing(self):
        return self.ionizing_simulations.has_edgeon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_contributions_ionizing(self):
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
        if self.do_earth: self.get_total_absorptions_earth()

        # Face-on
        if self.do_faceon: self.get_total_absorptions_faceon()

        # Edge-on
        if self.do_edgeon: self.get_total_absorptions_edgeon()

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
        #self.total_absorptions_earth = self.total_simulations.observed_dust_frame
        self.total_absorptions_earth = self.model.total_dust_luminosity_map_earth

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
        #self.total_absorptions_faceon = self.total_simulations.faceon_observed_dust_frame
        self.total_absorptions_faceon = self.model.total_dust_luminosity_map_faceon

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

        # Load
        self.total_absorptions_edgeon = Frame.from_file(self.total_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def calculate_total_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        #self.total_absorptions_edgeon = self.total_simulations.edgeon_observed_dust_sed
        self.total_absorptions_edgeon = self.model.total_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------

    def get_young_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the young stellar absorption maps ...")

        # Earth
        if self.do_earth: self.get_young_absorptions_earth()

        # Face-on
        if self.do_faceon: self.get_young_absorptions_faceon()

        # Edge-on
        if self.do_edgeon: self.get_young_absorptions_edgeon()

    # -----------------------------------------------------------------

    def get_young_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_young_absorptions_earth: self.load_young_absorptions_earth()

        # Calculate
        else: self.calculate_young_absorptions_earth()

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
        #self.young_absorptions_earth = self.young_simulations.observed_dust_frame
        self.young_absorptions_earth = self.model.young_dust_luminosity_map_earth

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
        #self.young_absorptions_faceon = self.young_simulations.faceon_observed_dust_frame
        self.young_absorptions_faceon = self.model.young_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    def get_young_absorptions_edgeon(self):

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

        self.young_absorptions_edgeon = Frame.from_file(self.young_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def calculate_young_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        #self.young_absorptions_edgeon = self.young_simulations.edgeon_observed_dust_frame
        self.young_absorptions_edgeon = self.model.young_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------

    def get_ionizing_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the ionizing stellar absorption maps ...")

        # Earth
        if self.do_earth: self.get_ionizing_absorptions_earth()

        # Face-on
        if self.do_faceon: self.get_ionizing_absorptions_faceon()

        # Edge-on
        if self.do_edgeon: self.get_ionizing_absorptions_edgeon()

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
        #self.ionizing_absorptions_earth = self.ionizing_simulations.observed_dust_frame
        self.ionizing_absorptions_earth = self.model.sfr_dust_luminosity_map_earth

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

    def calculate_ionizing_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        #self.ionizing_absorptions_faceon = self.ionizing_simulations.faceon_observed_dust_frame
        self.ionizing_absorptions_faceon = self.model.sfr_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    def get_ionizing_absorptions_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        # Load
        if self.has_ionizing_absorptions_edgeon: self.load_ionizing_absorptions_edgeon()

        # Create
        else: self.calculate_ionizing_absorptions_edgeon()

    # -----------------------------------------------------------------

    def load_ionizing_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        self.ionizing_absorptions_edgeon = Frame.from_file(self.ionizing_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def calculate_ionizing_absorptions_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        self.ionizing_absorptions_edgeon = self.model.sfr_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------

    def get_internal_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.get_internal_absorptions_earth()

        # Faceon
        if self.do_faceon: self.get_internal_absorptions_faceon()

        # Edge-on
        if self.do_edgeon: self.get_internal_absorptions_edgeon()

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

    def calculate_internal_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        self.internal_absorptions_earth = self.model.sfr_intrinsic_dust_luminosity_map_earth

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

        self.internal_absorptions_faceon = self.model.sfr_intrinsic_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    def get_internal_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_internal_absorptions_edgeon: self.load_internal_absorptions_edgeon()

        # Calculate
        else: self.calculate_internal_absorptions_edgeon()

    # -----------------------------------------------------------------

    def load_internal_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        self.internal_absorptions_edgeon = Frame.from_file(self.internal_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def calculate_internal_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        self.internal_absorptions_edgeon = self.model.sfr_intrinsic_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------

    def get_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:

            self.get_map_earth()
            self.get_map_earth_diffuse()

        # Face-on
        if self.do_faceon:

            self.get_map_faceon()
            self.get_map_faceon_diffuse()

        # Edge-on
        if self.do_edgeon:

            self.get_map_edgeon()
            self.get_map_edgeon_diffuse()

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        #return self.young_absorptions_earth + self.ionizing_absorptions_earth
        young, ionizing = uniformize(self.young_absorptions_earth, self.ionizing_absorptions_earth, distance=self.galaxy_distance)
        return young + ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_earth_diffuse(self):

        """
        This function ...
        :return:
        """

        #return self.unevolved_absorptions_earth - self.internal_absorptions_earth
        unevolved, internal = uniformize(self.unevolved_absorptions_earth, self.internal_absorptions_earth, distance=self.galaxy_distance)
        return unevolved - internal

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorptions_earth_diffuse(self):

        """
        This function ...
        :return:
        """

        #return self.total_absorptions_earth - self.internal_absorptions_earth
        total, internal = uniformize(self.total_absorptions_earth, self.internal_absorptions_earth, distance=self.galaxy_distance)
        return total - internal

    # -----------------------------------------------------------------

    def get_map_earth(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_map_earth: self.load_map_earth()

        # Create
        else: self.calculate_map_earth()

    # -----------------------------------------------------------------

    def load_map_earth(self):

        """
        This function ...
        :return:
        """

        self.map_earth = Frame.from_file(self.map_earth_path)

    # -----------------------------------------------------------------

    def calculate_map_earth(self):

        """
        Thisn function ...
        :return:
        """

        # Calculate
        #self.map_earth = self.unevolved_absorptions_earth_with_internal / self.total_absorptions_earth_with_internal
        self.map_earth = self.unevolved_absorptions_earth / self.total_absorptions_earth

    # -----------------------------------------------------------------

    def get_map_earth_diffuse(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_map_earth_diffuse: self.load_map_earth_diffuse()

        # Create
        else: self.calculate_map_earth_diffuse()

    # -----------------------------------------------------------------

    def load_map_earth_diffuse(self):

        """
        This function ...
        :return:
        """

        self.map_earth_diffuse = Frame.from_file(self.map_earth_diffuse_path)

    # -----------------------------------------------------------------

    def calculate_map_earth_diffuse(self):

        """
        This function ...
        :return:
        """

        self.map_earth_diffuse = self.unevolved_absorptions_earth_diffuse / self.total_absorptions_earth_diffuse

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        #print(self.young_absorptions_faceon.distance)
        #print(self.ionizing_absorptions_faceon.distance)
        #return self.young_absorptions_faceon + self.ionizing_absorptions_faceon
        young, ionizing = uniformize(self.young_absorptions_faceon, self.ionizing_absorptions_faceon, distance=self.galaxy_distance)
        return young + ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_faceon_diffuse(self):

        """
        Thisn function ...
        :return:
        """

        #print(self.unevolved_absorptions_faceon.distance)
        #print(self.internal_absorptions_faceon.distance)
        #return self.unevolved_absorptions_faceon - self.internal_absorptions_faceon
        unevolved, internal = uniformize(self.unevolved_absorptions_faceon, self.internal_absorptions_faceon)
        return unevolved - internal

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorptions_faceon_diffuse(self):

        """
        This function ...
        :return:
        """

        #return self.total_absorptions_faceon - self.internal_absorptions_faceon
        total, internal = uniformize(self.total_absorptions_faceon, self.internal_absorptions_faceon)
        return total - internal

    # -----------------------------------------------------------------

    def get_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_map_faceon: self.load_map_faceon()

        # Create
        else: self.calculate_map_faceon()

    # -----------------------------------------------------------------

    def load_map_faceon(self):

        """
        This function ...
        :return:
        """

        self.map_faceon = Frame.from_file(self.map_faceon_path)

    # -----------------------------------------------------------------

    def calculate_map_faceon(self):

        """
        Thisf unction ...
        :return:
        """

        # Calculate
        #self.map_faceon = self.unevolved_absorptions_faceon_with_internal / self.total_absorptions_faceon_with_internal
        self.map_faceon = self.unevolved_absorptions_faceon / self.total_absorptions_faceon

    # -----------------------------------------------------------------

    def get_map_faceon_diffuse(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_map_faceon_diffuse: self.load_map_faceon_diffuse()

        # Create
        else: self.calculate_map_faceon_diffuse()

    # -----------------------------------------------------------------

    def load_map_faceon_diffuse(self):

        """
        This function ...
        :return:
        """

        self.map_faceon_diffuse = Frame.from_file(self.map_faceon_diffuse_path)

    # -----------------------------------------------------------------

    def calculate_map_faceon_diffuse(self):

        """
        This function ...
        :return:
        """

        self.map_faceon_diffuse = self.unevolved_absorptions_faceon_diffuse / self.total_absorptions_faceon_diffuse

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        #return self.young_absorptions_edgeon + self.ionizing_absorptions_edgeon
        young, ionizing = uniformize(self.young_absorptions_edgeon, self.ionizing_absorptions_edgeon)
        return young + ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_edgeon_diffuse(self):

        """
        Thisf unction ...
        :return:
        """

        #return self.unevolved_absorptions_edgeon - self.internal_absorptions_edgeon
        unevolved, internal = uniformize(self.unevolved_absorptions_edgeon, self.internal_absorptions_edgeon)
        return unevolved - internal

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorptions_edgeon_diffuse(self):

        """
        This function ...
        :return:
        """

        #return self.total_absorptions_edgeon - self.internal_absorptions_edgeon
        total, internal = uniformize(self.total_absorptions_edgeon, self.internal_absorptions_edgeon)
        return total - internal

    # -----------------------------------------------------------------

    def get_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_map_edgeon: self.load_map_edgeon()

        # Create
        else: self.calculate_map_edgeon()

    # -----------------------------------------------------------------

    def load_map_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        self.map_edgeon = Frame.from_file(self.map_edgeon_path)

    # -----------------------------------------------------------------

    def calculate_map_edgeon(self):

        """
        This function ...
        :return:
        """

        #self.map_edgeon = self.unevolved_absorptions_edgeon_with_internal / self.total_absorptions_edgeon_with_internal
        self.map_edgeon = self.unevolved_absorptions_edgeon / self.total_absorptions_edgeon

    # -----------------------------------------------------------------

    def get_map_edgeon_diffuse(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_map_edgeon_diffuse: self.load_map_edgeon_diffuse()

        # Create
        else: self.calculate_map_edgeon_diffuse()

    # -----------------------------------------------------------------

    def load_map_edgeon_diffuse(self):

        """
        This function ...
        :return:
        """

        self.map_edgeon_diffuse = Frame.from_file(self.map_edgeon_diffuse_path)

    # -----------------------------------------------------------------

    def calculate_map_edgeon_diffuse(self):

        """
        This function ...
        :return:
        """

        self.map_edgeon_diffuse = self.unevolved_absorptions_edgeon_diffuse / self.total_absorptions_edgeon_diffuse

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

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

        # Write maps
        self.write_maps()

    # -----------------------------------------------------------------

    @lazyproperty
    def absorptions_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.projected_heating_path, "absorptions")

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.projected_heating_path, "maps")

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
        if self.do_earth: self.write_total_absorptions_earth()

        # Face-on
        if self.do_faceon: self.write_total_absorptions_faceon()

        # Edge-on
        if self.do_edgeon: self.write_total_absorptions_edgeon()

    # -----------------------------------------------------------------

    @property
    def total_absorptions_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "total_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.total_absorptions_earth_path)

    # -----------------------------------------------------------------

    def remove_total_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.total_absorptions_earth_path)

    # -----------------------------------------------------------------

    def write_total_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        self.total_absorptions_earth.saveto(self.total_absorptions_earth_path)

    # -----------------------------------------------------------------

    @property
    def total_absorptions_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "total_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.total_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def remove_total_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.total_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def write_total_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.total_absorptions_faceon.saveto(self.total_absorptions_faceon_path)

    # -----------------------------------------------------------------

    @property
    def total_absorptions_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "total_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.total_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def remove_total_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.total_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def write_total_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        self.total_absorptions_edgeon.saveto(self.total_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def write_young_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.write_young_absorptions_earth()

        # Faceon
        if self.do_faceon: self.write_young_absorptions_faceon()

        # Edgeon
        if self.do_edgeon: self.write_young_absorptions_edgeon()

    # -----------------------------------------------------------------

    @property
    def young_absorptions_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "young_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_absorptions_earth_path)

    # -----------------------------------------------------------------

    def remove_young_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.young_absorptions_earth_path)

    # -----------------------------------------------------------------

    def write_young_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        self.young_absorptions_earth.saveto(self.young_absorptions_earth_path)

    # -----------------------------------------------------------------

    @property
    def young_absorptions_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "young_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def remove_young_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.young_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def write_young_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        self.young_absorptions_faceon.saveto(self.young_absorptions_faceon_path)

    # -----------------------------------------------------------------

    @property
    def young_absorptions_edgeon_path(self):

        """
        This function ...
        """

        return fs.join(self.absorptions_path, "young_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def remove_young_absorptions_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        fs.remove_file(self.young_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def write_young_absorptions_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        self.young_absorptions_edgeon.saveto(self.young_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def write_ionizing_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.write_ionizing_absorptions_earth()

        # Face-on
        if self.do_faceon: self.write_ionizing_absorptions_faceon()

        # Edge-on
        if self.do_edgeon: self.write_ionizing_absorptions_edgeon()

    # -----------------------------------------------------------------

    @property
    def ionizing_absorptions_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "ionizing_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_ionizing_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ionizing_absorptions_earth_path)

    # -----------------------------------------------------------------

    def remove_ionizing_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        return fs.remove_file(self.ionizing_absorptions_earth_path)

    # -----------------------------------------------------------------

    def write_ionizing_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        self.ionizing_absorptions_earth.saveto(self.ionizing_absorptions_earth_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_absorptions_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "ionizing_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_ionizing_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ionizing_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def remove_ionizing_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.ionizing_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def write_ionizing_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        self.ionizing_absorptions_faceon.saveto(self.ionizing_absorptions_faceon_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_absorptions_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "ionizing_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_ionizing_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ionizing_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def remove_ionizing_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.ionizing_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def write_ionizing_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        self.ionizing_absorptions_edgeon.saveto(self.ionizing_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def write_internal_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.write_internal_absorptions_earth()

        # Face-on
        if self.do_faceon: self.write_internal_absorptions_faceon()

        # Edge-on
        if self.do_edgeon: self.write_internal_absorptions_edgeon()

    # -----------------------------------------------------------------

    @property
    def internal_absorptions_earth_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.absorptions_path, "internal_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_internal_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.internal_absorptions_earth_path)

    # -----------------------------------------------------------------

    def write_internal_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        self.internal_absorptions_earth.saveto(self.internal_absorptions_earth_path)

    # -----------------------------------------------------------------

    @property
    def internal_absorptions_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "internal_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_internal_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.internal_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def remove_internal_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.remove_file(self.internal_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def write_internal_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.internal_absorptions_faceon.saveto(self.internal_absorptions_faceon_path)

    # -----------------------------------------------------------------

    @property
    def internal_absorptions_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "internal_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_internal_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.internal_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def remove_internal_absorptions_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        fs.remove_file(self.internal_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def write_internal_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.internal_absorptions_edgeon.saveto(self.internal_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps of the heating fraction ...")

        # Earth
        if self.do_earth:

            self.write_map_earth()
            self.write_map_earth_diffuse()

        # Face-on
        if self.do_faceon:

            self.write_map_faceon()
            self.write_map_faceon_diffuse()

        # Edge-on
        if self.do_edgeon:

            self.write_map_edgeon()
            self.write_map_edgeon_diffuse()

    # -----------------------------------------------------------------

    @property
    def map_earth_path(self):

        """
        Thisn function ...
        :return:
        """

        return fs.join(self.maps_path, "earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_earth(self):

        """
        Thisnfunctino ...
        :return:
        """

        return fs.is_file(self.map_earth_path)

    # -----------------------------------------------------------------

    def remove_map_earth(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_earth_path)

    # -----------------------------------------------------------------

    def write_map_earth(self):

        """
        Thisf unction ...
        :return:
        """

        # Write
        self.map_earth.saveto(self.map_earth_path)

    # -----------------------------------------------------------------

    @property
    def map_earth_diffuse_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "earth_diffuse.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_earth_diffuse(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.map_earth_diffuse_path)

    # -----------------------------------------------------------------

    def write_map_earth_diffuse(self):

        """
        Thisfunction ...
        :return:
        """

        # Save
        self.map_earth_diffuse.saveto(self.map_earth_diffuse_path)

    # -----------------------------------------------------------------

    @property
    def map_faceon_path(self):

        """
        Thisk function ...
        :return:
        """

        return fs.join(self.maps_path, "faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.map_faceon_path)

    # -----------------------------------------------------------------

    def remove_map_faceon(self):

        """
        Thisfunction ...
        :return:
        """

        fs.remove_file(self.map_faceon_path)

    # -----------------------------------------------------------------

    def write_map_faceon(self):

        """
        This function ...
        :return:
        """

        self.map_faceon.saveto(self.map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def map_faceon_diffuse_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "faceon_diffuse.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_faceon_diffuse(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.map_faceon_diffuse_path)

    # -----------------------------------------------------------------

    def write_map_faceon_diffuse(self):

        """
        This function ...
        :return:
        """

        # Save
        self.map_faceon_diffuse.saveto(self.map_faceon_diffuse_path)

    # -----------------------------------------------------------------

    @property
    def map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_edgeon(self):

        """
        Thisnfunction ...
        :return:
        """

        return fs.is_file(self.map_edgeon_path)

    # -----------------------------------------------------------------

    def remove_map_edgeon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_edgeon_path)

    # -----------------------------------------------------------------

    def write_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Save
        self.map_edgeon.saveto(self.map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def map_edgeon_diffuse_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "edgeon_diffuse.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_edgeon_diffuse(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.map_edgeon_diffuse_path)

    # -----------------------------------------------------------------

    def write_map_edgeon_diffuse(self):

        """
        Thisfunction ...
        :return:
        """

        self.map_edgeon_diffuse.saveto(self.map_edgeon_diffuse_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Maps of the absorbed energy
        self.plot_absorption_maps()

        # Maps of the total heating fraction
        self.plot_maps()

    # -----------------------------------------------------------------

    def plot_absorption_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the total absorbed energy ...")

        # Earth
        if self.do_earth: self.plot_absorption_maps_earth()

        # Faceon
        if self.do_faceon: self.plot_absorption_maps_faceon()

        # Edgeon
        if self.do_edgeon: self.plot_absorption_maps_edgeon()

    # -----------------------------------------------------------------

    def plot_absorption_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the total absorbed energy from the earth projection ...")

        # Total
        self.plot_absorption_map_earth_total()

        # Young
        self.plot_absorption_map_earth_young()

        # Ionizing
        self.plot_absorption_map_earth_ionizing()

        # Internal
        self.plot_absorption_map_earth_internal()

        # Unevolved
        self.plot_absorption_map_earth_unevolved()

    # -----------------------------------------------------------------

    def plot_absorption_map_earth_total(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_earth_young(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_earth_ionizing(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_earth_internal(self):

        """
        Thisfunction ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_earth_unevolved(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the total absorbed energy from the face-on projection ...")

        # Total
        self.plot_absorption_map_faceon_total()

        # Young
        self.plot_absorption_map_faceon_young()

        # Ionizing
        self.plot_absorption_map_faceon_ionizing()

        # Internal
        self.plot_absorption_map_faceon_internal()

        # Unevolved
        self.plot_absorption_map_faceon_unevolved()

    # -----------------------------------------------------------------

    def plot_absorption_map_faceon_total(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_faceon_young(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_faceon_ionizing(self):

        """
        Thisf unction ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_faceon_internal(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_faceon_unevolved(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the total absorbed energy from the edge-on projection ...")

        # Total
        self.plot_absorption_map_edgeon_total()

        # Young
        self.plot_absorption_map_edgeon_young()

        # Ionizing
        self.plot_absorption_map_edgeon_ionizing()

        # Internal
        self.plot_absorption_map_edgeon_internal()

        # Unevolved
        self.plot_absorption_map_edgeon_unevolved()

    # -----------------------------------------------------------------

    def plot_absorption_map_edgeon_total(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_edgeon_young(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_edgeon_ionizing(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_edgeon_internal(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_absorption_map_edgeon_unevolved(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:
            self.plot_map_earth()
            self.plot_map_earth_diffuse()

        # Face-on
        if self.do_faceon:
            self.plot_map_faceon()
            self.plot_map_faceon_diffuse()

        # Edge-on
        if self.do_edgeon:
            self.plot_map_edgeon()
            self.plot_map_edgeon_diffuse()

    # -----------------------------------------------------------------

    def plot_map_earth(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_map_earth_diffuse(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_map_faceon(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_map_faceon_diffuse(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_map_edgeon(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_map_edgeon_diffuse(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------
