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

        # Maps of heating fraction
        self.map_earth = None
        self.map_faceon = None
        self.map_edgeon = None

        # Maps of heating fraction (diffuse)
        self.map_earth_diffuse = None
        self.map_faceon_diffuse = None
        self.map_edgeon_diffuse = None

        # Cubes of spectral heating (looking at dust emission)
        self.cube_earth = None
        self.cube_faceon = None
        self.cube_edgeon = None

        # Cubes of spectral heating (looking at dust absorption)
        self.cube_earth_absorption = None
        self.cube_faceon_absorption = None
        self.cube_edgeon_absorption = None

        # Curves of spectral heating
        self.curve_earth = None
        self.curve_earth_absorption = None
        self.curve_faceon = None
        self.curve_faceon_absorption = None
        self.curve_edgeon = None
        self.curve_edgeon_absorption = None

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

        # Fix the cubes
        self.fix_cubes()

        # Get the curves of the heating fraction per wavelength
        self.get_curves()

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

        return self.has_earth_cube_all
        #return self.has_earth_cube_contributions_all

    # -----------------------------------------------------------------

    @lazyproperty
    def do_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_faceon_cube_all
        #return self.has_faceon_cube_contributions_all

    # -----------------------------------------------------------------

    @lazyproperty
    def do_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_edgeon_cube_all
        #return self.has_edgeon_cube_contributions_all

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

    def fix_cube_emission(self, cube):

        """
        Thisfunction ...
        :param cube:
        :return:
        """

        # Truncate
        cube.truncate(min_wavelength=min_wavelength_emission)

        # Fix
        cube.replace_negatives_by_nans()
        cube.replace_infs_by_nans()

        # Replace
        cube.replace_by_nans_where_greater_than(1.1)
        cube.cutoff_greater(1.)

        # Interpolate nans
        cube.interpolate_nans(sigma=3.)

        # Set flag
        cube.metadata["fixed"] = True

    # -----------------------------------------------------------------

    def fix_cube_absorption(self, cube):

        """
        This function ...
        :param cube:
        :return:
        """

        # Truncate
        cube.truncate(max_wavelength=max_wavelength_absorption)

        # Fix
        cube.replace_negatives_by_nans()
        cube.replace_infs_by_nans()

        # Replace
        cube.replace_by_nans_where_greater_than(1.1)
        cube.cutoff_greater(1.)

        # Interpolate nans
        cube.interpolate_nans(sigma=3.)

        # Set flag
        cube.metadata["fixed"] = True

    # -----------------------------------------------------------------

    def get_cubes(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.get_cubes_earth()

        # Faceon
        if self.do_faceon: self.get_cubes_faceon()

        # Edgeon
        if self.do_edgeon: self.get_cubes_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_cubes_earth_emission(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_emission_cubes_earth

    # -----------------------------------------------------------------

    @property
    def do_cubes_earth_absorption(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_absorption_cubes_earth

    # -----------------------------------------------------------------

    def get_cubes_earth(self):

        """
        This function ...
        :return:
        """

        # Dust emission
        if self.do_cubes_earth_emission: self.get_cube_earth()

        # Dust absorption
        if self.do_cubes_earth_absorption: self.get_cube_earth_absorption()

    # -----------------------------------------------------------------

    @property
    def do_cubes_faceon_emission(self):

        """
        Thisfunction ...
        :return:
        """

        return self.has_dust_emission_cubes_faceon

    # -----------------------------------------------------------------

    @property
    def do_cubes_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        return self.has_dust_absorption_cubes_faceon

    # -----------------------------------------------------------------

    def get_cubes_faceon(self):

        """
        This function ...
        :return: 
        """

        # Dust emission
        if self.do_cubes_faceon_emission: self.get_cube_faceon()

        # Dust absorption
        if self.do_cubes_faceon_absorption: self.get_cube_faceon_absorption()

    # -----------------------------------------------------------------

    @property
    def do_cubes_edgeon_emission(self):

        """
        Thisfunction ...
        :return:
        """

        return self.has_dust_emission_cubes_edgeon

    # -----------------------------------------------------------------

    @property
    def do_cubes_edgeon_absorption(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_dust_absorption_cubes_edgeon

    # -----------------------------------------------------------------

    def get_cubes_edgeon(self):
        
        """
        This function ...
        :return:
        """

        # Dust emission
        if self.do_cubes_edgeon_emission: self.get_cube_edgeon()

        # Dust absorption
        if self.do_cubes_edgeon_absorption: self.get_cube_edgeon_absorption()

    # -----------------------------------------------------------------

    def get_cube_earth(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_cube_earth: self.load_cube_earth()

        # Create
        else: self.create_cube_earth()

    # -----------------------------------------------------------------

    def load_cube_earth(self):

        """
        Thisfunction ...
        :return:
        """

        self.cube_earth = DataCube.from_file(self.cube_earth_path)

    # -----------------------------------------------------------------

    @property
    def old_dust_emission_cube_earth(self):
        return self.model.old_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_cube_earth(self):
        return self.model.has_old_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def young_dust_emission_cube_earth(self):
        return self.model.young_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_cube_earth(self):
        return self.model.has_young_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_cube_earth(self):
        return self.model.sfr_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_cube_earth(self):
        return self.model.has_sfr_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_emission_cube_earth(self):
        return self.young_dust_emission_cube_earth + self.ionizing_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_cube_earth(self):
        return self.has_young_dust_emission_cube_earth and self.has_ionizing_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_cube_earth(self):
        return self.old_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_cube_earth(self):
        return self.has_old_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def total_dust_emission_cube_earth(self):
        return self.model.total_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_cube_earth(self):
        return self.model.has_total_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_cubes_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_unevolved_dust_emission_cube_earth and self.has_total_dust_emission_cube_earth and self.has_evolved_dust_emission_cube_earth

    # -----------------------------------------------------------------

    def create_cube_earth(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN DUST EMISSION WAVELENGTHS
        self.cube_earth = 0.5 * (self.unevolved_dust_emission_cube_earth + (self.total_dust_emission_cube_earth - self.evolved_dust_emission_cube_earth)) / self.total_dust_emission_cube_earth

    # -----------------------------------------------------------------

    def get_cube_earth_absorption(self):

        """
        Thisfunction ...
        :return:
        """

        # LOOKING IN ABSORPTION?

        # Load
        if self.has_cube_earth_absorption: self.load_cube_earth_absorption()

        # Create
        else: self.create_cube_earth_absorption()

    # -----------------------------------------------------------------

    def load_cube_earth_absorption(self):

        """
        This function ...
        :return:
        """

        # Load
        self.cube_earth_absorption = DataCube.from_file(self.cube_earth_absorption_path)

    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_cube_earth(self):
        return self.model.old_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_cube_earth(self):
        return self.model.has_old_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def young_dust_absorption_cube_earth(self):
        return self.model.young_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_young_dust_absorption_cube_earth(self):
        return self.model.has_young_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_cube_earth(self):
        return self.model.sfr_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_cube_earth(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_cube_earth(self):
        #return self.young_dust_absorption_cube_earth + self.ionizing_dust_absorption_cube_earth
        return self.model.unevolved_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_cube_earth(self):
        return self.model.has_unevolved_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def evolved_dust_absorption_cube_earth(self):
        return self.old_dust_absorption_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_absorption_cube_earth(self):
        return self.has_old_dust_absorption_cube_earth

    # -----------------------------------------------------------------

    @property
    def total_dust_absorption_cube_earth(self):
        return self.model.total_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_total_dust_absorption_cube_earth(self):
        return self.model.has_total_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_dust_absorption_cubes_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_unevolved_dust_absorption_cube_earth and self.has_total_dust_absorption_cube_earth and self.has_evolved_dust_absorption_cube_earth

    # -----------------------------------------------------------------

    def create_cube_earth_absorption(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN ABSORBED STELLAR WAVELENGTHS
        self.cube_earth_absorption = 0.5 * (self.unevolved_dust_absorption_cube_earth + (self.total_dust_absorption_cube_earth - self.evolved_dust_absorption_cube_earth)) / self.total_dust_absorption_cube_earth

    # -----------------------------------------------------------------

    def get_cube_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_cube_faceon: self.load_cube_faceon()

        # Create
        else: self.create_cube_faceon()

    # -----------------------------------------------------------------

    def load_cube_faceon(self):

        """
        This function ...
        :return:
        """

        self.cube_faceon = DataCube.from_file(self.cube_faceon_path)

    # -----------------------------------------------------------------

    @property
    def old_dust_emission_cube_faceon(self):
        return self.model.old_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_cube_faceon(self):
        return self.model.has_old_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def young_dust_emission_cube_faceon(self):
        return self.model.young_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_cube_faceon(self):
        return self.model.has_young_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_cube_faceon(self):
        return self.model.sfr_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_cube_faceon(self):
        return self.model.has_sfr_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def unevolved_dust_emission_cube_faceon(self):
        return self.young_dust_emission_cube_faceon + self.ionizing_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_cube_faceon(self):
        return self.has_young_dust_emission_cube_faceon and self.has_ionizing_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_cube_faceon(self):
        return self.old_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_cube_faceon(self):
        return self.has_old_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def total_dust_emission_cube_faceon(self):
        return self.model.total_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_cube_faceon(self):
        return self.model.has_total_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_cubes_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_unevolved_dust_emission_cube_faceon and self.has_total_dust_emission_cube_faceon and self.has_evolved_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    def create_cube_faceon(self):

        """
        THs function ...
        :return:
        """

        # LOOKING IN DUST EMISSION WAVELENGTHS
        self.cube_faceon = 0.5 * (self.unevolved_dust_emission_cube_faceon + (self.total_dust_emission_cube_faceon - self.evolved_dust_emission_cube_faceon)) / self.total_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    def get_cube_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN ABSORPTION

        # Load
        if self.has_cube_faceon_absorption: self.load_cube_faceon_absorption()

        # Create
        else: self.create_cube_faceon_absorption()

    # -----------------------------------------------------------------

    def load_cube_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        self.cube_faceon_absorption = DataCube.from_file(self.cube_faceon_absorption_path)

    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_cube_faceon(self):
        return self.model.old_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_cube_faceon(self):
        return self.model.has_old_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def young_dust_absorption_cube_faceon(self):
        return self.model.young_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_absorption_cube_faceon(self):
        return self.model.has_young_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_cube_faceon(self):
        return self.model.sfr_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_cube_faceon(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_cube_faceon(self):
        return self.young_dust_absorption_cube_faceon + self.ionizing_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_cube_faceon(self):
        return self.has_young_dust_absorption_cube_faceon and self.has_ionizing_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_absorption_cube_faceon(self):
        return self.old_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_absorption_cube_faceon(self):
        return self.has_old_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def total_dust_absorption_cube_faceon(self):
        return self.model.total_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_absorption_cube_faceon(self):
        return self.model.has_total_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_dust_absorption_cubes_faceon(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_unevolved_dust_absorption_cube_faceon and self.has_total_dust_absorption_cube_faceon and self.has_evolved_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    def create_cube_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        self.cube_faceon_absorption = 0.5 * (self.unevolved_dust_absorption_cube_faceon + (self.total_dust_absorption_cube_faceon - self.evolved_dust_absorption_cube_faceon)) / self.total_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    def get_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_cube_edgeon: self.load_cube_edgeon()

        # Create
        else: self.create_cube_edgeon()

    # -----------------------------------------------------------------

    def load_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        self.cube_edgeon = DataCube.from_file(self.cube_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def old_dust_emission_cube_edgeon(self):
        return self.model.old_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_cube_edgeon(self):
        return self.model.has_old_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def young_dust_emission_cube_edgeon(self):
        return self.model.young_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_cube_edgeon(self):
        return self.model.has_young_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_cube_edgeon(self):
        return self.model.sfr_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_cube_edgeon(self):
        return self.model.has_sfr_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_emission_cube_edgeon(self):
        return self.young_dust_emission_cube_edgeon + self.ionizing_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_cube_edgeon(self):
        return self.has_young_dust_emission_cube_edgeon and self.has_ionizing_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_cube_edgeon(self):
        return self.old_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_cube_edgeon(self):
        return self.has_old_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def total_dust_emission_cube_edgeon(self):
        return self.model.total_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_cube_edgeon(self):
        return self.model.has_total_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_cubes_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_unevolved_dust_emission_cube_edgeon and self.has_total_dust_emission_cube_edgeon and self.has_evolved_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    def create_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN DUST EMISSION WAVELENGTHS
        self.cube_edgeon = 0.5 * (self.unevolved_dust_emission_cube_edgeon + (self.total_dust_emission_cube_edgeon - self.evolved_dust_emission_cube_edgeon)) / self.total_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    def get_cube_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN ABSORPTION

        # Load
        if self.has_cube_edgeon_absorption: self.load_cube_edgeon_absorption()

        # Create
        else: self.create_cube_edgeon_absorption()

    # -----------------------------------------------------------------

    def load_cube_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        # Load
        self.cube_edgeon_absorption = DataCube.from_file(self.cube_edgeon_absorption_path)

    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_cube_edgeon(self):
        return self.model.old_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_cube_edgeon(self):
        return self.model.has_old_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def young_dust_absorption_cube_edgeon(self):
        return self.model.young_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_absorption_cube_edgeon(self):
        return self.model.has_young_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_cube_edgeon(self):
        return self.model.sfr_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_cube_edgeon(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_cube_edgeon(self):
        return self.young_dust_absorption_cube_edgeon + self.ionizing_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_cube_edgeon(self):
        return self.has_young_dust_absorption_cube_edgeon and self.has_ionizing_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_absorption_cube_edgeon(self):
        return self.old_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_absorption_cube_edgeon(self):
        return self.has_old_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def total_dust_absorption_cube_edgeon(self):
        return self.model.total_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_absorption_cube_edgeon(self):
        return self.model.has_total_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_dust_absorption_cubes_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_unevolved_dust_absorption_cube_edgeon and self.has_total_dust_absorption_cube_edgeon and self.has_evolved_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    def create_cube_edgeon_absorption(self):

        """
        Thisf unction ...
        :return:
        """

        self.cube_edgeon_absorption = 0.5 * (self.unevolved_dust_absorption_cube_edgeon + (self.total_dust_absorption_cube_edgeon - self.evolved_dust_absorption_cube_edgeon)) / self.total_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def do_fix_cube_earth(self):
        return "fixed" not in self.cube_earth.metadata or not self.cube_earth.metadata["fixed"]

    # -----------------------------------------------------------------

    @property
    def do_fix_cube_earth_absorption(self):
        return "fixed" not in self.cube_earth_absorption.metadata or not self.cube_earth_absorption.metadata["fixed"]

    # -----------------------------------------------------------------

    @property
    def do_fix_cube_faceon(self):
        return "fixed" not in self.cube_faceon.metadata or not self.cube_faceon.metadata["fixed"]

    # -----------------------------------------------------------------

    @property
    def do_fix_cube_faceon_absorption(self):
        return "fixed" not in self.cube_faceon_absorption.metadata or not self.cube_faceon_absorption.metadata["fixed"]

    # -----------------------------------------------------------------

    @property
    def do_fix_cube_edgeon(self):
        return "fixed" not in self.cube_edgeon.metadata or not self.cube_edgeon.metadata["fixed"]

    # -----------------------------------------------------------------

    @property
    def do_fix_cube_edgeon_absorption(self):
        return "fixed" not in self.cube_edgeon_absorption.metadata or not self.cube_edgeon_absorption.metadata["fixed"]

    # -----------------------------------------------------------------

    def fix_cubes(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:
            self.fix_cube_earth()
            self.fix_cube_earth_absorption()

        # Face-on
        if self.do_faceon:
            self.fix_cube_faceon()
            self.fix_cube_faceon_absorption()

        # Edge-on
        if self.do_edgeon:
            self.fix_cube_edgeon()
            self.fix_cube_edgeon_absorption()

    # -----------------------------------------------------------------

    def fix_cube_earth(self):

        """
        This function ...
        :return:
        """

        self.fix_cube_emission(self.cube_earth)

    # -----------------------------------------------------------------

    def fix_cube_earth_absorption(self):

        """
        This function ...
        :return:
        """

        self.fix_cube_absorption(self.cube_earth_absorption)

    # -----------------------------------------------------------------

    def fix_cube_faceon(self):

        """
        This function ...
        :return:
        """

        self.fix_cube_emission(self.cube_faceon)

    # -----------------------------------------------------------------

    def fix_cube_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        self.fix_cube_absorption(self.cube_faceon_absorption)

    # -----------------------------------------------------------------

    def fix_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        self.fix_cube_emission(self.cube_edgeon)

    # -----------------------------------------------------------------

    def fix_cube_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        self.fix_cube_absorption(self.cube_edgeon_absorption)

    # -----------------------------------------------------------------

    def get_curves(self):

        """
        This function ...
        :return:
        """

        if self.do_earth:
            self.get_curve_earth()
            self.get_curve_earth_absorption()

        if self.do_faceon:
            self.get_curve_faceon()
            self.get_curve_faceon_absorption()

        if self.do_edgeon:
            self.get_curve_edgeon()
            self.get_curve_edgeon_absorption()

    # -----------------------------------------------------------------

    def get_curve_earth(self):

        """
        Thisf unction ...
        :return:
        """

        if self.has_curve_earth: self.load_curve_earth()
        else: self.create_curve_earth()

    # -----------------------------------------------------------------

    def load_curve_earth(self):

        """
        This function ...
        :return:
        """

        self.curve_earth = WavelengthCurve.from_file(self.curve_earth_path)

    # -----------------------------------------------------------------

    def create_curve_earth(self):

        """
        This function ...
        :return:
        """

        self.curve_earth = self.cube_earth.global_curve("Funev_emission", measure="mean", description="Fraction of emitted energy by unevolved stars")

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_wavelengths(self):
        return self.curve_earth.wavelengths(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_values(self):
        return self.curve_earth.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_evolved(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values("Fev_emission",
                                                           self.curve_earth_wavelengths,
                                                           1. - self.curve_earth_values,
                                                           description="Fraction of emitted energy by evolved stars")

    # -----------------------------------------------------------------

    def get_curve_earth_absorption(self):

        """
        Thisf unction ...
        :return:
        """

        if self.has_curve_earth_absorption: self.load_curve_earth_absorption()
        else: self.create_curve_earth_absorption()

    # -----------------------------------------------------------------

    def load_curve_earth_absorption(self):

        """
        This function ...
        :return:
        """

        self.curve_earth_absorption = WavelengthCurve.from_file(self.curve_earth_absorption_path)

    # -----------------------------------------------------------------

    def create_curve_earth_absorption(self):

        """
        This function ..
        :return:
        """

        self.curve_earth_absorption = self.cube_earth_absorption.global_curve("Funev_absorption", measure="mean",
                                                       description="Fraction of absorbed energy by unevolved stars")

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_absorption_wavelengths(self):
        return self.curve_earth_absorption.wavelengths(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_absorption_values(self):
        return self.curve_earth_absorption.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_absorption_evolved(self):
        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values("Fev_absorption",
                                                           self.curve_earth_absorption_wavelengths,
                                                           1. - self.curve_earth_absorption_values,
                                                           description="Fraction of absorbed energy by evolved stars")

    # -----------------------------------------------------------------

    def get_curve_faceon(self):

        """
        This function ...
        :return:
        """

        if self.has_curve_faceon: self.load_curve_faceon()
        else: self.create_curve_faceon()

    # -----------------------------------------------------------------

    def load_curve_faceon(self):

        """
        This function ...
        :return:
        """

        self.curve_faceon = WavelengthCurve.from_file(self.curve_faceon_path)

    # -----------------------------------------------------------------

    def create_curve_faceon(self):

        """
        Thisf unction ...
        :return:
        """

        self.curve_faceon = self.cube_faceon.global_curve("Funev_emission_faceon", measure="mean",
                                             description="Fraction of emitted energy by unevolved stars (face-on)")

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_wavelengths(self):
        return self.curve_faceon.wavelengths(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_values(self):
        return self.curve_faceon.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_evolved(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values("Fev_emission_faceon",
                                                           self.curve_faceon_wavelengths,
                                                           1. - self.curve_faceon_values,
                                                           description="Fraction of emitted energy by evolved stars (face-on)")

    # -----------------------------------------------------------------

    def get_curve_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        if self.has_curve_faceon_absorption: self.load_curve_faceon_absorption()
        else: self.create_curve_faceon_absorption()

    # -----------------------------------------------------------------

    def load_curve_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        self.curve_faceon_absorption = WavelengthCurve.from_file(self.curve_faceon_absorption_path)

    # -----------------------------------------------------------------

    def create_curve_faceon_absorption(self):

        """
        Thisf unction ...
        :return:
        """

        self.curve_faceon_absorption = self.cube_faceon_absorption.global_curve("Funev_absorption_faceon", measure="mean",
                                                        description="Fraction of absorbed energy by unevolved stars (face-on)")

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_absorption_wavelengths(self):
        return self.curve_faceon_absorption.wavelengths(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_absorption_values(self):
        return self.curve_faceon_absorption.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_absorption_evolved(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values("Fev_absorption_faceon",
                                                           self.curve_faceon_absorption_wavelengths,
                                                           1. - self.curve_faceon_absorption_values,
                                                           description="Fraction of absorbed energy by evolved stars (face-on)")

    # -----------------------------------------------------------------

    def get_curve_edgeon(self):

        """
        This function ...
        :return:
        """

        if self.has_curve_edgeon: self.load_curve_edgeon()
        else: self.create_curve_edgeon()

    # -----------------------------------------------------------------

    def load_curve_edgeon(self):

        """
        Thisfunction ...
        :return:
        """

        self.curve_edgeon = WavelengthCurve.from_file(self.curve_edgeon_path)

    # -----------------------------------------------------------------

    def create_curve_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        self.curve_edgeon = self.cube_edgeon.global_curve("Funev_emission_edgeon", measure="mean",
                                             description="Fraction of emitted energy by unevolved stars (edge-on)")

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_wavelengths(self):
        return self.curve_edgeon.wavelengths(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_values(self):
        return self.curve_edgeon.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_evolved(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values("Fev_emission_edgeon", self.curve_edgeon_wavelengths,
                                                           1. - self.curve_edgeon_values,
                                                           description="Fraction of emitted energy by evolved stars (edge-on)")

    # -----------------------------------------------------------------

    def get_curve_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        if self.has_curve_edgeon_absorption: self.load_curve_edgeon_absorption()
        else: self.create_curve_edgeon_absorption()

    # -----------------------------------------------------------------

    def load_curve_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        self.curve_edgeon_absorption = WavelengthCurve.from_file(self.curve_edgeon_absorption_path)

    # -----------------------------------------------------------------

    def create_curve_edgeon_absorption(self):

        """
        Thisf unction ...
        :return:
        """

        self.curve_edgeon_absorption = self.cube_edgeon_absorption.global_curve("Funev_absorption_edgeon", measure="mean",
                                                        description="Fraction of absorbed energy by unevolved stars (edge-on)")

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_absorption_wavelengths(self):
        return self.curve_edgeon_absorption.wavelengths(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_absorption_values(self):
        return self.curve_edgeon_absorption.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_absorption_evolved(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values("Fev_absorption_edgeon",
                                                           self.curve_edgeon_absorption_wavelengths,
                                                           1. - self.curve_edgeon_absorption_values,
                                                           description="Fraction of absorbed energy by evolved stars (edge-on)")

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

        # Write cubes
        self.write_cubes()

        # Write the curves
        self.write_curves()

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

    @lazyproperty
    def cubes_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.projected_heating_path, "cubes")

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.projected_heating_path, "curves")

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

    def write_cubes(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:

            self.write_cube_earth()
            self.write_cube_earth_absorption()

        # Face-on
        if self.do_faceon:

            self.write_cube_faceon()
            self.write_cube_faceon_absorption()

        # Edge-on
        if self.do_edgeon:

            self.write_cube_edgeon()
            self.write_cube_edgeon_absorption()

    # -----------------------------------------------------------------

    @property
    def cube_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cubes_path, "earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_earth(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.cube_earth_path)

    # -----------------------------------------------------------------

    def remove_cube_earth(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.cube_earth_path)

    # -----------------------------------------------------------------

    def write_cube_earth(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_earth.saveto(self.cube_earth_path)

    # -----------------------------------------------------------------

    @property
    def cube_earth_absorption_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cubes_path, "earth_absorption.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_earth_absorption(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.is_file(self.cube_earth_absorption_path)

    # -----------------------------------------------------------------

    def remove_cube_earth_absorption(self):

        """
        Thisfunction ...
        :return:
        """

        fs.remove_file(self.cube_earth_absorption_path)

    # -----------------------------------------------------------------

    def write_cube_earth_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_earth_absorption.saveto(self.cube_earth_absorption_path)

    # -----------------------------------------------------------------

    @property
    def cube_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cubes_path, "faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.cube_faceon_path)

    # -----------------------------------------------------------------

    def remove_cube_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.cube_faceon_path)

    # -----------------------------------------------------------------

    def write_cube_faceon(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_faceon.saveto(self.cube_faceon_path)

    # -----------------------------------------------------------------

    @property
    def cube_faceon_absorption_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cubes_path, "faceon_absorption.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.cube_faceon_absorption_path)

    # -----------------------------------------------------------------

    def remove_cube_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.cube_faceon_absorption_path)

    # -----------------------------------------------------------------

    def write_cube_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_faceon_absorption.saveto(self.cube_faceon_absorption_path)

    # -----------------------------------------------------------------

    @property
    def cube_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cubes_path, "edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_edgeon(self):

        """
        THis function ...
        :return:
        """

        return fs.is_file(self.cube_edgeon_path)

    # -----------------------------------------------------------------

    def remove_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.cube_edgeon_path)

    # -----------------------------------------------------------------

    def write_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_edgeon.saveto(self.cube_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def cube_edgeon_absorption_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cubes_path, "edgeon_absorption.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_edgeon_absorption(self):

        """
        This fnuction ...
        :return:
        """

        return fs.is_file(self.cube_edgeon_absorption_path)

    # -----------------------------------------------------------------

    def write_cube_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_edgeon_absorption.saveto(self.cube_edgeon_absorption_path)

    # -----------------------------------------------------------------

    def write_curves(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:
            if not self.has_curve_earth: self.write_curve_earth()
            if not self.has_curve_earth_absorption: self.write_curve_earth_absorption()

        # Face-on
        if self.do_faceon:
            if not self.has_curve_faceon: self.write_curve_faceon()
            if not self.has_curve_faceon_absorption: self.write_curve_faceon_absorption()

        # Edge-on
        if self.do_edgeon:
            if not self.has_curve_edgeon: self.write_curve_edgeon()
            if not self.has_curve_edgeon_absorption: self.write_curve_edgeon_absorption()

    # -----------------------------------------------------------------

    @property
    def curve_earth_path(self):
        return fs.join(self.curves_path, "earth.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth(self):
        return fs.is_file(self.curve_earth_path)

    # -----------------------------------------------------------------

    def write_curve_earth(self):

        """
        This function ...
        :return:
        """

        # Save
        self.curve_earth.saveto(self.curve_earth_path)

    # -----------------------------------------------------------------

    @property
    def curve_earth_absorption_path(self):
        return fs.join(self.curves_path, "earth_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_absorption(self):
        return fs.is_file(self.curve_earth_absorption_path)

    # -----------------------------------------------------------------

    def write_curve_earth_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.curve_earth_absorption.saveto(self.curve_earth_absorption_path)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_path(self):
        return fs.join(self.curves_path, "faceon.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon(self):
        return fs.is_file(self.curve_faceon_path)

    # -----------------------------------------------------------------

    def write_curve_faceon(self):

        """
        This functino ...
        :return:
        """

        # Save
        self.curve_faceon.saveto(self.curve_faceon_path)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_absorption_path(self):
        return fs.join(self.curves_path, "faceon_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_absorption(self):
        return fs.is_file(self.curve_faceon_absorption_path)

    # -----------------------------------------------------------------

    def write_curve_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.curve_faceon_absorption.saveto(self.curve_faceon_absorption_path)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_path(self):
        return fs.join(self.curves_path, "edgeon.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon(self):
        return fs.is_file(self.curve_edgeon_path)

    # -----------------------------------------------------------------

    def write_curve_edgeon(self):

        """
        This function ...
        :return:
        """

        # Save
        self.curve_edgeon.saveto(self.curve_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_absorption_path(self):
        return fs.join(self.curves_path, "edgeon_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_absorption(self):
        return fs.is_file(self.curve_edgeon_absorption_path)

    # -----------------------------------------------------------------

    def write_curve_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.curve_edgeon_absorption.saveto(self.curve_edgeon_absorption_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Maps of spectral heating
        self.plot_spectral_maps()

        # Curves of spectral heating
        self.plot_curves()

    # -----------------------------------------------------------------

    def plot_spectral_maps(self):

        """
        This function ...
        :return:
        """

        if self.do_earth:
            self.plot_spectral_maps_emission_earth()
            self.plot_spectral_maps_absorption_earth()

        if self.do_faceon:
            self.plot_spectral_maps_emission_faceon()
            self.plot_spectral_maps_absorption_faceon()

        if self.do_edgeon:
            self.plot_spectral_maps_emission_edgeon()
            self.plot_spectral_maps_absorption_edgeon()

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_emission_earth(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_earth.frame_for_filter(fltr, convolve=False)

    # -----------------------------------------------------------------

    def get_spectral_map_emission_earth_path(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "earth_emission_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_earth(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_emission_earth_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_emission_earth(self):

        """
        Thsf unction ...
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust emission from the earth projection ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if map exists
            if self.has_spectral_map_emission_earth(fltr): continue

            # Get the map
            frame = self.get_spectral_map_emission_earth(fltr)

            # Get the path
            path = self.get_spectral_map_emission_earth_path(fltr)

            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_absorption_earth(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_earth_absorption.frame_for_filter(fltr, convolve=False)

    # -----------------------------------------------------------------

    def get_spectral_map_absorption_earth_path(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "earth_absorption_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_earth(self, fltr):

        """
        Thisn function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_absorption_earth_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_absorption_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust absorption from the earth projection ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if map exists
            if self.has_spectral_map_absorption_earth(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_earth(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_earth_path(fltr)

            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_emission_faceon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_faceon.frame_for_filter(fltr, convolve=False)

    # -----------------------------------------------------------------

    def get_spectral_map_emission_faceon_path(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "faceon_emission_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_faceon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_emission_faceon_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_emission_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust emission from the face-on projection ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if map exists
            if self.has_spectral_map_emission_faceon(fltr): continue

            # Get the map
            frame = self.get_spectral_map_emission_faceon(fltr)

            # Get the path
            path = self.get_spectral_map_emission_faceon_path(fltr)

            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_absorption_faceon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_faceon_absorption.frame_for_filter(fltr, convolve=False)

    # -----------------------------------------------------------------

    def get_spectral_map_absorption_faceon_path(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "faceon_absorption_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_faceon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_absorption_faceon_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_absorption_faceon(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust absorption from the face-on projection ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if map exists
            if self.has_spectral_map_absorption_faceon(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_faceon(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_faceon_path(fltr)

            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_emission_edgeon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_edgeon.frame_for_filter(fltr, convolve=False)

    # -----------------------------------------------------------------

    def get_spectral_map_emission_edgeon_path(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "edgeon_emission_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_edgeon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_emission_edgeon_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_emission_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust emission from the edge-on projection ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if the map exists
            if self.has_spectral_map_emission_edgeon(fltr): continue
            
            # Get the map
            frame = self.get_spectral_map_emission_edgeon(fltr)
            
            # Get the path
            path = self.get_spectral_map_emission_edgeon_path(fltr)
            
            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_absorption_edgeon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_edgeon_absorption.frame_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_spectral_map_absorption_edgeon_path(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "edgeon_absorption_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_edgeon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_absorption_edgeon_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_absorption_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust absorption from the edge-on projection ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if the map exists
            if self.has_spectral_map_absorption_edgeon(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_edgeon(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_edgeon_path(fltr)

            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    def plot_curves(self):

        """
        This function ...
        :return:
        """

        if self.do_earth:
            if not self.has_curve_earth_emission_plot: self.plot_curve_earth_emission()
            if not self.has_curve_earth_absorption_plot: self.plot_curve_earth_absorption()

        if self.do_faceon:
            if not self.has_curve_faceon_emission_plot: self.plot_curve_faceon_emission()
            if not self.has_curve_faceon_absorption_plot: self.plot_curve_faceon_absorption()

        if self.do_edgeon:
            if not self.has_curve_edgeon_emission_plot: self.plot_curve_edgeon_emission()
            if not self.has_curve_edgeon_absorption_plot: self.plot_curve_edgeon_absorption()

    # -----------------------------------------------------------------

    @property
    def curve_earth_emission_plot_path(self):
        return fs.join(self.curves_path, "earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_emission_plot(self):
        return fs.is_file(self.curve_earth_emission_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_earth_emission(self):

        """
        This function ...
        :return:
        """

        # Plot
        plotting.plot_curve(self.curve_earth, path=self.curve_earth_emission_plot_path)

    # -----------------------------------------------------------------

    @property
    def curve_earth_absorption_plot_path(self):
        return fs.join(self.curves_path, "earth_absorption.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_absorption_plot(self):
        return fs.is_file(self.curve_earth_absorption_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_earth_absorption(self):

        """
        This function ...
        :return:
        """

        # Plot
        plotting.plot_curve(self.curve_earth_absorption, path=self.curve_earth_absorption_plot_path)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_emission_plot_path(self):
        return fs.join(self.curves_path, "faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_emission_plot(self):
        return fs.is_file(self.curve_faceon_emission_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_faceon_emission(self):

        """
        This function ...
        :return:
        """

        # Plot
        plotting.plot_curve(self.curve_faceon, path=self.curve_faceon_emission_plot_path)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_absorption_plot_path(self):
        return fs.join(self.curves_path, "faceon_absorption.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_absorption_plot(self):
        return fs.is_file(self.curve_faceon_absorption_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        # Plot
        plotting.plot_curve(self.curve_faceon_absorption, path=self.curve_faceon_absorption_plot_path)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_emission_plot_path(self):
        return fs.join(self.curves_path, "edgeon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_emission_plot(self):
        return fs.is_file(self.curve_edgeon_emission_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_edgeon_emission(self):

        """
        Thisfunction ...
        :return:
        """

        # Plot
        plotting.plot_curve(self.curve_edgeon, path=self.curve_edgeon_emission_plot_path)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_absorption_plot_path(self):
        return fs.join(self.curves_path, "edgeon_absorption.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_absorption_plot(self):
        return fs.is_file(self.curve_edgeon_absorption_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        # Plot
        plotting.plot_curve(self.curve_edgeon_absorption, path=self.curve_edgeon_absorption_plot_path)

# -----------------------------------------------------------------
