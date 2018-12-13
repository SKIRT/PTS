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
from ....core.tools.utils import lazyproperty, lazyfileproperty
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
        super(ProjectedDustHeatingAnalyser, self).setup()

        # Check
        if (not self.do_earth) and (not self.do_faceon) and (not self.do_edgeon): raise RuntimeError("Cannot create any map: not enough simulation data is present")

    # -----------------------------------------------------------------

    @property
    def do_earth(self):
        return self.has_earth_cube_all

    # -----------------------------------------------------------------

    @lazyproperty
    def do_faceon(self):
        return self.has_faceon_cube_all

    # -----------------------------------------------------------------

    @lazyproperty
    def do_edgeon(self):
        return self.has_edgeon_cube_all

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
    # TOTAL ABSORPTION MAPS
    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "total_absorptions_earth_path", True, write=False)
    def total_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        return self.model.total_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "total_absorptions_faceon_path", True, write=False)
    def total_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.total_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "total_absorptions_edgeon_path", True, write=False)
    def total_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.total_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------
    # YOUNG ABSORPTION MAPS
    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "young_absorptions_earth_path", True, write=False)
    def young_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        return self.model.young_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "young_absorptions_faceon_path", True, write=False)
    def young_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.young_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "young_absorptions_edgeon_path", True, write=False)
    def young_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.young_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------
    # IONIZING ABSORPTION MAPS
    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "ionizing_absorptions_earth_path", True, write=False)
    def ionizing_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "ionizing_absorptions_faceon_path", True, write=False)
    def ionizing_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "ionizing_absorptions_edgeon_path", True, write=False)
    def ionizing_absorptions_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        return self.model.sfr_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------
    # INTERNAL ABSORPTION MAPS
    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "internal_absorptions_earth_path", True, write=False)
    def internal_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_intrinsic_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "internal_absorptions_faceon_path", True, write=False)
    def internal_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_intrinsic_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "internal_absorptions_edgeon_path", True, write=False)
    def internal_absorptions_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        return self.model.sfr_intrinsic_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------
    # UNEVOLVED ABSORPTION MAPS
    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        young, ionizing = uniformize(self.young_absorptions_earth, self.ionizing_absorptions_earth, distance=self.galaxy_distance)
        return young + ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_earth_diffuse(self):

        """
        This function ...
        :return:
        """

        unevolved, internal = uniformize(self.unevolved_absorptions_earth, self.internal_absorptions_earth, distance=self.galaxy_distance)
        return unevolved - internal

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        young, ionizing = uniformize(self.young_absorptions_faceon, self.ionizing_absorptions_faceon, distance=self.galaxy_distance)
        return young + ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_faceon_diffuse(self):

        """
        Thisn function ...
        :return:
        """

        unevolved, internal = uniformize(self.unevolved_absorptions_faceon, self.internal_absorptions_faceon)
        return unevolved - internal

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        young, ionizing = uniformize(self.young_absorptions_edgeon, self.ionizing_absorptions_edgeon)
        return young + ionizing

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorptions_edgeon_diffuse(self):

        """
        Thisf unction ...
        :return:
        """

        unevolved, internal = uniformize(self.unevolved_absorptions_edgeon, self.internal_absorptions_edgeon)
        return unevolved - internal

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorptions_edgeon_diffuse(self):

        """
        This function ...
        :return:
        """

        total, internal = uniformize(self.total_absorptions_edgeon, self.internal_absorptions_edgeon)
        return total - internal

    # -----------------------------------------------------------------
    # TOTAL DIFFUSE ABSORPTION MAPS
    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorptions_earth_diffuse(self):

        """
        This function ...
        :return:
        """

        total, internal = uniformize(self.total_absorptions_earth, self.internal_absorptions_earth, distance=self.galaxy_distance)
        return total - internal

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorptions_faceon_diffuse(self):

        """
        This function ...
        :return:
        """

        # return self.total_absorptions_faceon - self.internal_absorptions_faceon
        total, internal = uniformize(self.total_absorptions_faceon, self.internal_absorptions_faceon)
        return total - internal

    # -----------------------------------------------------------------
    # HEATING FRACTION MAPS
    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "map_earth_path", True, write=False)
    def map_earth(self):

        """
        Thisf unction ...
        :return:
        """

        return self.unevolved_absorptions_earth / self.total_absorptions_earth

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "map_earth_diffuse_path", True, write=False)
    def map_earth_diffuse(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_absorptions_earth_diffuse / self.total_absorptions_earth_diffuse

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "map_faceon_path", True, write=False)
    def map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_absorptions_faceon / self.total_absorptions_faceon

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "map_faceon_diffuse_path", True, write=False)
    def map_faceon_diffuse(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_absorptions_faceon_diffuse / self.total_absorptions_faceon_diffuse

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "map_edgeon_path", True, write=False)
    def map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_absorptions_edgeon / self.total_absorptions_edgeon

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "map_edgeon_diffuse_path", True, write=False)
    def map_edgeon_diffuse(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_absorptions_edgeon_diffuse / self.total_absorptions_edgeon_diffuse

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write absorption luminosity maps
        self.write_absorptions()

        # Write heating fraction maps
        self.write_maps()

    # -----------------------------------------------------------------

    @lazyproperty
    def absorptions_path(self):
        return fs.create_directory_in(self.projected_heating_path, "absorptions")

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_path(self):
        return fs.create_directory_in(self.projected_heating_path, "maps")

    # -----------------------------------------------------------------

    def write_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the absorption maps ...")

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
        return fs.join(self.absorptions_path, "total_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_absorptions_earth(self):
        return fs.is_file(self.total_absorptions_earth_path)

    # -----------------------------------------------------------------

    def remove_total_absorptions_earth(self):
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
        return fs.join(self.absorptions_path, "total_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_absorptions_faceon(self):
        return fs.is_file(self.total_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def remove_total_absorptions_faceon(self):
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
        return fs.join(self.absorptions_path, "total_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_absorptions_edgeon(self):
        return fs.is_file(self.total_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def remove_total_absorptions_edgeon(self):
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
        return fs.join(self.absorptions_path, "young_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_absorptions_earth(self):
        return fs.is_file(self.young_absorptions_earth_path)

    # -----------------------------------------------------------------

    def remove_young_absorptions_earth(self):
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
        return fs.join(self.absorptions_path, "young_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_absorptions_faceon(self):
        return fs.is_file(self.young_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def remove_young_absorptions_faceon(self):
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
        return fs.join(self.absorptions_path, "young_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_absorptions_edgeon(self):
        return fs.is_file(self.young_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def remove_young_absorptions_edgeon(self):
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
        return fs.join(self.absorptions_path, "ionizing_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_ionizing_absorptions_earth(self):
        return fs.is_file(self.ionizing_absorptions_earth_path)

    # -----------------------------------------------------------------

    def remove_ionizing_absorptions_earth(self):
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
        return fs.join(self.absorptions_path, "ionizing_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_ionizing_absorptions_faceon(self):
        return fs.is_file(self.ionizing_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def remove_ionizing_absorptions_faceon(self):
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
        return fs.join(self.absorptions_path, "ionizing_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_ionizing_absorptions_edgeon(self):
        return fs.is_file(self.ionizing_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def remove_ionizing_absorptions_edgeon(self):
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
        return fs.join(self.absorptions_path, "internal_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_internal_absorptions_earth(self):
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
        return fs.join(self.absorptions_path, "internal_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_internal_absorptions_faceon(self):
        return fs.is_file(self.internal_absorptions_faceon_path)

    # -----------------------------------------------------------------

    def remove_internal_absorptions_faceon(self):
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
        return fs.join(self.absorptions_path, "internal_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_internal_absorptions_edgeon(self):
        return fs.is_file(self.internal_absorptions_edgeon_path)

    # -----------------------------------------------------------------

    def remove_internal_absorptions_edgeon(self):
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
        if self.do_earth: self.write_maps_earth()

        # Face-on
        if self.do_faceon: self.write_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.write_maps_edgeon()

    # -----------------------------------------------------------------

    def write_maps_earth(self):

        """
        This function ...
        :param self:
        :return:
        """

        self.write_map_earth()
        self.write_map_earth_diffuse()

    # -----------------------------------------------------------------

    @property
    def map_earth_path(self):
        return fs.join(self.maps_path, "earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_earth(self):
        return fs.is_file(self.map_earth_path)

    # -----------------------------------------------------------------

    def remove_map_earth(self):
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
        return fs.join(self.maps_path, "earth_diffuse.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_earth_diffuse(self):
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

    def write_maps_faceon(self):

        """
        This function ...
        :param self:
        :return:
        """

        self.write_map_faceon()
        self.write_map_faceon_diffuse()

    # -----------------------------------------------------------------

    @property
    def map_faceon_path(self):
        return fs.join(self.maps_path, "faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_faceon(self):
        return fs.is_file(self.map_faceon_path)

    # -----------------------------------------------------------------

    def remove_map_faceon(self):
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
        return fs.join(self.maps_path, "faceon_diffuse.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_faceon_diffuse(self):
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

    def write_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        self.write_map_edgeon()
        self.write_map_edgeon_diffuse()

    # -----------------------------------------------------------------

    @property
    def map_edgeon_path(self):
        return fs.join(self.maps_path, "edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_edgeon(self):
        return fs.is_file(self.map_edgeon_path)

    # -----------------------------------------------------------------

    def remove_map_edgeon(self):
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
        return fs.join(self.maps_path, "edgeon_diffuse.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_edgeon_diffuse(self):
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

    @property
    def absorption_map_names(self):
        return ["total", "total_diffuse", "young", "ionizing", "internal", "unevolved", "unevolved_diffuse"]

    # -----------------------------------------------------------------

    def plot_absorption_map(self, frame, path):

        """
        This function ...
        :param frame:
        :param path:
        :return:
        """

        plotting.plot_frame(frame, path=path, colorbar=True)

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

    @property
    def absorption_maps_earth(self):
        return [self.total_absorptions_earth, self.total_absorptions_earth_diffuse, self.young_absorptions_earth, self.ionizing_absorptions_earth, self.internal_absorptions_earth, self.unevolved_absorptions_earth, self.unevolved_absorptions_earth_diffuse]

    # -----------------------------------------------------------------

    @property
    def absorption_map_paths_earth(self):
        return [fs.join(self.absorptions_path, name + "_earth.pdf") for name in self.absorption_map_names]

    # -----------------------------------------------------------------

    def plot_absorption_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the total absorbed energy from the earth projection ...")

        # Loop
        for frame, path in zip(self.absorption_maps_earth, self.absorption_map_paths_earth):

            # Check
            if fs.is_file(path): continue

            # Plot
            self.plot_absorption_map(frame, path)

    # -----------------------------------------------------------------

    @property
    def absorption_maps_faceon(self):
        return [self.total_absorptions_faceon, self.total_absorptions_faceon_diffuse, self.young_absorptions_faceon, self.ionizing_absorptions_faceon, self.internal_absorptions_faceon, self.unevolved_absorptions_faceon, self.unevolved_absorptions_faceon_diffuse]

    # -----------------------------------------------------------------

    @property
    def absorption_map_paths_faceon(self):
        return [fs.join(self.absorptions_path, name + "_faceon.pdf") for name in self.absorption_map_names]

    # -----------------------------------------------------------------

    def plot_absorption_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the total absorbed energy from the face-on projection ...")

        # Loop
        for frame, path in zip(self.absorption_maps_faceon, self.absorption_map_paths_faceon):

            # Check
            if fs.is_file(path): continue

            # Plot
            self.plot_absorption_map(frame, path)

    # -----------------------------------------------------------------

    @property
    def absorption_maps_edgeon(self):
        return [self.total_absorptions_edgeon, self.total_absorptions_edgeon_diffuse, self.young_absorptions_edgeon, self.ionizing_absorptions_edgeon, self.internal_absorptions_edgeon, self.unevolved_absorptions_edgeon, self.unevolved_absorptions_edgeon_diffuse]

    # -----------------------------------------------------------------

    @property
    def absorption_map_paths_edgeon(self):
        return [fs.join(self.absorptions_path, name + "_edgeon.pdf") for name in self.absorption_map_names]

    # -----------------------------------------------------------------

    def plot_absorption_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the total absorbed energy from the edge-on projection ...")

        # Loop
        for frame, path in zip(self.absorption_maps_edgeon, self.absorption_map_paths_edgeon):

            # Check
            if fs.is_file(path): continue

            # Plot
            self.plot_absorption_map(frame, path)

    # -----------------------------------------------------------------

    @property
    def heating_fraction_interval(self):
        return (-1,1,)

    # -----------------------------------------------------------------

    def plot_heating_fraction_map(self, frame, path):

        """
        This function ...
        :param frame:
        :param path:
        :return:
        """

        plotting.plot_map(frame, path=path, interval=self.heating_fraction_interval, colorbar=True)

    # -----------------------------------------------------------------

    @property
    def heating_map_names(self):
        return ["earth", "earth_diffuse", "faceon", "faceon_diffuse", "edgeon", "edgeon_diffuse"]

    # -----------------------------------------------------------------

    @property
    def heating_maps(self):
        return [self.map_earth, self.map_earth_diffuse, self.map_faceon, self.map_faceon_diffuse, self.map_edgeon, self.map_edgeon_diffuse]

    # -----------------------------------------------------------------

    @property
    def heating_map_paths(self):
        return [fs.join(self.maps_path, name + ".fits") for name in self.heating_map_names]

    # -----------------------------------------------------------------

    def plot_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the heating fraction maps ...")

        # Loop over the maps
        for frame, path in zip(self.heating_maps, self.heating_map_paths):

            # Check
            if fs.is_file(path): continue

            # Plot
            self.plot_heating_fraction_map(frame, path)

# -----------------------------------------------------------------
