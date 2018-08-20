#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating.spectral Contains the SpectralDustHeatingAnalyser class.

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

class SpectralDustHeatingAnalyser(DustHeatingAnalysisComponent):
    
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
        super(SpectralDustHeatingAnalyser, self).__init__(*args, **kwargs)

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
        super(SpectralDustHeatingAnalyser, self).setup()

    # -----------------------------------------------------------------

    @property
    def total_simulations(self):
        return self.model.total_simulations

    # -----------------------------------------------------------------

    @property
    def young_simulations(self):
        return self.model.young_simulations

    # -----------------------------------------------------------------

    @property
    def ionizing_simulations(self):
        return self.model.sfr_simulations

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
        return self.cube_earth is not None and ("fixed" not in self.cube_earth.metadata or not self.cube_earth.metadata["fixed"])

    # -----------------------------------------------------------------

    @property
    def do_fix_cube_earth_absorption(self):
        return self.cube_earth_absorption is not None and ("fixed" not in self.cube_earth_absorption.metadata or not self.cube_earth_absorption.metadata["fixed"])

    # -----------------------------------------------------------------

    @property
    def do_fix_cube_faceon(self):
        return self.cube_faceon is not None and ("fixed" not in self.cube_faceon.metadata or not self.cube_faceon.metadata["fixed"])

    # -----------------------------------------------------------------

    @property
    def do_fix_cube_faceon_absorption(self):
        return self.cube_faceon_absorption is not None and ("fixed" not in self.cube_faceon_absorption.metadata or not self.cube_faceon_absorption.metadata["fixed"])

    # -----------------------------------------------------------------

    @property
    def do_fix_cube_edgeon(self):
        return self.cube_edgeon is not None and ("fixed" not in self.cube_edgeon.metadata or not self.cube_edgeon.metadata["fixed"])

    # -----------------------------------------------------------------

    @property
    def do_fix_cube_edgeon_absorption(self):
        return self.cube_edgeon_absorption is not None and ("fixed" not in self.cube_edgeon_absorption.metadata or not self.cube_edgeon_absorption.metadata["fixed"])

    # -----------------------------------------------------------------

    def fix_cubes(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:
            if self.do_fix_cube_earth: self.fix_cube_earth()
            if self.do_fix_cube_earth_absorption: self.fix_cube_earth_absorption()

        # Face-on
        if self.do_faceon:
            if self.do_fix_cube_faceon: self.fix_cube_faceon()
            if self.do_fix_cube_faceon_absorption: self.fix_cube_faceon_absorption()

        # Edge-on
        if self.do_edgeon:
            if self.do_fix_cube_edgeon: self.fix_cube_edgeon()
            if self.do_fix_cube_edgeon_absorption: self.fix_cube_edgeon_absorption()

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

    @property
    def do_curve_earth_emission(self):
        return self.do_cubes_earth_emission

    # -----------------------------------------------------------------

    @property
    def do_curve_earth_absorption(self):
        return self.do_cubes_earth_absorption

    # -----------------------------------------------------------------

    @property
    def do_curve_faceon_emission(self):
        return self.do_cubes_faceon_emission

    # -----------------------------------------------------------------

    @property
    def do_curve_faceon_absorption(self):
        return self.do_cubes_faceon_absorption

    # -----------------------------------------------------------------

    @property
    def do_curve_edgeon_emission(self):
        return self.do_cubes_edgeon_emission

    # -----------------------------------------------------------------

    @property
    def do_curve_edgeon_absorption(self):
        return self.do_cubes_edgeon_absorption

    # -----------------------------------------------------------------

    def get_curves(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:
            if self.do_curve_earth_emission: self.get_curve_earth()
            if self.do_curve_earth_absorption: self.get_curve_earth_absorption()

        # Face-on
        if self.do_faceon:
            if self.do_curve_faceon_emission: self.get_curve_faceon()
            if self.do_curve_faceon_absorption: self.get_curve_faceon_absorption()

        # Edge-on
        if self.do_edgeon:
            if self.do_curve_edgeon_emission: self.get_curve_edgeon()
            if self.do_curve_edgeon_absorption: self.get_curve_edgeon_absorption()

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
        return self.curve_earth.wavelengths(asarray=True, unit="micron")

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
                                                           description="Fraction of emitted energy by evolved stars",
                                                           wavelength_unit="micron")

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
        return self.curve_earth_absorption.wavelengths(asarray=True, unit="micron")

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
                                                           description="Fraction of absorbed energy by evolved stars",
                                                           wavelength_unit="micron")

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
        return self.curve_faceon.wavelengths(asarray=True, unit="micron")

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
                                                           description="Fraction of emitted energy by evolved stars (face-on)",
                                                           wavelength_unit="micron")

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
        return self.curve_faceon_absorption.wavelengths(asarray=True, unit="micron")

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
                                                           description="Fraction of absorbed energy by evolved stars (face-on)",
                                                           wavelength_unit="micron")

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
        return self.curve_edgeon.wavelengths(asarray=True, unit="micron")

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
                                                           description="Fraction of emitted energy by evolved stars (edge-on)",
                                                           wavelength_unit="micron")

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
        return self.curve_edgeon_absorption.wavelengths(asarray=True, unit="micron")

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
                                                           description="Fraction of absorbed energy by evolved stars (edge-on)",
                                                           wavelength_unit="micron")

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

        # Write cubes
        self.write_cubes()

        # Write the curves
        self.write_curves()

    # -----------------------------------------------------------------

    @lazyproperty
    def cubes_path(self):
        return fs.create_directory_in(self.projected_heating_path, "cubes")

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_path(self):
        return fs.create_directory_in(self.projected_heating_path, "curves")

    # -----------------------------------------------------------------

    @property
    def do_write_cube_earth(self):
        return self.do_cubes_earth_emission and not self.has_cube_earth

    # -----------------------------------------------------------------

    @property
    def do_write_cube_earth_absorption(self):
        return self.do_cubes_earth_absorption and not self.has_cube_earth_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_cube_faceon(self):
        return self.do_cubes_faceon_emission and not self.has_cube_faceon

    # -----------------------------------------------------------------

    @property
    def do_write_cube_faceon_absorption(self):
        return self.do_cubes_faceon_absorption and not self.has_cube_faceon_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_cube_edgeon(self):
        return self.do_cubes_edgeon_emission and not self.has_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def do_write_cube_edgeon_absorption(self):
        return self.do_cubes_edgeon_absorption and not self.has_cube_edgeon_absorption

    # -----------------------------------------------------------------

    def write_cubes(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:

            if self.do_write_cube_earth: self.write_cube_earth()
            if self.do_write_cube_earth_absorption: self.write_cube_earth_absorption()

        # Face-on
        if self.do_faceon:

            if self.do_write_cube_faceon: self.write_cube_faceon()
            if self.do_write_cube_faceon_absorption: self.write_cube_faceon_absorption()

        # Edge-on
        if self.do_edgeon:

            if self.do_write_cube_edgeon: self.write_cube_edgeon()
            if self.do_write_cube_edgeon_absorption: self.write_cube_edgeon_absorption()

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

    @property
    def do_write_curve_earth(self):
        return self.do_curve_earth_emission and not self.has_curve_earth

    # -----------------------------------------------------------------

    @property
    def do_write_curve_earth_absorption(self):
        return self.do_curve_earth_absorption and not self.has_curve_earth_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_curve_faceon(self):
        return self.do_curve_faceon_emission and not self.has_curve_faceon

    # -----------------------------------------------------------------

    @property
    def do_write_curve_faceon_absorption(self):
        return self.do_curve_faceon_absorption and not self.has_curve_faceon_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_curve_edgeon(self):
        return self.do_curve_edgeon_emission and not self.has_curve_edgeon

    # -----------------------------------------------------------------

    @property
    def do_write_curve_edgeon_absorption(self):
        return self.do_curve_edgeon_absorption and not self.has_curve_edgeon_absorption

    # -----------------------------------------------------------------

    def write_curves(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:
            if self.do_write_curve_earth: self.write_curve_earth()
            if self.do_write_curve_earth_absorption: self.write_curve_earth_absorption()

        # Face-on
        if self.do_faceon:
            if self.do_write_curve_faceon: self.write_curve_faceon()
            if self.do_write_curve_faceon_absorption: self.write_curve_faceon_absorption()

        # Edge-on
        if self.do_edgeon:
            if self.do_write_curve_edgeon: self.write_curve_edgeon()
            if self.do_write_curve_edgeon_absorption: self.write_curve_edgeon_absorption()

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

        # Maps of the absorbed energy
        self.plot_absorption_maps()

        # Maps of the total heating fraction
        self.plot_maps()

        # Maps of spectral heating
        self.plot_spectral_maps()

        # Curves of spectral heating
        self.plot_curves()

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

    @property
    def do_plot_spectral_maps_emission_earth(self):
        return self.do_cubes_earth_emission

    # -----------------------------------------------------------------

    @property
    def do_plot_spectral_maps_absorption_earth(self):
        return self.do_cubes_earth_absorption

    # -----------------------------------------------------------------

    @property
    def do_plot_spectral_maps_emission_faceon(self):
        return self.do_cubes_faceon_emission

    # -----------------------------------------------------------------

    @property
    def do_plot_spectral_maps_absorption_faceon(self):
        return self.do_cubes_faceon_absorption

    # -----------------------------------------------------------------

    @property
    def do_plot_spectral_maps_emission_edgeon(self):
        return self.do_cubes_edgeon_emission

    # -----------------------------------------------------------------

    @property
    def do_plot_spectral_maps_absorption_edgeon(self):
        return self.do_cubes_edgeon_absorption

    # -----------------------------------------------------------------

    def plot_spectral_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:
            if self.do_plot_spectral_maps_emission_earth: self.plot_spectral_maps_emission_earth()
            if self.do_plot_spectral_maps_absorption_earth: self.plot_spectral_maps_absorption_earth()

        # Face-on
        if self.do_faceon:
            if self.do_plot_spectral_maps_emission_faceon: self.plot_spectral_maps_emission_faceon()
            if self.do_plot_spectral_maps_absorption_faceon: self.plot_spectral_maps_absorption_faceon()

        # Edge-on
        if self.do_edgeon:
            if self.do_plot_spectral_maps_emission_edgeon: self.plot_spectral_maps_emission_edgeon()
            if self.do_plot_spectral_maps_absorption_edgeon: self.plot_spectral_maps_absorption_edgeon()

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

    @property
    def do_plot_curve_earth_emission(self):
        return self.do_curve_earth_emission and not self.has_curve_earth_emission_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_earth_absorption(self):
        return self.do_curve_earth_absorption and not self.has_curve_earth_absorption_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_faceon_emission(self):
        return self.do_curve_faceon_emission and not self.has_curve_faceon_emission_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_faceon_absorption(self):
        return self.do_curve_faceon_absorption and not self.has_curve_faceon_absorption_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_edgeon_emission(self):
        return self.do_curve_edgeon_emission and not self.has_curve_edgeon_emission_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_edgeon_absorption(self):
        return self.do_curve_edgeon_absorption and not self.has_curve_edgeon_absorption_plot

    # -----------------------------------------------------------------

    def plot_curves(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:
            if self.do_plot_curve_earth_emission: self.plot_curve_earth_emission()
            if self.do_plot_curve_earth_absorption: self.plot_curve_earth_absorption()

        # Face-on
        if self.do_faceon:
            if self.do_plot_curve_faceon_emission: self.plot_curve_faceon_emission()
            if self.do_plot_curve_faceon_absorption: self.plot_curve_faceon_absorption()

        # Edge-on
        if self.do_edgeon:
            if self.do_plot_curve_edgeon_emission: self.plot_curve_edgeon_emission()
            if self.do_plot_curve_edgeon_absorption: self.plot_curve_edgeon_absorption()

    # -----------------------------------------------------------------

    @property
    def curve_earth_emission_plot_path(self):
        return fs.join(self.curves_path, "earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_emission_plot(self):
        return fs.is_file(self.curve_earth_emission_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_earth(self):

        """
        This function ...
        :return:
        """

        return {"unevolved": self.curve_earth, "evolved": self.curve_earth_evolved}

    # -----------------------------------------------------------------

    def plot_curve_earth_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust emission from the earth projection ...")

        # Plot
        plotting.plot_curves(self.curves_earth, path=self.curve_earth_emission_plot_path, xlog=True)

    # -----------------------------------------------------------------

    @property
    def curve_earth_absorption_plot_path(self):
        return fs.join(self.curves_path, "earth_absorption.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_absorption_plot(self):
        return fs.is_file(self.curve_earth_absorption_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_earth_absorption(self):

        """
        This function ...
        :return:
        """

        return {"unevolved": self.curve_earth_absorption, "evolved": self.curve_earth_absorption_evolved}

    # -----------------------------------------------------------------

    def plot_curve_earth_absorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust absorption from the earth projection ...")

        # Plot
        plotting.plot_curves(self.curves_earth_absorption, path=self.curve_earth_absorption_plot_path, xlog=True)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_emission_plot_path(self):
        return fs.join(self.curves_path, "faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_emission_plot(self):
        return fs.is_file(self.curve_faceon_emission_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_faceon(self):

        """
        This function ...
        :return:
        """

        return {"unevolved": self.curve_faceon, "evolved": self.curve_faceon_evolved}

    # -----------------------------------------------------------------

    def plot_curve_faceon_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust emission from the face-on projection ...")

        # Plot
        plotting.plot_curves(self.curves_faceon, path=self.curve_faceon_emission_plot_path, xlog=True)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_absorption_plot_path(self):
        return fs.join(self.curves_path, "faceon_absorption.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_absorption_plot(self):
        return fs.is_file(self.curve_faceon_absorption_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_faceon_absorption(self):

        """
        Thisf unction ...
        :return:
        """

        return {"unevolved": self.curve_faceon_absorption, "evolved": self.curve_faceon_absorption_evolved}

    # -----------------------------------------------------------------

    def plot_curve_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust absorption from the face-on projection ...")

        # Plot
        plotting.plot_curves(self.curves_faceon_absorption, path=self.curve_faceon_absorption_plot_path, xlog=True)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_emission_plot_path(self):
        return fs.join(self.curves_path, "edgeon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_emission_plot(self):
        return fs.is_file(self.curve_edgeon_emission_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_edgeon(self):

        """
        This function ...
        :return:
        """

        return {"unevolved": self.curve_edgeon, "evolved": self.curve_edgeon_evolved}

    # -----------------------------------------------------------------

    def plot_curve_edgeon_emission(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust emission from the edge-on projection ...")

        # Plot
        plotting.plot_curves(self.curves_edgeon, path=self.curve_edgeon_emission_plot_path, xlog=True)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_absorption_plot_path(self):
        return fs.join(self.curves_path, "edgeon_absorption.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_absorption_plot(self):
        return fs.is_file(self.curve_edgeon_absorption_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        return {"unevolved": self.curve_edgeon_absorption, "evolved": self.curve_edgeon_absorption_evolved}

    # -----------------------------------------------------------------

    def plot_curve_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust absorption from the edge-on projection ...")

        # Plot
        plotting.plot_curves(self.curves_edgeon_absorption, path=self.curve_edgeon_absorption_plot_path, xlog=True)

# -----------------------------------------------------------------
