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

        # The dust absorbed energy maps
        self.total_absorption = None
        self.total_absorption_faceon = None
        self.young_absorption = None
        self.young_absorption_faceon = None
        self.ionizing_absorption = None
        self.ionizing_absorption_faceon = None
        self.internal_absorption = None
        self.internal_absorption_faceon = None

        # Maps of heating
        self.map = None
        self.map_faceon = None

        # Cubes of spectral heating
        self.cube = None
        self.cube_faceon = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the absorption maps
        self.get_absorptions()

        # Get the maps of the heating fraction (earth and faceon)
        self.get_maps()

        # Get the cube of the heating fraction per wavelength (earth and faceon)
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
    def young_simulations(self):

        """
        This function ...
        :return:
        """

        return self.model.young_simulations

    # -----------------------------------------------------------------

    @property
    def ionizing_simulations(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_simulations

    # -----------------------------------------------------------------

    def get_absorptions(self):

        """
        This function ...
        :return:
        """

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

        # Earth
        self.get_total_absorption()

        # Face-on
        self.get_total_absorption_faceon()

    # -----------------------------------------------------------------

    def get_total_absorption(self):
        
        """
        This function ...
        :return: 
        """

        # Get total absorption map
        if self.has_total_absorption: self.load_total_absorption()

        # Calculate
        else: self.calculate_total_absorption()

    # -----------------------------------------------------------------

    def load_total_absorption(self):

        """
        Thisf unction ...
        :return:
        """

        # Load
        self.total_absorption = Frame.from_file(self.total_absorption)

    # -----------------------------------------------------------------

    def calculate_total_absorption(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.total_absorption = self.total_simulations.observed_dust_frame

    # -----------------------------------------------------------------

    def get_total_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_total_absorption_faceon: self.load_total_absorptions_faceon()

        # Calculate
        else: self.calculate_total_absorptions_faceon()

    # -----------------------------------------------------------------

    def load_total_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        self.total_absorption_faceon = Frame.from_file(self.total_absorption_faceon_path)

    # -----------------------------------------------------------------

    def calculate_total_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.total_absorption_faceon = self.total_simulations.faceon_observed_dust_sed

    # -----------------------------------------------------------------

    def get_young_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.get_young_absorption()

        # Face-on
        self.get_young_absorption_faceon()

    # -----------------------------------------------------------------

    def get_young_absorption(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_young_absorption: self.load_young_absorption()

        # Calculate
        else: self.calculate_young_absorption()

    # -----------------------------------------------------------------

    def load_young_absorption(self):

        """
        This function ...
        :return:
        """

        # Load
        self.young_absorption = Frame.from_file(self.young_absorption_path)

    # -----------------------------------------------------------------

    def calculate_young_absorption(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.young_absorption = self.young_simulations.observed_dust_frame

    # -----------------------------------------------------------------

    def get_young_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_young_absorption_faceon: self.load_young_absorption_faceon()

        # Calculate
        else: self.calculate_young_absorption_faceon()

    # -----------------------------------------------------------------

    def load_young_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        self.young_absorption_faceon = Frame.from_file(self.young_absorption_faceon_path)

    # -----------------------------------------------------------------

    def calculate_young_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.young_absorption_faceon = self.young_simulations.faceon_observed_dust_frame

    # -----------------------------------------------------------------

    def get_ionizing_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.get_ionizing_absorption()

        # Face-on
        self.get_ionizing_absorption_faceon()

    # -----------------------------------------------------------------

    def get_ionizing_absorption(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_ionizing_absorption: self.load_ionizing_absorption()

        # Calculate
        else: self.calculate_ionizing_absorption()

    # -----------------------------------------------------------------

    def load_ionizing_absorption(self):

        """
        This function ...
        :return:
        """

        # Load
        self.ionizing_absorption = Frame.from_file(self.ionizing_absorption_path)

    # -----------------------------------------------------------------

    def calculate_ionizing_absorption(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.ionizing_absorption = self.ionizing_simulations.observed_dust_frame

    # -----------------------------------------------------------------

    def get_ionizing_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_ionizing_absorption_faceon: self.load_ionizing_absorption_faceon()

        # Calculate
        else: self.calculate_ionizing_absorption_faceon()

    # -----------------------------------------------------------------

    def load_ionizing_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        self.ionizing_absorption_faceon = Frame.from_file(self.ionizing_absorption_faceon_path)

    # -----------------------------------------------------------------

    def calculate_ionizing_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        # Integrates the dust spectral cube
        self.ionizing_absorption_faceon = self.ionizing_simulations.faceon_observed_dust_frame

    # -----------------------------------------------------------------

    def get_internal_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.get_internal_absorption()

        # Faceon
        self.get_internal_absorption_faceon()

    # -----------------------------------------------------------------

    def get_internal_absorption(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_internal_absorption: self.load_internal_absorption()

        # Calculate
        else: self.calculate_internal_absorption()

    # -----------------------------------------------------------------

    def load_internal_absorption(self):

        """
        This function ...
        :return:
        """

        # Load
        self.internal_absorption = Frame.from_file(self.internal_absorption_path)

    # -----------------------------------------------------------------

    def calculate_internal_absorption(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get_internal_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        if self.has_internal_absorption_faceon: self.load_internal_absorption_faceon()

        # Calculate
        else: self.calculate_internal_absorption_faceon()

    # -----------------------------------------------------------------

    def load_internal_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        self.internal_absorption_faceon = Frame.from_file(self.internal_absorption_faceon_path)

    # -----------------------------------------------------------------

    def calculate_internal_absorption_faceon(self):

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
        self.get_earth_map()

        # Face-on
        self.get_faceon_map()

    # -----------------------------------------------------------------

    def get_earth_map(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def get_faceon_map(self):

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

        # Write
        self.write_total_absorption()

        # Write
        self.write_total_absorption_faceon()

    # -----------------------------------------------------------------

    @property
    def total_absorption_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "total.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_absorption(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.total_absorption_path)

    # -----------------------------------------------------------------

    def remove_total_absorption(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.total_absorption_path)

    # -----------------------------------------------------------------

    def write_total_absorption(self):

        """
        This function ...
        :return:
        """

        self.total_absorption.saveto(self.total_absorption_path)

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
    def young_absorption_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "young.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_absorption(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_absorption_path)

    # -----------------------------------------------------------------

    def remove_young_absorption(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.young_absorption_path)

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
