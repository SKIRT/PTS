#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.energy.projected Contains the ProjectedEnergyAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component import AnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.tools.utils import lazyproperty
from ....magic.core.frame import Frame

# -----------------------------------------------------------------

class ProjectedEnergyAnalyser(AnalysisComponent):
    
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
        super(ProjectedEnergyAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The analysis run
        self.analysis_run = None

        # The emissions
        self.emission_map_earth = None
        self.emission_map_faceon = None
        self.emission_map_edgeon = None

        # The absorptions
        self.absorption_map_earth = None
        self.absorption_map_faceon = None
        self.absorption_map_edgeon = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get the dust emission maps
        self.get_emissions()

        # Get the dust absorption maps
        self.get_absorptions()

        # Write
        self.write()

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ProjectedEnergyAnalyser, self).setup(**kwargs)

        # Load the run
        self.load_run()

    # -----------------------------------------------------------------

    def load_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the analysis run " + self.config.run + " ...")

        # Get the run
        self.analysis_run = self.get_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def model(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.model

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
    def observed_cube_dust(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.observed_cube_dust

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_frame(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube_dust.integrate()

    # -----------------------------------------------------------------

    @property
    def observed_cube_dust_faceon(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.faceon_observed_cube_dust

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_frame_faceon(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube_dust_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def observed_cube_dust_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.edgeon_observed_cube_dust

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_frame_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube_dust_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def observed_cube_transparent(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.observed_cube_transparent

    # -----------------------------------------------------------------

    @property
    def observed_cube_transparent_faceon(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.faceon_observed_cube_transparent

    # -----------------------------------------------------------------

    @property
    def observed_cube_transparent_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.edgeon_observed_cube_transparent

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_transparent_frame(self):

        """
        This function ....
        :return:
        """

        return self.observed_cube_transparent.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_transparent_frame_faceon(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube_transparent_faceon.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_transparent_frame_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube_transparent_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_frame(self):

        """
        This function ...
        :return:
        """

        return self.observed_transparent_frame

    # -----------------------------------------------------------------

    @property
    def observed_cube(self):

        """
        This function ..
        :return:
        """

        return self.total_simulations.observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_cube(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube - self.observed_cube_dust

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_frame(self):

        """
        This function ...
        :return:
        """

        return self.observed_stellar_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_frame_faceon(self):

        """
        This function ...
        :return:
        """

        return self.observed_transparent_frame_faceon

    # -----------------------------------------------------------------

    @property
    def observed_cube_faceon(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_cube_faceon(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube_faceon - self.observed_cube_dust_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_frame_faceon(self):

        """
        This function ...
        :return:
        """

        return self.observed_stellar_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_frame_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.observed_transparent_frame_edgeon

    # -----------------------------------------------------------------

    @property
    def observed_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.total_simulations.edgeon_observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube_edgeon - self.observed_cube_dust_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_frame_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.observed_stellar_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def absorbed_stellar_frame(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_frame - self.observed_stellar_frame

    # -----------------------------------------------------------------

    @lazyproperty
    def absorbed_stellar_frame_faceon(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_frame_faceon - self.observed_stellar_frame_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def absorbed_stellar_frame_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_frame_edgeon - self.observed_stellar_frame_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def projected_energy_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.energy_path, "projected")

    # -----------------------------------------------------------------

    @lazyproperty
    def emissions_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.projected_energy_path, "emissions")

    # -----------------------------------------------------------------

    @lazyproperty
    def absorptions_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.projected_energy_path, "absorptions")

    # -----------------------------------------------------------------

    def get_emissions(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.get_emissions_earth()

        # Faceon
        self.get_emissions_faceon()

        # Edgeon
        self.get_emissions_edgeon()

    # -----------------------------------------------------------------

    def get_emissions_earth(self):

        """
        This function ...
        :return:
        """

        if self.has_emission_map_earth: self.load_emissions_earth()
        else: self.calculate_emissions_earth()

    # -----------------------------------------------------------------

    def load_emissions_earth(self):

        """
        This function ...
        :return:
        """

        # Load
        self.emission_map_earth = Frame.from_file(self.emission_map_earth_path)

    # -----------------------------------------------------------------

    def calculate_emissions_earth(self):

        """
        This function ...
        :return:
        """

        # Get grame
        self.emission_map_earth = self.observed_dust_frame

    # -----------------------------------------------------------------

    def get_emissions_faceon(self):

        """
        This function ...
        :return:
        """

        if self.has_emission_map_faceon: self.load_emissions_faceon()
        else: self.calculate_emissions_faceon()

    # -----------------------------------------------------------------

    def load_emissions_faceon(self):

        """
        This function ...
        :return:
        """

        # Load
        self.emission_map_faceon = Frame.from_file(self.emission_map_faceon_path)

    # -----------------------------------------------------------------

    def calculate_emissions_faceon(self):

        """
        This function ...
        :return:
        """

        self.emission_map_faceon = self.observed_dust_frame_faceon

    # -----------------------------------------------------------------

    def get_emissions_edgeon(self):

        """
        This function ...
        :return:
        """

        if self.has_emission_map_edgeon: self.load_emissions_edgeon()
        else: self.calculate_emissions_edgeon()

    # -----------------------------------------------------------------

    def load_emissions_edgeon(self):

        """
        This function ...
        :return:
        """

        # Load
        self.emission_map_edgeon = Frame.from_file(self.emission_map_edgeon_path)

    # -----------------------------------------------------------------

    def calculate_emissions_edgeon(self):

        """
        This function ...
        :return:
        """

        self.emission_map_edgeon = self.observed_dust_frame_edgeon

    # -----------------------------------------------------------------

    def get_absorptions(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.get_absorptions_earth()

        # Face-on
        self.get_absorptions_faceon()

        # Edge-on
        self.get_absorptions_edgeon()

    # -----------------------------------------------------------------

    def get_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        if self.has_absorption_map_earth: self.load_absorptions_earth()
        else: self.calculate_absorptions_earth()

    # -----------------------------------------------------------------

    def load_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        self.absorption_map_earth = Frame.from_file(self.absorption_map_earth_path)

    # -----------------------------------------------------------------

    def calculate_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        self.absorption_map_earth = self.absorbed_stellar_frame

    # -----------------------------------------------------------------

    def get_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        if self.has_absorption_map_faceon: self.load_absorptions_faceon()
        else: self.calculate_absorptions_faceon()

    # -----------------------------------------------------------------

    def load_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        self.absorption_map_faceon = Frame.from_file(self.absorption_map_faceon_path)

    # -----------------------------------------------------------------

    def calculate_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        self.absorption_map_faceon = self.absorbed_stellar_frame_faceon

    # -----------------------------------------------------------------

    def get_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        if self.has_absorption_map_edgeon: self.load_absorptions_edgeon()
        else: self.calculate_absorptions_edgeon()

    # -----------------------------------------------------------------

    def load_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        self.absorption_map_edgeon = Frame.from_file(self.absorption_map_edgeon_path)

    # -----------------------------------------------------------------

    def calculate_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        self.absorption_map_edgeon = self.absorbed_stellar_frame_edgeon

    # -----------------------------------------------------------------

    def write(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Emissions
        self.write_emissions()

        # Absorptions
        self.write_absorptions()

    # -----------------------------------------------------------------

    def write_emissions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust emission maps ...")

        # Earth
        if self.do_write_emissions_earth: self.write_emissions_earth()

        # Face-on
        if self.do_write_emissions_faceon: self.write_emissions_faceon()

        # Edge-on
        if self.do_write_emissions_edgeon: self.write_emissions_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_emissions_earth(self):

        """
        This function ...
        :return:
        """

        return not self.has_emission_map_earth and self.emission_map_earth is not None

    # -----------------------------------------------------------------

    @property
    def emission_map_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.emissions_path, "earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_emission_map_earth(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.emission_map_earth_path)

    # -----------------------------------------------------------------

    def remove_emission_map_earth(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.emission_map_earth_path)

    # -----------------------------------------------------------------

    def write_emissions_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust emission map from the earth orientation ...")

        # Save
        self.emission_map_earth.saveto(self.emission_map_earth_path)

    # -----------------------------------------------------------------

    @property
    def do_write_emissions_faceon(self):

        """
        Thisn function ...
        :return:
        """

        return not self.has_emission_map_faceon and self.emission_map_faceon is not None

    # -----------------------------------------------------------------

    @property
    def emission_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.emissions_path, "faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_emission_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.emission_map_faceon_path)

    # -----------------------------------------------------------------

    def remove_emission_map_faceon(self):

        """
        Thisfunction ...
        :return:
        """

        fs.remove_file(self.emission_map_faceon_path)

    # -----------------------------------------------------------------

    def write_emissions_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust emission map from the faceon orientation ...")

        # Save
        self.emission_map_faceon.saveto(self.emission_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_emissions_edgeon(self):

        """
        This function ...
        :return:
        """

        return not self.has_emission_map_edgeon and self.emission_map_edgeon is not None

    # -----------------------------------------------------------------

    @property
    def emission_map_edgeon_path(self):

        """
        Thisnfunction ..
        :return:
        """

        return fs.join(self.emissions_path, "edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_emission_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.emission_map_edgeon_path)

    # -----------------------------------------------------------------

    def remove_emission_map_edgeon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.emission_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_emissions_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust emission map from the edge-on orientation ...")

        # Save
        self.emission_map_edgeon.saveto(self.emission_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust absorption maps ...")

        # Earth
        if self.do_write_absorptions_earth: self.write_absorptions_earth()

        # Face-on
        if self.do_write_absorptions_faceon: self.write_absorptions_faceon()

        # Edge-on
        if self.do_write_absorptions_edgeon: self.write_absorptions_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        return not self.has_absorption_map_earth and self.absorption_map_earth is not None

    # -----------------------------------------------------------------

    @property
    def absorption_map_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_absorption_map_earth(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.absorption_map_earth_path)

    # -----------------------------------------------------------------

    def remove_absorption_map_earth(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.absorption_map_earth_path)

    # -----------------------------------------------------------------

    def write_absorptions_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust absorption map from the earth projection ...")

        # Save
        self.absorption_map_earth.saveto(self.absorption_map_earth_path)

    # -----------------------------------------------------------------

    @property
    def do_write_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        return not self.has_absorption_map_faceon and self.absorption_map_faceon is not None

    # -----------------------------------------------------------------

    @property
    def absorption_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_absorption_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.absorption_map_faceon_path)

    # -----------------------------------------------------------------

    def remove_absorption_map_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.absorption_map_faceon_path)

    # -----------------------------------------------------------------

    def write_absorptions_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust absorption map from the faceon projection ...")

        # Save
        self.absorption_map_faceon.saveto(self.absorption_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        return not self.has_absorption_map_edgeon and self.absorption_map_edgeon is not None

    # -----------------------------------------------------------------

    @property
    def absorption_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.absorptions_path, "edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_absorption_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.absorption_map_edgeon_path)

    # -----------------------------------------------------------------

    def remove_absorption_map_edgeon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.absorption_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_absorptions_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust absorption map from the edgeon projection ...")

        # Save
        self.absorption_map_edgeon.saveto(self.absorption_map_edgeon_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
