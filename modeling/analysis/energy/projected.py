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
from ..component import AnalysisRunComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.tools.utils import lazyproperty
from ....magic.core.frame import Frame, regularize_frame, nan_value
from ....magic.tools.plotting import plot_frame_contours, plot_frame, plot_map
from ....core.tools.utils import lazyfileproperty

# -----------------------------------------------------------------

class ProjectedEnergyAnalyser(AnalysisRunComponent):
    
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

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Write
        self.write()

        # Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ProjectedEnergyAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @property
    def total_simulations(self):
        return self.model.total_simulations

    # -----------------------------------------------------------------

    @property
    def do_earth(self):
        return True

    # -----------------------------------------------------------------

    @lazyproperty
    def do_faceon(self):
        return self.total_simulations.has_faceon_observed_cube_orientation and self.total_simulations.has_other_observed_cube_contributions_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def do_edgeon(self):
        return self.total_simulations.has_edgeon_observed_cube_orientation and self.total_simulations.has_other_observed_cube_contributions_edgeon

    # -----------------------------------------------------------------

    @property
    def observed_cube_scattered(self):
        return self.total_simulations.observed_cube_scattered

    # -----------------------------------------------------------------

    @property
    def observed_cube_scattered_faceon(self):
        return self.total_simulations.faceon_observed_cube_scattered

    # -----------------------------------------------------------------

    @property
    def observed_cube_scattered_edgeon(self):
        return self.total_simulations.edgeon_observed_cube_scattered

    # -----------------------------------------------------------------

    @property
    def observed_cube_dust(self):
        return self.total_simulations.observed_cube_dust

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_frame(self):
        return self.observed_cube_dust.integrate()

    # -----------------------------------------------------------------

    @property
    def observed_cube_dust_faceon(self):
        return self.total_simulations.faceon_observed_cube_dust

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_frame_faceon(self):
        return self.observed_cube_dust_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def observed_cube_dust_edgeon(self):
        return self.total_simulations.edgeon_observed_cube_dust

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_frame_edgeon(self):
        return self.observed_cube_dust_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def observed_cube_transparent(self):
        return self.total_simulations.observed_cube_transparent

    # -----------------------------------------------------------------

    @property
    def observed_cube_transparent_faceon(self):
        return self.total_simulations.faceon_observed_cube_transparent

    # -----------------------------------------------------------------

    @property
    def observed_cube_transparent_edgeon(self):
        return self.total_simulations.edgeon_observed_cube_transparent

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_transparent_frame(self):
        return self.observed_cube_transparent.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_transparent_frame_faceon(self):
        return self.observed_cube_transparent_faceon.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_transparent_frame_edgeon(self):
        return self.observed_cube_transparent_edgeon.integrate()

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_cube(self):
        return self.observed_cube_transparent

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_frame(self):
        return self.observed_transparent_frame

    # -----------------------------------------------------------------

    @property
    def observed_cube(self):
        return self.total_simulations.observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_cube(self):
        return self.observed_cube - self.observed_cube_dust

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_frame(self):
        return self.observed_stellar_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_cube_faceon(self):
        return self.observed_cube_transparent_faceon

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_frame_faceon(self):
        return self.observed_transparent_frame_faceon

    # -----------------------------------------------------------------

    @property
    def observed_cube_faceon(self):
        return self.total_simulations.faceon_observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_cube_faceon(self):
        return self.observed_cube_faceon - self.observed_cube_dust_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_frame_faceon(self):
        return self.observed_stellar_cube_faceon.integrate()

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_cube_edgeon(self):
        return self.observed_cube_transparent_edgeon

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_frame_edgeon(self):
        return self.observed_transparent_frame_edgeon

    # -----------------------------------------------------------------

    @property
    def observed_cube_edgeon(self):
        return self.total_simulations.edgeon_observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_cube_edgeon(self):
        return self.observed_cube_edgeon - self.observed_cube_dust_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_frame_edgeon(self):
        return self.observed_stellar_cube_edgeon.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def absorbed_stellar_cube(self):
        return self.intrinsic_stellar_cube - self.observed_stellar_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def absorbed_stellar_frame(self):
        return self.intrinsic_stellar_frame - self.observed_stellar_frame

    # -----------------------------------------------------------------

    @lazyproperty
    def absorbed_stellar_cube_faceon(self):
        return self.intrinsic_stellar_cube_faceon - self.observed_stellar_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def absorbed_stellar_frame_faceon(self):
        return self.intrinsic_stellar_frame_faceon - self.observed_stellar_frame_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def absorbed_stellar_cube_edgeon(self):
        return self.intrinsic_stellar_cube_edgeon - self.observed_stellar_cube_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def absorbed_stellar_frame_edgeon(self):
        return self.intrinsic_stellar_frame_edgeon - self.observed_stellar_frame_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def projected_energy_path(self):
        return fs.create_directory_in(self.analysis_run.energy_path, "projected")

    # -----------------------------------------------------------------

    @lazyproperty
    def emissions_path(self):
        return fs.create_directory_in(self.projected_energy_path, "emissions")

    # -----------------------------------------------------------------

    @lazyproperty
    def absorptions_path(self):
        return fs.create_directory_in(self.projected_energy_path, "absorptions")

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths_path(self):
        return fs.create_directory_in(self.projected_energy_path, "wavelengths")

    # -----------------------------------------------------------------
    # EMISSION
    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "emission_map_earth_path", True, write=False)
    def emission_map_earth(self):

        """
        Thisf unction ...
        :return:
        """

        return self.observed_dust_frame

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "emission_map_faceon_path", True, write=False)
    def emission_map_faceon(self):

        """
        Thisfunction ...
        :return:
        """

        return self.observed_dust_frame_faceon

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "emission_map_edgeon_path", True, write=False)
    def emission_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.observed_dust_frame_edgeon

    # -----------------------------------------------------------------
    # ABSORPTION
    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "absorption_map_earth_path", True, write=False)
    def absorption_map_earth(self):

        """
        This function ...
        :return:
        """

        # Get
        frame = self.absorbed_stellar_frame

        # Remove negatives
        frame.replace_negatives_by_zeroes()

        # Return
        return frame

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "absorption_map_faceon_path", True, write=False)
    def absorption_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Get
        frame = self.absorbed_stellar_frame_faceon

        # Remove negatives
        frame.replace_negatives_by_zeroes()

        # Return
        return frame

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "absorption_map_edgeon_path", True, write=False)
    def absorption_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Get
        frame = self.absorbed_stellar_frame_edgeon

        # Remove negatives
        frame.replace_negatives_by_zeroes()

        # Return
        return frame

    # -----------------------------------------------------------------
    # ENERGY BALANCE
    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "map_earth_path", True, write=False)
    def map_earth(self):

        """
        This function ...
        :return:
        """

        # Make absolute map
        map_earth = self.emission_map_earth - self.absorption_map_earth

        # Make relative
        map_earth /= self.emission_map_earth

        # Return
        return map_earth

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "map_faceon_path", True, write=False)
    def map_faceon(self):

        """
        This function ...
        :return:
        """

        # Make absolute map
        map_faceon = self.emission_map_faceon - self.absorption_map_faceon

        # Make relative
        map_faceon /= self.emission_map_faceon

        # Return
        return map_faceon

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "map_edgeon_path", True, write=False)
    def map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Make absolute map
        map_edgeon = self.emission_map_edgeon - self.absorption_map_edgeon

        # Make relative
        map_edgeon /= self.emission_map_edgeon

        # Return
        return map_edgeon

    # -----------------------------------------------------------------
    # WAVELENGTH MAPS
    #   1. Scattered
    # -----------------------------------------------------------------

    @lazyproperty
    def mean_scattered_wavelengths_earth(self):
        return self.observed_cube_scattered.mean_wavelengths()

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "scattered_wavelength_map_earth_path", True, write=False)
    def scattered_wavelength_map_earth(self):

        """
        This function ...
        :return:
        """

        scat_wavelengths_earth = self.mean_scattered_wavelengths_earth
        nans = scat_wavelengths_earth.dilate_nans(radius=4, niterations=2)
        return scat_wavelengths_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_scattered_wavelengths_faceon(self):
        return self.observed_cube_scattered_faceon.mean_wavelengths()

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "scattered_wavelength_map_faceon_path", True, write=False)
    def scattered_wavelength_map_faceon(self):

        """
        Thisn function ...
        :return:
        """

        scat_wavelengths_faceon = self.mean_scattered_wavelengths_faceon
        return scat_wavelengths_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_scattered_wavelengths_edgeon(self):
        return self.observed_cube_scattered_edgeon.mean_wavelengths()

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "scattered_wavelength_map_edgeon_path", True, write=False)
    def scattered_wavelength_map_edgeon(self):

        """
        This function ...
        :return:
        """

        scat_wavelengths_edgeon = self.mean_scattered_wavelengths_edgeon
        return scat_wavelengths_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_absorbed_wavelengths_earth(self):
        return self.absorbed_stellar_cube.mean_wavelengths()

    # -----------------------------------------------------------------
    #   2. Absorbed
    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "absorbed_wavelength_map_earth_path", True, write=False)
    def absorbed_wavelength_map_earth(self):

        """
        This function ...
        :return:
        """

        abs_wavelengths_earth = self.mean_absorbed_wavelengths_earth

        negatives = abs_wavelengths_earth.negatives.disk_dilated(radius=3)
        negatives.fill_holes()
        mask = negatives.largest(connectivity=4).inverse()
        mask.disk_dilate(radius=3, niterations=2)
        abs_wavelengths_earth.apply_mask(mask, fill=nan_value)
        abs_wavelengths_earth.replace_negatives_by_nans()
        #self.abs_wavelengths_earth.cutoff_greater(1)
        abs_wavelengths_earth.replace_by_nans_where_greater_than(1)

        return abs_wavelengths_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_absorbed_wavelengths_faceon(self):
        return self.absorbed_stellar_cube_faceon.mean_wavelengths()

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "absorbed_wavelength_map_faceon_path", True, write=False)
    def absorbed_wavelength_map_faceon(self):

        """
        This function ...
        :return:
        """

        abs_wavelengths_faceon = self.mean_absorbed_wavelengths_faceon
        return abs_wavelengths_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_absorbed_wavelengths_edgeon(self):
        return self.absorbed_stellar_cube_edgeon.mean_wavelengths()

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "absorbed_wavelength_map_edgeon_path", True, write=False)
    def absorbed_wavelength_map_edgeon(self):

        """
        This function ...
        :return:
        """

        abs_wavelengths_edgeon = self.mean_absorbed_wavelengths_edgeon
        return abs_wavelengths_edgeon

    # -----------------------------------------------------------------
    #   3. Emitted
    # -----------------------------------------------------------------

    @lazyproperty
    def mean_emitted_wavelengths_earth(self):
        return self.observed_cube_dust.mean_wavelengths()

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "emitted_wavelength_map_earth_path", True, write=False)
    def emitted_wavelength_map_earth(self):

        """
        Thisnfunction ...
        :return:
        """

        em_wavelengths_earth = self.mean_emitted_wavelengths_earth
        em_wavelengths_earth.replace_by_nans_where_greater_than(200)
        em_wavelengths_earth.replace_by_nans_where_smaller_than(50)
        return em_wavelengths_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_emitted_wavelengths_faceon(self):
        return self.observed_cube_dust_faceon.mean_wavelengths()

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "emitted_wavelength_map_faceon_path", True, write=False)
    def emitted_wavelength_map_faceon(self):

        """
        This function ...
        :return:
        """

        em_wavelengths_faceon = self.mean_emitted_wavelengths_faceon
        return em_wavelengths_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_emitted_wavelengths_edgeon(self):
        return self.observed_cube_dust_edgeon.mean_wavelengths()

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "emitted_wavelength_map_edgeon_path", True, write=False)
    def emitted_wavelength_map_edgeon(self):

        """
        This function ...
        :return:
        """

        em_wavelengths_edgeon = self.mean_emitted_wavelengths_edgeon
        return em_wavelengths_edgeon

    # -----------------------------------------------------------------
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

        # Maps
        self.write_maps()

        # Wavelength maps
        self.write_wavelength_maps()

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
        return not self.has_emission_map_earth

    # -----------------------------------------------------------------

    @property
    def emission_map_earth_path(self):
        return fs.join(self.emissions_path, "earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_emission_map_earth(self):
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
        return not self.has_emission_map_faceon

    # -----------------------------------------------------------------

    @property
    def emission_map_faceon_path(self):
        return fs.join(self.emissions_path, "faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_emission_map_faceon(self):
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
        return not self.has_emission_map_edgeon

    # -----------------------------------------------------------------

    @property
    def emission_map_edgeon_path(self):
        return fs.join(self.emissions_path, "edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_emission_map_edgeon(self):
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
        return not self.has_absorption_map_earth

    # -----------------------------------------------------------------

    @property
    def absorption_map_earth_path(self):
        return fs.join(self.absorptions_path, "earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_absorption_map_earth(self):
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
        return not self.has_absorption_map_faceon

    # -----------------------------------------------------------------

    @property
    def absorption_map_faceon_path(self):
        return fs.join(self.absorptions_path, "faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_absorption_map_faceon(self):
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
        return not self.has_absorption_map_edgeon

    # -----------------------------------------------------------------

    @property
    def absorption_map_edgeon_path(self):
        return fs.join(self.absorptions_path, "edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_absorption_map_edgeon(self):
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

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_map_earth: self.write_map_earth()

        # Face-on
        if self.do_write_map_faceon: self.write_map_faceon()

        # Edge-on
        if self.do_write_map_edgeon: self.write_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_map_earth(self):
        return not self.has_map_earth

    # -----------------------------------------------------------------

    @property
    def map_earth_path(self):
        return fs.join(self.projected_energy_path, "earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_earth(self):
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
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the energy balance from the earth projection ...")

        # Save
        self.map_earth.saveto(self.map_earth_path)

    # -----------------------------------------------------------------

    @property
    def do_write_map_faceon(self):
        return not self.has_map_faceon

    # -----------------------------------------------------------------

    @property
    def map_faceon_path(self):
        return fs.join(self.projected_energy_path, "faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_faceon(self):
        return fs.is_file(self.map_faceon_path)

    # -----------------------------------------------------------------

    def remove_map_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_faceon_path)

    # -----------------------------------------------------------------

    def write_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the energy balance from the faceon projection ...")

        # Save
        self.map_faceon.saveto(self.map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_map_edgeon(self):
        return not self.has_map_edgeon

    # -----------------------------------------------------------------

    @property
    def map_edgeon_path(self):
        return fs.join(self.projected_energy_path, "edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_map_edgeon(self):
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

        # Inform the user
        log.info("Writing the map of the energy balance from the edgeon projection ...")

        # Save
        self.map_edgeon.saveto(self.map_edgeon_path)

    # -----------------------------------------------------------------

    def write_wavelength_maps(self):

        """
        This function ...
        :return:
        """

        # Scattered radiation
        self.write_scattered_wavelength_maps()

        # Absorbed radiation
        self.write_absorbed_wavelength_maps()

        # Emitted dust radiation
        self.write_emitted_wavelength_maps()

    # -----------------------------------------------------------------

    def write_scattered_wavelength_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_scat_wavelengths_earth: self.write_scattered_wavelength_map_earth()

        # Face-on
        if self.do_write_scat_wavelengths_faceon: self.write_scattered_wavelength_map_faceon()

        # Edge-on
        if self.do_write_scat_wavelengths_edgeon: self.write_scattered_wavelength_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_scat_wavelengths_earth(self):
        return not self.has_scattered_wavelength_map_earth

    # -----------------------------------------------------------------

    @property
    def scattered_wavelength_map_earth_path(self):
        return fs.join(self.wavelengths_path, "scattered_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_scattered_wavelength_map_earth(self):
        return fs.is_file(self.scattered_wavelength_map_earth_path)

    # -----------------------------------------------------------------

    def remove_scattered_wavelength_map_earth(self):

        """
        Thisfnuction ...
        :return:
        """

        fs.remove_file(self.scattered_wavelength_map_earth_path)

    # -----------------------------------------------------------------

    def write_scattered_wavelength_map_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the mean scattered wavelength from the earth projection ...")

        # Save
        self.scattered_wavelength_map_earth.saveto(self.scattered_wavelength_map_earth_path)

    # -----------------------------------------------------------------

    @property
    def do_write_scat_wavelengths_faceon(self):
        return not self.has_scattered_wavelength_map_faceon

    # -----------------------------------------------------------------

    @property
    def scattered_wavelength_map_faceon_path(self):
        return fs.join(self.wavelengths_path, "scattered_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_scattered_wavelength_map_faceon(self):
        return fs.is_file(self.scattered_wavelength_map_faceon_path)

    # -----------------------------------------------------------------

    def remove_scattered_wavelength_map_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.scattered_wavelength_map_faceon_path)

    # -----------------------------------------------------------------

    def write_scattered_wavelength_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the mean scattered wavelength from the faceon projection ...")

        # Save
        self.scattered_wavelength_map_faceon.saveto(self.scattered_wavelength_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_scat_wavelengths_edgeon(self):
        return not self.has_scattered_wavelength_map_edgeon

    # -----------------------------------------------------------------

    @property
    def scattered_wavelength_map_edgeon_path(self):
        return fs.join(self.wavelengths_path, "scattered_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_scattered_wavelength_map_edgeon(self):
        return fs.is_file(self.scattered_wavelength_map_edgeon_path)

    # -----------------------------------------------------------------

    def remove_scattered_wavelength_map_edgeon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.scattered_wavelength_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_scattered_wavelength_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the mean scattered wavelength from the edgeon projection ...")

        # Save
        self.scattered_wavelength_map_edgeon.saveto(self.scattered_wavelength_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_absorbed_wavelength_maps(self):

        """
        Thisf unction ...
        :return:
        """

        # Earth
        if self.do_write_abs_wavelengths_earth: self.write_absorbed_wavelength_map_earth()

        # Face-on
        if self.do_write_abs_wavelengths_faceon: self.write_absorbed_wavelength_map_faceon()

        # Edge-on
        if self.do_write_abs_wavelengths_edgeon: self.write_absorbed_wavelength_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_abs_wavelengths_earth(self):
        return not self.has_absorbed_wavelength_map_earth

    # -----------------------------------------------------------------

    @property
    def absorbed_wavelength_map_earth_path(self):
        return fs.join(self.wavelengths_path, "absorbed_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_absorbed_wavelength_map_earth(self):
        return fs.is_file(self.absorbed_wavelength_map_earth_path)

    # -----------------------------------------------------------------

    def remove_absorbed_wavelength_map_earth(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.absorbed_wavelength_map_earth_path)

    # -----------------------------------------------------------------

    def write_absorbed_wavelength_map_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the mean absorbed wavelength from the earth projection ...")

        # Save
        self.absorbed_wavelength_map_earth.saveto(self.absorbed_wavelength_map_earth_path)

    # -----------------------------------------------------------------

    @property
    def do_write_abs_wavelengths_faceon(self):
        return not self.has_absorbed_wavelength_map_faceon

    # -----------------------------------------------------------------

    @property
    def absorbed_wavelength_map_faceon_path(self):
        return fs.join(self.wavelengths_path, "absorbed_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_absorbed_wavelength_map_faceon(self):
        return fs.is_file(self.absorbed_wavelength_map_faceon_path)

    # -----------------------------------------------------------------

    def remove_absorbed_wavelength_map_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.absorbed_wavelength_map_faceon_path)

    # -----------------------------------------------------------------

    def write_absorbed_wavelength_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the mean absorbed wavelength from the faceon projection ...")

        # Save
        self.absorbed_wavelength_map_faceon.saveto(self.absorbed_wavelength_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_abs_wavelengths_edgeon(self):
        return not self.has_absorbed_wavelength_map_edgeon

    # -----------------------------------------------------------------

    @property
    def absorbed_wavelength_map_edgeon_path(self):
        return fs.join(self.wavelengths_path, "absorbed_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_absorbed_wavelength_map_edgeon(self):
        return fs.is_file(self.absorbed_wavelength_map_edgeon_path)

    # -----------------------------------------------------------------

    def remove_absorbed_wavelength_map_edgeon(self):

        """
        Thisfunction ...
        :return:
        """

        fs.remove_file(self.absorbed_wavelength_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_absorbed_wavelength_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of the mean absorbed wavelength from the edgeon projection ...")

        # Save
        self.absorbed_wavelength_map_edgeon.saveto(self.absorbed_wavelength_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_emitted_wavelength_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_em_wavelengths_earth: self.write_emitted_wavelength_map_earth()

        # Face-on
        if self.do_write_em_wavelengths_faceon: self.write_emitted_wavelength_map_faceon()

        # Edge-on
        if self.do_write_em_wavelengths_edgeon: self.write_emitted_wavelength_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_em_wavelengths_earth(self):
        return not self.has_emitted_wavelength_map_earth

    # -----------------------------------------------------------------

    @property
    def emitted_wavelength_map_earth_path(self):
        return fs.join(self.wavelengths_path, "emitted_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_emitted_wavelength_map_earth(self):
        return fs.is_file(self.emitted_wavelength_map_earth_path)

    # -----------------------------------------------------------------

    def remove_emitted_wavelength_map_earth(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.emitted_wavelength_map_earth_path)

    # -----------------------------------------------------------------

    def write_emitted_wavelength_map_earth(self):

        """
        Thisf unction ...
        :return:
        """

        self.emitted_wavelength_map_earth.saveto(self.emitted_wavelength_map_earth_path)

    # -----------------------------------------------------------------

    @property
    def do_write_em_wavelengths_faceon(self):
        return not self.has_emitted_wavelength_map_faceon

    # -----------------------------------------------------------------

    @property
    def emitted_wavelength_map_faceon_path(self):
        return fs.join(self.wavelengths_path, "emitted_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_emitted_wavelength_map_faceon(self):
        return fs.is_file(self.emitted_wavelength_map_faceon_path)

    # -----------------------------------------------------------------

    def remove_emitted_wavelength_map_faceon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.emitted_wavelength_map_faceon_path)

    # -----------------------------------------------------------------

    def write_emitted_wavelength_map_faceon(self):

        """
        This function ...
        :return:
        """

        self.emitted_wavelength_map_faceon.saveto(self.emitted_wavelength_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_em_wavelengths_edgeon(self):
        return not self.has_emitted_wavelength_map_edgeon

    # -----------------------------------------------------------------

    @property
    def emitted_wavelength_map_edgeon_path(self):
        return fs.join(self.wavelengths_path, "emitted_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_emitted_wavelength_map_edgeon(self):
        return fs.is_file(self.emitted_wavelength_map_edgeon_path)

    # -----------------------------------------------------------------

    def remove_emitted_wavelength_map_edgeon(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.emitted_wavelength_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_emitted_wavelength_map_edgeon(self):

        """
        This function ...
        :return:
        """

        self.emitted_wavelength_map_edgeon.saveto(self.emitted_wavelength_map_edgeon_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the maps
        self.plot_maps()

        # Plot the wavelength maps
        self.plot_wavelength_maps()

    # -----------------------------------------------------------------

    def plot_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_map_earth: self.plot_map_earth()

        # Faceon
        if self.do_plot_map_faceon: self.plot_map_faceon()

        # Edgeon
        if self.do_plot_map_edgeon: self.plot_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_plot_map_earth(self):
        return self.do_earth and not self.has_map_earth_plot and self.map_earth is not None

    # -----------------------------------------------------------------

    @property
    def map_earth_plot_path(self):
        return fs.join(self.projected_energy_path, "earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_map_earth_plot(self):
        return fs.is_file(self.map_earth_plot_path)

    # -----------------------------------------------------------------

    def remove_map_earth_plot(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_earth_plot_path)

    # -----------------------------------------------------------------

    def plot_map_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the earth map ...")

        # Plot
        #plot_frame_contours(self.map_earth, nlevels=5, path=self.map_earth_plot_path, plot_data=True, interval=[-1,1], data_cmap="RdBu")
        #plot_frame(self.map_earth, path=self.map_earth_plot_path, interval=(-1,1), cmap="RdBu", colorbar=True, scale="linear")
        #plot_frame(self.map_earth, path=self.map_earth_plot_path, interval=(-1.,1.), cmap="RdBu", colorbar=True, scale="linear")

        # Regularize
        #regularize_frame(self.map_earth, dilate_nans=True, dilation_radius=4, dilation_niterations=5, cutoff_above=2,
        #                 cutoff_below=-2, no_infs=True, replace_infs=nan_value)

        regularize_frame(self.map_earth, dilate_nans=True, dilation_radius=4, dilation_niterations=4, cutoff_above=1,
                        cutoff_below=0, no_infs=True, replace_infs=nan_value)

        # Plot
        #plot_frame_contours(self.map_earth, interval=(-2,2), data_cmap="RdBu", colorbar=True, nlevels=5, path=self.map_earth_plot_path, plot_data=True, single_colour="green")

        # Plot
        #plot_map(self.map_earth, interval=(-2,2), cmap="RdBu", path=self.map_earth_plot_path)
        plot_map(self.map_earth, interval=(0,1), path=self.map_earth_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_map_faceon(self):
        return self.do_faceon and not self.has_map_faceon_plot and self.map_faceon is not None

    # -----------------------------------------------------------------

    @property
    def map_faceon_plot_path(self):
        return fs.join(self.projected_energy_path, "faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_map_faceon_plot(self):
        return fs.is_file(self.map_faceon_plot_path)

    # -----------------------------------------------------------------

    def remove_map_faceon_plot(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_faceon_plot_path)

    # -----------------------------------------------------------------

    def plot_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the faceon map ...")

        # Plot
        #plot_frame_contours(self.map_faceon, nlevels=5, path=self.map_faceon_plot_path, plot_data=True, interval=[-1, 1], data_cmap="RdBu")
        plot_frame(self.map_faceon, path=self.map_faceon_plot_path, interval=(-1., 1.), cmap="RdBu")

    # -----------------------------------------------------------------

    @property
    def do_plot_map_edgeon(self):
        return self.do_edgeon and not self.has_map_edgeon_plot and self.map_edgeon is not None

    # -----------------------------------------------------------------

    @property
    def map_edgeon_plot_path(self):
        return fs.join(self.projected_energy_path, "edgeon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_map_edgeon_plot(self):
        return fs.is_file(self.map_edgeon_plot_path)

    # -----------------------------------------------------------------

    def remove_map_edgeon_plot(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.map_edgeon_plot_path)

    # -----------------------------------------------------------------

    def plot_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the edgeon map ...")

        # Plot
        #plot_frame_contours(self.map_edgeon, nlevels=5, path=self.map_edgeon_plot_path, plot_data=True, interval=[-1, 1], data_cmap="RdBu")
        plot_frame(self.map_edgeon, path=self.map_edgeon_plot_path, interval=[-1, 1], cmap="RdBu")

    # -----------------------------------------------------------------

    def plot_wavelength_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelength maps ...")

        # Scattered radiation
        self.plot_scattered_wavelength_maps()

        # Absorbed radiation
        self.plot_absorbed_wavelength_maps()

        # Emitted dust radiation
        self.plot_emitted_wavelength_maps()

    # -----------------------------------------------------------------

    def plot_scattered_wavelength_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_scat_wavelengths_earth: self.plot_scattered_wavelength_map_earth()

        # Faceon
        if self.do_plot_scat_wavelengths_faceon: self.plot_scattered_wavelength_map_faceon()

        # Edgeon
        if self.do_plot_scat_wavelengths_edgeon: self.plot_scattered_wavelength_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_plot_scat_wavelengths_earth(self):
        return self.do_earth and not self.has_scat_wavelengths_earth_plot and self.scat_wavelengths_earth is not None

    # -----------------------------------------------------------------

    @property
    def scat_wavelengths_earth_plot_path(self):
        return fs.join(self.wavelengths_path, "scattered_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_scat_wavelengths_earth_plot(self):
        return fs.is_file(self.scat_wavelengths_earth_plot_path)

    # -----------------------------------------------------------------

    def plot_scattered_wavelength_map_earth(self):

        """
        This function ...
        :return:
        """

        plot_map(self.scat_wavelengths_earth, path=self.scat_wavelengths_earth_plot_path, interval=(0.7,1))

    # -----------------------------------------------------------------

    @property
    def do_plot_scat_wavelengths_faceon(self):
        return self.do_faceon and not self.has_scat_wavelengths_faceon_plot and self.scat_wavelengths_faceon is not None

    # -----------------------------------------------------------------

    @property
    def scat_wavelengths_faceon_plot_path(self):
        return fs.join(self.wavelengths_path, "scattered_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_scat_wavelengths_faceon_plot(self):
        return fs.is_file(self.scat_wavelengths_faceon_plot_path)

    # -----------------------------------------------------------------

    def plot_scattered_wavelength_map_faceon(self):

        """
        This function ...
        :return:
        """

        plot_map(self.scat_wavelengths_faceon, path=self.scat_wavelengths_faceon_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_scat_wavelengths_edgeon(self):
        return self.do_edgeon and not self.has_scat_wavelengths_edgeon_plot and self.scat_wavelengths_edgeon is not None

    # -----------------------------------------------------------------

    @property
    def scat_wavelengths_edgeon_plot_path(self):
        return fs.join(self.wavelengths_path, "scattered_edgeon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_scat_wavelengths_edgeon_plot(self):
        return fs.is_file(self.scat_wavelengths_edgeon_plot_path)

    # -----------------------------------------------------------------

    def plot_scattered_wavelength_map_edgeon(self):

        """
        This function ...
        :return:
        """

        plot_map(self.scat_wavelengths_edgeon, path=self.scat_wavelengths_edgeon_plot_path, interval=(0.7,1), cmap="inferno_r")

    # -----------------------------------------------------------------

    def plot_absorbed_wavelength_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_abs_wavelengths_earth: self.plot_absorbed_wavelength_map_earth()

        # Faceon
        if self.do_plot_abs_wavelengths_faceon: self.plot_absorbed_wavelength_map_faceon()

        # Edgeon
        if self.do_plot_abs_wavelengths_edgeon: self.plot_absorbed_wavelength_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_plot_abs_wavelengths_earth(self):
        return self.do_earth and not self.has_abs_wavelengths_earth_plot and self.abs_wavelengths_earth is not None

    # -----------------------------------------------------------------

    @property
    def abs_wavelengths_earth_plot_path(self):
        return fs.join(self.wavelengths_path, "absorbed_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_abs_wavelengths_earth_plot(self):
        return fs.is_file(self.abs_wavelengths_earth_plot_path)

    # -----------------------------------------------------------------

    def plot_absorbed_wavelength_map_earth(self):

        """
        This function ...
        :return:
        """

        plot_map(self.abs_wavelengths_earth, path=self.abs_wavelengths_earth_plot_path, interval=(0.1,1), cmap="inferno_r")

    # -----------------------------------------------------------------

    @property
    def do_plot_abs_wavelengths_faceon(self):
        return self.do_faceon and not self.has_abs_wavelengths_faceon_plot and self.abs_wavelengths_faceon is not None

    # -----------------------------------------------------------------

    @property
    def abs_wavelengths_faceon_plot_path(self):
        return fs.join(self.wavelengths_path, "absorbed_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_abs_wavelengths_faceon_plot(self):
        return fs.is_file(self.abs_wavelengths_faceon_plot_path)

    # -----------------------------------------------------------------

    def plot_absorbed_wavelength_map_faceon(self):

        """
        This function ...
        :return:
        """

        plot_map(self.abs_wavelengths_faceon, path=self.abs_wavelengths_faceon_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_abs_wavelengths_edgeon(self):
        return self.do_edgeon and not self.has_abs_wavelengths_edgeon_plot and self.abs_wavelengths_edgeon is not None

    # -----------------------------------------------------------------

    @property
    def abs_wavelengths_edgeon_plot_path(self):
        return fs.join(self.wavelengths_path, "absorbed_edgeon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_abs_wavelengths_edgeon_plot(self):
        return fs.is_file(self.abs_wavelengths_edgeon_plot_path)

    # -----------------------------------------------------------------

    def plot_absorbed_wavelength_map_edgeon(self):

        """
        This function ...
        :return:
        """

        plot_map(self.abs_wavelengths_edgeon, path=self.abs_wavelengths_earth_plot_path)

    # -----------------------------------------------------------------

    def plot_emitted_wavelength_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_em_wavelengths_earth: self.plot_emitted_wavelength_map_earth()

        # Faceon
        if self.do_plot_em_wavelengths_faceon: self.plot_emitted_wavelength_map_faceon()

        # Edgeon
        if self.do_plot_em_wavelengths_edgeon: self.plot_emitted_wavelength_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_plot_em_wavelengths_earth(self):
        return self.do_earth and not self.has_em_wavelengths_earth_plot and self.em_wavelengths_earth is not None

    # -----------------------------------------------------------------

    @property
    def em_wavelengths_earth_plot_path(self):
        return fs.join(self.wavelengths_path, "emitted_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_em_wavelengths_earth_plot(self):
        return fs.is_file(self.em_wavelengths_earth_plot_path)

    # -----------------------------------------------------------------

    def plot_emitted_wavelength_map_earth(self):

        """
        This functino ...
        :return:
        """

        plot_map(self.em_wavelengths_earth, path=self.em_wavelengths_earth_plot_path, interval=(80,150), cmap="inferno_r", background_color="black")

    # -----------------------------------------------------------------

    @property
    def do_plot_em_wavelengths_faceon(self):
        return self.do_faceon and not self.has_em_wavelengths_faceon_plot and self.em_wavelengths_faceon is not None

    # -----------------------------------------------------------------

    @property
    def em_wavelengths_faceon_plot_path(self):
        return fs.join(self.wavelengths_path, "emitted_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_em_wavelengths_faceon_plot(self):
        return fs.is_file(self.em_wavelengths_faceon_plot_path)

    # -----------------------------------------------------------------

    def plot_emitted_wavelength_map_faceon(self):

        """
        This function ...
        :return:
        """

        plot_map(self.em_wavelengths_faceon, path=self.em_wavelengths_faceon_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_em_wavelengths_edgeon(self):
        return self.do_edgeon and not self.has_em_wavelengths_edgeon_plot and self.em_wavelengths_edgeon is not None

    # -----------------------------------------------------------------

    @property
    def em_wavelengths_edgeon_plot_path(self):
        return fs.join(self.wavelengths_path, "emitted_edgeon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_em_wavelengths_edgeon_plot(self):
        return fs.is_file(self.em_wavelengths_edgeon_plot_path)

    # -----------------------------------------------------------------

    def plot_emitted_wavelength_map_edgeon(self):

        """
        This function ...
        :return:
        """

        plot_map(self.em_wavelengths_edgeon, path=self.em_wavelengths_edgeon_plot_path)

# -----------------------------------------------------------------
