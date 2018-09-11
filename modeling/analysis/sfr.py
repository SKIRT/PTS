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
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import AnalysisComponent, AnalysisRunComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.core.frame import Frame
from ...magic.core.image import Image
from ...core.tools.utils import lazyproperty, lazyfileproperty
from ..core.data import Data3D
from ..projection.data import project_data
from ..core.model import oliver_stellar_mass, salim_fuv_to_sfr, kennicutt_evans_fuv_to_sfr
from ...core.units.parsing import parse_unit as u
from ...magic.tools.plotting import plot_map
from ...magic.core.list import uniformize

# -----------------------------------------------------------------

projected_name = "projected"
cell_name = "cell"

# -----------------------------------------------------------------

salim_name = "Salim"
ke_name = "K&E"
mappings_name = "Mappings"
mappings_ke_name = "Mappings and K&E"

# -----------------------------------------------------------------

methods = OrderedDict()
methods[salim_name] = "Salim conversion from intrinsic FUV luminosity of unevolved stars"
methods[ke_name] = "Kennicutt & Evans conversion from intrinsic FUV luminosity of unevolved stars"
methods[mappings_name] = "MAPPINGS SFR of ionizing stars"
methods[mappings_ke_name] = "MAPPINGS SFR of ionizing stars + Kennicutt & Evans conversion from intrinsic FUV luminosity of young stars"

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

        # Replot?
        if self.config.replot:
            self.config.replot_projected = True
            self.config.replot_cell_maps = True

        # Replot projected
        if self.config.replot_projected:
            self.config.replot_projected_sfr = True
            self.config.replot_projected_mass = True
            self.config.replot_projected_ssfr = True

        # SFR
        if self.config.replot_sfr:
            self.config.replot_projected_sfr = True

        # Mass
        if self.config.replot_mass:
            self.config.replot_projected_mass = True

        # sSFR
        if self.config.replot_ssfr:
            self.config.replot_projected_ssfr = True
            self.config.replot_cell_ssfr_maps = True

        # Replot SFR
        if self.config.replot_projected_sfr:

            # Salim
            if self.has_projected_sfr_salim_earth_map_plot: fs.remove_file(self.projected_sfr_salim_earth_map_plot_path)
            if self.has_projected_sfr_salim_faceon_map_plot: fs.remove_file(self.projected_sfr_salim_faceon_map_plot_path)

            # K&E
            if self.has_projected_sfr_ke_earth_map_plot: fs.remove_file(self.projected_sfr_ke_earth_map_plot_path)
            if self.has_projected_sfr_ke_faceon_map_plot: fs.remove_file(self.projected_sfr_ke_faceon_map_plot_path)

            # MAPPINGS
            if self.has_projected_sfr_mappings_earth_map_plot: fs.remove_file(self.projected_sfr_mappings_earth_map_plot_path)
            if self.has_projected_sfr_mappings_faceon_map_plot: fs.remove_file(self.projected_sfr_mappings_faceon_map_plot_path)

            # MAPPINGS + K&E
            if self.has_projected_sfr_mappings_ke_earth_map_plot: fs.remove_file(self.projected_sfr_mappings_ke_earth_map_plot_path)
            if self.has_projected_sfr_mappings_ke_faceon_map_plot: fs.remove_file(self.projected_sfr_mappings_ke_faceon_map_plot_path)

        # Replot mass
        if self.config.replot_projected_mass:

            if self.has_projected_mass_earth_map_plot: fs.remove_file(self.projected_mass_earth_map_plot_path)
            if self.has_projected_mass_faceon_map_plot: fs.remove_file(self.projected_mass_faceon_map_plot_path)
        
        # Replot sSFR
        if self.config.replot_projected_ssfr:

            # Salim
            if self.has_projected_ssfr_salim_earth_map_plot: fs.remove_file(self.projected_ssfr_salim_earth_map_plot_path)
            if self.has_projected_ssfr_salim_faceon_map_plot: fs.remove_file(self.projected_ssfr_salim_faceon_map_plot_path)

            # K&E
            if self.has_projected_ssfr_ke_earth_map_plot: fs.remove_file(self.projected_ssfr_ke_earth_map_plot_path)
            if self.has_projected_ssfr_ke_faceon_map_plot: fs.remove_file(self.projected_ssfr_ke_faceon_map_plot_path)

            # MAPPINGS
            if self.has_projected_ssfr_mappings_earth_map_plot: fs.remove_file(self.projected_ssfr_mappings_earth_map_plot_path)
            if self.has_projected_ssfr_mappings_faceon_map_plot: fs.remove_file(self.projected_ssfr_mappings_faceon_map_plot_path)

            # MAPPINGS + K&E
            if self.has_projected_ssfr_mappings_ke_earth_map_plot: fs.remove_file(self.projected_ssfr_mappings_ke_earth_map_plot_path)
            if self.has_projected_ssfr_mappings_ke_faceon_map_plot: fs.remove_file(self.projected_ssfr_mappings_ke_faceon_map_plot_path)

        # Replot cell maps
        if self.config.replot_cell_maps:

            #self.config.replot_cell_sfr_maps = True
            #self.config.replot_cell_mass_maps = True
            self.config.replot_cell_ssfr_maps = True

        # Cell sSFR maps
        if self.config.replot_cell_ssfr_maps:

            # Salim
            if self.has_cell_ssfr_salim_map_plot: fs.remove_file(self.cell_ssfr_salim_map_plot_path)

            # K&E
            if self.has_cell_ssfr_ke_map_plot: fs.remove_file(self.cell_ssfr_ke_map_plot_path)

            # MAPPINGS
            if self.has_cell_ssfr_mappings_map_plot: fs.remove_file(self.cell_ssfr_mappings_map_plot_path)

            # MAPPINGS + K&E
            if self.has_cell_ssfr_mappings_ke_map_plot: fs.remove_file(self.cell_ssfr_mappings_ke_map_plot_path)

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
    # SFR MAPS
    #   1. SALIM
    # -----------------------------------------------------------------

    @property
    def projected_sfr_salim_earth_path(self):
        return fs.join(self.projected_path, "sfr_salim_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_salim_earth(self):
        return fs.is_file(self.projected_sfr_salim_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_salim_earth_path", True, write=False)
    def sfr_salim_earth_map(self):

        """
        projected star formation rate from the earth projection
        :return:
        """

        return self.model.total_star_formation_rate_map_earth_salim

    # -----------------------------------------------------------------

    @property
    def projected_sfr_salim_faceon_path(self):
        return fs.join(self.projected_path, "sfr_salim_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_salim_faceon(self):
        return fs.is_file(self.projected_sfr_salim_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_salim_faceon_path", True, write=False)
    def sfr_salim_faceon_map(self):

        """
        projected star formation rate from the faceon projection
        :return:
        """

        return self.model.total_star_formation_rate_map_faceon_salim

    # -----------------------------------------------------------------
    #   2. KENNICUTT & EVANS
    # -----------------------------------------------------------------

    @property
    def projected_sfr_ke_earth_path(self):
        return fs.join(self.projected_path, "sfr_ke_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_ke_earth(self):
        return fs.is_file(self.projected_sfr_ke_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_ke_earth_path", True, write=False)
    def sfr_ke_earth_map(self):

        """
        Thisf unction ...
        :return:
        """

        return self.model.total_star_formation_rate_map_earth_ke

    # -----------------------------------------------------------------

    @property
    def projected_sfr_ke_faceon_path(self):
        return fs.join(self.projected_path, "sfr_ke_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_ke_faceon(self):
        return fs.is_file(self.projected_sfr_ke_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_ke_faceon_path", True, write=False)
    def sfr_ke_faceon_map(self):

        """
        This function ...
        :return:
        """

        return self.model.total_star_formation_rate_map_faceon_ke

    # -----------------------------------------------------------------
    #   3. MAPPINGS
    # -----------------------------------------------------------------

    @property
    def projected_sfr_mappings_earth_path(self):
        return fs.join(self.projected_path, "sfr_mappings_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_mappings_earth(self):
        return fs.is_file(self.projected_sfr_mappings_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_mappings_earth_path", True, write=False)
    def sfr_mappings_earth_map(self):

        """
        This function ...
        :return:
        """

        return self.model.star_formation_rate_map_earth

    # -----------------------------------------------------------------

    @property
    def projected_sfr_mappings_faceon_path(self):
        return fs.join(self.projected_path, "sfr_mappings_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_mappings_faceon(self):
        return fs.is_file(self.projected_sfr_mappings_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_mappings_faceon_path", True, write=False)
    def sfr_mappings_faceon_map(self):

        """
        This function ...
        :return:
        """

        return self.model.star_formation_rate_map_faceon

    # -----------------------------------------------------------------
    #   4. MAPPINGS + KENNICUTT & EVANS
    # -----------------------------------------------------------------

    @property
    def projected_sfr_mappings_ke_earth_path(self):
        return fs.join(self.projected_path, "sfr_mappings_ke_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_mappings_ke_earth(self):
        return fs.is_file(self.projected_sfr_mappings_ke_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_mappings_ke_earth_path", True, write=False)
    def sfr_mappings_ke_earth_map(self):

        """
        This function ...
        :return:
        """

        mappings, young_ke = uniformize(self.model.star_formation_rate_map_earth, self.model.young_star_formation_rate_map_earth_ke)
        return mappings + young_ke

    # -----------------------------------------------------------------

    @property
    def projected_sfr_mappings_ke_faceon_path(self):
        return fs.join(self.projected_path, "sfr_mappings_ke_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_mappings_ke_faceon(self):
        return fs.is_file(self.projected_sfr_mappings_ke_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_mappings_ke_faceon_path", True, write=False)
    def sfr_mappings_ke_faceon_map(self):

        """
        This function ...
        :return:
        """

        mappings, young_ke = uniformize(self.model.star_formation_rate_map_faceon, self.model.young_star_formation_rate_map_faceon_ke)
        return mappings + young_ke

    # -----------------------------------------------------------------
    # STELLAR MASS MAPS
    # -----------------------------------------------------------------

    @property
    def projected_mass_earth_path(self):
        return fs.join(self.projected_path, "mass_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_mass_earth(self):
        return fs.is_file(self.projected_mass_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_mass_earth_path", True, write=False)
    def stellar_mass_earth_map(self):

        """
        projected stellar mass from the earth projection
        :return:
        """

        return self.model.total_stellar_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def projected_mass_faceon_path(self):
        return fs.join(self.projected_path, "mass_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_mass_faceon(self):
        return fs.is_file(self.projected_mass_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_mass_faceon_path", True, write=False)
    def stellar_mass_faceon_map(self):

        """
        This function ...
        :return:
        """

        return self.model.total_stellar_mass_map_faceon

    # -----------------------------------------------------------------
    # sSFR MAPS
    #   1. SALIM
    # -----------------------------------------------------------------

    @property
    def projected_ssfr_salim_earth_path(self):
        return fs.join(self.projected_path, "ssfr_salim_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_salim_earth(self):
        return fs.is_file(self.projected_ssfr_salim_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_salim_earth_path", True, write=False)
    def ssfr_salim_earth_map(self):

        """
        projected specific star formation rate from the earth projection
        :return:
        """

        return self.model.total_ssfr_map_earth_salim

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_salim_faceon_path(self):
        return fs.join(self.projected_path, "ssfr_salim_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_salim_faceon(self):
        return fs.is_file(self.projected_ssfr_salim_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_salim_faceon_path", True, write=False)
    def ssfr_salim_faceon_map(self):

        """
        This function ...
        :return:
        """

        return self.model.total_ssfr_map_faceon_salim

    # -----------------------------------------------------------------
    #   2. KENNICUTT & EVANS
    # -----------------------------------------------------------------

    @property
    def projected_ssfr_ke_earth_path(self):
        return fs.join(self.projected_path, "ssfr_ke_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_ke_earth(self):
        return fs.is_file(self.projected_ssfr_ke_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_ke_earth_path", True, write=False)
    def ssfr_ke_earth_map(self):

        """
        This function ...
        :return:
        """

        return self.model.total_ssfr_map_earth_ke

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_ke_faceon_path(self):
        return fs.join(self.projected_path, "ssfr_ke_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_ke_faceon(self):
        return fs.is_file(self.projected_ssfr_ke_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_ke_faceon_path", True, write=False)
    def ssfr_ke_faceon_map(self):

        """
        This function ...
        :return:
        """

        return self.model.total_ssfr_map_faceon_ke

    # -----------------------------------------------------------------
    #   3. MAPPINGS
    # -----------------------------------------------------------------

    @property
    def projected_ssfr_mappings_earth_path(self):
        return fs.join(self.projected_path, "ssfr_mappings_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_mappings_earth(self):
        return fs.is_file(self.projected_ssfr_mappings_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_mappings_earth_path", True, write=False)
    def ssfr_mappings_earth_map(self):

        """
        This function ...
        :return:
        """

        sfr, stellar_mass = uniformize(self.sfr_mappings_earth_map, self.stellar_mass_earth_map, convert=False)
        return sfr / stellar_mass

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_mappings_faceon_path(self):
        return fs.join(self.projected_path, "ssfr_mappings_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_mappings_faceon(self):
        return fs.is_file(self.projected_ssfr_mappings_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_mappings_faceon_path", True, write=False)
    def ssfr_mappings_faceon_map(self):

        """
        This function ...
        :return:
        """

        sfr, stellar_mass = uniformize(self.sfr_mappings_faceon_map, self.stellar_mass_faceon_map, convert=False)
        return sfr / stellar_mass

    # -----------------------------------------------------------------
    #   4. MAPPINGS + KENNICUTT & EVANS
    # -----------------------------------------------------------------

    @property
    def projected_ssfr_mappings_ke_earth_path(self):
        return fs.join(self.projected_path, "ssfr_mappings_ke_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_mappings_ke_earth(self):
        return fs.is_file(self.projected_ssfr_mappings_ke_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_mappings_ke_earth_path", True, write=False)
    def ssfr_mappings_ke_earth_map(self):

        """
        This function ...
        :return:
        """

        sfr, stellar_mass = uniformize(self.sfr_mappings_ke_earth_map, self.stellar_mass_earth_map, convert=False)
        return sfr / stellar_mass

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_mappings_ke_faceon_path(self):
        return fs.join(self.projected_path, "ssfr_mappings_ke_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_mappings_ke_faceon(self):
        return fs.is_file(self.projected_ssfr_mappings_ke_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_mappings_ke_faceon_path", True, write=False)
    def ssfr_mappings_ke_faceon_map(self):

        """
        This function ...
        :return:
        """

        sfr, stellar_mass = uniformize(self.sfr_mappings_ke_faceon_map, self.stellar_mass_faceon_map, convert=False)
        return sfr / stellar_mass

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def sfr(self):
        return self.model.sfr

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_msun_yr(self):
        return self.sfr.to("Msun/yr").value

    # -----------------------------------------------------------------

    @property
    def length_unit(self):
        return "pc"

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

    @lazyfileproperty(Data3D, "cell_fuv_path", True, write=False)
    def fuv_data(self):

        """
        cell FUV luminosity data
        :return:
        """

        # Create the data
        return Data3D(self.fuv_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                      self.unevolved_cell_fuv_luminosities, length_unit=self.length_unit, unit=self.fuv_luminosity_unit,
                      description=self.fuv_description, distance=self.galaxy_distance, wavelength=self.fuv_wavelength)

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

    @property
    def i1_name(self):
        return "I1"

    # -----------------------------------------------------------------

    @property
    def i1_description(self):
        return "Intrinsic I1 luminosity of evolved stars"

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_i1_path", True, write=False)
    def i1_data(self):

        """
        cell I1 luminosity data
        :return:
        """

        return Data3D(self.i1_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                      self.old_cell_i1_luminosities, length_unit=self.length_unit, unit=self.i1_luminosity_unit,
                      description=self.i1_description, distance=self.galaxy_distance,
                      wavelength=self.i1_wavelength)

    # -----------------------------------------------------------------

    @property
    def cell_sfr_salim_path(self):
        return fs.join(self.cell_path, "sfr_salim.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_salim(self):
        return fs.is_file(self.cell_sfr_salim_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_sfr_salim_path", True, write=False)
    def sfr_salim_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell star formation rate (Salim) ...")

        # Calculate the SFR data from the intrinsic FUV data
        return salim_fuv_to_sfr(self.fuv_data)

    # -----------------------------------------------------------------

    @property
    def cell_sfr_ke_path(self):
        return fs.join(self.cell_path, "sfr_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_ke(self):
        return fs.is_file(self.cell_sfr_ke_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_sfr_ke_path", True, write=False)
    def sfr_ke_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell star formation rate (Kennicutt & Evans) ...")

        # Calculate
        return kennicutt_evans_fuv_to_sfr(self.fuv_data)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_mappings_values(self):
        return self.sfr_cell_normalized_mass * self.sfr_msun_yr

    # -----------------------------------------------------------------

    @property
    def cell_sfr_mappings_path(self):
        return fs.join(self.cell_path, "sfr_mappings.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_mappings(self):
        return fs.is_file(self.cell_sfr_mappings_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_sfr_mappings_path", True, write=False)
    def sfr_mappings_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell star formation rate (MAPPINGS) ...")

        # Create
        return Data3D(self.sfr_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                      self.sfr_mappings_values, length_unit=self.length_unit, unit="Msun/yr",
                      description=self.sfr_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def cell_sfr_mappings_ke_path(self):
        return fs.join(self.cell_path, "sfr_mappings_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_mappings_ke(self):
        return fs.is_file(self.cell_sfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_sfr_mappings_ke_path", True, write=False)
    def sfr_mappings_ke_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell star formation rate (MAPPINGS + Kennicutt & Evans) ...")

        # Calculate in Msun/yr
        sfr_values = self.sfr_mappings_values + kennicutt_evans_fuv_to_sfr(self.young_cell_fuv_luminosities, unit=self.fuv_luminosity_unit)

        # Create
        return Data3D(self.sfr_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                      sfr_values, length_unit=self.length_unit, unit="Msun/yr",
                      description=self.sfr_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_mass_path", True, write=False)
    def stellar_mass_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell stellar mass ...")

        # Calculate the stellar mass data from the I1 data
        return oliver_stellar_mass(self.i1_data, hubble_type=self.hubble_stage_type, hubble_subtype=self.hubble_stage_subtype)

    # -----------------------------------------------------------------

    @property
    def cell_sfrs_salim(self):
        return self.sfr_salim_data.values

    # -----------------------------------------------------------------

    @property
    def sfr_salim_unit(self):
        return self.sfr_salim_data.unit

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
    def cell_ssfrs_salim(self):
        return self.cell_sfrs_salim / self.cell_masses

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_salim_unit(self):
        return self.sfr_salim_unit / self.stellar_mass_unit

    # -----------------------------------------------------------------

    @property
    def ssfr_name(self):
        return "sSFR"

    # -----------------------------------------------------------------

    @property
    def ssfr_description(self):
        return "Specific star formation rate"

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_salim_path(self):
        return fs.join(self.cell_path, "ssfr_salim.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_salim(self):
        return fs.is_file(self.cell_ssfr_salim_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_ssfr_salim_path", True, write=False)
    def ssfr_salim_data(self):

        """
        cell specific star formation rate
        :return:
        """

        # Create the data
        return Data3D(self.ssfr_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates, self.cell_ssfrs_salim,
                      length_unit=self.length_unit, unit=self.ssfr_salim_unit, description=self.ssfr_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def cell_sfrs_ke(self):
        return self.sfr_ke_data.values

    # -----------------------------------------------------------------

    @property
    def sfr_ke_unit(self):
        return self.sfr_ke_data.unit

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfrs_ke(self):
        return self.cell_sfrs_ke / self.cell_masses

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_ke_unit(self):
        return self.sfr_ke_unit / self.stellar_mass_unit

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_ke_path(self):
        return fs.join(self.cell_path, "ssfr_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_ke(self):
        return fs.is_file(self.cell_ssfr_ke_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_ssfr_ke_path", True, write=False)
    def ssfr_ke_data(self):

        """
        This function ...
        :return:
        """

        # Create the data
        return Data3D(self.ssfr_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates, self.cell_ssfrs_ke,
                      length_unit=self.length_unit, unit=self.ssfr_ke_unit, description=self.ssfr_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def cell_sfrs_mappings(self):
        return self.sfr_mappings_data.values

    # -----------------------------------------------------------------

    @property
    def sfr_mappings_unit(self):
        return self.sfr_mappings_data.unit

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfrs_mappings(self):
        return self.cell_sfrs_mappings / self.cell_masses

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_mappings_unit(self):
        return self.sfr_mappings_unit / self.stellar_mass_unit

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_path(self):
        return fs.join(self.cell_path, "ssfr_mappings.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_mappings(self):
        return fs.is_file(self.cell_ssfr_mappings_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_ssfr_mappings_path", True, write=False)
    def ssfr_mappings_data(self):

        """
        This function ...
        :return:
        """

        # Create the data
        return Data3D(self.ssfr_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates,
                      self.cell_ssfrs_mappings, length_unit=self.length_unit, unit=self.ssfr_mappings_unit, description=self.ssfr_description,
                      distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def cell_sfrs_mappings_ke(self):
        return self.sfr_mappings_ke_data.values

    # -----------------------------------------------------------------

    @property
    def sfr_mappings_ke_unit(self):
        return self.sfr_mappings_ke_data.unit

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfrs_mappings_ke(self):
        return self.cell_sfrs_mappings_ke / self.cell_masses

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_mappings_ke_unit(self):
        return self.sfr_mappings_ke_unit / self.stellar_mass_unit

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_ke_path(self):
        return fs.join(self.cell_path, "ssfr_mappings_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_mappings_ke(self):
        return fs.is_file(self.cell_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_ssfr_mappings_ke_path", True, write=False)
    def ssfr_mappings_ke_data(self):

        """
        This function ...
        :return:
        """

        # Create the data
        return Data3D(self.ssfr_name, self.cell_x_coordinates, self.cell_y_coordinates, self.cell_z_coordinates, self.cell_ssfrs_mappings_ke,
                      length_unit=self.length_unit, unit=self.ssfr_mappings_ke_unit, description=self.ssfr_description, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # CELL SFR MAP
    # -----------------------------------------------------------------

    @property
    def sfr_name(self):
        return "SFR"

    # -----------------------------------------------------------------

    @property
    def sfr_description(self):
        return "Star formation rate"

    # -----------------------------------------------------------------

    @property
    def cell_sfr_salim_map_path(self):
        return fs.join(self.cell_path, "sfr_salim_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_salim_map(self):
        return fs.is_file(self.cell_sfr_salim_map_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "cell_sfr_salim_map_path", True, write=False)
    def sfr_salim_data_faceon_map(self):

        """
        map of the cell star formation rate
        :return:
        """

        return project_data(self.sfr_name, self.sfr_salim_data, self.faceon_projection, description=self.sfr_description)

    # -----------------------------------------------------------------

    @property
    def cell_sfr_ke_map_path(self):
        return fs.join(self.cell_path, "sfr_ke_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_ke_map(self):
        return fs.is_file(self.cell_sfr_ke_map_path)

    #-----------------------------------------------------------------

    @lazyfileproperty(Frame, "cell_sfr_ke_map_path", True, write=False)
    def sfr_ke_data_faceon_map(self):

        """
        This function ...
        :return:
        """

        return project_data(self.sfr_name, self.sfr_ke_data, self.faceon_projection, description=self.sfr_description)

    # -----------------------------------------------------------------

    @property
    def cell_sfr_mappings_ke_map_path(self):
        return fs.join(self.cell_path, "sfr_mappings_ke_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_mappings_ke_map(self):
        return fs.is_file(self.cell_sfr_mappings_ke_map_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "cell_sfr_mappings_ke_map_path", True, write=False)
    def sfr_mappings_ke_data_faceon_map(self):

        """
        This function ...
        :return:
        """

        return project_data(self.sfr_name, self.sfr_mappings_ke_data, self.faceon_projection, description=self.sfr_description)

    # -----------------------------------------------------------------
    # CELL STELLAR MASS MAP
    # -----------------------------------------------------------------

    @property
    def stellar_mass_name(self):
        return "Mstellar"

    # -----------------------------------------------------------------

    @property
    def stellar_mass_description(self):
        return "Stellar mass"

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "cell_mass_map_path", True, write=False)
    def stellar_mass_data_faceon_map(self):

        """
        map of the cell stellar mass
        :return:
        """

        return project_data(self.stellar_mass_name, self.stellar_mass_data, self.faceon_projection, description=self.stellar_mass_description)

    # -----------------------------------------------------------------
    # CELL SSFR MAP
    # -----------------------------------------------------------------

    @property
    def cell_ssfr_salim_map_path(self):
        return fs.join(self.cell_path, "ssfr_salim_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_salim_map(self):
        return fs.is_file(self.cell_ssfr_salim_map_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Image, "cell_ssfr_salim_map_path", True, write=False)
    def ssfr_salim_data_faceon_map(self):

        """
        map of the cell specific star formation rate
        """

        return project_data(self.ssfr_name, self.ssfr_salim_data, self.faceon_projection, description=self.ssfr_description, interpolate=True, as_image=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_data_faceon_frame(self):
        return self.ssfr_salim_data_faceon_map.primary

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_ke_map_path(self):
        return fs.join(self.cell_path, "ssfr_ke_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_ke_map(self):
        return fs.is_file(self.cell_ssfr_ke_map_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Image, "cell_ssfr_ke_map_path", True, write=False)
    def ssfr_ke_data_faceon_map(self):

        """
        This function ...
        :return:
        """

        return project_data(self.ssfr_name, self.ssfr_ke_data, self.faceon_projection, description=self.ssfr_description, interpolate=True, as_image=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_data_faceon_frame(self):
        return self.ssfr_ke_data_faceon_map.primary

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_map_path(self):
        return fs.join(self.cell_path, "ssfr_mappings_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_mappings_map(self):
        return fs.is_file(self.cell_ssfr_mappings_map_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "cell_ssfr_mappings_map_path", True, write=False)
    def ssfr_mappings_data_faceon_map(self):

        """
        This function ...
        :return:
        """

        return project_data(self.ssfr_name, self.ssfr_mappings_data, self.faceon_projection, description=self.ssfr_description, interpolate=True, as_image=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_data_faceon_frame(self):
        return self.ssfr_mappings_data_faceon_map.primary

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_ke_map_path(self):
        return fs.join(self.cell_path, "ssfr_mappings_ke_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_mappings_ke_map(self):
        return fs.is_file(self.cell_ssfr_mappings_ke_map_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "cell_ssfr_mappings_ke_map_path", True, write=False)
    def ssfr_mappings_ke_data_faceon_map(self):

        """
        This function ...
        :return:
        """

        return project_data(self.ssfr_name, self.ssfr_mappings_ke_data, self.faceon_projection, description=self.ssfr_description, interpolate=True, as_image=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_data_faceon_frame(self):
        return self.ssfr_mappings_ke_data_faceon_map.primary

    # -----------------------------------------------------------------
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
        if self.config.project: self.write_cell_maps()

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

    def write_projected_sfr(self):

        """
        This function ...
        :return:
        """

        # Salim
        self.write_projected_sfr_salim()

        # K&E
        self.write_projected_sfr_ke()

        # MAPPINGS
        self.write_projected_sfr_mappings()

        # MAPPINGS + K&E
        self.write_projected_sfr_mappings_ke()

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_salim_earth(self):
        return not self.has_projected_sfr_salim_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_salim_faceon(self):
        return not self.has_projected_sfr_salim_faceon

    # -----------------------------------------------------------------

    def write_projected_sfr_salim(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_projected_sfr_salim_earth: self.write_projected_sfr_salim_earth()

        # Faceon
        if self.do_write_projected_sfr_salim_faceon: self.write_projected_sfr_salim_faceon()

    # -----------------------------------------------------------------

    def write_projected_sfr_salim_earth(self):

        """
        Thisn function ...
        :return:
        """

        # Write
        self.sfr_salim_earth_map.saveto(self.projected_sfr_salim_earth_path)

    # -----------------------------------------------------------------

    def write_projected_sfr_salim_faceon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_salim_faceon_map.saveto(self.projected_sfr_salim_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_ke_earth(self):
        return not self.has_projected_sfr_ke_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_ke_faceon(self):
        return not self.has_projected_sfr_ke_faceon

    # -----------------------------------------------------------------

    def write_projected_sfr_ke(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_projected_sfr_ke_earth: self.write_projected_sfr_ke_earth()

        # Faceon
        if self.do_write_projected_sfr_ke_faceon: self.write_projected_sfr_ke_faceon()

    # -----------------------------------------------------------------

    def write_projected_sfr_ke_earth(self):

        """
        This function ...
        :return:
        """

        self.sfr_ke_earth_map.saveto(self.projected_sfr_ke_earth_path)

    # -----------------------------------------------------------------

    def write_projected_sfr_ke_faceon(self):

        """
        This function ...
        :return:
        """

        self.sfr_ke_faceon_map.saveto(self.projected_sfr_ke_faceon_path)

    # -----------------------------------------------------------------
    
    @property
    def do_write_projected_sfr_mappings_earth(self):
        return not self.has_projected_sfr_mappings_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_mappings_faceon(self):
        return not self.has_projected_sfr_mappings_faceon

    # -----------------------------------------------------------------

    def write_projected_sfr_mappings(self):
        
        """
        This function ...
        :return: 
        """
        
        # Earth
        if self.do_write_projected_sfr_mappings_earth: self.write_projected_sfr_mappings_earth()
        
        # Faceon
        if self.do_write_projected_sfr_mappings_faceon: self.write_projected_sfr_mappings_faceon()
        
    # -----------------------------------------------------------------

    def write_projected_sfr_mappings_earth(self):
        
        """
        This function ...
        :return: 
        """
        
        self.sfr_mappings_earth_map.saveto(self.projected_sfr_mappings_earth_path)
        
    # -----------------------------------------------------------------

    def write_projected_sfr_mappings_faceon(self):

        """
        This function ...
        :return: 
        """
        
        self.sfr_mappings_faceon_map.saveto(self.projected_sfr_mappings_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_mappings_ke_earth(self):
        return not self.has_projected_sfr_mappings_ke_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_mappings_ke_faceon(self):
        return not self.has_projected_sfr_mappings_ke_faceon

    # -----------------------------------------------------------------

    def write_projected_sfr_mappings_ke(self):

        """
        Thisf unction ...
        :return:
        """

        # Earth
        if self.do_write_projected_sfr_mappings_ke_earth: self.write_projected_sfr_mappings_ke_earth()

        # Faceon
        if self.do_write_projected_sfr_mappings_ke_faceon: self.write_projected_sfr_mappings_ke_faceon()

    # -----------------------------------------------------------------

    def write_projected_sfr_mappings_ke_earth(self):

        """
        This function ...
        :return:
        """

        self.sfr_mappings_ke_earth_map.saveto(self.projected_sfr_mappings_ke_earth_path)

    # -----------------------------------------------------------------

    def write_projected_sfr_mappings_ke_faceon(self):

        """
        This function ...
        :return:
        """

        self.sfr_mappings_ke_faceon_map.saveto(self.projected_sfr_mappings_ke_faceon_path)

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

    def write_projected_mass_earth(self):

        """
        This function ...
        :return:
        """

        # Write
        self.stellar_mass_earth_map.saveto(self.projected_mass_earth_path)

    # -----------------------------------------------------------------

    def write_projected_mass_faceon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.stellar_mass_faceon_map.saveto(self.projected_mass_faceon_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr(self):

        """
        This function ...
        :return:
        """

        # Salim
        self.write_projected_ssfr_salim()

        # K&E
        self.write_projected_ssfr_ke()

        # MAPPINGS
        self.write_projected_ssfr_mappings()

        # MAPPINGS + K&E
        self.write_projected_ssfr_mappings_ke()

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_salim_earth(self):
        return not self.has_projected_ssfr_salim_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_salim_faceon(self):
        return not self.has_projected_ssfr_salim_faceon

    # -----------------------------------------------------------------

    def write_projected_ssfr_salim(self):

        """
        Thisf unction ...
        :return:
        """

        # Earth
        if self.do_write_projected_ssfr_salim_earth: self.write_projected_ssfr_salim_earth()

        # Faceon
        if self.do_write_projected_ssfr_salim_faceon: self.write_projected_ssfr_salim_faceon()

    # -----------------------------------------------------------------

    def write_projected_ssfr_salim_earth(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_salim_earth_map.saveto(self.projected_ssfr_salim_earth_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_salim_faceon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_salim_faceon_map.saveto(self.projected_ssfr_salim_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_ke_earth(self):
        return not self.has_projected_ssfr_ke_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_ke_faceon(self):
        return not self.has_projected_ssfr_ke_faceon

    # -----------------------------------------------------------------

    def write_projected_ssfr_ke(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_projected_ssfr_ke_earth: self.write_projected_ssfr_ke_earth()

        # Faceon
        if self.do_write_projected_ssfr_ke_faceon: self.write_projected_ssfr_ke_faceon()

    # -----------------------------------------------------------------

    def write_projected_ssfr_ke_earth(self):

        """
        This function ...
        :return:
        """

        self.ssfr_ke_earth_map.saveto(self.projected_ssfr_ke_earth_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_ke_faceon(self):

        """
        This function ...
        :return:
        """

        self.ssfr_ke_faceon_map.saveto(self.projected_ssfr_ke_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_mappings_earth(self):
        return not self.has_projected_ssfr_mappings_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_mappings_faceon(self):
        return not self.has_projected_ssfr_mappings_faceon

    # -----------------------------------------------------------------

    def write_projected_ssfr_mappings(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_projected_ssfr_mappings_earth: self.write_projected_ssfr_mappings_earth()

        # Faceon
        if self.do_write_projected_ssfr_mappings_faceon: self.write_projected_ssfr_mappings_faceon()

    # -----------------------------------------------------------------

    def write_projected_ssfr_mappings_earth(self):

        """
        This function ...
        :return:
        """

        self.ssfr_mappings_earth_map.saveto(self.projected_ssfr_mappings_earth_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_mappings_faceon(self):

        """
        This function ...
        :return:
        """

        self.ssfr_mappings_faceon_map.saveto(self.projected_ssfr_mappings_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_mappings_ke_earth(self):
        return not self.has_projected_ssfr_mappings_ke_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_mappings_ke_faceon(self):
        return not self.has_projected_ssfr_mappings_ke_faceon

    # -----------------------------------------------------------------

    def write_projected_ssfr_mappings_ke(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_projected_ssfr_mappings_ke_earth: self.write_projected_ssfr_mappings_ke_earth()

        # Faceon
        if self.do_write_projected_ssfr_mappings_ke_faceon: self.write_projected_ssfr_mappings_ke_faceon()

    # -----------------------------------------------------------------

    def write_projected_ssfr_mappings_ke_earth(self):

        """
        This function ...
        :return:
        """

        self.ssfr_mappings_ke_earth_map.saveto(self.projected_ssfr_mappings_ke_earth_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_mappings_ke_faceon(self):

        """
        This function ...
        :return:
        """

        self.ssfr_mappings_ke_faceon_map.saveto(self.projected_ssfr_mappings_ke_faceon_path)

    # -----------------------------------------------------------------
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
    def do_write_cell_mass(self):
        return not self.has_cell_mass

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
        self.write_cell_sfr()

        # Stellar mass
        if self.do_write_cell_mass: self.write_cell_mass()

        # Specific star formation rate
        self.write_cell_ssfr()

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
    def do_write_cell_sfr_salim(self):
        return not self.has_cell_sfr_salim

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_ke(self):
        return not self.has_cell_sfr_ke

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_mappings(self):
        return not self.has_cell_sfr_mappings

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_mappings_ke(self):
        return not self.has_cell_sfr_mappings_ke

    # -----------------------------------------------------------------

    def write_cell_sfr(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Writing the cell star formation rate ...")

        # Salim
        if self.do_write_cell_sfr_salim: self.write_cell_sfr_salim()

        # K&E
        if self.do_write_cell_sfr_ke: self.write_cell_sfr_ke()

        # MAPPINGS
        if self.do_write_cell_sfr_mappings: self.write_cell_sfr_mappings()

        # MAPPINGS + K&E
        if self.do_write_cell_sfr_mappings_ke: self.write_cell_sfr_mappings_ke()

    # -----------------------------------------------------------------

    def write_cell_sfr_salim(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_salim_data.saveto(self.cell_sfr_salim_path)

    # -----------------------------------------------------------------

    def write_cell_sfr_ke(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_ke_data.saveto(self.cell_sfr_ke_path)

    # -----------------------------------------------------------------

    def write_cell_sfr_mappings(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_mappings_data.saveto(self.cell_sfr_mappings_path)

    # -----------------------------------------------------------------

    def write_cell_sfr_mappings_ke(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_mappings_ke_data.saveto(self.cell_sfr_mappings_ke_path)

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
    def do_write_cell_ssfr_salim(self):
        return not self.has_cell_ssfr_salim

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_ke(self):
        return not self.has_cell_ssfr_ke

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_mappings(self):
        return not self.has_cell_ssfr_mappings

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_mappings_ke(self):
        return not self.has_cell_ssfr_mappings_ke

    # -----------------------------------------------------------------

    def write_cell_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cell specific star formation rate ...")

        # Salim
        if self.do_write_cell_ssfr_salim: self.write_cell_ssfr_salim()

        # K&E
        if self.do_write_cell_ssfr_ke: self.write_cell_ssfr_ke()

        # MAPPINGS
        if self.do_write_cell_ssfr_mappings: self.write_cell_ssfr_mappings()

        # MAPPINGS + K&E
        if self.do_write_cell_ssfr_mappings_ke: self.write_cell_ssfr_mappings_ke()

    # -----------------------------------------------------------------

    def write_cell_ssfr_salim(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_salim_data.saveto(self.cell_ssfr_salim_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_ke(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_ke_data.saveto(self.cell_ssfr_ke_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_mappings(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_mappings_data.saveto(self.cell_ssfr_mappings_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_mappings_ke(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_mappings_ke_data.saveto(self.cell_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @property
    def do_write_cell_mass_map(self):
        return not self.has_cell_mass_map

    # -----------------------------------------------------------------

    def write_cell_maps(self):

        """
        Thsi function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps of the cell star formation rates ...")

        # Star formation rate
        # NO: MAP OF SFR DATA IS INCORRECT BECAUSE IT TAKES THE AVERAGE OF CELL VALUES
        #self.write_cell_sfr_maps()

        # Stellar mass
        # NO: MAP OF STELLAR MASS DATA IS INCORRECT BECAUSE IT TAKES THE AVERAGE OF CELL VALUES
        #if self.do_write_cell_mass_map: self.write_cell_mass_map()

        # Specific star formation rate
        self.write_cell_ssfr_maps()

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_salim_map(self):
        return not self.has_cell_sfr_salim_map

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_ke_map(self):
        return not self.has_cell_sfr_ke_map

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_mappings_ke_map(self):
        return not self.has_cell_sfr_mappings_ke_map

    # -----------------------------------------------------------------

    def write_cell_sfr_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps of the cell star formation rate ...")

        # Salim
        if self.do_write_cell_sfr_salim_map: self.write_cell_sfr_salim_map()

        # K&E
        if self.do_write_cell_sfr_ke_map: self.write_cell_sfr_ke_map()

        # MAPPINGS + K&E
        if self.do_write_cell_sfr_mappings_ke_map: self.write_cell_sfr_mappings_ke_map()

    # -----------------------------------------------------------------

    def write_cell_sfr_salim_map(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_salim_data_faceon_map.saveto(self.cell_sfr_salim_map_path)

    # -----------------------------------------------------------------

    def write_cell_sfr_ke_map(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_ke_data_faceon_map.saveto(self.cell_sfr_ke_map_path)

    # -----------------------------------------------------------------

    def write_cell_sfr_mappings_ke_map(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_mappings_ke_data_faceon_map.saveto(self.cell_sfr_mappings_ke_map_path)

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
    def do_write_cell_ssfr_salim_map(self):
        return not self.has_cell_ssfr_salim_map

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_ke_map(self):
        return not self.has_cell_ssfr_ke_map

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_mappings_map(self):
        return not self.has_cell_ssfr_mappings_map

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_mappings_ke_map(self):
        return not self.has_cell_ssfr_mappings_ke_map

    # -----------------------------------------------------------------

    def write_cell_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps of the cell specific star formation rate ...")

        # Salim
        if self.do_write_cell_ssfr_salim_map: self.write_cell_ssfr_salim_map()

        # K&E
        if self.do_write_cell_ssfr_ke_map: self.write_cell_ssfr_ke_map()

        # MAPPINGS
        if self.do_write_cell_ssfr_mappings_map: self.write_cell_ssfr_mappings_map()

        # MAPPINGS + K&E
        if self.do_write_cell_ssfr_mappings_ke_map: self.write_cell_ssfr_mappings_ke_map()

    # -----------------------------------------------------------------

    def write_cell_ssfr_salim_map(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_salim_data_faceon_map.saveto(self.cell_ssfr_salim_map_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_ke_map(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_ke_data_faceon_map.saveto(self.cell_ssfr_ke_map_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_mappings_map(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_mappings_data_faceon_map.saveto(self.cell_ssfr_mappings_map_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_mappings_ke_map(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_mappings_ke_data_faceon_map.saveto(self.cell_ssfr_mappings_ke_map_path)

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
        if self.config.project: self.plot_cell_maps()

    # -----------------------------------------------------------------

    @property
    def sfr_limits(self):
        #return (1e-6*self.sfr_msun_yr, 1e-3*self.sfr_msun_yr,) # not per solid angle
        return (1e-8*self.sfr_msun_yr, 1e-5*self.sfr_msun_yr,)

    # -----------------------------------------------------------------

    @property
    def ssfr_limits(self):
        #return (1e-13, 1e-9,)
        return (1e-14, 1e-10,)

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

    def plot_sfr_map(self, frame, path):

        """
        This function ...
        :param frame:
        :param path:
        :return:
        """

        # Create frame per pixel (angular) area
        converted = frame / frame.pixel_solid_angle.to("arcsec2").value

        # Plot
        plot_map(converted, path=path, cmap="inferno", colorbar=True, interval=self.sfr_limits, scale="log", background_color="black", title="Star formation rate (Msun / yr / arcsec2)")

    # -----------------------------------------------------------------

    def plot_mass_map(self, frame, path):

        """
        Thisn function ...
        :param frame:
        :param path:
        :return:
        """

        plot_map(frame, path=path, cmap="inferno", colorbar=True, scale="log", background_color="black", title="Stellar mass")

    # -----------------------------------------------------------------

    def plot_ssfr_map(self, frame, path):

        """
        This function ...
        :param frame:
        :param path:
        :return:
        """

        # Plot
        plot_map(frame, path=path, cmap="inferno", colorbar=True, interval=self.ssfr_limits, scale="log", background_color="black", title="Specific star formation rate (1/yr)") # bg color is for invalid pixels (NaN)

    # -----------------------------------------------------------------

    def plot_projected_sfr(self):

        """
        This function ...
        :return:
        """

        # Salim
        self.plot_projected_sfr_salim()

        # K&E
        self.plot_projected_sfr_ke()

        # MAPPINGS
        self.plot_projected_sfr_mappings()

        # MAPPINGS + K&E
        self.plot_projected_sfr_mappings_ke()

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_salim_earth_map(self):
        return not self.has_projected_sfr_salim_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_salim_faceon_map(self):
        return not self.has_projected_sfr_salim_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_sfr_salim(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_sfr_salim_earth_map: self.plot_projected_sfr_salim_earth()

        # Faceon
        if self.do_plot_projected_sfr_salim_faceon_map: self.plot_projected_sfr_salim_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_sfr_salim_earth_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_salim_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_salim_earth_map_plot(self):
        return fs.is_file(self.projected_sfr_salim_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_salim_earth(self):

        """
        This function ...
        :return:
        """

        # Plot
        self.plot_sfr_map(self.sfr_salim_earth_map, self.projected_sfr_salim_earth_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def projected_sfr_salim_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_salim_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_salim_faceon_map_plot(self):
        return fs.is_file(self.projected_sfr_salim_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_salim_faceon(self):

        """
        This function ...
        :return:
        """

        # Plot
        self.plot_sfr_map(self.sfr_salim_faceon_map, self.projected_sfr_salim_faceon_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_ke_earth(self):
        return not self.has_projected_sfr_ke_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_ke_faceon(self):
        return not self.has_projected_sfr_ke_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_sfr_ke(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_sfr_ke_earth: self.plot_projected_sfr_ke_earth()

        # Faceon
        if self.do_plot_projected_sfr_ke_faceon: self.plot_projected_sfr_ke_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_sfr_ke_earth_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_ke_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_ke_earth_map_plot(self):
        return fs.is_file(self.projected_sfr_ke_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_ke_earth(self):

        """
        This function ...
        :return:
        """

        self.plot_sfr_map(self.sfr_ke_earth_map, self.projected_sfr_ke_earth_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def projected_sfr_ke_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_ke_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_ke_faceon_map_plot(self):
        return fs.is_file(self.projected_sfr_ke_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_ke_faceon(self):

        """
        This function ...
        :return:
        """

        self.plot_sfr_map(self.sfr_ke_faceon_map, self.projected_sfr_ke_faceon_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_mappings_earth(self):
        return not self.has_projected_sfr_mappings_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_mappings_faceon(self):
        return not self.has_projected_sfr_mappings_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_sfr_mappings(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_sfr_mappings_earth: self.plot_projected_sfr_mappings_earth()

        # Faceon
        if self.do_plot_projected_sfr_mappings_faceon: self.plot_projected_sfr_mappings_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_sfr_mappings_earth_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_mappings_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_mappings_earth_map_plot(self):
        return fs.is_file(self.projected_sfr_mappings_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_mappings_earth(self):

        """
        This function ...
        :return:
        """

        self.plot_sfr_map(self.sfr_mappings_earth_map, self.projected_sfr_mappings_earth_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def projected_sfr_mappings_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_mappings_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_mappings_faceon_map_plot(self):
        return fs.is_file(self.projected_sfr_mappings_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_mappings_faceon(self):

        """
        This function ...
        :return:
        """

        self.plot_sfr_map(self.sfr_mappings_faceon_map, self.projected_sfr_mappings_faceon_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_mappings_ke_earth(self):
        return not self.has_projected_sfr_mappings_ke_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_mappings_ke_faceon(self):
        return not self.has_projected_sfr_mappings_ke_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_sfr_mappings_ke(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_sfr_mappings_ke_earth: self.plot_projected_sfr_mappings_ke_earth()

        # Faceon
        if self.do_plot_projected_sfr_mappings_ke_faceon: self.plot_projected_sfr_mappings_ke_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_sfr_mappings_ke_earth_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_mappings_ke_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_mappings_ke_earth_map_plot(self):
        return fs.is_file(self.projected_sfr_mappings_ke_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_mappings_ke_earth(self):

        """
        Thisf unction ...
        :return:
        """

        self.plot_sfr_map(self.sfr_mappings_ke_earth_map, self.projected_sfr_mappings_ke_earth_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def projected_sfr_mappings_ke_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_mappings_ke_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_mappings_ke_faceon_map_plot(self):
        return fs.is_file(self.projected_sfr_mappings_ke_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_mappings_ke_faceon(self):

        """
        This function ...
        :return:
        """

        self.plot_sfr_map(self.sfr_mappings_ke_faceon_map, self.projected_sfr_mappings_ke_faceon_map_plot_path)

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
        self.plot_mass_map(self.stellar_mass_earth_map, self.projected_mass_earth_map_plot_path)

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
        self.plot_mass_map(self.stellar_mass_faceon_map, self.projected_mass_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_ssfr(self):

        """
        This function ...
        :return:
        """

        # Salim
        self.plot_projected_ssfr_salim()

        # K&E
        self.plot_projected_ssfr_ke()

        # MAPPINGS
        self.plot_projected_ssfr_mappings()

        # MAPPINGS + K&E
        self.plot_projected_ssfr_mappings_ke()

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_ssfr_salim_earth(self):
        return not self.has_projected_ssfr_salim_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_ssfr_salim_faceon(self):
        return not self.has_projected_ssfr_salim_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_ssfr_salim(self):

        """
        Thisn function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_ssfr_salim_earth: self.plot_projected_ssfr_salim_earth()

        # Faceon
        if self.do_plot_projected_ssfr_salim_faceon: self.plot_projected_ssfr_salim_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_salim_earth_map_plot_path(self):
        return fs.join(self.projected_path, "ssfr_salim_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_salim_earth_map_plot(self):
        return fs.is_file(self.projected_ssfr_salim_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_ssfr_salim_earth(self):

        """
        This function ...
        :return:
        """

        # Plot
        self.plot_ssfr_map(self.ssfr_salim_earth_map, self.projected_ssfr_salim_earth_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_salim_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "ssfr_salim_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_salim_faceon_map_plot(self):
        return fs.is_file(self.projected_ssfr_salim_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_ssfr_salim_faceon(self):

        """
        This function ...
        :return:
        """

        # Plot
        self.plot_ssfr_map(self.ssfr_salim_faceon_map, self.projected_ssfr_salim_faceon_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_ssfr_ke_earth(self):
        return not self.has_projected_ssfr_ke_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_ssfr_ke_faceon(self):
        return not self.has_projected_ssfr_ke_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_ssfr_ke(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_ssfr_ke_earth: self.plot_projected_ssfr_ke_earth()

        # Faceon
        if self.do_plot_projected_ssfr_ke_faceon: self.plot_projected_ssfr_ke_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_ke_earth_map_plot_path(self):
        return fs.join(self.projected_path, "ssfr_ke_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_ke_earth_map_plot(self):
        return fs.is_file(self.projected_ssfr_ke_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_ssfr_ke_earth(self):

        """
        Thisn function ...
        :return:
        """

        # Plot
        self.plot_ssfr_map(self.ssfr_ke_earth_map, self.projected_ssfr_ke_earth_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_ke_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "ssfr_ke_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_ke_faceon_map_plot(self):
        return fs.is_file(self.projected_ssfr_ke_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_ssfr_ke_faceon(self):

        """
        This function ...
        :return:
        """

        # Plot
        self.plot_ssfr_map(self.ssfr_ke_faceon_map, self.projected_ssfr_ke_faceon_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_ssfr_mappings_earth(self):
        return not self.has_projected_ssfr_mappings_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_ssfr_mappings_faceon(self):
        return not self.has_projected_ssfr_mappings_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_ssfr_mappings(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_ssfr_mappings_earth: self.plot_projected_ssfr_mappings_earth()

        # Faceon
        if self.do_plot_projected_ssfr_mappings_faceon: self.plot_projected_ssfr_mappings_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_mappings_earth_map_plot_path(self):
        return fs.join(self.projected_path, "ssfr_mappings_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_mappings_earth_map_plot(self):
        return fs.is_file(self.projected_ssfr_mappings_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_ssfr_mappings_earth(self):

        """
        This function ...
        :return:
        """

        # Plot
        self.plot_ssfr_map(self.ssfr_mappings_earth_map, self.projected_ssfr_mappings_earth_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_mappings_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "ssfr_mappings_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_mappings_faceon_map_plot(self):
        return fs.is_file(self.projected_ssfr_mappings_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_ssfr_mappings_faceon(self):

        """
        This function ...
        :return:
        """

        # Plot
        self.plot_ssfr_map(self.ssfr_mappings_faceon_map, self.projected_ssfr_mappings_faceon_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_ssfr_mappings_ke_earth(self):
        return not self.has_projected_ssfr_mappings_ke_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_ssfr_mappings_ke_faceon(self):
        return not self.has_projected_ssfr_mappings_ke_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_ssfr_mappings_ke(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_ssfr_mappings_ke_earth: self.plot_projected_ssfr_mappings_ke_earth()

        # Faceon
        if self.do_plot_projected_ssfr_mappings_ke_faceon: self.plot_projected_ssfr_mappings_ke_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_mappings_ke_earth_map_plot_path(self):
        return fs.join(self.projected_path, "ssfr_mappings_ke_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_mappings_ke_earth_map_plot(self):
        return fs.is_file(self.projected_ssfr_mappings_ke_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_ssfr_mappings_ke_earth(self):

        """
        Thisn function ...
        :return:
        """

        # Plot
        self.plot_ssfr_map(self.ssfr_mappings_ke_earth_map, self.projected_ssfr_mappings_ke_earth_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_mappings_ke_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "ssfr_mappings_ke_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_mappings_ke_faceon_map_plot(self):
        return fs.is_file(self.projected_ssfr_mappings_ke_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_ssfr_mappings_ke_faceon(self):

        """
        This function ...
        :return:
        """

        self.plot_ssfr_map(self.ssfr_mappings_ke_faceon_map, self.projected_ssfr_mappings_ke_faceon_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_cell_mass_map(self):
        return not self.has_cell_mass_map_plot

    # -----------------------------------------------------------------

    def plot_cell_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the cell star formation rates ...")

        # Star formation rate
        # NO: MAP OF STELLAR MASS DATA IS INCORRECT BECAUSE IT TAKES THE AVERAGE OF CELL VALUES
        #self.plot_cell_sfr_maps()

        # Stellar mass
        # NO: MAP OF STELLAR MASS DATA IS INCORRECT BECAUSE IT TAKES THE AVERAGE OF CELL VALUES
        #if self.do_plot_cell_mass_map: self.plot_cell_mass_map()

        # Specific star formation rate
        self.plot_cell_ssfr_maps()

    # -----------------------------------------------------------------

    @property
    def do_plot_cell_sfr_salim_map(self):
        return not self.has_cell_sfr_salim_map_plot

    #-----------------------------------------------------------------

    @property
    def do_plot_cell_sfr_ke_map(self):
        return not self.has_cell_sfr_ke_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_cell_sfr_mappings_ke_map(self):
        return not self.has_cell_sfr_mappings_ke_map_plot

    # -----------------------------------------------------------------

    def plot_cell_sfr_maps(self):

        """
        This function ...
        :return:
        """

        # Salim
        if self.do_plot_cell_sfr_salim_map: self.plot_cell_sfr_salim_map()

        # K&E
        if self.do_plot_cell_sfr_ke_map: self.plot_cell_sfr_ke_map()

        # MAPPINGS + K&E
        if self.do_plot_cell_sfr_mappings_ke_map: self.plot_cell_sfr_mappings_ke_map()

    # -----------------------------------------------------------------

    @property
    def cell_sfr_salim_map_plot_path(self):
        return fs.join(self.cell_path, "sfr_salim_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_salim_map_plot(self):
        return fs.is_file(self.cell_sfr_salim_map_plot_path)

    # -----------------------------------------------------------------

    def plot_cell_sfr_salim_map(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_map(self.sfr_salim_data_faceon_map, path=self.cell_sfr_salim_map_plot_path, cmap="inferno", colorbar=True)

    # -----------------------------------------------------------------

    @property
    def cell_sfr_ke_map_plot_path(self):
        return fs.join(self.cell_path, "sfr_ke_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_ke_map_plot(self):
        return fs.is_file(self.cell_sfr_ke_map_plot_path)

    # -----------------------------------------------------------------

    def plot_cell_sfr_ke_map(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_map(self.sfr_ke_data_faceon_map, path=self.cell_sfr_ke_map_plot_path, cmap="inferno", colorbar=True)

    # -----------------------------------------------------------------

    @property
    def cell_sfr_mappings_ke_map_plot_path(self):
        return fs.join(self.cell_path, "sfr_mappings_ke_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_mappings_ke_map_plot(self):
        return fs.is_file(self.cell_sfr_mappings_ke_map_plot_path)

    # -----------------------------------------------------------------

    def plot_cell_sfr_mappings_ke_map(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_map(self.sfr_mappings_ke_data_faceon_map, path=self.cell_sfr_mappings_ke_map_plot_path, cmap="inferno", colorbar=True)

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
    def do_plot_cell_ssfr_salim_map(self):
        return not self.has_cell_ssfr_salim_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_cell_ssfr_ke_map(self):
        return not self.has_cell_ssfr_ke_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_cell_ssfr_mappings_map(self):
        return not self.has_cell_ssfr_mappings_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_cell_ssfr_mappings_ke_map(self):
        return not self.has_cell_ssfr_mappings_ke_map_plot

    # -----------------------------------------------------------------

    def plot_cell_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        # Salim
        if self.do_plot_cell_ssfr_salim_map: self.plot_cell_ssfr_salim_map()

        # K&E
        if self.do_plot_cell_ssfr_ke_map: self.plot_cell_ssfr_ke_map()

        # MAPPINGS
        if self.do_plot_cell_ssfr_mappings_map: self.plot_cell_ssfr_mappings_map()

        # MAPPINGS + K&E
        if self.do_plot_cell_ssfr_mappings_ke_map: self.plot_cell_ssfr_mappings_ke_map()

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_salim_map_plot_path(self):
        return fs.join(self.cell_path, "ssfr_salim_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_salim_map_plot(self):
        return fs.is_file(self.cell_ssfr_salim_map_plot_path)

    # -----------------------------------------------------------------

    def plot_cell_ssfr_salim_map(self):

        """
        This function ...
        :return:
        """

        # Plot
        #self.plot_ssfr_map(self.ssfr_salim_data_faceon_map, self.cell_ssfr_salim_map_plot_path)
        self.plot_ssfr_map(self.ssfr_salim_data_faceon_frame, self.cell_ssfr_salim_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_ke_map_plot_path(self):
        return fs.join(self.cell_path, "ssfr_ke_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_ke_map_plot(self):
        return fs.is_file(self.cell_ssfr_ke_map_plot_path)

    # -----------------------------------------------------------------

    def plot_cell_ssfr_ke_map(self):

        """
        This function ...
        :return:
        """

        # Plot
        #self.plot_ssfr_map(self.ssfr_ke_data_faceon_map, self.cell_ssfr_ke_map_plot_path)
        self.plot_ssfr_map(self.ssfr_ke_data_faceon_frame, self.cell_ssfr_ke_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_map_plot_path(self):
        return fs.join(self.cell_path, "ssfr_mappings_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_mappings_map_plot(self):
        return fs.is_file(self.cell_ssfr_mappings_map_plot_path)

    # -----------------------------------------------------------------

    def plot_cell_ssfr_mappings_map(self):

        """
        This function ...
        :return:
        """

        # Plot
        #self.plot_ssfr_map(self.ssfr_mappings_data_faceon_map, self.cell_ssfr_mappings_map_plot_path)
        self.plot_ssfr_map(self.ssfr_mappings_data_faceon_frame, self.cell_ssfr_mappings_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_ke_map_plot_path(self):
        return fs.join(self.cell_path, "ssfr_mappings_ke_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_mappings_ke_map_plot(self):
        return fs.is_file(self.cell_ssfr_mappings_ke_map_plot_path)

    # -----------------------------------------------------------------

    def plot_cell_ssfr_mappings_ke_map(self):

        """
        This function ...
        :return:
        """

        # Plot
        #self.plot_ssfr_map(self.ssfr_mappings_ke_data_faceon_map, self.cell_ssfr_mappings_ke_map_plot_path)
        self.plot_ssfr_map(self.ssfr_mappings_ke_data_faceon_frame, self.cell_ssfr_mappings_ke_map_plot_path)

# -----------------------------------------------------------------
