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
from .component import AnalysisRunComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.core.frame import Frame
from ...magic.core.image import Image
from ...core.tools.utils import lazyproperty, lazyfileproperty
from ..core.data import Data3D, SpectralData3D
from ..projection.data import project_data
from ..core.model import oliver_stellar_mass, salim_fuv_to_sfr, kennicutt_evans_fuv_to_sfr, kennicutt_tir_to_sfr, calzetti_24um_to_sfr
from ...core.units.parsing import parse_unit as u
from ...magic.tools.plotting import plot_map
from ...magic.core.list import uniformize
from ...core.filter.filter import parse_filter
from ...magic.tools.colours import make_colour_map
from ...core.basics.configuration import open_box
from ...core.data.sed import SED

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

    @lazyproperty
    def specific_luminosity_unit(self):
        return u("W/micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def frequency_luminosity_unit(self):
        return u("W/Hz")

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
    def projected_sfr_salim_earth_corrected_path(self):
        return fs.join(self.projected_path, "sfr_salim_earth_corrected.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_salim_earth_corrected(self):
        return fs.is_file(self.projected_sfr_salim_earth_corrected_path)

    # -----------------------------------------------------------------

    @property
    def intrinsic_fuv_luminosity_map_corrected_earth(self):
        fuv_lum, absorbed = uniformize(self.model.intrinsic_fuv_luminosity_map_earth, self.internal_absorbed_fuv_luminosity_map_earth, convolve=False, wavelength=self.fuv_wavelength)
        return fuv_lum + absorbed

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_salim_earth_corrected_path", True, write=False)
    def sfr_salim_earth_map_corrected(self):

        """
        This function ...
        :return:
        """

        return salim_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_corrected_earth, distance=self.galaxy_distance)

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

    @property
    def projected_sfr_salim_faceon_corrected_path(self):
        return fs.join(self.projected_path, "sfr_salim_faceon_corrected.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_salim_faceon_corrected(self):
        return fs.is_file(self.projected_sfr_salim_faceon_corrected_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity_map_corrected_faceon(self):
        fuv_lum, absorbed = uniformize(self.model.intrinsic_fuv_luminosity_map_faceon, self.internal_absorbed_fuv_luminosity_map_faceon, convolve=False, wavelength=self.fuv_wavelength)
        return fuv_lum + absorbed

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_salim_faceon_corrected_path", True, write=False)
    def sfr_salim_faceon_map_corrected(self):

        """
        This function ...
        :return:
        """

        return salim_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_corrected_faceon, distance=self.galaxy_distance)

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
    def projected_sfr_ke_earth_corrected_path(self):
        return fs.join(self.projected_path, "sfr_ke_earth_corrected.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_ke_earth_corrected(self):
        return fs.is_file(self.projected_sfr_ke_earth_corrected_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_ke_earth_corrected_path", True, write=False)
    def sfr_ke_earth_map_corrected(self):

        """
        This function ...
        :return:
        """

        return kennicutt_evans_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_corrected_earth, distance=self.galaxy_distance)

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

    @property
    def projected_sfr_ke_faceon_corrected_path(self):
        return fs.join(self.projected_path, "sfr_ke_faceon_corrected.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_ke_faceon_corrected(self):
        return fs.is_file(self.projected_sfr_ke_faceon_corrected_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_ke_faceon_corrected_path", True, write=False)
    def sfr_ke_faceon_map_corrected(self):

        """
        Thisk function ...
        :return:
        """

        return kennicutt_evans_fuv_to_sfr(self.intrinsic_fuv_luminosity_map_corrected_faceon, distance=self.galaxy_distance)

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
    #   5. TIR
    # -----------------------------------------------------------------

    @property
    def projected_sfr_tir_earth_path(self):
        return fs.join(self.projected_path, "sfr_tir_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_tir_earth(self):
        return fs.is_file(self.projected_sfr_tir_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_tir_earth_path", True, write=False)
    def sfr_tir_earth_map(self):

        """
        This function ...
        :return:
        """

        return self.model.total_star_formation_rate_map_earth_tir

    # -----------------------------------------------------------------

    @property
    def projected_sfr_tir_faceon_path(self):
        return fs.join(self.projected_path, "sfr_tir_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_tir_faceon(self):
        return fs.is_file(self.projected_sfr_tir_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_tir_faceon_path", True, write=False)
    def sfr_tir_faceon_map(self):

        """
        This function ...
        :return:
        """

        return self.model.total_star_formation_rate_map_faceon_tir

    # -----------------------------------------------------------------
    #   6. 24 micron
    # -----------------------------------------------------------------

    @property
    def projected_sfr_24um_earth_path(self):
        return fs.join(self.projected_path, "sfr_24um_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_24um_earth(self):
        return fs.is_file(self.projected_sfr_24um_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_24um_earth_path", True, write=False)
    def sfr_24um_earth_map(self):

        """
        This function ...
        :return:
        """

        return self.model.total_star_formation_rate_map_earth_24um

    # -----------------------------------------------------------------

    @property
    def projected_sfr_24um_faceon_path(self):
        return fs.join(self.projected_path, "sfr_24um_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_24um_faceon(self):
        return fs.is_file(self.projected_sfr_24um_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_sfr_24um_faceon_path", True, write=False)
    def sfr_24um_faceon_map(self):

        """
        This function ...
        :return:
        """

        return self.model.total_star_formation_rate_map_faceon_24um

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
    def projected_ssfr_salim_earth_corrected_path(self):
        return fs.join(self.projected_path, "ssfr_salim_earth_corrected.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_salim_earth_corrected(self):
        return fs.is_file(self.projected_ssfr_salim_earth_corrected_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_salim_earth_corrected_path", True, write=False)
    def ssfr_salim_earth_map_corrected(self):

        """
        This function ...
        :return:
        """

        sfr, old_stellar_mass = uniformize(self.sfr_salim_earth_map_corrected, self.model.old_stellar_mass_map_earth, convert=False)
        return sfr / old_stellar_mass

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

    @property
    def projected_ssfr_salim_faceon_corrected_path(self):
        return fs.join(self.projected_path, "ssfr_salim_faceon_corrected.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_salim_faceon_corrected(self):
        return fs.is_file(self.projected_ssfr_salim_faceon_corrected_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_salim_faceon_corrected_path", True, write=False)
    def ssfr_salim_faceon_map_corrected(self):

        """
        This function ...
        :return:
        """

        sfr, old_stellar_mass = uniformize(self.sfr_salim_faceon_map_corrected, self.model.old_stellar_mass_map_faceon, convert=False)
        return sfr / old_stellar_mass

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
    def projected_ssfr_ke_earth_corrected_path(self):
        return fs.join(self.projected_path, "ssfr_ke_earth_corrected.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_ke_earth_corrected(self):
        return fs.is_file(self.projected_ssfr_ke_earth_corrected_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_ke_earth_corrected_path", True, write=False)
    def ssfr_ke_earth_map_corrected(self):

        """
        This function ...
        :return:
        """

        sfr, old_stellar_mass = uniformize(self.sfr_ke_earth_map_corrected, self.model.old_stellar_mass_map_earth, convert=False)
        return sfr / old_stellar_mass

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

    @property
    def projected_ssfr_ke_faceon_corrected_path(self):
        return fs.join(self.projected_path, "ssfr_ke_faceon_corrected.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_ke_faceon_corrected(self):
        return fs.is_file(self.projected_ssfr_ke_faceon_corrected_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_ke_faceon_corrected_path", True, write=False)
    def ssfr_ke_faceon_map_corrected(self):

        """
        This function ...
        :return:
        """

        sfr, old_stellar_mass = uniformize(self.sfr_ke_faceon_map_corrected, self.model.old_stellar_mass_map_faceon, convert=False)
        return sfr / old_stellar_mass

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
    # FUV-H
    # -----------------------------------------------------------------

    @property
    def projected_ssfr_fuv_h_earth_path(self):
        return fs.join(self.projected_path, "ssfr_fuv_h_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_fuv_h_earth(self):
        return fs.is_file(self.projected_ssfr_fuv_h_earth_path)

    # -----------------------------------------------------------------

    @property
    def total_observed_earth_cube(self):
        return self.model.total_observed_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def total_observed_faceon_cube(self):
        return self.model.total_observed_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def total_observed_edgeon_cube(self):
        return self.model.total_observed_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def total_fuv_luminosity_map_earth(self):
        return self.total_observed_earth_cube.get_frame_for_wavelength(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_h_luminosity_map_earth(self):
        return self.total_observed_earth_cube.get_frame_for_wavelength(self.h_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_fuv_h_earth_path", True, write=False)
    def ssfr_fuv_h_earth_map(self):

        """
        This function ...
        :return:
        """

        # Get maps
        fuv = self.total_fuv_luminosity_map_earth
        h = self.total_h_luminosity_map_earth

        # Return the colour map
        return make_colour_map(fuv, h)

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_fuv_h_faceon_path(self):
        return fs.join(self.projected_path, "ssfr_fuv_h_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_fuv_h_faceon(self):
        return fs.is_file(self.projected_ssfr_fuv_h_faceon_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_fuv_luminosity_map_faceon(self):
        return self.total_observed_faceon_cube.get_frame_for_wavelength(self.fuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_h_luminosity_map_faceon(self):
        return self.total_observed_faceon_cube.get_frame_for_wavelength(self.h_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_fuv_h_faceon_path", True, write=False)
    def ssfr_fuv_h_faceon_map(self):

        """
        This function ...
        :return:
        """

        # Get maps
        fuv = self.total_fuv_luminosity_map_faceon
        h = self.total_h_luminosity_map_faceon

        # Return the colour map
        return make_colour_map(fuv, h)

    # -----------------------------------------------------------------
    # FUV-r
    # -----------------------------------------------------------------

    @property
    def projected_ssfr_fuv_r_earth_path(self):
        return fs.join(self.projected_path, "ssfr_fuv_r_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_fuv_r_earth(self):
        return fs.is_file(self.projected_ssfr_fuv_r_earth_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_r_luminosity_map_earth(self):
        return self.total_observed_earth_cube.get_frame_for_wavelength(self.r_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_fuv_r_earth_path", True, write=False)
    def ssfr_fuv_r_earth_map(self):

        """
        This function ...
        :return:
        """

        # Get maps
        fuv = self.total_fuv_luminosity_map_earth
        r = self.total_r_luminosity_map_earth

        # Return the colour map
        return make_colour_map(fuv, r)

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_fuv_r_faceon_path(self):
        return fs.join(self.projected_path, "ssfr_fuv_r_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_fuv_r_faceon(self):
        return fs.is_file(self.projected_ssfr_fuv_r_faceon_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_r_luminosity_map_faceon(self):
        return self.total_observed_faceon_cube.get_frame_for_wavelength(self.r_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_fuv_r_faceon_path", True, write=False)
    def ssfr_fuv_r_faceon_map(self):

        """
        This function ...
        :return:
        """

        # Get maps
        fuv = self.total_fuv_luminosity_map_faceon
        r = self.total_r_luminosity_map_faceon

        # Return the colour map
        return make_colour_map(fuv, r)

    # -----------------------------------------------------------------
    # NUV-H
    # -----------------------------------------------------------------

    @lazyproperty
    def total_nuv_luminosity_map_earth(self):
        return self.total_observed_earth_cube.get_frame_for_wavelength(self.nuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_nuv_h_earth_path(self):
        return fs.join(self.projected_path, "ssfr_nuv_h_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_nuv_h_earth(self):
        return fs.is_file(self.projected_ssfr_nuv_h_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_nuv_h_earth_path", True, write=False)
    def ssfr_nuv_h_earth_map(self):

        """
        This function ...
        :return:
        """

        # Get maps
        nuv = self.total_nuv_luminosity_map_earth
        h = self.total_h_luminosity_map_earth

        # Return the colour map
        return make_colour_map(nuv, h)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_nuv_luminosity_map_faceon(self):
        return self.total_observed_faceon_cube.get_frame_for_wavelength(self.nuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_nuv_h_faceon_path(self):
        return fs.join(self.projected_path, "ssfr_nuv_h_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_nuv_h_faceon(self):
        return fs.is_file(self.projected_ssfr_nuv_h_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_nuv_h_faceon_path", True, write=False)
    def ssfr_nuv_h_faceon_map(self):

        """
        This function ...
        :return:
        """

        # Get maps
        nuv = self.total_nuv_luminosity_map_faceon
        h = self.total_h_luminosity_map_faceon

        # Return
        return make_colour_map(nuv, h)

    # -----------------------------------------------------------------
    # NUV-r
    # -----------------------------------------------------------------

    @property
    def projected_ssfr_nuv_r_earth_path(self):
        return fs.join(self.projected_path, "ssfr_nuv_r_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_nuv_r_earth(self):
        return fs.is_file(self.projected_ssfr_nuv_r_earth_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_nuv_r_earth_path", True, write=False)
    def ssfr_nuv_r_earth_map(self):

        """
        This function ...
        :return:
        """

        # Get maps
        nuv = self.total_nuv_luminosity_map_earth
        r = self.total_r_luminosity_map_earth

        # Return
        return make_colour_map(nuv, r)

    # -----------------------------------------------------------------

    @property
    def projected_ssfr_nuv_r_faceon_path(self):
        return fs.join(self.projected_path, "ssfr_nuv_r_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_projected_ssfr_nuv_r_faceon(self):
        return fs.is_file(self.projected_ssfr_nuv_r_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "projected_ssfr_nuv_r_faceon_path", True, write=False)
    def ssfr_nuv_r_faceon_map(self):

        """
        This function ...
        :return:
        """

        # Get maps
        nuv = self.total_nuv_luminosity_map_faceon
        r = self.total_r_luminosity_map_faceon

        # Return
        return make_colour_map(nuv, r)

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
    def bulge_intrinsic_fuv_luminosity(self):
        return self.model.intrinsic_fuv_luminosity_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_intrinsic_fuv_luminosity_scalar(self):
        return self.bulge_intrinsic_fuv_luminosity.to(self.specific_luminosity_unit, wavelength=self.fuv_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @property
    def disk_intrinsic_fuv_luminosity(self):
        return self.model.intrinsic_fuv_luminosity_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_intrinsic_fuv_luminosity_scalar(self):
        return self.disk_intrinsic_fuv_luminosity.to(self.specific_luminosity_unit, wavelength=self.fuv_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @property
    def young_intrinsic_fuv_luminosity(self):
        return self.model.intrinsic_fuv_luminosity_young

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_fuv_luminosity_scalar(self):
        return self.young_intrinsic_fuv_luminosity.to(self.specific_luminosity_unit, wavelength=self.fuv_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @property
    def sfr_intrinsic_fuv_luminosity(self):
        return self.model.intrinsic_fuv_luminosity_sfr

    # -----------------------------------------------------------------

    @property
    def sfr_intrinsic_fuv_luminosity_corrected(self):
        uncorrected = self.sfr_intrinsic_fuv_luminosity.to(self.specific_luminosity_unit, distance=self.galaxy_distance, wavelength=self.fuv_wavelength)
        absorbed_internal = self.internal_absorbed_fuv_luminosity.to(self.specific_luminosity_unit, distance=self.galaxy_distance, wavelength=self.fuv_wavelength)
        return uncorrected + absorbed_internal

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_fuv_luminosity_scalar(self):
        return self.sfr_intrinsic_fuv_luminosity.to(self.specific_luminosity_unit, wavelength=self.fuv_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_fuv_luminosity_scalar_corrected(self):
        return self.sfr_intrinsic_fuv_luminosity_corrected.to(self.specific_luminosity_unit, wavelength=self.fuv_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @property
    def bulge_intrinsic_i1_luminosity(self):
        return self.model.intrinsic_i1_luminosity_old_bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_intrinsic_i1_luminosity_scalar(self):
        return self.bulge_intrinsic_i1_luminosity.to(self.specific_luminosity_unit, wavelength=self.i1_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @property
    def disk_intrinsic_i1_luminosity(self):
        return self.model.intrinsic_i1_luminosity_old_disk

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_intrinsic_i1_luminosity_scalar(self):
        return self.disk_intrinsic_i1_luminosity.to(self.specific_luminosity_unit, wavelength=self.i1_wavelength, distance=self.galaxy_distance).value

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
    def sfr_cell_fuv_luminosities_corrected(self):
        return self.sfr_cell_normalized_mass * self.sfr_intrinsic_fuv_luminosity_scalar_corrected

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_cell_fuv_luminosities(self):
        return self.young_cell_fuv_luminosities + self.sfr_cell_fuv_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_cell_fuv_luminosities_corrected(self):
        return self.young_cell_fuv_luminosities + self.sfr_cell_fuv_luminosities_corrected

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_cell_fuv_luminosities(self):
        return self.bulge_cell_normalized_mass * self.bulge_intrinsic_fuv_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_cell_fuv_luminosities(self):
        return self.disk_cell_normalized_mass * self.disk_intrinsic_fuv_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def old_cell_fuv_luminosities(self):
        return self.bulge_cell_fuv_luminosities + self.disk_cell_fuv_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def total_cell_fuv_luminosities(self):
        return self.old_cell_fuv_luminosities + self.unevolved_cell_fuv_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_frequency_luminosity_conversion_factor(self):
        return self.specific_luminosity_unit.conversion_factor(self.frequency_luminosity_unit, wavelength=self.fuv_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_cell_fuv_frequency_luminosities(self):
        return self.total_cell_fuv_luminosities * self.fuv_frequency_luminosity_conversion_factor

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

        # Create the data with external xyz
        return Data3D.from_values(self.fuv_name, self.unevolved_cell_fuv_luminosities, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname, length_unit=self.length_unit, description=self.fuv_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.specific_luminosity_unit,
                                  distance=self.galaxy_distance, wavelength=self.fuv_wavelength)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_fuv_corrected_path", True, write=False)
    def fuv_data_corrected(self):

        """
        This function ...
        :return:
        """

        # Create the data with external xyz
        return Data3D.from_values(self.fuv_name, self.unevolved_cell_fuv_luminosities_corrected, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname, length_unit=self.length_unit, description=self.fuv_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.specific_luminosity_unit,
                                  distance=self.galaxy_distance, wavelength=self.fuv_wavelength)

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

        # Create the data with external xyz
        return Data3D.from_values(self.i1_name, self.old_cell_i1_luminosities,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.i1_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.specific_luminosity_unit,
                                  distance=self.galaxy_distance, wavelength=self.i1_wavelength)

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
    def cell_sfr_salim_corrected_path(self):
        return fs.join(self.cell_path, "sfr_salim_corrected.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_salim_corrected(self):
        return fs.is_file(self.cell_sfr_salim_corrected_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_sfr_salim_corrected_path", True, write=False)
    def sfr_salim_data_corrected(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell star formation rate (Salim, corrected) ...")

        # Calculate
        return salim_fuv_to_sfr(self.fuv_data_corrected)

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

    @property
    def cell_sfr_ke_corrected_path(self):
        return fs.join(self.cell_path, "sfr_ke_corrected.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_ke_corrected(self):
        return fs.is_file(self.cell_sfr_ke_corrected_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_sfr_ke_corrected_path", True, write=False)
    def sfr_ke_data_corrected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell star formation rate (Kennicutt & Evans, corrected) ...")

        # Calculate
        return kennicutt_evans_fuv_to_sfr(self.fuv_data_corrected)

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

    @property
    def sfr_unit(self):
        return "Msun/yr"

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_sfr_mappings_path", True, write=False)
    def sfr_mappings_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell star formation rate (MAPPINGS) ...")

        # Create the data with external xyz
        return Data3D.from_values(self.sfr_name, self.sfr_mappings_values,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.sfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.sfr_unit,
                                  distance=self.galaxy_distance)

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
        sfr_values = self.sfr_mappings_values + kennicutt_evans_fuv_to_sfr(self.young_cell_fuv_luminosities, unit=self.specific_luminosity_unit)

        # Create the data with external xyz
        return Data3D.from_values(self.sfr_name, sfr_values,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.sfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.sfr_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def dust_luminosity_unit(self):
        return self.total_contribution_absorption_unit

    # -----------------------------------------------------------------

    @property
    def total_cell_dust_luminosities(self):
        return self.total_contribution_absorption_luminosities # total absorbed energy = total emitted energy
        # integrating the dust emission spectral data is too expensive

    # -----------------------------------------------------------------

    @property
    def cell_sfr_tir_path(self):
        return fs.join(self.cell_path, "sfr_tir.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_tir(self):
        return fs.is_file(self.cell_sfr_tir_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_sfr_tir_path", True, write=False)
    def sfr_tir_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell star formation rate (TIR) ...")

        # Calculate
        sfr_values = kennicutt_tir_to_sfr(self.total_cell_dust_luminosities, unit=self.dust_luminosity_unit)

        # Create the data with external xyz
        return Data3D.from_values(self.sfr_name, sfr_values, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.sfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit="Msun/yr",
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # ABSORPTION
    # -----------------------------------------------------------------

    @property
    def absorption_path(self):
        return self.analysis_run.absorption_path

    # -----------------------------------------------------------------

    @property
    def absorption_sfr_path(self):
        return fs.join(self.absorption_path, "sfr")

    # -----------------------------------------------------------------
    #  Internal
    # -----------------------------------------------------------------

    @property
    def sfr_absorption_properties_path(self):
        return fs.join(self.absorption_path, "sfr.txt")

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_properties(self):
        return open_box(self.sfr_absorption_properties_path)

    # -----------------------------------------------------------------

    @property
    def internal_absorbed_luminosity(self):
        return self.sfr_absorption_properties.internal.absorbed

    # -----------------------------------------------------------------

    @property
    def internal_absorbed_fuv_luminosity(self):
        return self.sfr_absorption_properties.internal.absorbed_fuv

    # -----------------------------------------------------------------

    @property
    def internal_absorption_sed_path(self):
        return fs.join(self.absorption_sfr_path, "absorption_internal.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def internal_absorption_sed(self):
        return SED.from_file(self.internal_absorption_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def internal_absorbed_fuv_luminosity_map_earth(self):
        return self.model.sfr_map_earth.normalized(to=self.sfr_intrinsic_fuv_luminosity_corrected)

    # -----------------------------------------------------------------

    @lazyproperty
    def internal_absorbed_fuv_luminosity_map_faceon(self):
        return self.model.sfr_map_faceon.normalized(to=self.sfr_intrinsic_fuv_luminosity_corrected)

    # -----------------------------------------------------------------

    @lazyproperty
    def internal_absorbed_fuv_luminosity_map_edgeon(self):
        return self.model.sfr_map_edgeon.normalized(to=self.sfr_intrinsic_fuv_luminosity_corrected)

    # -----------------------------------------------------------------
    #  Total
    # -----------------------------------------------------------------

    @property
    def total_emission_spectral_data_path(self):
        return fs.join(self.absorption_path, "total_emission.dat")

    # -----------------------------------------------------------------

    @property
    def total_absorption_spectral_data_path(self):
        return fs.join(self.absorption_path, "total_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_emission_spectral_data(self):
        return fs.is_file(self.total_emission_spectral_data_path)

    # -----------------------------------------------------------------

    @property
    def has_total_absorption_spectral_data(self):
        return fs.is_file(self.total_absorption_spectral_data_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_emission_spectral_data(self):
        return SpectralData3D.from_file(self.total_emission_spectral_data_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_spectral_data(self):
        return SpectralData3D.from_file(self.total_absorption_spectral_data_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_cell_24um_luminosities(self):
        return self.total_emission_spectral_data.get_array_for_wavelength(self.model.mips24_wavelength)

    # -----------------------------------------------------------------

    @property
    def mips24_luminosity_unit(self):
        return self.total_emission_spectral_data.unit

    # -----------------------------------------------------------------

    @property
    def cell_sfr_24um_path(self):
        return fs.join(self.cell_path, "sfr_24um.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_sfr_24um(self):
        return fs.is_file(self.cell_sfr_24um_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_sfr_24um_path", True, write=False)
    def sfr_24um_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the cell star formation rate (24 micron) ...")

        # Calculate
        sfr_values = calzetti_24um_to_sfr(self.total_cell_24um_luminosities, unit=self.mips24_luminosity_unit)

        # Create the data with external xyz
        return Data3D.from_values(self.sfr_name, sfr_values, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.sfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit="Msun/yr",
                                  distance=self.galaxy_distance)

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
    def cell_sfrs_salim_corrected(self):
        return self.sfr_salim_data_corrected.values

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
    def cell_ssfrs_salim_corrected(self):
        return self.cell_sfrs_salim_corrected / self.cell_masses

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

        # Create the data with external xyz
        return Data3D.from_values(self.ssfr_name, self.cell_ssfrs_salim,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.ssfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.ssfr_salim_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_salim_corrected_path(self):
        return fs.join(self.cell_path, "ssfr_salim_corrected.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_salim_corrected(self):
        return fs.is_file(self.cell_ssfr_salim_corrected_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_ssfr_salim_corrected_path", True, write=False)
    def ssfr_salim_data_corrected(self):

        """
        This function ...
        :return:
        """

        # Create the data with external xyz
        return Data3D.from_values(self.ssfr_name, self.cell_ssfrs_salim_corrected,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.ssfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.ssfr_salim_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def cell_sfrs_ke(self):
        return self.sfr_ke_data.values

    # -----------------------------------------------------------------

    @property
    def cell_sfrs_ke_corrected(self):
        return self.sfr_ke_data_corrected.values

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
    def cell_ssfrs_ke_corrected(self):
        return self.cell_sfrs_ke_corrected / self.cell_masses

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

        # Create the data with external xyz
        return Data3D.from_values(self.ssfr_name, self.cell_ssfrs_ke,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.ssfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.ssfr_ke_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_ke_corrected_path(self):
        return fs.join(self.cell_path, "ssfr_ke_corrected.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_ke_corrected(self):
        return fs.is_file(self.cell_ssfr_ke_corrected_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_ssfr_ke_corrected_path", True, write=False)
    def ssfr_ke_data_corrected(self):

        """
        This function ...
        :return:
        """

        # Create the data with external xyz
        return Data3D.from_values(self.ssfr_name, self.cell_ssfrs_ke_corrected,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.ssfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.ssfr_ke_unit,
                                  distance=self.galaxy_distance)

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

        # Create the data with external xyz
        return Data3D.from_values(self.ssfr_name, self.cell_ssfrs_mappings,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.ssfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.ssfr_mappings_unit,
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

        # Create the data with external xyz
        return Data3D.from_values(self.ssfr_name, self.cell_ssfrs_mappings_ke,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.ssfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.ssfr_mappings_ke_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # H luminosities
    # -----------------------------------------------------------------

    @lazyproperty
    def h_wavelength(self):
        return parse_filter("2MASS H").wavelength

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_intrinsic_h_luminosity(self):
        return self.model.bulge_simulations.intrinsic_photometry_at(self.h_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_intrinsic_h_luminosity_scalar(self):
        return self.bulge_intrinsic_h_luminosity.to(self.specific_luminosity_unit, wavelength=self.h_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_intrinsic_h_luminosity(self):
        return self.model.disk_simulations.intrinsic_photometry_at(self.h_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_intrinsic_h_luminosity_scalar(self):
        return self.disk_intrinsic_h_luminosity.to(self.specific_luminosity_unit, wavelength=self.h_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_h_luminosity(self):
        return self.model.young_simulations.intrinsic_photometry_at(self.h_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_h_luminosity_scalar(self):
        return self.young_intrinsic_h_luminosity.to(self.specific_luminosity_unit, wavelength=self.h_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_h_luminosity(self):
        return self.model.sfr_simulations.intrinsic_photometry_at(self.h_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_h_luminosity_scalar(self):
        return self.sfr_intrinsic_h_luminosity.to(self.specific_luminosity_unit, wavelength=self.h_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_cell_h_luminosities(self):
        return self.bulge_cell_normalized_mass * self.bulge_intrinsic_h_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_cell_h_luminosities(self):
        return self.disk_cell_normalized_mass * self.disk_intrinsic_h_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def old_cell_h_luminosities(self):
        return self.bulge_cell_h_luminosities + self.disk_cell_h_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def young_cell_h_luminosities(self):
        return self.young_cell_normalized_mass * self.young_intrinsic_h_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_cell_h_luminosities(self):
        return self.sfr_cell_normalized_mass * self.sfr_intrinsic_h_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_cell_h_luminosities(self):
        return self.young_cell_h_luminosities + self.sfr_cell_h_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def total_cell_h_luminosities(self):
        return self.old_cell_h_luminosities + self.unevolved_cell_h_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def h_frequency_luminosity_conversion_factor(self):
        return self.specific_luminosity_unit.conversion_factor(self.frequency_luminosity_unit, wavelength=self.h_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_cell_h_frequency_luminosities(self):
        return self.total_cell_h_luminosities * self.h_frequency_luminosity_conversion_factor

    # -----------------------------------------------------------------
    # R luminosities
    # -----------------------------------------------------------------

    @lazyproperty
    def r_wavelength(self):
        return parse_filter("SDSS r").wavelength

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_intrinsic_r_luminosity(self):
        return self.model.bulge_simulations.intrinsic_photometry_at(self.r_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_intrinsic_r_luminosity_scalar(self):
        return self.bulge_intrinsic_r_luminosity.to(self.specific_luminosity_unit, wavelength=self.r_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_intrinsic_r_luminosity(self):
        return self.model.disk_simulations.intrinsic_photometry_at(self.r_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_intrinsic_r_luminosity_scalar(self):
        return self.disk_intrinsic_r_luminosity.to(self.specific_luminosity_unit, wavelength=self.r_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_r_luminosity(self):
        return self.model.young_simulations.intrinsic_photometry_at(self.r_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_r_luminosity_scalar(self):
        return self.young_intrinsic_r_luminosity.to(self.specific_luminosity_unit, wavelength=self.r_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_r_luminosity(self):
        return self.model.sfr_simulations.intrinsic_photometry_at(self.r_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_r_luminosity_scalar(self):
        return self.sfr_intrinsic_r_luminosity.to(self.specific_luminosity_unit, wavelength=self.r_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_cell_r_luminosities(self):
        return self.bulge_cell_normalized_mass * self.bulge_intrinsic_r_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_cell_r_luminosities(self):
        return self.disk_cell_normalized_mass * self.disk_intrinsic_r_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def old_cell_r_luminosities(self):
        return self.bulge_cell_r_luminosities + self.disk_cell_r_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def young_cell_r_luminosities(self):
        return self.young_cell_normalized_mass * self.young_intrinsic_r_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_cell_r_luminosities(self):
        return self.sfr_cell_normalized_mass * self.sfr_intrinsic_r_luminosity_scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_cell_r_luminosities(self):
        return self.young_cell_r_luminosities + self.sfr_cell_r_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def total_cell_r_luminosities(self):
        return self.old_cell_r_luminosities + self.unevolved_cell_r_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def r_frequency_luminosity_conversion_factor(self):
        return self.specific_luminosity_unit.conversion_factor(self.frequency_luminosity_unit, wavelength=self.r_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_cell_r_frequency_luminosities(self):
        return self.total_cell_r_luminosities * self.r_frequency_luminosity_conversion_factor

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def cell_ssfr_fuv_h_path(self):
        return fs.join(self.cell_path, "ssfr_fuv_h.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_fuv_h(self):
        return fs.is_file(self.cell_ssfr_fuv_h_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_ssfr_fuv_h_path", True, write=False)
    def ssfr_fuv_h_data(self):

        """
        Thisf unction ...
        :return:
        """

        # Calculate the colours
        fuv_h = -2.5 * np.log10(self.total_cell_fuv_frequency_luminosities / self.total_cell_h_frequency_luminosities)

        # Create the data
        return Data3D.from_values(self.ssfr_name, fuv_h, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.ssfr_description,
                                  xyz_filepath = self.cell_coordinates_filepath, distance = self.galaxy_distance) # no unit

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_fuv_r_path(self):
        return fs.join(self.cell_path, "ssfr_fuv_r.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_fuv_r(self):
        return fs.is_file(self.cell_ssfr_fuv_r_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_ssfr_fuv_r_path", True, write=False)
    def ssfr_fuv_r_data(self):

        """
        This function ...
        :return:
        """

        # Calculate the colours
        fuv_r = -2.5 * np.log10(self.total_cell_fuv_frequency_luminosities / self.total_cell_r_frequency_luminosities)

        # Create the data
        return Data3D.from_values(self.ssfr_name, fuv_r, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.ssfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, distance=self.galaxy_distance)  # no unit

    # -----------------------------------------------------------------
    # NUV luminosities
    #   Bulge
    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_intrinsic_nuv_luminosity(self):
        return self.model.bulge_simulations.intrinsic_photometry_at(self.nuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_intrinsic_nuv_luminosity_scalar(self):
        return self.bulge_intrinsic_nuv_luminosity.to(self.specific_luminosity_unit, wavelength=self.nuv_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_cell_nuv_luminosities(self):
        return self.bulge_cell_normalized_mass * self.bulge_intrinsic_nuv_luminosity_scalar

    # -----------------------------------------------------------------
    #   Disk
    # -----------------------------------------------------------------

    @lazyproperty
    def disk_intrinsic_nuv_luminosity(self):
        return self.model.disk_simulations.intrinsic_photometry_at(self.nuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_intrinsic_nuv_luminosity_scalar(self):
        return self.disk_intrinsic_nuv_luminosity.to(self.specific_luminosity_unit, wavelength=self.nuv_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_cell_nuv_luminosities(self):
        return self.disk_cell_normalized_mass * self.disk_intrinsic_nuv_luminosity_scalar

    # -----------------------------------------------------------------
    #   Old
    # -----------------------------------------------------------------

    @lazyproperty
    def old_cell_nuv_luminosities(self):
        return self.bulge_cell_nuv_luminosities + self.disk_cell_nuv_luminosities

    # -----------------------------------------------------------------
    #   Young
    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_nuv_luminosity(self):
        return self.model.young_simulations.intrinsic_photometry_at(self.nuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_intrinsic_nuv_luminosity_scalar(self):
        return self.young_intrinsic_nuv_luminosity.to(self.specific_luminosity_unit, wavelength=self.nuv_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def young_cell_nuv_luminosities(self):
        return self.young_cell_normalized_mass * self.young_intrinsic_nuv_luminosity_scalar

    # -----------------------------------------------------------------
    #   SFR
    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_nuv_luminosity(self):
        return self.model.sfr_simulations.intrinsic_photometry_at(self.nuv_wavelength, interpolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_intrinsic_nuv_luminosity_scalar(self):
        return self.sfr_intrinsic_nuv_luminosity.to(self.specific_luminosity_unit, wavelength=self.nuv_wavelength, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_cell_nuv_luminosities(self):
        return self.sfr_cell_normalized_mass * self.sfr_intrinsic_nuv_luminosity_scalar

    # -----------------------------------------------------------------
    #   Unevolved
    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_cell_nuv_luminosities(self):
        return self.sfr_cell_normalized_mass * self.sfr_intrinsic_nuv_luminosity_scalar

    # -----------------------------------------------------------------
    #   Total
    # -----------------------------------------------------------------

    @lazyproperty
    def total_cell_nuv_luminosities(self):
        return self.old_cell_nuv_luminosities + self.unevolved_cell_nuv_luminosities

    # -----------------------------------------------------------------

    @lazyproperty
    def nuv_frequency_luminosity_conversion_factor(self):
        return self.specific_luminosity_unit.conversion_factor(self.frequency_luminosity_unit, wavelength=self.nuv_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_cell_nuv_frequency_luminosities(self):
        return self.total_cell_nuv_luminosities * self.nuv_frequency_luminosity_conversion_factor

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def cell_ssfr_nuv_h_path(self):
        return fs.join(self.cell_path, "ssfr_nuv_h.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_nuv_h(self):
        return fs.is_file(self.cell_ssfr_nuv_h_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_ssfr_nuv_h_path", True, write=False)
    def ssfr_nuv_h_data(self):

        """
        This function ...
        :return:
        """

        # Calculate the colours
        nuv_h = -2.5 * np.log10(self.total_cell_nuv_frequency_luminosities / self.total_cell_h_frequency_luminosities)

        # Create the data
        return Data3D.from_values(self.ssfr_name, nuv_h, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.ssfr_description,
                                  xyz_filepath = self.cell_coordinates_filepath, distance = self.galaxy_distance) # no unit

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_nuv_r_path(self):
        return fs.join(self.cell_path, "ssfr_nuv_r.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_nuv_r(self):
        return fs.is_file(self.cell_ssfr_nuv_r_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "cell_ssfr_nuv_r_path", True, write=False)
    def ssfr_nuv_r_data(self):

        """
        This function ...
        :return:
        """

        # Calculate the colours
        nuv_r = -2.5 * np.log10(self.total_cell_nuv_frequency_luminosities / self.total_cell_r_frequency_luminosities)

        # Create the data
        return Data3D.from_values(self.ssfr_name, nuv_r, self.cell_x_coordinates_colname,
                                  self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.ssfr_description,
                                  xyz_filepath=self.cell_coordinates_filepath, distance=self.galaxy_distance)  # no unit

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

        # TIR
        self.write_projected_sfr_tir()

        # 24 micron
        self.write_projected_sfr_24um()

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_salim_earth(self):
        return not self.has_projected_sfr_salim_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_salim_faceon(self):
        return not self.has_projected_sfr_salim_faceon

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_salim_earth_corrected(self):
       return not self.has_projected_sfr_salim_earth_corrected

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_salim_faceon_corrected(self):
        return not self.has_projected_sfr_salim_faceon_corrected

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

        # Corrected earth
        if self.do_write_projected_sfr_salim_earth_corrected: self.write_projected_sfr_salim_earth_corrected()

        # Corrected faceon
        if self.do_write_projected_sfr_salim_faceon_corrected: self.write_projected_sfr_salim_faceon_corrected()

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

    def write_projected_sfr_salim_earth_corrected(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_salim_earth_map_corrected.saveto(self.projected_sfr_salim_earth_corrected_path)

    # -----------------------------------------------------------------

    def write_projected_sfr_salim_faceon_corrected(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_salim_faceon_map_corrected.saveto(self.projected_sfr_salim_faceon_corrected_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_ke_earth(self):
        return not self.has_projected_sfr_ke_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_ke_faceon(self):
        return not self.has_projected_sfr_ke_faceon

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_ke_earth_corrected(self):
        return not self.has_projected_sfr_ke_earth_corrected

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_ke_faceon_corrected(self):
        return not self.has_projected_sfr_ke_faceon_corrected

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

        # Corrected earth
        if self.do_write_projected_sfr_ke_earth_corrected: self.write_projected_sfr_ke_earth_corrected()

        # Corrected faceon
        if self.do_write_projected_sfr_ke_faceon_corrected: self.write_projected_sfr_ke_faceon_corrected()

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

    def write_projected_sfr_ke_earth_corrected(self):

        """
        This function ...
        :return:
        """

        self.sfr_ke_earth_map_corrected.saveto(self.projected_sfr_ke_earth_corrected_path)

    # -----------------------------------------------------------------

    def write_projected_sfr_ke_faceon_corrected(self):

        """
        This function ...
        :return:
        """

        self.sfr_ke_faceon_map_corrected.saveto(self.projected_sfr_ke_faceon_corrected_path)

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
    def do_write_projected_sfr_tir_earth(self):
        return not self.has_projected_sfr_tir_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_tir_faceon(self):
        return not self.has_projected_sfr_tir_faceon

    # -----------------------------------------------------------------

    def write_projected_sfr_tir(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_projected_sfr_tir_earth: self.write_projected_sfr_tir_earth()

        # Faceon
        if self.do_write_projected_sfr_tir_faceon: self.write_projected_sfr_tir_faceon()

    # -----------------------------------------------------------------

    def write_projected_sfr_tir_earth(self):

        """
        This function ...
        :return:
        """

        self.sfr_tir_earth_map.saveto(self.projected_sfr_tir_earth_path)

    # -----------------------------------------------------------------

    def write_projected_sfr_tir_faceon(self):

        """
        This function ...
        :return:
        """

        self.sfr_tir_faceon_map.saveto(self.projected_sfr_tir_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_24um_earth(self):
        return not self.has_projected_sfr_24um_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_sfr_24um_faceon(self):
        return not self.has_projected_sfr_24um_faceon

    # -----------------------------------------------------------------

    def write_projected_sfr_24um(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_projected_sfr_24um_earth: self.write_projected_sfr_24um_earth()

        # Faceon
        if self.do_write_projected_sfr_24um_faceon: self.write_projected_sfr_24um_faceon()

    # -----------------------------------------------------------------

    def write_projected_sfr_24um_earth(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_24um_earth_map.saveto(self.projected_sfr_24um_earth_path)

    # -----------------------------------------------------------------

    def write_projected_sfr_24um_faceon(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_24um_faceon_map.saveto(self.projected_sfr_24um_faceon_path)

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

        # FUV-H
        self.write_projected_ssfr_fuv_h()

        # FUV-r
        self.write_projected_ssfr_fuv_r()

        # NUV-H
        self.write_projected_ssfr_nuv_h()

        # NUV-r
        self.write_projected_ssfr_nuv_r()

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_salim_earth(self):
        return not self.has_projected_ssfr_salim_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_salim_faceon(self):
        return not self.has_projected_ssfr_salim_faceon

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_salim_earth_corrected(self):
        return not self.has_projected_ssfr_salim_earth_corrected

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_salim_faceon_corrected(self):
        return not self.has_projected_ssfr_salim_faceon_corrected

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

        # Corrected earth
        if self.do_write_projected_ssfr_salim_earth_corrected: self.write_projected_ssfr_salim_earth_corrected()

        # Corrected faceon
        if self.do_write_projected_ssfr_salim_faceon_corrected: self.write_projected_ssfr_salim_faceon_corrected()

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

    def write_projected_ssfr_salim_earth_corrected(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_salim_earth_map_corrected.saveto(self.projected_ssfr_salim_earth_corrected_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_salim_faceon_corrected(self):

        """
        This fucntion ...
        :return:
        """

        # Write
        self.ssfr_salim_faceon_map_corrected.saveto(self.projected_ssfr_salim_faceon_corrected_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_ke_earth(self):
        return not self.has_projected_ssfr_ke_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_ke_faceon(self):
        return not self.has_projected_ssfr_ke_faceon

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_ke_earth_corrected(self):
        return not self.has_projected_ssfr_ke_earth_corrected

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_ke_faceon_corrected(self):
        return not self.has_projected_ssfr_ke_faceon_corrected

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

        # Corrected earth
        if self.do_write_projected_ssfr_ke_earth_corrected: self.write_projected_ssfr_ke_earth_corrected()

        # Corrected faceon
        if self.do_write_projected_ssfr_ke_faceon_corrected: self.write_projected_ssfr_ke_faceon_corrected()

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

    def write_projected_ssfr_ke_earth_corrected(self):

        """
        This function ...
        :return:
        """

        self.ssfr_ke_earth_map_corrected.saveto(self.projected_ssfr_ke_earth_corrected_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_ke_faceon_corrected(self):

        """
        This function ...
        :return:
        """

        self.ssfr_ke_faceon_map_corrected.saveto(self.projected_ssfr_ke_faceon_corrected_path)

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

    @property
    def do_write_projected_ssfr_fuv_h_earth(self):
        return not self.has_projected_ssfr_fuv_h_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_fuv_h_faceon(self):
        return not self.has_projected_ssfr_fuv_h_faceon

    # -----------------------------------------------------------------

    def write_projected_ssfr_fuv_h(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_write_projected_ssfr_fuv_h_earth: self.write_projected_ssfr_fuv_h_earth()

        # Faceon
        if self.do_write_projected_ssfr_fuv_h_faceon: self.write_projected_ssfr_fuv_h_faceon()

    # -----------------------------------------------------------------

    def write_projected_ssfr_fuv_h_earth(self):

        """
        This function ...
        :return:
        """

        self.ssfr_fuv_h_earth_map.saveto(self.projected_ssfr_fuv_h_earth_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_fuv_h_faceon(self):

        """
        This function ...
        :return:
        """

        self.ssfr_fuv_h_faceon_map.saveto(self.projected_ssfr_fuv_h_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_fuv_r_earth(self):
        return not self.has_projected_ssfr_fuv_r_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_fuv_r_faceon(self):
        return not self.has_projected_ssfr_fuv_r_faceon

    # -----------------------------------------------------------------

    def write_projected_ssfr_fuv_r(self):

        """
        Thisf unction ...
        :return:
        """

        # Earth
        if self.do_write_projected_ssfr_fuv_r_earth: self.write_projected_ssfr_fuv_r_earth()

        # Faceon
        if self.do_write_projected_ssfr_fuv_r_faceon: self.write_projected_ssfr_fuv_r_faceon()

    # -----------------------------------------------------------------

    def write_projected_ssfr_fuv_r_earth(self):

        """
        This function ...
        :return:
        """

        self.ssfr_fuv_r_earth_map.saveto(self.projected_ssfr_fuv_r_earth_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_fuv_r_faceon(self):

        """
        This function ...
        :return:
        """

        self.ssfr_fuv_r_faceon_map.saveto(self.projected_ssfr_fuv_r_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_nuv_h_earth(self):
        return not self.has_projected_ssfr_nuv_h_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_nuv_h_faceon(self):
        return not self.has_projected_ssfr_nuv_h_faceon

    # -----------------------------------------------------------------

    def write_projected_ssfr_nuv_h(self):

        """
        This function ...
        :return:
        """

        # Earth
        self.write_projected_ssfr_nuv_h_earth()

        # Faceon
        self.write_projected_ssfr_nuv_h_faceon()

    # -----------------------------------------------------------------

    def write_projected_ssfr_nuv_h_earth(self):

        """
        This function ...
        :return:
        """

        self.ssfr_nuv_h_earth_map.saveto(self.projected_ssfr_nuv_h_earth_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_nuv_h_faceon(self):

        """
        This function ...
        :return:
        """

        self.ssfr_nuv_h_faceon_map.saveto(self.projected_ssfr_nuv_h_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_nuv_r_earth(self):
        return not self.has_projected_ssfr_nuv_r_earth

    # -----------------------------------------------------------------

    @property
    def do_write_projected_ssfr_nuv_r_faceon(self):
        return not self.has_projected_ssfr_nuv_r_faceon

    # -----------------------------------------------------------------

    def write_projected_ssfr_nuv_r(self):

        """
        Thisn function ...
        :return:
        """

        # Earth
        self.write_projected_ssfr_nuv_r_earth()

        # Faceon
        self.write_projected_ssfr_nuv_r_faceon()

    # -----------------------------------------------------------------

    def write_projected_ssfr_nuv_r_earth(self):

        """
        This function ...
        :return:
        """

        self.ssfr_nuv_r_earth_map.saveto(self.projected_ssfr_nuv_r_earth_path)

    # -----------------------------------------------------------------

    def write_projected_ssfr_nuv_r_faceon(self):

        """
        Thins function ...
        :return:
        """

        self.ssfr_nuv_r_faceon_map.saveto(self.projected_ssfr_nuv_r_faceon_path)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def do_write_cell_fuv(self):
        return not self.has_cell_fuv

    # -----------------------------------------------------------------

    @property
    def do_write_cell_fuv_corrected(self):
        return not self.has_cell_fuv_corrected

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

        # Intrinsic FUV luminosity, corrected for internal absorption
        if self.do_write_cell_fuv_corrected: self.write_cell_fuv_corrected()

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
        log.info("Writing the cell intrinsic FUV luminosities ...")

        # Write
        self.fuv_data.saveto(self.cell_fuv_path)

    # -----------------------------------------------------------------

    @property
    def cell_fuv_corrected_path(self):
        return fs.join(self.cell_path, "fuv_corrected.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_fuv_corrected(self):
        return fs.is_file(self.cell_fuv_corrected_path)

    # -----------------------------------------------------------------

    def write_cell_fuv_corrected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cell intrinsic FUV luminosities (corrected) ...")

        # Write
        self.fuv_data_corrected.saveto(self.cell_fuv_corrected_path)

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
    def do_write_cell_sfr_salim_corrected(self):
        return not self.has_cell_sfr_salim_corrected

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_ke(self):
        return not self.has_cell_sfr_ke

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_ke_corrected(self):
        return not self.has_cell_sfr_ke_corrected

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_mappings(self):
        return not self.has_cell_sfr_mappings

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_mappings_ke(self):
        return not self.has_cell_sfr_mappings_ke

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_tir(self):
        return not self.has_cell_sfr_tir

    # -----------------------------------------------------------------

    @property
    def do_write_cell_sfr_24um(self):
        return not self.has_cell_sfr_24um

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

        # Salim, corrected
        if self.do_write_cell_sfr_salim_corrected: self.write_cell_sfr_salim_corrected()

        # K&E
        if self.do_write_cell_sfr_ke: self.write_cell_sfr_ke()

        # K&E, corrected
        if self.do_write_cell_sfr_ke_corrected: self.write_cell_sfr_ke_corrected()

        # MAPPINGS
        if self.do_write_cell_sfr_mappings: self.write_cell_sfr_mappings()

        # MAPPINGS + K&E
        if self.do_write_cell_sfr_mappings_ke: self.write_cell_sfr_mappings_ke()

        # TIR
        if self.do_write_cell_sfr_tir: self.write_cell_sfr_tir()

        # 24 micron
        if self.do_write_cell_sfr_24um: self.write_cell_sfr_24um()

    # -----------------------------------------------------------------

    def write_cell_sfr_salim(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_salim_data.saveto(self.cell_sfr_salim_path)

    # -----------------------------------------------------------------

    def write_cell_sfr_salim_corrected(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_salim_data_corrected.saveto(self.cell_sfr_salim_corrected_path)

    # -----------------------------------------------------------------

    def write_cell_sfr_ke(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_ke_data.saveto(self.cell_sfr_ke_path)

    # -----------------------------------------------------------------

    def write_cell_sfr_ke_corrected(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_ke_data_corrected.saveto(self.cell_sfr_ke_corrected_path)

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

    def write_cell_sfr_tir(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_tir_data.saveto(self.cell_sfr_tir_path)

    # -----------------------------------------------------------------

    def write_cell_sfr_24um(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_24um_data.saveto(self.cell_sfr_24um_path)

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
    def do_write_cell_ssfr_salim_corrected(self):
        return not self.has_cell_ssfr_salim_corrected

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_ke(self):
        return not self.has_cell_ssfr_ke

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_ke_corrected(self):
        return not self.has_cell_ssfr_ke_corrected

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_mappings(self):
        return not self.has_cell_ssfr_mappings

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_mappings_ke(self):
        return not self.has_cell_ssfr_mappings_ke

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_fuv_h(self):
        return not self.has_cell_ssfr_fuv_h

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_fuv_r(self):
        return not self.has_cell_ssfr_fuv_r

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_nuv_h(self):
        return not self.has_cell_ssfr_nuv_h

    # -----------------------------------------------------------------

    @property
    def do_write_cell_ssfr_nuv_r(self):
        return not self.has_cell_ssfr_nuv_r

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

        # Salim, corrected
        if self.do_write_cell_ssfr_salim_corrected: self.write_cell_ssfr_salim_corrected()

        # K&E
        if self.do_write_cell_ssfr_ke: self.write_cell_ssfr_ke()

        # K&E, corrected
        if self.do_write_cell_ssfr_ke_corrected: self.write_cell_ssfr_ke_corrected()

        # MAPPINGS
        if self.do_write_cell_ssfr_mappings: self.write_cell_ssfr_mappings()

        # MAPPINGS + K&E
        if self.do_write_cell_ssfr_mappings_ke: self.write_cell_ssfr_mappings_ke()

        # FUV-H
        if self.do_write_cell_ssfr_fuv_h: self.write_cell_ssfr_fuv_h()

        # FUV-R
        if self.do_write_cell_ssfr_fuv_r: self.write_cell_ssfr_fuv_r()

        # NUV-H
        if self.do_write_cell_ssfr_nuv_h: self.write_cell_ssfr_nuv_h()

        # NUV-R
        if self.do_write_cell_ssfr_nuv_r: self.write_cell_ssfr_nuv_r()

    # -----------------------------------------------------------------

    def write_cell_ssfr_salim(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_salim_data.saveto(self.cell_ssfr_salim_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_salim_corrected(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_salim_data_corrected.saveto(self.cell_ssfr_salim_corrected_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_ke(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_ke_data.saveto(self.cell_ssfr_ke_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_ke_corrected(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_ke_data_corrected.saveto(self.cell_ssfr_ke_corrected_path)

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

    def write_cell_ssfr_fuv_h(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_fuv_h_data.saveto(self.cell_ssfr_fuv_h_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_fuv_r(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_fuv_r_data.saveto(self.cell_ssfr_fuv_r_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_nuv_h(self):

        """
        This function ...
        :return:
        """

        # Write
        self.ssfr_nuv_h_data.saveto(self.cell_ssfr_nuv_h_path)

    # -----------------------------------------------------------------

    def write_cell_ssfr_nuv_r(self):

        """
        Thins function ...
        :return:
        """

        # Write
        self.ssfr_nuv_r_data.saveto(self.cell_ssfr_nuv_r_path)

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

        # TIR
        self.plot_projected_sfr_tir()

        # 24um
        self.plot_projected_sfr_24um()

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
    def do_plot_projected_sfr_tir_earth(self):
        return not self.has_projected_sfr_tir_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_tir_faceon(self):
        return not self.has_projected_sfr_tir_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_sfr_tir(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_sfr_tir_earth: self.plot_projected_sfr_tir_earth()

        # Faceon
        if self.do_plot_projected_sfr_tir_faceon: self.plot_projected_sfr_tir_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_sfr_tir_earth_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_tir_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_tir_earth_map_plot(self):
        return fs.is_file(self.projected_sfr_tir_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_tir_earth(self):

        """
        This function ...
        :return:
        """

        self.plot_sfr_map(self.sfr_tir_earth_map, self.projected_sfr_tir_earth_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def projected_sfr_tir_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_tir_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_tir_faceon_map_plot(self):
        return fs.is_file(self.projected_sfr_tir_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_tir_faceon(self):

        """
        This function ...
        :return:
        """

        self.plot_sfr_map(self.sfr_tir_faceon_map, self.projected_sfr_tir_faceon_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_24um_earth(self):
        return not self.has_projected_sfr_24um_earth_map_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_projected_sfr_24um_faceon(self):
        return not self.has_projected_sfr_24um_faceon_map_plot

    # -----------------------------------------------------------------

    def plot_projected_sfr_24um(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_plot_projected_sfr_24um_earth: self.plot_projected_sfr_24um_earth()

        # Faceon
        if self.do_plot_projected_sfr_24um_faceon: self.plot_projected_sfr_24um_faceon()

    # -----------------------------------------------------------------

    @property
    def projected_sfr_24um_earth_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_24um_earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_24um_earth_map_plot(self):
        return fs.is_file(self.projected_sfr_24um_earth_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_24um_earth(self):

        """
        This function ...
        :return:
        """

        self.plot_sfr_map(self.sfr_24um_earth_map, self.projected_sfr_24um_earth_map_plot_path)

    # -----------------------------------------------------------------

    @property
    def projected_sfr_24um_faceon_map_plot_path(self):
        return fs.join(self.projected_path, "sfr_24um_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_projected_sfr_24um_faceon_map_plot(self):
        return fs.is_file(self.projected_sfr_24um_faceon_map_plot_path)

    # -----------------------------------------------------------------

    def plot_projected_sfr_24um_faceon(self):

        """
        This function ...
        :return:
        """

        self.plot_sfr_map(self.sfr_24um_faceon_map, self.projected_sfr_24um_faceon_map_plot_path)

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
