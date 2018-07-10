#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.properties Contains the PropertiesAnalyser class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import AnalysisRunComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty
from ...core.tools.serialization import write_dict
from ...core.basics.containers import DefaultOrderedDict
from ...magic.core.frame import Frame
from ...magic.tools.plotting import plot_map

from .component import bol_map_name, intr_stellar_map_name, obs_stellar_map_name, diffuse_dust_map_name, dust_map_name
from .component import scattered_map_name, absorbed_diffuse_map_name, fabs_diffuse_map_name, fabs_map_name, stellar_mass_map_name, ssfr_map_name
from .component import attenuated_map_name, direct_map_name, sfr_map_name, i1_map_name, intr_i1_map_name, fuv_map_name
from .component import intr_fuv_map_name, dust_mass_map_name, stellar_lum_map_name, intr_dust_map_name
from .component import diffuse_mass_map_name, mass_map_name, earth_name, faceon_name, edgeon_name

# -----------------------------------------------------------------

class PropertiesAnalyser(AnalysisRunComponent):
    
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
        super(PropertiesAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The maps
        self.total_maps = DefaultOrderedDict(OrderedDict)
        self.bulge_maps = DefaultOrderedDict(OrderedDict)
        self.disk_maps = DefaultOrderedDict(OrderedDict)
        self.old_maps = DefaultOrderedDict(OrderedDict)
        self.young_maps = DefaultOrderedDict(OrderedDict)
        self.sfr_maps = DefaultOrderedDict(OrderedDict)
        self.unevolved_maps = DefaultOrderedDict(OrderedDict)
        self.dust_maps = DefaultOrderedDict(OrderedDict)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get the maps
        self.get_maps()

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

        # Call the setup of the base class
        super(PropertiesAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @property
    def model(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.model

    # -----------------------------------------------------------------

    def get_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the maps ...")

        # Total
        self.get_total_maps()

        # Bulge
        self.get_bulge_maps()

        # Disk
        self.get_disk_maps()

        # Old
        self.get_old_maps()

        # Young
        self.get_young_maps()

        # SFR
        self.get_sfr_maps()

        # Unevolved
        self.get_unevolved_maps()

        # Dust
        self.get_dust_maps()

    # -----------------------------------------------------------------

    def get_total_maps(self):

        """
        Thisf unction ...
        :return:
        """

        # Earth
        if self.do_earth: self.get_total_maps_earth()

        # Face-on
        if self.do_faceon: self.get_total_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.get_total_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def total_earth_maps(self):

        """
        This function ...
        :return:
        """

        return self.total_maps[earth_name]

    # -----------------------------------------------------------------

    def get_total_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the total maps from the earth projection ...")

        # Bolometric luminosity
        if self.has_total_earth_map(bol_map_name): self.total_maps[earth_name][bol_map_name] = self.load_total_earth_map(bol_map_name)
        elif self.model.has_total_bolometric_luminosity_map_earth: self.total_maps[earth_name][bol_map_name] = self.model.total_bolometric_luminosity_map_earth

        # Intrinsic stellar luminosity (transparent luminosity)
        if self.has_total_earth_map(intr_stellar_map_name): self.total_maps[earth_name][intr_stellar_map_name] = self.load_total_earth_map(intr_stellar_map_name)
        elif self.model.has_total_intrinsic_stellar_luminosity_map_earth: self.total_maps[earth_name][intr_stellar_map_name] = self.model.total_intrinsic_stellar_luminosity_map_earth

        # Observed stellar luminosity
        if self.has_total_earth_map(obs_stellar_map_name): self.total_maps[earth_name][obs_stellar_map_name] = self.load_total_earth_map(obs_stellar_map_name)
        elif self.model.has_total_observed_stellar_luminosity_map_earth: self.total_maps[earth_name][obs_stellar_map_name] = self.model.total_observed_stellar_luminosity_map_earth

        # Diffuse dust emission luminosity
        if self.has_total_earth_map(diffuse_dust_map_name): self.total_maps[earth_name][diffuse_dust_map_name] = self.load_total_earth_map(diffuse_dust_map_name)
        elif self.model.has_total_diffuse_dust_luminosity_map_earth: self.total_maps[earth_name][diffuse_dust_map_name] = self.model.total_diffuse_dust_luminosity_map_earth

        # Dust emission luminosity
        if self.has_total_earth_map(dust_map_name): self.total_maps[earth_name][dust_map_name] = self.load_total_earth_map(dust_map_name)
        elif self.model.has_total_dust_luminosity_map_earth: self.total_maps[earth_name][dust_map_name] = self.model.total_dust_luminosity_map_earth

        # Scattered stellar luminosity
        if self.has_total_earth_map(scattered_map_name): self.total_maps[earth_name][scattered_map_name] = self.load_total_earth_map(scattered_map_name)
        elif self.model.has_total_scattered_stellar_luminosity_map_earth: self.total_maps[earth_name][scattered_map_name] = self.model.total_scattered_stellar_luminosity_map_earth

        # Absorbed stellar luminosity (by diffuse dust) (extinction)
        if self.has_total_earth_map(absorbed_diffuse_map_name): self.total_maps[earth_name][absorbed_diffuse_map_name] = self.load_total_earth_map(absorbed_diffuse_map_name)
        elif self.model.has_total_absorbed_diffuse_stellar_luminosity_map_earth: self.total_maps[earth_name][absorbed_diffuse_map_name] = self.model.total_absorbed_diffuse_stellar_luminosity_map_earth

        # Fraction of energy absorbed by DIFFUSE dust
        if self.has_total_earth_map(fabs_diffuse_map_name): self.total_maps[earth_name][fabs_diffuse_map_name] = self.load_total_earth_map(fabs_diffuse_map_name)
        elif self.model.has_total_fabs_diffuse_map_earth: self.total_maps[earth_name][fabs_diffuse_map_name] = self.model.total_fabs_diffuse_map_earth

        # Fraction of energy absorbed by dust
        if self.has_total_earth_map(fabs_map_name): self.total_maps[earth_name][fabs_map_name] = self.load_total_earth_map(fabs_map_name)
        elif self.model.has_total_fabs_map_earth: self.total_maps[earth_name][fabs_map_name] = self.model.total_fabs_map_earth

        # Attenuated stellar luminosity (attenuation)
        if self.has_total_earth_map(attenuated_map_name): self.total_maps[earth_name][attenuated_map_name] = self.load_total_earth_map(attenuated_map_name)
        elif self.model.has_total_attenuated_stellar_luminosity_map_earth: self.total_maps[earth_name][attenuated_map_name] = self.model.total_attenuated_stellar_luminosity_map_earth

        # Direct luminosity
        if self.has_total_earth_map(direct_map_name): self.total_maps[earth_name][direct_map_name] = self.load_total_earth_map(direct_map_name)
        elif self.model.has_total_direct_stellar_luminosity_map_earth: self.total_maps[earth_name][direct_map_name] = self.model.total_direct_stellar_luminosity_map_earth

        # Star formation rate
        if self.has_total_earth_map(sfr_map_name): self.total_maps[earth_name][sfr_map_name] = self.load_total_earth_map(sfr_map_name)
        elif self.model.has_total_star_formation_rate_map_earth: self.total_maps[earth_name][sfr_map_name] = self.model.total_star_formation_rate_map_earth

        # Stellar mass
        if self.has_total_earth_map(stellar_mass_map_name): self.total_maps[earth_name][stellar_mass_map_name] = self.load_total_earth_map(stellar_mass_map_name)
        elif self.model.has_total_stellar_mass_map_earth: self.total_maps[earth_name][stellar_mass_map_name] = self.model.total_stellar_mass_map_earth

        # Specific star formation rate
        if self.has_total_earth_map(ssfr_map_name): self.total_maps[earth_name][ssfr_map_name] = self.load_total_earth_map(ssfr_map_name)
        elif self.model.has_total_ssfr_map_earth: self.total_maps[earth_name][ssfr_map_name] = self.model.total_ssfr_map_earth

    # -----------------------------------------------------------------

    @property
    def total_faceon_maps(self):

        """
        Thisn function ...
        :return:
        """

        return self.total_maps[faceon_name]

    # -----------------------------------------------------------------

    def get_total_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the total maps from the face-on projection ...")

        # Bolometric luminosity
        if self.has_total_faceon_map(bol_map_name): self.total_maps[faceon_name][bol_map_name] = self.load_total_faceon_map(bol_map_name)
        elif self.model.has_total_bolometric_luminosity_map_faceon: self.total_maps[faceon_name][bol_map_name] = self.model.total_bolometric_luminosity_map_faceon

        # Intrinsic stellar luminosity (transparent luminosity)
        if self.has_total_faceon_map(intr_stellar_map_name): self.total_maps[faceon_name][intr_stellar_map_name] = self.load_total_faceon_map(intr_stellar_map_name)
        elif self.model.has_total_intrinsic_stellar_luminosity_map_faceon: self.total_maps[faceon_name][intr_stellar_map_name] = self.model.total_intrinsic_stellar_luminosity_map_faceon

        # Observed stellar luminosity
        if self.has_total_faceon_map(obs_stellar_map_name): self.total_maps[faceon_name][obs_stellar_map_name] = self.load_total_faceon_map(obs_stellar_map_name)
        elif self.model.has_total_observed_stellar_luminosity_map_faceon: self.total_maps[faceon_name][obs_stellar_map_name] = self.model.total_observed_stellar_luminosity_map_faceon

        # Diffuse dust emission luminosity
        if self.has_total_faceon_map(diffuse_dust_map_name): self.total_maps[faceon_name][diffuse_dust_map_name] = self.load_total_faceon_map(diffuse_dust_map_name)
        elif self.model.has_total_diffuse_dust_luminosity_map_faceon: self.total_maps[faceon_name][diffuse_dust_map_name] = self.model.total_diffuse_dust_luminosity_map_faceon

        # Dust emission luminosity
        if self.has_total_faceon_map(dust_map_name): self.total_maps[faceon_name][dust_map_name] = self.load_total_faceon_map(dust_map_name)
        elif self.model.has_total_dust_luminosity_map_faceon: self.total_maps[faceon_name][dust_map_name] = self.model.total_dust_luminosity_map_faceon

        # Scattered stellar luminosity
        if self.has_total_faceon_map(scattered_map_name): self.total_maps[faceon_name][scattered_map_name] = self.load_total_faceon_map(scattered_map_name)
        elif self.model.has_total_scattered_stellar_luminosity_map_faceon: self.total_maps[faceon_name][scattered_map_name] = self.model.total_scattered_stellar_luminosity_map_faceon

        # Absorbed stellar luminosity (by diffuse dust) (extinction)
        if self.has_total_faceon_map(absorbed_diffuse_map_name): self.total_maps[faceon_name][absorbed_diffuse_map_name] = self.load_total_faceon_map(absorbed_diffuse_map_name)
        elif self.model.has_total_absorbed_diffuse_stellar_luminosity_map_faceon: self.total_maps[faceon_name][absorbed_diffuse_map_name] = self.model.total_absorbed_diffuse_stellar_luminosity_map_faceon

        # Fraction of energy absorbed by DIFFUSE dust
        if self.has_total_faceon_map(fabs_diffuse_map_name): self.total_maps[faceon_name][fabs_diffuse_map_name] = self.load_total_faceon_map(fabs_diffuse_map_name)
        elif self.model.has_total_fabs_diffuse_map_faceon: self.total_maps[faceon_name][fabs_diffuse_map_name] = self.model.total_fabs_diffuse_map_faceon

        # Fraction of energy absorbed by dust
        if self.has_total_faceon_map(fabs_map_name): self.total_maps[faceon_name][fabs_map_name] = self.load_total_faceon_map(fabs_map_name)
        elif self.model.has_total_fabs_map_faceon: self.total_maps[faceon_name][fabs_map_name] = self.model.total_fabs_map_faceon

        # Attenuated stellar luminosity (attenuation)
        if self.has_total_faceon_map(attenuated_map_name): self.total_maps[faceon_name][attenuated_map_name] = self.load_total_faceon_map(attenuated_map_name)
        elif self.model.has_total_attenuated_stellar_luminosity_map_faceon: self.total_maps[faceon_name][attenuated_map_name] = self.model.total_attenuated_stellar_luminosity_map_faceon

        # Direct luminosity
        if self.has_total_faceon_map(direct_map_name): self.total_maps[faceon_name][direct_map_name] = self.load_total_faceon_map(direct_map_name)
        elif self.model.has_total_direct_stellar_luminosity_map_faceon: self.total_maps[faceon_name][direct_map_name] = self.model.total_direct_stellar_luminosity_map_faceon

        # Star formation rate
        if self.has_total_faceon_map(sfr_map_name): self.total_maps[faceon_name][sfr_map_name] = self.load_total_faceon_map(sfr_map_name)
        elif self.model.has_total_star_formation_rate_map_faceon: self.total_maps[faceon_name][sfr_map_name] = self.model.total_star_formation_rate_map_faceon

        # Stellar mass
        if self.has_total_faceon_map(stellar_mass_map_name): self.total_maps[faceon_name][stellar_mass_map_name] = self.load_total_faceon_map(stellar_mass_map_name)
        elif self.model.has_total_stellar_mass_map_faceon: self.total_maps[faceon_name][stellar_mass_map_name] = self.model.total_stellar_mass_map_faceon

        # Specific star formation rate
        if self.has_total_faceon_map(ssfr_map_name): self.total_maps[faceon_name][ssfr_map_name] = self.load_total_faceon_map(ssfr_map_name)
        elif self.model.has_total_ssfr_map_faceon: self.total_maps[faceon_name][ssfr_map_name] = self.model.total_ssfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def total_edgeon_maps(self):

        """
        This function ...
        :return:
        """

        return self.total_maps[edgeon_name]

    # -----------------------------------------------------------------

    def get_total_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the total maps from the edge-on projection ...")

        # Bolometric luminosity
        if self.has_total_edgeon_map(bol_map_name): self.total_maps[edgeon_name][bol_map_name] = self.load_total_edgeon_map(bol_map_name)
        elif self.model.has_total_bolometric_luminosity_map_edgeon: self.total_maps[edgeon_name][bol_map_name] = self.model.total_bolometric_luminosity_map_edgeon

        # Intrinsic stellar luminosity (transparent luminosity)
        if self.has_total_edgeon_map(intr_stellar_map_name): self.total_maps[edgeon_name][intr_stellar_map_name] = self.load_total_edgeon_map(intr_stellar_map_name)
        elif self.model.has_total_intrinsic_stellar_luminosity_map_edgeon: self.total_maps[edgeon_name][intr_stellar_map_name] = self.model.total_intrinsic_stellar_luminosity_map_edgeon

        # Observed stellar luminosity
        if self.has_total_edgeon_map(obs_stellar_map_name): self.total_maps[edgeon_name][obs_stellar_map_name] = self.load_total_edgeon_map(obs_stellar_map_name)
        elif self.model.has_total_observed_stellar_luminosity_map_edgeon: self.total_maps[edgeon_name][obs_stellar_map_name] = self.model.total_observed_stellar_luminosity_map_edgeon

        # Diffuse dust emission luminosity
        if self.has_total_edgeon_map(diffuse_dust_map_name): self.total_maps[edgeon_name][diffuse_dust_map_name] = self.load_total_edgeon_map(diffuse_dust_map_name)
        elif self.model.has_total_diffuse_dust_luminosity_map_edgeon: self.total_maps[edgeon_name][diffuse_dust_map_name] = self.model.total_diffuse_dust_luminosity_map_edgeon

        # Dust emission luminosity
        if self.has_total_edgeon_map(dust_map_name): self.total_maps[edgeon_name][dust_map_name] = self.load_total_edgeon_map(dust_map_name)
        elif self.model.has_total_dust_luminosity_map_edgeon: self.total_maps[edgeon_name][dust_map_name] = self.model.total_dust_luminosity_map_edgeon

        # Scattered stellar luminosity
        if self.has_total_edgeon_map(scattered_map_name): self.total_maps[edgeon_name][scattered_map_name] = self.load_total_edgeon_map(scattered_map_name)
        elif self.model.has_total_scattered_stellar_luminosity_map_edgeon: self.total_maps[edgeon_name][scattered_map_name] = self.model.total_scattered_stellar_luminosity_map_edgeon

        # Absorbed stellar luminosity (by diffuse dust) (extinction)
        if self.has_total_edgeon_map(absorbed_diffuse_map_name): self.total_maps[edgeon_name][absorbed_diffuse_map_name] = self.load_total_edgeon_map(absorbed_diffuse_map_name)
        elif self.model.has_total_absorbed_diffuse_stellar_luminosity_map_edgeon: self.total_maps[edgeon_name][absorbed_diffuse_map_name] = self.model.total_absorbed_diffuse_stellar_luminosity_map_edgeon

        # Fraction of energy absorbed by DIFFUSE dust
        if self.has_total_edgeon_map(fabs_diffuse_map_name): self.total_maps[edgeon_name][fabs_diffuse_map_name] = self.load_total_edgeon_map(fabs_diffuse_map_name)
        elif self.model.has_total_fabs_diffuse_map_edgeon: self.total_maps[edgeon_name][fabs_diffuse_map_name] = self.model.total_fabs_diffuse_map_edgeon

        # Fraction of energy absorbed by dust
        if self.has_total_edgeon_map(fabs_map_name): self.total_maps[edgeon_name][fabs_map_name] = self.load_total_edgeon_map(fabs_map_name)
        elif self.model.has_total_fabs_map_edgeon: self.total_maps[edgeon_name][fabs_map_name] = self.model.total_fabs_map_edgeon

        # Attenuated stellar luminosity (attenuation)
        if self.has_total_edgeon_map(attenuated_map_name): self.total_maps[edgeon_name][attenuated_map_name] = self.load_total_edgeon_map(attenuated_map_name)
        elif self.model.has_total_attenuated_stellar_luminosity_map_edgeon: self.total_maps[edgeon_name][attenuated_map_name] = self.model.total_attenuated_stellar_luminosity_map_edgeon

        # Direct luminosity
        if self.has_total_edgeon_map(direct_map_name): self.total_maps[edgeon_name][direct_map_name] = self.load_total_edgeon_map(direct_map_name)
        elif self.model.has_total_direct_stellar_luminosity_map_edgeon: self.total_maps[edgeon_name][direct_map_name] = self.model.total_direct_stellar_luminosity_map_edgeon

        # Star formation rate
        if self.has_total_edgeon_map(sfr_map_name): self.total_maps[edgeon_name][sfr_map_name] = self.load_total_edgeon_map(sfr_map_name)
        elif self.model.has_total_star_formation_rate_map_edgeon: self.total_maps[edgeon_name][sfr_map_name] = self.model.total_star_formation_rate_map_edgeon

        # Stellar mass
        if self.has_total_edgeon_map(stellar_mass_map_name): self.total_maps[edgeon_name][stellar_mass_map_name] = self.load_total_edgeon_map(stellar_mass_map_name)
        elif self.model.has_total_stellar_mass_map_edgeon: self.total_maps[edgeon_name][stellar_mass_map_name] = self.model.total_stellar_mass_map_edgeon

        # Specific star formation rate
        if self.has_total_edgeon_map(ssfr_map_name): self.total_maps[edgeon_name][ssfr_map_name] = self.load_total_edgeon_map(ssfr_map_name)
        elif self.model.has_total_ssfr_map_edgeon: self.total_maps[edgeon_name][ssfr_map_name] = self.model.total_ssfr_map_edgeon

    # -----------------------------------------------------------------

    def get_bulge_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.get_bulge_maps_earth()

        # Face-on
        if self.do_faceon: self.get_bulge_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.get_bulge_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def bulge_earth_maps(self):

        """
        This function ...
        :return:
        """

        return self.bulge_maps[earth_name]

    # -----------------------------------------------------------------

    def get_bulge_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Bolometric luminosity
        if self.has_bulge_earth_map(bol_map_name): self.bulge_maps[earth_name][bol_map_name] = self.load_bulge_earth_map(bol_map_name)
        elif self.model.has_old_bulge_bolometric_luminosity_map_earth: self.bulge_maps[earth_name][bol_map_name] = self.model.old_bulge_bolometric_luminosity_map_earth

        # Direct
        if self.has_bulge_earth_map(direct_map_name): self.bulge_maps[earth_name][direct_map_name] = self.load_bulge_earth_map(direct_map_name)
        elif self.model.has_old_bulge_direct_stellar_luminosity_map_earth: self.bulge_maps[earth_name][direct_map_name] = self.model.old_bulge_direct_stellar_luminosity_map_earth

        # (observed) I1 lum
        if self.has_bulge_earth_map(i1_map_name): self.bulge_maps[earth_name][i1_map_name] = self.load_bulge_earth_map(i1_map_name)
        elif self.model.has_old_bulge_i1_luminosity_map_earth: self.bulge_maps[earth_name][i1_map_name] = self.model.old_bulge_i1_luminosity_map_earth

        # Intrinsic I1
        if self.has_bulge_earth_map(intr_i1_map_name): self.bulge_maps[earth_name][intr_i1_map_name] = self.load_bulge_earth_map(intr_i1_map_name)
        elif self.model.has_old_bulge_intrinsic_i1_luminosity_map_earth: self.bulge_maps[earth_name][intr_i1_map_name] = self.model.old_bulge_intrinsic_i1_luminosity_map_earth

        # Dust luminosity
        if self.has_bulge_earth_map(dust_map_name): self.bulge_maps[earth_name][dust_map_name] = self.load_bulge_earth_map(dust_map_name)
        elif self.model.has_old_bulge_dust_luminosity_map_earth: self.bulge_maps[earth_name][dust_map_name] = self.model.old_bulge_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def bulge_faceon_maps(self):

        """
        This function ...
        :return:
        """

        return self.bulge_maps[faceon_name]

    # -----------------------------------------------------------------

    def get_bulge_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Bolometric luminosity
        if self.has_bulge_faceon_map(bol_map_name): self.bulge_maps[faceon_name][bol_map_name] = self.load_bulge_faceon_map(bol_map_name)
        elif self.model.has_old_bulge_bolometric_luminosity_map_faceon: self.bulge_maps[faceon_name][bol_map_name] = self.model.old_bulge_bolometric_luminosity_map_faceon

        # Direct
        if self.has_bulge_faceon_map(direct_map_name): self.bulge_maps[faceon_name][direct_map_name] = self.load_bulge_faceon_map(direct_map_name)
        elif self.model.has_old_bulge_direct_stellar_luminosity_map_faceon: self.bulge_maps[faceon_name][direct_map_name] = self.model.old_bulge_direct_stellar_luminosity_map_faceon

        # (observed) I1 lum
        if self.has_bulge_faceon_map(i1_map_name): self.bulge_maps[faceon_name][i1_map_name] = self.load_bulge_faceon_map(i1_map_name)
        elif self.model.has_old_bulge_i1_luminosity_map_faceon: self.bulge_maps[faceon_name][i1_map_name] = self.model.old_bulge_i1_luminosity_map_faceon

        # Intrinsic I1
        if self.has_bulge_faceon_map(intr_i1_map_name): self.bulge_maps[faceon_name][intr_i1_map_name] = self.load_bulge_faceon_map(intr_i1_map_name)
        elif self.model.has_old_bulge_intrinsic_i1_luminosity_map_faceon: self.bulge_maps[faceon_name][intr_i1_map_name] = self.model.old_bulge_intrinsic_i1_luminosity_map_faceon

        # Dust luminosity
        if self.has_bulge_faceon_map(dust_map_name): self.bulge_maps[faceon_name][dust_map_name] = self.load_bulge_faceon_map(dust_map_name)
        elif self.model.has_old_bulge_dust_luminosity_map_faceon: self.bulge_maps[faceon_name][dust_map_name] = self.model.old_bulge_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def bulge_edgeon_maps(self):

        """
        This function ...
        :return:
        """

        return self.bulge_maps[edgeon_name]

    # -----------------------------------------------------------------

    def get_bulge_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Bolometric luminosity
        if self.has_bulge_edgeon_map(bol_map_name): self.bulge_maps[edgeon_name][bol_map_name] = self.load_bulge_edgeon_map(bol_map_name)
        elif self.model.has_old_bulge_bolometric_luminosity_map_edgeon: self.bulge_maps[edgeon_name][bol_map_name] = self.model.old_bulge_bolometric_luminosity_map_edgeon

        # Direct
        if self.has_bulge_edgeon_map(direct_map_name): self.bulge_maps[edgeon_name][direct_map_name] = self.load_bulge_edgeon_map(direct_map_name)
        elif self.model.has_old_bulge_direct_stellar_luminosity_map_edgeon: self.bulge_maps[edgeon_name][direct_map_name] = self.model.old_bulge_direct_stellar_luminosity_map_edgeon

        # (observed) I1 lum
        if self.has_bulge_edgeon_map(i1_map_name): self.bulge_maps[edgeon_name][i1_map_name] = self.load_bulge_edgeon_map(i1_map_name)
        elif self.model.has_old_bulge_i1_luminosity_map_edgeon: self.bulge_maps[edgeon_name][i1_map_name] = self.model.old_bulge_i1_luminosity_map_edgeon

        # Intrinsic I1
        if self.has_bulge_edgeon_map(intr_i1_map_name): self.bulge_maps[edgeon_name][intr_i1_map_name] = self.load_bulge_edgeon_map(intr_i1_map_name)
        elif self.model.has_old_bulge_intrinsic_i1_luminosity_map_edgeon: self.bulge_maps[edgeon_name][intr_i1_map_name] = self.model.old_bulge_intrinsic_i1_luminosity_map_edgeon

        # Dust luminosity
        if self.has_bulge_edgeon_map(dust_map_name): self.bulge_maps[edgeon_name][dust_map_name] = self.load_bulge_edgeon_map(dust_map_name)
        elif self.model.has_old_bulge_dust_luminosity_map_edgeon: self.bulge_maps[edgeon_name][dust_map_name] = self.model.old_bulge_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------

    def get_disk_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.get_disk_maps_earth()

        # Face-on
        if self.do_faceon: self.get_disk_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.get_disk_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def disk_earth_maps(self):

        """
        This function ...
        :return:
        """

        return self.disk_maps[earth_name]

    # -----------------------------------------------------------------

    def get_disk_maps_earth(self):

        """
        This unction ...
        :return:
        """

        # Bolometric
        if self.has_disk_earth_map(bol_map_name): self.disk_maps[earth_name][bol_map_name] = self.load_disk_earth_map(bol_map_name)
        elif self.model.has_old_disk_bolometric_luminosity_map_earth: self.disk_maps[earth_name][bol_map_name] = self.model.old_disk_bolometric_luminosity_map_earth

        # Direct
        if self.has_disk_earth_map(direct_map_name): self.disk_maps[earth_name][direct_map_name] = self.load_disk_earth_map(direct_map_name)
        elif self.model.has_old_disk_direct_stellar_luminosity_map_earth: self.disk_maps[earth_name][direct_map_name] = self.model.old_disk_direct_stellar_luminosity_map_earth

        # (observed) I1
        if self.has_disk_earth_map(i1_map_name): self.disk_maps[earth_name][i1_map_name] = self.load_disk_earth_map(i1_map_name)
        elif self.model.has_old_disk_i1_luminosity_map_earth: self.disk_maps[earth_name][i1_map_name] = self.model.old_disk_i1_luminosity_map_earth

        # Intrinsic I1
        if self.has_disk_earth_map(intr_i1_map_name): self.disk_maps[earth_name][intr_i1_map_name] = self.load_disk_earth_map(intr_i1_map_name)
        elif self.model.has_old_disk_intrinsic_i1_luminosity_map_earth: self.disk_maps[earth_name][intr_i1_map_name] = self.model.old_disk_intrinsic_i1_luminosity_map_earth

        # Dust luminosity
        if self.has_disk_earth_map(dust_map_name): self.disk_maps[earth_name][dust_map_name] = self.load_disk_earth_map(dust_map_name)
        elif self.model.has_old_disk_dust_luminosity_map_earth: self.disk_maps[earth_name][dust_map_name] = self.model.old_disk_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def disk_faceon_maps(self):

        """
        This function ...
        :return:
        """

        return self.disk_maps[faceon_name]

    # -----------------------------------------------------------------

    def get_disk_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_disk_faceon_map(bol_map_name): self.disk_maps[faceon_name][bol_map_name] = self.load_disk_faceon_map(bol_map_name)
        elif self.model.has_old_disk_bolometric_luminosity_map_faceon: self.disk_maps[faceon_name][bol_map_name] = self.model.old_disk_bolometric_luminosity_map_faceon

        # Direct
        if self.has_disk_faceon_map(direct_map_name): self.disk_maps[faceon_name][direct_map_name] = self.load_disk_faceon_map(direct_map_name)
        elif self.model.has_old_disk_direct_stellar_luminosity_map_faceon: self.disk_maps[faceon_name][direct_map_name] = self.model.old_disk_direct_stellar_luminosity_map_faceon

        # (observed) I1
        if self.has_disk_faceon_map(i1_map_name): self.disk_maps[faceon_name][i1_map_name] = self.load_disk_faceon_map(i1_map_name)
        elif self.model.has_old_disk_i1_luminosity_map_faceon: self.disk_maps[faceon_name][i1_map_name] = self.model.old_disk_i1_luminosity_map_faceon

        # Intrinsic I1
        if self.has_disk_faceon_map(intr_i1_map_name): self.disk_maps[faceon_name][intr_i1_map_name] = self.load_disk_faceon_map(intr_i1_map_name)
        elif self.model.has_old_disk_intrinsic_i1_luminosity_map_faceon: self.disk_maps[faceon_name][intr_i1_map_name] = self.model.old_disk_intrinsic_i1_luminosity_map_faceon

        # Dust luminosity
        if self.has_disk_faceon_map(dust_map_name): self.disk_maps[faceon_name][dust_map_name] = self.load_disk_faceon_map(dust_map_name)
        elif self.model.has_old_disk_dust_luminosity_map_faceon: self.disk_maps[faceon_name][dust_map_name] = self.model.old_disk_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def disk_edgeon_maps(self):

        """
        This function ...
        :return:
        """

        return self.disk_maps[edgeon_name]

    # -----------------------------------------------------------------

    def get_disk_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_disk_edgeon_map(bol_map_name): self.disk_maps[edgeon_name][bol_map_name] = self.load_disk_edgeon_map(bol_map_name)
        elif self.model.has_old_disk_bolometric_luminosity_map_edgeon: self.disk_maps[edgeon_name][bol_map_name] = self.model.old_disk_bolometric_luminosity_map_edgeon

        # Direct
        if self.has_disk_edgeon_map(direct_map_name): self.disk_maps[edgeon_name][direct_map_name] = self.load_disk_edgeon_map(direct_map_name)
        elif self.model.has_old_disk_direct_stellar_luminosity_map_edgeon: self.disk_maps[edgeon_name][direct_map_name] = self.model.old_disk_direct_stellar_luminosity_map_edgeon

        # (observed) I1
        if self.has_disk_edgeon_map(i1_map_name): self.disk_maps[edgeon_name][i1_map_name] = self.load_disk_edgeon_map(i1_map_name)
        elif self.model.has_old_disk_i1_luminosity_map_edgeon: self.disk_maps[edgeon_name][i1_map_name] = self.model.old_disk_i1_luminosity_map_edgeon

        # Intrinsic I1
        if self.has_disk_edgeon_map(intr_i1_map_name): self.disk_maps[edgeon_name][intr_i1_map_name] = self.load_disk_edgeon_map(intr_i1_map_name)
        elif self.model.has_old_disk_intrinsic_i1_luminosity_map_edgeon: self.disk_maps[edgeon_name][intr_i1_map_name] = self.model.old_disk_intrinsic_i1_luminosity_map_edgeon

        # Dust luminosity
        if self.has_disk_edgeon_map(dust_map_name): self.disk_maps[edgeon_name][dust_map_name] = self.load_disk_edgeon_map(dust_map_name)
        elif self.model.has_old_disk_dust_luminosity_map_edgeon: self.disk_maps[edgeon_name][dust_map_name] = self.model.old_disk_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------

    def get_old_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.get_old_maps_earth()

        # Face-on
        if self.do_faceon: self.get_old_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.get_old_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def old_earth_maps(self):

        """
        This function ...
        :return:
        """

        return self.old_maps[earth_name]

    # -----------------------------------------------------------------

    def get_old_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_old_earth_map(bol_map_name): self.old_maps[earth_name][bol_map_name] = self.load_old_earth_map(bol_map_name)
        elif self.model.has_old_bolometric_luminosity_map_earth: self.old_maps[earth_name][bol_map_name] = self.model.old_bolometric_luminosity_map_earth

        # Direct
        if self.has_old_earth_map(direct_map_name): self.old_maps[earth_name][direct_map_name] = self.load_old_earth_map(direct_map_name)
        elif self.model.has_old_direct_stellar_luminosity_map_earth: self.old_maps[earth_name][direct_map_name] = self.model.old_direct_stellar_luminosity_map_earth

        # (observed) I1
        if self.has_old_earth_map(i1_map_name): self.old_maps[earth_name][i1_map_name] = self.load_old_earth_map(i1_map_name)
        elif self.model.has_old_i1_luminosity_map_earth: self.old_maps[earth_name][i1_map_name] = self.model.old_i1_luminosity_map_earth

        # Intrinsic I1
        if self.has_old_earth_map(intr_i1_map_name): self.old_maps[earth_name][intr_i1_map_name] = self.load_old_earth_map(intr_i1_map_name)
        elif self.model.has_old_intrinsic_i1_luminosity_map_earth: self.old_maps[earth_name][intr_i1_map_name] = self.model.old_intrinsic_i1_luminosity_map_earth

        # Dust luminosity
        if self.has_old_earth_map(dust_map_name): self.old_maps[earth_name][dust_map_name] = self.load_old_earth_map(dust_map_name)
        elif self.model.has_old_dust_luminosity_map_earth: self.old_maps[earth_name][dust_map_name] = self.model.old_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def old_faceon_maps(self):

        """
        This function ...
        :return:
        """

        return self.old_maps[faceon_name]

    # -----------------------------------------------------------------

    def get_old_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_old_faceon_map(bol_map_name): self.old_maps[faceon_name][bol_map_name] = self.load_old_faceon_map(bol_map_name)
        elif self.model.has_old_bolometric_luminosity_map_faceon: self.old_maps[faceon_name][bol_map_name] = self.model.old_bolometric_luminosity_map_faceon

        # Direct
        if self.has_old_faceon_map(direct_map_name): self.old_maps[faceon_name][direct_map_name] = self.load_old_faceon_map(direct_map_name)
        elif self.model.has_old_direct_stellar_luminosity_map_faceon: self.old_maps[faceon_name][direct_map_name] = self.model.old_direct_stellar_luminosity_map_faceon

        # (observed) I1
        if self.has_old_faceon_map(i1_map_name): self.old_maps[faceon_name][i1_map_name] = self.load_old_faceon_map(i1_map_name)
        elif self.model.has_old_i1_luminosity_map_faceon: self.old_maps[faceon_name][i1_map_name] = self.model.old_i1_luminosity_map_faceon

        # Intrinsic I1
        if self.has_old_faceon_map(intr_i1_map_name): self.old_maps[faceon_name][intr_i1_map_name] = self.load_old_faceon_map(intr_i1_map_name)
        elif self.model.has_old_intrinsic_i1_luminosity_map_faceon: self.old_maps[faceon_name][intr_i1_map_name] = self.model.old_intrinsic_i1_luminosity_map_faceon

        # Dust luminosity
        if self.has_old_faceon_map(dust_map_name): self.old_maps[faceon_name][dust_map_name] = self.load_old_faceon_map(dust_map_name)
        elif self.model.has_old_dust_luminosity_map_faceon: self.old_maps[faceon_name][dust_map_name] = self.model.old_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def old_edgeon_maps(self):

        """
        This function ...
        :return:
        """

        return self.old_maps[edgeon_name]

    # -----------------------------------------------------------------

    def get_old_maps_edgeon(self):

        """
        This ufnction ...
        :return:
        """

        # Bolometric
        if self.has_old_edgeon_map(bol_map_name): self.old_maps[edgeon_name][bol_map_name] = self.load_old_edgeon_map(bol_map_name)
        elif self.model.has_old_bolometric_luminosity_map_edgeon: self.old_maps[edgeon_name][bol_map_name] = self.model.old_bolometric_luminosity_map_edgeon

        # Direct
        if self.has_old_edgeon_map(direct_map_name): self.old_maps[edgeon_name][direct_map_name] = self.load_old_edgeon_map(direct_map_name)
        elif self.model.has_old_direct_stellar_luminosity_map_edgeon: self.old_maps[edgeon_name][direct_map_name] = self.model.old_direct_stellar_luminosity_map_edgeon

        # (observed) I1
        if self.has_old_edgeon_map(i1_map_name): self.old_maps[edgeon_name][i1_map_name] = self.load_old_edgeon_map(i1_map_name)
        elif self.model.has_old_i1_luminosity_map_edgeon: self.old_maps[edgeon_name][i1_map_name] = self.model.old_i1_luminosity_map_edgeon

        # Intrinsic I1
        if self.has_old_edgeon_map(intr_i1_map_name): self.old_maps[edgeon_name][intr_i1_map_name] = self.load_old_edgeon_map(intr_i1_map_name)
        elif self.model.has_old_intrinsic_i1_luminosity_map_edgeon: self.old_maps[edgeon_name][intr_i1_map_name] = self.model.old_intrinsic_i1_luminosity_map_edgeon

        # Dust luminosity
        if self.has_old_edgeon_map(dust_map_name): self.old_maps[edgeon_name][dust_map_name] = self.load_old_edgeon_map(dust_map_name)
        elif self.model.has_old_dust_luminosity_map_edgeon: self.old_maps[edgeon_name][dust_map_name] = self.model.old_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------

    def get_young_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.get_young_maps_earth()

        # Face-on
        if self.do_faceon: self.get_young_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.get_young_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def young_earth_maps(self):

        """
        This function ...
        :return:
        """

        return self.young_maps[earth_name]

    # -----------------------------------------------------------------

    def get_young_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_young_earth_map(bol_map_name): self.young_maps[earth_name][bol_map_name] = self.load_young_earth_map(bol_map_name)
        elif self.model.has_young_bolometric_luminosity_map_earth: self.young_maps[earth_name][bol_map_name] = self.model.young_bolometric_luminosity_map_earth

        # Direct
        if self.has_young_earth_map(direct_map_name): self.young_maps[earth_name][direct_map_name] = self.load_young_earth_map(direct_map_name)
        elif self.model.has_young_direct_stellar_luminosity_map_earth: self.young_maps[earth_name][direct_map_name] = self.model.young_direct_stellar_luminosity_map_earth

        # (observed) FUV
        if self.has_young_earth_map(fuv_map_name): self.young_maps[earth_name][fuv_map_name] = self.load_young_earth_map(fuv_map_name)
        elif self.model.has_young_fuv_luminosity_map_earth: self.young_maps[earth_name][fuv_map_name] = self.model.young_fuv_luminosity_map_earth

        # Intrinsic FUV
        if self.has_young_earth_map(intr_fuv_map_name): self.young_maps[earth_name][intr_fuv_map_name] = self.load_young_earth_map(intr_fuv_map_name)
        elif self.model.has_young_intrinsic_fuv_luminosity_map_earth: self.young_maps[earth_name][intr_fuv_map_name] = self.model.young_intrinsic_fuv_luminosity_map_earth

        # Dust luminosity
        if self.has_young_earth_map(dust_map_name): self.young_maps[earth_name][dust_map_name] = self.load_young_earth_map(dust_map_name)
        elif self.model.has_young_dust_luminosity_map_earth: self.young_maps[earth_name][dust_map_name] = self.model.young_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def young_faceon_maps(self):

        """
        This function ...
        :return:
        """

        return self.young_maps[faceon_name]

    # -----------------------------------------------------------------

    def get_young_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_young_faceon_map(bol_map_name): self.young_maps[faceon_name][bol_map_name] = self.load_young_faceon_map(bol_map_name)
        elif self.model.has_young_bolometric_luminosity_map_faceon: self.young_maps[faceon_name][bol_map_name] = self.model.young_bolometric_luminosity_map_faceon

        # Direct
        if self.has_young_faceon_map(direct_map_name): self.young_maps[faceon_name][direct_map_name] = self.load_young_faceon_map(direct_map_name)
        elif self.model.has_young_direct_stellar_luminosity_map_faceon: self.young_maps[faceon_name][direct_map_name] = self.model.young_direct_stellar_luminosity_map_faceon

        # (observed) FUV
        if self.has_young_faceon_map(fuv_map_name): self.young_maps[faceon_name][fuv_map_name] = self.load_young_faceon_map(fuv_map_name)
        elif self.model.has_young_fuv_luminosity_map_faceon: self.young_maps[faceon_name][fuv_map_name] = self.model.young_fuv_luminosity_map_faceon

        # Intrinsic FUV
        if self.has_young_faceon_map(intr_fuv_map_name): self.young_maps[faceon_name][intr_fuv_map_name] = self.load_young_faceon_map(intr_fuv_map_name)
        elif self.model.has_young_intrinsic_fuv_luminosity_map_faceon: self.young_maps[faceon_name][intr_fuv_map_name] = self.model.young_intrinsic_fuv_luminosity_map_faceon

        # Dust luminosity
        if self.has_young_faceon_map(dust_map_name): self.young_maps[faceon_name][dust_map_name] = self.load_young_faceon_map(dust_map_name)
        elif self.model.has_young_dust_luminosity_map_faceon: self.young_maps[faceon_name][dust_map_name] = self.model.young_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def young_edgeon_maps(self):

        """
        This function ...
        :return:
        """

        return self.young_maps[edgeon_name]

    # -----------------------------------------------------------------

    def get_young_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_young_edgeon_map(bol_map_name): self.young_maps[edgeon_name][bol_map_name] = self.load_young_edgeon_map(bol_map_name)
        elif self.model.has_young_bolometric_luminosity_map_edgeon: self.young_maps[edgeon_name][bol_map_name] = self.model.young_bolometric_luminosity_map_edgeon

        # Direct
        if self.has_young_edgeon_map(direct_map_name): self.young_maps[edgeon_name][direct_map_name] = self.load_young_edgeon_map(direct_map_name)
        elif self.model.has_young_direct_stellar_luminosity_map_edgeon: self.young_maps[edgeon_name][direct_map_name] = self.model.young_direct_stellar_luminosity_map_edgeon

        # (observed) FUV
        if self.has_young_edgeon_map(fuv_map_name): self.young_maps[edgeon_name][fuv_map_name] = self.load_young_edgeon_map(fuv_map_name)
        elif self.model.has_young_fuv_luminosity_map_edgeon: self.young_maps[edgeon_name][fuv_map_name] = self.model.young_fuv_luminosity_map_edgeon

        # Intrinsic FUV
        if self.has_young_edgeon_map(intr_fuv_map_name): self.young_maps[edgeon_name][intr_fuv_map_name] = self.load_young_edgeon_map(intr_fuv_map_name)
        elif self.model.has_young_intrinsic_fuv_luminosity_map_edgeon: self.young_maps[edgeon_name][intr_fuv_map_name] = self.model.young_intrinsic_fuv_luminosity_map_edgeon

        # Dust luminosity
        if self.has_young_edgeon_map(dust_map_name): self.young_maps[edgeon_name][dust_map_name] = self.load_young_edgeon_map(dust_map_name)
        elif self.model.has_young_dust_luminosity_map_edgeon: self.young_maps[edgeon_name][dust_map_name] = self.model.young_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------

    def get_sfr_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.get_sfr_maps_earth()

        # Face-on
        if self.do_faceon: self.get_sfr_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.get_sfr_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def sfr_earth_maps(self):

        """
        This function ...
        :return:
        """

        return self.sfr_maps[earth_name]

    # -----------------------------------------------------------------

    def get_sfr_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_sfr_earth_map(bol_map_name): self.sfr_maps[earth_name][bol_map_name] = self.load_sfr_earth_map(bol_map_name)
        elif self.model.has_sfr_bolometric_luminosity_map_earth: self.sfr_maps[earth_name][bol_map_name] = self.model.sfr_bolometric_luminosity_map_earth

        # Direct
        if self.has_sfr_earth_map(direct_map_name): self.sfr_maps[earth_name][direct_map_name] = self.load_sfr_earth_map(direct_map_name)
        elif self.model.has_sfr_direct_stellar_luminosity_map_earth: self.sfr_maps[earth_name][direct_map_name] = self.model.sfr_direct_stellar_luminosity_map_earth

        # (observed) FUV
        if self.has_sfr_earth_map(fuv_map_name): self.sfr_maps[earth_name][fuv_map_name] = self.load_sfr_earth_map(fuv_map_name)
        elif self.model.has_sfr_fuv_luminosity_map_earth: self.sfr_maps[earth_name][fuv_map_name] = self.model.sfr_fuv_luminosity_map_earth

        # Intrinsic FUV
        if self.has_sfr_earth_map(intr_fuv_map_name): self.sfr_maps[earth_name][intr_fuv_map_name] = self.load_sfr_earth_map(intr_fuv_map_name)
        elif self.model.has_sfr_intrinsic_fuv_luminosity_map_earth: self.sfr_maps[earth_name][intr_fuv_map_name] = self.model.sfr_intrinsic_fuv_luminosity_map_earth

        # SFR
        if self.has_sfr_earth_map(sfr_map_name): self.sfr_maps[earth_name][sfr_map_name] = self.load_sfr_earth_map(sfr_map_name)
        elif self.model.has_star_formation_rate_map_earth: self.sfr_maps[earth_name][sfr_map_name] = self.model.star_formation_rate_map_earth

        # Dust mass
        if self.has_sfr_earth_map(dust_mass_map_name): self.sfr_maps[earth_name][dust_mass_map_name] = self.load_sfr_earth_map(dust_mass_map_name)
        elif self.model.has_sfr_dust_mass_map_earth: self.sfr_maps[earth_name][dust_mass_map_name] = self.model.sfr_dust_mass_map_earth

        # Stellar bolometric luminosity
        if self.has_sfr_earth_map(stellar_lum_map_name): self.sfr_maps[earth_name][stellar_lum_map_name] = self.load_sfr_earth_map(stellar_lum_map_name)
        elif self.model.has_sfr_stellar_luminosity_map_earth: self.sfr_maps[earth_name][stellar_lum_map_name] = self.model.sfr_stellar_luminosity_map_earth

        # Intrinsic dust luminosity
        if self.has_sfr_earth_map(intr_dust_map_name): self.sfr_maps[earth_name][intr_dust_map_name] = self.load_sfr_earth_map(intr_dust_map_name)
        elif self.model.has_sfr_intrinsic_dust_luminosity_map_earth: self.sfr_maps[earth_name][intr_dust_map_name] = self.model.sfr_intrinsic_dust_luminosity_map_earth

        # Dust bolometric luminosity
        if self.has_sfr_earth_map(dust_map_name): self.sfr_maps[earth_name][dust_map_name] = self.load_sfr_earth_map(dust_map_name)
        elif self.model.has_sfr_dust_luminosity_map_earth: self.sfr_maps[earth_name][dust_map_name] = self.model.sfr_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def sfr_faceon_maps(self):

        """
        This function ...
        :return:
        """

        return self.sfr_maps[faceon_name]

    # -----------------------------------------------------------------

    def get_sfr_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_sfr_faceon_map(bol_map_name): self.sfr_maps[faceon_name][bol_map_name] = self.load_sfr_faceon_map(bol_map_name)
        elif self.model.has_sfr_bolometric_luminosity_map_faceon: self.sfr_maps[faceon_name][bol_map_name] = self.model.sfr_bolometric_luminosity_map_faceon

        # Direct
        if self.has_sfr_faceon_map(direct_map_name): self.sfr_maps[faceon_name][direct_map_name] = self.load_sfr_faceon_map(direct_map_name)
        elif self.model.has_sfr_direct_stellar_luminosity_map_faceon: self.sfr_maps[faceon_name][direct_map_name] = self.model.sfr_direct_stellar_luminosity_map_faceon

        # (observed) FUV
        if self.has_sfr_faceon_map(fuv_map_name): self.sfr_maps[faceon_name][fuv_map_name] = self.load_sfr_faceon_map(fuv_map_name)
        elif self.model.has_sfr_fuv_luminosity_map_faceon: self.sfr_maps[faceon_name][fuv_map_name] = self.model.sfr_fuv_luminosity_map_faceon

        # Intrinsic FUV
        if self.has_sfr_faceon_map(intr_fuv_map_name): self.sfr_maps[faceon_name][intr_fuv_map_name] = self.load_sfr_faceon_map(intr_fuv_map_name)
        elif self.model.has_sfr_intrinsic_fuv_luminosity_map_faceon: self.sfr_maps[faceon_name][intr_fuv_map_name] = self.model.sfr_intrinsic_fuv_luminosity_map_faceon

        # SFR
        if self.has_sfr_faceon_map(sfr_map_name): self.sfr_maps[faceon_name][sfr_map_name] = self.load_sfr_faceon_map(sfr_map_name)
        elif self.model.has_star_formation_rate_map_faceon: self.sfr_maps[faceon_name][sfr_map_name] = self.model.star_formation_rate_map_faceon

        # Dust mass
        if self.has_sfr_faceon_map(dust_mass_map_name): self.sfr_maps[faceon_name][dust_mass_map_name] = self.load_sfr_faceon_map(dust_mass_map_name)
        elif self.model.has_sfr_dust_mass_map_faceon: self.sfr_maps[faceon_name][dust_mass_map_name] = self.model.sfr_dust_mass_map_faceon

        # Stellar bolometric luminosity
        if self.has_sfr_faceon_map(stellar_lum_map_name): self.sfr_maps[faceon_name][stellar_lum_map_name] = self.load_sfr_faceon_map(stellar_lum_map_name)
        elif self.model.has_sfr_stellar_luminosity_map_faceon: self.sfr_maps[faceon_name][stellar_lum_map_name] = self.model.sfr_stellar_luminosity_map_faceon

        # Intrinsic dust luminosity
        if self.has_sfr_faceon_map(intr_dust_map_name): self.sfr_maps[faceon_name][intr_dust_map_name] = self.load_sfr_faceon_map(intr_dust_map_name)
        elif self.model.has_sfr_intrinsic_dust_luminosity_map_faceon: self.sfr_maps[faceon_name][intr_dust_map_name] = self.model.sfr_intrinsic_dust_luminosity_map_faceon

        # Dust bolometric luminosity
        if self.has_sfr_faceon_map(dust_map_name): self.sfr_maps[faceon_name][dust_map_name] = self.load_sfr_faceon_map(dust_map_name)
        elif self.model.has_sfr_dust_luminosity_map_faceon: self.sfr_maps[faceon_name][dust_map_name] = self.model.sfr_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def sfr_edgeon_maps(self):

        """
        This function ...
        :return:
        """

        return self.sfr_maps[edgeon_name]

    # -----------------------------------------------------------------

    def get_sfr_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_sfr_edgeon_map(bol_map_name): self.sfr_maps[edgeon_name][bol_map_name] = self.load_sfr_edgeon_map(bol_map_name)
        elif self.model.has_sfr_bolometric_luminosity_map_edgeon: self.sfr_maps[edgeon_name][bol_map_name] = self.model.sfr_bolometric_luminosity_map_edgeon

        # Direct
        if self.has_sfr_edgeon_map(direct_map_name): self.sfr_maps[edgeon_name][direct_map_name] = self.load_sfr_edgeon_map(direct_map_name)
        elif self.model.has_sfr_direct_stellar_luminosity_map_edgeon: self.sfr_maps[edgeon_name][direct_map_name] = self.model.sfr_direct_stellar_luminosity_map_edgeon

        # (observed) FUV
        if self.has_sfr_edgeon_map(fuv_map_name): self.sfr_maps[edgeon_name][fuv_map_name] = self.load_sfr_edgeon_map(fuv_map_name)
        elif self.model.has_sfr_fuv_luminosity_map_edgeon: self.sfr_maps[edgeon_name][fuv_map_name] = self.model.sfr_fuv_luminosity_map_edgeon

        # Intrinsic FUV
        if self.has_sfr_edgeon_map(intr_fuv_map_name): self.sfr_maps[edgeon_name][intr_fuv_map_name] = self.load_sfr_edgeon_map(intr_fuv_map_name)
        elif self.model.has_sfr_intrinsic_fuv_luminosity_map_edgeon: self.sfr_maps[edgeon_name][intr_fuv_map_name] = self.model.sfr_intrinsic_fuv_luminosity_map_edgeon

        # SFR
        if self.has_sfr_edgeon_map(sfr_map_name): self.sfr_maps[edgeon_name][sfr_map_name] = self.load_sfr_edgeon_map(sfr_map_name)
        elif self.model.has_star_formation_rate_map_edgeon: self.sfr_maps[edgeon_name][sfr_map_name] = self.model.star_formation_rate_map_edgeon

        # Dust mass
        if self.has_sfr_edgeon_map(dust_mass_map_name): self.sfr_maps[edgeon_name][dust_mass_map_name] = self.load_sfr_edgeon_map(dust_mass_map_name)
        elif self.model.has_sfr_dust_mass_map_edgeon: self.sfr_maps[edgeon_name][dust_mass_map_name] = self.model.sfr_dust_mass_map_edgeon

        # Stellar bolometric luminosity
        if self.has_sfr_edgeon_map(stellar_lum_map_name): self.sfr_maps[edgeon_name][stellar_lum_map_name] = self.load_sfr_edgeon_map(stellar_lum_map_name)
        elif self.model.has_sfr_stellar_luminosity_map_edgeon: self.sfr_maps[edgeon_name][stellar_lum_map_name] = self.model.sfr_stellar_luminosity_map_edgeon

        # Intrinsic dust luminosity
        if self.has_sfr_edgeon_map(intr_dust_map_name): self.sfr_maps[edgeon_name][intr_dust_map_name] = self.load_sfr_edgeon_map(intr_dust_map_name)
        elif self.model.has_sfr_intrinsic_dust_luminosity_map_edgeon: self.sfr_maps[edgeon_name][intr_dust_map_name] = self.model.sfr_intrinsic_dust_luminosity_map_edgeon

        # Dust bolometric luminosity
        if self.has_sfr_edgeon_map(dust_map_name): self.sfr_maps[edgeon_name][dust_map_name] = self.load_sfr_edgeon_map(dust_map_name)
        elif self.model.has_sfr_dust_luminosity_map_edgeon: self.sfr_maps[edgeon_name][dust_map_name] = self.model.sfr_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------

    def get_unevolved_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.get_unevolved_maps_earth()

        # Face-on
        if self.do_faceon: self.get_unevolved_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.get_unevolved_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def unevolved_earth_maps(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_maps[earth_name]

    # -----------------------------------------------------------------

    def get_unevolved_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_unevolved_earth_map(bol_map_name): self.unevolved_maps[earth_name][bol_map_name] = self.load_unevolved_earth_map(bol_map_name)
        elif self.model.has_unevolved_bolometric_luminosity_map_earth: self.unevolved_maps[earth_name][bol_map_name] = self.model.unevolved_bolometric_luminosity_map_earth

        # Direct
        if self.has_unevolved_earth_map(direct_map_name): self.unevolved_maps[earth_name][direct_map_name] = self.load_unevolved_earth_map(direct_map_name)
        elif self.model.has_unevolved_direct_stellar_luminosity_map_earth: self.unevolved_maps[earth_name][direct_map_name] = self.model.unevolved_direct_stellar_luminosity_map_earth

        # FUV
        if self.has_unevolved_earth_map(fuv_map_name): self.unevolved_maps[earth_name][fuv_map_name] = self.load_unevolved_earth_map(direct_map_name)
        elif self.model.has_unevolved_fuv_luminosity_map_earth: self.unevolved_maps[earth_name][fuv_map_name] = self.model.unevolved_fuv_luminosity_map_earth

        # Intrinsic FUV
        if self.has_unevolved_earth_map(intr_fuv_map_name): self.unevolved_maps[earth_name][intr_fuv_map_name] = self.load_unevolved_earth_map(intr_fuv_map_name)
        elif self.model.has_unevolved_intrinsic_fuv_luminosity_map_earth: self.unevolved_maps[earth_name][intr_fuv_map_name] = self.model.unevolved_intrinsic_fuv_luminosity_map_earth

        # SFR
        if self.has_unevolved_earth_map(sfr_map_name): self.unevolved_maps[earth_name][sfr_map_name] = self.load_unevolved_earth_map(sfr_map_name)
        elif self.model.has_unevolved_star_formation_rate_map_earth: self.unevolved_maps[earth_name][sfr_map_name] = self.model.unevolved_star_formation_rate_map_earth

        # Dust luminosity
        if self.has_unevolved_earth_map(dust_map_name): self.unevolved_maps[earth_name][dust_map_name] = self.load_unevolved_earth_map(dust_map_name)
        elif self.model.has_unevolved_dust_luminosity_map_earth: self.unevolved_maps[earth_name][dust_map_name] = self.model.unevolved_dust_luminosity_map_earth

    # -----------------------------------------------------------------

    @property
    def unevolved_faceon_maps(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_maps[faceon_name]

    # -----------------------------------------------------------------

    def get_unevolved_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_unevolved_faceon_map(bol_map_name): self.unevolved_maps[faceon_name][bol_map_name] = self.load_unevolved_faceon_map(bol_map_name)
        elif self.model.has_unevolved_bolometric_luminosity_map_faceon: self.unevolved_maps[faceon_name][bol_map_name] = self.model.unevolved_bolometric_luminosity_map_faceon

        # Direct
        if self.has_unevolved_faceon_map(direct_map_name): self.unevolved_maps[faceon_name][direct_map_name] = self.load_unevolved_faceon_map(direct_map_name)
        elif self.model.has_unevolved_direct_stellar_luminosity_map_faceon: self.unevolved_maps[faceon_name][direct_map_name] = self.model.unevolved_direct_stellar_luminosity_map_faceon

        # FUV
        if self.has_unevolved_faceon_map(fuv_map_name): self.unevolved_maps[faceon_name][fuv_map_name] = self.load_unevolved_faceon_map(direct_map_name)
        elif self.model.has_unevolved_fuv_luminosity_map_faceon: self.unevolved_maps[faceon_name][fuv_map_name] = self.model.unevolved_fuv_luminosity_map_faceon

        # Intrinsic FUV
        if self.has_unevolved_faceon_map(intr_fuv_map_name): self.unevolved_maps[faceon_name][intr_fuv_map_name] = self.load_unevolved_faceon_map(intr_fuv_map_name)
        elif self.model.has_unevolved_intrinsic_fuv_luminosity_map_faceon: self.unevolved_maps[faceon_name][intr_fuv_map_name] = self.model.unevolved_intrinsic_fuv_luminosity_map_faceon

        # SFR
        if self.has_unevolved_faceon_map(sfr_map_name): self.unevolved_maps[faceon_name][sfr_map_name] = self.load_unevolved_faceon_map(sfr_map_name)
        elif self.model.has_unevolved_star_formation_rate_map_faceon: self.unevolved_maps[faceon_name][sfr_map_name] = self.model.unevolved_star_formation_rate_map_faceon

        # Dust luminosity
        if self.has_unevolved_faceon_map(dust_map_name): self.unevolved_maps[faceon_name][dust_map_name] = self.load_unevolved_faceon_map(dust_map_name)
        elif self.model.has_unevolved_dust_luminosity_map_faceon: self.unevolved_maps[faceon_name][dust_map_name] = self.model.unevolved_dust_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def unevolved_edgeon_maps(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_maps[edgeon_name]

    # -----------------------------------------------------------------

    def get_unevolved_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Bolometric
        if self.has_unevolved_edgeon_map(bol_map_name): self.unevolved_maps[edgeon_name][bol_map_name] = self.load_unevolved_edgeon_map(bol_map_name)
        elif self.model.has_unevolved_bolometric_luminosity_map_edgeon: self.unevolved_maps[edgeon_name][bol_map_name] = self.model.unevolved_bolometric_luminosity_map_edgeon

        # Direct
        if self.has_unevolved_edgeon_map(direct_map_name): self.unevolved_maps[edgeon_name][direct_map_name] = self.load_unevolved_edgeon_map(direct_map_name)
        elif self.model.has_unevolved_direct_stellar_luminosity_map_edgeon: self.unevolved_maps[edgeon_name][direct_map_name] = self.model.unevolved_direct_stellar_luminosity_map_edgeon

        # FUV
        if self.has_unevolved_edgeon_map(fuv_map_name): self.unevolved_maps[edgeon_name][fuv_map_name] = self.load_unevolved_edgeon_map(direct_map_name)
        elif self.model.has_unevolved_fuv_luminosity_map_edgeon: self.unevolved_maps[edgeon_name][fuv_map_name] = self.model.unevolved_fuv_luminosity_map_edgeon

        # Intrinsic FUV
        if self.has_unevolved_edgeon_map(intr_fuv_map_name): self.unevolved_maps[edgeon_name][intr_fuv_map_name] = self.load_unevolved_edgeon_map(intr_fuv_map_name)
        elif self.model.has_unevolved_intrinsic_fuv_luminosity_map_edgeon: self.unevolved_maps[edgeon_name][intr_fuv_map_name] = self.model.unevolved_intrinsic_fuv_luminosity_map_edgeon

        # SFR
        if self.has_unevolved_edgeon_map(sfr_map_name): self.unevolved_maps[edgeon_name][sfr_map_name] = self.load_unevolved_edgeon_map(sfr_map_name)
        elif self.model.has_unevolved_star_formation_rate_map_edgeon: self.unevolved_maps[edgeon_name][sfr_map_name] = self.model.unevolved_star_formation_rate_map_edgeon

        # Dust luminosity
        if self.has_unevolved_edgeon_map(dust_map_name): self.unevolved_maps[edgeon_name][dust_map_name] = self.load_unevolved_edgeon_map(dust_map_name)
        elif self.model.has_unevolved_dust_luminosity_map_edgeon: self.unevolved_maps[edgeon_name][dust_map_name] = self.model.unevolved_dust_luminosity_map_edgeon

    # -----------------------------------------------------------------

    def get_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.get_dust_maps_earth()

        # Face-on
        if self.do_faceon: self.get_dust_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.get_dust_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def dust_earth_maps(self):

        """
        This function ...
        :return:
        """

        return self.dust_maps[earth_name]

    # -----------------------------------------------------------------

    def get_dust_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Dust mass
        if self.has_dust_earth_map(diffuse_mass_map_name): self.dust_maps[earth_name][diffuse_mass_map_name] = self.load_dust_earth_map(diffuse_mass_map_name)
        elif self.model.has_diffuse_dust_mass_map_earth: self.dust_maps[earth_name][diffuse_mass_map_name] = self.model.diffuse_dust_mass_map_earth

        # Total dust mass
        if self.has_dust_earth_map(mass_map_name): self.dust_maps[earth_name][mass_map_name] = self.load_dust_earth_map(mass_map_name)
        elif self.model.has_dust_mass_map_earth: self.dust_maps[earth_name][mass_map_name] = self.model.dust_mass_map_earth

    # -----------------------------------------------------------------

    @property
    def dust_faceon_maps(self):

        """
        This function ...
        :return:
        """

        return self.dust_maps[faceon_name]

    # -----------------------------------------------------------------

    def get_dust_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Dust mass
        if self.has_dust_faceon_map(diffuse_mass_map_name): self.dust_maps[faceon_name][diffuse_mass_map_name] = self.load_dust_faceon_map(diffuse_mass_map_name)
        elif self.model.has_diffuse_dust_mass_map_faceon: self.dust_maps[faceon_name][diffuse_mass_map_name] = self.model.diffuse_dust_mass_map_faceon

        # Total dust mass
        if self.has_dust_faceon_map(mass_map_name): self.dust_maps[faceon_name][mass_map_name] = self.load_dust_faceon_map(mass_map_name)
        elif self.model.has_dust_mass_map_faceon: self.dust_maps[faceon_name][mass_map_name] = self.model.dust_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def dust_edgeon_maps(self):

        """
        This function ...
        :return:
        """

        return self.dust_maps[edgeon_name]

    # -----------------------------------------------------------------

    def get_dust_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Dust mass
        if self.has_dust_edgeon_map(diffuse_mass_map_name): self.dust_maps[edgeon_name][diffuse_mass_map_name] = self.load_dust_edgeon_map(diffuse_mass_map_name)
        elif self.model.has_diffuse_dust_mass_map_edgeon: self.dust_maps[edgeon_name][diffuse_mass_map_name] = self.model.diffuse_dust_mass_map_edgeon

        # Total dust mass
        if self.has_dust_edgeon_map(mass_map_name): self.dust_maps[edgeon_name][mass_map_name] = self.load_dust_edgeon_map(mass_map_name)
        elif self.model.has_dust_mass_map_edgeon: self.dust_maps[edgeon_name][mass_map_name] = self.model.dust_mass_map_edgeon

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the parameters
        self.write_parameters()

        # Write the maps
        self.write_maps()

    # -----------------------------------------------------------------

    @lazyproperty
    def parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_path, "parameters")

    # -----------------------------------------------------------------

    def write_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the parameters ...")

        # Intrinsic parameters
        self.write_intrinsic_parameters()

        # Derived parameter values of total model
        self.write_total_parameters()

        # Derived parameter values of bulge
        self.write_bulge_parameters()

        # Derived parameter values of disk
        self.write_disk_parameters()

        # Derived parameter values of old stellar component
        self.write_old_parameters()

        # Derived parameter values of young stellar component
        self.write_young_parameters()

        # Derived parameter values of SFR component
        self.write_sfr_parameters()

        # Derived parameter values of unevolved components
        self.write_unevolved_parameters()

        # Derived parameter values of dust component
        self.write_dust_parameters()

    # -----------------------------------------------------------------

    @property
    def intrinsic_parameters(self):

        """
        This function ...
        :return:
        """

        return self.model.parameter_values

    # -----------------------------------------------------------------

    @property
    def intrinsic_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.parameters_path, "intrinsic.dat")

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_parameters(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.intrinsic_parameters_path)

    # -----------------------------------------------------------------

    def write_intrinsic_parameters(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Writing the intrinsic model parameters ...")

        # Write
        write_dict(self.intrinsic_parameters, self.intrinsic_parameters_path)

    # -----------------------------------------------------------------

    @property
    def total_parameters(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_total

    # -----------------------------------------------------------------

    @property
    def total_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.parameters_path, "total.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_parameters(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.total_parameters_path)

    # -----------------------------------------------------------------

    def write_total_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the derived parameters of the total model ...")

        # Write
        write_dict(self.total_parameters, self.total_parameters_path)

    # -----------------------------------------------------------------

    @property
    def bulge_parameters(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_bulge

    # -----------------------------------------------------------------

    @property
    def bulge_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.parameters_path, "bulge.dat")

    # -----------------------------------------------------------------

    @property
    def has_bulge_parameters(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.bulge_parameters_path)

    # -----------------------------------------------------------------

    def write_bulge_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the derived parameters of the old stellar bulge ...")

        # Write
        write_dict(self.bulge_parameters, self.bulge_parameters_path)

    # -----------------------------------------------------------------

    @property
    def disk_parameters(self):

        """
        Thisn function ...
        :return:
        """

        return self.model.derived_parameter_values_disk

    # -----------------------------------------------------------------

    @property
    def disk_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.parameters_path, "disk.dat")

    # -----------------------------------------------------------------

    @property
    def has_disk_parameters(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.disk_parameters_path)

    # -----------------------------------------------------------------

    def write_disk_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the derived parameters of the old stellar disk ...")

        # Write
        write_dict(self.disk_parameters, self.disk_parameters_path)

    # -----------------------------------------------------------------

    @property
    def old_parameters(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_old

    # -----------------------------------------------------------------

    @property
    def old_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.parameters_path, "old.dat")

    # -----------------------------------------------------------------

    @property
    def has_old_parameters(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_parameters_path)

    # -----------------------------------------------------------------

    def write_old_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the derived parameters of the old stellar components ...")

        # Write
        write_dict(self.old_parameters, self.old_parameters_path)

    # -----------------------------------------------------------------

    @property
    def young_parameters(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_young

    # -----------------------------------------------------------------

    @property
    def young_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.parameters_path, "young.dat")

    # -----------------------------------------------------------------

    @property
    def has_young_parameters(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_parameters_path)

    # -----------------------------------------------------------------

    def write_young_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the derived parameters of the young stellar component ...")

        # Write
        write_dict(self.young_parameters, self.young_parameters_path)

    # -----------------------------------------------------------------

    @property
    def sfr_parameters(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_sfr

    # -----------------------------------------------------------------

    @property
    def sfr_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.parameters_path, "sfr.dat")

    # -----------------------------------------------------------------

    @property
    def has_sfr_parameters(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_parameters_path)

    # -----------------------------------------------------------------

    def write_sfr_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the derived parameters of the SFR stellar component ...")

        # Write
        write_dict(self.sfr_parameters, self.sfr_parameters_path)

    # -----------------------------------------------------------------

    @property
    def unevolved_parameters(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_unevolved

    # -----------------------------------------------------------------

    @property
    def unevolved_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.parameters_path, "unevolved.dat")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_parameters(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.unevolved_parameters_path)

    # -----------------------------------------------------------------

    def write_unevolved_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the derived parameters of the unevolved stellar components ...")

        # Write
        write_dict(self.unevolved_parameters, self.unevolved_parameters_path)

    # -----------------------------------------------------------------

    @property
    def dust_parameters(self):

        """
        This function ...
        :return:
        """

        return self.model.derived_parameter_values_dust

    # -----------------------------------------------------------------

    @property
    def dust_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.parameters_path, "dust.dat")

    # -----------------------------------------------------------------

    @property
    def has_dust_parameters(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.dust_parameters_path)

    # -----------------------------------------------------------------

    def write_dust_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the derived parameters of the dust component ...")

        # Write
        write_dict(self.dust_parameters, self.dust_parameters_path)

    # -----------------------------------------------------------------

    @property
    def maps_total_earth_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_total_earth_path

    # -----------------------------------------------------------------

    @property
    def maps_total_faceon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_total_faceon_path

    # -----------------------------------------------------------------

    @property
    def maps_total_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_total_edgeon_path

    # -----------------------------------------------------------------

    @property
    def maps_bulge_earth_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_bulge_earth_path

    # -----------------------------------------------------------------

    @property
    def maps_bulge_faceon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_bulge_faceon_path

    # -----------------------------------------------------------------

    @property
    def maps_bulge_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_bulge_edgeon_path

    # -----------------------------------------------------------------

    @property
    def maps_disk_earth_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_disk_earth_path

    # -----------------------------------------------------------------

    @property
    def maps_disk_faceon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_disk_faceon_path

    # -----------------------------------------------------------------

    @property
    def maps_disk_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_disk_edgeon_path

    # -----------------------------------------------------------------

    @property
    def maps_old_earth_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_old_earth_path

    # -----------------------------------------------------------------

    @property
    def maps_old_faceon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_old_faceon_path

    # -----------------------------------------------------------------

    @property
    def maps_old_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_old_edgeon_path

    # -----------------------------------------------------------------

    @property
    def maps_young_earth_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_young_earth_path

    # -----------------------------------------------------------------

    @property
    def maps_young_faceon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_young_faceon_path

    # -----------------------------------------------------------------

    @property
    def maps_young_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_young_edgeon_path

    # -----------------------------------------------------------------

    @property
    def maps_sfr_earth_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_sfr_earth_path

    # -----------------------------------------------------------------

    @property
    def maps_sfr_faceon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_sfr_faceon_path

    # -----------------------------------------------------------------

    @property
    def maps_sfr_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_sfr_edgeon_path

    # -----------------------------------------------------------------

    @property
    def maps_unevolved_earth_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_unevolved_earth_path

    # -----------------------------------------------------------------

    @property
    def maps_unevolved_faceon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_unevolved_faceon_path

    # -----------------------------------------------------------------

    @property
    def maps_unevolved_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_unevolved_edgeon_path

    # -----------------------------------------------------------------

    @property
    def maps_dust_earth_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_dust_earth_path

    # -----------------------------------------------------------------

    @property
    def maps_dust_faceon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_dust_faceon_path

    # -----------------------------------------------------------------

    @property
    def maps_dust_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return self.properties_maps_dust_edgeon_path

    # -----------------------------------------------------------------

    @property
    def do_earth(self):

        """
        This function ...
        :return:
        """

        return self.config.earth

    # -----------------------------------------------------------------

    @property
    def do_faceon(self):

        """
        This function ...
        :return:
        """

        return self.config.faceon

    # -----------------------------------------------------------------

    @property
    def do_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.config.edgeon

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Total
        self.write_total_maps()

        # Bulge
        self.write_bulge_maps()

        # Disk
        self.write_disk_maps()

        # Old
        self.write_old_maps()

        # Young
        self.write_young_maps()

        # SFR
        self.write_sfr_maps()

        # Unevolved
        self.write_unevolved_maps()

        # Dust
        self.write_dust_maps()

    # -----------------------------------------------------------------

    def write_total_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the total maps ...")

        # Earth
        if self.do_earth: self.write_total_maps_earth()

        # Face-on
        if self.do_faceon: self.write_total_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.write_total_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def total_earth_map_names(self):

        """
        This function ...
        :return:
        """

        return self.total_maps[earth_name].keys()

    # -----------------------------------------------------------------

    def get_total_earth_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_total_earth_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_total_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_total_earth_map_path(name))

    # -----------------------------------------------------------------

    def load_total_earth_map(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_total_earth_map_path(name))

    # -----------------------------------------------------------------

    def remove_total_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_total_earth_map_path(name))

    # -----------------------------------------------------------------

    def write_total_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the total maps in the earth projection ...")

        # Loop over the map names
        for name in self.total_earth_map_names:

            # Has to be written?
            if self.has_total_earth_map(name): continue

            # Determine the path
            path = self.get_total_earth_map_path(name)

            # Save
            self.total_maps[earth_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def total_faceon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.total_maps[faceon_name].keys()

    # -----------------------------------------------------------------

    def get_total_faceon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_total_faceon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_total_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_total_faceon_map_path(name))

    # -----------------------------------------------------------------

    def load_total_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_total_faceon_map_path(name))

    # -----------------------------------------------------------------

    def remove_total_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_total_faceon_map_path(name))

    # -----------------------------------------------------------------

    def write_total_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the total maps in the faceon projection ...")

        # Loop over the map names
        for name in self.total_faceon_map_names:

            # Has to be written?
            if self.has_total_faceon_map(name): continue

            # Determine the path
            path = self.get_total_faceon_map_path(name)

            # Save
            self.total_maps[faceon_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def total_edgeon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.total_maps[edgeon_name].keys()

    # -----------------------------------------------------------------

    def get_total_edgeon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_total_edgeon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_total_edgeon_map(self, name):

        """
        Thisnf unction ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_total_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def load_total_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_total_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def remove_total_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_total_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def write_total_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the total maps in the edgeon projection ...")

        # Loop over the map names
        for name in self.total_edgeon_map_names:

            # Has to be written?
            if self.has_total_edgeon_map(name): continue

            # Determine the path
            path = self.get_total_edgeon_map_path(name)

            # Save
            self.total_maps[edgeon_name][name].saveto(path)

    # -----------------------------------------------------------------

    def write_bulge_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the bulge maps ...")

        # Earth
        if self.do_earth: self.write_bulge_maps_earth()

        # Face-on
        if self.do_faceon: self.write_bulge_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.write_bulge_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def bulge_earth_map_names(self):

        """
        This function ...
        :return:
        """

        return self.bulge_maps[earth_name].keys()

    # -----------------------------------------------------------------

    def get_bulge_earth_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_bulge_earth_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_bulge_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_bulge_earth_map_path(name))

    # -----------------------------------------------------------------

    def load_bulge_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_bulge_earth_map_path(name))

    # -----------------------------------------------------------------

    def remove_bulge_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_bulge_earth_map_path(name))

    # -----------------------------------------------------------------

    def write_bulge_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the bulge maps in the earth projection ...")

        # Loop over the map names
        for name in self.bulge_earth_map_names:

            # Write?
            if self.has_bulge_earth_map(name): continue

            # Get path
            path = self.get_bulge_earth_map_path(name)

            # Save
            self.bulge_maps[earth_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def bulge_faceon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.bulge_maps[faceon_name].keys()

    # -----------------------------------------------------------------

    def get_bulge_faceon_map_path(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_bulge_faceon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_bulge_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_bulge_faceon_map_path(name))

    # -----------------------------------------------------------------

    def load_bulge_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_bulge_faceon_map_path(name))

    # -----------------------------------------------------------------

    def remove_bulge_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_bulge_faceon_map_path(name))

    # -----------------------------------------------------------------

    def write_bulge_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the bulge maps in the face-on projection ...")

        # Loop over the map names
        for name in self.bulge_faceon_map_names:

            # Write
            if self.has_bulge_faceon_map(name): continue

            # Get path
            path = self.get_bulge_faceon_map_path(name)

            # Save
            self.bulge_maps[faceon_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def bulge_edgeon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.bulge_maps[edgeon_name].keys()

    # -----------------------------------------------------------------

    def get_bulge_edgeon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_bulge_edgeon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_bulge_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_bulge_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def load_bulge_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_bulge_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def remove_bulge_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_bulge_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def write_bulge_maps_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Writing the bulge maps in the edge-on projection ...")

        # Loop over the the map names
        for name in self.bulge_edgeon_map_names:

            # Write?
            if self.has_bulge_edgeon_map(name): continue

            # Get path
            path = self.get_bulge_edgeon_map_path(name)

            # Save
            self.bulge_maps[edgeon_name][name].saveto(path)

    # -----------------------------------------------------------------

    def write_disk_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the disk maps ...")

        # Earth
        if self.do_earth: self.write_disk_maps_earth()

        # Face-on
        if self.do_faceon: self.write_disk_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.write_disk_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def disk_earth_map_names(self):

        """
        Thisn function ...
        :return:
        """

        return self.disk_maps[earth_name].keys()

    # -----------------------------------------------------------------

    def get_disk_earth_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_disk_earth_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_disk_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_disk_earth_map_path(name))

    # -----------------------------------------------------------------

    def load_disk_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_disk_earth_map_path(name))

    # -----------------------------------------------------------------

    def remove_disk_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_disk_earth_map_path(name))

    # -----------------------------------------------------------------

    def write_disk_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the disk maps in the earth projection ...")

        # Loop over the map names
        for name in self.disk_earth_map_names:

            # Write?
            if self.has_disk_earth_map(name): continue

            # Get path
            path = self.get_disk_earth_map_path(name)

            # Save
            self.disk_maps[earth_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def disk_faceon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.disk_maps[faceon_name].keys()

    # -----------------------------------------------------------------

    def get_disk_faceon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_disk_faceon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_disk_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_disk_faceon_map_path(name))

    # -----------------------------------------------------------------

    def load_disk_faceon_map(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_disk_faceon_map_path(name))

    # -----------------------------------------------------------------

    def remove_disk_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_disk_faceon_map_path(name))

    # -----------------------------------------------------------------

    def write_disk_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the disk maps in the faceon projection ...")

        # Loop over the map names
        for name in self.disk_faceon_map_names:

            # Write?
            if self.has_disk_faceon_map(name): continue

            # Get path
            path = self.get_disk_faceon_map_path(name)

            # Save
            self.disk_maps[faceon_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def disk_edgeon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.disk_maps[edgeon_name].keys()

    # -----------------------------------------------------------------

    def get_disk_edgeon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_disk_edgeon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_disk_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_disk_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def load_disk_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_disk_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def remove_disk_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_disk_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def write_disk_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the disk maps in the edgeon projection ...")

        # Loop over the map names
        for name in self.disk_edgeon_map_names:

            # Write?
            if self.has_disk_edgeon_map(name): continue

            # Get path
            path = self.get_disk_edgeon_map_path(name)

            # Save
            self.disk_maps[edgeon_name][name].saveto(path)

    # -----------------------------------------------------------------

    def write_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the old maps ...")

        # Earth
        if self.do_earth: self.write_old_maps_earth()

        # Face-on
        if self.do_faceon: self.write_old_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.write_old_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def old_earth_map_names(self):

        """
        This function ...
        :return:
        """

        return self.old_maps[earth_name].keys()

    # -----------------------------------------------------------------

    def get_old_earth_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_old_earth_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_old_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_old_earth_map_path(name))

    # -----------------------------------------------------------------

    def load_old_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_old_earth_map_path(name))

    # -----------------------------------------------------------------

    def remove_old_earth_map(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_old_earth_map_path(name))

    # -----------------------------------------------------------------

    def write_old_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the old maps in the earth projection ...")

        # Loop over the map names
        for name in self.old_earth_map_names:

            # Write
            if self.has_old_earth_map(name): continue

            # Get path
            path = self.get_old_earth_map_path(name)

            # Save
            self.old_maps[earth_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def old_faceon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.old_maps[faceon_name].keys()

    # -----------------------------------------------------------------

    def get_old_faceon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_old_faceon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_old_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_old_faceon_map_path(name))

    # -----------------------------------------------------------------

    def load_old_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_old_faceon_map_path(name))

    # -----------------------------------------------------------------

    def remove_old_faceon_map(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_old_faceon_map_path(name))

    # -----------------------------------------------------------------

    def write_old_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the old maps in the face-on projection ...")

        # Loop over the map names
        for name in self.old_faceon_map_names:

            # Write
            if self.has_old_faceon_map(name): continue

            # Get path
            path = self.get_old_faceon_map_path(name)

            # Save
            self.old_maps[faceon_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def old_edgeon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.old_maps[edgeon_name].keys()

    # -----------------------------------------------------------------

    def get_old_edgeon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_old_edgeon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_old_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_old_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def load_old_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_old_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def remove_old_edgeon_map(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_old_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def write_old_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the old maps in the edge-on projection ...")

        # Loop over the map names
        for name in self.old_edgeon_map_names:

            # Write
            if self.has_old_edgeon_map(name): continue

            # Get path
            path = self.get_old_edgeon_map_path(name)

            # Save
            self.old_maps[edgeon_name][name].saveto(path)

    # -----------------------------------------------------------------

    def write_young_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth: self.write_young_maps_earth()

        # Face-on
        if self.do_faceon: self.write_young_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.write_young_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def young_earth_map_names(self):

        """
        This function ...
        :return:
        """

        return self.young_maps[earth_name].keys()

    # -----------------------------------------------------------------

    def get_young_earth_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_young_earth_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_young_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_young_earth_map_path(name))

    # -----------------------------------------------------------------

    def load_young_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_young_earth_map_path(name))

    # -----------------------------------------------------------------

    def remove_young_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_young_earth_map_path(name))

    # -----------------------------------------------------------------

    def write_young_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the young maps in the earth projection ...")

        # Loop over the map names
        for name in self.young_earth_map_names:

            # Write?
            if self.has_young_earth_map(name): continue

            # Get the path
            path = self.get_young_earth_map_path(name)

            # Save
            self.young_maps[earth_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def young_faceon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.young_maps[faceon_name].keys()

    # -----------------------------------------------------------------

    def get_young_faceon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_young_faceon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_young_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_young_faceon_map_path(name))

    # -----------------------------------------------------------------

    def load_young_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_young_faceon_map_path(name))

    # -----------------------------------------------------------------

    def remove_young_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_young_faceon_map_path(name))

    # -----------------------------------------------------------------

    def write_young_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the young maps in the faceon projection ...")

        # Loop over the map names
        for name in self.young_faceon_map_names:

            # Write?
            if self.has_young_faceon_map(name): continue

            # Get the path
            path = self.get_young_faceon_map_path(name)

            # Save
            self.young_maps[faceon_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def young_edgeon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.young_maps[edgeon_name].keys()

    # -----------------------------------------------------------------

    def get_young_edgeon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_young_edgeon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_young_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_young_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def load_young_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_young_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def remove_young_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_young_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def write_young_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the young maps in the edgeon projection ...")

        # Loop over the map names
        for name in self.young_edgeon_map_names:

            # Write?
            if self.has_young_edgeon_map(name): continue

            # Get the path
            path = self.get_young_edgeon_map_path(name)

            # Save
            self.young_maps[edgeon_name][name].saveto(path)

    # -----------------------------------------------------------------

    def write_sfr_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SFR maps ...")

        # Earth
        if self.do_earth: self.write_sfr_maps_earth()

        # Face-on
        if self.do_faceon: self.write_sfr_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.write_sfr_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def sfr_earth_map_names(self):

        """
        This function ...
        :return:
        """

        return self.sfr_maps[earth_name].keys()

    # -----------------------------------------------------------------

    def get_sfr_earth_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_sfr_earth_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_sfr_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_sfr_earth_map_path(name))

    # -----------------------------------------------------------------

    def load_sfr_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_sfr_earth_map_path(name))

    # -----------------------------------------------------------------

    def remove_sfr_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_sfr_earth_map_path(name))

    # -----------------------------------------------------------------

    def write_sfr_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SFR maps in the earth projection ...")

        # Loop over the map names
        for name in self.sfr_earth_map_names:

            # Write?
            if self.has_sfr_earth_map(name): continue

            # Get path
            path = self.get_sfr_earth_map_path(name)

            # Save
            self.sfr_maps[earth_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def sfr_faceon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.sfr_maps[faceon_name].keys()

    # -----------------------------------------------------------------

    def get_sfr_faceon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_sfr_faceon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_sfr_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_sfr_faceon_map_path(name))

    # -----------------------------------------------------------------

    def load_sfr_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_sfr_faceon_map_path(name))

    # -----------------------------------------------------------------

    def remove_sfr_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_sfr_faceon_map_path(name))

    # -----------------------------------------------------------------

    def write_sfr_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SFR maps in the faceon projection ...")

        # Loop over the map names
        for name in self.sfr_faceon_map_names:

            # Write?
            if self.has_sfr_faceon_map(name): continue

            # Get path
            path = self.get_sfr_faceon_map_path(name)

            # Save
            self.sfr_maps[faceon_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def sfr_edgeon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.sfr_maps[edgeon_name].keys()

    # -----------------------------------------------------------------

    def get_sfr_edgeon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_sfr_edgeon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_sfr_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_sfr_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def load_sfr_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_sfr_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def remove_sfr_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_sfr_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def write_sfr_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SFR maps in the edgeon projection ...")

        # Loop over the map names
        for name in self.sfr_edgeon_map_names:

            # Write?
            if self.has_sfr_edgeon_map(name): continue

            # Get path
            path = self.get_sfr_edgeon_map_path(name)

            # Save
            self.sfr_maps[edgeon_name][name].saveto(path)

    # -----------------------------------------------------------------

    def write_unevolved_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the unevolved maps ...")

        # Earth
        if self.do_earth: self.write_unevolved_maps_earth()

        # Face-on
        if self.do_faceon: self.write_unevolved_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.write_unevolved_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def unevolved_earth_map_names(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_maps[earth_name].keys()

    # -----------------------------------------------------------------

    def get_unevolved_earth_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_unevolved_earth_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_unevolved_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_unevolved_earth_map_path(name))

    # -----------------------------------------------------------------

    def load_unevolved_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_unevolved_earth_map_path(name))

    # -----------------------------------------------------------------

    def remove_unevolved_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_unevolved_earth_map_path(name))

    # -----------------------------------------------------------------

    def write_unevolved_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the unevolved maps in the earth projection ...")

        # Loop over the map names
        for name in self.unevolved_earth_map_names:

            # Write?
            if self.has_unevolved_earth_map(name): continue

            # Get path
            path = self.get_unevolved_earth_map_path(name)

            # Save
            self.unevolved_maps[earth_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def unevolved_faceon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_maps[faceon_name].keys()

    # -----------------------------------------------------------------

    def get_unevolved_faceon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_unevolved_faceon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_unevolved_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_unevolved_faceon_map_path(name))

    # -----------------------------------------------------------------

    def load_unevolved_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_unevolved_faceon_map_path(name))

    # -----------------------------------------------------------------

    def remove_unevolved_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_unevolved_faceon_map_path(name))

    # -----------------------------------------------------------------

    def write_unevolved_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the unevolved maps in the faceon projection ...")

        # Loop over the map names
        for name in self.unevolved_faceon_map_names:

            # Write?
            if self.has_unevolved_faceon_map(name): continue

            # Get path
            path = self.get_unevolved_faceon_map_path(name)

            # Save
            self.unevolved_maps[faceon_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def unevolved_edgeon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_maps[edgeon_name].keys()

    # -----------------------------------------------------------------

    def get_unevolved_edgeon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_unevolved_edgeon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_unevolved_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_unevolved_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def load_unevolved_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_unevolved_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def remove_unevolved_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_unevolved_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def write_unevolved_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the unevolved maps in the edgeon projection ...")

        # Loop over the map names
        for name in self.unevolved_edgeon_map_names:

            # Write?
            if self.has_unevolved_edgeon_map(name): continue

            # Get path
            path = self.get_unevolved_edgeon_map_path(name)

            # Save
            self.unevolved_maps[edgeon_name][name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps ...")

        # Earth
        if self.do_earth: self.write_dust_maps_earth()

        # Face-on
        if self.do_faceon: self.write_dust_maps_faceon()

        # Edge-on
        if self.do_edgeon: self.write_dust_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def dust_earth_map_names(self):

        """
        This function ...
        :return:
        """

        return self.dust_maps[earth_name].keys()

    # -----------------------------------------------------------------

    def get_dust_earth_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_dust_earth_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_dust_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_dust_earth_map_path(name))

    # -----------------------------------------------------------------

    def load_dust_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_dust_earth_map_path(name))

    # -----------------------------------------------------------------

    def remove_dust_earth_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_dust_earth_map_path(name))

    # -----------------------------------------------------------------

    def write_dust_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps in the earth projection ...")

        # Loop over the map names
        for name in self.dust_earth_map_names:

            # Write?
            if self.has_dust_earth_map(name): continue

            # Get path
            path = self.get_dust_earth_map_path(name)

            # Save
            self.dust_maps[earth_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def dust_faceon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.dust_maps[faceon_name].keys()

    # -----------------------------------------------------------------

    def get_dust_faceon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_dust_faceon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_dust_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_dust_faceon_map_path(name))

    # -----------------------------------------------------------------

    def load_dust_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_dust_faceon_map_path(name))

    # -----------------------------------------------------------------

    def remove_dust_faceon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_dust_faceon_map_path(name))

    # -----------------------------------------------------------------

    def write_dust_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps in the faceon projection ...")

        # Loop over the map names
        for name in self.dust_faceon_map_names:

            # Write?
            if self.has_dust_faceon_map(name): continue

            # Get path
            path = self.get_dust_faceon_map_path(name)

            # Save
            self.dust_maps[faceon_name][name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def dust_edgeon_map_names(self):

        """
        This function ...
        :return:
        """

        return self.dust_maps[edgeon_name].keys()

    # -----------------------------------------------------------------

    def get_dust_edgeon_map_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_dust_edgeon_path, name + ".fits")

    # -----------------------------------------------------------------

    def has_dust_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_dust_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def load_dust_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return Frame.from_file(self.get_dust_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def remove_dust_edgeon_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fs.remove_file(self.get_dust_edgeon_map_path(name))

    # -----------------------------------------------------------------

    def write_dust_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps in the edgeon projection ...")

        # Loop over the map names
        for name in self.dust_edgeon_map_names:

            # Write?
            if self.has_dust_edgeon_map(name): continue

            # Get path
            path = self.get_dust_edgeon_map_path(name)

            # Save
            self.dust_maps[edgeon_name][name].saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Maps
        self.plot_maps()

    # -----------------------------------------------------------------

    def plot_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps ...")

        # Total
        self.plot_total_maps()

        # Bulge
        self.plot_bulge_maps()

        # Disk
        self.plot_disk_maps()

        # Old
        self.plot_old_maps()

        # Young
        self.plot_young_maps()

        # SFR
        self.plot_sfr_maps()

        # Unevolved
        self.plot_unevolved_maps()

        # Dust
        self.plot_dust_maps()

    # -----------------------------------------------------------------

    def plot_total_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the total model ...")

        # Earth
        if self.do_earth: self.plot_total_maps_earth()

        # Faceon
        if self.do_faceon: self.plot_total_maps_faceon()

        # Edgeon
        if self.do_edgeon: self.plot_total_maps_edgeon()

    # -----------------------------------------------------------------

    def get_total_earth_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_total_earth_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_total_earth_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_total_earth_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_total_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.total_earth_map_names:

            # Needs plotting?
            if self.has_total_earth_map_plot(name): continue

            # Get path
            path = self.get_total_earth_map_plot_path(name)

            # Plot
            self.plot_map(self.total_earth_maps[name], path)

    # -----------------------------------------------------------------

    def get_total_faceon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_total_faceon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_total_faceon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_total_faceon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_total_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.total_faceon_map_names:

            # Needs plotting?
            if self.has_total_faceon_map_plot(name): continue

            # Get path
            path = self.get_total_faceon_map_plot_path(name)

            # Plot
            self.plot_map(self.total_faceon_maps[name], path)

    # -----------------------------------------------------------------

    def get_total_edgeon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_total_edgeon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_total_edgeon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_total_edgeon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_total_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.total_edgeon_map_names:

            # Needs plotting?
            if self.has_total_edgeon_map_plot(name): continue

            # Get path
            path = self.get_total_edgeon_map_plot_path(name)

            # Plot
            self.plot_map(self.total_edgeon_maps[name], path)

    # -----------------------------------------------------------------

    def plot_bulge_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the old stellar bulge component ...")

        # Earth
        if self.do_earth: self.plot_bulge_maps_earth()

        # Faceon
        if self.do_faceon: self.plot_bulge_maps_faceon()

        # Edgeon
        if self.do_edgeon: self.plot_bulge_maps_edgeon()

    # -----------------------------------------------------------------

    def get_bulge_earth_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_bulge_earth_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_bulge_earth_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_bulge_earth_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_bulge_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.bulge_earth_map_names:

            # Needs plotting?
            if self.has_bulge_earth_map_plot(name): continue

            # Get path
            path = self.get_bulge_earth_map_plot_path(name)

            # Plot
            self.plot_map(self.bulge_earth_maps[name], path)

    # -----------------------------------------------------------------

    def get_bulge_faceon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_bulge_faceon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_bulge_faceon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_bulge_faceon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_bulge_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.bulge_faceon_map_names:

            # Needs plotting?
            if self.has_bulge_faceon_map_plot(name): continue

            # Get path
            path = self.get_bulge_faceon_map_plot_path(name)

            # Plot
            self.plot_map(self.bulge_faceon_maps[name], path)

    # -----------------------------------------------------------------

    def get_bulge_edgeon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_bulge_edgeon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_bulge_edgeon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_bulge_edgeon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_bulge_maps_edgeon(self):

        """
        Thisfunction ...
        :return:
        """

        # Loop over the map names
        for name in self.bulge_edgeon_map_names:

            # Needs plotting?
            if self.has_bulge_edgeon_map_plot(name): continue

            # Get path
            path = self.get_bulge_edgeon_map_plot_path(name)

            # Plot
            self.plot_map(self.bulge_edgeon_maps[name], path)

    # -----------------------------------------------------------------

    def plot_disk_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the old stellar disk component ...")

        # Earth
        if self.do_earth: self.plot_disk_maps_earth()

        # Faceon
        if self.do_faceon: self.plot_disk_maps_faceon()

        # Edgeon
        if self.do_edgeon: self.plot_disk_maps_edgeon()

    # -----------------------------------------------------------------

    def get_disk_earth_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_disk_earth_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_disk_earth_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_disk_earth_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_disk_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.disk_earth_map_names:

            # Needs plotting?
            if self.has_disk_earth_map_plot(name): continue

            # Get path
            path = self.get_disk_earth_map_plot_path(name)

            # Plot
            self.plot_map(self.disk_earth_maps[name], path)

    # -----------------------------------------------------------------

    def get_disk_faceon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_disk_faceon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_disk_faceon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_disk_faceon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_disk_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.disk_faceon_map_names:

            # Needs plotting?
            if self.has_disk_faceon_map_plot(name): continue

            # Get path
            path = self.get_disk_faceon_map_plot_path(name)

            # Plot
            self.plot_map(self.disk_faceon_maps[name], path)

    # -----------------------------------------------------------------

    def get_disk_edgeon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_disk_edgeon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_disk_edgeon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_disk_edgeon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_disk_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.disk_edgeon_map_names:

            # Needs plotting?
            if self.has_disk_edgeon_map_plot(name): continue

            # Get path
            path = self.get_disk_edgeon_map_plot_path(name)

            # Plot
            self.plot_map(self.disk_edgeon_maps[name], path)

    # -----------------------------------------------------------------

    def plot_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the old stellar component ...")

        # Earth
        if self.do_earth: self.plot_old_maps_earth()

        # Faceon
        if self.do_faceon: self.plot_old_maps_faceon()

        # Edgeon
        if self.do_edgeon: self.plot_old_maps_edgeon()

    # -----------------------------------------------------------------

    def get_old_earth_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_old_earth_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_old_earth_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_old_earth_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_old_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.old_earth_map_names:

            # Needs plotting?
            if self.has_old_earth_map_plot(name): continue

            # Get path
            path = self.get_old_earth_map_plot_path(name)

            # Plot
            self.plot_map(self.old_earth_maps[name], path)

    # -----------------------------------------------------------------

    def get_old_faceon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_old_faceon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_old_faceon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_old_faceon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_old_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.old_faceon_map_names:

            # Needs plotting?
            if self.has_old_faceon_map_plot(name): continue

            # Get path
            path = self.get_old_faceon_map_plot_path(name)

            # Plot
            self.plot_map(self.old_faceon_maps[name], path)

    # -----------------------------------------------------------------

    def get_old_edgeon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_old_edgeon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_old_edgeon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_old_edgeon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_old_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.old_edgeon_map_names:

            # Needs plotting?
            if self.has_old_edgeon_map_plot(name): continue

            # Get path
            path = self.get_old_edgeon_map_plot_path(name)

            # Plot
            self.plot_map(self.old_edgeon_maps[name], path)

    # -----------------------------------------------------------------

    def plot_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the young stellar component ...")

        # Earth
        if self.do_earth: self.plot_young_maps_earth()

        # Faceon
        if self.do_faceon: self.plot_young_maps_faceon()

        # Edgeon
        if self.do_edgeon: self.plot_young_maps_edgeon()

    # -----------------------------------------------------------------

    def get_young_earth_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_young_earth_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_young_earth_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_young_earth_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_young_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.young_earth_map_names:

            # Needs plotting?
            if self.has_young_earth_map_plot(name): continue

            # Get path
            path = self.get_young_earth_map_plot_path(name)

            # Plot
            self.plot_map(self.young_earth_maps[name], path)

    # -----------------------------------------------------------------

    def get_young_faceon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_young_faceon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_young_faceon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_young_faceon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_young_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.young_faceon_map_names:

            # Needs plotting?
            if self.has_young_faceon_map_plot(name): continue

            # Get path
            path = self.get_young_faceon_map_plot_path(name)

            # Plot
            self.plot_map(self.young_faceon_maps[name], path)

    # -----------------------------------------------------------------

    def get_young_edgeon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_young_edgeon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_young_edgeon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_young_edgeon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_young_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Loop over the map names
        for name in self.young_edgeon_map_names:

            # Needs plotting?
            if self.has_young_edgeon_map_plot(name): continue

            # Get path
            path = self.get_young_edgeon_map_plot_path(name)

            # Plot
            self.plot_map(self.young_edgeon_maps[name], path)

    # -----------------------------------------------------------------

    def plot_sfr_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the SFR component ...")

        # Earth
        if self.do_earth: self.plot_sfr_maps_earth()

        # Faceon
        if self.do_faceon: self.plot_sfr_maps_faceon()

        # Edgeon
        if self.do_edgeon: self.plot_sfr_maps_edgeon()

    # -----------------------------------------------------------------

    def get_sfr_earth_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_sfr_earth_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_sfr_earth_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_sfr_earth_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_sfr_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.sfr_earth_map_names:

            # Needs plotting?
            if self.has_sfr_earth_map_plot(name): continue

            # Get path
            path = self.get_sfr_earth_map_plot_path(name)

            # Plot
            self.plot_map(self.sfr_earth_maps[name], path)

    # -----------------------------------------------------------------

    def get_sfr_faceon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_sfr_faceon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_sfr_faceon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_sfr_faceon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_sfr_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.sfr_faceon_map_names:

            # Needs plotting?
            if self.has_sfr_faceon_map_plot(name): continue

            # Get path
            path = self.get_sfr_faceon_map_plot_path(name)

            # Plot
            self.plot_map(self.sfr_faceon_maps[name], path)

    # -----------------------------------------------------------------

    def get_sfr_edgeon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_sfr_edgeon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_sfr_edgeon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_sfr_edgeon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_sfr_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.sfr_edgeon_map_names:

            # Needs plotting?
            if self.has_sfr_edgeon_map_plot(name): continue

            # Get path
            path = self.get_sfr_edgeon_map_plot_path(name)

            # Plot
            self.plot_map(self.sfr_edgeon_maps[name], path)

    # -----------------------------------------------------------------

    def plot_unevolved_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the unevolved stellar component ...")

        # Earth
        if self.do_earth: self.plot_unevolved_maps_earth()

        # Faceon
        if self.do_faceon: self.plot_unevolved_maps_faceon()

        # Edgeon
        if self.do_edgeon: self.plot_unevolved_maps_edgeon()

    # -----------------------------------------------------------------

    def get_unevolved_earth_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_unevolved_earth_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_unevolved_earth_map_plot(self, name):

        """
        This function ..
        :param name:
        :return:
        """

        return fs.is_file(self.get_unevolved_earth_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_unevolved_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.unevolved_earth_map_names:

            # Needs plotting?
            if self.has_unevolved_earth_map_plot(name): continue

            # Get path
            path = self.get_unevolved_earth_map_plot_path(name)

            # Plot
            self.plot_map(self.unevolved_earth_maps[name], path)

    # -----------------------------------------------------------------

    def get_unevolved_faceon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_unevolved_faceon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_unevolved_faceon_map_plot(self, name):

        """
        This function ..
        :param name:
        :return:
        """

        return fs.is_file(self.get_unevolved_faceon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_unevolved_maps_faceon(self):

        """
        Thisn function ...
        :return:
        """

        # Loop over the maps
        for name in self.unevolved_faceon_map_names:

            # Needs plotting?
            if self.has_unevolved_faceon_map_plot(name): continue

            # Get path
            path = self.get_unevolved_faceon_map_plot_path(name)

            # Plot
            self.plot_map(self.unevolved_faceon_maps[name], path)

    # -----------------------------------------------------------------

    def get_unevolved_edgeon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_unevolved_edgeon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_unevolved_edgeon_map_plot(self, name):

        """
        This function ..
        :param name:
        :return:
        """

        return fs.is_file(self.get_unevolved_edgeon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_unevolved_maps_edgeon(self):

        """
        This function ...
        :return:
        """
        
        # Loop over the maps
        for name in self.unevolved_edgeon_map_names:

            # Needs plotting?
            if self.has_unevolved_edgeon_map_plot(name): continue

            # Get path
            path = self.get_unevolved_edgeon_map_plot_path(name)

            # Plot
            self.plot_map(self.unevolved_edgeon_maps[name], path)

    # -----------------------------------------------------------------

    def plot_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the dust component ...")

        # Earth
        if self.do_earth: self.plot_dust_maps_earth()

        # Faceon
        if self.do_faceon: self.plot_dust_maps_faceon()

        # Edgeon
        if self.do_edgeon: self.plot_dust_maps_edgeon()

    # -----------------------------------------------------------------

    def get_dust_earth_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_dust_earth_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_dust_earth_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_dust_earth_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_dust_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_earth_map_names:

            # Needs plotting
            if self.has_dust_earth_map_plot(name): continue

            # Get path
            path = self.get_dust_earth_map_plot_path(name)

            # Plot
            self.plot_map(self.dust_earth_maps[name], path)

    # -----------------------------------------------------------------

    def get_dust_faceon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_dust_faceon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_dust_faceon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_dust_faceon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_dust_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_faceon_map_names:

            # Needs plotting
            if self.has_dust_faceon_map_plot(name): continue

            # Get path
            path = self.get_dust_faceon_map_plot_path(name)

            # Plot
            self.plot_map(self.dust_faceon_maps[name], path)

    # -----------------------------------------------------------------

    def get_dust_edgeon_map_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.maps_dust_edgeon_path, name + ".pdf")

    # -----------------------------------------------------------------

    def has_dust_edgeon_map_plot(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.is_file(self.get_dust_edgeon_map_plot_path(name))

    # -----------------------------------------------------------------

    def plot_dust_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_edgeon_map_names:

            # Needs plotting
            if self.has_dust_edgeon_map_plot(name): continue

            # Get path
            path = self.get_dust_edgeon_map_plot_path(name)

            # Plot
            self.plot_map(self.dust_edgeon_maps[name], path)

    # -----------------------------------------------------------------

    def plot_map(self, frame, path):

        """
        This function ...
        :param frame:
        :param path:
        :return:
        """

        plot_map(frame, path=path, cmap="inferno", colorbar=True)

# -----------------------------------------------------------------
