#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.component Contains the AnalysisComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from .context import AnalysisContext
from ...core.tools.utils import lazyproperty
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

bol_map_name = "bol"
intr_stellar_map_name = "intr_stellar" # intrinsic stellar (bol) luminosity (transparent)
obs_stellar_map_name = "obs_stellar" # observed stellar (bol) luminosity
diffuse_dust_map_name = "diffuse_dust"
dust_map_name = "dust" # dust (bol) luminosity
#dust_with_internal_map_name = "dust_with_internal" # dust (bol) luminosity + internal dust (MAPPINGS)
scattered_map_name = "scattered" # scattered stellar luminosity
absorbed_diffuse_map_name = "absorbed_diffuse"
#absorbed_map_name = "absorbed" # absorbed stellar luminosity
#absorbed_with_internal_map_name = "absorbed_with_internal" # absorbed stellar luminosity + internal absorption (MAPPINGS)
fabs_diffuse_map_name = "fabs_diffuse"
fabs_map_name = "fabs"
attenuated_map_name = "attenuated" # attenuated stellar luminosity
direct_map_name = "direct" # direct stellar luminosity
stellar_mass_map_name = "stellar_mass"  # stellar mass
ssfr_map_name = "ssfr" # specific star formation rate

i1_map_name = "i1"
intr_i1_map_name = "intr_i1"

fuv_map_name = "fuv"
intr_fuv_map_name = "intr_fuv"

sfr_map_name = "sfr"
dust_mass_map_name = "dust_mass"
stellar_lum_map_name = "stellar_lum"
#dust_lum_map_name = "dust_lum"
intr_dust_map_name = "intr_dust"
diffuse_mass_map_name = "diffuse_mass"
mass_map_name = "mass"
#total_mass_map_name = "total_mass"
lum_map_name = "lum"
total_lum_map_name = "total_lum"

# -----------------------------------------------------------------

earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

class AnalysisComponent(GalaxyModelingComponent):
    
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
        super(AnalysisComponent, self).__init__(*args, **kwargs)

        # The analysis context
        self.context = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisComponent, self).setup()

        # Create the analysis context
        self.context = AnalysisContext(self.analysis_path)

    # -----------------------------------------------------------------

    @property
    def timing_table_path(self):

        """
        This function ...
        :return:
        """

        return self.context.timing_table_path

    # -----------------------------------------------------------------

    @property
    def memory_table_path(self):

        """
        This function ...
        :return:
        """

        return self.context.memory_table_path

    # -----------------------------------------------------------------

    @property
    def timing_table(self):

        """
        This function ...
        :return:
        """

        return self.context.timing_table

    # -----------------------------------------------------------------

    @property
    def memory_table(self):

        """
        This function ...
        :return:
        """

        return self.context.memory_table

    # -----------------------------------------------------------------

    @property
    def cached_table(self):

        """
        This fucntion ...
        :return:
        """

        return self.context.cached_table

    # -----------------------------------------------------------------

    @property
    def analysis_run_names(self):

        """
        This function ...
        :return:
        """

        return self.context.analysis_run_names

    # -----------------------------------------------------------------

    def get_run_path(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        return self.context.get_run_path(run_name)

    # -----------------------------------------------------------------

    def get_run_info(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        return self.context.get_run_info(run_name)

    # -----------------------------------------------------------------

    def get_run(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        return self.context.get_run(run_name)

# -----------------------------------------------------------------

class AnalysisRunComponent(AnalysisComponent):

    """
    This class ...
    """

    @lazyproperty
    def analysis_run(self):

        """
        This function ...
        :return:
        """

        return self.get_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.model_name

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
    def properties_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.properties_path

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_path, "maps")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_total_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_path, "total")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_total_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_total_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_total_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_total_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_total_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_total_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_bulge_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_path, "bulge")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_bulge_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_bulge_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_bulge_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_bulge_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_bulge_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_bulge_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_disk_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_path, "disk")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_disk_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_disk_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_disk_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_disk_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_disk_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_disk_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_old_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_path, "old")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_old_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_old_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_old_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_old_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_old_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_old_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_young_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_path, "young")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_young_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_young_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_young_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_young_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_young_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_young_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_sfr_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_path, "sfr")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_sfr_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_sfr_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_sfr_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_sfr_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_sfr_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_sfr_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_unevolved_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_path, "unevolved")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_unevolved_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_unevolved_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_unevolved_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_unevolved_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_unevolved_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_unevolved_path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_dust_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_path, "dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_dust_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_dust_path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_dust_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_dust_path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def properties_maps_dust_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_maps_dust_path, edgeon_name)

    # -----------------------------------------------------------------

    def get_total_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric luminosity
        if which == bol_map_name:

            if orientation == earth_name: return self.model.total_bolometric_luminosity_map_earth
            elif orientation == faceon_name: return self.model.total_bolometric_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.total_bolometric_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Intrinsic stellar luminosity (transparent luminosity)
        if which == intr_stellar_map_name:

            if orientation == earth_name: return self.model.total_intrinsic_stellar_luminosity_map_earth
            elif orientation == faceon_name: return self.model.total_intrinsic_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.total_intrinsic_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Observed stellar luminosity
        elif which == obs_stellar_map_name:

            if orientation == earth_name: return self.model.total_observed_stellar_luminosity_map_earth
            elif orientation == faceon_name: return self.model.total_observed_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.total_observed_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Diffuse dust emission luminosity
        elif which == diffuse_dust_map_name:

            if orientation == earth_name: return self.model.total_diffuse_dust_luminosity_map_earth
            elif orientation == faceon_name: return self.model.total_diffuse_dust_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.total_diffuse_dust_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Dust emission luminosity
        elif which == dust_map_name:

            if orientation == earth_name: return self.model.total_dust_luminosity_map_earth
            elif orientation == faceon_name: return self.model.total_dust_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.total_dust_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Scattered stellar luminosity
        elif which == scattered_map_name:

            if orientation == earth_name: return self.model.total_scattered_stellar_luminosity_map_earth
            elif orientation == faceon_name: return self.model.total_scattered_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.total_scattered_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Absorbed stellar luminosity (by diffuse dust) (extinction)
        # absorbed = transparent - observed stellar (= observed - dust = direct + scattered)
        elif which == absorbed_diffuse_map_name:

            if orientation == earth_name: return self.model.total_absorbed_diffuse_stellar_luminosity_map_earth
            elif orientation == faceon_name: return self.model.total_absorbed_diffuse_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.total_absorbed_diffuse_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Absorbed stellar luminosity (extinction)
        # CUBE INFORMATION IS NOT AVAILABLE, SO MAP IS NOT USEFUL (IS JUST THE SAME AS DUST EMISSION MAP)
        #elif which == absorbed_map_name:

        # Fraction of energy absorbed by DIFFUSE dust
        elif which == fabs_diffuse_map_name:

            if orientation == earth_name: return self.model.total_fabs_diffuse_map_earth
            elif orientation == faceon_name: return self.model.total_fabs_diffuse_map_faceon
            elif orientation == edgeon_name: return self.model.total_fabs_diffuse_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Fraction of energy absorbed by dust
        elif which == fabs_map_name:

            if orientation == earth_name: return self.model.total_fabs_map_earth
            elif orientation == faceon_name: return self.model.total_fabs_map_faceon
            elif orientation == edgeon_name: return self.model.total_fabs_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Attenuated stellar luminosity (attenuation)
        elif which == attenuated_map_name: # attenuated = transparent - direct stellar

            if orientation == earth_name: return self.model.total_attenuated_stellar_luminosity_map_earth
            elif orientation == faceon_name: return self.model.total_attenuated_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.total_attenuated_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Direct luminosity
        elif which == direct_map_name:

            if orientation == earth_name: return self.model.total_direct_stellar_luminosity_map_earth
            elif orientation == faceon_name: return self.model.total_direct_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.total_direct_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Star formation rate
        elif which == sfr_map_name:

            if orientation == earth_name: return self.model.total_star_formation_rate_map_earth
            elif orientation == faceon_name: return self.model.total_star_formation_rate_map_faceon
            elif orientation == edgeon_name: return self.model.total_star_formation_rate_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Stellar mass
        elif which == stellar_mass_map_name:

            if orientation == earth_name: return self.model.total_stellar_mass_map_earth
            elif orientation == faceon_name: return self.model.total_stellar_mass_map_faceon
            elif orientation == edgeon_name: return self.model.total_stellar_mass_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Specific star formation rate
        elif which == ssfr_map_name:

            if orientation == earth_name: return self.model.total_ssfr_map_earth
            elif orientation == faceon_name: return self.model.total_ssfr_map_faceon
            elif orientation == edgeon_name: return self.model.total_ssfr_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------

    def get_bulge_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric luminosity
        if which == bol_map_name:

            if orientation == earth_name: return self.model.old_bulge_bolometric_luminosity_map
            elif orientation == faceon_name: return self.model.old_bulge_bolometric_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_bulge_bolometric_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Direct
        elif which == direct_map_name:

            if orientation == earth_name: return self.model.old_bulge_direct_stellar_luminosity_map
            elif orientation == faceon_name: return self.model.old_bulge_direct_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_bulge_direct_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation +"'")

        # (observed) I1 lum
        elif which == i1_map_name:

            if orientation == earth_name: return self.model.old_bulge_i1_luminosity_map
            elif orientation == faceon_name: return self.model.old_bulge_i1_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_bulge_i1_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Intrinsic I1
        elif which == intr_i1_map_name:

            if orientation == earth_name: return self.model.old_bulge_intrinsic_i1_luminosity_map
            elif orientation == faceon_name: return self.model.old_bulge_intrinsic_i1_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_bulge_intrinsic_i1_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Dust luminosity
        elif which == dust_map_name:

            if orientation == earth_name: return self.model.old_bulge_dust_luminosity_map
            elif orientation == faceon_name: return self.model.old_bulge_dust_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_bulge_dust_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------

    def get_disk_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric
        if which == bol_map_name:

            if orientation == earth_name: return self.model.old_disk_bolometric_luminosity_map
            elif orientation == faceon_name: return self.model.old_disk_bolometric_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_disk_bolometric_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Direct
        if which == direct_map_name:

            if orientation == earth_name: return self.model.old_disk_direct_stellar_luminosity_map
            elif orientation == faceon_name: return self.model.old_disk_direct_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_disk_direct_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # (observed) I1
        elif which == i1_map_name:

            if orientation == earth_name: return self.model.old_disk_i1_luminosity_map
            elif orientation == faceon_name: return self.model.old_disk_i1_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_disk_i1_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Intrinsic I1
        elif which == intr_i1_map_name:

            if orientation == earth_name: return self.model.old_disk_intrinsic_i1_luminosity_map
            elif orientation == faceon_name: return self.model.old_disk_intrinsic_i1_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_disk_intrinsic_i1_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Dust luminosity
        elif which == dust_map_name:

            if orientation == earth_name: return self.model.old_disk_dust_luminosity_map
            elif orientation == faceon_name: return self.model.old_disk_dust_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_disk_dust_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------

    def get_old_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric
        if which == bol_map_name:

            if orientation == earth_name: return self.model.old_bolometric_luminosity_map
            elif orientation == faceon_name: return self.model.old_bolometric_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_bolometric_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Direct
        if which == direct_map_name:

            if orientation == earth_name: return self.model.old_direct_stellar_luminosity_map
            elif orientation == faceon_name: return self.model.old_direct_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_direct_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # (observed) I1
        elif which == i1_map_name:

            if orientation == earth_name: return self.model.old_i1_luminosity_map
            elif orientation == faceon_name: return self.model.old_i1_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_i1_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Intrinsic I1
        elif which == intr_i1_map_name:

            if orientation == earth_name: return self.model.old_intrinsic_i1_luminosity_map
            elif orientation == faceon_name: return self.model.old_intrinsic_i1_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_intrinsic_i1_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Dust luminosity
        elif which == dust_map_name:

            if orientation == earth_name: return self.model.old_dust_luminosity_map
            elif orientation == faceon_name: return self.model.old_bulge_dust_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.old_bulge_dust_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    def get_young_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric
        if which == bol_map_name:

            if orientation == earth_name: return self.model.young_bolometric_luminosity_map
            elif orientation == faceon_name: return self.model.young_bolometric_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.young_bolometric_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Direct
        elif which == direct_map_name:

            if orientation == earth_name: return self.model.young_direct_stellar_luminosity_map
            elif orientation == faceon_name: return self.model.young_direct_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.young_direct_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # (observed) FUV
        elif which == fuv_map_name:

            if orientation == earth_name: return self.model.young_fuv_luminosity_map
            elif orientation == faceon_name: return self.model.young_fuv_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.young_fuv_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Intrinsic FUV
        elif which == intr_fuv_map_name:

            if orientation == earth_name: return self.model.young_intrinsic_fuv_luminosity_map
            elif orientation == faceon_name: return self.model.young_intrinsic_fuv_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.young_intrinsic_fuv_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Dust luminosity
        elif which == dust_map_name:

            if orientation == earth_name: return self.model.young_dust_luminosity_map
            elif orientation == faceon_name: return self.model.young_dust_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.young_dust_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------

    def get_sfr_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric
        if which == bol_map_name:

            if orientation == earth_name: return self.model.sfr_bolometric_luminosity_map
            elif orientation == faceon_name: return self.model.sfr_bolometric_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.sfr_bolometric_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Direct
        elif which == direct_map_name:

            if orientation == earth_name: return self.model.sfr_direct_stellar_luminosity_map
            elif orientation == faceon_name: return self.model.sfr_direct_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.sfr_direct_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # (observed) FUV
        elif which == fuv_map_name:

            if orientation == earth_name: return self.model.sfr_fuv_luminosity_map
            elif orientation == faceon_name: return self.model.sfr_fuv_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.sfr_fuv_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Intrinsic FUV
        elif which == intr_fuv_map_name:

            if orientation == earth_name: return self.model.sfr_intrinsic_fuv_luminosity_map
            elif orientation == faceon_name: return self.model.sfr_intrinsic_fuv_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.sfr_intrinsic_fuv_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # SFR
        elif which == sfr_map_name:

            if orientation == earth_name: return self.model.star_formation_rate_map
            elif orientation == faceon_name: return self.model.star_formation_rate_map_faceon
            elif orientation == edgeon_name: return self.model.star_formation_rate_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Dust mass
        elif which == dust_mass_map_name:

            if orientation == earth_name: return self.model.sfr_dust_mass_map
            elif orientation == faceon_name: return self.model.sfr_dust_mass_map_faceon
            elif orientation == edgeon_name: return self.model.sfr_dust_mass_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Stellar bolometric luminosity
        elif which == stellar_lum_map_name:

            if orientation == earth_name: return self.model.sfr_stellar_luminosity_map
            elif orientation == faceon_name: return self.model.sfr_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.sfr_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Intrinsic dust luminosity
        elif which == intr_dust_map_name:

            if orientation == earth_name: return self.model.sfr_intrinsic_dust_luminosity_map
            elif orientation == faceon_name: return self.model.sfr_intrinsic_dust_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.sfr_intrinsic_dust_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Dust bolometric luminosity
        elif which == dust_map_name:

            if orientation == earth_name: return self.model.sfr_dust_luminosity_map
            elif orientation == faceon_name: return self.model.sfr_dust_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.sfr_dust_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------

    def get_unevolved_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Bolometric
        if which == bol_map_name:

            if orientation == earth_name: return self.model.unevolved_bolometric_luminosity_map
            elif orientation == faceon_name: return self.model.unevolved_bolometric_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.unevolved_bolometric_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Direct
        elif which == direct_map_name:

            if orientation == earth_name: return self.model.unevolved_direct_stellar_luminosity_map
            elif orientation == faceon_name: return self.model.unevolved_direct_stellar_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.unevolved_direct_stellar_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # FUV
        elif which == fuv_map_name:

            if orientation == earth_name: return self.model.unevolved_fuv_luminosity_map
            elif orientation == faceon_name: return self.model.unevolved_fuv_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.unevolved_fuv_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Intrinsic FUV
        elif which == intr_fuv_map_name:

            if orientation == earth_name: return self.model.unevolved_intrinsic_fuv_luminosity_map
            elif orientation == faceon_name: return self.model.unevolved_intrinsic_fuv_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.unevolved_intrinsic_fuv_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # SFR
        elif which == sfr_map_name:

            if orientation == earth_name: return self.model.unevolved_star_formation_rate_map
            elif orientation == faceon_name: return self.model.unevolved_star_formation_rate_map_faceon
            elif orientation == edgeon_name: return self.model.unevolved_star_formation_rate_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Dust luminosity
        elif which == dust_map_name:

            if orientation == earth_name: return self.model.unevolved_dust_luminosity_map
            elif orientation == faceon_name: return self.model.unevolved_dust_luminosity_map_faceon
            elif orientation == edgeon_name: return self.model.unevolved_dust_luminosity_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

    # -----------------------------------------------------------------

    def get_dust_map(self, which, orientation=earth_name):

        """
        This function ...
        :param which:
        :param orientation:
        :return:
        """

        # Dust mass
        if which == diffuse_mass_map_name:

            if orientation == earth_name: return self.model.diffuse_dust_mass_map
            elif orientation == faceon_name: return self.model.diffuse_dust_mass_map_faceon
            elif orientation == edgeon_name: return self.model.diffuse_dust_mass_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Total dust mass
        elif which == mass_map_name:

            if orientation == earth_name: return self.model.dust_mass_map
            elif orientation == faceon_name: return self.model.dust_mass_map_faceon
            elif orientation == edgeon_name: return self.model.dust_mass_map_edgeon
            else: raise ValueError("Invalid orientation: '" + orientation + "'")

        # Invalid
        else: raise ValueError("Invalid argument: '" + which + "'")

# -----------------------------------------------------------------
