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

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty
from ...core.tools.serialization import write_dict

# -----------------------------------------------------------------

class PropertiesAnalyser(AnalysisComponent):
    
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

        # The analysis run
        self.analysis_run = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup of the base class
        super(PropertiesAnalyser, self).setup(**kwargs)

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
    def properties_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.properties_path

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

    @lazyproperty
    def maps_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.properties_path, "maps")

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.maps_path, "earth")

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.maps_path, "faceon")

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.maps_path, "edgeon")

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
        self.write_total_maps_earth()

        # Face-on
        self.write_total_maps_faceon()

        # Edge-on
        self.write_total_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_total_bol_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_total_bolometric_luminosity_map and not self.has_total_bol_map

    # -----------------------------------------------------------------

    @property
    def total_bol_map(self):

        """
        This function ...
        :return:
        """

        return self.model.total_bolometric_luminosity_map

    # -----------------------------------------------------------------

    @property
    def total_bol_map_path(self):

        """
        Thisn function ...
        :return:
        """

        return fs.join(self.maps_earth_path, "total_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_bol_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.total_bol_map_path)

    # -----------------------------------------------------------------

    def write_total_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the total maps in the earth projection ...")

        # Bolometric
        if self.do_write_total_bol_map: self.total_bol_map.saveto(self.total_bol_map_path)

    # -----------------------------------------------------------------

    def write_total_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the total maps in the faceon projection ...")

        # Bolometric
        if self.has_total_bol_map_faceon: self.total_bol_map_faceon.saveto(self.total_bol_map_faceon_path)

    # -----------------------------------------------------------------

    def write_total_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the total maps in the edgeon projection ...")

        # Bolometric
        if self.has_total_bol_map_edgeon: self.total_bol_map_edgeon.saveto(self.total_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_bulge_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the bulge maps ...")

        # Earth
        self.write_bulge_maps_earth()

        # Face-on
        self.write_bulge_maps_faceon()

        # Edge-on
        self.write_bulge_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_bulge_i1_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_bulge_i1_luminosity_map and not self.has_bulge_i1_map

    # -----------------------------------------------------------------

    @property
    def bulge_i1_map(self):

        """
        This function ...
        :return:
        """

        return self.model.bulge_i1_luminosity_map

    # -----------------------------------------------------------------

    @property
    def bulge_i1_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_earth_path, "bulge_i1.fits")

    # -----------------------------------------------------------------

    @property
    def has_bulge_i1_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.bulge_i1_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_bulge_bol_map(self):

        """
        Thisf unction ...
        :return:
        """

        return self.model.has_bulge_bolometric_luminosity_map and not self.has_bulge_bol_map

    # -----------------------------------------------------------------

    @property
    def bulge_bol_map(self):

        """
        This function ...
        :return:
        """

        return self.model.bulge_bolometric_luminosity_map

    # -----------------------------------------------------------------

    @property
    def bulge_bol_map_path(self):

        """
        Thisn function ...
        :return:
        """

        return fs.join(self.maps_earth_path, "bulge_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_bulge_bol_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.bulge_bol_map_path)

    # -----------------------------------------------------------------

    def write_bulge_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the bulge maps in the earth projection ...")

        # I1
        if self.do_write_bulge_i1_map: self.bulge_i1_map.saveto(self.bulge_i1_map_path)

        # Bol
        if self.do_write_bulge_bol_map: self.bulge_bol_map.saveto(self.bulge_bol_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_bulge_i1_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_bulge_i1_luminosity_map_faceon and not self.has_bulge_i1_map_faceon

    # -----------------------------------------------------------------

    @property
    def bulge_i1_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.bulge_i1_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def bulge_i1_map_faceon_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "bulge_i1.fits")

    # -----------------------------------------------------------------

    @property
    def has_bulge_i1_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.bulge_i1_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_bulge_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_bulge_bolometric_luminosity_map_faceon and not self.has_bulge_bol_map_faceon

    # -----------------------------------------------------------------

    @property
    def bulge_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.bulge_bolometric_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def bulge_bol_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "bulge_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_bulge_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.bulge_bol_map_faceon_path)

    # -----------------------------------------------------------------

    def write_bulge_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the bulge maps in the face-on projection ...")

        # I1
        if self.do_write_bulge_i1_map_faceon: self.bulge_i1_map_faceon.saveto(self.bulge_i1_map_faceon_path)

        # Bol
        if self.do_write_bulge_bol_map_faceon: self.bulge_bol_map_faceon.saveto(self.bulge_bol_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_bulge_i1_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_bulge_i1_luminosity_map_edgeon and not self.has_bulge_i1_map_edgeon

    # -----------------------------------------------------------------

    @property
    def bulge_i1_map_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        return self.model.bulge_i1_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def bulge_i1_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "bulge_i1.fits")

    # -----------------------------------------------------------------

    @property
    def has_bulge_i1_map_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        return fs.is_file(self.bulge_i1_map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_bulge_bol_map_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        return self.model.has_bulge_bolometric_luminosity_map_edgeon and not self.has_bulge_bol_map_edgeon

    # -----------------------------------------------------------------

    @property
    def bulge_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.bulge_bolometric_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def bulge_bol_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "bulge_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_bulge_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.bulge_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_bulge_maps_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Writing the bulge maps in the edge-on projection ...")

        # I1
        if self.do_write_bulge_i1_map_edgeon: self.bulge_i1_map_edgeon.saveto(self.bulge_i1_map_edgeon_path)

        # Bol
        if self.do_write_bulge_bol_map_edgeon: self.bulge_bol_map_edgeon.saveto(self.bulge_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def has_disk_i1_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.disk_i1_map_path)

    # -----------------------------------------------------------------

    @property
    def disk_i1_map(self):

        """
        This function ...
        :return:
        """

        return self.model.old_disk_i1_luminosity_map

    # -----------------------------------------------------------------

    @property
    def disk_i1_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_earth_path, "disk_i1.fits")

    # -----------------------------------------------------------------

    @property
    def has_disk_bol_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.disk_bol_map_path)

    # -----------------------------------------------------------------

    @property
    def disk_bol_map(self):

        """
        This function ...
        :return:
        """

        return self.model.old_disk_bolometric_luminosity_map

    # -----------------------------------------------------------------

    @property
    def disk_bol_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "disk_bol.fits")

    # -----------------------------------------------------------------

    def write_disk_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the disk maps ...")

        # Earth
        self.write_disk_maps_earth()

        # Face-on
        self.write_disk_maps_faceon()

        # Edge-on
        self.write_disk_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_disk_i1_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_old_disk_i1_luminosity_map and not self.has_disk_i1_map

    # -----------------------------------------------------------------

    @property
    def do_write_disk_bol_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_old_disk_bolometric_luminosity_map and not self.has_disk_bol_map

    # -----------------------------------------------------------------

    def write_disk_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the disk maps in the earth projection ...")

        # I1 lum
        if self.do_write_disk_i1_map: self.disk_i1_map.saveto(self.disk_i1_map_path)

        # Bol lum
        if self.do_write_disk_bol_map: self.disk_bol_map.saveto(self.disk_bol_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_disk_i1_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_old_disk_i1_luminosity_map and not self.has_disk_i1_map_faceon

    # -----------------------------------------------------------------

    @property
    def disk_i1_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.old_disk_i1_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def disk_i1_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "disk_i1.fits")

    # -----------------------------------------------------------------

    @property
    def has_disk_i1_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.disk_i1_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_disk_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_old_disk_bolometric_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def disk_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.disk_bolometric_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def disk_bol_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "disk_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_disk_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.disk_bol_map_faceon_path)

    # -----------------------------------------------------------------

    def write_disk_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the disk maps in the faceon projection ...")

        # I1 lum
        if self.do_write_disk_i1_map_faceon: self.disk_i1_map_faceon.saveto(self.disk_i1_map_faceon_path)

        # Bol lum
        if self.do_write_disk_bol_map_faceon: self.disk_bol_map_faceon.saveto(self.disk_bol_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_disk_i1_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_disk_i1_luminosity_map_edgeon and not self.has_disk_i1_map_edgeon

    # -----------------------------------------------------------------

    @property
    def disk_i1_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.disk_i1_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def disk_i1_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "disk_i1.fits")

    # -----------------------------------------------------------------

    @property
    def has_disk_i1_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.disk_i1_map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_disk_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_disk_bol_luminosity_map_edgeon and not self.has_disk_bol_map_edgeon

    # -----------------------------------------------------------------

    @property
    def disk_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.disk_bol_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def disk_bol_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "disk_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_disk_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.disk_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_disk_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the disk maps in the edgeon projection ...")

        # I1 lum
        if self.do_write_disk_i1_map_edgeon: self.disk_i1_map_edgeon.saveto(self.disk_i1_map_edgeon_path)

        # Bol lum
        if self.do_write_disk_bol_map_edgeon: self.disk_bol_map_edgeon.saveto(self.disk_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the old maps ...")

        # Earth
        self.write_old_maps_earth()

        # Face-on
        self.write_old_maps_faceon()

        # Edge-on
        self.write_old_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_old_i1_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_old_i1_luminosity_map and not self.has_old_i1_map

    # -----------------------------------------------------------------

    @property
    def old_i1_map(self):

        """
        This function ...
        :return:
        """

        return self.model.old_i1_luminosity_map

    # -----------------------------------------------------------------

    @property
    def old_i1_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_earth_path, "old_i1.fits")

    # -----------------------------------------------------------------

    @property
    def has_old_i1_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_i1_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_old_bol_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_old_bolometric_luminosity_map

    # -----------------------------------------------------------------

    @property
    def old_bol_map(self):

        """
        Thisf unction ...
        :return:
        """

        return self.model.old_bolometric_luminosity_map

    # -----------------------------------------------------------------

    @property
    def old_bol_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_earth_path, "old_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_old_bol_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_bol_map_path)

    # -----------------------------------------------------------------

    def write_old_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the old maps in the earth projection ...")

        # I1
        if self.do_write_old_i1_map: self.old_i1_map.saveto(self.old_i1_map_path)

        # Bol
        if self.do_write_old_bol_map: self.old_bol_map.saveto(self.old_bol_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_old_i1_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_old_i1_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def old_i1_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.old_i1_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def old_i1_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "old_i1.fits")

    # -----------------------------------------------------------------

    @property
    def has_old_i1_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_i1_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_old_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_old_bolometric_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def old_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.old_bolometric_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def old_bol_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "old_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_old_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_bol_map_faceon_path)

    # -----------------------------------------------------------------

    def write_old_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the old maps in the face-on projection ...")

        # I1
        if self.do_write_old_i1_map_faceon: self.old_i1_map_faceon.saveto(self.old_i1_map_faceon_path)

        # Bol
        if self.do_write_old_bol_map_faceon: self.old_bol_map_faceon.saveto(self.old_bol_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_old_i1_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_old_i1_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_i1_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.old_i1_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_i1_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "old_i1.fits")

    # -----------------------------------------------------------------

    @property
    def has_old_i1_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_i1_map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_old_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_old_bolometric_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.old_bolometric_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def old_bol_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "old_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_old_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_old_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the old maps in the edge-on projection ...")

        # I1
        if self.do_write_old_i1_map_edgeon: self.old_i1_map_edgeon.saveto(self.old_i1_map_edgeon_path)

        # Bol
        if self.do_write_old_bol_map_edgeon: self.old_bol_map_edgeon.saveto(self.old_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the young maps ...")

        # Earth
        self.write_young_maps_earth()

        # Face-on
        self.write_young_maps_faceon()

        # Edge-on
        self.write_young_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_young_fuv_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_young_fuv_luminosity_map and not self.has_young_fuv_map

    # -----------------------------------------------------------------

    @property
    def young_fuv_map(self):

        """
        This function ...
        :return:
        """

        return self.model.young_fuv_luminosity_map

    # -----------------------------------------------------------------

    @property
    def young_fuv_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "young_fuv.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_fuv_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_fuv_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_young_bol_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_young_bolometric_luminosity_map and not self.has_young_bol_map

    # -----------------------------------------------------------------

    @property
    def young_bol_map(self):

        """
        This function ...
        :return:
        """

        return self.model.young_bolometric_luminosity_map

    # -----------------------------------------------------------------

    @property
    def young_bol_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "young_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_bol_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_bol_map_path)

    # -----------------------------------------------------------------

    def write_young_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the young maps in the earth projection ...")

        # FUV lum
        if self.do_write_young_fuv_map: self.young_fuv_map.saveto(self.young_fuv_map_path)

        # Bol lum
        if self.do_write_young_bol_map: self.young_bol_map.saveto(self.young_bol_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_young_fuv_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_young_fuv_luminosity_map_faceon and not self.has_young_fuv_map_faceon

    # -----------------------------------------------------------------

    @property
    def young_fuv_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.young_fuv_luminosity_map

    # -----------------------------------------------------------------

    @property
    def young_fuv_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "young_fuv.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_fuv_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_fuv_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_young_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_young_bolometric_luminosity_map_faceon and not self.has_young_bol_map_faceon

    # -----------------------------------------------------------------

    @property
    def young_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.young_bolometric_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def young_bol_map_faceon_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "young_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_bol_map_faceon_path)

    # -----------------------------------------------------------------

    def write_young_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the young maps in the faceon projection ...")

        # FUV lum
        if self.do_write_young_fuv_map_faceon: self.young_fuv_map_faceon.saveto(self.young_fuv_map_faceon_path)

        # Bol lum
        if self.do_write_young_bol_map_faceon: self.young_bol_map_faceon.saveto(self.young_bol_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_young_fuv_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_young_fuv_luminosity_map_edgeon and not self.has_young_fuv_map_edgeon

    # -----------------------------------------------------------------

    @property
    def young_fuv_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.young_fuv_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def young_fuv_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "young_fuv.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_fuv_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_fuv_map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_young_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_young_bolometric_luminosity_map_edgeon and not self.has_young_bol_map_edgeon

    # -----------------------------------------------------------------

    @property
    def young_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.young_bolometric_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def young_bol_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "young_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_young_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_young_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the young maps in the edgeon projection ...")

        # FUV lum
        if self.do_write_young_fuv_map_edgeon: self.young_fuv_map_edgeon.saveto(self.young_fuv_map_edgeon_path)

        # Bol lum
        if self.do_write_young_bol_map_edgeon: self.young_bol_map_edgeon.saveto(self.young_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_fuv_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_sfr_fuv_luminosity_map and not self.has_sfr_fuv_map

    # -----------------------------------------------------------------

    @property
    def sfr_fuv_map(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_fuv_luminosity_map

    # -----------------------------------------------------------------

    @property
    def sfr_fuv_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "sfr_fuv.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_fuv_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_fuv_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_bol_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_sfr_bolometric_luminosity_map and not self.has_sfr_bol_map

    # -----------------------------------------------------------------

    @property
    def sfr_bol_map(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_bolometric_luminosity_map

    # -----------------------------------------------------------------

    @property
    def sfr_bol_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "sfr_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_bol_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_bol_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_star_formation_rate_map and not self.has_sfr_map

    # -----------------------------------------------------------------

    @property
    def sfr_map(self):

        """
        This function ...
        :return:
        """

        return self.model.star_formation_rate_map

    # -----------------------------------------------------------------

    @property
    def sfr_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "sfr.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_sfr_dust_mass_map and not self.has_sfr_dust_mass_map

    # -----------------------------------------------------------------

    @property
    def sfr_dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_dust_mass_map

    # -----------------------------------------------------------------

    @property
    def sfr_dust_mass_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_earth_path, "sfr_dust_mass.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_dust_mass_map_path)

    # -----------------------------------------------------------------

    def write_sfr_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SFR maps ...")

        # Earth
        self.write_sfr_maps_earth()

        # Face-on
        self.write_sfr_maps_faceon()

        # Edge-on
        self.write_sfr_maps_edgeon()

    # -----------------------------------------------------------------

    def write_sfr_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SFR maps in the earth projection ...")

        # FUV lum
        if self.do_write_sfr_fuv_map: self.sfr_fuv_map.saveto(self.sfr_fuv_map_path)

        # Bol lum
        if self.do_write_sfr_bol_map: self.sfr_bol_map.saveto(self.sfr_bol_map_path)

        # Star formation rate
        if self.do_write_sfr_map: self.sfr_map.saveto(self.sfr_map_path)

        # Dust mass
        if self.do_write_sfr_dust_mass_map: self.sfr_dust_mass_map.saveto(self.sfr_dust_mass_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_fuv_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_sfr_fuv_luminosity_map_faceon and not self.has_sfr_fuv_map_faceon

    # -----------------------------------------------------------------

    @property
    def sfr_fuv_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def sfr_fuv_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "sfr_fuv.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_fuv_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_fuv_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_sfr_bolometric_luminosity_map_faceon and not self.has_sfr_bol_map_faceon

    # -----------------------------------------------------------------

    @property
    def sfr_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_bolometric_luminosity_map

    # -----------------------------------------------------------------

    @property
    def sfr_bol_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "sfr_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_bol_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_star_formation_rate_map_faceon and not self.has_sfr_map_faceon

    # -----------------------------------------------------------------

    @property
    def sfr_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.star_formation_rate_map_faceon

    # -----------------------------------------------------------------

    @property
    def sfr_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "sfr.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_map_faceon(self):

        """
        Thois function ...
        :return:
        """

        return fs.is_file(self.sfr_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_sfr_dust_mass_map_faceon and not self.has_sfr_dust_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def sfr_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_dust_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def sfr_dust_mass_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "sfr_dust_mass.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_dust_mass_map_faceon_path)

    # -----------------------------------------------------------------

    def write_sfr_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SFR maps in the faceon projection ...")

        # FUV lum
        if self.do_write_sfr_fuv_map_faceon: self.sfr_fuv_map_faceon.saveto(self.sfr_fuv_map_faceon_path)

        # Bol lum
        if self.do_write_sfr_bol_map_faceon: self.sfr_bol_map_faceon.saveto(self.sfr_bol_map_faceon_path)

        # Star formation rate
        if self.do_write_sfr_map_faceon: self.sfr_map_faceon.saveto(self.sfr_map_faceon_path)

        # Dust mass
        if self.do_write_sfr_dust_mass_map_faceon: self.sfr_dust_mass_map_faceon.saveto(self.sfr_dust_mass_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_fuv_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_sfr_fuv_luminosity_map_edgeon and not self.has_sfr_fuv_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_fuv_map_edgeon(self):

        """
        This function ...
        :return:
        """
        return self.model.sfr_fuv_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_fuv_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "sfr_fuv.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_fuv_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_fuv_map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_sfr_bolometric_luminosity_map_edgeon and not self.has_sfr_bol_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_bolometric_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_bol_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "sfr_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_star_formation_rate_map_edgeon and not self.has_sfr_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.star_formation_rate_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "sfr.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_sfr_dust_mass_map_edgeon and not self.has_sfr_dust_mass_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.sfr_dust_mass_map_edgeon

    # -----------------------------------------------------------------

    @property
    def sfr_dust_mass_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "sfr_dust_mass.fits")

    # -----------------------------------------------------------------

    @property
    def has_sfr_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sfr_dust_mass_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_sfr_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SFR maps in the edgeon projection ...")

        # FUV lum
        if self.do_write_sfr_fuv_map_edgeon: self.sfr_fuv_map_edgeon.saveto(self.sfr_fuv_map_edgeon_path)

        # Bol lum
        if self.do_write_sfr_bol_map_edgeon: self.sfr_bol_map_edgeon.saveto(self.sfr_bol_map_edgeon_path)

        # Star formation rate
        if self.do_write_sfr_map_edgeon: self.sfr_map_edgeon.saveto(self.sfr_map_edgeon_path)

        # Dust mass
        if self.do_write_sfr_dust_mass_map_edgeon: self.sfr_dust_mass_map_edgeon.saveto(self.sfr_dust_mass_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_unevolved_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the unevolved maps ...")

        # Earth
        self.write_unevolved_maps_earth()

        # Face-on
        self.write_unevolved_maps_faceon()

        # Edge-on
        self.write_unevolved_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_unevolved_fuv_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_unevolved_fuv_luminosity_map and not self.has_unevolved_fuv_map

    # -----------------------------------------------------------------

    @property
    def unevolved_fuv_map(self):

        """
        This function ...
        :return:
        """

        return self.model.unevolved_fuv_luminosity_map

    # -----------------------------------------------------------------

    @property
    def unevolved_fuv_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_earth_path, "unevolved_fuv.fits")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_fuv_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.unevolved_fuv_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_unevolved_bol_map(self):

        """
        Thisf unction ...
        :return:
        """

        return self.model.has_unevolved_bolometric_luminosity_map and not self.has_unevolved_bol_map

    # -----------------------------------------------------------------

    @property
    def unevolved_bol_map(self):

        """
        This function ...
        :return:
        """

        return self.model.unevolved_bolometric_luminosity_map

    # -----------------------------------------------------------------

    @property
    def unevolved_bol_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_earth_path, "unevolved_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_bol_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.unevolved_bol_map_path)

    # -----------------------------------------------------------------

    def write_unevolved_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the unevolved maps in the earth projection ...")

        # FUV lum
        if self.do_write_unevolved_fuv_map: self.unevolved_fuv_map.saveto(self.unevolved_fuv_map_path)

        # Bol lum
        if self.do_write_unevolved_bol_map: self.unevolved_bol_map.saveto(self.unevolved_bol_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_unevolved_fuv_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_unevolved_fuv_luminosity_map_faceon and not self.has_unevolved_fuv_map_faceon

    # -----------------------------------------------------------------

    @property
    def unevolved_fuv_map_faceon(self):

        """
        Thisn function ...
        :return:
        """

        return self.model.unevolved_fuv_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def unevolved_fuv_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "unevolved_fuv.fits")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_fuv_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.unevolved_fuv_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_unevolved_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_unevolved_bolometric_luminosity_map_faceon and not self.has_unevolved_bol_map_faceon

    # -----------------------------------------------------------------

    @property
    def unevolved_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.unevolved_bolometric_luminosity_map_faceon

    # -----------------------------------------------------------------

    @property
    def unevolved_bol_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "unevolved_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_bol_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.unevolved_bol_map_faceon_path)

    # -----------------------------------------------------------------

    def write_unevolved_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the unevolved maps in the faceon projection ...")

        # FUV lum
        if self.do_write_unevolved_fuv_map_faceon: self.unevolved_fuv_map_faceon.saveto(self.unevolved_fuv_map_faceon_path)

        # Bol lum
        if self.do_write_unevolved_bol_map_faceon: self.unevolved_bol_map_faceon.saveto(self.unevolved_bol_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_unevolved_fuv_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_unevolved_fuv_luminosity_map_edgeon and not self.has_unevolved_fuv_map_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_fuv_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.unevolved_fuv_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_fuv_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "unevolved_fuv.fits")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_fuv_map_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        return fs.is_file(self.unevolved_fuv_map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_unevolved_bol_map_edgeon(self):

        """
        Thisn function ...
        :return:
        """

        return self.model.has_unevolved_bolometric_luminosity_map_edgeon and not self.has_unevolved_bol_map_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.unevolved_bolometric_luminosity_map_edgeon

    # -----------------------------------------------------------------

    @property
    def unevolved_bol_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "unevolved_bol.fits")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_bol_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.unevolved_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_unevolved_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the unevolved maps in the edgeon projection ...")

        # FUV lum
        if self.do_write_unevolved_fuv_map_edgeon: self.unevolved_fuv_map_edgeon.saveto(self.unevolved_fuv_map_edgeon_path)

        # Bol lum
        if self.do_write_unevolved_bol_map_edgeon: self.unevolved_bol_map_edgeon.saveto(self.unevolved_bol_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps ...")

        # Earth
        self.write_dust_maps_earth()

        # Face-on
        self.write_dust_maps_faceon()

        # Edge-on
        self.write_dust_maps_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_dust_mass_map and not self.has_dust_mass_map

    # -----------------------------------------------------------------

    @property
    def dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return self.model.dust_mass_map

    # -----------------------------------------------------------------

    @property
    def dust_mass_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_earth_path, "dust_mass.fits")

    # -----------------------------------------------------------------

    @property
    def has_dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.dust_mass_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_total_dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return self.model.has_total_dust_mass_map and not self.has_total_dust_mass_map

    # -----------------------------------------------------------------

    @property
    def total_dust_mass_map(self):

        """
        This function ...
        :return:
        """

        return self.model.total_dust_mass_map

    # -----------------------------------------------------------------

    @property
    def total_dust_mass_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_path, "total_dust_mass.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_dust_mass_map(self):

        """
        Thisn function ...
        :return:
        """

        return fs.is_file(self.total_dust_mass_map_path)

    # -----------------------------------------------------------------

    def write_dust_maps_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps in the earth projection ...")

        # Dust mass
        if self.do_write_dust_mass_map: self.dust_mass_map.saveto(self.dust_mass_map_path)

        # Total dust mass
        if self.do_write_total_dust_mass_map: self.total_dust_mass_map.saveto(self.total_dust_mass_map_path)

    # -----------------------------------------------------------------

    @property
    def do_write_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_dust_mass_map_faceon and not self.has_dust_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.dust_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def dust_mass_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "dust_mass_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.dust_mass_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_total_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_total_dust_mass_map_faceon and not self.has_total_dust_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def total_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.model.total_dust_mass_map_faceon

    # -----------------------------------------------------------------

    @property
    def total_dust_mass_map_faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_faceon_path, "total_dust_mass.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_dust_mass_map_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.total_dust_mass_map_faceon_path)

    # -----------------------------------------------------------------

    def write_dust_maps_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps in the faceon projection ...")

        # Dust mass
        if self.do_write_dust_mass_map_faceon: self.dust_mass_map_faceon.saveto(self.dust_mass_map_faceon_path)

        # Total dust mass
        if self.do_write_total_dust_mass_map_faceon: self.total_dust_mass_map_faceon.saveto(self.total_dust_mass_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_dust_mass_map_edgeon and not self.has_dust_mass_map_edgeon

    # -----------------------------------------------------------------

    @property
    def dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.dust_mass_map_edgeon

    # -----------------------------------------------------------------

    @property
    def dust_mass_map_edgeon_path(self):

        """
        This functino ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "dust_mass.fits")

    # -----------------------------------------------------------------

    @property
    def has_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.dust_mass_map_edgeon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_total_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.has_total_dust_mass_map_edgeon and not self.has_total_dust_mass_map_edgeon

    # -----------------------------------------------------------------

    @property
    def total_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.model.total_dust_mass_map_edgeon

    # -----------------------------------------------------------------

    @property
    def total_dust_mass_map_edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_edgeon_path, "total_dust_mass.fits")

    # -----------------------------------------------------------------

    @property
    def has_total_dust_mass_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.total_dust_mass_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_dust_maps_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps in the edgeon projection ...")

        # Dust mass
        if self.do_write_dust_mass_map_edgeon: self.dust_mass_map_edgeon.saveto(self.dust_mass_map_edgeon_path)

        # Total dust mass
        if self.do_write_total_dust_mass_map_edgeon: self.total_dust_mass_map_edgeon.saveto(self.total_dust_mass_map_edgeon_path)

# -----------------------------------------------------------------
