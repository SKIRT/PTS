#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.attenuation Contains the AttenuationMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import MapsComponent
from ...magic.maps.attenuation.cortese import CorteseAttenuationMapsMaker
from ...magic.maps.attenuation.buat import BuatAttenuationMapsMaker
from ...core.basics.configuration import prompt_string_list

# -----------------------------------------------------------------

methods = ["cortese", "buat"]

# -----------------------------------------------------------------

class AttenuationMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(AttenuationMapMaker, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self.maps_attenuation_path

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Setup
        self.setup(**kwargs)

        # 2. Cortese
        self.make_cortese_attenuation_maps()

        # 3. Buat
        self.make_buat_attenuation_maps()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def make_cortese_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making Cortese attenuation maps ...")

        # Define method name
        method_name = "cortese"

        # Create the map maker
        maker = CorteseAttenuationMapsMaker()

        # Get the input
        fuv = self.get_frame_for_filter(self.fuv_filter)

        # Get maps and origins
        tirs = self.get_tir_maps(flatten=True)
        ssfrs = self.get_ssfr_maps(flatten=True)

        # Select only certain TIR maps
        if self.config.select_tir: tirs = select_maps(tirs, "TIR maps")

        # Select only certain sSFR maps
        if self.config.select_ssfr: ssfrs = select_maps(ssfrs, "sSFR maps")

        # Origins
        tirs_origins = self.get_tir_origins(flatten=True)
        ssfrs_origins = self.get_ssfr_origins(flatten=True)

        # Methods
        tirs_methods = self.get_tir_methods(flatten=True)
        ssfrs_methods = self.get_ssfr_methods(flatten=True)

        # Get current maps
        current = self.get_current_maps_method(method_name)

        # Run the map maker
        maker.run(fuv=fuv, tirs=tirs, ssfrs=ssfrs, tirs_origins=tirs_origins, ssfrs_origins=ssfrs_origins, tirs_methods=tirs_methods, ssfrs_methods=ssfrs_methods, method_name=method_name, maps=current)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

        # Set the extra maps
        self.extra_maps[method_name] = maker.tirtofuvs

    # -----------------------------------------------------------------

    def make_buat_attenuation_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making Buat attenuation maps ...")

        # Define method name
        method_name = "buat"

        # Create the map maker
        maker = BuatAttenuationMapsMaker()

        # Get the input
        fuv = self.get_frame_for_filter(self.fuv_filter)
        nuv = self.get_frame_for_filter(self.nuv_filter)
        tirs = self.get_tir_maps(flatten=True)
        tirs_origins = self.get_tir_origins(flatten=True)
        tirs_methods = self.get_tir_methods(flatten=True)

        # Select only certain TIR maps
        if self.config.select_tir: tirs = select_maps(tirs, "TIR maps")

        #print("TIRS methods", tirs_mesthods)

        # Get current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Run the map maker
        maker.run(fuv=fuv, nuv=nuv, tirs=tirs, tirs_origins=tirs_origins, tirs_methods=tirs_methods, method_name=method_name, maps=current)

        #print("Maker methods", maker.methods.keys())
        #print("keys", maker.maps.keys())

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

        # Set the extra maps
        self.extra_maps[method_name] = maker.tirtofuvs

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps
        self.write_maps()

        # Write the origins
        self.write_origins()

        # Write the methods
        self.write_methods()

        # Write the extra maps
        self.write_extra_maps()

# -----------------------------------------------------------------

def select_maps(maps, title):

    """
    This function ...
    :param maps:
    :param title:
    :return:
    """

    # Select names interactively
    names = prompt_string_list("names", title, choices=maps.keys())

    # New maps
    new_maps = dict()

    # Get selection
    for name in names: new_maps[name] = maps[name]

    # Return the selected maps
    return new_maps

# -----------------------------------------------------------------
