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
from ...core.basics.containers import create_subdict

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

        # Define name for extra maps
        self.extra_maps_name = "TIRtoFUV"

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self.maps_attenuation_path

    # -----------------------------------------------------------------

    @property
    def make_cortese(self):

        """
        This function ...
        :return:
        """

        return "cortese" in self.config.methods

    # -----------------------------------------------------------------

    @property
    def make_buat(self):

        """
        This function ...
        :return:
        """

        return "buat" in self.config.methods

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
        if self.make_cortese: self.make_cortese_attenuation_maps()

        # 3. Buat
        if self.make_buat: self.make_buat_attenuation_maps()

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

        # Plot
        maker.config.plot = self.config.plot

        # Get the input
        fuv = self.get_frame_for_filter(self.fuv_filter)

        # Get maps and origins
        tirs = self.get_tir_maps(flatten=True, methods=self.config.tir_methods)
        ssfrs = self.get_ssfr_maps(flatten=True)

        # Origins
        tirs_origins = self.get_tir_origins(flatten=True, methods=self.config.tir_methods)
        ssfrs_origins = self.get_ssfr_origins(flatten=True)

        # Methods
        tirs_methods = self.get_tir_methods(flatten=True, methods=self.config.tir_methods)
        ssfrs_methods = self.get_ssfr_methods(flatten=True)

        # Get only certain TIR maps
        if self.config.tirs is not None:
            tirs = create_subdict(tirs, self.config.tirs)
            tirs_origins = create_subdict(tirs_origins, self.config.tirs)
            tirs_methods = create_subdict(tirs_methods, self.config.tirs)

        # Get only certain sSFR maps
        if self.config.ssfrs is not None:
            ssfrs = create_subdict(ssfrs, self.config.ssfrs)
            ssfrs_origins = create_subdict(ssfrs_origins, self.config.ssfrs)
            ssfrs_methods = create_subdict(ssfrs_methods, self.config.ssfrs)

        # Select only certain TIR maps
        if self.config.select_tir:
            tirs, tir_names = select_maps(tirs, "TIR maps to create Cortese attenuation maps", return_names=True)
            tirs_origins = create_subdict(tirs_origins, tir_names)
            tirs_methods = create_subdict(tirs_methods, tir_names)

        # Select only certain sSFR maps
        if self.config.select_ssfr:
            ssfrs, ssfr_names = select_maps(ssfrs, "sSFR maps to create Cortese attenuation maps", return_names=True)
            ssfrs_origins = create_subdict(ssfrs_origins, ssfr_names)
            ssfrs_methods = create_subdict(ssfrs_methods, ssfr_names)

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

        # Plot
        #maker.config.plot = self.config.plot

        # Get FUV and NUV maps
        fuv = self.get_frame_for_filter(self.fuv_filter)
        nuv = self.get_frame_for_filter(self.nuv_filter)

        # Get TIR maps, origins and methods
        tirs = self.get_tir_maps(flatten=True, methods=self.config.tir_methods)
        tirs_origins = self.get_tir_origins(flatten=True, methods=self.config.tir_methods)
        tirs_methods = self.get_tir_methods(flatten=True, methods=self.config.tir_methods)

        # Get only certain TIR maps
        if self.config.tirs is not None:
            tirs = create_subdict(tirs, self.config.tirs)
            tirs_origins = create_subdict(tirs_origins, self.config.tirs)
            tirs_methods = create_subdict(tirs_methods, self.config.tirs)

        # Select only certain TIR maps
        if self.config.select_tir:
            tirs, tir_names = select_maps(tirs, "TIR maps to create Buat attenuation maps", return_names=True)
            tirs_origins = create_subdict(tirs_origins, tir_names)
            tirs_methods = create_subdict(tirs_methods, tir_names)

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

def select_maps(maps, title, return_names=False):

    """
    This function ...
    :param maps:
    :param title:
    :param return_names:
    :return:
    """

    # Select names interactively
    names = prompt_string_list("names", title, choices=maps.keys())

    # New maps
    new_maps = dict()

    # Get selection
    for name in names: new_maps[name] = maps[name]

    # Return the selected maps
    if return_names: return new_maps, names
    else: return new_maps

# -----------------------------------------------------------------
