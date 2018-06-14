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
from .component import MapMakingComponent
from ...magic.maps.attenuation.cortese import CorteseAttenuationMapsMaker
from ...magic.maps.attenuation.buat import BuatAttenuationMapsMaker

# -----------------------------------------------------------------

methods = ["cortese", "buat"]

# -----------------------------------------------------------------

class AttenuationMapMaker(MapMakingComponent):

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

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Cortese
        if self.make_cortese: self.make_cortese_attenuation_maps()

        # 3. Buat
        if self.make_buat: self.make_buat_attenuation_maps()

        # 4. Write
        self.write()

        # 5. Plot
        if self.config.plot: self.plot()

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
        maker.config.plot = self.config.debug_plots

        # Get the input
        fuv = self.get_frame_for_filter(self.fuv_filter)

        # Get TIR maps
        tirs, tirs_origins, tirs_methods, tirs_nans = self.select_tir_maps(self.config.tirs, methods=self.config.tir_methods,
                                                                           prompt=self.config.select_tir, title="TIR maps to create Cortese attenuation maps")

        # Get sSFR maps
        ssfrs, ssfrs_origins, ssfrs_methods, ssfrs_nans = self.select_ssfr_maps(self.config.ssfrs, prompt=self.config.select_ssfr,
                                                                                title="sSFR maps to create Cortese attenuation maps", method="colours")

        # Get current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Get current extra maps
        if self.config.remake : current_extra = dict()
        else: current_extra = self.get_current_extra_maps_method(method_name)

        # Run the map maker
        maker.run(fuv=fuv, tirs=tirs, ssfrs=ssfrs, tirs_origins=tirs_origins, ssfrs_origins=ssfrs_origins,
                  tirs_methods=tirs_methods, tirs_nans=tirs_nans, ssfrs_methods=ssfrs_methods, ssfrs_nans=ssfrs_nans,
                  method_name=method_name, maps=current, region_of_interest=self.truncation_ellipse, tir_to_fuvs=current_extra)

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
        #maker.config.plot = self.config.debug_plots

        # Get FUV and NUV maps
        fuv = self.get_frame_for_filter(self.fuv_filter)
        nuv = self.get_frame_for_filter(self.nuv_filter)

        # Get TIR maps, origins, methods and NaN maps
        tirs, tirs_origins, tirs_methods, tirs_nans = self.select_tir_maps(self.config.tirs, methods=self.config.tir_methods,
                                                                           prompt=self.config.select_tir,
                                                                           title="TIR maps to create Buat attenuation maps")

        # Get current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Get current extra maps
        if self.config.remake: current_extra = dict()
        else: current_extra = self.get_current_extra_maps_method(method_name)

        # Run the map maker
        maker.run(fuv=fuv, nuv=nuv, tirs=tirs, tirs_origins=tirs_origins, tirs_methods=tirs_methods,
                  tirs_nans=tirs_nans, method_name=method_name, maps=current, region_of_interest=self.truncation_ellipse,
                  tir_to_fuvs=current_extra)

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

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot
        self.plot_attenuation()

# -----------------------------------------------------------------
