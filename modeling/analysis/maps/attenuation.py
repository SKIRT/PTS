#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.attenuation Contains the AttenuationMapsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsAnalysisComponent
from ....core.basics.log import log
from ....magic.maps.attenuation.cortese import CorteseAttenuationMapsMaker
from ....magic.maps.attenuation.buat import BuatAttenuationMapsMaker

# -----------------------------------------------------------------

class AttenuationMapsAnalyser(MapsAnalysisComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(AttenuationMapsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # 2. Make the maps
        self.make_maps()

        # 3. Write
        self.write()

        # 4. Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AttenuationMapsAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making attenuation maps ...")

        # 1. Cortese
        self.make_cortese_attenuation_maps()

        # 2. Buat
        self.make_buat_attenuation_maps()

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

        # Origins
        tirs_origins = self.get_tir_origins(flatten=True)
        ssfrs_origins = self.get_ssfr_origins(flatten=True)

        # Methods
        tirs_methods = self.get_tir_methods(flatten=True)
        ssfrs_methods = self.get_ssfr_methods(flatten=True)

        # Get current maps
        current = self.get_current_maps_method(method_name)

        # Run the map maker
        maker.run(fuv=fuv, tirs=tirs, ssfrs=ssfrs, tirs_origins=tirs_origins, ssfrs_origins=ssfrs_origins,
                  tirs_methods=tirs_methods, ssfrs_methods=ssfrs_methods, method_name=method_name, maps=current)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def make_buat_attenuation_maps(self):

        """
        This fucntion ...
        :return:
        """

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

        # print("TIRS methods", tirs_methods)

        # Get current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Run the map maker
        maker.run(fuv=fuv, nuv=nuv, tirs=tirs, tirs_origins=tirs_origins, tirs_methods=tirs_methods,
                  method_name=method_name, maps=current)

        # print("Maker methods", maker.methods.keys())
        # print("keys", maker.maps.keys())

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return:
        """

        return self.attenuation_path

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the colour maps
        self.write_maps()

        # Write origins
        self.write_origins()

        # Write the methods
        self.write_methods()

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

# -----------------------------------------------------------------
