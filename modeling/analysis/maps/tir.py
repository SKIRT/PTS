#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.tir Contains the TIRMapsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsAnalysisComponent
from ....core.basics.log import log
from ....magic.maps.tir.single import SingleBandTIRMapMaker
from ....magic.maps.tir.multi import MultiBandTIRMapMaker

# -----------------------------------------------------------------

class TIRMapsAnalyser(MapsAnalysisComponent):

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
        super(TIRMapsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Make the maps
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TIRMapsAnalyser, self).setup(**kwargs)

        # Load the analysis run
        self.load_run()

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making TIR map ...")

        # 2. Make maps based on a single band
        self.make_tir_maps_single()

        # 3. Make maps based on multiple bands
        self.make_tir_maps_multi()

    # -----------------------------------------------------------------

    def make_tir_maps_single(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Making maps based on a single band ...")

        # Set the method name
        method_name = "single"

        # Create
        maker = SingleBandTIRMapMaker()

        # Get current
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Run
        frames, errors = self.load_data_singleband()
        maker.run(frames=frames, errors=errors, maps=current, method_name=method_name)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def make_tir_maps_multi(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making maps based on multiple bands ...")

        # Set method name
        method_name = "multi"

        # Create
        maker = MultiBandTIRMapMaker()

        # Get current
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Run
        frames, errors = self.load_data_multiband()
        maker.run(frames=frames, errors=errors, maps=current, method_name=method_name)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------



# -----------------------------------------------------------------
