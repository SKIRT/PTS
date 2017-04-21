#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.attenuation.attenuation Contains the AttenuationMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ..component import MapsComponent
from .cortese import CorteseAttenuationMapMaker

# -----------------------------------------------------------------

class AttenuationMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(AttenuationMapMaker, self).__init__(config, interactive)

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

        # 3. Write
        self.write()

    # -----------------------------------------------------------------

    def make_cortese_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making Cortese attenuation maps ...")

        # Create the map maker
        maker = CorteseAttenuationMapMaker()

        # Run the map maker
        maker.run()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
