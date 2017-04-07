#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.tir.tir Contains the TIRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component import MapsComponent
from ....core.tools.logging import log
#from .blackbody import BlackBodyDustMapMaker
#from .emission import EmissionDustMapMaker
#from .buat import BuatDustMapMaker
#from .cortese import CorteseDustMapMaker
#from ....core.tools import filesystem as fs
#from .tirtofuv import TIRtoFUVMapMaker
from ....magic.core.image import Image

# -----------------------------------------------------------------

class TIRMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(TIRMapMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The TIR to FUV ratio (in log)
        self.log_tir_to_fuv = None

        # The dust maps (and error maps)
        self.maps = dict()
        self.error_maps = dict()

        # The best dust map
        self.map = None

        # The image of significance masks
        self.significance = Image()

        # The cutoff mask
        self.cutoff_mask = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Calculate the significance masks
        self.calculate_significance()

        # 2. Make the TIR to FUV map
        if self.config.make_buat or self.config.make_cortese: self.make_tir_to_fuv()

        

        # Make the cutoff mask
        self.make_cutoff_mask()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DustMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    

    # -----------------------------------------------------------------

    def write_significance_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the significance masks ...")

        # Write
        self.significance.saveto(self.dust_significance_path)

    # -----------------------------------------------------------------

    def write_cutoff_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cutoff mask ...")

        # Write
        self.cutoff_mask.saveto(self.dust_cutoff_path)

# -----------------------------------------------------------------
