#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

# Import standard modules


# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ..core.datacube import DataCube
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

class BatchImagePreparer(Configurable):
    
    """
    This class ...
    """
    
    def __init__(self, config=None):
        
        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(BatchImagePreparer, self).__init__(config)

        # Frames / or do we want a datacube ...
        self.frames = []

    # -----------------------------------------------------------------

    def add_frame(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        self.frames = []

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        self.setup()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------
