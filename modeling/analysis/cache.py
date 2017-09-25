#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.cache Contains the AnalysisRunsCacher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.basics.log import log

# -----------------------------------------------------------------

class AnalysisRunsCacher(AnalysisComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(AnalysisComponent, self).__init__(*args, **kwargs)

        # -- Attributes --

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param simulation:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 6. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisRunsCacher, self).setup()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
