#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.attenuation.component Contains the AttenuationAnalysisComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ..component import AnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log

# -----------------------------------------------------------------

class AttenuationAnalysisComponent(AnalysisComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(AttenuationAnalysisComponent, self).__init__(*args, **kwargs)

        # The analysis run
        self.analysis_run = None

    # -----------------------------------------------------------------

    @property
    def attenuation_curve_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.attenuation_path, "curve")

    # -----------------------------------------------------------------

    @property
    def attenuation_map_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.attenuation_path, "map")

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AttenuationAnalysisComponent, self).setup(**kwargs)

        # Load the analysis run
        self.load_run()

    # -----------------------------------------------------------------

    def load_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the analysis run " + self.config.run + " ...")

        # Get the run
        self.analysis_run = self.get_run(self.config.run)

# -----------------------------------------------------------------
