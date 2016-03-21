#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.component Contains the FittingComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem

# -----------------------------------------------------------------

class FittingComponent(ModelingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(FittingComponent, self).__init__(config)

        # -- Attributes --

        # The path to the fit/in directory
        self.fit_in_path = None

        # The path to the fit/out directory
        self.fit_out_path = None

        # The path to the ski file
        self.fit_ski_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingComponent, self).setup()

        # Set the output path
        self.config.output_path = self.fit_path

        # Set the path to the fit/in path
        self.fit_in_path = filesystem.join(self.fit_path, "in")

        # Set the path to the fit/out path
        self.fit_out_path = filesystem.join(self.fit_path, "out")

        # Create the fit/in and fit/out directories
        filesystem.create_directories([self.fit_in_path, self.fit_out_path])

        # Determine the path to the ski file
        self.fit_ski_path = filesystem.join(self.fit_path, self.galaxy_name + ".ski")

# -----------------------------------------------------------------
