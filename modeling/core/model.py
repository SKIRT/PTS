#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.model Contains the Model class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ..core.mappings import Mappings

# -----------------------------------------------------------------

class Model(object):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        self.simulation_name = None
        self.chi_squared = None
        self.parameter_values = None

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values.keys() if self.parameter_values is not None else None

    # -----------------------------------------------------------------

    @property
    def mappings(self):

        """
        This function ...
        :return:
        """

        # Get the relevant parameters
        metallicity = self.parameter_values["metallicity"]
        compactness = self.parameter_values["sfr_compactness"]
        pressure = self.parameter_values["sfr_pressure"]
        covering_factor = self.parameter_values["sfr_covering"]
        #sfr = self.parameter_values[""] # has to be derived from the MAPPINGS SED from the FUV luminosity
        sfr = 1.

        # Create the MAPPINGS template
        mappings = Mappings(metallicity, compactness, pressure, covering_factor, sfr)

        # Return the mappings template
        return mappings

# -----------------------------------------------------------------
