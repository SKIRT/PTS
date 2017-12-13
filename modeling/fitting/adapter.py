#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.adapter Contains the FittingAdapter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.basics.log import log

# -----------------------------------------------------------------

# Default magnitude ranges
magnitude_ranges = dict()
magnitude_ranges["fuv_young"] = 2
magnitude_ranges["dust_mass"] = 1
magnitude_ranges["fuv_ionizing"] = 3

# -----------------------------------------------------------------

# Default percentual deviation from the initial value
percentual_ranges = dict()
percentual_ranges["distance"] = 20
percentual_ranges["ionizing_scaleheight"] = 25
percentual_ranges["sfr_compactness"] = 25
percentual_ranges["old_scaleheight"] = 25
percentual_ranges["position_angle"] = 10
percentual_ranges["metallicity"] = 50
percentual_ranges["young_scaleheight"] = 25
percentual_ranges["sfr_covering"] = 25
percentual_ranges["dust_scaleheight"] = 25
percentual_ranges["i1_old"] = 10
percentual_ranges["sfr_pressure"] = 25
percentual_ranges["inclination"] = 15

# -----------------------------------------------------------------

default_magnitude_range = 2
default_percentual_range = 20

# -----------------------------------------------------------------

class FittingAdapter(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(FittingAdapter, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The fitting run
        self.fitting_run = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)



        # 13. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingConfigurer, self).setup()

    # -----------------------------------------------------------------



    # -----------------------------------------------------------------

    

# -----------------------------------------------------------------
