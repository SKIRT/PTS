#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.generation Contains the Generation class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules

# -----------------------------------------------------------------

class GenerationInfo(object):

    """
    This function ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Set info
        self.name = kwargs.pop("name", None)
        self.index = kwargs.pop("index", None)
        self.method = kwargs.pop("method", None)
        self.wavelength_grid_level = kwargs.pop("wavelength_grid_level", None)
        #self.dust_grid_level = kwargs.pop("dust_grid_level", None)
        self.model_representation = kwargs.pop("model_representation", None)
        self.nsimulations = kwargs.pop("nsimulations", None)
        self.npackages = kwargs.pop("npackages", None)
        self.selfabsorption = kwargs.pop("selfabsorption", None)
        self.transient_heating = kwargs.pop("transient_heating", None)

        self.path = kwargs.pop("path", None)
        self.individuals_table_path = kwargs.pop("individuals_table_path", None)
        self.parameters_table_path = kwargs.pop("parameters_table_path", None)
        self.chi_squared_table_path = kwargs.pop("chi_squared_table_path", None)

# -----------------------------------------------------------------

class Generation(object):
    
    """
    This class...
    """
    
    def __init__(self):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        #super(FittingComponent, self).__init__(config)

        # -- Attributes --

        pass

# -----------------------------------------------------------------
