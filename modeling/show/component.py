#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.show.component Contains the ShowComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class ShowComponent(GalaxyModelingComponent):
    
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
        super(ShowComponent, self).__init__(*args, **kwargs)

        # -- Attributes --

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ShowComponent, self).setup(**kwargs)

# -----------------------------------------------------------------
