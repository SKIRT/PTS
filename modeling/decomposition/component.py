#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.decomposition.component Contains the DecompositionComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class DecompositionComponent(GalaxyModelingComponent):
    
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
        super(DecompositionComponent, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The path to the components/2D directory
        self.components_2d_path = None

        # The path to the components/residuals directory
        self.components_residuals_path = None

        # The paths to the different galaxy projections
        self.earth_projection_path = None
        self.edgeon_projection_path = None
        self.faceon_projection_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DecompositionComponent, self).setup(**kwargs)

        # Set ...
        self.components_2d_path = fs.create_directory_in(self.components_path, "2D")

        # Set ...
        self.components_residuals_path = fs.create_directory_in(self.components_path, "residuals")

        # The paths to the different galaxy projections
        self.earth_projection_path = fs.join(self.components_projections_path, "earth.proj")
        self.edgeon_projection_path = fs.join(self.components_projections_path, "edgeon.proj")
        self.faceon_projection_path = fs.join(self.components_projections_path, "faceon.proj")

# -----------------------------------------------------------------
