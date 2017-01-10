#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.decomposition.s4g Contains the S4GDecomposer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import DecompositionComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...magic.misc.s4g import S4G

# -----------------------------------------------------------------

# Online table url
s4g_decomposition_table_link = "http://www.oulu.fi/astronomy/S4G_PIPELINE4/s4g_p4_table8.dat"

# Local table path
local_table_path = fs.join(introspection.pts_dat_dir("modeling"), "s4g", "s4g_p4_table8.dat")

# -----------------------------------------------------------------

class S4GDecomposer(DecompositionComponent):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(S4GDecomposer, self).__init__(config)

        # The S4G interface
        self.s4g = S4G()

        # The path for the S4G decomposition models
        self.components_2d_s4g_path = None

        # The dictionary of components
        self.components = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Inform the user
        log.info("Getting the structural galaxy parameters from the S4G catalog ...")

        # Get the components
        self.get_components()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(S4GDecomposer, self).setup()

        # Set the path
        self.components_2d_s4g_path = fs.create_directory_in(self.components_2d_path, "S4G")

        # Set the galaxy name
        self.s4g.config.galaxy_name = self.galaxy_name

        # Don't show results on console
        self.s4g.config.show = False

    # -----------------------------------------------------------------

    def get_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galaxy components from S4G ...")

        # Run
        self.s4g.run()

        # Get the components
        self.components = self.s4g.components

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the components
        self.write_components()

    # -----------------------------------------------------------------

    def write_components(self):

        """
        This function ...
        :return:
        """

        # Loop over the components
        for name in self.components:

            # Determine the path
            path = fs.join(self.components_2d_s4g_path, name + ".mod")

            # Save the model
            self.components[name].saveto(path)

# -----------------------------------------------------------------
