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
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...magic.services.s4g import S4G

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

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(S4GDecomposer, self).__init__(*args, **kwargs)

        # The S4G interface
        self.s4g = S4G()

        # The path for the S4G decomposition models
        self.components_2d_s4g_path = None

        # The dictionary of components
        self.components = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the components
        self.get_components()

        # 3. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(S4GDecomposer, self).setup(**kwargs)

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
        self.s4g.run(inclination=self.galaxy_inclination)

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
