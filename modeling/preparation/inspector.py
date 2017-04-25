#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.inspector Contains the PreparationInspector class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PreparationComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log

# -----------------------------------------------------------------

class PreparationInspector(PreparationComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(PreparationInspector, self).__init__(config, interactive)

        # Maps of the significance levels
        self.significance_maps = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Inspect errors
        self.inspect_significance()

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
        super(PreparationInspector, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def inspect_significance(self):

        """
        This function ...
        :return:
        """

        # Load all data
        for name in self.dataset.names:

            # Debugging
            log.debug("Calculating significance map for the " + name + " image ...")

            # Add the level map to the dictionary
            self.significance_maps[name] = self.dataset.get_significance(name, self.config.levels)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the uer
        log.info("Writing ...")

        # Write significance maps
        self.write_significance_maps()

    # -----------------------------------------------------------------

    def write_significance_maps(self):

        """
        This function ...
        :return:
        """

        # Loop over the images
        for name in self.significance_maps:

            # Determine the path
            path = fs.join(self.inspect_path, name + "_significance.fits")

            # Save the error level map
            self.significance_maps[name].saveto(path)

# -----------------------------------------------------------------
