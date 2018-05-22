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
from ...core.basics.log import log
from ...core.tools import time

# -----------------------------------------------------------------

class PreparationInspector(PreparationComponent):
    
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
        super(PreparationInspector, self).__init__(*args, **kwargs)

        # Maps of the significance levels
        self.significance_maps = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inspect paths
        self.inspect_paths()

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

        # Set the output path
        directory_name = time.unique_name(self.command_name())
        self.config.output = fs.create_directory_in(self.inspect_path, directory_name)

    # -----------------------------------------------------------------

    def inspect_paths(self):
        
        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Inspecting the image paths ...")

        #self.get_prep_path()

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

        # Infomr the user
        log.info("Writing significance maps ...")

        # Loop over the images
        for name in self.significance_maps:

            # Determine the path
            path = self.output_path_file(name + "_significance.fits")

            # Save the error level map
            self.significance_maps[name].saveto(path)

# -----------------------------------------------------------------
