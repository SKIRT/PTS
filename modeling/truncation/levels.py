#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.truncation.levels Contains the SignificanceLevelsSetter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import TruncationComponent
from ...core.basics.log import log
from ...core.tools.serialization import write_dict
from ...core.basics.configuration import prompt_real

# -----------------------------------------------------------------

class SignificanceLevelsSetter(TruncationComponent):
    
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
        super(SignificanceLevelsSetter, self).__init__(*args, **kwargs)

        # The levels
        self.levels = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Prompt
        if not self.has_levels: self.prompt()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SignificanceLevelsSetter, self).setup(**kwargs)

        # If levels are passed as an argument
        if self.config.levels is not None:

            # Initialize levels dict
            self.levels = dict()

            # Loop over the different image names
            for name in self.prep_names:

                # Check whether level is specified
                if name not in self.config.levels: raise ValueError("The significance level for the '" + name + "' image is not specified")

                # Set the level
                self.levels[name] = self.config.levels[name]

    # -----------------------------------------------------------------

    @property
    def has_levels(self):

        """
        This function ...
        :return:
        """

        return self.levels is not None

    # -----------------------------------------------------------------

    def prompt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the significance level for each galaxy image ...")

        # Initialize
        self.levels = dict()

        # Loop over the galaxy images
        for name in self.prep_names:

            # Set properties
            level_name = name.replace(" ", "_") + "_level"
            description = "significance level for the '" + name + "' image"

            # Prompt
            level = prompt_real(level_name, description, required=True)

            # Set the level
            self.levels[name] = level

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ellipse
        self.write_levels()

    # -----------------------------------------------------------------

    def write_levels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the significance levels ...")

        # Write the dictionary of significance levels
        write_dict(self.levels, self.significance_levels_path)

# -----------------------------------------------------------------
