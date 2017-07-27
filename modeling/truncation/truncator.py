#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.truncation.truncator Contains the Truncator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import TruncationComponent
from ...core.tools.logging import log

# -----------------------------------------------------------------

class Truncator(TruncationComponent):
    
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
        super(Truncator, self).__init__(*args, **kwargs)

        # --- Attributes ---

        # The truncation ellipse
        self.ellipse = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)



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
        super(Truncator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ellipse
        self.write_ellipse()

    # -----------------------------------------------------------------

    def write_ellipse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the truncation ellipse ...")

        # Save the ellipse
        self.ellipse.saveto(self.truncation_ellipse_path)

# -----------------------------------------------------------------
