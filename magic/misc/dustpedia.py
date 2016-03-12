#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.dustpedia Contains the DustPedia class, which provides an interface to the DustPedia database.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import requests

# Import the relevant PTS classes and modules
from ...core.tools.logging import log

# -----------------------------------------------------------------

# login: http://dustpedia.astro.noa.gr/Account/Login
# data: http://dustpedia.astro.noa.gr/Data

# -----------------------------------------------------------------

class DustPedia(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get_images(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        pass

# -----------------------------------------------------------------
