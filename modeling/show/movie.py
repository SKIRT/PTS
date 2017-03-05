#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.show.movie Contains the MovieMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import ShowComponent

# -----------------------------------------------------------------

class MovieMaker(ShowComponent):

    """
    This class makes a movie of a radiative transfer model for a galaxy
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        super(MovieMaker, self).__init__(config, interactive)

# -----------------------------------------------------------------
