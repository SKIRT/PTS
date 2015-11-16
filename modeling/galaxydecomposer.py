#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy import log
import astropy.logger

# Import the relevant PTS modules
from ..tools import configuration

# *****************************************************************

class GalaxyDecomposer(object):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        ## Configuration

        self.config = configuration.set("galaxydecomposer", config)

        ## Logging

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

    # *****************************************************************

    def run(self):

        """
        This function ...
        :return:
        """

        pass

# *****************************************************************
