#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.launcher This module can be used to launch SKIRT/FitSKIRT simulations in a convenient way

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import inspect
import pxssh

# Import astronomical modules
from astropy import log
import astropy.logger

# Import the relevant PTS modules
from pts import configuration

# -----------------------------------------------------------------

class SkirtRemote(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        ## Configuration

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "skirtremote.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config = configuration.open(config, default_config)

        ## Logging

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

        ## Attributes

        # Create the SSH interface
        self.ssh = pxssh.pxssh()

    # -----------------------------------------------------------------

    def login(self):

        """
        This function ...
        :return:
        """

        # Connect to the remote host
        self.ssh.login(self.config.host, self.config.user, self.config.password)

    # -----------------------------------------------------------------

    def submit(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def status(self):

        """
        This function ..
        :return:
        """

        pass

# -----------------------------------------------------------------