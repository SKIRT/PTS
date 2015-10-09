#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# *****************************************************************

class WCSFinder(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        if config is None:

            directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))

            # Load the default configurations for the star remover
            config_path = os.path.join(directory, "config", "skyextractor.cfg")
            self.config = Config(file(config_path))

        else: self.config = config
        
    # *****************************************************************

    def run(self, frame):

        """
        This function ...
        :return:
        """

        pass
        
# *****************************************************************
