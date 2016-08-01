#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.data.sed Contains the SEDFetcher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .component import DataComponent
from ...core.tools.logging import log
from ...dustpedia.core.database import DustPediaDatabase, get_account

# -----------------------------------------------------------------

class SEDFetcher(DataComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(SEDFetcher, self).__init__(config)

        # The DustPedia database
        self.database = DustPediaDatabase()

        # The SED
        self.sed = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 3. Fetch the SED
        self.fetch_sed()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        super(SEDFetcher, self).setup()

        # Get username and password for the DustPedia database
        if self.config.database.username is not None:
            username = self.config.database.username
            password = self.config.database.password
        else: username, password = get_account()

        # Login to the DustPedia database
        self.database.login(username, password)

    # -----------------------------------------------------------------

    def fetch_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the SED ...")

        # Get the SED
        self.sed = self.database.get_sed(self.ngc_id_nospaces)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        self.write_sed()

    # -----------------------------------------------------------------

    def write_sed(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
