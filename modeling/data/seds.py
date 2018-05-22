#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.data.sedfetching Contains the SEDFetcher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from .component import DataComponent
from ...magic.services.seds import SEDFetcher as _SEDFetcher

# -----------------------------------------------------------------

# INTERESTING LINK FOR HALPHA IMAGES AND FLUXES for CERTAIN GALAXIES: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+AS/137/495&-to=3

# -----------------------------------------------------------------

class SEDFetcher(DataComponent):

    """
    This class ...
    """
    
    def __init__(self, *args, **kwargs):
    
        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SEDFetcher, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The SED fetcher
        self.fetcher = _SEDFetcher()

        # The SEDs
        self.seds = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        """

        # 2. Get the SEDs
        self.get_seds()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SEDFetcher, self).setup(**kwargs)

        # Configure the SED fetcher
        self.fetcher.config.galaxy_name = self.galaxy_name
        self.fetcher.config.catalogs = self.config.catalogs
        self.fetcher.config.list = False
        self.fetcher.config.write = False

    # -----------------------------------------------------------------

    def get_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the SEDs ...")

        # Run the SED fetcher
        self.fetcher.run(ngc_name=self.ngc_name)

        # Get the seds
        self.seds = self.fetcher.seds

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write out the SEDs
        self.write_seds()

    # -----------------------------------------------------------------

    def write_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SEDs ...")

        # Loop over the different SEDs
        for label in self.seds:

            # Debugging info
            log.debug("Writing " + label + " SED ...")

            # Determine the path to the new SED file
            sed_path = fs.join(self.data_seds_path, label + ".dat")

            # Save the SED at the specified location
            self.seds[label].saveto(sed_path)

# -----------------------------------------------------------------
