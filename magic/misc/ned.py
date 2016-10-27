#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.ned Contains the NED class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import astronomical modules
from astroquery.ned import Ned

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from pts.core.basics.filter import Filter
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

class NED(Configurable):

    """
    This class ..
    """

    def __init__(self, config=None):

        """
        This function ...
        :param config:
        """

        # Call the constructor of the base class
        super(NED, self).__init__(config)

        # Image info
        self.images = defaultdict(list)

        # Images of unknown filters
        self.unknown = []

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        self.find()

        self.list()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        super(NED, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def find(self):

        """
        This function ...
        :return:
        """

        # Get the list
        urls = Ned.get_image_list(self.config.galaxy)

        images = []

        # Print the list
        for url in urls:

            # Get the name
            name = fs.strip_extension(fs.strip_extension(fs.name(url))) # strip both the .gz as the .fits extension

            # Get the bibcode
            try: bibcode = url.split("img/")[1].split("/")[0]
            except IndexError: bibcode = None

            if ":" in name:

                splitted = name.split(":")

                if splitted[0].startswith("NGC_"):
                    band = splitted[0].split("NGC_")[1][5:]
                    try:
                        filter = Filter.from_string(band)
                        splitted = [self.config.galaxy, None, band, splitted[1]]
                    except: pass

                if len(splitted) == 3:

                    splitted = [self.config.galaxy, None, splitted[1], splitted[2]]

                elif len(splitted) == 2:

                    info_and_band = splitted[0].split("NGC_")[1][5:]
                    splitted = [self.config.galaxy, None, info_and_band, splitted[1]]

                galaxy_name = splitted[0]
                unknown = splitted[1]
                band = splitted[2]
                source = splitted[3]

                try:
                    year = int(source[-4:])
                    if year < 1985: continue
                except ValueError: year = None

                images.append((band, year, bibcode, url))

            elif "_" in name:

                splitted = name.split("_")

                band = splitted[-1]

                images.append((band, None, bibcode, url))

            elif "." in name:

                splitted = name.split(".")

                galaxy_name = splitted[0]

                images.append((None, None, bibcode, url))

        # Print
        for band, year, bibcode, url in images:

            if band is None: fltr = None
            elif "Ha" in band or "H-alpha" in band or "H_alph" in band: fltr = Filter.from_string("Ha")
            else:

                try: fltr = Filter.from_string(band)
                except ValueError: fltr = None

            #print(fltr, year, bibcode, url)

            if fltr is None:

                self.unknown.append((bibcode, year, url))

            else:

                fltrstring = str(fltr)

                # Add to the images dictionary
                self.images[fltrstring].append((bibcode, year, url))

    # -----------------------------------------------------------------

    def list(self):

        """
        This function ...
        :return:
        """

        if self.config.filter is not None:

            if str(self.config.filter) not in self.images: log.error("No images for this filter are found")

            else:

                log.success("Found images:")
                log.info("")

                for bibcode, year, url in self.images[str(self.config.filter)]:

                    name = fs.name(url)

                    log.info(" - " + name)

                log.info("")

        else:

            # List known
            self.list_filters()

            # List unknown
            if self.config.list_unknown: self.list_unknown()

    # -----------------------------------------------------------------

    def list_filters(self):

        """
        This function ...
        :return:
        """

        # Loop over the filters
        for fltrstring in self.images:

            #if self.config.filter is not None and str(self.config.filter) != fltrstring: continue

            log.success(fltrstring + ":")
            log.info("")

            for bibcode, year, url in self.images[fltrstring]:

                name = fs.name(url)

                log.info(" - " + name)

            log.info("")

    # -----------------------------------------------------------------

    def list_unknown(self):

        """
        This function ...
        :return:
        """

        log.error("Unknown filters:")
        log.info("")

        for bibcode, year, url in self.unknown:

            name = fs.name(url)

            log.info(" - " + name)

        log.info("")

# -----------------------------------------------------------------
