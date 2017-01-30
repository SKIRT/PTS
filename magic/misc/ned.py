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
import re
from collections import defaultdict

import sys
reload(sys)
sys.setdefaultencoding('utf-8')

# Import astronomical modules
from astroquery.ned import Ned
from astroquery import nasa_ads as ads

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.basics.filter import parse_filter, BroadBandFilter
from ...core.basics.configurable import Configurable
from ...core.tools import formatting as fmt

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

        # A regular expression object that strips away special unicode characters, used on the remote console output
        self.ansi_escape = re.compile(r'\x1b[^m]*m')

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Find the images
        self.find()

        # 3. List the images
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
                        filter = parse_filter(band)
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
            elif "Ha" in band or "H-alpha" in band or "H_alph" in band: fltr = BroadBandFilter("Ha")
            else:

                try: fltr = parse_filter(band)
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

                print(fmt.green + fmt.bold + "Found images: (" + str(len(self.images[str(self.config.filter)])) + ")" + fmt.reset)
                print("")

                for bibcode, year, url in self.images[str(self.config.filter)]:

                    name = fs.name(url)

                    results = ads.ADS.query_simple(bibcode)
                    authors = results["authors"][0]

                    if len(authors) == 1: authorstring = authors[0]
                    elif len(authors) == 2: authorstring = "and".join(authors)
                    else: authorstring = authors[0] + " et al."

                    title = results["title"][0][0]
                    journal = results["journal"][0][0].split(",")[0]
                    citations = results["citations"][0][0]

                    print(fmt.underlined + name + fmt.reset)
                    print("")
                    if year is not None: print(" * year:", year)
                    print(" * title:", title)
                    print(" * journal:", journal)
                    print(" * citations:", citations)
                    print(" * authors:", authorstring)
                    print(" * url:", url)
                    print("")

        else:

            # List known
            self.list_filters()

            # List unknown
            if self.config.unknown: self.list_unknown()

    # -----------------------------------------------------------------

    def list_filters(self):

        """
        This function ...
        :return:
        """

        # Loop over the filters
        for fltrstring in self.sorted_filter_names:

            print(fmt.bold + fmt.green + fltrstring + ": (" + str(len(self.images[fltrstring])) + ")" + fmt.reset)
            print("")

            for bibcode, year, url in self.images[fltrstring]:

                name = fs.name(url)

                results = ads.ADS.query_simple(bibcode)
                authors = results["authors"][0]

                if len(authors) == 1: authorstring = authors[0]
                elif len(authors) == 2: authorstring = "and".join(authors)
                else: authorstring = authors[0] + " et al."

                title = results["title"][0][0]
                journal = results["journal"][0][0].split(",")[0]
                citations = results["citations"][0][0]

                print(fmt.underlined + name + fmt.reset)
                print("")
                if year is not None: print(" * year:", year)
                print(" * title:", title)
                print(" * journal:", journal)
                print(" * citations:", citations)
                print(" * authors:", authorstring)
                print(" * url:", url)
                print("")

    # -----------------------------------------------------------------

    def list_unknown(self):

        """
        This function ...
        :return:
        """

        print(fmt.red + fmt.bold + "Unknown filters: (" + str(len(self.unknown)) + ")" + fmt.reset)
        print("")

        for bibcode, year, url in self.unknown:

            name = fs.name(url)

            results = ads.ADS.query_simple(bibcode)

            if len(results) > 0:

                authors = results["authors"][0]

                if len(authors) == 1: authorstring = authors[0]
                elif len(authors) == 2: authorstring = "and".join(authors)
                else: authorstring = authors[0] + " et al."

                title = results["title"][0][0]
                journal = results["journal"][0][0].split(",")[0]
                citations = results["citations"][0][0]

            else: title = journal = citations = authorstring = None

            #print(type(title))
            #print(str(title))
            #print(self.ansi_escape.sub('', str(title)).replace('\x1b[K', '').split("\r\n")[1:-1])
            #print(title.replace("\u2014", "")) if title is not None else print("")

            #print(title.astype('U')) if title is not None else print("")

            #print(title.encode('utf-8')) if title is not None else print("")

            print(fmt.underlined + name + fmt.reset)
            print("")
            if year is not None: print(" * year:", year)
            if title is not None: print(" * title:", title)
            if journal is not None: print(" * journal:", journal)
            if citations is not None: print(" * citations:", citations)
            if authorstring is not None: print(" * authors:", authorstring)
            print(" * url:", url)
            print("")

    # -----------------------------------------------------------------

    @property
    def sorted_filter_names(self):

        """
        This function ...
        :return:
        """

        names = sorted(self.images.keys(), key=lambda key: parse_filter(key).pivot)
        return names

# -----------------------------------------------------------------
