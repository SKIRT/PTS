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
from collections import OrderedDict
import requests
from lxml import html

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

        r = requests.get("http://dustpedia.astro.noa.gr/Data?GalaxyName="+galaxy_name+"&tLow=&tHigh=&vLow=&vHigh=&inclLow=&inclHigh=&d25Low=&d25High=&SearchButton=Search")

        page_as_string = r.content

        tree = html.fromstring(page_as_string)

        tables = [ e for e in tree.iter() if e.tag == 'table']
        table = tables[-1]

        table_rows = [ e for e in table.iter() if e.tag == 'tr']
        column_headings =[ e.text_content() for e in table_rows[0].iter() if e.tag == 'th']

        my_results = []
        for row in table_rows[1:]:
            cell_content = [ e.text_content() for e in row.iter() if e.tag == 'td']
            temp_dict = OrderedDict()
            for numb, cell in enumerate(cell_content):
                if numb == 0:
                    temp_dict['row_label'] = cell.strip()
                else:
                    dict_key = column_headings[numb]
                    temp_dict[dict_key] = cell

            my_results.append(temp_dict)

        #return my_results

        #my_results list should only have one entry # = one row of the table (because one galaxy is queried)

        # the orderedDict that represents the row is structured like this:

        # dict["row_label"] = entry of the first column
        # dict["a"] = entry of the second column with title "a"
        # ... for as many columns as there are
        # Here: only 2 columns, "a" = "Galaxy Fits Files"

        splitted = my_results[0]["Galaxy Fits Files"].split("\r\n")

        fits_list = [entry.strip() for entry in splitted if entry.strip()]

        # The list of FITS file names for to the specified galaxy
        return fits_list

# -----------------------------------------------------------------
