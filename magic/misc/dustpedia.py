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
import urllib
from collections import OrderedDict
import requests
from lxml import html

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import tables
from ..core.image import Image
from ..core.frame import Frame
from ...core.tools import filesystem

# -----------------------------------------------------------------

# The base link
base_link = "http://dustpedia.astro.noa.gr"

# login: http://dustpedia.astro.noa.gr/Account/Login
login_link = "http://dustpedia.astro.noa.gr/Account/Login"

# data: http://dustpedia.astro.noa.gr/Data
data_link = "http://dustpedia.astro.noa.gr/Data"

# user
user_link = "http://dustpedia.astro.noa.gr/Account/UserProfile"

# -----------------------------------------------------------------

class DustPedia(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Create the session
        self.session = requests.session()

    # -----------------------------------------------------------------

    def login(self, username, password):

        """
        This function ...
        :param username:
        :param password:
        :return:
        """

        r = self.session.get(user_link)

        p = self.session.post(login_link, {'UserName': username, 'password': password})

        # Check login
        r = self.session.get(user_link)
        assert username in r.content

    # -----------------------------------------------------------------

    def __del__(self):

        """
        The destructor ...
        :return:
        """

        # Close the session
        self.session.close()

    # -----------------------------------------------------------------

    def get_image_links(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        r = self.session.get(data_link + "?GalaxyName="+galaxy_name+"&tLow=&tHigh=&vLow=&vHigh=&inclLow=&inclHigh=&d25Low=&d25High=&SearchButton=Search")

        page_as_string = r.content

        tree = html.fromstring(page_as_string)

        tables = [ e for e in tree.iter() if e.tag == 'table']
        table = tables[-1]

        table_rows = [ e for e in table.iter() if e.tag == 'tr']
        column_headings =[ e.text_content() for e in table_rows[0].iter() if e.tag == 'th']

        galaxy_info = None
        image_links = []

        for row in table_rows[1:]:

            column_index = 0

            #row_content = []

            for e in row.iter():

                if e.tag != "td": continue

                if column_index == 0:

                    galaxy_info = e.text_content()

                else:

                    #print("CHILDREN")
                    #for ee in e.iterchildren(): print(ee)
                    #print()
                    #print("DESCENDANTS")
                    #for ee in e.iterdescendants(): print(ee)
                    #print()
                    #print("LINKS")
                    for ee in e.iterlinks():

                        link = base_link + ee[2]

                        image_links.append(link)

                    #print()
                    #print("SIBLINGS")
                    #for ee in e.itersiblings(): print(ee)

                column_index += 1

        return image_links

    # -----------------------------------------------------------------

    def get_image_names(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Example link: http://dustpedia.astro.noa.gr/Data/GetImage?imageName=NGC3031_Planck_10600.fits&instrument=Planck

        names = []

        for link in self.get_image_links(galaxy_name):

            name = link.split("imageName=")[1].split("&instrument")[0]
            names.append(name)

        return names

    # -----------------------------------------------------------------

    def get_image(self, galaxy_name, image_name):

        """
        This function ...
        :param galaxy_name:
        :param image_name:
        :return:
        """

        get_link = None

        for link in self.get_image_links(galaxy_name):

            link_name = link.split("imageName=")[1].split("&instrument")[0]

            if link_name == image_name:

                get_link = link
                break

        local_path = filesystem.join(filesystem.home(), "test.fits")

        self.download_image(get_link, local_path)

        # Open the image
        frame = Frame.from_file(local_path)

        # Remove the file
        filesystem.remove_file(local_path)

        # Return the image frame
        return frame

    # -----------------------------------------------------------------

    def get_galaxy_info(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        names = ["Name", "RA", "DEC", "Hubble Stage", "Hubble Type", "V", "D25", "Inclination"]

        r = self.session.get(data_link + "?GalaxyName="+galaxy_name+"&tLow=&tHigh=&vLow=&vHigh=&inclLow=&inclHigh=&d25Low=&d25High=&SearchButton=Search")

        page_as_string = r.content

        tree = html.fromstring(page_as_string)

        table_list = [ e for e in tree.iter() if e.tag == 'table']
        table = table_list[-1]

        table_rows = [ e for e in table.iter() if e.tag == 'tr']

        galaxy_info = None

        for row in table_rows[1:]:

            column_index = 0

            for e in row.iter():

                if e.tag != "td": continue

                if column_index == 0:

                    #for ee in e.iterchildren(): print(ee.text_content())
                    #for ee in e.iterdescendants(): print(ee.text_content())
                    #for ee in e.itersiblings(): print(ee.text_content())

                    galaxy_info = e.text_content()

                    #print(galaxy_info)

                column_index += 1

        splitted = galaxy_info.split("\r\n")

        lines = [split.strip() for split in splitted if split.strip()]

        name = None
        ra = None
        dec = None
        stage = None
        type = None
        v = None
        d25 = None
        i = None

        for line in lines:

            if "Name" in line: name = line.split(": ")[1]
            elif "RA(2000)" in line: ra = float(line.split(": ")[1])
            elif "DEC(2000)" in line: dec = float(line.split(": ")[1])
            elif "Hubble Stage(T)" in line: stage = float(line.split(": ")[1])
            elif "Hubble Type: Sab" in line: type = line.split(": ")[1]
            elif "V (km/s)" in line: v = float(line.split(": ")[1])
            elif "D25 (arcmin)" in line: d25 = float(line.split(": ")[1])
            elif "Inclination (deg.)" in line: i = float(line.split(": ")[1])

        data = [[name], [ra], [dec], [stage], [type], [v], [d25], [i]]

        table = tables.new(data, names)

        return table

    # -----------------------------------------------------------------

    def download_image(self, link, local_path):

        """
        This function ...
        :param link:
        :param local_path:
        :return:
        """

        # NOTE the stream=True parameter
        r = self.session.get(link, stream=True)

        with open(local_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)
                    #f.flush() # commented by recommendation from J.F.Sebastian

# -----------------------------------------------------------------
