#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.database Contains the DustPediaDatabase class,
#  which provides an interface to the DustPedia database.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import tempfile
import requests
from lxml import html

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import tables
from ...magic.core.frame import Frame
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from .sample import DustPediaSample
from ...core.filter.filter import parse_filter
from ...core.tools import network
from ...core.tools import progress

# -----------------------------------------------------------------

# The base link
base_link = "http://dustpedia.astro.noa.gr"

# login: http://dustpedia.astro.noa.gr/Account/Login
login_link = "http://dustpedia.astro.noa.gr/Account/Login"

# data: http://dustpedia.astro.noa.gr/Data
data_link = "http://dustpedia.astro.noa.gr/Data"

# user
user_link = "http://dustpedia.astro.noa.gr/Account/UserProfile"

# print preview
print_preview_link = "http://dustpedia.astro.noa.gr/Data/GalaxiesPrintView"

# -----------------------------------------------------------------

class DustPediaDatabase(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Determine the path to a temporary directory
        self.temp_path = tempfile.gettempdir()

        # Create the session
        self.session = requests.session()

        # A flag that states whether we are connected
        self.connected = False

        # DustPedia sample
        self.sample = DustPediaSample()

    # -----------------------------------------------------------------

    def login(self, username, password):

        """
        This function ...
        :param username:
        :param password:
        :return:
        """

        # Inform the user
        log.info("Logging in to the DustPedia database ...")

        r = self.session.get(user_link)
        p = self.session.post(login_link, {'UserName': username, 'password': password})

        # Check login
        r = self.session.get(user_link)

        # Check whether the login was succesful
        self.connected = username in r.content

        # If the login failed, raise an error
        if not self.connected: raise RuntimeError("Login failed")
        else: log.success("Succesfully connected to the DustPedia database")

    # -----------------------------------------------------------------

    def logout(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Logging out from the DustPedia database ...")

        # Disconnect
        if self.connected: self.session.close()

    # -----------------------------------------------------------------

    def __del__(self):

        """
        The destructor ...
        :return:
        """

        # Log out from the database
        self.logout()

    # -----------------------------------------------------------------

    def get_galaxy_names(self, parameters):

        """
        This function ...
        :param parameters:
        :return:
        """

        link = page_link_from_parameters(parameters)

        r = self.session.get(link)
        r = self.session.get(print_preview_link)

        page_as_string = r.content

        tree = html.fromstring(page_as_string)

        tables = [e for e in tree.iter() if e.tag == 'table']
        table = tables[-1]

        table_rows = [e for e in table.iter() if e.tag == 'tr']
        column_headings = [e.text_content() for e in table_rows[0].iter() if e.tag == 'th']

        galaxy_names = []

        for row in table_rows[1:]:

            column_index = 0

            for e in row.iter():

                if e.tag != "td": continue

                galaxy_info = e.text_content()

                galaxy_name = galaxy_info.split("Name: ")[1].split(" \r\n")[0]
                galaxy_names.append(galaxy_name)
                break

        return galaxy_names

    # -----------------------------------------------------------------

    def get_galaxies(self, parameters):

        """
        This function ...
        :param parameters:
        :return:
        """

        name_column = []
        ra_column = []
        dec_column = []
        stage_column = []
        type_column = []
        v_column = []
        d25_column = []
        i_column = []

        link = page_link_from_parameters(parameters)

        r = self.session.get(link)
        r = self.session.get(print_preview_link)

        page_as_string = r.content

        tree = html.fromstring(page_as_string)

        page_tables = [e for e in tree.iter() if e.tag == 'table']
        table = page_tables[-1]

        table_rows = [e for e in table.iter() if e.tag == 'tr']
        column_headings = [e.text_content() for e in table_rows[0].iter() if e.tag == 'th']

        info_boxes = []

        for row in table_rows[1:]:

            column_index = 0

            for e in row.iter():

                if e.tag != "td": continue

                galaxy_info = e.text_content()

                #galaxy_name = galaxy_info.split("Name: ")[1].split(" \r\n")[0]
                #galaxy_names.append(galaxy_name)
                info_boxes.append(galaxy_info)
                break

        for box in info_boxes:

            lines = box.split("\r\n")

            name = None
            ra = None
            dec = None
            stage = None
            type = None
            v = None
            d25 = None
            i = None

            for line in lines:

                if "Name" in line: name = line.split(": ")[1].strip()
                elif "RA(2000)" in line: ra = float(line.split(": ")[1])
                elif "DEC(2000)" in line: dec = float(line.split(": ")[1])
                elif "Hubble Stage(T)" in line: stage = float(line.split(": ")[1])
                elif "Hubble Type: Sab" in line: type = line.split(": ")[1].strip()
                elif "V (km/s)" in line: v = float(line.split(": ")[1])
                elif "D25 (arcmin)" in line: d25 = float(line.split(": ")[1])
                elif "Inclination (deg.)" in line: i = float(line.split(": ")[1])

            if add_hubble_type(type, parameters):

                name_column.append(name)
                ra_column.append(ra)
                dec_column.append(dec)
                stage_column.append(stage)
                type_column.append(type)
                v_column.append(v)
                d25_column.append(d25)
                i_column.append(i)

        # Create the table
        names = ["Name", "RA", "DEC", "Hubble stage", "Hubble type", "V", "D25", "Inclination"]
        data = [name_column, ra_column, dec_column, stage_column, type_column, v_column, d25_column, i_column]
        table = tables.new(data, names)

        return table

    # -----------------------------------------------------------------

    def get_image_urls(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Inform the user
        log.info("Getting the URLs of the available images for galaxy '" + galaxy_name + "' ...")

        # Go to the page
        r = self.session.get(data_link + "?GalaxyName="+galaxy_name+"&tLow=&tHigh=&vLow=&vHigh=&inclLow=&inclHigh=&d25Low=&d25High=&SearchButton=Search")

        page_as_string = r.content

        tree = html.fromstring(page_as_string)

        tables = [ e for e in tree.iter() if e.tag == 'table']
        table = tables[-1]

        table_rows = [ e for e in table.iter() if e.tag == 'tr']
        column_headings = [ e.text_content() for e in table_rows[0].iter() if e.tag == 'th']

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

        # Inform the user
        log.info("Getting the names of the images that are available for galaxy '" + galaxy_name + "' ...")

        # Example link: http://dustpedia.astro.noa.gr/Data/GetImage?imageName=NGC3031_Planck_10600.fits&instrument=Planck

        names = []

        for url in self.get_image_urls(galaxy_name):

            name = url.split("imageName=")[1].split("&instrument")[0]
            names.append(name)

        return names

    # -----------------------------------------------------------------

    def get_image_names_and_urls(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Initialize a dictionary for the urls
        urls = dict()

        # Loop over all the urls
        for url in self.get_image_urls(galaxy_name):

            name = url.split("imageName=")[1].split("&instrument")[0]
            urls[name] = url

        # Return the dictionary of urls
        return urls

    # -----------------------------------------------------------------

    def get_image_filters(self, galaxy_name):

        """
        This function ...
        :return:
        """

        galaxy_name = self.sample.get_name(galaxy_name)

        names = self.get_image_names(galaxy_name)

        filters = []

        for name in names:

            if "Error" in name: continue

            # Skip DSS
            if "DSS" in name and "SDSS" not in name: continue

            # Get the filter
            fltr_string = name.split(galaxy_name + "_")[1].split(".fits")[0]
            fltr = parse_filter(fltr_string)

            # Add the filter
            filters.append(fltr)

        # Return the list of filters
        return filters

    # -----------------------------------------------------------------

    def get_image(self, galaxy_name, image_name):

        """
        This function ...
        :param galaxy_name:
        :param image_name:
        :return:
        """

        # Inform the user
        log.info("Getting the image '" + image_name + "' for galaxy '" + galaxy_name + "' ...")

        # Determine a temporary path for the image file
        local_path = fs.join(self.temp_path, image_name)

        # Download the image to the temporary directory
        self.download_image(galaxy_name, image_name, local_path)

        # Open the image
        frame = Frame.from_file(local_path)

        # Remove the file
        fs.remove_file(local_path)

        # Return the image frame
        return frame

    # -----------------------------------------------------------------

    def download_image(self, galaxy_name, image_name, path):

        """
        This function ...
        :param galaxy_name:
        :param image_name:
        :param path:
        :return:
        """

        # Inform the user
        log.info("Downloading the image '" + image_name + "' for galaxy '" + galaxy_name + "' to '" + path + " ...")

        get_link = None

        for url in self.get_image_urls(galaxy_name):

            link_name = url.split("imageName=")[1].split("&instrument")[0]

            if link_name == image_name:
                get_link = url
                break

        # Download
        network.download_file(get_link, path, progress_bar=log.is_debug(), stream=True, session=self.session)

    # -----------------------------------------------------------------

    def download_images(self, galaxy_name, path):

        """
        This function ...
        :param galaxy_name:
        :param path: directory
        :return:
        """

        # Inform the user
        log.info("Downloading all images for galaxy '" + galaxy_name + "' to '" + path + " ...")

        # Loop over the image URLS found for this galaxy
        for url in self.get_image_urls(galaxy_name):

            # Determine path
            image_name = url.split("imageName=")[1].split("&instrument")[0]
            image_path = fs.join(path, image_name)

            # Download this image
            network.download_file(url, image_path, progress_bar=log.is_debug(), stream=True, session=self.session)

    # -----------------------------------------------------------------

    def get_galaxy_info(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Inform the user
        log.info("Getting general information about galaxy '" + galaxy_name + "' ...")

        # Set the names of the table columns
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

    def download_photometry(self, dir_path):

        """
        This function ...
        :param dir_path:
        :return:
        """

        # main photometry: http://dustpedia.astro.noa.gr/Content/tempFiles/DustPedia_Aperture_Photometry.csv
        # IRAS: http://dustpedia.astro.noa.gr/Content/tempFiles/DustPedia_IRAS_SCANPI.csv
        # Planck: http://dustpedia.astro.noa.gr/Content/tempFiles/DustPedia_Planck_CCS2.csv
        # Release notes: http://dustpedia.astro.noa.gr/Content/tempFiles/Photometry_Notes.pdf

        urls = ["http://dustpedia.astro.noa.gr/Content/tempFiles/DustPedia_Aperture_Photometry.csv",
                "http://dustpedia.astro.noa.gr/Content/tempFiles/DustPedia_IRAS_SCANPI.csv",
                "http://dustpedia.astro.noa.gr/Content/tempFiles/DustPedia_Planck_CCS2.csv"]

        # Download the photometry files
        network.download_files(urls, dir_path)

# -----------------------------------------------------------------

def page_link_from_parameters(parameters):

    """
    This function ...
    :param parameters:
    :return:
    """

    if "T" in parameters:
        tlow = str(parameters["T"][0]) if parameters["T"][0] is not None else ""
        thigh = str(parameters["T"][1]) if parameters["T"][1] is not None else ""
    else:
        tlow = thigh = ""

    if "V" in parameters:
        vlow = str(parameters["V"][0]) if parameters["V"][0] is not None else ""
        vhigh = str(parameters["V"][1]) if parameters["V"][1] is not None else ""
    else:
        vlow = vhigh = ""

    if "inclination" in parameters:
        incllow = str(parameters["inclination"][0]) if parameters["inclination"][0] is not None else ""
        inclhigh = str(parameters["inclination"][1]) if parameters["inclination"][1] is not None else ""
    else:
        incllow = inclhigh = ""

    if "D25" in parameters:
        d25low = str(parameters["D25"][0]) if parameters["D25"][0] is not None else ""
        d25high = str(parameters["D25"][1]) if parameters["D25"][1] is not None else ""
    else:
        d25low = d25high = ""

    return data_link + "?GalaxyName=&tLow=" + tlow + "&tHigh=" + thigh + "&vLow=" + vlow + \
            "&vHigh=" + vhigh + "&inclLow=" + incllow + "&inclHigh=" + inclhigh + "&d25Low=" + d25low + \
            "&d25High=" + d25high + "&SearchButton=Search"

# -----------------------------------------------------------------

def get_account():

    """
    This function ...
    :return:
    """

    return introspection.get_account("dustpedia")

# -----------------------------------------------------------------

def add_hubble_type(hubble_type, parameters):

    """
    This function ...
    :param hubble_type:
    :param parameters:
    :return:
    """

    if "Hubble type" in parameters:

        if isinstance(parameters["Hubble type"], basestring): return hubble_type == parameters["Hubble type"]
        else: return hubble_type in parameters["Hubble type"]

    else: return True

# -----------------------------------------------------------------
