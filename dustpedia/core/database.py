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
import requests
from collections import OrderedDict
from lxml import html
import numpy as np

# Import astronomical modules
from astropy.io.fits import getheader
from astropy.units import Unit
from astropy.coordinates import Angle, SkyCoord

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools import tables
from ...magic.core.frame import Frame
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.filter.filter import parse_filter
from ...core.tools import network
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...core.tools import time
from ...core.tools import types
from ...magic.core.list import FrameList, NamedFrameList
from ...magic.core.list import CoordinateSystemList, NamedCoordinateSystemList
from ...core.tools import archive
from ...core.units.parsing import parse_unit as u
from ...core.basics.containers import DefaultOrderedDict
from ...core.tools.utils import lazyproperty
from ...core.basics.map import Map
from ...core.tools.utils import memoize_method

# -----------------------------------------------------------------

# The base link
base_link = "http://dustpedia.astro.noa.gr"

# login: http://dustpedia.astro.noa.gr/Account/Login
login_link = "http://dustpedia.astro.noa.gr/Account/Login"

# data: http://dustpedia.astro.noa.gr/Data
data_link = "http://dustpedia.astro.noa.gr/Data"

# MBB: http://dustpedia.astro.noa.gr/MBB
mbb_link = "http://dustpedia.astro.noa.gr/MBB"

# user
user_link = "http://dustpedia.astro.noa.gr/Account/UserProfile"

# print preview
print_preview_link = "http://dustpedia.astro.noa.gr/Data/GalaxiesPrintView"

# Account page
account_link = "http://dustpedia.astro.noa.gr/Account/UserProfile"

# -----------------------------------------------------------------

# http://dustpedia.astro.noa.gr/Content/tempFiles/mbb/dustpedia_mbb_results.csv

#all_mmb_results_url = "http://dustpedia.astro.noa.gr/Content/tempFiles/mbb/dustpedia_mbb_results.csv"
all_mmb_results_url = "http://dustpedia.astro.noa.gr/Content/tempFiles/mbb/dustpedia_mbb_results_v2.dat"

# emissivity of κλ=8.52 x (250/λ)1.855 cm2/gr adapted to the THEMIS model
# A bootstrap analysis was used to calculate the uncertainties in dust temperatures and masses.

# -----------------------------------------------------------------

all_cigale_results_url = "http://dustpedia.astro.noa.gr/Content/tempFiles/cigale/dustpedia_cigale_results_v5.dat"
cigale_results_filename = fs.name(all_cigale_results_url)

# -----------------------------------------------------------------

# The path to the PTS kernels directory
if not fs.is_directory(introspection.pts_ext_dir): fs.create_directory(introspection.pts_ext_dir)
dustpedia_path = fs.join(introspection.pts_ext_dir, "dustpedia")
if not fs.is_directory(dustpedia_path): fs.create_directory(dustpedia_path)

# -----------------------------------------------------------------

def get_mbb_dust_mass(galaxy_name):

    """
    This function ...
    :param galaxy_name: 
    :return: 
    """

    username, password = get_account()

    database = DustPediaDatabase()
    database.login(username, password)

    dust_mass = database.get_dust_mass_black_body(galaxy_name)
    return dust_mass

# -----------------------------------------------------------------

def get_cigale_dust_mass(galaxy_name):

    """
    This function ...
    :param galaxy_name:
    :return:
    """

    username, password = get_account()

    database = DustPediaDatabase()
    database.login(username, password)

    dust_mass = database.get_dust_mass_cigale(galaxy_name)
    return dust_mass

# -----------------------------------------------------------------

def has_cigale_parameters(galaxy_name):

    """
    This function ...
    :param galaxy_name:
    :return:
    """

    # Check whether has
    database = DustPediaDatabase()
    return database.get_cigale_parameters(galaxy_name)

# -----------------------------------------------------------------

def get_cigale_parameters(galaxy_name):

    """
    This function ...
    :param galaxy_name:
    :return:
    """

    username, password = get_account()
    database = DustPediaDatabase()

    # Log-in
    try: database.login(username, password)
    except requests.ConnectionError:
        log.warning("The database in unavailable. Check your network connection.")

    # Still try to get the parameters
    return database.get_cigale_parameters(galaxy_name)

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
        self.temp_path = introspection.create_temp_dir(time.unique_name("database"))

        # Create the session
        self.session = requests.session()

        # A flag that states whether we are connected
        self.connected = False

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

    @property
    def full_user_name(self):

        """
        This function ...
        :return:
        """

        if not self.connected: raise ValueError("Must be logged in")

        # Get account page text
        r = self.session.get(account_link)
        page_as_string = r.content

        # Get the name
        full_name = page_as_string.split("    Welcome ")[1].split(",")[0]

        # Return the name
        return full_name

    # -----------------------------------------------------------------

    def reset(self, username, password):

        """
        This fucntion ...
        :return:
        """

        # Logout
        self.logout()

        # Create new session
        self.session = requests.session()

        # Login
        self.login(username, password)

    # -----------------------------------------------------------------

    def logout(self):

        """
        This function ...
        :return:
        """

        # Disconnect
        if self.connected:

            # Inform the user
            log.info("Logging out from the DustPedia database ...")

            # Close the session
            self.session.close()

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

        # Inform the user
        log.info("Getting the galaxy names ...")

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

        # Return the table
        return table

    # -----------------------------------------------------------------

    @memoize_method
    def get_image_urls(self, galaxy_name, error_maps=True):

        """
        This function ...
        :param galaxy_name:
        :param error_maps:
        :return:
        """

        #print(self.session.__getstate__())

        # Inform the user
        log.info("Getting the URLs of the available images for galaxy '" + galaxy_name + "' ...")

        # Go to the page
        formatted_galaxy_name = galaxy_name.replace("+", "%2B")
        url = data_link + "?GalaxyName="+formatted_galaxy_name+"&tLow=&tHigh=&vLow=&vHigh=&inclLow=&inclHigh=&d25Low=&d25High=&SearchButton=Search"
        r = self.session.get(url)

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

                        name = fs.name(link)
                        if "_Error" in name and not error_maps: continue

                        image_links.append(link)

                    #print()
                    #print("SIBLINGS")
                    #for ee in e.itersiblings(): print(ee)

                column_index += 1

        # Request other page
        #r = self.session.get(user_link)

        #self.session.prepare_request()

        #print(image_links)

        # Filter out galaxies that match to the specified name (e.g. NGC1351 and NGC1351A)
        galaxy_image_links = []
        for link in image_links:
            image_name = link.split("imageName=")[1].split("&instrument")[0]
            image_galaxy_name = image_name.split("_")[0]
            #print(image_galaxy_name, formatted_galaxy_name)
            if image_galaxy_name != formatted_galaxy_name: continue
            galaxy_image_links.append(link)

        #print(galaxy_image_links)

        # Return
        #return image_links
        return galaxy_image_links

    # -----------------------------------------------------------------

    @memoize_method
    def get_image_names(self, galaxy_name, error_maps=True):

        """
        This function ...
        :param galaxy_name:
        :param error_maps:
        :return:
        """

        # Inform the user
        log.info("Getting the names of the images that are available for galaxy '" + galaxy_name + "' ...")

        # Example link: http://dustpedia.astro.noa.gr/Data/GetImage?imageName=NGC3031_Planck_10600.fits&instrument=Planck

        names = []

        for url in self.get_image_urls(galaxy_name):

            name = url.split("imageName=")[1].split("&instrument")[0]
            if "_Error" in name and not error_maps: continue
            names.append(name)

        return names

    # -----------------------------------------------------------------

    @memoize_method
    def get_image_names_and_urls(self, galaxy_name, error_maps=True):

        """
        This function ...
        :param galaxy_name:
        :param error_maps:
        :return:
        """

        # Initialize a dictionary for the urls
        urls = dict()

        # Loop over all the urls
        for url in self.get_image_urls(galaxy_name):

            # Get the filename
            name = url.split("imageName=")[1].split("&instrument")[0]

            # Skip error frames if requested
            if "_Error" in name and not error_maps: continue

            # Check whether this file is a compressed image
            # Get the bare name (without .gz)
            if archive.is_archive(name): name = archive.bare_name(name)

            # Add an entry to the dictionary
            urls[name] = url

        # Return the dictionary of urls
        return urls

    # -----------------------------------------------------------------

    @memoize_method
    def get_image_name_for_filter(self, galaxy_name, fltr):

        """
        Thisf ucntion ...
        :param galaxy_name: 
        :param fltr: 
        :return: 
        """

        # Convert into fltr
        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        # Loop over all the images
        names = self.get_image_names(galaxy_name, error_maps=False)
        for name in names:

            # Skip DSS
            if "DSS" in name and "SDSS" not in name: continue

            # Get the filter
            fltr_string = name.split(galaxy_name + "_")[1].split(".fits")[0]
            name_fltr = parse_filter(fltr_string)

            # Match
            if fltr == name_fltr: return name

        # No match found
        return None

    # -----------------------------------------------------------------

    @memoize_method
    def get_image_names_and_filters(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Initialize dictionary
        filters = OrderedDict()

        # Loop over all the images
        names = self.get_image_names(galaxy_name, error_maps=False)
        for name in names:

            # Skip DSS
            if "DSS" in name and "SDSS" not in name: continue

            # Get the filter
            fltr_string = name.split(galaxy_name + "_")[1].split(".fits")[0]
            fltr = parse_filter(fltr_string)

            # Set the filter
            filters[name] = fltr

        # Return the dictionary of filters
        return filters

    # -----------------------------------------------------------------

    @memoize_method
    def get_image_names_and_filters_per_observatory(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Initialize dictionary
        observatories = DefaultOrderedDict(dict)

        # Get names and filters
        filters = self.get_image_names_and_filters(galaxy_name)

        # Loop over the filters
        for name in filters:

            # Get the filter
            fltr = filters[name]

            # Add to appropriate part of the dictionary
            observatories[fltr.observatory][name] = fltr

        # Return the dictionary
        return observatories

    # -----------------------------------------------------------------

    @memoize_method
    def get_image_filters(self, galaxy_name):

        """
        This function ...
        :return:
        """

        # Initialize list
        filters = []

        formatted_galaxy_name = galaxy_name.replace("+", "%2B")

        # Loop over the image names
        names = self.get_image_names(galaxy_name, error_maps=False)
        for name in names:

            #print(name, galaxy_name, formatted_galaxy_name)

            # Get the filter
            fltr_string = name.split(formatted_galaxy_name + "_")[1].split(".fits")[0]

            # Skip DSS
            if "DSS" in fltr_string and "SDSS" not in fltr_string: continue

            # Get the filter
            try: fltr = parse_filter(fltr_string)
            except ValueError:
                log.warning("Invalid filter string: '" + fltr_string + "': skipping image ...")
                continue

            # Add the filter
            filters.append(fltr)

        # Return the list of filters
        return filters

    # -----------------------------------------------------------------

    @memoize_method
    def get_image_filters_per_observatory(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Initialize
        observatories = DefaultOrderedDict(list)

        # Get filters
        filters = self.get_image_filters(galaxy_name)

        # Loop over the filters
        for fltr in filters:

            # Add to dict
            observatories[fltr.observatory].append(fltr)

        # Return
        return observatories

    # -----------------------------------------------------------------

    @memoize_method
    def get_image_filters_per_instrument(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Initialize
        instruments = DefaultOrderedDict(list)

        # Get filters
        filters = self.get_image_filters(galaxy_name)

        # Loop over the filters
        for fltr in filters:

            # Add to dict
            instruments[fltr.instrument].append(fltr)

        # Return
        return instruments

    # -----------------------------------------------------------------

    @memoize_method
    def get_observatories(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return self.get_image_filters_per_observatory(galaxy_name).keys()

    # -----------------------------------------------------------------

    @memoize_method
    def get_instruments(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return self.get_image_filters_per_instrument(galaxy_name).keys()

    # -----------------------------------------------------------------

    @memoize_method
    def has_galex(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return "GALEX" in self.get_observatories(galaxy_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_sdss(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return "SDSS" in self.get_instruments(galaxy_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_2mass(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return "2MASS" in self.get_observatories(galaxy_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_spitzer(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return "Spitzer" in self.get_observatories(galaxy_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_irac(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return "IRAC" in self.get_instruments(galaxy_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_mips(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return "MIPS" in self.get_instruments(galaxy_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_wise(self, galaxy_name):

        """
        Thisf unction ...
        :param galaxy_name:
        :return:
        """

        return "WISE" in self.get_observatories(galaxy_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_pacs(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return "Pacs" in self.get_instruments(galaxy_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_spire(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return "SPIRE" in self.get_instruments(galaxy_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_herschel(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return "Herschel" in self.get_observatories(galaxy_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_planck(self, galaxy_name):

        """
        Thisf unction ...
        :param galaxy_name:
        :return:
        """

        return "Planck" in self.get_observatories(galaxy_name)

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

    def get_frame(self, galaxy_name, image_name):

        """
        This function ...
        :param galaxy_name: 
        :param image_name: 
        :return: 
        """

        return self.get_image(galaxy_name, image_name)

    # -----------------------------------------------------------------

    def get_image_for_filter(self, galaxy_name, fltr):

        """
        THis function ...
        :param galaxy_name: 
        :param fltr: 
        :return: 
        """

        # Get the name
        name = self.get_image_name_for_filter(galaxy_name, fltr)

        # Return the image
        return self.get_image(galaxy_name, name)

    # -----------------------------------------------------------------

    def get_frame_for_filter(self, galaxy_name, fltr):

        """
        This function ...
        :param galaxy_name: 
        :param fltr: 
        :return: 
        """

        return self.get_image_for_filter(galaxy_name, fltr)

    # -----------------------------------------------------------------

    def get_images_for_filters(self, galaxy_name, filters):

        """
        This function ...
        :param galaxy_name:
        :param filters: 
        :return: 
        """

        images = []

        # Loop over the filters
        for fltr in filters: images.append(self.get_image_for_filter(galaxy_name, fltr))

        # Return the images
        return images

    # -----------------------------------------------------------------

    def get_frames_for_filters(self, galaxy_name, filters):

        """
        This function ...
        :param galaxy_name: 
        :param filters: 
        :return: 
        """

        return self.get_images_for_filters(galaxy_name, filters)

    # -----------------------------------------------------------------

    def get_framelist_for_filters(self, galaxy_name, filters, named=False):

        """
        This function ...
        :param galaxy_name: 
        :param filters: 
        :param named:
        :return: 
        """

        # Initialize list
        if named: frames = NamedFrameList()
        else: frames = FrameList()

        # Loop over the filters
        for fltr in filters:

            # Get image name
            name = self.get_image_name_for_filter(galaxy_name, fltr)

            # Get the frame
            frame = self.get_frame(galaxy_name, name)

            # Add the frame
            if named: frames.append(frame, name)
            else: frames.append(frame, fltr)

        # Return the frame list
        return frames

    # -----------------------------------------------------------------

    def get_header(self, galaxy_name, image_name):

        """
        This function ...
        :param galaxy_name:
        :param image_name:
        :return:
        """

        # Inform the user
        log.info("Getting the header for the '" + image_name + "' for galaxy '" + galaxy_name + "' ...")

        # Determine a temporary path for the image file
        local_path = fs.join(self.temp_path, image_name)

        # Download the image to the temporary directory
        self.download_image(galaxy_name, image_name, local_path)

        # Load the header
        header = getheader(local_path)

        # Return the header
        return header

    # -----------------------------------------------------------------

    def get_header_for_filter(self, galaxy_name, fltr):

        """
        This function ...
        :param galaxy_name: 
        :param fltr:
        :return: 
        """

        # Get name
        name = self.get_image_name_for_filter(galaxy_name, fltr)

        # Return the header
        return self.get_header(galaxy_name, name)

    # -----------------------------------------------------------------

    def get_headers_for_filters(self, galaxy_name, filters):

        """
        This function ...
        :param galaxy_name: 
        :param filters: 
        :return: 
        """

        headers = []
        for fltr in filters: headers.append(self.get_header_for_filter(galaxy_name, fltr))
        return headers

    # -----------------------------------------------------------------

    def get_wcs(self, galaxy_name, image_name):

        """
        This function ...
        :param galaxy_name:
        :param image_name:
        :return:
        """

        # Inform the user
        log.info("Getting the coordinate system for the '" + image_name + "' for galaxy '" + galaxy_name + "' ...")

        # Get the header
        header = self.get_header(galaxy_name, image_name)

        # Create and return the coordinate system
        return CoordinateSystem(header=header)

    # -----------------------------------------------------------------

    def get_coordinate_system(self, galaxy_name, image_name):

        """
        This function ...
        :param galaxy_name: 
        :param image_name: 
        :return: 
        """

        return self.get_wcs(galaxy_name, image_name)

    # -----------------------------------------------------------------

    def get_coordinate_system_for_filter(self, galaxy_name, fltr):

        """
        This function ...
        :param galaxy_name: 
        :param fltr: 
        :return: 
        """

        # Get name
        name = self.get_image_name_for_filter(galaxy_name, fltr)

        # Return the coordinate system
        return self.get_coordinate_system(galaxy_name, name)

    # -----------------------------------------------------------------

    def get_coordinate_systems_for_filters(self, galaxy_name, filters):

        """
        This function ...
        :param galaxy_name: 
        :param filters: 
        :return: 
        """

        coordinate_systems = []
        for fltr in filters: coordinate_systems.append(self.get_coordinate_system_for_filter(galaxy_name, fltr))
        return coordinate_systems

    # -----------------------------------------------------------------

    def get_coordinate_system_list_for_filter(self, galaxy_name, filters, named=False):

        """
        This function ...
        :param galaxy_name: 
        :param filters: 
        :param named:
        :return: 
        """

        if named: coordinate_systems = NamedCoordinateSystemList()
        else: coordinate_systems = CoordinateSystemList()

        # Loop over the filters
        for fltr in filters:

            # Get name
            name = self.get_image_name_for_filter(galaxy_name, fltr)

            # Get coordinate system
            wcs = self.get_coordinate_system(galaxy_name, name)

            # Add
            if named: coordinate_systems.append(name, wcs)
            else: coordinate_systems.append(wcs, fltr=fltr)

        # Return
        return coordinate_systems

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
        network.download_file(get_link, path, progress_bar=log.is_debug, stream=True, session=self.session)

    # -----------------------------------------------------------------

    def download_image_from_url(self, url, path):

        """
        This function ...
        :param url:
        :param path:
        :return:
        """

        # Download
        filepath = network.download_file(url, path, progress_bar=log.is_debug, stream=True, session=self.session)

        # Check if compressed
        if archive.is_archive(filepath): # based on extension

            # Not compressed
            if not archive.is_compressed(filepath): # NOT ACTUALLY COMPRESSED

                # Warning
                log.warning("Image is not compressed but extension is '" + fs.get_extension(filepath) + "': renaming ...")

                # Determine uncompressed name
                name = archive.bare_name(filepath)

                # Rename the file
                fs.rename_file_path(filepath, name)

            # Compressed
            else:

                # Debugging
                log.debug("Decompressing ...")

                # Decompress
                archive.decompress_file_in_place(filepath, remove=True)

        # Get the filename
        # filename = url.split("imageName=")[1].split("&")[0]

        # OLD CODE: SOMETIMES, THE IMAGE NAME WOULD SAY .GZ BUT NOT ACTUALLY COMPRESSED
        # # Download archive, decompress
        # if archive.is_archive(filename):
        #
        #     # Download to temporary path
        #     compressed_path = fs.join(self.temp_path, filename)
        #
        #     # Remove potential file with same name
        #     if fs.is_file(compressed_path): fs.remove_file(compressed_path)
        #
        #     # Debugging
        #     log.debug("Downloading ...")
        #
        #     # Download compressed file
        #     network.download_file(url, compressed_path, progress_bar=log.is_debug, stream=True, session=self.session)
        #
        #     # Debugging
        #     log.debug("Decompressing ...")
        #
        #     # Decompress
        #     archive.decompress_file(compressed_path, path)
        #
        # # Regular download
        # else: network.download_file(url, path, progress_bar=log.is_debug, stream=True, session=self.session)

    # -----------------------------------------------------------------

    def download_images(self, galaxy_name, path, error_maps=True, instruments=None, not_instruments=None):

        """
        This function ...
        :param galaxy_name:
        :param path: directory
        :param error_maps:
        :param instruments:
        :param not_instruments:
        :return:
        """

        # Inform the user
        log.info("Downloading all images for galaxy '" + galaxy_name + "' to '" + path + " ...")

        # Loop over the image URLS found for this galaxy
        for url in self.get_image_urls(galaxy_name, error_maps=error_maps):

            # Determine instrument
            instrument = url.split("&instrument=")[1].strip()
            if instruments is not None and instrument not in instruments: continue
            if not_instruments is not None and instrument in not_instruments: continue

            # Determine path
            image_name = url.split("imageName=")[1].split("&instrument")[0]
            image_path = fs.join(path, image_name)

            # Download this image
            #network.download_file(url, image_path, progress_bar=log.is_debug, stream=True, session=self.session)

            # Download (and decompress)
            self.download_image_from_url(url, image_path)

    # -----------------------------------------------------------------

    @memoize_method
    def has_galaxy(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Get page
        formatted_galaxy_name = galaxy_name.replace("+", "%2B")
        url = data_link + "?GalaxyName=" + formatted_galaxy_name + "&tLow=&tHigh=&vLow=&vHigh=&inclLow=&inclHigh=&d25Low=&d25High=&SearchButton=Search"
        # print(url)
        r = self.session.get(url)
        page_as_string = r.content

        tree = html.fromstring(page_as_string)

        table_list = [e for e in tree.iter() if e.tag == 'table']
        # print(len(table_list), table_list)

        # Return whether there are tables on the page
        return len(table_list) > 0

    # -----------------------------------------------------------------

    @memoize_method
    def get_galaxy_info(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Inform the user
        log.info("Getting general information about galaxy '" + galaxy_name + "' ...")

        # Get page
        formatted_galaxy_name = galaxy_name.replace("+", "%2B")
        url = data_link + "?GalaxyName="+formatted_galaxy_name+"&tLow=&tHigh=&vLow=&vHigh=&inclLow=&inclHigh=&d25Low=&d25High=&SearchButton=Search"
        #print(url)
        r = self.session.get(url)
        page_as_string = r.content

        tree = html.fromstring(page_as_string)

        table_list = [ e for e in tree.iter() if e.tag == 'table']
        #print(len(table_list), table_list)

        # No result for this galaxy
        if len(table_list) == 0: return None

        table = table_list[-1]

        table_rows = [ e for e in table.iter() if e.tag == 'tr']

        galaxy_info = None

        #print([row.text_content() for row in table_rows])
        #print("nrows", len(table_rows))

        # Get appropriate rows
        if len(table_rows) == 2: rows = table_rows[1:]
        #elif len(table_rows) == 3: rows = [table_rows[0]]
        else: rows = [table_rows[1]]

        for row in rows:

            column_index = 0

            for e in row.iter():

                if e.tag != "td": continue

                if column_index == 0:

                    #for ee in e.iterchildren(): print(ee.text_content())
                    #for ee in e.iterdescendants(): print(ee.text_content())
                    #for ee in e.itersiblings(): print(ee.text_content())

                    galaxy_info = e.text_content()
                    break

                    #print(galaxy_info)

                column_index += 1

        #print(galaxy_info)

        # Get the lines with the galaxy info that are not empty
        splitted = galaxy_info.split("\r\n")
        lines = [split.strip() for split in splitted if split.strip()]

        #print(lines)

        # Initialize variables
        name = None
        ra = None
        dec = None
        stage = None
        type = None
        v = None
        d25 = None
        i = None

        # Loop over the lines
        for line in lines:

            # Galaxy name
            if "Name" in line: name = line.split(": ")[1]

            # Right ascension
            elif "RA(2000)" in line: ra = float(line.split(": ")[1])

            # Declination
            elif "DEC(2000)" in line: dec = float(line.split(": ")[1])

            # Hubble stage
            elif "Hubble Stage(T)" in line:

                splitted = line.split(": ")
                if len(splitted) > 1 and splitted[1].strip() != "": stage = float(splitted[1])
                else: stage = None

            # Hubble stage
            elif "Hubble Type" in line:

                splitted = line.split(": ")
                if len(splitted) > 1 and splitted[1].strip() != "": type = splitted[1].strip()
                else: type = None

            # Velocity
            elif "V (km/s)" in line:

                splitted = line.split(": ")
                if len(splitted) > 1 and splitted[1].strip() != "": v = float(splitted[1])
                else: v = None

            # D25
            elif "D25 (arcmin)" in line:

                splitted = line.split(": ")
                if len(splitted) > 1 and splitted[1].strip() != "": d25 = float(splitted[1])
                else: d25 = None

            # Inclination angle
            elif "Inclination (deg.)" in line:

                splitted = line.split(": ")
                if len(splitted) > 1 and splitted[1].strip() != "": i = float(splitted[1])
                else: i = None

        # Create mapping
        info = Map()
        info.name = name
        info.position = SkyCoord(ra=ra, dec=dec, unit="deg")
        info.stage = stage
        info.type = type
        info.velocity = v * u("km/s") if v is not None else None
        info.d25 = Angle(d25, "arcmin") if d25 is not None else None
        info.inclination = Angle(i, "deg") if i is not None else None

        # Return the info
        return info

    # -----------------------------------------------------------------

    @memoize_method
    def get_galaxy_info_table(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Get the info
        info = self.get_galaxy_info(galaxy_name)

        # Set the names of the table columns
        names = ["Name", "RA", "DEC", "Hubble Stage", "Hubble Type", "V", "D25", "Inclination"]
        data = [[info.name], [info.position.ra.degree], [info.position.dec.degree], [info.stage], [info.type], [info.velocity.to("km/s").value], [info.d25.to("arcmin").value], [info.inclination.degree]]

        # Create and return the table
        table = tables.new(data, names)
        return table

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_black_body_parameters(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return: 
        """

        # Inform the user
        log.info("Getting general information about galaxy '" + galaxy_name + "' ...")

        # http://dustpedia.astro.noa.gr/MBB?GalaxyName=NGC3031&tLow=&tHigh=&vLow=&vHigh=&inclLow=&inclHigh=&d25Low=&d25High=&SearchButton=Search

        formatted_galaxy_name = galaxy_name.replace("+", "%2B")
        r = self.session.get(mbb_link + "?GalaxyName=" + formatted_galaxy_name + "&tLow=&tHigh=&vLow=&vHigh=&inclLow=&inclHigh=&d25Low=&d25High=&SearchButton=Search")

        page_as_string = r.content

        tree = html.fromstring(page_as_string)

        table_list = [e for e in tree.iter() if e.tag == 'table']
        table = table_list[-1]

        table_rows = [e for e in table.iter() if e.tag == 'tr']

        galaxy_info = None

        for row in table_rows[1:]:

            column_index = 0

            for e in row.iter():

                if e.tag != "td": continue

                if column_index == 0:

                    # for ee in e.iterchildren(): print(ee.text_content())
                    # for ee in e.iterdescendants(): print(ee.text_content())
                    # for ee in e.itersiblings(): print(ee.text_content())

                    galaxy_info = e.text_content()

                    # print(galaxy_info)

                column_index += 1

        splitted = galaxy_info.split("\r\n")

        lines = [split.strip() for split in splitted if split.strip()]

        #return lines

        # Dust Temperature (K): 22.6±0.7
        # Dust Mass (M_sun): 4900000±1000000
        # Dust Luminosity (L_sun): 2.10E+09

        temperature = None
        temperature_error = None
        mass = None
        mass_error = None
        luminosity = None
        luminosity_error = None

        #for index in range(len(lines)):
        index = 0
        while index < len(lines):

            line = lines[index]

            if "Dust Temperature" in line:

                next_line = lines[index+1]
                valuestr, errorstr = next_line.split("&plusmn")
                value = float(valuestr)
                error = float(errorstr)

                temperature = value * u("K")
                temperature_error = error * u("K")

                index += 1

            elif "Dust Mass" in line:

                next_line = lines[index+1]
                valuestr, errorstr = next_line.split("&plusmn")
                value = float(valuestr)
                error = float(errorstr)

                mass = value * u("Msun")
                mass_error = error * u("Msun")

                index += 1

            elif "Dust Luminosity" in line:

                next_line = lines[index+1]

                #valuestr, errorstr = next_line.split("&plusmn")
                #value = float(valuestr)
                #error = float(errorstr)

                value = float(next_line)
                error = None

                luminosity = value * u("Lsun")
                #luminosity_error = error * u("Lsun")
                luminosity_error = None

            index += 1

        # Return the parameters
        return mass, mass_error, temperature, temperature_error, luminosity, luminosity_error

    # -----------------------------------------------------------------

    def download_dust_black_body_plot(self, galaxy_name, path):

        """
        This function ...
        :param galaxy_name: 
        :param path:
        :return: 
        """

        # http://dustpedia.astro.noa.gr/Content/Dustpedia_SEDs_THEMIS/NGC3031.png

        url = "http://dustpedia.astro.noa.gr/Content/Dustpedia_SEDs_THEMIS/" + galaxy_name + ".png"

        # Download
        filepath = network.download_file(url, path, session=self.session, progress_bar=log.is_debug)

        # Return the filepath
        return filepath

    # -----------------------------------------------------------------

    def show_dust_black_body_plot(self, galaxy_name):

        """
        This function ...
        :param galaxy_name: 
        :return: 
        """

        # Get
        filepath = self.download_dust_black_body_plot(galaxy_name, self.temp_path)

        # Open the file
        fs.open_file(filepath)

    # -----------------------------------------------------------------

    def download_dust_black_body_table(self, path):

        """
        This function ... 
        :param path: 
        :return: 
        """

        filepath = network.download_file(all_mmb_results_url, path, session=self.session, progress_bar=log.is_debug)
        return filepath

    # -----------------------------------------------------------------

    def get_dust_black_body_table(self):

        """
        This function ...
        :return: 
        """

        filepath = self.download_dust_black_body_table(self.temp_path)
        return tables.from_file(filepath, format="ascii")

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_black_body_table(self):

        """
        This function ...
        :return:
        """

        return self.get_dust_black_body_table()

    # -----------------------------------------------------------------

    @memoize_method
    def get_black_body_galaxy_index(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = tables.find_index(self.dust_black_body_table, galaxy_name, "Name")
        return index

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_mass_black_body(self, galaxy_name):

        """
        This function ....
        :param galaxy_name:
        :return:
        """

        index = self.get_black_body_galaxy_index(galaxy_name)
        value = self.dust_black_body_table["Mdust__Mo"][index] * Unit("Msun")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_mass_error_black_body(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_black_body_galaxy_index(galaxy_name)
        value = self.dust_black_body_table["Mdust_err"][index] * Unit("Msun")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_temperature_black_body(self, galaxy_name):

        """
        This function ....
        :param galaxy_name:
        :return:
        """

        index = self.get_black_body_galaxy_index(galaxy_name)
        value = self.dust_black_body_table["Tdust__K"][index] * Unit("K")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_temperature_error_black_body(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_black_body_galaxy_index(galaxy_name)
        value = self.dust_black_body_table["Tdust_err"][index] * Unit("K")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_luminosity_black_body(self, galaxy_name):

        """
        This fucntion ...
        :param galaxy_name:
        :return:
        """

        index = self.get_black_body_galaxy_index(galaxy_name)
        value = self.dust_black_body_table["Ldust__Lo"][index] * Unit("Lsun")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_luminosity_error_black_body(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_black_body_galaxy_index(galaxy_name)
        value = self.dust_black_body_table["Ldust_err"][index] * Unit("Lsun")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_chi_squared_black_body(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_black_body_galaxy_index(galaxy_name)
        value = self.dust_black_body_table["nchi2"][index]
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def has_dust_black_body_table_parameters(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Get the index of the galaxy in the black body table
        index = self.get_black_body_galaxy_index(galaxy_name)
        if index is None: return False

        # Get first parameter value
        temperature = self.dust_black_body_table["Tdust__K"][index]
        return not np.isnan(temperature)

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_black_body_table_parameters(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Get the index of the galaxy in the black body table
        index = self.get_black_body_galaxy_index(galaxy_name)

        # Get the values
        dust_temperature = self.dust_black_body_table["Tdust__K"][index] * Unit("K")
        dust_temperature_error = self.dust_black_body_table["Tdust_err"][index] * Unit("K")
        dust_luminosity = self.dust_black_body_table["Ldust__Lo"][index] * Unit("Lsun")
        dust_luminosity_error = self.dust_black_body_table["Ldust_err"][index] * Unit("Lsun")
        dust_mass = self.dust_black_body_table["Mdust__Mo"][index] * Unit("Msun")
        dust_mass_error = self.dust_black_body_table["Mdust_err"][index] * Unit("Msun")
        chi_squared = self.dust_black_body_table["nchi2"][index]

        # Create the parameters dictionary
        parameters = OrderedDict()
        parameters["dust_temperature"] = dust_temperature
        parameters["dust_temperature_error"] = dust_temperature_error
        parameters["dust_luminosity"] = dust_luminosity
        parameters["dust_luminosity_error"] = dust_luminosity_error
        parameters["dust_mass"] = dust_mass
        parameters["dust_mass_error"] = dust_mass_error
        parameters["chi_squared"] = chi_squared

        # Return the parameters dictionary
        return parameters

    # -----------------------------------------------------------------

    def download_cigale_table(self, path):

        """
        This function ....
        :param path:
        :return:
        """

        filepath = network.download_file(all_cigale_results_url, path, session=self.session, progress_bar=log.is_debug)
        return filepath

    # -----------------------------------------------------------------

    def get_cigale_table(self):

        """
        This function ...
        :return:
        """

        filepath = fs.join(dustpedia_path, cigale_results_filename)
        if not fs.is_file(filepath): filepath = self.download_cigale_table(dustpedia_path)

        # Load the table and return it
        return tables.from_file(filepath, format="ascii")

    # -----------------------------------------------------------------

    @lazyproperty
    def cigale_table(self):

        """
        This function ...
        :return:
        """

        return self.get_cigale_table()

    # -----------------------------------------------------------------

    @memoize_method
    def get_cigale_galaxy_index(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Get the galaxy index
        index = tables.find_index(self.cigale_table, galaxy_name, "Name")
        return index

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_mass_cigale(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["Mdust__Mo"][index] * Unit("Msun")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_mass_error_cigale(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["Mdust_err"][index] * Unit("Msun")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_sfr_cigale(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["SFR__Mo_per_yr"][index] * Unit("Msun/yr")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_sfr_error_cigale(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["SFR_err"][index] * Unit("Msun/yr")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_stellar_mass_cigale(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["Mstar__Mo"][index] * Unit("Msun")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_stellar_mass_error_cigale(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["Mstar_err"][index] * Unit("Msun")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_stellar_luminosity_cigale(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["Lstar__W"][index] * Unit("W")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_stellar_luminosity_error_cigale(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["Lstar_err"][index] * Unit("W")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_luminosity_cigale(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["Ldust__W"][index] * Unit("W")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_luminosity_error_cigale(self, galaxy_name):

        """
        THis function ....
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["Ldust_err"][index] * Unit("W")
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_fuv_attenuation_cigale(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["FUV_att"][index]
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def get_fuv_attenuation_error_cigale(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")
        value = self.cigale_table["FUV_att_err"][index]
        return value

    # -----------------------------------------------------------------

    @memoize_method
    def has_cigale_parameters(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: return False

        # Get first parameter value
        dust_mass = self.cigale_table["Mdust__Mo"][index]
        return not np.isnan(dust_mass)

    # -----------------------------------------------------------------

    @memoize_method
    def get_cigale_parameters(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Get the index of the galaxy in the Cigale results table
        index = self.get_cigale_galaxy_index(galaxy_name)
        if index is None: raise ValueError("Invalid galaxy name")

        # Get the parameters
        dust_mass = self.cigale_table["Mdust__Mo"][index] * Unit("Msun")
        dust_mass_error = self.cigale_table["Mdust_err"][index] * Unit("Msun")
        stellar_mass = self.cigale_table["Mstar__Mo"][index] * Unit("Msun")
        stellar_mass_error = self.cigale_table["Mstar_err"][index] * Unit("Msun")
        sfr = self.cigale_table["SFR__Mo_per_yr"][index] * Unit("Msun/yr")
        sfr_error = self.cigale_table["SFR_err"][index] * Unit("Msun/yr")
        stellar_luminosity = self.cigale_table["Lstar__W"][index] * Unit("W")
        stellar_luminosity_error = self.cigale_table["Lstar_err"][index] * Unit("W")
        dust_luminosity = self.cigale_table["Ldust__W"][index] * Unit("W")
        dust_luminosity_error = self.cigale_table["Ldust_err"][index] * Unit("W")
        fuv_attenuation = self.cigale_table["FUV_att"][index]
        fuv_attenuation_error = self.cigale_table["FUV_att_err"][index]

        # Create the parameters dictionary
        parameters = OrderedDict()
        parameters["sfr"] = sfr
        parameters["sfr_error"] = sfr_error
        parameters["stellar_mass"] = stellar_mass
        parameters["stellar_mass_error"] = stellar_mass_error
        parameters["dust_mass"] = dust_mass
        parameters["dust_mass_error"] = dust_mass_error
        parameters["stellar_luminosity"] = stellar_luminosity
        parameters["stellar_luminosity_error"] = stellar_luminosity_error
        parameters["dust_luminosity"] = dust_luminosity
        parameters["dust_luminosity_error"] = dust_luminosity_error
        parameters["fuv_attenuation"] = fuv_attenuation
        parameters["fuv_attenuation_error"] = fuv_attenuation_error

        # Return the parameters dictionary
        return parameters

    # -----------------------------------------------------------------

    def get_photometry_cutouts_url(self, galaxy_name):

        """
        This fucntion ...
        :param galaxy_name:
        :return:
        """

        # http://dustpedia.astro.noa.gr/Data/GetImage?imageName=ESO097-013_Thumbnail_Grid.png&mode=photometry
        url = "http://dustpedia.astro.noa.gr/Data/GetImage?imageName=" + galaxy_name + "_Thumbnail_Grid.png&mode=photometry"
        return url

    # -----------------------------------------------------------------

    def download_photometry_cutouts(self, galaxy_name, dir_path):

        """
        This function ...
        :param galaxy_name:
        :param dir_path:
        :return:
        """

        # Inform the user
        log.info("Downloading the photometry cutouts for galaxy '" + galaxy_name + "' ...")

        url = self.get_photometry_cutouts_url(galaxy_name)
        return network.download_file(url, dir_path, session=self.session, progress_bar=log.is_debug)

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
