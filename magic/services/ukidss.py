#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.ukidss Contains the UKIDSS class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astroquery.ukidss import Ukidss
from astroquery.vizier import Vizier

# Import the relevant PTS classes and modules
from ...core.tools import sequences
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.filter.filter import parse_filter
from ..tools import catalogs
from ..basics.coordinate import SkyCoordinate
from ...core.tools import network
from ...core.units.parsing import parse_unit as u
from ...core.tools import tables

# -----------------------------------------------------------------

bands = ['all', 'J', 'H', 'K', 'H2', 'Z', 'Y', 'Br']

# -----------------------------------------------------------------

class UKIDSS(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        #self.ukidss = Ukidss(username='xyz', password='secret', community='your_community', database='UKIDSSDR8PLUS', programme_id='GPS')
        self.database = Ukidss(database="UKIDSSDR10PLUS", programme_id="all")

        # The Vizier querying object
        self.vizier = Vizier()
        self.vizier.ROW_LIMIT = -1

    # -----------------------------------------------------------------

    @property
    def catalogs(self):

        """
        This function ...
        :return:
        """

        short = Ukidss.list_catalogs(style="short")
        long = Ukidss.list_catalogs(style="long")
        return sequences.zip_into_dict(short, long)

    # -----------------------------------------------------------------

    @property
    def databases(self):

        """
        This function ...
        :return:
        """

        return Ukidss.list_databases()

    # -----------------------------------------------------------------

    def get_position(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        result = self.vizier.query_object(galaxy_name, catalog=["VII/237"])
        table = result[0]

        if len(table) > 1: raise ValueError("Ambiguous result")

        # Get RA and DEC
        ra = table["_RAJ2000"][0]
        dec = table["_DEJ2000"][0]

        # Create sky coordinate and return it
        return SkyCoordinate(ra=ra, dec=dec, unit="deg")

    # -----------------------------------------------------------------

    def get_d25(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        result = self.vizier.query_object(galaxy_name, catalog=["VII/237"])
        table = result[0]

        if len(table) > 1: raise ValueError("Ambiguous result")

        # Calculate the diameter
        diameter = np.power(10.0, table["logD25"][0]) * 0.1 * u("arcmin") if table["logD25"][0] else None

        # Return
        return diameter

    # -----------------------------------------------------------------

    def get_r25(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        result = self.vizier.query_object(galaxy_name, catalog=["VII/237"])
        table = result[0]

        if len(table) > 1: raise ValueError("Ambiguous result")

        # Calculate the ratio
        ratio = np.power(10.0, table["logR25"][0]) if table["logR25"][0] else None

        # Return the ratio
        return ratio

    # -----------------------------------------------------------------

    def get_hyperleda_name(self, galaxy_name): # gets name = HYPERLEDA name, and checks whether in DustPedia sample

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Get the HYPERLEDA name
        objname = catalogs.get_hyperleda_name(galaxy_name)
        return objname

    # -----------------------------------------------------------------

    def get_image_urls(self, galaxy_name, band="all"):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :return:
        """

        # image_urls = Ukidss.get_image_list(coord.SkyCoord(ra=83.633083, dec=22.0145, unit=(u.deg, u.deg), frame='icrs'), frame_type='interleave', programme_id="GCS", waveband="K", radius=20*u.arcmin)
        # get_image_list(self, coordinates, waveband='all', frame_type='stack',
        #               image_width=1 * u.arcmin, image_height=None,
        #               radius=None, database='UKIDSSDR7PLUS',
        #               programme_id='all', get_query_payload=False):

        #radius = self.get_d25(galaxy_name)
        radius = None
        width = self.get_d25(galaxy_name)
        urls = self.database.get_image_list(galaxy_name, waveband=band, frame_type="all", radius=radius, image_width=width)
        return urls

    # -----------------------------------------------------------------

    def get_image_table(self, galaxy_name, band="all"):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :return:
        """

        urls = self.get_image_urls(galaxy_name, band)

        column_names = ["Name", "Path", "Band", "mfid", "extNo", "lx", "hx", "ly", "hy", "rf", "flip", "uniq", "xpos", "ypos", "ra", "dec"]

        names = []
        paths = []
        bands = []
        mfids = []
        extnos = []
        lxs = []
        hxs = []
        lys = []
        hys = []
        rfs = []
        flips = []
        uniqs = []
        xposs = []
        yposs = []
        ras = []
        decs = []

        for url in urls:

            path = url.split("?file=")[1].split("&")[0]
            name = fs.name(path)
            names.append(name)
            paths.append(path)

            band = url.split("&band=")[1].split("&")[0]
            bands.append(band)

            mfid = url.split("&mfid=")[1].split("&")[0]
            mfids.append(mfid)

            extno = url.split("&extNo=")[1].split("&")[0]
            extnos.append(extno)

            lx = url.split("&lx=")[1].split("&")[0]
            lxs.append(lx)

            hx = url.split("&hx=")[1].split("&")[0]
            hxs.append(hx)

            ly = url.split("&ly=")[1].split("&")[0]
            lys.append(ly)

            hy = url.split("&hy=")[1].split("&")[0]
            hys.append(ly)

            rf = url.split("&rf=")[1].split("&")[0]
            rfs.append(rf)

            flip = url.split("&flip=")[1].split("&")[0]
            flips.append(flip)

            uniq = url.split("&uniq=")[1].split("&")[0]
            uniqs.append(uniq)

            xpos = url.split("&xpos=")[1].split("&")[0]
            xposs.append(xpos)

            ypos = url.split("&ypos=")[1].split("&")[0]
            yposs.append(ypos)

            ra = url.split("&ra=")[1].split("&")[0]
            ras.append(ra)

            dec = url.split("&dec=")[1].split("&")[0]
            decs.append(dec)

        data = [names, paths, bands, mfids, extnos, lxs, hxs, lys, hys, rfs, flips, uniqs, xposs, yposs, ras, decs]
        table = tables.new(data, column_names)
        return table

    # -----------------------------------------------------------------

    def get_image_paths(self, galaxy_name, band="all"):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :return:
        """

        urls = self.get_image_urls(galaxy_name, band=band)
        paths = [url.split("?file=")[1].split("&")[0] for url in urls]
        return paths

    # -----------------------------------------------------------------

    def get_image_names(self, galaxy_name, band="all"):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :return:
        """

        return [fs.name(path) for path in self.get_image_paths(galaxy_name, band)]

    # -----------------------------------------------------------------

    def get_image_filter_names(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        filter_names = ["UKIDSS " + url.split("&band=")[1].split("&")[0] for url in self.get_image_urls(galaxy_name)]
        return list(set(filter_names))

    # -----------------------------------------------------------------

    def get_image_filters(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        return [parse_filter(name) for name in self.get_image_filter_names(galaxy_name)]

    # -----------------------------------------------------------------

    def get_image_names_and_urls(self, galaxy_name, band="all"):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get_image_names_and_filters(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def download_images(self, galaxy_name, path, band="all"):

        """
        This function ...
        :param galaxy_name:
        :param path:
        :param band:
        :return:
        """

        # Get urls
        urls = self.get_image_urls(galaxy_name, band)

        # Download each
        for url in urls:
            image_path = url.split("?file=")[1].split("&")[0]
            name = fs.name(image_path)
            self.download_image(url, name, path)

    # -----------------------------------------------------------------

    def download_image(self, url, name, path):

        """
        This function ...
        :param url:
        :param name:
        :param path:
        :return:
        """

        # Inform the user
        log.info("Downloading the image '" + url + "' to '" + path + " ...")

        # Download
        network.download_file(url, path, new_name=name, progress_bar=log.is_debug, stream=True)

# -----------------------------------------------------------------
