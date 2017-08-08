#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.properties Contains the DustPediaProperties class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools import tables
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.filter.filter import parse_filter
from ...core.basics.map import Map
from ...core.units.parsing import parse_unit as u
from .database import DustPediaDatabase
from .sample import DustPediaSample
from ...core.basics.table import SmartTable
from ...magic.basics.pixelscale import Pixelscale
from ...core.tools import types

# -----------------------------------------------------------------

# Load properties
properties_path = fs.join(introspection.pts_dat_dir("dustpedia"), "properties.dat")

# -----------------------------------------------------------------

calibration_references = dict()
calibration_references["a"] = ["http://dx.doi.org/10.1086/520512"]
calibration_references["b"] = ["https://www.sdss3.org/dr12/scope.php"]
calibration_references["c"] = ["http://dx.doi.org/10.1086/374362"]
calibration_references["d"] = ["http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html"]
calibration_references["e"] = ["https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/17/#_Toc410728305"]
calibration_references["f"] = ["https://irsa.ipac.caltech.edu/data/SPITZER/docs/mips/mipsinstrumenthandbook/42/#_Toc288032317"]
calibration_references["g"] = ["http://herschel.esac.esa.int/twiki/bin/view/Public/PacsCalibrationWeb"]
calibration_references["h"] = ["http://herschel.esac.esa.int/twiki/bin/view/Public/SpireCalibrationWeb"]
calibration_references["i"] = ["http://arxiv.org/abs/1502.01587"]
calibration_references["j"] = ["https://arxiv.org/abs/1505.08022"]
calibration_references["k"] = ["http://dx.doi.org/10.1086/376841", "http://dx.doi.org/10.1086/427938"]

# -----------------------------------------------------------------

reference_for_instrument = dict()
reference_for_instrument["GALEX"] = "a"
reference_for_instrument["SDSS"] = "b"
reference_for_instrument["2MASS"] = "c"
reference_for_instrument["WISE"] = "d"
reference_for_instrument["IRAC"] = "e"
reference_for_instrument["Mips"] = "f"
reference_for_instrument["Pacs"] = "g"
reference_for_instrument["SPIRE"] = "h"
reference_for_instrument["HFI"] = "i"
reference_for_instrument["LFI"] = "j"
reference_for_instrument["IRAS"] = "k"

# -----------------------------------------------------------------

def get_fwhm(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    properties = DustPediaProperties()
    return properties.get_fwhm(fltr)

# -----------------------------------------------------------------

class DustPediaProperties(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Initialize the data structure
        self.data = OrderedDict()

        # Load the data
        self.load_data()

    # -----------------------------------------------------------------

    def load_data(self):

        """
        THis function ...
        :return: 
        """

        properties = tables.from_file(properties_path, format="ascii.commented_header")

        # Convert properties in dictionary per filter
        for index in range(len(properties)):

            facility = properties["Facility"][index]
            band = properties["Band"][index]

            # Skip DSS
            if facility == "DSS1" or facility == "DSS2": continue

            filterstring = facility + " " + band
            fltr = parse_filter(filterstring)

            # Get the properties
            pixelscale = properties["Pixelsize"][index] * u("arcsec")
            fwhm = properties["FWHM"][index] * u("arcsec")
            calibration = properties["Calibration error"][index]

            # Set properties
            self.data[fltr] = Map(pixelscale=pixelscale, fwhm=fwhm, calibration=calibration)

    # -----------------------------------------------------------------

    def has_filter(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        return fltr in self.data

    # -----------------------------------------------------------------

    def create_table(self, galaxy_name, username=None, password=None):

        """
        This function ...
        :param galaxy_name: 
        :param username:
        :param password:
        :return: 
        """

        # Inform the user
        log.info("Creating a table with the data properties for galaxy '" + galaxy_name + "' ...")

        # Get DustPedia name
        sample = DustPediaSample()
        galaxy_name = sample.get_name(galaxy_name)

        # Login to DataBase
        database = DustPediaDatabase()
        if username is not None:
            if password is None: raise ValueError("Password is not specified")
            database.login(username, password)

        # Create table
        table = SmartTable()
        table.add_column_info("Observatory/telescope", str, None, "Observatory of telescope")
        table.add_column_info("Instrument", str, None, "Instrument")
        table.add_column_info("Band", str, None, "Band")
        table.add_column_info("Wavelength", float, "micron", "Filter effective wavelength")
        table.add_column_info("Pixelscale", float, "arcsec", "Image pixelscale")
        table.add_column_info("FWHM of PSF", float, "arcsec", "FWHM of the PSF")
        table.add_column_info("Calibration uncertainty (%)", float, None, "Percentual calibration uncertainty")
        table.add_column_info("Calibration uncertainty", float, "mag", "Calibration uncertaintity in magnitudes")

        # Loop over the images for M81
        filters = database.get_image_names_and_filters(galaxy_name)
        for name in filters:

            # Get the filter
            fltr = filters[name]

            # Get filter properties
            observatory = fltr.observatory
            instrument = fltr.instrument
            band = fltr.band
            wavelength = fltr.wavelength

            if observatory == "SDSS": observatory = "APO"

            # Get properties
            pixelscale = self.data[fltr].pixelscale
            fwhm = self.data[fltr].fwhm
            percentual_calibration = self.data[fltr].calibration

            # Get magnitude calibration
            calibration = calibration_magnitude_for_filter(fltr)

            # Set row values
            values = [observatory, instrument, band, wavelength, pixelscale, fwhm, percentual_calibration, calibration]

            # Add entry to the table
            table.add_row(values)

        # Sort the table
        table.sort("Wavelength")

        # Return the table
        return table

    # -----------------------------------------------------------------

    def get_fwhm(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        fwhm = self.data[fltr].fwhm
        return fwhm

    # -----------------------------------------------------------------

    @property
    def fwhms(self):

        """
        This function ...
        :return: 
        """

        fwhms = dict()
        for fltr in self.data:
            fwhms[fltr] = self.data[fltr].fwhm
        return fwhms

    # -----------------------------------------------------------------

    def get_calibration_error_percentual(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        percentual_calibration = self.data[fltr].calibration
        return percentual_calibration

    # -----------------------------------------------------------------

    def get_calibration_error_relative(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        percentual_calibration = self.data[fltr].calibration
        return 0.01 * percentual_calibration

    # -----------------------------------------------------------------

    def get_calibration_error_reference(self, fltr):

        """
        Thisf ucntion ...
        :param fltr: 
        :return: 
        """

        instrument = fltr.instrument
        if instrument not in reference_for_instrument: return None
        letter = reference_for_instrument[instrument]
        return calibration_references[letter]

    # -----------------------------------------------------------------

    def get_calibration_error_magnitude(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        # Get magnitude calibration
        calibration = calibration_magnitude_for_filter(fltr)
        return calibration

    # -----------------------------------------------------------------

    def has_calibration_error_magnitude(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        magn = calibration_magnitude_for_filter(fltr)
        return magn is not None

    # -----------------------------------------------------------------

    def get_pixelsize(self, fltr):

        """
        This function ...
        :param fltr
        :return: 
        """

        pixelsize = self.data[fltr].pixelscale
        return pixelsize

    # -----------------------------------------------------------------

    def get_pixelscale(self, fltr):

        """
        This function ...
        :param fltr
        :return: 
        """

        pixelsize = self.get_pixelsize(fltr)
        return Pixelscale(pixelsize, pixelsize)

# -----------------------------------------------------------------

calibration_magnitude = dict()
calibration_magnitude["GALEX FUV"] = 0.05
calibration_magnitude["GALEX NUV"] = 0.03
calibration_magnitude["2MASS J"] = 0.03
calibration_magnitude["2MASS H"] = 0.03
calibration_magnitude["2MASS Ks"] = 0.03

# -----------------------------------------------------------------

def calibration_magnitude_for_filter(fltr):

    """
    This function ...
    :param fltr: 
    :return: 
    """

    for key in calibration_magnitude:
        key_filter = parse_filter(key)
        if key_filter == fltr: return calibration_magnitude[key] * u("mag")
    return None

# -----------------------------------------------------------------