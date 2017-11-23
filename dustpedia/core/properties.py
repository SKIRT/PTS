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
import numpy as np
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

def has_fwhm(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    properties = DustPediaProperties()
    return properties.has_filter(fltr)

# -----------------------------------------------------------------

def has_calibration_error(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    properties = DustPediaProperties()
    return properties.has_filter(fltr)

# -----------------------------------------------------------------

def get_calibration_error(fltr, flux, errorbar=False):

    """
    Thisn function ...
    :param fltr:
    :param flux: can be array (Frame) or single value
    :param errorbar:
    :return:
    """

    from ...core.basics.errorbar import ErrorBar

    # Load DustPedia properties
    properties = DustPediaProperties()

    # Calculate calibration error frame
    if properties.has_calibration_error_magnitude(fltr):

        magnitude = properties.get_calibration_error_magnitude(fltr)

        # Real value or quantity
        if types.is_real_type(flux) or types.is_quantity(flux): calibration_error = calculate_calibration_error(flux, magnitude, errorbar=errorbar)

        # Assume it's a frame
        else:
            if errorbar: raise ValueError("Cannot return as errorbar")
            calibration_error = get_calibration_uncertainty_frame_from_magnitude(flux, magnitude)

    # From percentual calibration error
    else:

        # Calculate calibration errors with percentage
        fraction = properties.get_calibration_error_relative(fltr)
        calibration_error = flux * fraction

        # As error bar
        if errorbar: calibration_error = ErrorBar(calibration_error)

    # Return the calibration error
    return calibration_error

# -----------------------------------------------------------------

def calculate_calibration_error(flux, calibration_magn, errorbar=True):

    """
    This function ...
    :param flux:
    :param calibration_magn:
    :param errorbar:
    :return:
    """

    from ...core.basics.errorbar import ErrorBar
    from ...core.units.utils import jansky_to_ab, ab_to_jansky

    if hasattr(calibration_magn, "unit"): mag_error = calibration_magn.to("mag").value
    else: mag_error = float(calibration_magn)

    # Inform the user
    #log.info("Calculating the calibration error for the " + name + " image ...")

    # Get the calibration error
    #calibration_error = CalibrationError.from_filter(self.frames[name].filter)
    #mag_error = calibration_error.value

    # Convert the total flux to AB magnitude
    flux_mag = jansky_to_ab(flux)

    a = flux_mag - mag_error
    b = flux_mag + mag_error

    #log.debug("Flux value: " + str(flux))
    #log.debug("a magnitude: " + str(a))
    #log.debug("b magnitude: " + str(b))

    # Convert a and b to Jy
    a = ab_to_jansky(a)
    b = ab_to_jansky(b)

    # Calculate lower and upper limit of the flux
    c = a - flux
    d = b - flux

    #log.debug("a value: " + str(a))
    #log.debug("b value: " + str(b))
    #log.debug("c value: " + str(c))
    #log.debug("d value: " + str(d))

    # Create the error bar
    if errorbar: error = ErrorBar(d, c)

    # No errorbar: max value
    else: error = max(abs(d), abs(c))

    # Return the calibration error
    return error

# -----------------------------------------------------------------

def get_calibration_uncertainty_frame_from_magnitude(image, calibration_magn):

    """
    This fucntion ...
    :param image:
    :param calibration_magn:
    :return:
    """

    from ...magic.core.frame import Frame
    from ...magic.basics.mask import Mask as oldMask
    from ...core.units.utils import jansky_to_ab, ab_mag_zero_point

    # Convert the frame into AB magnitudes
    invalid = oldMask.is_zero_or_less(image)
    ab_frame = jansky_to_ab(image)
    # Set infinites to zero
    ab_frame[invalid] = 0.0

    # -----------------------------------------------------------------

    # The calibration uncertainty in AB magnitude
    #mag_error = calibration_magn.value
    if hasattr(calibration_magn, "unit"): mag_error = calibration_magn.to("mag").value
    else: mag_error = float(calibration_magn)

    # a = image[mag] - mag_error
    a = ab_frame - mag_error

    # b = image[mag] + mag_error
    b = ab_frame + mag_error

    # Convert a and b to Jy
    a = ab_mag_zero_point.to("Jy").value * np.power(10.0, -2. / 5. * a)
    b = ab_mag_zero_point.to("Jy").value * np.power(10.0, -2. / 5. * b)

    # c = a[Jy] - image[Jy]
    # c = a - jansky_frame
    c = a - image

    # d = image[Jy] - b[Jy]
    # d = jansky_frame - b
    d = image - b

    # ----------------------------------------------------------------- BELOW: if frame was not already in Jy

    # Calibration errors = max(c, d)
    # calibration_errors = np.maximum(c, d)  # element-wise maxima
    # calibration_errors[invalid] = 0.0 # set zero where AB magnitude could not be calculated

    # relative_calibration_errors = calibration_errors / jansky_frame
    # relative_calibration_errors[invalid] = 0.0

    # Check that there are no infinities or nans in the result
    # assert not np.any(np.isinf(relative_calibration_errors)) and not np.any(np.isnan(relative_calibration_errors))

    # The actual calibration errors in the same unit as the data
    # calibration_frame = self.image.frames.primary * relative_calibration_errors

    # -----------------------------------------------------------------

    calibration_frame = np.maximum(c, d)  # element-wise maxima
    calibration_frame[invalid] = 0.0  # set zero where AB magnitude could not be calculated

    new_invalid = oldMask.is_nan(calibration_frame) + oldMask.is_inf(calibration_frame)
    calibration_frame[new_invalid] = 0.0

    # Check that there are no infinities or nans in the result
    assert not np.any(np.isinf(calibration_frame)) and not np.any(np.isnan(calibration_frame))

    # Make frame from numpy array
    calibration_frame = Frame(calibration_frame)

    # Return the frame
    return calibration_frame

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

        # Loop over the images
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

    def create_tables(self, galaxy_name, username, password):

        """
        This function ...
        :param galaxy_name:
        :param username:
        :param password:
        :return:
        """

        # Inform the user
        log.info("Creating tables with the data properties for galaxy '" + galaxy_name + "' ...")

        # Get DustPedia name
        sample = DustPediaSample()
        galaxy_name = sample.get_name(galaxy_name)

        # Login to DataBase
        database = DustPediaDatabase()
        if username is not None:
            if password is None: raise ValueError("Password is not specified")
            database.login(username, password)

        # Loop over the images
        observatories = database.get_image_names_and_filters_per_observatory(galaxy_name)

        # Dict for tables
        tables = OrderedDict()

        # Loop over the observatories
        for observatory in observatories:

            # Create table
            table = SmartTable()
            table.add_column_info("Band", str, None, "Band")
            table.add_column_info("Wavelength", float, "micron", "Filter effective wavelength")
            table.add_column_info("Pixelscale", float, "arcsec", "Image pixelscale")
            table.add_column_info("FWHM of PSF", float, "arcsec", "FWHM of the PSF")
            table.add_column_info("Calibration uncertainty (%)", float, None, "Percentual calibration uncertainty")
            table.add_column_info("Calibration uncertainty", float, "mag", "Calibration uncertaintity in magnitudes")

            # Loop over the filters
            for name in observatories[observatory]:

                # Get the filter
                fltr = observatories[observatory][name]

                # Get filter properties
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
                values = [band, wavelength, pixelscale, fwhm, percentual_calibration, calibration]

                # Add entry to the table
                table.add_row(values)

            # Sort the table
            table.sort("Wavelength")

            # Add the table
            tables[observatory] = table

        # Return the dictionary of tables
        return tables

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