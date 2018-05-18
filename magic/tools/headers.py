#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.headers Contains functions for extracting information from FITS headers.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import warnings
import numpy as np

# Import astronomical modules
from astropy import coordinates
from astropy.io.fits import getheader

# Import the relevant PTS classes and modules
from ...core.filter.filter import parse_filter
from ..basics.coordinatesystem import CoordinateSystem
from ...core.basics.log import log
from ..basics.pixelscale import Pixelscale, angular_or_physical_pixelscale
from ...core.units.unit import PhotometricUnit
from ...core.units.parsing import parse_unit as u
from ...core.units.utils import interpret_physical_type
from ...core.tools import types

# -----------------------------------------------------------------

def get_header(path):

    """
    Thisf unction ...
    :param path:
    :return:
    """

    return getheader(path)

# -----------------------------------------------------------------

def flattened(header):

    """
    This function ...
    :param header:
    :return:
    """

    flat_header = copy.deepcopy(header)
    flat_header["NAXIS"] = 2
    if "NAXIS3" in flat_header: del flat_header["NAXIS3"]
    for key in flat_header:
        if "PLANE" in key: del flat_header[key]

    return flat_header

# -----------------------------------------------------------------

def get_pixelscale(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Search for different keywords indicating the pixelscale
    for keyword in ("PIXSCALE", "SECPIX", "PFOV", "PLTSCALE"):

        if keyword in header:

            # Get the scale and try to get the unit
            scale = header[keyword]
            if scale == "N/A": continue

            scale = float(scale)

            #print(scale, type(scale))
            # Composite unit
            #if isinstance(scale, CompositeUnit):
            #    scale = parse_quantity(str(scale))

            # If unit is not defined
            if not hasattr(scale, "unit"):

                try: unit = header.comments[keyword].split("[")[1].split("]")[0]
                except IndexError: unit = None

                # Parse the unit with Astropy
                if unit is not None:

                    unit = unit.replace("asec", "arcsec")
                    if not (unit.endswith("pixel") or unit.endswith("pix")): unit = unit + "/pix"
                    try: unit = u(unit)
                    except ValueError: unit = None

                # Debugging
                log.debug("pixelscale found in " + str(keyword) + " keyword = " + str(scale))
                log.debug("unit for the pixelscale = " + str(unit))

                # If no unit is found, guess that it's arcseconds / pixel ...
                if unit is None:
                    warnings.warn("Unit of the pixelscale is not defined in the header, assuming arcsec/pix ...")
                    unit = u("arcsec/pix")
                scale = scale * unit

            # Return the scale
            return angular_or_physical_pixelscale(scale)

    # Search for the 'PXSCAL1' and 'PXSCAL2' keywords
    for keyword_combination in (("PXSCAL1", "PXSCAL2"), ("XPIXSIZE", "YPIXSIZE")):

        if keyword_combination[0] in header and keyword_combination[1] in header:

            scale1 = float(header[keyword_combination[0]])
            if scale1 == "N/A": continue
            try: unit1 = header.comments[keyword_combination[0]].split("[")[1].split("]")[0]
            except IndexError: unit1 = None

            scale2 = float(header[keyword_combination[1]])
            if scale2 == "N/A": continue
            try: unit2 = header.comments[keyword_combination[1]].split("[")[1].split("]")[0]
            except IndexError: unit2 = None

            # Parse the unit with Astropy
            if unit1 is not None:

                unit1 = unit1.replace("asec", "arcsec")
                if not (unit1.endswith("pixel") or unit1.endswith("pix")): unit1 = unit1 + "/pix"
                try: unit1 = u(unit1)
                except ValueError: unit1 = None

            if unit2 is not None:

                unit2 = unit2.replace("asec", "arcsec")
                if not (unit2.endswith("pixel") or unit2.endswith("pix")): unit2 = unit2 + "/pix"
                try: unit2 = u(unit2)
                except ValueError: unit2 = None

            # Debugging
            log.debug("pixelscale found in PXSCAL1 and PXSCAL2 keywords = (" + str(scale1) + "," +str(scale2) + ")")
            log.debug("unit for the pixelscale = (" + str(unit1) + "," + str(unit2) + ")")

            # If no unit is found, guess that it's arcseconds / pixel ...
            if unit1 is None:
                warnings.warn("Unit of the pixelscale in axis1 is not defined in the header, assuming arcsec/pix ...")
                unit1 = u("arcsec/pix")
            if unit2 is None:
                warnings.warn("Unit of the pixelscale in axis2 is not defined in the header, assuming arcsec/pix ...")
                unit2 = u("arcsec/pix")

            # Set quantities
            scale1 = scale1 * unit1
            scale2 = scale2 * unit2

            # Return the pixelscale
            return angular_or_physical_pixelscale(scale1, scale2)

    # If none of the above keywords were found, return None
    else: return None

# -----------------------------------------------------------------

def get_psf_filter(header):

    """
    This function ...
    :param header: 
    :return: 
    """

    if "PSFFLTR" not in header: return None
    else: return parse_filter(header["PSFFLTR"])

# -----------------------------------------------------------------

def get_smoothing_factor(header):

    """
    This function ...
    :param header:
    :return:
    """

    if "SMOOTHF" not in header: return 1.
    else: return float(header["SMOOTHF"])

# -----------------------------------------------------------------

def get_wavelength(header):

    """
    This function ...
    :param header:
    :return:
    """

    # # Get the wavelength
    if "WAVELEN" in header: wavelength = get_quantity(header["WAVELEN"], "micron")
    elif "WVLNGTH" in header: wavelength = get_quantity(header["WVLNGTH"], "micron")
    else: wavelength = None
    return wavelength

# -----------------------------------------------------------------

def get_filter(name, header=None):

    """
    This function ...
    :param name:
    :param header:
    :return:
    """

    filterid = name.lower()
    channel = None
    wavelength = None
    frequency = None

    if "kernel" in filterid:
        log.debug("The image represents a kernel, so no filter will be set")
        return None

    # Get information from the header
    if header is not None:

        # HEADER EXPLICITLY SAYS THAT THERE IS NO FILTER!
        if "FILTER" in header and header["FILTER"] == "n/a": return None

        # Get information regarding the telescope and instrument
        if "TELESCOP" in header: filterid += " " + get_string(header["TELESCOP"]).lower()
        if "INSTRUME" in header: filterid += " " + get_string(header["INSTRUME"]).lower()
        if "ORIGIN" in header: filterid += " " + get_string(header["ORIGIN"]).lower()
        if "OBSERVAT" in header: filterid += " " + get_string(header["OBSERVAT"]).lower()

        # Look for explicit telescope or instrument specifications
        if "OBSERVAT" in header:

            value = header["OBSERVAT"]
            if "CTIO" in value:

                # CTIO FILTER
                if "FILTER" in header: filter_name = header["FILTER"]
                else: filter_name = ""

                # Search for FNAMES
                fltrname = get_filter_name_from_fnamei_keys(header)
                if fltrname is not None: filter_name += " " + fltrname

                # Create filterstring
                filterstring = "CTIO " + filter_name
                return parse_filter(filterstring)

        if "TELESCOP" in header:

            value = header["TELESCOP"]
            if "CTIO" in value:

                # CTIO FILTER
                if "FILTER" in header: filter_name = header["FILTER"]
                else: filter_name = ""

                # SEARCH FOR FNAMES
                fltrname = get_filter_name_from_fnamei_keys(header)
                if fltrname is not None: filter_name += " " + fltrname

                # Create filterstring
                filterstring = "CTIO " + filter_name
                return parse_filter(filterstring)

        fltr = get_filter_from_filteri_keys(header)
        if fltr is not None: return fltr

        # FILTERS
        if "FILTERS" in header:
            try:
                #print(header["FILTERS"])
                filter = parse_filter(header["FILTERS"].upper())
                return filter
            except ValueError: pass

        # Get a name describing the filter
        if "FILTER" in header:
            # One-letter filters are too ambiguous!
            if len(header["FILTER"]) > 1:
                try:
                    filter = parse_filter(header["FILTER"].upper())
                    return filter
                except ValueError: pass
            filterid += " " + get_string(header['FILTER']).lower()
        if "FLTRNM" in header:
            try:
                filter = parse_filter(header["FLTRNM"].upper())
                return filter
            except ValueError: pass
            filterid += " " + get_string(header['FLTRNM']).lower()

        # Get information about the channel number
        if "CHNLNUM" in header: channel = get_int(header["CHNLNUM"])
        elif "BAND" in header: channel = get_int(header["BAND"])
        else: channel = None

        # Get the wavelength
        wavelength = get_wavelength(header)

        # Get the frequency
        if "FREQ" in header: frequency = get_quantity(header["FREQ"], "GHz")
        elif "FREQUENCY" in header: frequency = get_quantity(header["FREQUENCY"], "GHz")
        else: frequency = None

    # Debug information
    log.debug("filterid = " + str(filterid))
    log.debug("channel = " + str(channel))
    log.debug("wavelength = " + str(wavelength))
    log.debug("frequency = " + str(frequency))

    final_filter_name = None

    # -- UV --

    # GALEX
    if "fuv" in filterid: final_filter_name = "GALEX FUV"
    elif "nuv" in filterid: final_filter_name = "GALEX NUV"

    # SWIFT
    elif "uw2" in filterid: final_filter_name = "SWIFT W2"
    elif "um2" in filterid: final_filter_name = "SWIFT M2"
    elif "uw1" in filterid: final_filter_name = "SWIFT W1"
    elif "swift" in filterid or "uvot" in filterid:

        if "w2" in filterid: final_filter_name = "SWIFT W2"
        elif "m2" in filterid: final_filter_name = "SWIFT M2"
        elif "w1" in filterid: final_filter_name = "SWIFT W1"
        else: log.warning("Could not determine which SWIFT UVOT filter was used for this image")

    # TODO: support other UV instruments

    # -- Optical --

    # SDSS
    elif "sdss" in filterid:

        if "-u" in filterid: final_filter_name = "SDSS u"
        elif "-g" in filterid: final_filter_name = "SDSS g"
        elif "-r" in filterid: final_filter_name = "SDSS r"
        elif "-i" in filterid: final_filter_name = "SDSS i"
        elif "-z" in filterid: final_filter_name = "SDSS z"
        else:

            if "sdss u" in filterid: final_filter_name = "SDSS u"
            elif "sdss g" in filterid: final_filter_name = "SDSS g"
            elif "sdss r" in filterid: final_filter_name = "SDSS r"
            elif "sdss i" in filterid: final_filter_name = "SDSS i"
            elif "sdss z" in filterid: final_filter_name = "SDSS z"
            else:

                if "u" in name: final_filter_name = "SDSS u"
                elif "g" in name: final_filter_name = "SDSS g"
                elif "r" in name: final_filter_name = "SDSS r"
                elif "i" in name: final_filter_name = "SDSS i"
                elif "z" in name: final_filter_name = "SDSS z"
                else: log.warning("Could not determine which SDSS filter was used for this image")

    # R band // not good; H alpha image was also identified as R band ...
    #elif "r" in filterid and "kpno" in filterid: return Filter("KPNO.Mosaic.R")

    # TODO: support other optical instruments

    # -- IR --

    # 2MASS filters
    elif "2mass" in filterid:

        if "h" in filterid: final_filter_name = "2MASS H"
        elif "j" in filterid: final_filter_name = "2MASS J"
        elif "k" in filterid: final_filter_name = "2MASS Ks"
        else: log.warning("Could not determine which 2MASS filter was used for this image")

    # UKIDSS/UKIRT filters
    elif "ukidss" in filterid:
    
        if "h" in filterid: final_filter_name = "UKIDSS H"
        elif "j" in filterid: final_filter_name = "UKIDSS J"
        elif "k" in filterid: final_filter_name = "UKIDSS K"
        else: log.warning("Could not determine which UKIDSS filter was used for this image")
    


    # IRAC filters
    elif "irac" in filterid:

        if "3.6" in filterid or "i1" in filterid: final_filter_name = "IRAC I1"
        elif "4.5" in filterid or "i2" in filterid: final_filter_name = "IRAC I2"
        elif "5.8" in filterid or "i3" in filterid: final_filter_name = "IRAC I3"
        elif "8.0" in filterid or "i4" in filterid: final_filter_name = "IRAC I4"
        else:  # Look at the channel number

            if channel is not None:

                if channel == 1: final_filter_name = "IRAC I1"
                elif channel == 2: final_filter_name = "IRAC I2"
                elif channel == 3: final_filter_name = "IRAC I3"
                elif channel == 4: final_filter_name = "IRAC I4"
                else: log.warning("Could not determine which IRAC filter was used for this image")

            elif wavelength is not None:

                if np.isclose(wavelength.to("micron").value, 3.6, rtol=0.05): final_filter_name = "IRAC I1"
                elif np.isclose(wavelength.to("micron").value, 4.5, rtol=0.05): final_filter_name = "IRAC I2"
                elif np.isclose(wavelength.to("micron").value, 5.8, rtol=0.05): final_filter_name = "IRAC I3"
                elif np.isclose(wavelength.to("micron").value, 8.0, rtol=0.05): final_filter_name = "IRAC I4"
                else: log.warning("Could not determine which IRAC filter was used for this image")

            else: log.warning("Could not determine which IRAC filter was used for this image")

    # WISE filters
    elif "wise" in filterid:

        if "w1" in filterid: final_filter_name = "WISE W1"
        elif "w2" in filterid: final_filter_name = "WISE W2"
        elif "w3" in filterid: final_filter_name = "WISE W3"
        elif "w4" in filterid: final_filter_name = "WISE W4"
        else:

            if channel is not None:

                if channel == 1: final_filter_name = "WISE W1"
                elif channel == 2: final_filter_name = "WISE W2"
                elif channel == 3: final_filter_name = "WISE W3"
                elif channel == 4: final_filter_name = "WISE W4"
                else: log.warning("Could not determine which WISE filter was used for this image")

            elif wavelength is not None:

                if np.isclose(wavelength.to("micron").value, 3.4, rtol=0.05): final_filter_name = "WISE W1"
                elif np.isclose(wavelength.to("micron").value, 4.6, rtol=0.05): final_filter_name = "WISE W2"
                elif np.isclose(wavelength.to("micron").value, 12., rtol=0.05): final_filter_name = "WISE W3"
                elif np.isclose(wavelength.to("micron").value, 22., rtol=0.05): final_filter_name = "WISE W4"
                else: log.warning("Could not determine which WISE filter was used for this image")

            elif "3.4" in filterid: final_filter_name = "WISE W1"
            elif "4.6" in filterid: final_filter_name = "WISE W2"
            elif "11.6" in filterid or "12" in filterid: final_filter_name = "WISE W3"
            elif "22.1" in filterid or "22" in filterid: final_filter_name = "WISE W4"
            else: log.warning("Could not determine which WISE filter was used for this image")

    # MIPS filters
    elif "mips" in filterid:

        if "24" in filterid: final_filter_name = "MIPS 24"
        elif "70" in filterid: final_filter_name = "MIPS 70"
        elif "160" in filterid: final_filter_name = "MIPS 160"
        else:

            if wavelength is not None:

                if np.isclose(wavelength.to("micron").value, 24., rtol=0.05): final_filter_name = "MIPS 24"
                elif np.isclose(wavelength.to("micron").value, 70., rtol=0.05): final_filter_name = "MIPS 70"
                elif np.isclose(wavelength.to("micron").value, 160., rtol=0.05): final_filter_name = "MIPS 160"
                else: log.warning("Could not determine which MIPS filter was used for this image")

            else: log.warning("Could not determine which MIPS filter was used for this image")

    # Spitzer bands (IRAC and MIPS but "irac" and "mips" are not in filterid)
    elif "spitzer" in filterid:

        # Look for IRAC wavelength
        if "3.6" in filterid or "i1" in filterid: final_filter_name = "IRAC I1"
        elif "4.5" in filterid or "i2" in filterid: final_filter_name = "IRAC I2"
        elif "5.8" in filterid or "i3" in filterid: final_filter_name = "IRAC I3"
        elif "8.0" in filterid or "i4" in filterid: final_filter_name = "IRAC I4"

        # Look for MIPS wavelengths
        elif "24" in filterid: final_filter_name = "MIPS 24"
        elif "70" in filterid: final_filter_name = "MIPS 70"
        elif "160" in filterid: final_filter_name = "MIPS 160"

        elif wavelength is not None:

            if np.isclose(wavelength.to("micron").value, 3.6, rtol=0.05): final_filter_name = "IRAC I1"
            elif np.isclose(wavelength.to("micron").value, 4.5, rtol=0.05): final_filter_name = "IRAC I2"
            elif np.isclose(wavelength.to("micron").value, 5.8, rtol=0.05): final_filter_name = "IRAC I3"
            elif np.isclose(wavelength.to("micron").value, 8.0, rtol=0.05): final_filter_name = "IRAC I4"
            elif np.isclose(wavelength.to("micron").value, 24., rtol=0.05): final_filter_name = "MIPS 24"
            elif np.isclose(wavelength.to("micron").value, 70., rtol=0.05): final_filter_name = "MIPS 70"
            elif np.isclose(wavelength.to("micron").value, 160., rtol=0.05): final_filter_name = "MIPS 160"
            else: log.warning("Could not determine which Spitzer filter was used for this image")

        else: log.warning("Could not determine which Spitzer filter was used for this image")

    # PACS filters
    elif "pacs" in filterid:

        if '70' in filterid or 'blue' in filterid: final_filter_name = "Pacs blue"
        elif '100' in filterid or 'green' in filterid: final_filter_name = "Pacs green"
        elif '160' in filterid or 'red' in filterid: final_filter_name = "Pacs red"
        else: log.warning("Could not determine which PACS filter was used for this image")

    # SPIRE filters
    elif "spire" in filterid:

        if "psw" in filterid or "250" in filterid: final_filter_name = "SPIRE PSW"
        elif "pmw" in filterid or "350" in filterid: final_filter_name = "SPIRE PMW"
        elif "plw" in filterid or "500" in filterid: final_filter_name = "SPIRE PLW"
        else:

            if channel is not None:

                if channel == 1: final_filter_name = "SPIRE PSW"
                elif channel == 2: final_filter_name = "SPIRE PMW"
                elif channel == 3: final_filter_name = "SPIRE PLW"
                else: log.warning("Could not determine which SPIRE filter was used for this image")

            elif wavelength is not None:

                if np.isclose(wavelength.to("micron").value, 250., rtol=0.05): final_filter_name = "SPIRE PSW"
                elif np.isclose(wavelength.to("micron").value, 350, rtol=0.05): final_filter_name = "SPIRE PMW"
                elif np.isclose(wavelength.to("micron").value, 500., rtol=0.05): final_filter_name = "SPIRE PLW"
                else: log.warning("Could not determine which SPIRE filter was used for this image")

    # -- H alpha --
    #elif "alpha" in filterid or "6561" in filterid or "656_1" in filterid: final_filter_name = "656_1"
    elif "alpha" in filterid or "6561" in filterid or "656_1" in filterid: final_filter_name = "Halpha"
    elif "ha" in filterid and "kpno" in filterid: final_filter_name = "Halpha"

    # Planck
    elif "planck" in filterid:

        if "350" in filterid: final_filter_name = "Planck 350"
        elif "550" in filterid: final_filter_name = "Planck 550"
        elif "850" in filterid: final_filter_name = "Planck 850"
        elif "1380" in filterid: final_filter_name = "Planck 1380"
        elif "2100" in filterid: final_filter_name = "Planck 2100"
        elif "3000" in filterid: final_filter_name = "Planck 3000"
        elif "4260" in filterid: final_filter_name = "Planck 4260"
        elif "6810" in filterid: final_filter_name = "Planck 6810"
        elif "10600" in filterid: final_filter_name = "Planck 10600"
        else:

            if "30" in filterid: final_filter_name = "Planck 30"
            elif "44" in filterid: final_filter_name = "Planck 44"
            elif "70" in filterid: final_filter_name = "Planck 70"
            elif "100" in filterid: final_filter_name = "Planck 100"
            elif "143" in filterid: final_filter_name = "Planck 143"
            elif "217" in filterid: final_filter_name = "Planck 217"
            elif "353" in filterid: final_filter_name = "Planck 353"
            elif "545" in filterid: final_filter_name = "Planck 545"
            elif "857" in filterid: final_filter_name = "Planck 857"
            else:

                # Wavelength
                if wavelength is not None:

                    if np.isclose(wavelength.to("micron").value, 350., rtol=0.05): final_filter_name = "Planck 350"
                    elif np.isclose(wavelength.to("micron").value, 550., rtol=0.05): final_filter_name = "Planck 550"
                    elif np.isclose(wavelength.to("micron").value, 850., rtol=0.05): final_filter_name = "Planck 850"
                    elif np.isclose(wavelength.to("micron").value, 1380., rtol=0.05): final_filter_name = "Planck 1380"
                    elif np.isclose(wavelength.to("micron").value, 2100., rtol=0.05): final_filter_name = "Planck 2100"
                    elif np.isclose(wavelength.to("micron").value, 3000., rtol=0.05): final_filter_name = "Planck 3000"
                    elif np.isclose(wavelength.to("micron").value, 4260., rtol=0.05): final_filter_name = "Planck 4260"
                    elif np.isclose(wavelength.to("micron").value, 6810, rtol=0.05): final_filter_name = "Planck 6810"
                    elif np.isclose(wavelength.to("micron").value, 10600, rtol=0.05): final_filter_name = "Planck 10600"
                    else: log.warning("Could not determine which Planck filter was used for this image")

                # Frequency
                elif frequency is not None:

                    if np.isclose(frequency.to("GHz").value, 30., rtol=0.05): final_filter_name = "Planck 30"
                    elif np.isclose(frequency.to("GHz").value, 44., rtol=0.05): final_filter_name = "Planck 44"
                    elif np.isclose(frequency.to("GHz").value, 70., rtol=0.05): final_filter_name = "Planck 70"
                    elif np.isclose(frequency.to("GHz").value, 100, rtol=0.05): final_filter_name = "Planck 100"
                    elif np.isclose(frequency.to("GHz").value, 143, rtol=0.05): final_filter_name = "Planck 143"
                    elif np.isclose(frequency.to("GHz").value, 217, rtol=0.05): final_filter_name = "Planck 217"
                    elif np.isclose(frequency.to("GHz").value, 353, rtol=0.05): final_filter_name = "Planck 353"
                    elif np.isclose(frequency.to("GHz").value, 545, rtol=0.05): final_filter_name = "Planck 545"
                    elif np.isclose(frequency.to("GHz").value, 857, rtol=0.05): final_filter_name = "Planck 857"
                    else: log.warning("Could not determine which Planck filter was used for this image")

    # No filter name could be composed
    if final_filter_name is None:

        # If wavelength was found
        if wavelength is not None:

            # Warning
            log.warning("Filter could not be identified, but wavelength is " + str(wavelength))

            # Set limits
            value = wavelength.to("micron").value
            five_percent = 0.05 * value
            lower = value - five_percent
            upper = value + five_percent

            log.warning("Setting the filter to a uniform filter between the wavelengths " + str(lower) + " micron and " + str(upper) + " micron")

            if "FILTER" in header: name = header["FILTER"].split("  / ")[0].replace(" ", "")
            else: name = None

            # Create a custom filter around the wavelength
            fltr = parse_filter((lower, upper), name=name)

        # Wavelength could not be found
        else:
            log.warning("Filter or wavelength could not be identified")
            fltr = None

    # Filter name could be composed
    else:

        # Create the filter
        fltr = parse_filter(final_filter_name)

        # Inform the user
        log.debug("Filter was identified as " + str(fltr))

    # Create and return a Filter object
    return fltr

# -----------------------------------------------------------------

def get_unit(header, density=False, brightness=False, density_strict=False, brightness_strict=False):

    """
    This function ...
    :param header:
    :param density:
    :param brightness:
    :param density_strict:
    :param brightness_strict:
    :return:
    """

    # Check whether physical type is defined
    if "PHYSTYPE" in header:

        physical_type = header["PHYSTYPE"]
        density_interpreted, brightness_interpreted = interpret_physical_type(physical_type)

        # CHECK DENSITY
        if density_strict and density_interpreted != density:

            if density: raise ValueError("Density_strict is True and density is True but found density = False from header physical type specification")
            else: raise ValueError("Density_strict is True and density is False but found density = True from header physical type specification")

        # DENSITY OK
        density = density_interpreted

        # CHECK BRIGHTNESS
        if brightness_strict and brightness_interpreted != brightness:

            if brightness: raise ValueError("Brightness_strict is True and brightness is True but found brightness = False from header physical type specification")
            else: raise ValueError("Brightness_strict is True and brightness is False but found brightness = True from header physical type specification")

        # BRIGHTNESS OK
        brightness = brightness_interpreted

        # Set strict, because we found specification in header
        density_strict = brightness_strict = True

    #else: density = brightness = density_strict = brightness_strict = False # these are the defaults

    unit = None

    # Look for different keywords
    for keyword in ("BUNIT", "SIGUNIT", "ZUNITS"):

        # Skip
        if keyword not in header: continue

        value = header[keyword].split("   / ")[0].rstrip()

        try:
            # Unit found
            unit = PhotometricUnit(value, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict)
            break

        # Parsing failed
        except ValueError: continue

    # Return the unit
    return unit

# -----------------------------------------------------------------

def get_fwhm(header):

    """
    This function ...
    :param header:
    :return:
    """

    fwhm = None

    for keyword in ["FWHM"]:

        if keyword not in header: continue

        # Get the FWHM
        fwhm = get_quantity(header["FWHM"], default_unit="arcsec")

    # Return the FWHM
    return fwhm

# -----------------------------------------------------------------

def get_distance(header):
    
    """
    This function ...
    :param header: 
    :return: 
    """

    distance = None

    for keyword in ["DISTANCE"]:

        if keyword not in header: continue

        # Get the distance
        distance = get_quantity(header["DISTANCE"], default_unit="Mpc")

    # Return the distance
    return distance

# -----------------------------------------------------------------

def get_zero_point(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Loop over all keys in the header
    for key in header:

        if "MAGZP" in key: return header[key]

    # If no keyword is found that states the zero-point, return None
    return None

# -----------------------------------------------------------------

def is_sky_subtracted(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Initially, set the boolean to False
    subtracted = False

    # Find the keyword
    if 'BACK_SUB' in header: subtracted = header['BACK_SUB']

    # Check type
    if not isinstance(subtracted, bool): raise ValueError("Invalid value for 'BACK_SUB'")

    # Return the boolean value
    return subtracted

# -----------------------------------------------------------------

def is_source_extracted(header):

    """
    This function ...
    :param header:
    :return:
    """

    extracted = False

    # Find the keyword
    if "SRC_EXTR" in header: extracted = header["SRC_EXTR"]

    # Check type
    if not isinstance(extracted, bool): raise ValueError("Invalid value for 'SRC_EXTR'")

    # Return the boolean value
    return extracted

# -----------------------------------------------------------------

def is_extinction_corrected(header):

    """
    This function ...
    :param header:
    :return:
    """

    corrected = False

    # Find the keyword
    if "EXT_CORR" in header: corrected = header["EXT_CORR"]

    # Check type
    if not isinstance(corrected, bool): raise ValueError("Invalid value for 'EXT_CORR'")

    # Return the boolean value
    return corrected

# -----------------------------------------------------------------

def get_number_of_frames(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Initially, set the boolean to False
    nframes = 1

    if 'NAXIS' in header:

        # If there are 3 axes, get the size of the third
        if header['NAXIS'] == 3: nframes = header['NAXIS3']

    # Return the boolean value
    return nframes

# -----------------------------------------------------------------

def get_frame_name_and_description(header, i, always_call_first_primary=True, absolute_index_names=True, frame_names=None,
                                   mask_names=None, segments_names=None):

    """
    This function ...
    :param header:
    :param i:
    :param always_call_first_primary:
    :param absolute_index_names:
    :param frame_names:
    :param mask_names:
    :param segments_names:
    :return:
    """

    planeX = "PLANE" + str(i)

    # Return the description
    if planeX in header: description = header[planeX]
    else: description = None

    plane_type = "frame"

    # FITS file created by AstroMagic
    if description is not None and "[" in description and "]" in description:

        name = description.split(" [")[0]
        plane_type = description.split("[")[1].split("]")[0]

    elif i == 0 and always_call_first_primary:

        # Get the name of this frame, but the first frame always gets the name 'primary' unless the
        # 'always_call_first_primary' flag is disabled

        description = "the primary signal map"
        name = "primary"

    elif description is not None:

        # Convert spaces to underscores and ignore things between parentheses
        name = description.split("(")[0].rstrip(" ").replace(" ", "_")

        # If the frame name contains 'error', use the standard name "errors" for this frame
        if 'error' in name: name = "errors"

    else: ## description is None

        description = str(i) + "th plane of the image"

        if absolute_index_names: name = "frame" + str(i)
        else:
            # Determine the next index for the frame based on the current frame indices
            from ...core.tools import numbers
            if frame_names is None: raise ValueError("Frame names must be defined")
            current_indices = [int(frame_name.split("frame")[1]) for frame_name in frame_names if frame_name.startswith("frame")]
            relative_index = numbers.lowest_missing_integer(current_indices)
            name = "frame" + str(relative_index)

    # Return the name and description
    return name, description, plane_type

# -----------------------------------------------------------------

def get_frame_description(header, i):

    """
    This function ...
    :param header:
    :param i:
    :return:
    """

    planeX = "PLANE" + str(i)

    # Return the description
    if planeX in header: return header[planeX]
    else: return None

# -----------------------------------------------------------------

def get_frame_index(header, name):

    """
    This function ...
    """

    for key in header:

        # Skip keys not ...
        if not "PLANE" in key: continue

        if header[key] == name: return int(key.split("PLANE")[1])

    return None

# -----------------------------------------------------------------

def get_frame_name(description):

    """
    This function ...
    :param description:
    :return:
    """

    # Convert spaces to underscores and ignore things between parentheses
    name = description.split("(")[0].rstrip(" ").replace(" ", "_")

    # If the frame name contains 'error', use the standard name "errors" for this frame
    if 'error' in name: name = "errors"

    # Return the frame name
    return name

# -----------------------------------------------------------------

def check_header_matches_image(image, header):

    """
    This function ...
    :param image:
    :param header:
    :return:
    """

    wcs = CoordinateSystem(header)

    # wcs.naxis attributes are deprecated, so we perform this check conditionally
    if ((hasattr(wcs,'naxis1') and hasattr(wcs,'naxis2')) and not
            (wcs.naxis1 == image.shape[1] and wcs.naxis2 == image.shape[0])):
        raise Exception("Image shape must match header shape.")

# -----------------------------------------------------------------

def get_pixel_mapping(header1, header2):

    """
    This function determines the mapping from pixel coordinates in header1 to pixel coordinates in header2 (the
    reference header). It takes the following arguments:
    :param header1:
    :param header2:
    :return: a NumPy array describing a grid of y,x pixel locations in the input header's pixel units but the output
    header's world units. It raises a TypeError if neither header is not a Header or WCS instance, and a
    NotImplementedError if the CTYPE in the header is not recognized.
    """

    # Get the WCS from the two headers
    wcs1 = CoordinateSystem(header1)
    wcs2 = CoordinateSystem(header2)

    # Convert the coordinates
    if not all([w1 == w2 for w1,w2 in zip(wcs1.wcs.ctype,wcs2.wcs.ctype)]):
        allowed_coords = ('GLON', 'GLAT', 'RA', 'DEC')
        if all([(any(word in w1 for word in allowed_coords) and
                 any(word in w2 for word in allowed_coords))
                for w1,w2 in zip(wcs1.wcs.ctype,wcs2.wcs.ctype)]):
            csys1 = ctype_to_csys(wcs1.wcs)
            csys2 = ctype_to_csys(wcs2.wcs)
            convert_coordinates = True
        else:
            # do unit conversions
            raise NotImplementedError("Unit conversions between {0} and {1} have not yet been implemented.".format(wcs1.wcs.ctype,wcs2.wcs.ctype))

    else: convert_coordinates = False

    # sigh... why does numpy use matrix convention?  Makes everything so
    # much harder...
    # WCS has naxis attributes because it is loaded with
    # _load_wcs_from_header
    outshape = [wcs2.naxis2,wcs2.naxis1]

    yy2,xx2 = np.indices(outshape)

    # get the world coordinates of the output image
    lon2,lat2 = wcs2.wcs_pix2world(xx2, yy2, 0)

    # Alternative
    #x = np.arange(wcs2.naxis1)
    #y = np.arange(wcs2.naxis2)
    #X, Y = np.meshgrid(x, y)
    #lon2, lat2 = wcs2.wcs_pix2world(X, Y, 0)

    if convert_coordinates:

        # Transform the world coordinates from the output image into the coordinate
        # system of the input image
        C2 = coordinates.SkyCoord(lon2, lat2, unit="deg", frame=csys2)
        C1 = C2.transform_to(csys1)
        lon2,lat2 = C1.spherical.lon.deg,C1.spherical.lat.deg

    xx1,yy1 = wcs1.wcs_world2pix(lon2, lat2, 0)
    grid = np.array([yy1.reshape(outshape),xx1.reshape(outshape)])

    # Return the grid
    return grid

# -----------------------------------------------------------------

def ctype_to_csys(wcs):

    """
    This function ...
    :param wcs:
    :return:
    """

    ctype = wcs.ctype[0]
    if 'RA' in ctype or 'DEC' in ctype:
        if wcs.equinox == 2000 or wcs.equinox == 2000.:
            return 'fk5'
        elif wcs.equinox == 1950 or wcs.equinox == 1950.:
            return 'fk4'
        else:
            raise NotImplementedError("Non-fk4/fk5 equinoxes are not allowed")
    elif 'GLON' in ctype or 'GLAT' in ctype:
        return 'galactic'

# -----------------------------------------------------------------

def get_float(entry):

    """
    This function ...
    :param entry:
    :return:
    """

    try: return float(entry)
    except ValueError:
        value = entry.split("   / ")[0].rstrip()
        return float(value)

# -----------------------------------------------------------------

def get_quantity(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    if types.is_string_type(entry): value = entry.split("   / ")[0].rstrip()
    else: value = entry

    # Try parsing
    try:

        num_value = float(value)
        if default_unit is None: raise RuntimeError("Default unit is not provided")
        unit = u(default_unit)

    except ValueError:

        #floats = re.findall("[-+]?\d*\.\d+|\d+", value)
        #assert len(floats) == 1
        #num_value = floats[0]

        #unit_description = value.split(str(floats[0])[-1:])[1].rstrip()

        #print("unit:::", unit_description)

        #unit = u.Unit(unit_description)

        composite_unit = u(value)
        if default_unit is None: raise RuntimeError("Default unit is not provided")
        num_value = composite_unit.to(default_unit)
        unit = u(default_unit)

    # Return quantity
    return num_value * unit

# -----------------------------------------------------------------

def get_string(entry):

    """
    This function ...
    :param entry:
    :return:
    """

    value = entry.split("   / ")[0].rstrip()
    return value

# -----------------------------------------------------------------

def get_int(entry):

    """
    This function ...
    :param entry:
    :return:
    """

    try: return int(entry)
    except ValueError:
        value = entry.split("   / ")[0].rstrip()
        return int(value)

# -----------------------------------------------------------------

def get_filter_from_filteri_keys(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Find FILTERI keywords
    for filter_index in range(5):

        filter_index_keyword = "FILTER" + str(filter_index)
        if filter_index_keyword in header:
            try:
                # print(header[filter_index_keyword])
                filter = parse_filter(header[filter_index_keyword].upper())
                return filter
            except ValueError:
                pass

    return None

# -----------------------------------------------------------------

def get_filter_name_from_fnamei_keys(header):

    """
    This function ...
    :param header:
    :return:
    """

    for filter_index in range(5):

        filter_index_keyword = "FNAME" + str(filter_index)
        if filter_index_keyword in header:
            if "diaphragm" in header[filter_index_keyword].lower(): continue # just a diaphragm
            return header[filter_index_keyword]

    return None

# -----------------------------------------------------------------
