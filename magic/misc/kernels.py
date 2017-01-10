#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.kernels Contains the AnianoKernels class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import urllib
import requests
from lxml import html
from collections import OrderedDict

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ..core.kernel import ConvolutionKernel
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.tools import archive
from ...core.basics.filter import Filter

# -----------------------------------------------------------------

# The path to the PTS kernels directory
kernels_path = fs.join(introspection.pts_user_dir, "kernels")

# -----------------------------------------------------------------

# The link to the Aniano high-resolution kernel files and the Aniano PSF files
aniano_kernels_highres_link = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2012_May/Kernels_fits_Files/Hi_Resolution/"
aniano_kernels_lowres_link = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2012_May/Kernels_fits_Files/Low_Resolution/"
aniano_psf_files_link = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2012_May/PSF_fits_Files/"

# -----------------------------------------------------------------

# Reference: Common-Resolution Convolution Kernels for Space- and Ground-Based Telescopes (G. Aniano et. al)
fwhms = {"GALEX_FUV": 4.48 * Unit("arcsec"),
         "GALEX_NUV": 5.05 * Unit("arcsec"),
         "IRAC_3.6": 1.90 * Unit("arcsec"),
         "IRAC_4.5": 1.81 * Unit("arcsec"),
         "IRAC_5.8": 2.11 * Unit("arcsec"),
         "IRAC_8.0": 2.82 * Unit("arcsec"),
         "WISE_FRAME_3.4": 5.79 * Unit("arcsec"),
         "WISE_FRAME_4.6": 6.37 * Unit("arcsec"),
         "WISE_FRAME_11.6": 6.60 * Unit("arcsec"),
         "WISE_FRAME_22.1": 11.89 * Unit("arcsec"),
         "MIPS_24": 6.43 * Unit("arcsec"),
         "MIPS_70": 18.74 * Unit("arcsec"),
         "MIPS_160": 38.78 * Unit("arcsec"),
         "PACS_70": 5.67 * Unit("arcsec"),
         "PACS_100": 7.04 * Unit("arcsec"),
         "PACS_160": 11.18 * Unit("arcsec"),
         "SPIRE_250": 18.15 * Unit("arcsec"),
         "SPIRE_350": 24.88 * Unit("arcsec"),
         "SPIRE_500": 36.09 * Unit("arcsec")}

# -----------------------------------------------------------------

variable_fwhms = ["SDSS u", "SDSS g", "SDSS r", "SDSS i", "SDSS z"]

# -----------------------------------------------------------------

aniano_names = {"UVOT UVM2": "Gauss_02.5",
                "UVOT UVW1": "Gauss_02.5",
                "UVOT UVW2": "Gauss_03.0",
                "GALEX FUV": "GALEX_FUV",
                "GALEX NUV": "GALEX_NUV",
                "SDSS u": "BiGauss_02.0",  # FWHM is actually variable
                "SDSS g": "BiGauss_02.0",  # FWHM is actually variable
                "SDSS r": "BiGauss_02.0",  # FWHM is actually variable
                "Mosaic Halpha": "Gauss_03.0",
                "SDSS i": "BiGauss_02.0",  # FWHM is actually variable
                "SDSS z": "BiGauss_02.0",  # FWHM is actually variable
                "2MASS J": "Gauss_03.5",
                "2MASS H": "Gauss_03.0",
                "2MASS Ks": "Gauss_03.5",
                "WISE W1": "WISE_FRAME_3.4",
                "IRAC I1": "IRAC_3.6",
                "IRAC I2": "IRAC_4.5",
                "WISE W2": "WISE_FRAME_4.6",
                "IRAC I3": "IRAC_5.8",
                "IRAC I4": "IRAC_8.0",
                "WISE W3": "WISE_FRAME_11.6",
                "WISE W4": "WISE_FRAME_22.1",
                "MIPS 24mu": "MIPS_24",
                "MIPS 70mu": None,
                "MIPS 160mu": None,
                "Pacs blue": "PACS_70",
                "Pacs red": "PACS_160",
                "SPIRE PSW": None,
                "SPIRE PMW": None,
                "SPIRE PLW": None}

# -----------------------------------------------------------------

class AnianoKernels(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # The path to the directory where the Aniano kernels are saved (to be reused)
        self.kernels_path = fs.join(kernels_path, "aniano")

    # -----------------------------------------------------------------

    def get_kernel(self, from_filter, to_filter, high_res=True, fwhm=None):

        """
        This function ...
        :param from_filter:
        :param to_filter:
        :param high_res:
        :param fwhm:
        :return:
        """

        # Get the local path to the kernel (will be downloaded if necessary)
        kernel_path, to_psf_name = self.get_kernel_path(from_filter, to_filter, high_res=high_res, fwhm=fwhm, return_name=True)

        # Load the kernel frame
        kernel = ConvolutionKernel.from_file(kernel_path)

        # Get the FWHM of the kernel (should have been done already!)
        if kernel.fwhm is None:

            # Find the appropriate FWHM
            if "Gauss" in to_psf_name or "Moffet" in to_psf_name: fwhm = float(to_psf_name.split("_")[1]) * Unit("arcsec")
            elif to_psf_name in fwhms: fwhm = fwhms[to_psf_name]
            else: fwhm = None

            # Set the FWHM of the kernel
            kernel.fwhm = fwhm

        # Return the kernel
        return kernel

    # -----------------------------------------------------------------

    def get_kernel_path(self, from_filter, to_filter, high_res=True, from_fwhm=None, to_fwhm=None, return_name=False, from_model="BiGauss", to_model="BiGauss"):

        """
        This function ...
        :param from_filter:
        :param to_filter:
        :param high_res:
        :param from_fwhm:
        :param to_fwhm:
        :param return_name:
        :param from_model: BiGauss is for SDSS
        :param to_model: BiGauss is for SDSS
        :return:
        """

        if isinstance(from_filter, basestring): from_filter = Filter.from_string(from_filter)
        if isinstance(to_filter, basestring): to_filter = Filter.from_string(to_filter)

        # For variable FWHM of input image
        if str(from_filter) in variable_fwhms: # Is SDSS

            # Check whether FWHM is specified
            if from_fwhm is None: raise ValueError("When convolving an SDSS image, the FWHM of that image must be specified")

            # Determine aniano name for the FWHM
            fwhm_arcsec = from_fwhm.to("arcsec").value
            aniano_name = from_model + "_0" + closest_half_integer_string(fwhm_arcsec) # BiGauss is for SDSS
            log.warning("The convolution kernel will be based on the FWHM of the original image, which is specified as " + str(from_fwhm) + ". Please check that this value is sensible. The aniano PSF that is taken for this FWHM is " + aniano_name)
            from_psf_name = aniano_name

        # Determine from_psf_name
        else: from_psf_name = aniano_names[str(from_filter)]

        # For variable FWHM of image with target resolution
        if str(to_filter) in variable_fwhms: # Is SDSS

            # Check whether FWHM is specified
            if to_fwhm is None: raise ValueError("When convolving to the resolution of a SDSS image, the FWHM of that image must be specified")

            # Determine aniano name for the FWHM
            fwhm_arcsec = to_fwhm.to("arcsec").value
            aniano_name = to_model + "_0" + closest_half_integer_string(fwhm_arcsec)  # BiGauss is for SDSS
            log.warning("The convolution kernel will be based on the FWHM of the target image, which is specified as " + str(to_fwhm) + ". Please check that this value is sensible. The aniano PSF that is taken for this FWHM is " + aniano_name)
            to_psf_name = aniano_name

        # Determine to_psf_name
        else: to_psf_name = aniano_names[str(to_filter)]

        # Determine the path to the kernel file
        if high_res: kernel_file_basename = "Kernel_HiRes_" + from_psf_name + "_to_" + to_psf_name
        else: kernel_file_basename = "Kernel_LoRes_" + from_psf_name + "_to_" + to_psf_name
        kernel_file_path = fs.join(self.kernels_path, kernel_file_basename + ".fits")

        # Download the kernel if it is not present
        if not fs.is_file(kernel_file_path):

            # Download the kernel
            self.download_kernel(kernel_file_basename)

            # Find the appropriate FWHM
            if "Gauss" in to_psf_name or "Moffet" in to_psf_name: fwhm = float(to_psf_name.split("_")[1]) * Unit("arcsec")
            elif to_psf_name in fwhms: fwhm = fwhms[to_psf_name]
            else: fwhm = None

            # Set the FWHM of the kernel
            kernel = ConvolutionKernel.from_file(kernel_file_path, fwhm=fwhm)
            #kernel.fwhm = fwhm
            kernel.saveto(kernel_file_path)

        # Return
        if return_name: return kernel_file_path, to_psf_name
        else: return kernel_file_path # Return the local kernel path

    # -----------------------------------------------------------------

    def get_psf(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Get the local path to the PSF file (will be downloaded if necessary)
        psf_path, psf_name = self.get_psf_path(fltr, return_name=True)

        # Load the PSF frame
        psf = ConvolutionKernel.from_file(psf_path)

        # Get the FWHM of the PSF
        if "Gauss" in psf_name or "Moffet" in psf_name: fwhm = float(psf_name.split("_")[1]) * Unit("arcsec")
        elif psf_name in fwhms: fwhm = fwhms[psf_name]
        else: fwhm = None

        # Set the FWHM of the PSF
        psf.fwhm = fwhm

        # Return the PSF
        return psf

    # -----------------------------------------------------------------

    def get_psf_path(self, fltr, return_name=False):

        """
        This function ...
        :param fltr:
        :param return_name:
        :return:
        """

        # Determine the aniano name for the PSF
        psf_name = aniano_names[str(fltr)]

        # Determine the path to the PSF file
        basename = "PSF_" + psf_name
        psf_file_path = fs.join(self.kernels_path, basename + ".fits")

        # Download the PSF file if it is not present
        if not fs.is_file(psf_file_path):

            # Download the PSF
            self.download_psf(basename)

            # Get the FWHM of the PSF
            if "Gauss" in psf_name or "Moffet" in psf_name: fwhm = float(psf_name.split("_")[1]) * Unit("arcsec")
            elif psf_name in fwhms: fwhm = fwhms[psf_name]
            else: fwhm = None

            # Set the FWHM of the PSF
            psf = ConvolutionKernel.from_file(psf_file_path, fwhm=fwhm)
            psf.saveto(psf_file_path)

        if return_name: return psf_file_path, psf_name
        else: return psf_file_path # Return the local PSF path

    # -----------------------------------------------------------------

    def get_fwhm(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Determine the aniano name for the PSF
        psf_name = aniano_names[str(fltr)]

        # Get the FWHM of the PSF
        if "Gauss" in psf_name or "Moffet" in psf_name: fwhm = float(psf_name.split("_")[1]) * Unit("arcsec")
        elif psf_name in fwhms: fwhm = fwhms[psf_name]
        else: fwhm = None

        # Return the FWHM
        return fwhm

    # -----------------------------------------------------------------

    @property
    def present_kernels(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.kernels_path, extension="fits", startswith="Kernel", returns="name")

    # -----------------------------------------------------------------

    @property
    def present_psfs(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.kernels_path, extension="fits", startswith="PSF", returns="name")

    # -----------------------------------------------------------------

    @property
    def online_kernels(self):

        """
        This function ...
        :return:
        """

        # From http://stackoverflow.com/questions/17709058/parsing-html-data-into-python-list-for-manipulation

        page_as_string = requests.get(aniano_kernels_highres_link).content

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

        # ...

        return my_results

    # -----------------------------------------------------------------

    @property
    def online_psfs(self):

        """
        This function ...
        :return:
        """

        return None

    # -----------------------------------------------------------------

    def download_kernel(self, kernel_basename):

        """
        This function ...
        :param kernel_basename:
        :return:
        """

        kernel_fitsname = kernel_basename + ".fits"
        kernel_gzname = kernel_basename + ".fits.gz"

        # Determine the link to the online kernel file
        link = aniano_kernels_highres_link if "HiRes" in kernel_basename else aniano_kernels_lowres_link
        kernel_link = link + kernel_gzname

        gz_path = fs.join(self.kernels_path, kernel_gzname)
        fits_path = fs.join(self.kernels_path, kernel_fitsname)

        # Inform the user
        log.info("Downloading kernel " + kernel_basename + " from " + link + " ...")

        # Download the kernel
        urllib.urlretrieve(kernel_link, gz_path)

        # Inform the user
        log.info("Decompressing kernel file ...")

        # Decompress the kernel FITS file
        archive.decompress_gz(gz_path, fits_path)

        # Remove the fits .gz file
        fs.remove_file(gz_path)

    # -----------------------------------------------------------------

    def download_psf(self, psf_basename):

        """
        This function ...
        :param psf_basename:
        :return:
        """

        psf_fitsname = psf_basename + ".fits"
        psf_gzname = psf_basename + ".fits.gz"

        # Determine the link to the online PSF file
        psf_link = aniano_psf_files_link + psf_gzname

        gz_path = fs.join(self.kernels_path, psf_gzname)
        fits_path = fs.join(self.kernels_path, psf_fitsname)

        # Inform the user
        log.info("Downloading PSF file " + psf_basename + " from " + aniano_psf_files_link + " ...")

        # Download the file
        urllib.urlretrieve(psf_link, gz_path)

        # Inform the user
        log.info("Decompressing PSF file ...")

        # Decompress the PSF FITS file
        archive.decompress_gz(gz_path, fits_path)

        # Remove the fits.gz file
        fs.remove_file(gz_path)

# -----------------------------------------------------------------

def closest_half_integer(number):

    """
    Round a number to the closest half integer.
    >>> closest_half_integer(1.3)
    1.5
    >>> closest_half_integer(2.6)
    2.5
    >>> closest_half_integer(3.0)
    3.0
    >>> closest_half_integer(4.1)
    4.0
    """

    return round(number * 2) / 2.0

# -----------------------------------------------------------------

def closest_half_integer_string(number):

    """
    This function ...
    :param number:
    :return:
    """

    value = closest_half_integer(number)
    return "{0:.1f}".format(round(value, 1))

# -----------------------------------------------------------------

# Reference: http://pypher.readthedocs.org/en/latest/

# -----------------------------------------------------------------

class PypherKernels(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # The path to the directory where the generated Pypher kernels are saved (to be reused)
        self.path = fs.join(kernels_path, "pypher")

# -----------------------------------------------------------------
