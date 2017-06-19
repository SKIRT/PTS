#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.convolution.aniano Contains the AnianoKernels class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import urllib
import requests
from lxml import html
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..core.kernel import ConvolutionKernel
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.tools import archive
from ...core.filter.filter import parse_filter
from ...core.units.parsing import parse_unit
from .kernels import Kernels, kernels_path, get_fwhm, has_variable_fwhm
from ...core.tools import types

# -----------------------------------------------------------------

# The link to the Aniano high-resolution kernel files and the Aniano PSF files
aniano_kernels_highres_link_2012 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2012/Kernels_fits_Files/Hi_Resolution/"
aniano_kernels_lowres_link_2012 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2012/Kernels_fits_Files/Low_Resolution/"
aniano_psf_files_link_2012 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2012/PSF_fits_Files/"

# 2017 !
#aniano_kernels_highres_link = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2017/Kernels_fits_Files/Hi_Resolution/"
#aniano_kernels_lowres_link = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2017/Kernels_fits_Files/Low_Resolution/"
#aniano_psf_files_link = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2017/PSF_fits_Files/"

aniano_kernels_highres_link = aniano_kernels_highres_link_2012
aniano_kernels_lowres_link = aniano_kernels_lowres_link_2012
aniano_psf_files_link = aniano_psf_files_link_2012

# -----------------------------------------------------------------

aniano_names = {"UVOT UVM2": "Gauss_02.5",
                "UVOT UVW1": "Gauss_02.5",
                "UVOT UVW2": "Gauss_03.0",
                "GALEX FUV": "GALEX_FUV",
                "GALEX NUV": "GALEX_NUV",
                "SDSS u": "BiGauss_02.0",  # FWHM is actually variable
                "SDSS g": "BiGauss_02.0",  # FWHM is actually variable
                "SDSS r": "BiGauss_02.0",  # FWHM is actually variable
                "Halpha": "Gauss_03.0",
                "SDSS i": "BiGauss_02.0",  # FWHM is actually variable
                "SDSS z": "BiGauss_02.0",  # FWHM is actually variable
                "2MASS J": "Gauss_03.5",   # FWHM is actually variable
                "2MASS H": "Gauss_03.0",   # FWHM is actually variable
                "2MASS Ks": "Gauss_03.5",  # FWHM is actually variable
                "WISE W1": "WISE_FRAME_3.4",
                "IRAC I1": "IRAC_3.6",
                "IRAC I2": "IRAC_4.5",
                "WISE W2": "WISE_FRAME_4.6",
                "IRAC I3": "IRAC_5.8",
                "IRAC I4": "IRAC_8.0",
                "WISE W3": "WISE_FRAME_11.6",
                "WISE W4": "WISE_FRAME_22.1",
                "MIPS 24mu": "MIPS_24",
                "MIPS 70mu": "MIPS_70",
                "MIPS 160mu": "MIPS_160",
                "Pacs blue": "PACS_70",
                "Pacs red": "PACS_160",
                "SPIRE PSW": "SPIRE_250",
                "SPIRE PMW": "SPIRE_350",
                "SPIRE PLW": "SPIRE_500"}

# -----------------------------------------------------------------

filter_names = {"GALEX_FUV": "GALEX FUV",
                 "GALEX_NUV": "GALEX NUV",
                 "IRAC_3.6": "IRAC I1",
                 "IRAC_4.5": "IRAC I2",
                 "IRAC_5.8": "IRAC I3",
                 "IRAC_8.0": "IRAC I4",
                 "WISE_FRAME_3.4": "WISE W1",
                 "WISE_FRAME_4.6": "WISE W2",
                 "WISE_FRAME_11.6": "WISE W3",
                 "WISE_FRAME_22.1": "WISE W4",
                 "MIPS_24": "Mips 24mu",
                 "MIPS_70": "Mips 70mu",
                 "MIPS_160": "Mips 160mu",
                 "PACS_70": "Pacs blue",
                 "PACS_100": "Pacs green",
                 "PACS_160": "Pacs red",
                 "SPIRE_250": "SPIRE PSW",
                 "SPIRE_350": "SPIRE PMW",
                 "SPIRE_500": "SPIRE PLW"}

# -----------------------------------------------------------------

def get_filters_and_aniano_names_lists():

    """
    This function ...
    :return:
    """

    filters = []
    names = []

    for filter_name, aniano_name in aniano_names.iteritems():

        filters.append(parse_filter(filter_name))
        names.append(aniano_name)

    # Return the lists
    return filters, names

# -----------------------------------------------------------------

def get_aniano_name_from_filter(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    filters, names = get_filters_and_aniano_names_lists()
    index = filters.index(fltr)
    return names[index]

# -----------------------------------------------------------------

def get_filter_from_aniano_name(name):

    """
    This function ...
    :param name:
    :return:
    """

    return parse_filter(filter_names[name])

# -----------------------------------------------------------------

def get_fwhm_for_aniano_name(psf_name):

    """
    This function ...
    :param psf_name:
    :return:
    """

    # Find the appropriate FWHM
    if "Gauss" in psf_name or "Moffet" in psf_name: fwhm = float(psf_name.split("_")[1]) * parse_unit("arcsec")
    elif psf_name in filter_names:
        fltr = get_filter_from_aniano_name(psf_name)
        fwhm = get_fwhm(fltr)
    else: fwhm = None

    return fwhm

# -----------------------------------------------------------------

class AnianoKernels(Kernels):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Call the constructor of the base class
        super(AnianoKernels, self).__init__()

        # The path to the directory where the Aniano kernels are saved (to be reused)
        self.kernels_path = fs.create_directory_in(kernels_path, "aniano")

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return: 
        """

        return [parse_filter(key) for key in aniano_names.keys()]

    # -----------------------------------------------------------------

    @property
    def filter_names(self):

        """
        This function ...
        :return: 
        """

        return [str(fltr) for fltr in self.filters]

    # -----------------------------------------------------------------

    def has_filter(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        return fltr in self.filters

    # -----------------------------------------------------------------

    def has_kernel_for_filters(self, from_filter, to_filter):

        """
        This function ...
        :param from_filter: 
        :param to_filter: 
        :return: 
        """

        return self.has_filter(from_filter) and self.has_filter(to_filter)

    # -----------------------------------------------------------------

    def get_kernel(self, from_filter, to_filter, high_res=True, from_fwhm=None, to_fwhm=None):

        """
        This function ...
        :param from_filter:
        :param to_filter:
        :param high_res:
        :param from_fwhm:
        :param to_fwhm:
        :return:
        """

        # Get the local path to the kernel (will be downloaded if necessary)
        kernel_path, to_psf_name = self.get_kernel_path(from_filter, to_filter, high_res=high_res, from_fwhm=from_fwhm, to_fwhm=to_fwhm, return_name=True)

        # Load the kernel frame
        kernel = ConvolutionKernel.from_file(kernel_path, from_filter=from_filter, to_filter=to_filter)

        # Get the FWHM of the kernel (should have been done already!)
        if kernel.fwhm is None:

            # Find the appropriate FWHM
            fwhm = get_fwhm_for_aniano_name(to_psf_name)
            if fwhm is None: fwhm = to_fwhm

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

        # Convert filter strings into filters
        if types.is_string_type(from_filter): from_filter = parse_filter(from_filter)
        if types.is_string_type(to_filter): to_filter = parse_filter(to_filter)

        # For variable FWHM of input image
        if has_variable_fwhm(from_filter): # Is SDSS

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
        if has_variable_fwhm(to_filter): # Is SDSS or 2MASS

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
            fwhm = get_fwhm_for_aniano_name(to_psf_name)

            # Set the FWHM of the kernel
            kernel = ConvolutionKernel.from_file(kernel_file_path, fwhm=fwhm, from_filter=from_filter, to_filter=to_filter)
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

        # Parse filter
        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        # Get the local path to the PSF file (will be downloaded if necessary)
        psf_path, psf_name = self.get_psf_path(fltr, return_name=True)

        # Load the PSF frame
        psf = ConvolutionKernel.from_file(psf_path, to_filter=fltr)

        # Get the FWHM of the PSF
        fwhm = get_fwhm_for_aniano_name(psf_name)

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
            fwhm = get_fwhm_for_aniano_name(psf_name)

            # Set the FWHM of the PSF
            psf = ConvolutionKernel.from_file(psf_file_path, fwhm=fwhm, to_filter=fltr)
            psf.saveto(psf_file_path)

        if return_name: return psf_file_path, psf_name
        else: return psf_file_path # Return the local PSF path

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
