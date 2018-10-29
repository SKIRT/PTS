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
import requests
from lxml import html
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..core.kernel import ConvolutionKernel
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.tools import archive
from ...core.filter.filter import parse_filter
from ...core.units.parsing import parse_unit
from .kernels import Kernels, kernels_path, get_fwhm, has_variable_fwhm
from ...core.tools import types
from ...core.tools import network
from ..core import fits
from ...core.filter.broad import BroadBandFilter
from ...core.filter.narrow import NarrowBandFilter
from ...core.tools import strings
from ...core.tools.numbers import nearest_half_integer, nearest_integer

# -----------------------------------------------------------------

# The link to the Aniano high-resolution kernel files and the Aniano PSF files
aniano_kernels_highres_link_2012 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2012/Kernels_fits_Files/Hi_Resolution/"
aniano_kernels_lowres_link_2012 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2012/Kernels_fits_Files/Low_Resolution/"
aniano_psf_files_link_2012 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2012/PSF_fits_Files/"

# 2017 !: REMOVED!!!
#aniano_kernels_highres_link_2017 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2017/Kernels_fits_Files/Hi_Resolution/"
#aniano_kernels_lowres_link_2017 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2017/Kernels_fits_Files/Low_Resolution/"
#aniano_psf_files_link_2017 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2017/PSF_FITS_Files/"

# 2018: should contain everything
aniano_kernels_highres_link_2018 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2018/Kernels_FITS_Files/Hi_Resolution/"
aniano_kernels_lowres_link_2018 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2018/Kernels_FITS_Files/Low_Resolution/"
aniano_psf_files_link_2018 = "http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2018/PSF_FITS_Files/"

# -----------------------------------------------------------------

# aniano_kernels_highres_link = aniano_kernels_highres_link_2012
# aniano_kernels_lowres_link = aniano_kernels_lowres_link_2012
# aniano_psf_files_link = aniano_psf_files_link_2012

aniano_kernels_highres_link = aniano_kernels_highres_link_2018
aniano_kernels_lowres_link = aniano_kernels_lowres_link_2018
aniano_psf_files_link = aniano_psf_files_link_2018

# -----------------------------------------------------------------

aniano_names = {"UVOT UVM2": "Gauss_02.5",
                "UVOT UVW1": "Gauss_02.5",
                "UVOT UVW2": "Gauss_03.0",
                "GALEX FUV": "GALEX_FUV",
                "GALEX NUV": "GALEX_NUV",
                "SDSS u": "BiGauss_02.0",  # FWHM is actually variable
                "SDSS g": "BiGauss_02.0",  # FWHM is actually variable
                "SDSS r": "BiGauss_02.0",  # FWHM is actually variable
                "Halpha": "Gauss_03.0", # FWHM is actually variable
                "Ha": "Gauss_03.0", # FWHM is actually variable
                "SDSS i": "BiGauss_02.0",  # FWHM is actually variable
                "SDSS z": "BiGauss_02.0",  # FWHM is actually variable
                "2MASS J": "Gauss_03.5",   # FWHM is actually variable
                "2MASS H": "Gauss_03.0",   # FWHM is actually variable
                "2MASS Ks": "Gauss_03.5",  # FWHM is actually variable
                "UKIDSS J": "BiGauss_01.0",  # FWHM is actually variable
                "UKIDSS H": "BiGauss_01.0",  # FWHM is actually variable
                "UKIDSS K": "BiGauss_01.0",  # FWHM is actually variable
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
                "Pacs green": "PACS_100",
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
        return [parse_filter(key) for key in aniano_names.keys()]

    # -----------------------------------------------------------------

    @property
    def filter_names(self):
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

        #print(from_filter)
        #print(to_filter)
        return self.has_filter(from_filter) and self.has_filter(to_filter)

    # -----------------------------------------------------------------

    def has_psf_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr: 
        :return: 
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        # Check whether the filter string is in the Aniano names
        return str(fltr) in aniano_names

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
        kernel_path, to_psf_name = self.get_kernel_path(from_filter, to_filter, high_res=high_res, from_fwhm=from_fwhm, to_fwhm=to_fwhm, return_name=True, check_valid=False)

        # Load the kernel frame
        try: kernel = ConvolutionKernel.from_file(kernel_path, from_filter=from_filter, to_filter=to_filter)
        except fits.DamagedFITSFileError:
            log.warning("The kernel image was probably damaged. Removing it and downloading it again ...")
            fs.remove_file(kernel_path)
            kernel_path, to_psf_name = self.get_kernel_path(from_filter, to_filter, high_res=high_res, from_fwhm=from_fwhm, to_fwhm=to_fwhm, return_name=True, check_valid=False)
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

    def get_kernel_basename(self, from_filter, to_filter, high_res=True, from_fwhm=None, to_fwhm=None, from_model=None,
                            to_model=None, return_name=False):

        """
        This function ...
        :param from_filter:
        :param to_filter:
        :param high_res:
        :param from_fwhm:
        :param to_fwhm:
        :param from_model:
        :param to_model:
        :param return_name:
        :return:
        """

        # Convert filter strings into filters
        if types.is_string_type(from_filter): from_filter = parse_filter(from_filter)
        if types.is_string_type(to_filter): to_filter = parse_filter(to_filter)

        # For variable FWHM of input image
        if has_variable_fwhm(from_filter):

            # Check whether FWHM is specified
            if from_fwhm is None: raise ValueError("The FWHM of that image must be specified")

            # Set model
            # ARE THERE OTHER FILTERS FOR WHICH IT SHOULD BE BIGAUSS?
            #if isinstance(from_filter, BroadBandFilter) and from_filter.is_sdss: from_model = "BiGauss"
            #elif isinstance(from_filter, NarrowBandFilter) and from_filter.is_halpha: from_model = "BiGauss"
            #elif from_model is None: from_model = "Gauss" # default is Gauss

            # Determine aniano name for the FWHM
            fwhm_arcsec = from_fwhm.to("arcsec").value
            #aniano_name = from_model + "_0" + closest_half_integer_string(fwhm_arcsec)

            # Determine Gauss or BiGauss
            if from_model is None:
                if fwhm_arcsec <= 2.5: from_model = "BiGauss"
                else: fwhm_arcsec = "Gauss"

            # Under 10
            if int(round(fwhm_arcsec)) < 10: aniano_name = from_model + "_0" + closest_half_integer_string(fwhm_arcsec)

            # Above 10
            #else: aniano_name = from_model + "_" + closest_half_integer_string(fwhm_arcsec)
            else: aniano_name = from_model + "_" + closest_integer_string(fwhm_arcsec)

            log.warning("The convolution kernel will be based on the FWHM of the original image, which is specified as " + str(from_fwhm) + ". Please check that this value is sensible. The aniano PSF that is taken for this FWHM is " + aniano_name)
            from_psf_name = aniano_name

        # Determine from_psf_name
        else: from_psf_name = aniano_names[str(from_filter)]

        # For variable FWHM of image with target resolution
        if has_variable_fwhm(to_filter):

            # Check whether FWHM is specified
            if to_fwhm is None: raise ValueError("The FWHM of that image must be specified")

            # Set model
            # ARE THERE OTHER FILTERS FOR WHICH IT SHOULD BE BIGAUSS?
            #if isinstance(to_filter, BroadBandFilter) and to_filter.is_sdss: to_model = "BiGauss"
            #elif isinstance(to_filter, NarrowBandFilter) and to_filter.is_halpha: to_model = "BiGauss"
            #elif to_model is None: to_model = "Gauss" # default is Gauss

            # Determine aniano name for the FWHM
            fwhm_arcsec = to_fwhm.to("arcsec").value
            #aniano_name = to_model + "_0" + closest_half_integer_string(fwhm_arcsec)

            # Determine Gauss or BiGauss
            if to_model is None:
                if fwhm_arcsec <= 2.5: to_model = "BiGauss"
                else: to_model = "Gauss"

            # Under 10
            if int(round(fwhm_arcsec)) < 10: aniano_name = to_model + "_0" + closest_half_integer_string(fwhm_arcsec)

            # Above 10
            #else: aniano_name = to_model + "_" + closest_half_integer_string(fwhm_arcsec)
            else: aniano_name = to_model + "_" + closest_integer_string(fwhm_arcsec)

            log.warning("The convolution kernel will be based on the FWHM of the target image, which is specified as " + str(to_fwhm) + ". Please check that this value is sensible. The aniano PSF that is taken for this FWHM is " + aniano_name)
            to_psf_name = aniano_name

        # Determine to_psf_name
        else: to_psf_name = aniano_names[str(to_filter)]

        # Determine the path to the kernel file
        if high_res: kernel_file_basename = "Kernel_HiRes_" + from_psf_name + "_to_" + to_psf_name
        else: kernel_file_basename = "Kernel_LoRes_" + from_psf_name + "_to_" + to_psf_name

        # Return
        if return_name: return kernel_file_basename, to_psf_name
        else: return kernel_file_basename

    # -----------------------------------------------------------------

    def get_kernel_path(self, from_filter, to_filter, high_res=True, from_fwhm=None, to_fwhm=None,
                        return_name=False, from_model=None, to_model=None, check_valid=True):

        """
        This function ...
        :param from_filter:
        :param to_filter:
        :param high_res:
        :param from_fwhm:
        :param to_fwhm:
        :param return_name:
        :param from_model:
        :param to_model:
        :param check_valid:
        :return:
        """

        # Get the kernel file basename
        kernel_file_basename, to_psf_name = self.get_kernel_basename(from_filter, to_filter, high_res=high_res, from_fwhm=from_fwhm,
                                                        to_fwhm=to_fwhm, from_model=from_model, to_model=to_model, return_name=True)

        # Determine kernel path
        kernel_file_path = fs.join(self.kernels_path, kernel_file_basename + ".fits")

        # Download the kernel if it is not present
        if not fs.is_file(kernel_file_path):

            # Download the kernel
            self.download_kernel(kernel_file_basename)

            # Initialize the kernel file
            self.initialize_kernel_file(to_psf_name, kernel_file_path, from_filter, to_filter)

        # CHeck whether the file is OK
        if check_valid and not fits.is_valid(kernel_file_path):

            # Give warning
            log.warning("The kernel file is damaged. Removing it and downloading it again ...")

            # Remove damaged file
            fs.remove_file(kernel_file_path)

            # Download the PSF
            self.download_kernel(kernel_file_basename)

            # Initialize the kernel file
            self.initialize_kernel_file(to_psf_name, kernel_file_path, from_filter, to_filter)

        # Return
        if return_name: return kernel_file_path, to_psf_name
        else: return kernel_file_path # Return the local kernel path

    # -----------------------------------------------------------------

    def initialize_kernel_file(self, to_psf_name, filepath, from_filter, to_filter):

        """
        This function ...
        :param kernel_basename:
        :param to_psf_name:
        :param filepath:
        :param from_filter:
        :param to_filter:
        :return:
        """

        # Debugging
        log.debug("Initializing the kernel file from the '" + str(from_filter) + "' to the '" + str(to_filter) + "' filter to '" + filepath + "' ...")

        # Find the appropriate FWHM
        fwhm = get_fwhm_for_aniano_name(to_psf_name)
        if fwhm is None: log.warning("The FWHM of the kernel is undefined")

        # Set the FWHM of the kernel
        kernel = ConvolutionKernel.from_file(filepath, fwhm=fwhm, from_filter=from_filter, to_filter=to_filter)
        kernel.saveto(filepath)

    # -----------------------------------------------------------------

    def get_psf(self, fltr, fwhm=None):

        """
        This function ...
        :param fltr:
        :param fwhm:
        :return:
        """

        # Parse filter
        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        # Get the local path to the PSF file (will be downloaded if necessary)
        psf_path, psf_name = self.get_psf_path(fltr, return_name=True, check_valid=False, fwhm=fwhm)

        # Load the PSF frame
        try: psf = ConvolutionKernel.from_file(psf_path, to_filter=fltr)
        except fits.DamagedFITSFileError:
            log.warning("The kernel image was probably damaged. Removing it and downloading it again ...")
            fs.remove_file(psf_path)
            psf_path, psf_name = self.get_psf_path(fltr, return_name=True, check_valid=False, fwhm=fwhm)
            psf = ConvolutionKernel.from_file(psf_path, to_filter=fltr)

        # Get the FWHM of the PSF
        fwhm_aniano = get_fwhm_for_aniano_name(psf_name)
        if fwhm is None: fwhm = fwhm_aniano
        else:
            rdiff = (fwhm.to("arcsec").value - fwhm_aniano.to("arcsec").value) / fwhm_aniano.to("arcsec").value
            if rdiff > 0.05: raise ValueError("The specified FWHM is more than 5% off compared to the value corresponding to the PSF '" + psf_name + "'")

        # Set the FWHM of the PSF
        psf.fwhm = fwhm

        # Return the PSF
        return psf

    # -----------------------------------------------------------------

    def get_psf_basename(self, fltr, return_name=False):

        """
        This function ...
        :param fltr:
        :param return_name:
        :return:
        """

        # Convert to filter
        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        # Determine the aniano name for the PSF
        psf_name = aniano_names[str(fltr)]

        # Determine the path to the PSF file
        if aniano_psf_files_link == aniano_psf_files_link_2018: basename = "PSF_Original_" + psf_name
        elif aniano_psf_files_link == aniano_psf_files_link_2012: basename = "PSF_" + psf_name
        else: raise ValueError("Invalid aniano PSF files URL '" + aniano_psf_files_link + "'")

        # Return
        if return_name: return basename, psf_name
        else: return basename

    # -----------------------------------------------------------------

    def get_psf_path(self, fltr, return_name=False, check_valid=True, fwhm=None):

        """
        This function ...
        :param fltr:
        :param return_name:
        :param check_valid:
        :param fwhm:
        :return:
        """

        # Get the PSF basename
        basename, psf_name = self.get_psf_basename(fltr, return_name=True)

        # Set path
        psf_file_path = fs.join(self.kernels_path, basename + ".fits")
        #print(psf_file_path)

        # Determine the (potential) compressed filepath
        compressed_psf_file_path = fs.join(self.kernels_path, basename + ".fits.gz")

        # Download the PSF file if it is not present
        if not fs.is_file(psf_file_path):

            # The compressed file is present
            if fs.is_file(compressed_psf_file_path):

                #print(compressed_psf_file_path)

                # Check whether the file is indeed compressed
                if archive.is_compressed(compressed_psf_file_path):

                    # Initialize the PSF file
                    self.initialize_compressed_psf_file(psf_name, compressed_psf_file_path, fltr, fwhm=fwhm)

                # Not actually compressed, just rename
                else:

                    # Rename
                    fs.rename_file_path(compressed_psf_file_path, basename + ".fits")

                    # Initialize the file
                    self.initialize_psf_file(psf_name, psf_file_path, fltr, fwhm=fwhm)

            # No file is present
            else:

                # Download the PSF
                self.download_psf(basename)

                # Initialize the PSF file
                self.initialize_psf_file(psf_name, psf_file_path, fltr, fwhm=fwhm)

        # Check whether maybe the regular FITS file is actually a compressed file
        elif archive.is_compressed(psf_file_path):

            # Initialize the PSF file
            self.initialize_compressed_psf_file(psf_name, psf_file_path, fltr, fwhm=fwhm)

        # CHeck whether the file is OK
        if check_valid and not fits.is_valid(psf_file_path):

            # Give warning
            log.warning("The PSF file is damaged. Removing it and downloading it again ...")

            # Remove damaged file
            fs.remove_file(psf_file_path)

            # Download the PSF
            self.download_psf(basename)

            # Initialize the PSF file
            self.initialize_psf_file(psf_name, psf_file_path, fltr, fwhm=fwhm)

        # Return
        if return_name: return psf_file_path, psf_name
        else: return psf_file_path # Return the local PSF path

    # -----------------------------------------------------------------

    def initialize_compressed_psf_file(self, psf_name, filepath, fltr, fwhm=None):

        """
        This function ...
        :param psf_name:
        :param filepath:
        :param fltr:
        :param fwhm:
        :return:
        """

        # Debugging
        log.debug("Decompressing the '" + psf_name + "' PSF file ...")

        # Fix extension
        filepath = archive.fix_extension(filepath)

        # Unpack
        fits_path = archive.decompress_file_in_place(filepath, remove=True)

        # Initialize psf file
        self.initialize_psf_file(psf_name, fits_path, fltr, fwhm=fwhm)

    # -----------------------------------------------------------------

    def initialize_psf_file(self, psf_name, filepath, fltr, fwhm=None):

        """
        This function ...
        :param psf_name:
        :param filepath:
        :param fltr:
        :param fwhm:
        :return:
        """

        # Debugging
        log.debug("Initializing the '" + psf_name + "' PSF file to '" + filepath + "' ...")

        # Get the FWHM of the PSF
        fwhm_aniano = get_fwhm_for_aniano_name(psf_name)
        if fwhm is None: fwhm = fwhm_aniano
        else:
            rdiff = (fwhm.to("arcsec").value - fwhm_aniano.to("arcsec").value) / fwhm_aniano.to("arcsec").value
            if rdiff > 0.05: raise ValueError("The specified FWHM is more than 5% off compared to the value corresponding to the PSF '" + psf_name + "'")

        # Set the FWHM of the PSF
        psf = ConvolutionKernel.from_file(filepath, fwhm=fwhm, to_filter=fltr)
        psf.saveto(filepath)

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

        # NEW: Check whether the link exists!!
        #if not network.exists(kernel_link): raise ValueError("The '" + kernel_basename + "' does not exist!")

        # Determine filepath
        gz_path = fs.join(self.kernels_path, kernel_gzname)
        fits_path = fs.join(self.kernels_path, kernel_fitsname)

        # Inform the user
        log.info("Downloading kernel " + kernel_basename + " from " + link + " ...")

        # CHECK
        if not network.is_available(kernel_link):

            # Was 2018
            if aniano_kernels_highres_link == aniano_kernels_highres_link_2018:

                log.warning("The 2018 " + kernel_basename + " kernel file is not present, trying the 2012 version ...")
                new_aniano_kernels_highres_link = aniano_kernels_highres_link_2012

            # Was 2012
            elif aniano_kernels_highres_link == aniano_kernels_highres_link_2012:

                log.warning("The 2012 " + kernel_basename + " kernel file is not present, trying the 2018 version ...")
                new_aniano_kernels_highres_link = aniano_kernels_highres_link_2018

            # Invalid
            else: raise ValueError("Unknown Aniano High-Res kernels URL '" + aniano_kernels_highres_link + "'")

            # TRY AGAIN
            kernel_link = new_aniano_kernels_highres_link + kernel_gzname

            # CHECK AGAIN
            if not network.is_available(kernel_link):

                # Try finding Gauss instead of BiGauss kernel
                if "to_BiGauss" in kernel_link:

                    # Replace with Gauss
                    new_kernel_link = strings.replace_last(kernel_link, "BiGauss", "Gauss")

                    # Check
                    if not network.is_available(new_kernel_link): raise ValueError("Cannot find the " + kernel_basename + " kernel file")

                    # Use
                    else:
                        log.warning("The '" + kernel_basename + "' kernel does not exist: taking the '" + strings.replace_last(kernel_basename, "BiGauss", "Gauss") + "' kernel instead ...")
                        kernel_link = new_kernel_link

                # Not found
                else: raise ValueError("Cannot find the " + kernel_basename + " kernel file")

        # Download the kernel
        network.download_file(kernel_link, gz_path, progress_bar=log.is_debug)

        # Inform the user
        log.info("Decompressing kernel file ...")

        # Decompress the PSF FITS file
        if archive.is_compressed(gz_path):

            # Decompress the kernel FITS file
            archive.decompress_gz(gz_path, fits_path)

            # Remove the fits .gz file
            fs.remove_file(gz_path)

        # Remove the gz extension
        else: fs.remove_extension(gz_path)

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

        # CHECK
        #print(aniano_psf_files_link, network.is_available(aniano_psf_files_link))
        #print(psf_link, network.is_available(psf_link))
        if not network.is_available(psf_link):

            # Was 2018
            if aniano_psf_files_link == aniano_psf_files_link_2018:

                log.warning("The 2018 " + psf_basename + " PSF file is not present, trying the 2012 version ...")
                new_aniano_psf_files_link = aniano_psf_files_link_2012
                new_psf_gzname = psf_gzname.replace("Original_", "")

            # Was 2012
            elif aniano_psf_files_link == aniano_psf_files_link_2012:

                log.warning("The 2012 " + psf_basename + " PSF file is not present, trying the 2018 version ...")
                new_aniano_psf_files_link = aniano_psf_files_link_2018
                new_psf_gzname = psf_gzname.replace("PSF_", "PSF_Original_")

            # Invalid
            else: raise ValueError("Unknown Aniano PSF files URL '" + aniano_psf_files_link + "'")

            # TRY AGAIN
            psf_link = new_aniano_psf_files_link + new_psf_gzname
            #print(psf_link)

            # CHECK AGAIN
            if not network.is_available(psf_link):

                # Try to find Gauss instead of BiGauss
                if "BiGauss" in psf_link:

                    # Replace with Gauss
                    new_psf_link = strings.replace(psf_link, "BiGauss", "Gauss")

                    # Check
                    if not network.is_available(new_psf_link): raise ValueError("Cannot find the " + psf_basename + " PSF file")

                    # Use
                    else:
                        log.warning("The '" + psf_basename + "' PSF file does not exist: taking the '" + strings.replace(psf_basename, "BiGauss", "Gauss") + "' PSF instead ...")
                        psf_link = new_psf_link

                # Not found
                else: raise ValueError("Cannot find the " + psf_basename + " PSF file")

        # Download the file
        network.download_file(psf_link, gz_path, progress_bar=log.is_debug)

        # Inform the user
        log.info("Decompressing PSF file if necessary ...")

        # Decompress the PSF FITS file
        if archive.is_compressed(gz_path):

            # Make decompressed file
            archive.decompress_gz(gz_path, fits_path)

            # Remove the fits.gz file
            fs.remove_file(gz_path)

        # Remove the gz extension
        else: fs.remove_extension(gz_path)

# -----------------------------------------------------------------

def closest_half_integer_string(number):

    """
    This function ...
    :param number:
    :return:
    """

    value = nearest_half_integer(number)
    return "{0:.1f}".format(round(value, 1))

# -----------------------------------------------------------------

def closest_integer_string(number):

    """
    Thisn function ...
    :param number:
    :return:
    """

    value = nearest_integer(number)
    return repr(value)

# -----------------------------------------------------------------
