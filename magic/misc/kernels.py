#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.kernels Contains the AnianoKernels class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import gzip
import shutil
import urllib
import requests
from lxml import html
from collections import OrderedDict

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ..core.frame import Frame
from ...core.tools import inspection, filesystem
from ...core.tools.logging import log

# -----------------------------------------------------------------

# The path to the PTS kernels directory
kernels_path = filesystem.join(inspection.pts_user_dir, "kernels")

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
        self.kernels_path = filesystem.join(kernels_path, "aniano")

    # -----------------------------------------------------------------

    def get_kernel(self, from_instrument, to_instrument, high_res=True):

        """
        This function ...
        :param from_instrument:
        :param to_instrument:
        :param high_res:
        :return:
        """

        # Get the local path to the kernel (will be downloaded if necessary)
        kernel_path = self.get_kernel_path(from_instrument, to_instrument, high_res=high_res)

        # Load the kernel frame
        kernel = Frame.from_file(kernel_path)

        # Get the FWHM of the kernel (should have been done already!)
        if kernel.fwhm is None:

            # Find the appropriate FWHM
            if "Gauss" in to_instrument or "Moffet" in to_instrument: fwhm = float(to_instrument.split("_")[1]) * Unit("arcsec")
            elif to_instrument in fwhms: fwhm = fwhms[to_instrument]
            else: fwhm = None

            # Set the FWHM of the kernel
            kernel.fwhm = fwhm

        # Return the kernel
        return kernel

    # -----------------------------------------------------------------

    def get_kernel_path(self, from_instrument, to_instrument, high_res=True):

        """
        This function ...
        :param from_instrument:
        :param to_instrument:
        :param high_res:
        :return:
        """

        # TODO: allow strings that are not literally the strings appearing in the Aniano FITS file names
        from_psf_name = from_instrument
        to_psf_name = to_instrument

        # Determine the path to the kernel file
        if high_res: kernel_file_basename = "Kernel_HiRes_" + from_psf_name + "_to_" + to_psf_name
        else: kernel_file_basename = "Kernel_LoRes_" + from_psf_name + "_to_" + to_psf_name
        kernel_file_path = filesystem.join(self.kernels_path, kernel_file_basename + ".fits")

        # Download the kernel if it is not present
        if not filesystem.is_file(kernel_file_path):

            # Download the kernel
            self.download_kernel(kernel_file_basename)

            # Find the appropriate FWHM
            if "Gauss" in to_instrument or "Moffet" in to_instrument: fwhm = float(to_instrument.split("_")[1]) * Unit("arcsec")
            elif to_instrument in fwhms: fwhm = fwhms[to_instrument]
            else: fwhm = None

            # Set the FWHM of the kernel
            kernel = Frame.from_file(kernel_file_path)
            kernel.fwhm = fwhm
            kernel.save(kernel_file_path)

        # Return the local kernel path
        return kernel_file_path

    # -----------------------------------------------------------------

    def get_psf(self, instrument):

        """
        This function ...
        :param instrument:
        :return:
        """

        # Get the local path to the PSF file (will be downloaded if necessary)
        psf_path = self.get_psf_path(instrument)

        # Load the PSF frame
        psf = Frame.from_file(psf_path)

        # Get the FWHM of the PSF
        if "Gauss" in instrument or "Moffet" in instrument: fwhm = float(instrument.split("_")[1]) * Unit("arcsec")
        elif instrument in fwhms: fwhm = fwhms[instrument]
        else: fwhm = None

        # Set the FWHM of the PSF
        psf.fwhm = fwhm

        # Return the PSF
        return psf

    # -----------------------------------------------------------------

    def get_psf_path(self, instrument):

        """
        This function ...
        :param instrument:
        :return:
        """

        # TODO: allow strings that are not literally the strings appearing in the Aniano FITS file names
        psf_name = instrument

        # Determine the path to the PSF file
        basename = "PSF_" + psf_name
        psf_file_path = filesystem.join(self.kernels_path, basename + ".fits")

        # Download the PSF file if it is not present
        if not filesystem.is_file(psf_file_path): self.download_psf(basename)

        # Return the local PSF path
        return psf_file_path

    # -----------------------------------------------------------------

    @property
    def present_kernels(self):

        """
        This function ...
        :return:
        """

        return filesystem.files_in_path(self.kernels_path, extension="fits", startswith="Kernel", returns="name")

    # -----------------------------------------------------------------

    @property
    def present_psfs(self):

        """
        This function ...
        :return:
        """

        return filesystem.files_in_path(self.kernels_path, extension="fits", startswith="PSF", returns="name")

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

        gz_path = filesystem.join(self.kernels_path, kernel_gzname)
        fits_path = filesystem.join(self.kernels_path, kernel_fitsname)

        # Inform the user
        log.info("Downloading kernel " + kernel_basename + " from " + link + " ...")

        # Download the kernel
        urllib.urlretrieve(kernel_link, gz_path)

        # Inform the user
        log.info("Unzipping kernel ...")

        # Unzip the kernel FITS file
        with gzip.open(gz_path, 'rb') as f_in:
            with open(fits_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Remove the fits.gz file
        filesystem.remove_file(gz_path)

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

        gz_path = filesystem.join(self.kernels_path, psf_gzname)
        fits_path = filesystem.join(self.kernels_path, psf_fitsname)

        # Inform the user
        log.info("Downloading PSF file " + psf_basename + " from " + aniano_psf_files_link + " ...")

        # Download the file
        urllib.urlretrieve(psf_link, gz_path)

        # Inform the user
        log.info("Unzipping PSF file ...")

        # Unzip the FITS file
        with gzip.open(gz_path, 'rb') as f_in:
            with open(fits_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Remove the fits.gz file
        filesystem.remove_file(gz_path)

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
        self.path = filesystem.join(kernels_path, "pypher")

# -----------------------------------------------------------------
