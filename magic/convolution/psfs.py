#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.psfs Contains classes that can be used to access PSFs.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ..core.kernel import ConvolutionKernel
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.tools import network

# -----------------------------------------------------------------

if not fs.is_directory(introspection.pts_ext_dir): fs.create_directory(introspection.pts_ext_dir)
psfs_path = fs.join(introspection.pts_ext_dir, "psfs")
if not fs.is_directory(psfs_path): fs.create_directory(psfs_path)

# -----------------------------------------------------------------

pacs_fwhms = [5.67 * Unit("arcsec"), 7.04 * Unit("arcsec"), 11.18 * Unit("arcsec")]
spire_fwhms = [18.15 * Unit("arcsec"), 24.88 * Unit("arcsec"), 36.09 * Unit("arcsec")]

# -----------------------------------------------------------------

class HerschelPSFs(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        herschel_psfs_path = fs.join(psfs_path, "Herschel")
        if not fs.is_directory(herschel_psfs_path): fs.create_directory(herschel_psfs_path)

        self.path = herschel_psfs_path

        # The path to the directory where the Aniano kernels are saved (to be reused)
        self.spire_path = fs.create_directory_in(herschel_psfs_path, "spire")
        self.pacs_path = fs.create_directory_in(herschel_psfs_path, "pacs")

    # -----------------------------------------------------------------

    def get_spire_psf(self, band, one_arcsec=False):

        """
        This function ...
        :param band:
        :param one_arcsec:
        :return:
        """

        if band == 1 or band == "PSW":

            if one_arcsec: psf_path = fs.join(self.spire_path, "0x5000241aL_PSW_bgmod10_1arcsec.fits")
            else: psf_path = fs.join(self.spire_path, "0x5000241aL_PSW_bgmod10_6arcsec.fits")
            fwhm = spire_fwhms[0]

        elif band == 2 or band == "PMW":

            if one_arcsec: psf_path = fs.join(self.spire_path, "0x5000241aL_PMW_bgmod10_1arcsec.fits")
            else: psf_path = fs.join(self.spire_path, "0x5000241aL_PMW_bgmod10_10arcsec.fits")
            fwhm = spire_fwhms[1]

        elif band == 3 or band == "PLW":

            if one_arcsec: psf_path = fs.join(self.spire_path, "0x5000241aL_PLW_bgmod10_1arcsec.fits")
            else: psf_path = fs.join(self.spire_path, "0x5000241aL_PLW_bgmod10_14arcsec.fits")
            fwhm = spire_fwhms[2]

        else: raise ValueError("Invalid option for 'band'")

        # Load the PSF frame
        psf = ConvolutionKernel.from_file(psf_path, fwhm=fwhm)

        # Set the FWHM of the PSF
        #if psf.fwhm is None: psf.fwhm = fwhm

        # Return the PSF frame
        return psf

    # -----------------------------------------------------------------

    def get_pacs_psf(self, band):

        """
        This function ...
        :param band:
        :return:
        """

        if band == 1 or band == "70" or band == "70mu":

            psf_path = fs.join()
            fwhm = pacs_fwhms[0]

        elif band == 2 or band == "100" or band == "100mu":

            psf_path = fs.join()
            fwhm = pacs_fwhms[1]

        elif band == 3 or band == "160" or band == "160mu":

            psf_path = fs.join()
            fwhm = pacs_fwhms[2]

        else: raise ValueError("Invalid option for 'band'")

        # Load the PSF frame
        psf = ConvolutionKernel.from_file(psf_path, fwhm=fwhm)

        # Set the FWHM of the PSF
        #if psf.fwhm is None: psf.fwhm = fwhm

        # Return the PSF frame
        return psf

# -----------------------------------------------------------------

# URLs
core_psfs_url = "http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/psfprf/apex_core_irac_PRFs.tar"
extended_psfs_url = "http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/psfprf/Cryo.IRAC.Extended.PSF.5X.070704.tar.gz"
extended_psfs_warm_url = "http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/psfprf/Warm.IRAC.Extended.PSF.5X.091113.tar.gz"

# -----------------------------------------------------------------

class IRACPSFs(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        irac_psfs_path = fs.join(psfs_path, "IRAC")
        if not fs.is_directory(irac_psfs_path): fs.create_directory(irac_psfs_path)

        self.path = irac_psfs_path

        # The path to the directory where the Aniano kernels are saved (to be reused)
        self.core_path = fs.create_directory_in(irac_psfs_path, "core")
        self.extended_path = fs.create_directory_in(irac_psfs_path, "extended")

    # -----------------------------------------------------------------

    def download_core_psfs(self):

        """
        This function ...
        :return:
        """

        # Download and decompress file
        network.download_and_decompress_file(core_psfs_url, self.path)

    # -----------------------------------------------------------------

    def download_extended_psfs(self):

        """
        This function ...
        :return:
        """

        # Download and decompress file
        network.download_and_decompress_file(extended_psfs_url, self.path)

    # -----------------------------------------------------------------

    def download_warm_extended_psfs(self):

        """
        This function ...
        :return:
        """

        # Download and decompress file
        network.download_and_decompress_file(extended_psfs_warm_url, self.path)

    # -----------------------------------------------------------------

    def get_psf_path(self, band, extended=True, warm=False):

        """
        This function ...
        :param band: 1, 2, 3, or 4
        :param extended:
        :param warm:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get_psf(self, band, extended=True, warm=False):

        """
        This function ...
        :param band: 1, 2, 3, or 4
        :param extended:
        :param warm:
        :return:
        """

        pass

# -----------------------------------------------------------------
