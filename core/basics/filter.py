#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.filter Working with wavelength filters (or bands).
#
# An instance of the Filter class in this module represents a particular wavelength band, including its response or
# transmission curve, and allows integrating a given spectrum over the band. A filter instance can be constructed from
# one of the provided resource files describing specific instruments or standard wavelength bands. Alternatively,
# a filter can be created with a uniform transmission curve over a certain wavelength range.

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import sys
import types
import numpy as np
from scipy.interpolate import interp1d
from lxml import etree

# -----------------------------------------------------------------
#  Filter class
# -----------------------------------------------------------------

## An instance of the Filter class represents a particular wavelength band, including its response or transmission
# curve and some basic properties such as its mean and pivot wavelengths. The class provides a function to
# integrate a given spectrum over the band. A filter instance can be constructed by name from one of
# the provided resource files describing specific instruments or standard wavelength bands. Alternatively,
# a filter can be created with a uniform transmission curve over a certain wavelength range.
#
# The precise formalae involved in the integration over the filter and the calculation of the pivot wavelength
# depend on whether the instrument counts photons (photon counter) or measures energy (bolometer).
#
class Filter:
    # ---------- Constructing -------------------------------------

    ## The constructor constructs a new Filter instance in one of two ways, depending on the type of the argument.
    #
    # If \em filterspec is a tuple with two numbers, a bolometer-type filter is constructed with a uniform
    # transmission curve in the indicated range. The (min,max) wavelength values must be expressed in micron.
    #
    # If \em filterspec is a string, the constructor locates a "$PYTHONPATH/dat/filter" directory, and searches
    # in that directory for a VOTable resource file with a name that contains the specified string. If there is
    # exactly one such file, its contents is loaded into the new Filter instance. If there is no such file, or
    # if the specified string matches multiple file names, the constructor raises an error.
    #
    # The VOTable resource files do not properly specify the type of filter (photon counter or bolometer).
    # The constructor classifies filters based on the instrument name in the filter identifier.
    # While this simple heuristic works for now, it is not guaranteed to be correct for future cases.
    #
    # The tables below list the filters available at the time of writing. More filter definitions can be
    # downloaded from the filter profile service web site http://svo2.cab.inta-csic.es/theory/fps.
    #
    #| Filter-spec | \f$\lambda_\text{ctr}\f$ | \f$\lambda_\text{piv}\f$ | Type | Description
    #|------- -----|--------------------------|--------------------------|------|------------
    #| 2MASS/2MASS.J | 1.239 | 1.239 | Pho | 2MASS J
    #| 2MASS/2MASS.H | 1.649 | 1.649 | Pho | 2MASS H
    #| 2MASS/2MASS.Ks | 2.163 | 2.164 | Pho | 2MASS Ks
    #| GALEX/GALEX.FUV | 0.1526 | 0.1535 | Pho | GALEX FUV
    #| GALEX/GALEX.NUV | 0.2329 | 0.2301 | Pho | GALEX NUV
    #| Herschel/Pacs.blue | 71.33 | 70.77 | Bol | Herschel Pacs blue filter
    #| Herschel/Pacs.green | 102.3 | 100.8 | Bol | Herschel Pacs green filter
    #| Herschel/Pacs.red | 166.1 | 161.9 | Bol | Herschel Pacs red filter
    #| Herschel/SPIRE.PSW_ext | 257.6 | 252.5 | Bol | Herschel SPIRE PSW filter (extended sources)
    #| Herschel/SPIRE.PMW_ext | 357.5 | 354.3 | Bol | Herschel SPIRE PMW filter (extended sources)
    #| Herschel/SPIRE.PLW_ext | 518.4 | 515.4 | Bol | Herschel SPIRE PLW filter (extended sources)
    #| IRAS/IRAS.12mu | 11.43 | 11.41 | Pho | IRAS 12 micron
    #| IRAS/IRAS.25mu | 23.97 | 23.61 | Pho | IRAS 25 micron
    #| IRAS/IRAS.60mu | 61.88 | 60.41 | Pho | IRAS 60 micron
    #| IRAS/IRAS.100mu | 101 | 101.1 | Pho | IRAS 100 micron
    #| Misc/MCPS.U | 0.3666 | 0.3646 | Pho | MCPS Johnson U
    #| Misc/MCPS.B | 0.4406 | 0.4434 | Pho | MCPS Johnson B
    #| Misc/MCPS.V | 0.5426 | 0.5493 | Pho | MCPS Johnson V
    #| Misc/MCPS.I | 0.8643 | 0.8739 | Pho | MCPS Johnson I
    #| SLOAN/SDSS.u | 0.3565 | 0.3557 | Pho | SDSS u
    #| SLOAN/SDSS.g | 0.47 | 0.4702 | Pho | SDSS g
    #| SLOAN/SDSS.r | 0.6174 | 0.6176 | Pho | SDSS r
    #| SLOAN/SDSS.i | 0.7534 | 0.749 | Pho | SDSS i
    #| SLOAN/SDSS.z | 0.8782 | 0.8947 | Pho | SDSS z
    #| Spitzer/IRAC.I1 | 3.547 | 3.551 | Pho | IRAC I1
    #| Spitzer/IRAC.I2 | 4.502 | 4.496 | Pho | IRAC I2
    #| Spitzer/IRAC.I3 | 5.716 | 5.724 | Pho | IRAC I3
    #| Spitzer/IRAC.I4 | 7.856 | 7.884 | Pho | IRAC I4
    #| Spitzer/MIPS.24mu | 23.47 | 23.59 | Bol | MIPS 24 microns
    #| Spitzer/MIPS.70mu | 70.52 | 70.89 | Bol | MIPS 70 microns
    #| Spitzer/MIPS.160mu | 156.9 | 155.4 | Bol | MIPS 160 microns
    #| UKIRT/UKIDSS.Z | 0.8825 | 0.8826 | Pho | UKIDSS Z
    #| UKIRT/UKIDSS.Y | 1.03 | 1.031 | Pho | UKIDSS Y
    #| UKIRT/UKIDSS.J | 1.249 | 1.25 | Pho | UKIDSS J
    #| UKIRT/UKIDSS.H | 1.638 | 1.635 | Pho | UKIDSS H
    #| UKIRT/UKIDSS.K | 2.206 | 2.206 | Pho | UKIDSS K
    #| WISE/WISE.W1 | 3.466 | 3.39 | Pho | WISE W1 filter
    #| WISE/WISE.W2 | 4.644 | 4.641 | Pho | WISE W2 filter
    #| WISE/WISE.W3 | 13.22 | 12.57 | Pho | WISE W3 filter
    #| WISE/WISE.W4 | 22.22 | 22.31 | Pho | WISE W4 filter
    #
    def __init__(self, filterspec):

        # string --> load from appropriate resource file
        if isinstance(filterspec, types.StringTypes):
            for pythondir in sys.path:
                filterdir = os.path.join(pythondir, "pts", "core", "dat", "filters")
                if os.path.isdir(filterdir):
                    filterfiles = filter(lambda fn: fn.endswith(".xml") and filterspec in fn, os.listdir(filterdir))
                    if len(filterfiles) > 1: raise ValueError("filter spec " + filterspec + " is ambiguous")
                    if len(filterfiles) < 1: raise ValueError("no filter found with spec " + filterspec)

                    # load the XML tree
                    tree = etree.parse(os.path.join(filterdir, filterfiles[0]))

                    # verify the wavelength unit to be Angstrom
                    unit = tree.xpath("//RESOURCE/PARAM[@name='WavelengthUnit'][1]/@value")[0]
                    if unit!='Angstrom': raise ValueError("VOTable uses unsupported unit: " + unit)

                    # load some basic properties (converting from Angstrom to micron)
                    self._WavelengthMin = 1e-4*float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthMin'][1]/@value")[0])
                    self._WavelengthMax = 1e-4*float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthMax'][1]/@value")[0])
                    self._WavelengthCen = 1e-4*float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthCen'][1]/@value")[0])
                    self._FilterID = tree.xpath("//RESOURCE/PARAM[@name='filterID'][1]/@value")[0]
                    self._Description = tree.xpath("//RESOURCE/PARAM[@name='Description'][1]/@value")[0]
                    self._Description = self._Description.replace("&#956;m", "micron")

                    # load the transmission table (converting wavelengths from Angstrom to micron)
                    values = np.array(tree.xpath("//RESOURCE/TABLE/DATA/TABLEDATA[1]/TR/TD/text()"), dtype=float)
                    if len(values)<4: raise ValueError("transmission table not found in filter definition")
                    self._Wavelengths,self._Transmission = np.reshape(values, (-1,2)).T
                    self._Wavelengths *= 1e-4

                    # determine the filter type (there seems to be no better heuristic than using the instrument name)
                    self._PhotonCounter = not any(["/"+x in self._FilterID.lower() for x in ("mips","pacs","spire")])

                    # calculate the pivot wavelength
                    if self._PhotonCounter:
                        integral1 = np.trapz(x=self._Wavelengths, y=self._Transmission*self._Wavelengths)
                        integral2 = np.trapz(x=self._Wavelengths, y=self._Transmission/self._Wavelengths)
                    else:
                        integral1 = np.trapz(x=self._Wavelengths, y=self._Transmission)
                        integral2 = np.trapz(x=self._Wavelengths, y=self._Transmission/(self._Wavelengths**2))
                    self._IntegratedTransmission = integral1
                    self._WavelengthPivot = np.sqrt(integral1/integral2)

                    # report success
                    return

            raise ValueError("filter resource path not found")

        # range --> construct ad hoc uniform bolometer
        else:
            self._WavelengthMin, self._WavelengthMax = map(float,filterspec)
            self._WavelengthCen = 0.5 * (self._WavelengthMin + self._WavelengthMax)
            self._FilterID = "Uniform_[{},{}]".format(self._WavelengthMin,self._WavelengthMax)
            self._Description = "Uniform filter in range [{},{}]".format(self._WavelengthMin,self._WavelengthMax)
            self._Wavelengths = np.array((self._WavelengthMin,self._WavelengthMax))
            self._Transmission = np.ones((2,))
            self._PhotonCounter = False
            self._IntegratedTransmission = self._WavelengthMax - self._WavelengthMin
            self._WavelengthPivot = np.sqrt(self._WavelengthMin * self._WavelengthMax)

    # ---------- Retrieving information -------------------------------

    ## This function returns a unique identifier for the filter.
    def filterID(self):
        return self._FilterID

    ## This function returns a human-readable description for the filter.
    def description(self):
        return self._Description

    ## This function returns the minimum wavelength for the filter, in micron.
    def minwavelength(self):
        return self._WavelengthMin

    ## This function returns the maximum wavelength for the filter, in micron.
    def maxwavelength(self):
        return self._WavelengthMax

    ## This function returns the center wavelength for the filter, in micron. The center wavelength is
    # defined as the wavelength halfway between the two points for which filter response or transmission
    # (depending on the filter type) is half maximum.
    def centerwavelength(self):
        return self._WavelengthCen

    ## This function returns the pivot wavelength for the filter, in micron. The pivot wavelength is defined
    # as the wavelength that connects the filter-averaged wavelength and frequency-style fluxes through
    # \f$\left<F_\nu\right> = \left<F_\lambda\right>\lambda_\text{pivot}^2/c\f$. The value depends
    # on the filter type. For a photon counter with response curve \f$R(\lambda)\f$,
    # \f[ \lambda_\mathrm{pivot} = \sqrt{ \frac{ \int\lambda R(\lambda) \,\mathrm{d}\lambda }
    #     {  \int R(\lambda) \,\mathrm{d}\lambda/\lambda } }. \f]
    # For a bolometer with transmission curve \f$T(\lambda)\f$,
    # \f[ \lambda_\mathrm{pivot} = \sqrt{ \frac{ \int T(\lambda) \,\mathrm{d}\lambda }
    #     { \int T(\lambda) \,\mathrm{d}\lambda/\lambda^2 } }. \f]
    def pivotwavelength(self):
        return self._WavelengthPivot

    # ---------- Integrating --------------------------------------

    ## This function calculates and returns the filter-averaged value \f$\left<F_\lambda\right>\f$ for a given
    # spectral energy distribution \f$F_\lambda(\lambda)\f$. The calculation depends
    # on the filter type. For a photon counter with response curve \f$R(\lambda)\f$,
    # \f[ \left<F_\lambda\right> = \frac{ \int\lambda F_\lambda(\lambda)R(\lambda) \,\mathrm{d}\lambda }
    #     { \int\lambda R(\lambda) \,\mathrm{d}\lambda }. \f]
    # For a bolometer with transmission curve \f$T(\lambda)\f$,
    # \f[ \left<F_\lambda\right> = \frac{ \int F_\lambda(\lambda)T(\lambda) \,\mathrm{d}\lambda }
    #     { \int T(\lambda) \,\mathrm{d}\lambda }. \f]
    #
    # The quantities \f$F_\lambda(\lambda)\f$ must be expressed per unit of wavelength (and \em not, for example,
    # per unit of frequency). The resulting \f$\left<F_\lambda\right>\f$ has the same units as the input distribition.
    # \f$F_\lambda(\lambda)\f$ can be expressed in any units (as long as it is per unit of wavelength) and it can
    # represent various quantities; for example a flux density, a surface density, or a luminosity density.
    #
    # The function accepts two arguments:
    # - \em wavelengths: a numpy array specifying the wavelengths \f$\lambda_\ell\f$, in micron, in increasing order,
    #   on which the spectral energy distribution is sampled. The integration is performed on a wavelength grid that
    #   combines the grid points given here with the grid points on which the filter response or transmission curve
    #   is defined.
    # - \em densities: a numpy array specifying the spectral energy distribution(s) \f$F_\lambda(\lambda_\ell)\f$
    #   per unit of wavelength. This can be an array with the same length as \em wavelengths, or a multi-dimensional
    #   array where the last dimension has the same length as \em wavelengths.
    #   The returned result will have the shape of \em densities minus the last (or only) dimension.
    def convolve(self, wavelengths, densities):
        # define short names for the involved wavelength grids
        wa = wavelengths
        wb = self._Wavelengths

        # create a combined wavelength grid, restricted to the overlapping interval
        w1 = wa[ (wa>=wb[0]) & (wa<=wb[-1]) ]
        w2 = wb[ (wb>=wa[0]) & (wb<=wa[-1]) ]
        w = np.unique(np.hstack((w1,w2)))
        if len(w) < 2: return 0

        # log-log interpolate SED and transmission on the combined wavelength grid
        # (use scipy interpolation function for SED because np.interp does not support broadcasting)
        F = np.exp(interp1d(np.log(wa), _log(densities), copy=False, bounds_error=False, fill_value=0.)(np.log(w)))
        T = np.exp(np.interp(np.log(w), np.log(wb), _log(self._Transmission), left=0., right=0.))

        # perform the integration
        if self._PhotonCounter:
            return np.trapz(x=w, y=w*F*T) / self._IntegratedTransmission
        else:
            return np.trapz(x=w, y=F*T) / self._IntegratedTransmission

## This private helper function returns the natural logarithm for positive values, and a large negative number
# (but not infinity) for zero or negative values.
def _log(X):
    zeromask = X<=0
    logX = np.empty(X.shape)
    logX[zeromask] = -750.  # the smallest (in magnitude) negative value x for which np.exp(x) returns zero
    logX[~zeromask] = np.log(X[~zeromask])
    return logX

# -----------------------------------------------------------------
