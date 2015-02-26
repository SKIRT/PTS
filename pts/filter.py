#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.filter Working with wavelength filters (or bands).
#
# An instance of the Filter class in this module represents a particular wavelength band, including its transmission
# curve, and allows integrating a given spectrum over the band. A filter instance can be constructed from one of
# the provided resource files describing specific instruments or standard wavelength bands. Alternatively,
# a filter can be created with a uniform transmission curve over a certain wavelength range.
#

# -----------------------------------------------------------------

import os
import os.path
import sys
import types
import numpy as np
from lxml import etree

# -----------------------------------------------------------------
#  Filter class
# -----------------------------------------------------------------

## An instance of the Filter class represents a particular wavelength band, including its transmission curve
# and some basic properties such as its mean and effective wavelengths. The class provides a function to
# integrate a given spectrum over the band. A filter instance can be constructed by name from one of
# the provided resource files describing specific instruments or standard wavelength bands. Alternatively,
# a filter can be created with a uniform transmission curve over a certain wavelength range.
#
class Filter:
    # ---------- Constructing -----------------------------

    ## The constructor constructs a new Filter instance in one of two ways, depending on the type of the argument.
    #
    # If \em filterspec is a string, the constructor locates a "$PYTHONPATH/dat/filter" directory, and searches
    # in that directory for a VOTable resource file with a name that contains the specified string. If there is
    # exactly one such file, its contents is loaded into the new Filter instance. If there is no such file, or
    # if the specified string matches multiple file names, the constructor raises an error. The table below
    # lists the filters available at the time of writing. More filter definitions can be downloaded from
    # the filter profile service web site http://svo2.cab.inta-csic.es/theory/fps .
    #
    #| Filterspec | Description
    #|--------|------------
    #| GALEX.FUV | GALEX FUV
    #| GALEX.NUV | GALEX NUV
    #| Pacs.blue | Herschel Pacs blue filter
    #| Pacs.green | Herschel Pacs green filter
    #| Pacs.red | Herschel Pacs red filter
    #| SPIRE.PLW | Herschel SPIRE PLW filter (extended sources)
    #| SPIRE.PMW | Herschel SPIRE PMW filter (extended sources)
    #| SPIRE.PSW | Herschel SPIRE PSW filter (extended sources)
    #| IRAS.100 | IRAS 100 micron
    #| IRAS.12 | IRAS 12 micron
    #| IRAS.25 | IRAS 25 micron
    #| IRAS.60 | IRAS 60 micron
    #| MCPS.B | MCPS Johnson B
    #| MCPS.I | MCPS Johnson I
    #| MCPS.U | MCPS Johnson U
    #| MCPS.V | MCPS Johnson V
    #| SDSS.g | SDSS g
    #| SDSS.i | SDSS i
    #| SDSS.r | SDSS r
    #| SDSS.u | SDSS u
    #| SDSS.z | SDSS z
    #| UKIDSS.H | UKIDSS H
    #| UKIDSS.J | UKIDSS J
    #| UKIDSS.K | UKIDSS K
    #| UKIDSS.Y | UKIDSS Y
    #| UKIDSS.Z | UKIDSS Z
    #
    # If \em filterspec is a tuple with two numbers, a filter is constructed with a uniform transmission curve
    # in the indicated range. The (min,max) wavelength values must be expressed in micron.
    #
    def __init__(self, filterspec):

        # string --> load from appropriate resource file
        if isinstance(filterspec, types.StringTypes):
            for pythondir in sys.path:
                filterdir = os.path.join(pythondir, "dat", "filters")
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
                    self._WavelengthMean = 1e-4*float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthMean'][1]/@value")[0])
                    self._FilterID = tree.xpath("//RESOURCE/PARAM[@name='filterID'][1]/@value")[0]
                    self._Description = tree.xpath("//RESOURCE/PARAM[@name='Description'][1]/@value")[0]
                    self._Description = self._Description.replace("&#956;m", "micron")

                    # load the transmission table (converting wavelengths from Angstrom to micron)
                    values = np.array(tree.xpath("//RESOURCE/TABLE/DATA/TABLEDATA[1]/TR/TD/text()"), dtype=float)
                    if len(values)<4: raise ValueError("transmission table not found in filter definition")
                    self._wavelengths,self._transmission = np.reshape(values, (-1,2)).T
                    self._wavelengths *= 1e-4

                    # calculate the effective width since the value in the file seems to be nonsense
                    self._WidthEff = np.trapz(x=self._wavelengths, y=self._transmission)
                    return

            raise ValueError("filter resource path not found")

        # range --> construct ad hoc
        else:
            self._WavelengthMin, self._WavelengthMax = map(float,filterspec)
            self._WavelengthCen = 0.5 * (self._WavelengthMin + self._WavelengthMax)
            self._WavelengthMean = self._WavelengthCen
            self._WidthEff = self._WavelengthMax - self._WavelengthMin
            self._Description = "Uniform filter in range [{},{}]".format(self._WavelengthMin,self._WavelengthMax)
            self._FilterID = "Uniform_[{},{}]".format(self._WavelengthMin,self._WavelengthMax)
            self._wavelengths = np.array((self._WavelengthMin,self._WavelengthMax))
            self._transmission = np.ones((2,))

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
    # defined as the wavelength halfway between the two points for which filter transmission is half maximum.
    def centerwavelength(self):
        return self._WavelengthCen

    ## This function returns the mean wavelength for the filter, in micron. The mean wavelength is
    # defined as \f$\int \lambda\,T(\lambda)\,\mathrm{d}\lambda\; / \int T(\lambda)\,\mathrm{d}\lambda\f$
    # where \f$T(\lambda)\f$ is the filter's transmission curve.
    def meanwavelength(self):
        return self._WavelengthMean

    ## This function returns the effective width for the filter, in micron. The effective width is defined as
    #  \f$\int T(\lambda)\,\mathrm{d}\lambda\f$ where \f$T(\lambda)\f$ is the filter's transmission curve.
    # The effective width is equivalent to the horizontal size of a rectangle with height equal to 100% transmission
    # and with the same area as the one covered by the filter transmission curve.
    def effectivewidth(self):
        return self._WidthEff

# -----------------------------------------------------------------
