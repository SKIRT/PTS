#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.skifile Reading and updating a SKIRT parameter file.
#
# An instance of the SkiFile class in this module allows reading from and updating an existing ski file.

# -----------------------------------------------------------------

import os.path
from datetime import datetime
from lxml import etree

# -----------------------------------------------------------------
#  Utility functions
# -----------------------------------------------------------------

## This private helper function converts a string containing a decimal number followed by a units string
# to an appropriate number. The current implementation requires that the units string equals "micron"
# and returns the number unchanged.
def _tonumber(stringwithunits):
    segments = stringwithunits.split()
    if len(segments) != 2: raise ValueError("Invalid number string")
    if segments[1] != "micron": raise ValueError("Unsupported units")
    return float(segments[0])

# -----------------------------------------------------------------
#  SkiFile class
# -----------------------------------------------------------------

## An instance of the SkiFile class represents a particular existing SKIRT parameter file (\em ski file).
# There are functions to read and/or update certain information in the ski file, such as obtaining or setting
# the value of a particular parameter. The intention is to encapsulate any knowledge about the ski file format
# and structure within this class, concentrating the update pain if and when that format changes.
# Consequently the public functions in this class are quite high-level, and specific rather than generic.
#
# Updates made to a SkiFile instance do \em not affect the underlying file; use the saveto() function to save
# the updated contents of a SkiFile instance to another file (or to replace the original file if so desired).
#
# A SkiFile class instance is always constructed from an existing ski file; creating a new ski file from scratch
# is not supported. To create a new ski file, start SKIRT in interactive mode (without any arguments).
#
class SkiFile:
    # ---------- Constructing and saving -----------------------------

    ## The constructor loads the contents of the specified ski file into a new SkiFile instance.
    # The filename \em must end with ".ski" or with "_parameters.xml".
    #
    def __init__(self, filepath):
        if not filepath.lower().endswith((".ski","_parameters.xml")):
            raise ValueError("Invalid filename extension for ski file")
        # load the XML tree (remove blank text to avoid confusing the pretty printer when saving)
        self.tree = etree.parse(os.path.expanduser(filepath), parser=etree.XMLParser(remove_blank_text=True))

    ## This function saves the (possibly updated) contents of the SkiFile instance into the specified file.
    # The filename \em must end with ".ski". Saving to and thus replacing the ski file from which this
    # SkiFile instance was originally constructed is allowed, but often not the intention.
    def saveto(self, filepath):
        if not filepath.lower().endswith(".ski"):
            raise ValueError("Invalid filename extension for ski file")
        # update the producer and time attributes on the root element
        root = self.tree.getroot()
        root.set("producer","Python Toolkit for SKIRT (SkiFile class)")
        root.set("time", datetime.now().strftime("%Y-%m-%dT%H:%M:%S"))
        # serialize the XML tree
        outfile = open(os.path.expanduser(filepath), "wb")
        outfile.write(etree.tostring(self.tree, encoding="UTF-8", xml_declaration=True, pretty_print=True))
        outfile.close()

    # ---------- Retrieving information -------------------------------

    ## This function returns a list of the wavelengths specified in the ski file for an oligochromatic simulation,
    # in micron. If the ski file specifies a panchromatic simulation, the function returns an empty list.
    # The current implementation requires that the wavelengths in the ski file are specified in micron.
    def wavelengths(self):
        # get the value of the wavelengths attribute on the OligoWavelengthGrid element (as a list of query results)
        results = self.tree.xpath("//OligoWavelengthGrid/@wavelengths")
        # if not found, return an empty list
        if len(results) != 1: return [ ]
        # split the first result in separate strings, check the units, and extract the numbers
        return [_tonumber(s) for s in results[0].split(",")]

    ## This function returns the first instrument's distance in parsec. The current implementation requires that
    # the distance attribute in the ski file specifies one of the following unit strings: pc, kpc, Mpc.
    def instrumentdistance(self):
        # get the first instrument element
        instruments = self.tree.xpath("//instruments/*[1]")
        if len(instruments) != 1: raise ValueError("No instruments in ski file")
        # get the distance and the unit string
        distance, unit = instruments[0].get("distance").split()
        # convert to pc
        units = { 'pc': 1.,  'kpc': 1e3,  'Mpc': 1e6 }
        return float(distance)*units[unit]

    ## This function returns the shape of the first instrument's frame, in pixels.
    def instrumentshape(self):
        # get the first instrument element
        instruments = self.tree.xpath("//instruments/*[1]")
        if len(instruments) != 1: raise ValueError("No instruments in ski file")
        # get its shape
        return ( int(instruments[0].get("pixelsX")), int(instruments[0].get("pixelsY")) )

    ## This function returns a list of instrument names, in order of occurrence in the ski file.
    def instrumentnames(self):
        # get the instrument elements
        instruments = self.tree.xpath("//instruments/*")
        # return their names
        return [ instr.get("instrumentName") for instr in instruments ]

    # ---------- Updating information ---------------------------------

    ## This function applies an XSLT transform to the ski file if an XPath condition evaluates to true.
    # The first argument is a string specifying an XPath 1.0 expression to be evaluated in the context of the XML
    # document representing the ski file; the expression value is converted to boolean according to XPath semantics.
    # If the value is true, the XSLT 1.0 transform specified in the second argument is applied to the XML document,
    # and the result replaces the original document. The second argument is a string containing one or more
    # \<xsl:template\> elements that specify the changes to be applied to the document. The \<xsl:stylesheet\>
    # element and the identity template are automatically added and must not be contained in the argument string.
    # The function returns true if the transform was applied, and false if it was not (i.e. the document is unchanged).
    def transformif(self, condition, templates):
        needed = self.tree.xpath("boolean(" + condition + ")")
        if needed:
            prefix  = '''<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
                           <xsl:template match="@*|node()">
                             <xsl:copy>
                               <xsl:apply-templates select="@*|node()"/>
                             </xsl:copy>
                           </xsl:template>'''
            postfix = '''</xsl:stylesheet>'''
            transform = etree.XSLT(etree.XML(prefix + templates + postfix))
            self.tree = transform(self.tree)
        return needed

    ## This function sets the number of photon packages on the MonteCarloSimulation element in the ski file
    # to the specified value
    def setpackages(self, number):
        # get the MonteCarloSimulation element
        elems = self.tree.xpath("//OligoMonteCarloSimulation | //PanMonteCarloSimulation")
        if len(elems) != 1: raise ValueError("No MonteCarloSimulation in ski file")
        # set the attribute value
        elems[0].set("packages", str(number))

    ## This function replaces any instruments in the ski file by a new list of perspective instruments
    # corresponding to the movie frames defined in the specified list. The instruments are named "0",
    # "1", "2"... corresponding to the zero-based frame index in the list. Each frame is given as a tuple
    # containing the following information: viewport shape (in pixels), viewport size, viewport position,
    # crosshair position, upwards position, and focal length (all in world coordinates, expressed in the
    # default units for length in the target ski file).
    # The components of each item are grouped in tuples, so the structure of the complete list is:
    # [ ((Nx,Ny),(Sx,Sy),(Vx,Vy,Vz),(Cx,Cy,Cz),(Ux,Uy,Uz),Fe) , ... ]
    def setperspectiveinstruments(self, frames):
        # get the instruments element
        parents = self.tree.xpath("//instruments")
        if len(parents) == 0: raise ValueError("No 'instruments' element in ski file")
        if len(parents) > 1: raise ValueError("Multiple 'instruments' elements in ski file")
        parent = parents[0]
        # remove the old instruments
        for instrument in parent.getchildren():
            parent.remove(instrument)
        # add a new instrument for each frame
        index = 0
        for pixels,size,view,cross,up,focal in frames:
            attrs = { "instrumentName" : str(index),
                      "pixelsX" : str(pixels[0]), "pixelsY" : str(pixels[1]), "width" : str(size[0]),
                      "viewX" : str(view[0]),  "viewY" : str(view[1]), "viewZ" : str(view[2]),
                      "crossX" : str(cross[0]), "crossY" : str(cross[1]), "crossZ" : str(cross[2]),
                      "upX" : str(up[0]), "upY" : str(up[1]),  "upZ" : str(up[2]), "focal" : str(focal) }
            parent.append(parent.makeelement("PerspectiveInstrument", attrs))
            index += 1

    ## This function sets the filename attribute of the SPHStellarComp element to the specified value.
    def setstarfile(self, filename):
        # get the SPHStellarComp element
        elems = self.tree.xpath("//SPHStellarComp")
        if len(elems) != 1: raise ValueError("No SPHStellarComp in ski file")
        # set the attribute value
        elems[0].set("filename", filename)

    ## This function sets the filename attribute of the SPHStarburstComp element to the specified value.
    def sethiifile(self, filename):
        # get the SPHStarburstComp element
        elems = self.tree.xpath("//SPHStarburstComp")
        if len(elems) != 1: raise ValueError("No SPHStarburstComp in ski file")
        # set the attribute value
        elems[0].set("filename", filename)

    ## This function sets the filename attribute of the SPHDustDistribution element to the specified value.
    def setgasfile(self, filename):
        # get the SPHDustDistribution element
        elems = self.tree.xpath("//SPHDustDistribution")
        if len(elems) != 1: raise ValueError("No SPHDustDistribution in ski file")
        # set the attribute value
        elems[0].set("filename", filename)

    ## This function sets any extentX, extentY and extentZ attributes to the specified value (converted to a string),
    # regardless of the element in which such attributes reside.
    def setextent(self, value):
        strvalue = str(value)
        for attr in self.tree.xpath("//*/@extentX"): attr.getparent().set("extentX", strvalue)
        for attr in self.tree.xpath("//*/@extentY"): attr.getparent().set("extentY", strvalue)
        for attr in self.tree.xpath("//*/@extentZ"): attr.getparent().set("extentZ", strvalue)

# -----------------------------------------------------------------
