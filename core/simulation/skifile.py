#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.skifile Reading and updating a SKIRT parameter file.
#
# An instance of the SkiFile class in this module allows reading from and updating an existing ski file.

# -----------------------------------------------------------------

# Import standard modules
import os.path
from datetime import datetime
from lxml import etree
from numpy import arctan

# Import the relevant PTS classes and modules
from .units import SkirtUnits
from ..basics.filter import Filter
from ..tools import archive as arch

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

        # Set the path to the ski file
        self.path = os.path.expanduser(filepath)

        # load the XML tree (remove blank text to avoid confusing the pretty printer when saving)
        self.tree = etree.parse(arch.opentext(self.path), parser=etree.XMLParser(remove_blank_text=True))

        # Replace path by the full, absolute path
        self.path = os.path.abspath(self.path)

    ## This function saves the (possibly updated) contents of the SkiFile instance into the specified file.
    # The filename \em must end with ".ski". Saving to and thus replacing the ski file from which this
    # SkiFile instance was originally constructed is allowed, but often not the intention.
    def saveto(self, filepath):
        if not filepath.lower().endswith(".ski"):
            raise ValueError("Invalid filename extension for ski file")
        # update the producer and time attributes on the root element
        root = self.tree.getroot()
        root.set("producer", "Python Toolkit for SKIRT (SkiFile class)")
        root.set("time", datetime.now().strftime("%Y-%m-%dT%H:%M:%S"))
        # serialize the XML tree
        outfile = open(os.path.expanduser(filepath), "wb")
        outfile.write(etree.tostring(self.tree, encoding="UTF-8", xml_declaration=True, pretty_print=True))
        outfile.close()

    ## This function saves the ski file to the original path
    def save(self): self.saveto(self.path)

    # ---------- Retrieving information -------------------------------

    ## This function returns a SkirtUnits object initialized with the SKIRT unit system ('SI', 'stellar', or
    # 'extragalactic') and the flux style ('neutral', 'wavelength' or 'frequency') specified in the ski file.
    def units(self):
        unitelements = self.tree.xpath("//units/*[1]")
        if len(unitelements) == 1:
            unitsystem = unitelements[0].tag
            fluxstyle = unitelements[0].get("fluxOutputStyle", default='neutral')
        else:
            unitsystem = 'extragalactic'
            fluxstyle = 'neutral'
        return SkirtUnits(unitsystem, fluxstyle)

    ## This function returns the number of wavelengths for oligochromatic or panchromatic simulations
    def nwavelengths(self):
        # Try to get the list of wavelengths from the ski file
        wavelengths = self.wavelengths()
        # If the list is not empty, retun its size
        if wavelengths: return len(wavelengths)
        # If the list is empty, the ski file either represents a panchromatic simulation (and we can get the
        # number of points directly from the tree) or a FileWavelengthGrid is used (in which case we raise an error)
        entry = self.tree.xpath("//wavelengthGrid/*[1]")[0]
        if entry.tag == 'FileWavelengthGrid':
            raise ValueError("The number of wavelengths is not defined within the ski file. Call wavelengthsfile().")
        else:
            return int(entry.get("points"))

    ## This function returns the name of the wavelengths file that is used for the simulation, if any
    def wavelengthsfile(self):
        entry = self.tree.xpath("//FileWavelengthGrid")
        if entry: return entry[0].get("filename")
        else: return None

    ## This function returns the number of photon packages per wavelength
    def packages(self):
        # Get the MonteCarloSimulation element
        elems = self.tree.xpath("//OligoMonteCarloSimulation | //PanMonteCarloSimulation")
        if len(elems) != 1: raise ValueError("No MonteCarloSimulation in ski file")
        # Get the number of packages
        return int(float(elems[0].get("packages")))

    ## This function returns the number of dust cells
    def ncells(self):

        xpoints = self.nxcells()
        ypoints = 1
        zpoints = 1

        try:
            ypoints = self.nycells()
        except ValueError: pass

        try:
            zpoints = self.nzcells()
        except ValueError: pass

        # Return the total number of dust cells
        return xpoints*ypoints*zpoints

    ## This function returns the number of dust cells in the x direction
    def nxcells(self):
        try:
            xpoints = int(self.tree.xpath("//meshX/*")[0].get("numBins"))
        except TypeError:
            raise ValueError("The number of dust cels is not defined within the ski file")
        return xpoints

    ## This function returns the number of dust cells in the y direction
    def nycells(self):
        try:
            ypoints = int(self.tree.xpath("//meshY/*")[0].get("numBins"))
        except TypeError:
            raise ValueError("The dimension of the dust grid is lower than 2")
        return ypoints

    ## This function returns the number of dust cells in the z direction
    def nzcells(self):
        try:
            zpoints = int(self.tree.xpath("//meshZ/*")[0].get("numBins"))
        except TypeError:
            raise ValueError("The dimension of the dust grid is lower than 3")
        return zpoints

    ## This function returns the dimension of the dust grid
    def dimension(self):
        # Try to find the number of points in the y direction
        try:
            int(self.tree.xpath("//dustGridStructure/*[1]")[0].get("pointsY"))
        except TypeError:
            return 1

        # Try to find the number of points in the z direction
        try:
            int(self.tree.xpath("//dustGridStructure/*[1]")[0].get("pointsZ"))
        except TypeError:
            return 2

        # If finding the number of ypoints and zpoints succeeded, the grid is 3-dimensional
        return 3

    ## This function returns the number of dust components
    def ncomponents(self):
        components = self.tree.xpath("//CompDustDistribution/components/*")
        return int(len(components))

    ## This function returns the number of dust library items
    def nlibitems(self):
        dustlib = self.tree.xpath("//dustLib/*")[0]
        if dustlib.tag == "AllCellsDustLib":
            return self.ncells()
        elif dustlib.tag == "Dim2DustLib":
            return dustlib.attrib["pointsTemperature"] * dustlib.attrib["pointsWavelength"]
        elif dustlib.tag == "Dim1DustLib":
            return dustlib.attrib["entries"]

    ## This function returns the number of dust populations (from all dust mixes combined)
    def npopulations(self):
        npops = 0
        # For each dust mix
        for dustmix in self.tree.xpath("//mix/*[1]"):
            if dustmix.tag in ["InterstellarDustMix", "Benchmark1DDustMix", "Benchmark2DDustMix", "DraineLiDustMix"]:
                npops += 1
            elif dustmix.tag == "TrustDustMix":
                npops += int(dustmix.attrib["graphitePops"])
                npops += int(dustmix.attrib["silicatePops"])
                npops += int(dustmix.attrib["PAHPops"])
            elif dustmix.tag == "ConfigurableDustMix":
                npops += len(self.tree.xpath("//ConfigurableDustMix/populations/*"))
        return npops

    ## This function returns the number of simple instruments
    def nsimpleinstruments(self):
        return len(self.tree.xpath("//SimpleInstrument"))

    ## This function returns the number of full instruments
    def nfullinstruments(self):
        return len(self.tree.xpath("//FullInstrument"))

    ## This function returns whether transient heating is enabled
    def transientheating(self):
        return len(self.tree.xpath("//TransientDustEmissivity"))

    ## This function returns whether dust emission is enabled
    def dustemission(self):
        return len(self.tree.xpath("//dustEmissivity"))

    @property
    def emission_boost(self):
        try:
            pandustsystem = self.tree.xpath("//PanDustSystem")[0]
            return float(pandustsystem.attrib["emissionBoost"])
        except:
            raise ValueError("Not a panchromatic simulation")

    ## This function returns whether dust selfabsorption is enabled
    def dustselfabsorption(self):
        try:
            pandustsystem = self.tree.xpath("//PanDustSystem")[0]
            return (pandustsystem.attrib["selfAbsorption"] == "true")
        except:
            return False

    ## This function returns the number of pixels for each of the instruments
    def npixels(self, nwavelengths=None):
        pixels = []
        nwavelengths = nwavelengths if nwavelengths is not None else self.nwavelengths()
        instruments = self.tree.xpath("//instruments/*")
        for instrument in instruments:
            type = instrument.tag
            name = instrument.attrib["instrumentName"]
            datacube = int(instrument.attrib["pixelsX"])*int(instrument.attrib["pixelsY"])*nwavelengths
            if type == "SimpleInstrument":
                pixels.append([name, type, datacube])
            elif type == "FullInstrument":
                scattlevels = int(instrument.attrib["scatteringLevels"])
                scattering = scattlevels + 1 if scattlevels > 0 else 0
                dustemission = 1 if self.dustemission() else 0
                npixels = datacube * (3 + scattering + dustemission)
                pixels.append([name, type, npixels])
        return pixels

    ## This function returns a list of the wavelengths specified in the ski file for an oligochromatic simulation,
    # in micron. If the ski file specifies a panchromatic simulation, the function returns an empty list.
    # The current implementation requires that the wavelengths in the ski file are specified in micron.
    def wavelengths(self):
        # get the value of the wavelengths attribute on the OligoWavelengthGrid element (as a list of query results)
        results = self.tree.xpath("//OligoWavelengthGrid/@wavelengths")
        # if not found, return an empty list
        if len(results) != 1: return []
        # split the first result in separate strings, extract the numbers using the appropriate units
        units = self.units()
        return [units.convert(s,to_unit='micron',quantity='wavelength') for s in results[0].split(",")]

    ## This function returns the first instrument's distance, in the specified units (default is 'pc').
    def instrumentdistance(self, unit='pc'):
        # get the first instrument element
        instruments = self.tree.xpath("//instruments/*[1]")
        if len(instruments) != 1: raise ValueError("No instruments in ski file")
        # get the distance including the unit string
        distance = instruments[0].get("distance")
        # convert to requested units
        return self.units().convert(distance, to_unit=unit, quantity='distance')

    ## This function returns the shape of the first instrument's frame, in pixels.
    def instrumentshape(self):
        # get the first instrument element
        instruments = self.tree.xpath("//instruments/*[1]")
        if len(instruments) != 1: raise ValueError("No instruments in ski file")
        # get its shape
        return ( int(instruments[0].get("pixelsX")), int(instruments[0].get("pixelsY")) )

    ## This function returns the angular area (in sr) of a single pixel in the first instrument's frame.
    def angularpixelarea(self):
        # get the first instrument element
        instruments = self.tree.xpath("//instruments/*[1]")
        if len(instruments) != 1: raise ValueError("No instruments in ski file")
        instrument = instruments[0]
        # get the distance in m
        d = self.units().convert(instrument.get("distance"), to_unit='m', quantity='distance')
        # get the field of view in m
        fovx = self.units().convert(instrument.get("fieldOfViewX"), to_unit='m', quantity='length')
        fovy = self.units().convert(instrument.get("fieldOfViewY"), to_unit='m', quantity='length')
        # get the number of pixels
        nx = int(instrument.get("pixelsX"))
        ny = int(instrument.get("pixelsY"))
        # calculate the angular pixel area
        sx = 2 * arctan(fovx / nx / d / 2)
        sy = 2 * arctan(fovy / ny / d / 2)
        return sx * sy

    ## This function returns a list of instrument names, in order of occurrence in the ski file.
    def instrumentnames(self):
        # get the instrument elements
        instruments = self.tree.xpath("//instruments/*")
        # return their names
        return [ instr.get("instrumentName") for instr in instruments ]

    ## This function returns true if the ski file specifies a DustLib with a StaggeredAssigner, false otherwise.
    def staggered(self):
        # get any StaggeredAssigner elements as a DustLib child
        staggered = self.tree.xpath("//dustLib/*/assigner/StaggeredAssigner")
        return len(staggered) > 0

    ## This function returns the dust fraction specified in an SPHDustDistribution,
    # or 0 if the element or the attribute are not present.
    def dustfraction(self):
        # get the value of the relevant attribute on the SPHDustDistribution element (as a list of query results)
        results = self.tree.xpath("//SPHDustDistribution/@dustFraction")
        # if not found, return zero
        if len(results) != 1: return 0
        # convert the first result
        return float(results[0])

    ## This function returns the maximum gas temperature specified in an SPHDustDistribution, in Kelvin,
    # or 0 if the element or the attribute are not present.
    def maximumtemperature(self):
        # get the value of the relevant attribute on the SPHDustDistribution element (as a list of query results)
        results = self.tree.xpath("//SPHDustDistribution/@maximumTemperature")
        # if not found, return zero
        if len(results) != 1: return 0
        # extract the number from the first result, assuming units of K
        return float(results[0].split()[0])

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

    ## This function sets the number of wavelengths
    def setnwavelengths(self, number):
        elems = self.tree.xpath("//wavelengthGrid/*[1]")
        elems[0].set("points", str(number))

    ## This function sets the number of dust cells in the x direction
    def setxdustcells(self, number):
        self.tree.xpath("//dustGridStructure/*[1]")[0].set("pointsX", str(number))

    ## This function sets the number of dust cells in the y direction
    def setydustcells(self, number):
        try:
            self.tree.xpath("//dustGridStructure/*[1]")[0].set("pointsY", str(number))
        except TypeError:
            raise ValueError("The dimension of the dust grid is lower than 2")

    ## This function sets the number of dust cells in the z direction
    def setzdustcells(self, number):
        try:
            self.tree.xpath("//dustGridStructure/*[1]")[0].set("pointsZ", str(number))
        except TypeError:
            raise ValueError("The dimension of the dust grid is lower than 3")

    ## This function increases the number of photon packages by a certain factor
    def increasepackages(self, factor):
        # Set the increased number of packages
        self.setpackages(self.packages()*factor)

    ## This function increases the number of dust cells by a certain factor
    def increasedustcells(self, factor):
        # Get the dimension of the dust grid
        dimension = self.dimension()
        # Set the increased number of dust cells in the x direction
        self.setxdustcells(int(round(self.nxcells() * factor**(1 / float(dimension)))))
        # Set the increased number of dust cells in the y direction
        if dimension > 1: self.setydustcells(int(round(self.nycells() * factor**(1 / float(dimension)))))
        # Set the increased number of dust cells in the z direction
        if dimension > 2: self.setzdustcells(int(round(self.nzcells() * factor**(1 / float(dimension)))))

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
        elems = self.tree.xpath("//SPHStellarComp[./sedFamily/BruzualCharlotSEDFamily]")
        if len(elems) != 1: raise ValueError("No SPHStellarComp with BruzualCharlotSEDFamily in ski file")
        # set the attribute value
        elems[0].set("filename", filename)

    ## This function sets the filename attribute of the SPHStarburstComp element to the specified value.
    def sethiifile(self, filename):
        # get the SPHStarburstComp element
        elems = self.tree.xpath("//SPHStellarComp[./sedFamily/MappingsSEDFamily]")
        if len(elems) != 1: raise ValueError("No SPHStellarComp with MappingsSEDFamily in ski file")
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

    ## This function returns the stellar system
    def get_stellar_system(self):

        return self.get_unique_base_element("stellarSystem")

    ## This function returns the dust system
    def get_dust_system(self):

        return self.get_unique_base_element("dustSystem")

    ## This function returns the list of stellar components
    def get_stellar_components(self, include_comments=False):

        # Get the stellar system
        stellar_system = self.get_stellar_system()

        # Get the 'components' element
        stellar_components_parents = stellar_system.xpath("components")

        # Check if only one 'components' element is present
        if len(stellar_components_parents) == 0: raise ValueError("Stellar system is not composed of components")
        elif len(stellar_components_parents) > 1: raise ValueError("Invalid ski file: multiple 'components' objects within stellar system")
        stellar_components = stellar_components_parents[0]

        # Return the stellar components as a list
        if include_comments: return stellar_components.getchildren()
        else: return [component for component in stellar_components.getchildren() if component.tag is not etree.Comment]

    ## This function returns the dust distribution
    def get_dust_distribution(self):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Return the dust distribution
        return get_unique_element(dust_system, "dustDistribution")

    ## This function returns the list of dust components
    def get_dust_components(self, include_comments=False):

        # Get the dust distribution
        dust_distribution = self.get_dust_distribution()

        # Check whether the dust distribution is a CompDustDistribution
        if not dust_distribution.tag == "CompDustDistribution": raise ValueError("Dust distribution is not composed of components")

        # Get the 'components' element
        dust_components_parents = dust_distribution.xpath("components")

        # Check if only one 'components' element is present
        if len(dust_components_parents) == 0: raise ValueError("Dust distribution is not composed of components")
        elif len(dust_components_parents) > 1: raise ValueError("Invalid ski file: multiple 'components' objects within dust distribution")
        dust_components = dust_components_parents[0]

        # Return the dust components as a list
        if include_comments: return dust_components.getchildren()
        else: return [component for component in dust_components.getchildren() if component.tag is not etree.Comment]

    def get_stellar_component(self, component_id):

        """
        This function returns the stellar component which is preceeded by a comment matching the description
        :param id: index or description
        :return:
        """

        # The component identifier is an integer number -> index of stellar components
        if isinstance(component_id, int):

            # Get all the stellar components (without comments)
            components = self.get_stellar_components()

            # Return the stellar component with the specified index
            return components[component_id]

        # The component identifier is a string -> get stellar component based on description
        elif isinstance(component_id, basestring):

            # Get the stellar components
            components = self.get_stellar_components(include_comments=True)

            # Loop over the different components
            for child in components:

                if child.tag is etree.Comment and child.text.strip() == component_id:

                    # Return the child element right after the comment element
                    return child.getnext()

            # If no match is found, give an error
            raise ValueError("No stellar component found with description '" + component_id + "'")

        # Invalid component id
        else: raise ValueError("Invalid component identifier (should be integer or string)")

    def get_dust_component(self, component_id):

        """
        This function returns the dust component which is preceeded by a comment matching the description
        :param component_id: index or description
        :return:
        """

        # The component identifier is an integer number -> index of dust components
        if isinstance(component_id, int):

            # Get all the dust components (without comments)
            components = self.get_dust_components()

            # Return the dust component with the specified index
            return components[component_id]

        # The component identifier is a string -> get dust component based on description
        elif isinstance(component_id, basestring):

            # Get the dust components
            components = self.get_dust_components(include_comments=True)

            # Loop over the different components
            for child in components:

                if child.tag is etree.Comment and child.text.strip() == component_id:

                    # Return the child element right after the comment element
                    return child.getnext()

            # If no match is found, give an error
            raise ValueError("No dust component found with description '" + component_id + "'")

        # Invalid component id
        else: raise ValueError("Invalid component identifier (should be integer or string)")

    def get_stellar_component_normalization(self, component_id):

        """
        This function returns the normalization element of the specified stellar component
        :param component_id:
        :return:
        """

        # Get the stellar component
        stellar_component = self.get_stellar_component(component_id)

        # Get normalization of this component
        return get_unique_element(stellar_component, "normalization")

    def get_stellar_component_luminosity(self, component_id):

        """
        This function returns the luminosity for a specified stellar component
        :param component_id:
        :return:
        """

        # Get the stellar component normalization of the component
        normalization = self.get_stellar_component_normalization(component_id)

        # Check the type of the normalization
        if normalization.tag == "BolLuminosityStellarCompNormalization":

            # Return the total luminosity and None for the band
            return get_quantity(normalization, "luminosity", default_unit="Lsun"), None

        elif normalization.tag == "LuminosityStellarCompNormalization":

            # Return the luminosity and the corresponding band
            return get_quantity(normalization, "luminosity", default_unit="Lsun"), Filter.from_string(normalization.get("band"))

        elif normalization.tag == "SpectralLuminosityStellarCompNormalization":

            # The (spectral) luminosity
            luminosity = get_quantity(normalization, "luminosity")

            # The wavelength
            wavelength = get_quantity(normalization, "wavelength")

            # Return the luminosity and the wavelength as quantities
            return luminosity, wavelength

    def set_stellar_component_luminosity(self, component_id, luminosity, filter_or_wavelength=None):

        """
        This function sets the luminosity for a specified stellar component
        :param component_id:
        :param luminosity:
        :param filter_or_wavelength:
        :return:
        """

        # Get the stellar component normalization of the component
        normalization = self.get_stellar_component_normalization(component_id)

        # No filter or wavelength is defined, use BolLuminosityStellarCompNormalization
        if filter_or_wavelength is None:

            # Get element that holds the normalization class
            parent = normalization.getparent()

            # Remove the old normalization
            parent.remove(normalization)

            # Make and add the new normalization element
            attrs = {"luminosity" : str(luminosity)}
            parent.append(parent.makeelement("BolLuminosityStellarCompNormalization", attrs))

        # Filter is defined, use LuminosityStellarCompNormalization
        elif isinstance(filter_or_wavelength, Filter):

            # Get element that holds the normalization class
            parent = normalization.getparent()

            # Remove the old normalization
            parent.remove(normalization)

            # Make and add the new normalization element
            attrs = {"luminosity": str(luminosity), "band": filter_or_wavelength.skirt_description}
            parent.append(parent.makeelement("LuminosityStellarCompNormalization", attrs))

        # Wavelength is defined as an Astropy quantity, use SpectralLuminosityStellarCompNormalization
        elif filter_or_wavelength.__class__.__name__ == "Quantity":

            # Get element that holds the normalization class
            parent = normalization.getparent()

            # Remove the old normalization
            parent.remove(normalization)

            # Make and add the new normalization element
            attrs = {"luminosity": str(luminosity), "wavelength": str(filter_or_wavelength)}
            parent.append(parent.makeelement("SpectralLuminosityStellarCompNormalization", attrs))

        # Invalid filter or wavelength argument
        else: raise ValueError("Invalid filter or wavelength")

    def get_dust_component_normalization(self, component_id):

        """
        This function returns the normalization element of the specified dust component
        :param component_id:
        :return:
        """

        # Get the dust component
        dust_component = self.get_dust_component(component_id)

        # Return the normalization
        return get_unique_element(dust_component, "normalization")

    def get_dust_component_mass(self, component_id):

        """
        This function returns the dust mass for a specified dust component
        :param component_id:
        :return:
        """

        # Get the dust component normalization of the component
        normalization = self.get_dust_component_normalization(component_id)

        # Check if the normalization is of type 'DustMassDustCompNormalization'
        if not normalization.tag == "DustMassDustCompNormalization": raise ValueError("Dust component normalization is not of type 'DustMassDustCompNormalization")

        # Get the dust mass and return it as a quantity
        return get_quantity(normalization, "dustMass")

    def set_dust_component_mass(self, component_id, mass):

        """
        This function sets the dust mass to a specific value, for a specified dust component (is none is given, the first dust component)
        :param component_id: None, index (integer) or description (string)
        :param mass: the dust mass as a quantity
        :return:
        """

        # Get the dust component normalization of the component
        normalization = self.get_dust_component_normalization(component_id)

        # Check if the normalization is of type 'DustMassDustCompNormalization'
        if not normalization.tag == "DustMassDustCompNormalization": raise ValueError("Dust component normalization is not of type 'DustMassDustCompNormalization")

        # Set the new dust mass
        normalization.set("dustMass", str(mass.to("Msun")) + " Msun")

    def get_wavelength_grid(self):

        """
        This function returns the wavelength grid
        :return:
        """

        # Get the wavelength grid
        return self.get_unique_base_element("wavelengthGrid")

    def set_file_wavelength_grid(self, filename):

        """
        This function sets the wavelength grid to a file
        :param filename:
        :return:
        """

        # Get the wavelength grid
        wavelength_grid = self.get_wavelength_grid()

        # Get the parent
        parent = wavelength_grid.getparent()

        # Remove the old wavelength grid
        parent.remove(wavelength_grid)

        # Make and add the new wavelength grid
        attrs = {"filename": filename}
        parent.append(parent.makeelement("FileWavelengthGrid", attrs))

    def get_stellar_component_geometry(self, component_id):

        """
        This function ...
        :param component_id:
        :return:
        """

        # Get the stellar component
        stellar_component = self.get_stellar_component(component_id)

        # Return the geometry element of the stellar component
        return get_unique_element(stellar_component, "geometry")

    def set_stellar_component_fits_geometry(self, component_id, filename, pixelscale, position_angle, inclination, x_size, y_size, x_center, y_center, scale_height):

        """
        This function ...
        :param component_id:
        :param filename:
        :param pixelscale:
        :param position_angle:
        :param inclination:
        :param x_size:
        :param y_size:
        :param x_center:
        :param y_center:
        :param scale_height:
        :return:
        """

        # Get the stellar component geometry
        geometry = self.get_stellar_component_geometry(component_id)

        # Get the parent
        parent = geometry.getparent()

        # Remove the old geometry
        parent.remove(geometry)

        # Create and add the new geometry
        attrs = {"filename": filename, "pixelScale": str(pixelscale), "positionAngle": str(position_angle.to("deg")) + " deg",
                 "inclination": str(position_angle.to("deg")) + " deg", "xelements": str(x_size), "yelements": str(y_size),
                 "xcenter": str(x_center), "ycenter": str(y_center)}
        parent.append(parent.makeelement("ReadFitsGeometry", attrs))

    def set_stellar_component_sersic_geometry(self, component_id, index, radius, y_flattening=1, z_flattening=1):

        """
        This function ...
        :param component_id:
        :param index:
        :param radius:
        :param y_flattening:
        :param z_flattening:
        :return:
        """

        # Get the stellar component geometry
        geometry = self.get_stellar_component_geometry(component_id)

        # Get the parent
        parent = geometry.getparent()

        # Remove the old geometry
        parent.remove(geometry)

        # Create and add the new geometry
        attrs = {"yFlattening": str(y_flattening), "zFlattening": str(z_flattening)}
        new_geometry = parent.makeelement("TriaxialGeometryDecorator", attrs)

        # Add sersic profile to the geometry
        attrs = {"index": str(index), "radius": str(radius)}
        new_geometry.append(new_geometry.makeelement("SersicGeometry", attrs))

        # Add the new geometry
        parent.append(new_geometry)

    def set_stellar_component_expdisk_geometry(self, component_id, radial_scale, axial_scale, radial_truncation=0, axial_truncation=0, inner_radius=0):

        """
        This function ...
        :param component_id:
        :param radial_scale:
        :param axial_scale:
        :param radial_truncation:
        :param axial_truncation:
        :param inner_radius:
        :return:
        """

        # Get the stellar component geometry
        geometry = self.get_stellar_component_geometry(component_id)

        # Get the parent
        parent = geometry.getparent()

        # Remove the old geometry
        parent.remove(geometry)

        # Create and add the new exponential disk geometry
        attrs = {"radialScale": str(radial_scale), "axialScale": str(axial_scale), "radialTrunc": str(radial_truncation), "axialTrunc": str(axial_truncation), "innerRadius": str(inner_radius)}
        parent.append(parent.makeelement("ExpDiskGeometry", attrs))

    def get_stellar_component_sed(self, component_id):

        """
        This function ...
        :param component_id:
        :return:
        """

        # Get the stellar component
        component = self.get_stellar_component(component_id)

        # Get the SED element
        return get_unique_element(component, "sed")

    def set_stellar_component_sed(self, component_id, template, age, metallicity):

        """
        This function ...
        :param component_id:
        :param template:
        :param age:
        :param metallicity:
        :return:
        """

        # The name of the template class in SKIRT
        template_class = template + "SED"

        # Get the stellar component SED
        sed = self.get_stellar_component_sed(component_id)

        # Get the parent
        parent = sed.getparent()

        # Remove the old SED element
        parent.remove(sed)

        # Create and add the new geometry
        attrs = {"age": str(age), "metallicity": str(metallicity)}
        parent.append(parent.makeelement(template_class, attrs))

    def get_dust_lib(self):

        """
        This function ...
        :return:
        """

        # Get the dust system
        dust_system = self.get_dust_system()

        # Return the dust lib element
        return get_unique_element(dust_system, "dustLib")

    def set_allcells_dust_lib(self):

        """
        This function ...
        :return:
        """

        # Get the dust lib
        lib = self.get_dust_lib()

        # Get the parent
        parent = lib.getparent()

        # Remove the old DustLib element
        parent.remove(lib)

        # Create and add the new library
        parent.append(parent.makeelement("AllCellsDustLib"), {})

    def set_2d_dust_lib(self, temperature_points=25, wavelength_points=10):

        """
        This function ...
        :param temperature_points:
        :param wavelength_points:
        :return:
        """

        # Get the dust lib
        lib = self.get_dust_lib()

        # Get the parent
        parent = lib.getparent()

        # Remove the old DustLib element
        parent.remove(lib)

        # Create and add the new library
        attrs = {"pointsTemperature": str(temperature_points), "pointsWavelength": str(wavelength_points)}
        parent.append(parent.makeelement("Dim2DustLib"), attrs)

    def set_1d_dust_lib(self, points):

        """
        This function ...
        :param points:
        :return:
        """

        # Get the dust lib
        lib = self.get_dust_lib()

        # Get the parent
        parent = lib.getparent()

        # Remove the old DustLib element
        parent.remove(lib)

        # Create and add the new library
        attrs = {"entries": str(points)}
        parent.append(parent.makeelement("Dim1DustLib"), attrs)

    def to_dict(self):

        """
        This function returns the ski file structure as a nested dictionary
        :return:
        """

        return recursive_dict(self.tree.getroot())

    def to_json(self):

        """
        This function returns the ski file to a json format
        :return:
        """

        import json
        return json.dumps(self.to_dict())

    def get_unique_base_element(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return get_unique_element(self.tree.getroot(), "//"+name)

# -----------------------------------------------------------------

def get_unique_element(element, name):

    """
    This function ...
    :param element:
    :param name:
    :return:
    """

    # Get child element of the given element
    parents = element.xpath(name)

    # Check if only one child element is present
    if len(parents) == 0: raise ValueError("Invalid ski file: no '" + name + "' elements within '" + element.tag + "'")
    elif len(parents) > 1: raise ValueError("Invalid ski file: multiple '" + name + "' elements within '" + element.tag + "'")
    parents = parents[0]

    # Check if only one child object is present
    if len(parents) == 0: raise ValueError("Invalid ski file: no '" + name + "' elements within '" + element.tag + "'")
    elif len(parents) > 1: raise ValueError("Invalid ski file: multiple '" + name + "' elements within '" + element.tag + "'")
    child = parents[0]

    # Return the child element
    return child

# -----------------------------------------------------------------

def get_quantity(element, name, default_unit=None):

    """
    This function ...
    :param element:
    :param name:
    :param default_unit:
    :return:
    """

    from astropy.units import Unit

    splitted = element.get(name).split()
    value = float(splitted[0])
    try: unit = splitted[1]
    except IndexError: unit = default_unit

    # Create a quantity object
    if unit is not None: value = value * Unit(unit)
    return value

# -----------------------------------------------------------------

def recursive_dict(element):
    return element.tag, dict(map(recursive_dict, element)) or element.text

# -----------------------------------------------------------------
