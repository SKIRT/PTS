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
import copy
from collections import OrderedDict
from datetime import datetime
from lxml import etree
from numpy import arctan
import warnings

# Import the relevant PTS classes and modules
from .units import SkirtUnits
from ..filter.filter import parse_filter, Filter
from ..tools import archive as arch
from ..tools import filesystem as fs
from ..tools.stringify import str_from_bool, str_from_angle
from ..tools import xml
from .input import SimulationInput
from ..tools import types

# -----------------------------------------------------------------

# Define parameters in the ski file that are actually quantities but should be entered as scalar values for SKIRT
fake_quantities = dict()
fake_quantities[("BolLuminosityStellarCompNormalization", "luminosity")] = "Lsun"

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
class SkiFile7:
    # ---------- Constructing and saving -----------------------------

    ## The constructor loads the contents of the specified ski file into a new SkiFile instance.
    # The filename \em must end with ".ski" or with "_parameters.xml".
    #
    def __init__(self, filepath=None, tree=None):

        # If filepath is passed
        if filepath is not None:

            # Check
            if tree is not None: raise ValueError("Cannot define both filepath and tree")
            if not filepath.lower().endswith((".ski","_parameters.xml")): raise ValueError("Invalid filename extension for ski file")

            # Set the path to the ski file
            self.path = os.path.expanduser(filepath)

            # Load the XML tree (remove blank text to avoid confusing the pretty printer when saving)
            lines = arch.get_lines(self.path)
            root = etree.fromstringlist(lines, parser=etree.XMLParser(remove_blank_text=True))
            self.tree = etree.ElementTree(root)

            # Replace path by the full, absolute path
            self.path = os.path.abspath(self.path)

        # If tree is passed
        elif tree is not None:

            self.tree = tree
            self.path = None

        # Missing input
        else: raise ValueError("Either filepath or tree must be passed to the constructor")

    ## Open a ski file from a path
    @classmethod
    def from_file(cls, path):
        return cls(filepath=path)

    ## Open a ski file from a remote path (and the remote instance)
    @classmethod
    def from_remote_file(cls, path, remote):

        import StringIO
        #output = StringIO.StringIO()
        #for line in remote.read_lines(path): output.write(line + "\n") # DOESN'T WORK??
        output = StringIO.StringIO(remote.get_text(path)) # WORKS!!

        # Load the XML tree (remove blank text to avoid confusing the pretty printer when saving)
        #tree = etree.fromstring(contents, parser=etree.XMLParser(remove_blank_text=True)) # doesn't work, cannot acces getroot()??
        #tree = etree.parse(remote.read_lines(path, add_sep=True), parser=etree.XMLParser(remove_blank_text=True)) # cannot parse from generator
        tree = etree.parse(output, parser=etree.XMLParser(remove_blank_text=True))

        # Create ski file from the tree
        return cls(tree=tree)

    ## This function converts the tree into a string
    def __str__(self):
        string = etree.tostring(self.tree, encoding="UTF-8", xml_declaration=True, pretty_print=False, with_tail=False)
        lines = string.split(">")
        newstring = ">\n".join(lines)
        return newstring

    ## This function saves the (possibly updated) contents of the SkiFile instance into the specified file.
    # The filename \em must end with ".ski". Saving to and thus replacing the ski file from which this
    # SkiFile instance was originally constructed is allowed, but often not the intention.
    def saveto(self, filepath, update_path=True, fix=False):
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

        # Update the ski file path
        if update_path: self.path = filepath

        # Fix, if requested
        if fix:
            from .skifile import fix_ski_file
            fix_ski_file(filepath)

    ## This function saves the ski file to the original path
    def save(self, fix=False): self.saveto(self.path, fix=fix)

    ## This function returns a copy (a deep copy) of this ski file
    def copy(self):
        ski = copy.deepcopy(self)
        ski.path = None # set the path to None so this copy won't be involuntarily saved over the original file
        return ski

    ## This function returns the ski contents as a string
    def to_string(self):
        return etree.tostring(self.tree)

    ## This function returns the ski contents as a list of strings (lines)
    def to_lines(self):
        return etree.tostringlist(self.tree)

    # ---------- Retrieving information -------------------------------

    @property
    def root(self):
        return self.tree.getroot()

    ## This function returns the simulation element
    def get_simulation(self):
        return self.tree.getroot().getchildren()[0]

    ## This property returns the simulation class (PanMonteCarloSimulation or OligoMonteCarloSimulation)
    @property
    def simulation_class(self):
        return self.get_simulation().tag

    ## This property returns the simulation type
    @property
    def simulation_type(self):
        if self.oligochromatic(): return "oligo"
        elif self.panchromatic(): return "pan"
        else: raise RuntimeError("Something went wrong")

    ## This property returns whether the ski file required input
    @property
    def needs_input(self): return len(self.input_files) > 0

    # This property returns the names of the input files
    @property
    def input_files(self):

        elements = xml.get_all_elements(self.tree.getroot())

        filenames = []

        for el in elements:
            for label in el.attrib:
                value = el.attrib[label]
                if "." in value:
                    before, after = value.rsplit(".", 1)
                    if after[0].isdigit() or after[0] == " ": continue
                    assert after[0].isalpha() # equivalent
                    filenames.append(value)
        return filenames

    # This property returns a list of (element, attribute_name) tuples
    @property
    def input_file_elements_and_attributes(self):

        elements = xml.get_all_elements(self.tree.getroot())

        result = []

        for el in elements:
            for label in el.attrib:
                value = el.attrib[label]
                if "." in value:
                    before, after = value.split(".")
                    if after[0].isdigit() or after[0] == " ": continue
                    assert after[0].isalpha() # equivalent

                    result.append((el, label))

        return result

    ## This function returns the paths to the input files
    def input_paths(self, input_path=None, working_directory=None, ignore_wavelength_grid=False):

        # Get the input file names
        filenames = self.input_files

        # Get the name of the wavelength grid file
        wavelengths_filename = self.wavelengthsfile()

        # Check whether the input path has been specified
        if input_path is None:

            # Initilize the input setting as a list
            input_paths = []

            # Loop over the files
            for filename in filenames:

                # Skip if wavelength grid
                if wavelengths_filename is not None and ignore_wavelength_grid and filename == wavelengths_filename: continue

                # Check whether the file is present in the working directory
                if not fs.contains_file(working_directory, filename): raise IOError("The input file " + filename + " is not present in the working directory")
                else: input_paths.append(fs.join(working_directory, filename))

            # Return the input paths
            return input_paths

        # Construct simulation input
        elif types.is_sequence(input_path): input_path = SimulationInput(*input_path)
        elif types.is_string_type(input_path): input_path = SimulationInput(input_path)
        elif isinstance(input_path, SimulationInput): pass
        else: raise ValueError("Invalid value for 'input_path': " + str(input_path))

        # Initialize a list to contain the file paths
        input_paths = []

        # Loop over the filenames as defined in the ski file
        for filename in filenames:

            # Skip if wavelength grid
            if wavelengths_filename is not None and ignore_wavelength_grid and filename == wavelengths_filename: continue

            # Check if in simulation input
            if filename not in input_path: raise ValueError("The simulation input specification (" + str(input_path) + ") does not contain the file '" + filename + "' needed for ski file '" + self.prefix + "'")

            # Get the path for this file
            path = input_path[filename]

            # Add the input path
            input_paths.append(path)

        # Return the input paths
        return input_paths

    ## This function allows to change the input filename/filepath into another filename/filepath
    def change_input_filename(self, old_filename, new_filename):

        result = self.input_file_elements_and_attributes
        for element, label in result:

            # Check
            filename = self.get_value(element, label)
            if filename != old_filename: continue

            # Set
            self.set_value(element, label, new_filename)

    ## This property gives the simulation prefix
    @property
    def prefix(self):
        if self.path is None: raise RuntimeError("Ski file has no path, so prefix cannot be determined")
        return os.path.basename(self.path).split(".ski")[0]

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

    @property
    def uses_wavelength_file(self):
        grid = self.get_wavelength_grid()
        return grid.tag == "FileWavelengthGrid"

    ## This function returns the number of wavelengths, either defined by the ski file, or from the input wavelengths file
    def get_nwavelengths(self, input_path=None):
        if self.uses_wavelength_file:
            if input_path is None: raise ValueError("Wavelengths are defined in a file but input path was not specified")
            return self.nwavelengthsfile(input_path)
        else: return self.nwavelengths()

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
            raise ValueError("The number of wavelengths is not defined within the ski file. Call nwavelengthsfile().")
        else:
            return int(entry.get("points"))

    ## This function returns the name of the wavelengths file, if any
    def wavelengthsfilename(self):
        # If this ski file contains a file wavelength grid
        entry = self.tree.xpath("//FileWavelengthGrid")
        if entry: return self.get_value(entry[0], "filename")

    ## This function returns the path of the wavelengths file that is used for the simulation, if any
    def wavelengthsfile(self, input_path=None):

        # If this ski file contains a file wavelength grid
        entry = self.tree.xpath("//FileWavelengthGrid")
        if entry:

            filename = self.get_value(entry[0], "filename")

            # Simulation input is specified
            if input_path is not None:

                # TODO: use find_input_filepath (from simulationinput module) here

                # List of file paths
                if types.is_sequence(input_path):

                    for path in input_path:
                        if os.path.basename(path) == filename:
                            wavelengths_path = path
                            break
                    else: raise ValueError("The list of input paths does not contain the path to the wavelengths file")

                # Directory path
                elif types.is_string_type(input_path): wavelengths_path = os.path.join(input_path, filename)

                # Simulation input object
                elif isinstance(input_path, SimulationInput):

                    # Check whether present in simulation input
                    if filename not in input_path: raise ValueError("The file '" + filename + "' with the input wavelengths could not be found within the simulation input specification")

                    # Otherwise, set the path
                    wavelengths_path = input_path[filename]

                # Dictionary
                elif types.is_dictionary(input_path):

                    #input_path = SimulationInput(**input_path)

                    # Check whether present in simulation input
                    if filename not in input_path: raise ValueError("The file '" + filename + "' with the input wavelengths could not be found within the simulation input specification")

                    # Otherwise, set the path
                    wavelengths_path = input_path[filename]

                # Invalid
                else: raise ValueError("Invalid value for 'input_path': '" + str(input_path) + "'")

            # Input path is not specified
            else: wavelengths_path = filename

            # Return the file path
            return wavelengths_path

        # No file wavelength grid in this ski file
        else: return None

    ## This function returns whether the wavelengths file could be found
    def find_wavelengthsfile(self, input_path):
        path = self.wavelengthsfile(input_path)
        return fs.is_file(path)

    ## This function returns the number of wavelength points as defined in the wavelengths file
    def nwavelengthsfile(self, input_path):
        path = self.wavelengthsfile(input_path)
        with open(path, 'r') as f: first_line = f.readline()
        nwavelengths = int(first_line.split("\n")[0])
        return nwavelengths

    ## This function
    def treegridfile(self, input_path=None):
        # If this ski file contains a file tree dust grid
        entry = self.tree.xpath("//FileTreeDustGrid")
        if entry:

            filename = self.get_value(entry[0], "filename")

            # Simulation input is specified
            if input_path is not None:

                # Find the file
                from .input import find_input_filepath
                tree_path = find_input_filepath(filename, input_path)

            # Input path is not specified
            else: tree_path = filename

            # Return the file path
            return tree_path

        # No file tree dust grid in this ski file
        else: return None

    ## This function returns the number of photon packages per wavelength
    def packages(self):
        # Get the MonteCarloSimulation element
        elems = self.tree.xpath("//OligoMonteCarloSimulation | //PanMonteCarloSimulation")
        if len(elems) != 1: raise ValueError("No MonteCarloSimulation in ski file")
        # Get the number of packages
        return int(float(elems[0].get("packages")))

    ## This function looks for elements with a given name
    def find_elements(self, name):
        return self.tree.xpath("//"+name+"/*")

    ## This function looks for a single element with a given name
    def find_element(self, name):
        elements = self.find_elements(name)
        if len(elements) > 1: raise ValueError("Ambigious result")
        elif len(elements) == 0: return None
        else: return elements[0]

    ## This function checks whether an element with the specified name is present in the hierarchy
    def has_element(self, name):
        return self.find_element(name) is not None

    ## This property returns the dimension of the dust grid
    @property
    def griddimension(self):
        if self.has_element("meshX") and self.has_element("meshY") and self.has_element("meshZ"): return 3
        elif self.has_element("meshR") and self.has_element("meshZ"): return 2
        elif self.has_element("meshR"): return 1
        else: raise ValueError("The grid dimension could not be determined")

    ## This function returns the number of dust cells, either defined by the ski file, or from the input grid file
    def get_ncells(self, input_path=None):
        if self.filetreegrid():
            if input_path is None: raise ValueError("Grid is defined in a file but input path was not specified")
            return self.ncellsfiletree(input_path)
        elif self.treegrid(): raise ValueError("The number of dust cells cannot be determined for a regular tree grid")
        else: return self.ncells()

    ## This function returns the number of dust cells
    def ncells(self):

        if self.griddimension == 1: return self.nrcells()
        elif self.griddimension == 2:

            rpoints = self.nrcells()
            zpoints = self.nzcells()

            return rpoints * zpoints

        else:

            xpoints = self.nxcells()
            ypoints = self.nycells()
            zpoints = self.nzcells()

            # Return the total number of dust cells
            return xpoints * ypoints * zpoints

    ## This function returns the number of dust cells in the x direction
    def nxcells(self):
        try:
            xpoints = int(self.tree.xpath("//meshX/*")[0].get("numBins"))
        except (TypeError, IndexError):
            raise ValueError("The number of dust cells is not defined within the ski file")
        return xpoints

    ## This function returns the number of dust cells in the y direction
    def nycells(self):
        try:
            ypoints = int(self.tree.xpath("//meshY/*")[0].get("numBins"))
        except (TypeError, IndexError):
            raise ValueError("The dimension of the dust grid is lower than 2")
        return ypoints

    ## This function returns the number of dust cells in the z direction
    def nzcells(self):
        try:
            zpoints = int(self.tree.xpath("//meshZ/*")[0].get("numBins"))
        except (TypeError, IndexError):
            raise ValueError("The dimension of the dust grid is lower than 3")
        return zpoints

    ## This function returns the number of dust cells in the radial direction
    def nrcells(self):
        return int(self.tree.xpath("//meshR/*")[0].get("numBins"))

    ## This function returns the grid type
    def gridtype(self):
        return self.get_dust_grid().tag

    ## This function returns True if a tree dust grid is used and False otherwise
    def treegrid(self):
        return "Tree" in self.gridtype()

    ## This function returns True if a file tree dust grid is used and False otherwise
    def filetreegrid(self):
        return self.gridtype() == "FileTreeDustGrid"

    ## This property returns the filename of the filetree
    @property
    def filetreegrid_name(self):
        grid = self.get_dust_grid()
        return grid.get("filename")

    ## This function
    def filetreegrid_path(self, input_path):

        filename = self.filetreegrid_name

        # Simulation input is specified
        if input_path is not None:

            # List of file paths
            if types.is_sequence(input_path):

                for path in input_path:
                    if os.path.basename(path) == filename:
                        filepath = path
                        break
                else: raise ValueError("The list of input paths does not contain the path to the file tree grid file")

            # Directory path
            elif types.is_string_type(input_path): filepath = os.path.join(input_path, filename)

            # Simulation input object
            elif isinstance(input_path, SimulationInput):

                # Check whether present in simulation input
                if filename not in input_path: raise ValueError("The file '" + filename + "' with the tree could not be found within the simulation input specification")

                # Otherwise, set the path
                filepath = input_path[filename]

            # Dictionary
            elif types.is_dictionary(input_path):

                # Check whether present in simulation input
                if filename not in input_path: raise ValueError("The file '" + filename + "' with the tree could not be found within the simulation input specification")

                # Otherwise, set the path
                filepath = input_path[filename]

            # Invalid
            else: raise ValueError("Invalid value for 'input_path': '" + str(input_path) + "'")

        # Input path is not specified
        else: filepath = filename

        # Return
        return filepath

    ## This function returns the file tree
    def filetreegrid_tree(self, input_path):
        filepath = self.filetreegrid_path(input_path)
        from .tree import DustGridTree
        return DustGridTree.from_file(filepath)

    ## This function returns the number of dust cells for a filetree grid
    def ncellsfiletree(self, input_path):
        #tree = self.filetreegrid_tree(input_path)
        #return tree.nleaves
        filepath = self.filetreegrid_path(input_path)
        from .tree import get_nleaves
        return get_nleaves(filepath)

    ## This function returns True when a tree dust grid is used that is not a file tree dust grid
    def treegrid_notfile(self):
        return self.treegrid() and not self.filetreegrid()

    ## This function returns the minimum level of the tree
    def tree_min_level(self):
        return int(self.get_dust_grid().attrib["minLevel"])

    ## This function returns the maximum level of the tree
    def tree_max_level(self):
        return int(self.get_dust_grid().attrib["maxLevel"])

    ## This function returns the search method for the tree
    def tree_search_method(self):
        return self.get_dust_grid().attrib["searchMethod"]

    ## This function returns the tree sample count
    def tree_sample_count(self):
        return int(self.get_dust_grid().attrib["sampleCount"])

    ## This function returns the maximum optical depth of the tree
    def tree_max_optical_depth(self):
        return float(self.get_dust_grid().attrib["maxOpticalDepth"])

    ## This function returns the maximum mass fraction of the tree
    def tree_max_mass_fraction(self):
        return float(self.get_dust_grid().attrib["maxMassFraction"])

    # This function returns the maximum density dispersion of the tree
    def tree_max_dens_disp(self):
        return float(self.get_dust_grid().attrib["maxDensDispFraction"])

    ## This function returns the number of dust components
    def ncomponents(self):
        components = self.tree.xpath("//CompDustDistribution/components/*")
        return int(len(components))

    ## This function returns the dust lib type
    def dustlib_type(self):
        return self.get_dust_lib().tag

    ## This function returns the dust lib dimension
    def dustlib_dimension(self):
        if self.dustlib_type == "AllCellsDustLib": return 3
        elif self.dustlib_type == "Dim2DustLib": return 2
        else: return 1

    ## This function returns the number of dust library items
    def nlibitems(self, ncells=None):
        dustlib = self.get_dust_lib()
        if dustlib.tag == "AllCellsDustLib":
            if ncells is not None: return ncells
            else: return self.ncells()
        elif dustlib.tag == "Dim2DustLib":
            temppoints = dustlib.attrib["pointsTemperature"] if "pointsTemperature" in dustlib.attrib else 25
            wavelengthpoints = dustlib.attrib["pointsWavelength"] if "pointsWavelength" in dustlib.attrib else 10
            return temppoints * wavelengthpoints
        elif dustlib.tag == "Dim1DustLib":
            return int(dustlib.attrib["entries"])

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
            elif dustmix.tag == "ThemisDustMix":
                npops += int(dustmix.attrib["hydrocarbonPops"])
                npops += int(dustmix.attrib["silicatePops"])
            elif dustmix.tag == "ZubkoDustMix":
                npops += int(dustmix.attrib["graphitePops"])
                npops += int(dustmix.attrib["silicatePops"])
                npops += int(dustmix.attrib["PAHPops"])*2
        return npops

    ## This function returns the number of simple instruments
    def nsimpleinstruments(self):
        return len(self.tree.xpath("//SimpleInstrument"))

    ## This function returns the number of full instruments
    def nfullinstruments(self):
        return len(self.tree.xpath("//FullInstrument"))

    ## This function returns whether transient heating is enabled
    def transientheating(self):
        return len(self.tree.xpath("//TransientDustEmissivity")) > 0

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

    def enable_selfabsorption(self):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Check if the dust system is of type 'PanDustSystem'
        if dust_system.tag != "PanDustSystem": raise ValueError("Not a panchromatic simulation")

        # Enable dust self-absorption
        self.set_value(dust_system, "selfAbsorption", "true")

    def disable_selfabsorption(self):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Check if the dust system is of type 'PanDustSystem'
        if dust_system.tag != "PanDustSystem": raise ValueError("Not a panchromatic simulation")

        # Disable dust self-absorption
        self.set_value(dust_system, "selfAbsorption", "false")

    def enable_all_dust_system_writing_options(self):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Loop over all elements of the dust system
        for element in dust_system.getiterator():

            # Check if any of the settings of this element is a writing option
            for setting_name, setting_value in element.items():

                # Skip settings that are not writing settings
                if not setting_name.startswith("write"): continue

                # Set the setting to true
                self.set_value(element, setting_name, "true")

    def set_write_convergence(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeConvergence' setting to true
        self.set_value(dust_system, "writeConvergence", str_from_bool(value, lower=True))

    def set_write_density(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeDensity' setting to true
        self.set_value(dust_system, "writeDensity", str_from_bool(value, lower=True))

    def set_write_depth_map(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeDepthMap' setting to true
        self.set_value(dust_system, "writeDepthMap", str_from_bool(value, lower=True))

    def set_write_quality(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeQuality' setting to true
        self.set_value(dust_system, "writeQuality", str_from_bool(value, lower=True))

    def set_write_cell_properties(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeCellProperties' setting to true
        self.set_value(dust_system, "writeCellProperties", str_from_bool(value, lower=True))

    def set_write_stellar_density(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeStellarDensity' setting to true
        self.set_value(dust_system, "writeStellarDensity", str_from_bool(value, lower=True))

    def set_write_cells_crossed(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeCellsCrossed' setting to true
        self.set_value(dust_system, "writeCellsCrossed", str_from_bool(value, lower=True))

    def set_write_emissivity(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeEmissivity' setting to true
        self.set_value(dust_system, "writeEmissivity", str_from_bool(value, lower=True))

    def set_write_temperature(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeTemperature' setting to true
        self.set_value(dust_system, "writeTemperature", str_from_bool(value, lower=True))

    def set_write_isrf(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeISRF' setting to true
        self.set_value(dust_system, "writeISRF", str_from_bool(value, lower=True))

    def set_write_absorption(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeAbsorption' setting to true
        self.set_value(dust_system, "writeAbsorption", str_from_bool(value, lower=True))

    def set_write_spectral_absorption(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeSpectralAbsorption' setting to true
        self.set_value(dust_system, "writeSpectralAbsorption", str_from_bool(value, lower=True))

    def set_write_spectral_emission(self, value=True):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Set the 'writeSpectralEmission' setting to true
        self.set_value(dust_system, "writeSpectralEmission", str_from_bool(value, lower=True))

    def set_write_grid(self, value=True):

        # Get the dust grid
        grid = self.get_dust_grid()

        # Set the 'writeGrid' setting to true
        self.set_value(grid, "writeGrid", str_from_bool(value, lower=True))

    def set_write_grid_tree(self, value=True):

        # Get the dust grid
        if not self.has_tree_dust_grid: raise ValueError("Cannot set this option because no tree dust grid is set")

        # Get the dust grid
        grid = self.get_dust_grid()

        # Set the flag
        self.set_value(grid, "writeTree", str_from_bool(value, lower=True))

    def disable_all_dust_system_writing_options(self):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Loop over all elements of the dust system
        for element in dust_system.getiterator():

            # Check if any of the settings of this element is a writing option
            for setting_name, setting_value in element.items():

                # Skip settings that are not writing settings
                if not setting_name.startswith("write"): continue

                # Set the setting to true
                self.set_value(element, setting_name, "false")

    def enable_all_writing_options(self):

        # Loop over all elements in the tree
        for element in self.tree.getiterator():

            # Check if any of the settings of this element is a writing option
            for setting_name, setting_value in element.items():

                # Skip settings that are not writing settings
                if not setting_name.startswith("write"): continue

                # Set the setting to true
                self.set_value(element, setting_name, "true")

    def disable_all_writing_options(self):

        # Loop over all elements in the tree
        for element in self.tree.getiterator():

            # Check if any of the settings of this element is a writing option
            for setting_name, setting_value in element.items():

                # Skip settings that are not writing settings
                if not setting_name.startswith("write"): continue

                # Set the setting to false
                self.set_value(element, setting_name, "false")

    ## This function returns the number of pixels for each of the instruments
    def npixels(self, nwavelengths=None):
        pixels = []
        nwavelengths = nwavelengths if nwavelengths is not None else self.nwavelengths()
        instruments = self.tree.xpath("//instruments/*")
        for instrument in instruments:
            type = instrument.tag
            name = instrument.attrib["instrumentName"]
            if type == "SimpleInstrument" or type == "FrameInstrument":
                datacube = int(instrument.attrib["pixelsX"]) * int(instrument.attrib["pixelsY"]) * nwavelengths
                pixels.append([name, type, datacube])
            elif type == "FullInstrument":
                datacube = int(instrument.attrib["pixelsX"]) * int(instrument.attrib["pixelsY"]) * nwavelengths
                try: scattlevels = int(instrument.attrib["scatteringLevels"])
                except KeyError: scattlevels = 0
                scattering = scattlevels + 1 if scattlevels > 0 else 0
                dustemission = 1 if self.dustemission() else 0
                npixels = datacube * (3 + scattering + dustemission)
                pixels.append([name, type, npixels])
            elif type == "SEDInstrument":
                pixels.append([name, type, nwavelengths])
        return pixels

    ## This function returns the total number of spatial pixels from all the instruments (an SED counts for one pixel)
    def nspatialpixels(self):
        npixels = 0
        # Loop over the instruments
        for instrument in self.get_instruments():
            # Count SEDInstrument as one pixel, get the number of x and y pixels for other types of instrument and calculate the area
            instrument_type = instrument.tag
            if instrument_type == "SEDInstrument": npixels += 1
            else: npixels += int(instrument.attrib["pixelsX"]) * int(instrument.attrib["pixelsY"])
        # Return the total amount of spatial pixels
        return npixels

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
        return [units.convert(s, to_unit='micron', quantity='wavelength') for s in results[0].split(",")]

    ## This property returns the wavelengths (for an oligochromatic simulation) as quantities
    @property
    def wavelength_list(self):
        from ..units.parsing import parse_unit as u
        return [wavelength * u("micron") for wavelength in self.wavelengths()]

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
        # get its shape (for SKIRT7 and SKIRT8)
        return ( int(instruments[0].get("pixelsX", instruments[0].get("numPixelsX"))),
                 int(instruments[0].get("pixelsY", instruments[0].get("numPixelsY"))) )

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
        # get the number of pixels (for SKIRT7 and SKIRT8)
        nx = int(instrument.get("pixelsX", instrument.get("numPixelsX")))
        ny = int(instrument.get("pixelsY", instrument.get("numPixelsY")))
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

    ## This function returns whether the ski file describes a oligochromatic simulation
    def oligochromatic(self):
        elems = self.tree.xpath("//OligoMonteCarloSimulation")
        return len(elems) > 0

    ## This function returns whether the ski file describes a panchromatic simulation
    def panchromatic(self):
        elems = self.tree.xpath("//PanMonteCarloSimulation")
        return len(elems) > 0

    ## This function converts the ski file to a ski file that describes an oligochromatic simulation
    def to_oligochromatic(self, wavelengths):

        if self.oligochromatic(): warnings.warn("The simulation is already oligochromatic")
        else:

            from ..units.stringify import represent_quantity

            simulation = self.tree.xpath("//PanMonteCarloSimulation")[0]
            simulation.tag = "OligoMonteCarloSimulation"

            # Remove the old wavelength grid
            wavelength_grid = self.get_wavelength_grid()
            parent = wavelength_grid.getparent()
            self.set_value(parent, "type", "OligoWavelengthGrid")
            parent.remove(wavelength_grid)

            # If the wavelength is a quantity, make a list of the one wavelength
            if hasattr(wavelengths, "unit"): wavelengths = [wavelengths]

            # Make the oligochromatic wavelength grid
            attrs = {"wavelengths": ", ".join(map(represent_quantity, wavelengths))}
            parent.append(parent.makeelement("OligoWavelengthGrid", attrs))

            # Adapt the stellar components
            components = self.get_stellar_components()
            for component in components:

                luminosities = ", ".join(["1"] * len(wavelengths))

                component.tag = "OligoStellarComp"
                self.set_value(component, "luminosities", luminosities)

                for child in component.getchildren():
                    if child.tag == "sed" or child.tag == "normalization": component.remove(child)

            # Adapt the dust system
            dust_system = self.get_dust_system()
            parent = dust_system.getparent()
            self.set_value(parent, "type", "OligoDustSystem")
            dust_system.tag = "OligoDustSystem"

            # Remove dust system settings
            if "writeAbsorption" in dust_system.attrib: dust_system.attrib.pop("writeAbsorption")
            dust_system.attrib.pop("writeISRF")
            dust_system.attrib.pop("writeTemperature")
            dust_system.attrib.pop("writeEmissivity")
            dust_system.attrib.pop("selfAbsorption")
            dust_system.attrib.pop("emissionBoost")
            if "cycles" in dust_system.attrib: dust_system.attrib.pop("cycles")
            if "emissionBias" in dust_system.attrib: dust_system.attrib.pop("emissionBias")

            for child in dust_system.getchildren():
                if child.tag == "dustEmissivity" or child.tag == "dustLib": dust_system.remove(child)

    ## This function converts the ski file to a ski file that describes a panchromatic simulation
    def to_panchromatic(self):

        if self.panchromatic(): warnings.warn("The simulation is already panchromatic")
        else:

            simulation = self.tree.xpath("//OligoMonteCarloSimulation")[0]
            simulation.tag = "PanMonteCarloSimulation"

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
        elems[0].set("packages", str(int(number)))

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
        dimension = self.griddimension
        # Set the increased number of dust cells in the x direction
        self.setxdustcells(int(round(self.nxcells() * factor**(1 / float(dimension)))))
        # Set the increased number of dust cells in the y direction
        if dimension > 1: self.setydustcells(int(round(self.nycells() * factor**(1 / float(dimension)))))
        # Set the increased number of dust cells in the z direction
        if dimension > 2: self.setzdustcells(int(round(self.nzcells() * factor**(1 / float(dimension)))))

    ## This function sets the maximum mass fraction of the tree dust grid in the ski file
    # to the specified value
    def setmaxmassfraction(self, number):
        # get the tree dust grid element
        elems = self.tree.xpath("//BinTreeDustGrid | //OctTreeDustGrid")
        if len(elems) != 1: raise ValueError("No tree dust grid in ski file")
        # set the attribute value
        elems[0].set("maxMassFraction", str(number))

    ## This function sets the extent (i.e. half-size in three dimensions) of the tree dust grid in the ski file
    # to the specified value, which must be an astropy quantity expressed in a length unit
    def setdustextent(self, extent):
        # get the tree dust grid element
        elems = self.tree.xpath("//BinTreeDustGrid | //OctTreeDustGrid")
        if len(elems) != 1: raise ValueError("No tree dust grid in ski file")
        # set the attribute value
        self.set_quantity(elems[0], "minX", -extent)
        self.set_quantity(elems[0], "minY", -extent)
        self.set_quantity(elems[0], "minZ", -extent)
        self.set_quantity(elems[0], "maxX", extent)
        self.set_quantity(elems[0], "maxY", extent)
        self.set_quantity(elems[0], "maxZ", extent)

    ## This function sets the dust fraction of the SPH dust distribution in the ski file
    # to the specified value
    def setdustfraction(self, number):
        # get the tree dust grid element
        elems = self.tree.xpath("//SPHDustDistribution")
        if len(elems) != 1: raise ValueError("No SPHDustDistribution in ski file")
        # set the attribute value
        elems[0].set("dustFraction", str(number))

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

    @property
    def has_stellar_system(self):
        try:
            system = self.get_stellar_system()
            return True
        except ValueError: return False

    ## This function returns the dust system
    def get_dust_system(self):
        return self.get_unique_base_element("dustSystem")

    @property
    def has_dust_system(self):
        try:
            system = self.get_dust_system()
            return True
        except ValueError: return False

    ## This function removes the complete dust system
    def remove_dust_system(self):
        dust_system = self.get_dust_system()
        parent = dust_system.getparent()
        parent.getparent().remove(parent)

    ## THis function removes the complete stellar system
    def remove_stellar_system(self):
        stellar_system = self.get_stellar_system()
        parent = stellar_system.getparent()
        parent.getparent().remove(parent)

    ## This property returns the number of stellar components
    @property
    def nstellar_components(self):
        return len(self.get_stellar_component_ids())

    ##
    def get_stellar_components_object(self):

        # Get the stellar system
        stellar_system = self.get_stellar_system()

        # Get the 'components' element
        stellar_components_parents = stellar_system.xpath("components")

        # Check if only one 'components' element is present
        if len(stellar_components_parents) == 0:
            raise ValueError("Stellar system is not composed of components")
        elif len(stellar_components_parents) > 1:
            raise ValueError("Invalid ski file: multiple 'components' objects within stellar system")
        stellar_components = stellar_components_parents[0]

        return stellar_components

    ## This function returns the list of stellar components
    def get_stellar_components(self, include_comments=False):

        stellar_components = self.get_stellar_components_object()

        # Return the stellar components as a list
        if include_comments: return stellar_components.getchildren()
        else: return [component for component in stellar_components.getchildren() if component.tag is not etree.Comment]

    ## This function returns the list of stellar component names
    def get_stellar_component_ids(self):

        # Initialize a list to contain the component ids
        ids = []

        # Get the list of stellar components
        components = self.get_stellar_components(include_comments=True)

        # keep track of the number of actual components
        number_of_components = 0

        # Loop over the components (also includes the comments)
        i = 0
        while i < len(components):

            if components[i].tag is etree.Comment:

                ids.append(components[i].text.strip())
                i += 2 # skip the next component -> it is the component corresponding to this comment

            # No name -> add the index of this component as the ID
            else:
                ids.append(number_of_components)
                i += 1

            # Increment the number of components
            number_of_components += 1

        # Return the list of names
        return ids

    ## This function checks whether a stellar component with a certain ID is present
    def has_stellar_component(self, component_id):
        return component_id in self.get_stellar_component_ids()

    ## This function returns the dust distribution
    def get_dust_distribution(self):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Return the dust distribution
        return xml.get_unique_element(dust_system, "dustDistribution")

    ## This fucntion
    def get_dust_components_object(self):

        # Get the dust distribution
        dust_distribution = self.get_dust_distribution()

        # Check whether the dust distribution is a CompDustDistribution
        if not dust_distribution.tag == "CompDustDistribution": raise ValueError(
            "Dust distribution is not composed of components")

        # Get the 'components' element
        dust_components_parents = dust_distribution.xpath("components")

        # Check if only one 'components' element is present
        if len(dust_components_parents) == 0:
            raise ValueError("Dust distribution is not composed of components")
        elif len(dust_components_parents) > 1:
            raise ValueError("Invalid ski file: multiple 'components' objects within dust distribution")

        dust_components = dust_components_parents[0]
        return dust_components

    ## This function returns the list of dust components
    def get_dust_components(self, include_comments=False):

        dust_components = self.get_dust_components_object()

        # Return the dust components as a list
        if include_comments: return dust_components.getchildren()
        else: return [component for component in dust_components.getchildren() if component.tag is not etree.Comment]

    ## This functions returns a list with the ids of the different dust components (the id is a name if this is defined
    #  for the component, otherwise it is the index of the component)
    def get_dust_component_ids(self):

        # Initialize a list to contain the component ids
        ids = []

        # Get the list of dust components
        components = self.get_dust_components(include_comments=True)

        # keep track of the number of actual components
        number_of_components = 0

        # Loop over the components (also includes the comments)
        i = 0
        while i < len(components):

            if components[i].tag is etree.Comment:

                ids.append(components[i].text.strip())
                i += 2 # skip the next component -> it is the component corresponding to this comment

            # No name -> add the index of this component as the ID
            else:
                ids.append(number_of_components)
                i += 1

            # Increment the number of components
            number_of_components += 1

        # Return the list of names
        return ids

    @property
    def ndust_components(self):
        return len(self.get_dust_component_ids())

    ## This functions checks whether a dust component with a certain ID exists
    def has_dust_component(self, component_id):
        return component_id in self.get_dust_component_ids()

    ## This function returns the stellar component that is recognized by the specified id (index or name)
    def get_stellar_component(self, component_id):

        # The component identifier is an integer number -> index of stellar components
        if types.is_integer_type(component_id):

            # Get all the stellar components (without comments)
            components = self.get_stellar_components()

            # Return the stellar component with the specified index
            return components[component_id]

        # The component identifier is a string -> get stellar component based on description
        elif types.is_string_type(component_id):

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

    ## This function returns the dust component that is recognized by the specified id (index or name)
    def get_dust_component(self, component_id):

        # The component identifier is an integer number -> index of dust components
        if types.is_integer_type(component_id):

            # Get all the dust components (without comments)
            components = self.get_dust_components()

            # Return the dust component with the specified index
            return components[component_id]

        # The component identifier is a string -> get dust component based on description
        elif types.is_string_type(component_id):

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

    ## This functions removes the stellar component with the specified ID
    def remove_stellar_component(self, component_id):

        # Get the stellar component with the specified ID
        component = self.get_stellar_component(component_id)

        # Get the previous item
        previous = component.getprevious()

        # Get the parent
        parent = component.getparent()

        # Check whether the previous item is a comment
        if previous.tag is etree.Comment:

            # If the comment states the component ID, remove it
            if previous.text.strip() == component_id: parent.remove(previous)

            # If the comment preceeding the component does not have the name of that component (it must by definition),
            # something strange is going on ...
            else: raise ValueError("Something is wrong with the ski file")

        # Remove the stellar component
        parent.remove(component)

    ## This function removes the dust component with the specified ID
    def remove_dust_component(self, component_id):

        # Get the dust component with the specified ID
        component = self.get_dust_component(component_id)

        # Get the previous item
        previous = component.getprevious()

        # Get the parent
        parent = component.getparent()

        # Check whether the previous item is a comment
        if previous.tag is etree.Comment:

            # If the comment states the component ID, remove it
            if previous.text.strip() == component_id: parent.remove(previous)

            # If the comment preceeding the component does not have the name of that component (it must by definition),
            # something strange is going on ...
            else: raise ValueError("Something is wrong with the ski file")

        # Remove the dust component
        parent.remove(component)

    ## This function removes the stellar components except for the component(s) with the specified ID(s)
    def remove_stellar_components_except(self, component_ids):

        if types.is_string_type(component_ids): component_ids = [component_ids]

        # Loop over the stellar component IDs
        for id_i in self.get_stellar_component_ids():

            # Skip IDs that are specified by the user
            if id_i in component_ids: continue

            # Remove all other stellar components
            self.remove_stellar_component(id_i)

    ## This function removes the dust components except for the component(s) with the specified ID(s)
    def remove_dust_components_except(self, component_ids):

        if types.is_string_type(component_ids): component_ids = [component_ids]

        # Loop over the dust component IDs
        for id_i in self.get_dust_component_ids():

            # Skip IDs that are specified by the user
            if id_i in component_ids: continue

            # Remove all other dust components
            self.remove_dust_component(id_i)

    ## This function removes all stellar components
    def remove_all_stellar_components(self):

        # Loop over the stellar component IDs
        for id_i in self.get_stellar_component_ids(): self.remove_stellar_component(id_i)

    ## This function removes all dust components
    def remove_all_dust_components(self):

        # Loop over the dust component IDs
        for id_i in self.get_dust_component_ids(): self.remove_dust_component(id_i)

    ## This function creates a new stellar component
    def create_new_stellar_component(self, component_id=None, geometry=None, geometry_type=None, geometry_properties=None,
                                     sed_type=None, sed_properties=None, normalization_type=None,
                                     normalization_properties=None, luminosities=None, sed_template=None, age=None,
                                     filter_or_wavelength=None, luminosity=None, metallicity=None, compactness=None,
                                     pressure=None, covering_factor=None):

        # Panchromatic simulation
        if self.panchromatic():

            # Call the implementation
            self._create_new_stellar_component_pan(component_id=component_id, geometry=geometry, geometry_type=geometry_type,
                                                   geometry_properties=geometry_properties, sed_type=sed_type, sed_properties=sed_properties,
                                                   normalization_type=normalization_type, normalization_properties=normalization_properties, sed_template=sed_template, age=age,
                                                   filter_or_wavelength=filter_or_wavelength, luminosity=luminosity, metallicity=metallicity, compactness=compactness,
                                                   pressure=pressure, covering_factor=covering_factor)

        # Oligochromatic simulation
        else:

            # Check for too much input
            if sed_type is not None or sed_properties is not None: raise ValueError("Cannot specify 'sed_type' or 'sed_properties' for oligochromatic simulations")
            if normalization_type is not None or normalization_properties is not None: raise ValueError("Cannot specify 'normalization_type' or 'normalization_properties' for oligochromatic simulations")

            # Check that luminosities are passed
            if luminosities is None: raise ValueError("Luminosities must be passed")

            # Create
            self._create_new_stellar_component_oligo(component_id=component_id, geometry=geometry, geometry_type=geometry_type,
                                                    geometry_properties=geometry_properties, luminosities=luminosities)

    ## This function is the implementation for adding a new stellar component in panchromatic ski files
    def _create_new_stellar_component_pan(self, component_id=None, geometry=None, geometry_type=None, geometry_properties=None,
                                     sed_type=None, sed_properties=None, normalization_type=None,
                                     normalization_properties=None, sed_template=None, age=None,
                                     filter_or_wavelength=None, luminosity=None, metallicity=None, compactness=None,
                                     pressure=None, covering_factor=None):

        # Get the stellar system
        stellar_system = self.get_stellar_system()

        # Get the 'components' element
        stellar_components_parent = xml.get_unique_element_direct(stellar_system, "components")

        # Create the stellar component
        stellar_component = stellar_components_parent.makeelement("PanStellarComp", {})

        # Create children
        geometry_parent = stellar_component.makeelement("geometry", {"type":"Geometry"})
        sed_parent = stellar_component.makeelement("sed", {"type":"StellarSED"})
        normalization_parent = stellar_component.makeelement("normalization", {"type":"StellarCompNormalization"})

        # Create geometry
        if geometry_type is not None:
            if geometry is not None: raise ValueError("Cannot specify 'geometry' and 'geometry_type'")
            if geometry_properties is None: raise ValueError("Geometry properties must be defined")
            geometry = self.create_element(geometry_type, geometry_properties)
            geometry_parent.append(geometry)

        # Create SED
        if sed_type is not None:
            if sed_properties is None: sed_properties = {}
            sed = self.create_element(sed_type, sed_properties)
            sed_parent.append(sed)

        # Create normalization
        if normalization_type is not None:
            if normalization_properties is None: normalization_properties = {}
            normalization = self.create_element(normalization_type, normalization_properties)
            normalization_parent.append(normalization)

        # Set geometry, sed and normalization to the stellar component
        stellar_component.append(geometry_parent)
        stellar_component.append(sed_parent)
        stellar_component.append(normalization_parent)

        # Add the component ID, if possible
        if component_id is not None:
            comment = etree.Comment(" " + component_id + " ")
            stellar_components_parent.append(comment)

        # Add the new stellar component
        stellar_components_parent.append(stellar_component)

        # Set component ID
        if component_id is None: component_id = self.nstellar_components - 1

        # Set geometry
        if geometry is not None: self.set_stellar_component_geometry(component_id, geometry)

        # Set SED
        if sed_template is not None:

            # Set MAPPINGS or other template
            if sed_template == "Mappings": self.set_stellar_component_mappingssed(component_id, metallicity, compactness, pressure, covering_factor)
            else: self.set_stellar_component_sed(component_id, sed_template, age, metallicity)

        # Set normalization based on filter or wavelength
        if filter_or_wavelength is not None and types.is_quantity(filter_or_wavelength): self.set_stellar_component_normalization_wavelength(component_id, filter_or_wavelength) # actually this step is not needed?
        if luminosity is not None:
            #if filter_or_wavelength is None: raise ValueError("If luminosity is passed, filter or wavelength should be specified") # not true: bolometric is also possible
            self.set_stellar_component_luminosity(component_id, luminosity, filter_or_wavelength=filter_or_wavelength)

    ## This function is the implementation for adding a new stellar component in oligochromatic skifiles
    def _create_new_stellar_component_oligo(self, component_id=None, geometry=None, geometry_type=None,
                                            geometry_properties=None, luminosities=None):

        from ..tools.stringify import tostr

        # Check the number of luminosities
        if not self.uses_wavelength_file and len(luminosities) != self.nwavelengths(): raise ValueError("The number of luminosities must match the number of wavelengths")

        # Get the stellar system
        stellar_system = self.get_stellar_system()

        # Get the 'components' element
        stellar_components_parent = xml.get_unique_element_direct(stellar_system, "components")

        # Set luminositites property for oligochomratic stellar component
        stellar_comp_properties = dict()
        stellar_comp_properties["luminosities"] = ", ".join(map(tostr, luminosities))

        # Create the stellar component
        stellar_component = stellar_components_parent.makeelement("OligoStellarComp", stellar_comp_properties)

        # Create geometry child
        geometry_parent = stellar_component.makeelement("geometry", {"type": "Geometry"})

        # Create geometry
        if geometry_type is not None:
            if geometry is not None: raise ValueError("Cannot specify 'geometry' and 'geometry_type'")
            if geometry_properties is None: raise ValueError("Geometry properties must be defined")
            geometry = self.create_element(geometry_type, geometry_properties)
            geometry_parent.append(geometry)

        # Set geometry to the stellar component
        stellar_component.append(geometry_parent)

        # Add the component ID, if possible
        if component_id is not None:
            comment = etree.Comment(" " + component_id + " ")
            stellar_components_parent.append(comment)

        # Add the new stellar component
        stellar_components_parent.append(stellar_component)

        # Set component ID
        if component_id is None: component_id = self.nstellar_components - 1

        # Set geometry
        if geometry is not None: self.set_stellar_component_geometry(component_id, geometry)

    ## This function creates a new dust component
    def create_new_dust_component(self, component_id=None, geometry=None, geometry_type=None, geometry_properties=None,
                                  mix=None, mix_type=None, mix_properties=None, normalization_type=None,
                                  normalization_value=None, normalization_properties=None, mass=None,
                                  hydrocarbon_pops=25, silicate_pops=25, # for THEMIS
                                  graphite_populations=7, silicate_populations=7, pah_populations=5, write_mix=True, # for Zubko
                                  write_mean_mix=True, write_size=True):

        # Get the dust distribution
        dust_distribution = self.get_dust_distribution()

        # Check dust distribution type
        if dust_distribution.tag != "CompDustDistribution": raise ValueError("The dust distribution is not a 'CompDustDistribution', so adding components is not possible")

        # Get the 'components' element
        dust_components_parent = xml.get_unique_element_direct(dust_distribution, "components")

        # Create the new dut component
        dust_component = dust_components_parent.makeelement("DustComp", {})

        # Create children
        geometry_parent = dust_component.makeelement("geometry", {"type": "Geometry"})
        mix_parent = dust_component.makeelement("mix", {"type": "DustMix"})
        normalization_parent = dust_component.makeelement("normalization", {"type": "DustCompNormalization"})

        # Create geometry
        if geometry_type is not None:
            if geometry is not None: raise ValueError("Cannot specify 'geometry' and 'geometry_type'")
            if geometry_properties is None: raise ValueError("Geometry properties must be defined")
            geometry = self.create_element(geometry_type, geometry_properties)
            geometry_parent.append(geometry)

        # Create mix
        if mix_type is not None:
            if mix is not None: raise ValueError("Cannot specify 'mix' and 'mix_type'")
            if mix_properties is None: raise ValueError("Mix properties must be defined")
            mix = self.create_element(mix_type, mix_properties)
            mix_parent.append(mix)

        # Create normalization
        if normalization_type is not None:
            if normalization_value is not None: raise ValueError("Cannot specify 'normalization_value' and 'normalization_type'")
            if normalization_properties is None: raise ValueError("Normalization properties must be defined")
            normalization = self.create_element(normalization_type, normalization_properties)
            normalization_parent.append(normalization)

        # Set geometry, mix and normalization to the dust component
        dust_component.append(geometry_parent)
        dust_component.append(mix_parent)
        dust_component.append(normalization_parent)

        # Add the component ID, if possible
        if component_id is not None:
            comment = etree.Comment(" " + component_id + " ")
            dust_components_parent.append(comment)

        # Add the new dust component
        dust_components_parent.append(dust_component)

        # Set the component ID
        if component_id is None: component_id = self.ndust_components - 1

        # Set geometry
        if geometry is not None: self.set_dust_component_geometry(component_id, geometry)

        # Set mix
        if mix is not None:
            if mix == "themis": self.set_dust_component_themis_mix(component_id, hydrocarbon_pops=hydrocarbon_pops, silicate_pops=silicate_pops, write_mix=write_mix, write_mean_mix=write_mean_mix, write_size=write_size)
            elif mix == "zubko": self.set_dust_component_zubko_mix(component_id, graphite_populations=graphite_populations, silicate_populations=silicate_populations, pah_populations=pah_populations, write_mix=write_mix, write_mean_mix=write_mean_mix, write_size=write_size)
            else: raise ValueError("Invalid mix: '" + str(mix) + "': must be 'themis' or 'zubko' (for now)")

        # Set normalization
        if normalization_value is not None: self.set_dust_component_normalization(component_id, normalization_value)

        # Set dust mass
        if mass is not None:
            #self.set_dust_component_mass(component_id, mass) # only works for existing normalization element
            self.set_dust_component_normalization(component_id, mass)

    ## This function returns all properties of the stellar component with the specified id
    def get_stellar_component_properties(self, component_id):

        # Get the stellar component
        stellar_component = self.get_stellar_component(component_id)

        # Get the properties
        return xml.get_properties(stellar_component)

    ## This function returns all properties of the stellar component with the specified id
    def get_dust_component_properties(self, component_id):

        # Get the dust component
        dust_component = self.get_dust_component(component_id)

        # Get the properties
        return xml.get_properties(dust_component)

    ## This functions returns the normalization of the stellar component with the specified id
    def get_stellar_component_normalization(self, component_id):

        # Get the stellar component
        stellar_component = self.get_stellar_component(component_id)

        # Get normalization of this component
        return xml.get_unique_element(stellar_component, "normalization")

    ## This function sets the wavelength for the spectral luminosity normalization
    def set_stellar_component_normalization_wavelength(self, component_id, wavelength):

        # Get the normalization element
        normalization = self.get_stellar_component_normalization(component_id)

        # Set the wavelength
        # element, name, value, default_unit=None
        self.set_quantity(normalization, "wavelength", wavelength)

    ## This function sets the spectral luminosity for the spectral luminosity normalization
    def set_stellar_component_normalization_spectral_luminosity(self, component_id, luminosity):

        # Get the normalization element
        normalization = self.get_stellar_component_normalization(component_id)

        # Set the wavelength
        # element, name, value, default_unit=None
        self.set_quantity(normalization, "luminosity", luminosity)

    ## This function returns the wavelength or filter for normalization of the specified stellar component
    def get_stellar_component_normalization_wavelength_or_filter(self, component_id):

        normalization = self.get_stellar_component_normalization(component_id)

        if normalization.tag == "BolLuminosityStellarCompNormalization": return None
        elif normalization.tag == "LuminosityStellarCompNormalization": return parse_filter(normalization.get("band"))
        elif normalization.tag == "SpectralLuminosityStellarCompNormalization": return self.get_quantity(normalization, "wavelength")
        else: raise ValueError("Unrecognized stellar component normalization: " + normalization.tag)

    ## This function returns the luminosity of the stellar component with the specified id,
    #   - if the normalization is by bolometric luminosity, returns (luminosity [as Astropy quantity], None)
    #   - if the normalization is by luminosity in a specific band, returns (luminosity [as Astropy quantity], Filter object)
    #   - if the normalization is by spectral luminosity at a specific wavelength, returns (spectral luminosity [as Astropy quantity], wavelength [as Astropy quantity])
    def get_stellar_component_luminosity(self, component_id, return_wavelength=True):

        # Get the stellar component normalization of the component
        normalization = self.get_stellar_component_normalization(component_id)

        # Check the type of the normalization
        if normalization.tag == "BolLuminosityStellarCompNormalization":

            # Return the total luminosity and None for the band
            return self.get_quantity(normalization, "luminosity", default_unit="Lsun"), None

        elif normalization.tag == "LuminosityStellarCompNormalization":

            # Return the luminosity and the corresponding band
            return self.get_quantity(normalization, "luminosity"), parse_filter(normalization.get("band"))

        elif normalization.tag == "SpectralLuminosityStellarCompNormalization":

            # The (spectral) luminosity
            luminosity = self.get_quantity(normalization, "luminosity")

            # The wavelength
            wavelength = self.get_quantity(normalization, "wavelength")

            # Return the luminosity and the wavelength as quantities
            if return_wavelength: return luminosity, wavelength
            else: return luminosity

    ## For oligochromatic simulations
    def set_stellar_component_luminosities(self, component_id, luminosities):

        # Get the stellar component normalization of the component
        component = self.get_stellar_component(component_id)

        # Set the 'luminosities' attribute
        component.set("luminosities", " ".join(map(str, luminosities)))

    ## This function removes the current stellar component normalization, and returns the parent
    def remove_stellar_component_normalization(self, component_id):

        # Try getting the current normalization
        try:

            # Get the stellar component normalization
            normalization = self.get_stellar_component_normalization(component_id)

            # Get the parent
            parent = normalization.getparent()

            # Remove the old normalization
            parent.remove(normalization)

        # No normalization yet
        except ValueError:

            # Get the stellar component
            stellar_component = self.get_stellar_component(component_id)

            # Get the 'normalization' element
            try: parent = xml.get_unique_element_direct(stellar_component, "normalization")
            except ValueError:
                parent = stellar_component.makeelement("normalization", {"type": "StellarCompNormalization"})
                stellar_component.append(parent)

        # Return the 'normalization' parent
        return parent

    ## This function sets the luminosity of the stellar component with the specified id,
    #  - if filter_or_wavelength is None, the specified luminosity [as Astropy quantity] is interpreted as a bolometric luminosity
    #  - if filter_or_wavelength is a Filter instance, the luminosity [as Astropy quantity] is interpreted as the luminosity in the corresponding band
    #  - if filter_or_wavelength is a wavelength [as an Astropy quantity], the luminosity should be the spectral luminosity [as Astropy quantity] at that wavelength
    def set_stellar_component_luminosity(self, component_id, luminosity, filter_or_wavelength=None):

        from ..units.stringify import represent_quantity

        # Remove stellar component normalization, return the parent
        parent = self.remove_stellar_component_normalization(component_id)

        # No filter or wavelength is defined, use BolLuminosityStellarCompNormalization
        if filter_or_wavelength is None:

            # Make and add the new normalization element
            attrs = {"luminosity" : luminosity.to("Lsun").value}
            parent.append(parent.makeelement("BolLuminosityStellarCompNormalization", attrs))

        # Filter is defined, use LuminosityStellarCompNormalization
        elif isinstance(filter_or_wavelength, Filter):

            # Make and add the new normalization element
            attrs = {"luminosity": represent_quantity(luminosity), "band": filter_or_wavelength.skirt_description}
            parent.append(parent.makeelement("LuminosityStellarCompNormalization", attrs))

        # Wavelength is defined as an Astropy quantity, use SpectralLuminosityStellarCompNormalization
        elif filter_or_wavelength.__class__.__name__ == "Quantity":

            # Make and add the new normalization element
            attrs = {"luminosity": represent_quantity(luminosity), "wavelength": represent_quantity(filter_or_wavelength)}
            parent.append(parent.makeelement("SpectralLuminosityStellarCompNormalization", attrs))

        # Invalid filter or wavelength argument
        else: raise ValueError("Invalid filter or wavelength")

    ## This function returns the normalization of the dust component with the specified id
    def get_dust_component_normalization(self, component_id):

        # Get the dust component
        dust_component = self.get_dust_component(component_id)

        # Return the normalization
        return xml.get_unique_element(dust_component, "normalization")

    ## This function returns the dust mix for the dust component with the specified id
    def get_dust_component_mix(self, component_id):

        # Get the dust component
        dust_component = self.get_dust_component(component_id)

        # Return the dust mix
        return xml.get_unique_element(dust_component, "mix")

    ## This function removes the current dust component mix, and returns the parent
    def remove_dust_component_mix(self, component_id):

        # Try getting the current mix
        try:

            # Get the dust component mix
            mix = self.get_dust_component_mix(component_id)

            # Get the parent
            parent = mix.getparent()

            # Remove the old mix
            parent.remove(mix)

        # No geometry yet
        except ValueError:

            # Get the dust component
            dust_component = self.get_dust_component(component_id)

            # Get the 'mix' element
            try: parent = xml.get_unique_element_direct(dust_component, "mix")
            except ValueError:

                parent = dust_component.makeelement("mix", {"type": "DustMix"})
                dust_component.append(parent)

        # Return the 'mix' parent
        return parent

    ## This functions sets a THEMIS dust mix model for the dust component with the specified id
    def set_dust_component_themis_mix(self, component_id, hydrocarbon_pops=25, silicate_pops=25, write_mix=True, write_mean_mix=True, write_size=True):

        # Remove current mix, return the parent
        parent = self.remove_dust_component_mix(component_id)

        # Make and add the new mix
        attrs = {"writeMix": str_from_bool(write_mix, lower=True), "writeMeanMix": str_from_bool(write_mean_mix, lower=True),
                 "writeSize": str_from_bool(write_size, lower=True), "hydrocarbonPops": str(hydrocarbon_pops),
                 "silicatePops": str(silicate_pops)}
        parent.append(parent.makeelement("ThemisDustMix", attrs))

    ## This function sets a Zubko dust mix model for the dust component with the specified id
    def set_dust_component_zubko_mix(self, component_id, graphite_populations=7, silicate_populations=7, pah_populations=5, write_mix=True, write_mean_mix=True, write_size=True):

        # Remove current mix, return the parent
        parent = self.remove_dust_component_mix(component_id)

        # Make and add the new mix
        attrs = {"writeMix": str_from_bool(write_mix, lower=True), "writeMeanMix": str_from_bool(write_mean_mix, lower=True),
                 "writeSize": str_from_bool(write_size, lower=True), "graphitePops": str(graphite_populations), "silicatePops": str(silicate_populations), "PAHPops": str(pah_populations)}
        parent.append(parent.makeelement("ZubkoDustMix", attrs))

    ## This function returns the mass of the dust component with the specified id, as an Astropy quantity
    def get_dust_component_mass(self, component_id):

        # Get the dust component normalization of the component
        normalization = self.get_dust_component_normalization(component_id)

        # Check if the normalization is of type 'DustMassDustCompNormalization'
        if not normalization.tag == "DustMassDustCompNormalization": raise ValueError("Dust component normalization is not of type 'DustMassDustCompNormalization")

        # Get the dust mass and return it as a quantity
        return self.get_quantity(normalization, "dustMass")

    ## This function removes the current dust component normalization, and returns the parent
    def remove_dust_component_normalization(self, component_id):

        # Try getting the current normalization
        try:

            # Get the dust component normalization
            normalization = self.get_dust_component_normalization(component_id)

            # Get the parent
            parent = normalization.getparent()

            # Remove the old normalization
            parent.remove(normalization)

        # No normalization yet
        except ValueError:

            # Get the dust component
            dust_component = self.get_dust_component(component_id)

            # Get the 'normalization' element
            try: parent = xml.get_unique_element_direct(dust_component, "normalization")
            except ValueError:

                parent = dust_component.makeelement("normalization", {"type": "DustCompNormalization"})
                dust_component.append(parent)

        # Return the 'normalization' parent
        return parent

    ## This function sets the dust component normalization
    def set_dust_component_normalization(self, component_id, value, wavelength=None, direction=None):

        # Get parent
        parent = self.remove_dust_component_normalization(component_id)

        from ..units.helper import is_mass, parse_quantity
        from ..units.stringify import represent_quantity

        # If mass is given: dust mass normalization
        if is_mass(value):

            # Get the mass as a quantity
            mass = parse_quantity(value)

            # Create normalization
            attrs = {"dustMass": represent_quantity(mass)}
            normalization = parent.makeelement("DustMassDustCompNormalization", attrs)

            # Add the normalization
            parent.append(normalization)

        # Other:
        else: raise NotImplementedError("Only dust mass normalization has been implemented")

        #<Type name="RadialDustCompNormalization" concrete="true" base="DustCompNormalization" title="normalization by defining the radial optical depth at some wavelength">
        #<Type name="FaceOnDustCompNormalization" concrete="true" base="DustCompNormalization" title="normalization by defining the face-on optical depth at some wavelength">
        #<Type name="EdgeOnDustCompNormalization" concrete="true" base="DustCompNormalization" title="normalization by defining the edge-on optical depth at some wavelength">
        #<Type name="XDustCompNormalization" concrete="true" base="DustCompNormalization" title="normalization by defining the X-axis optical depth at some wavelength">
        #<Type name="YDustCompNormalization" concrete="true" base="DustCompNormalization" title="normalization by defining the Y-axis optical depth at some wavelength">
        #<Type name="ZDustCompNormalization" concrete="true" base="DustCompNormalization" title="normalization by defining the Z-axis optical depth at some wavelength">
        #<Type name="CompDustDistribution" concrete="true" base="DustDistribution" title="a dust distribution composed of various dust components">

    ## This function sets the mass of the dust component with the specified id. The mass should be an Astropy quantity.
    def set_dust_component_mass(self, component_id, mass):

        # Get the dust component normalization of the component
        normalization = self.get_dust_component_normalization(component_id)

        # Check if the normalization is of type 'DustMassDustCompNormalization'
        if not normalization.tag == "DustMassDustCompNormalization": raise ValueError("Dust component normalization is not of type 'DustMassDustCompNormalization")

        # Set the new dust mass
        self.set_quantity(normalization, "dustMass", mass)

    ## This function returns the wavelength grid
    def get_wavelength_grid(self):
        # Get the wavelength grid
        return self.get_unique_base_element("wavelengthGrid")

    def set_minwavelength(self, value):
        self.set_quantity(self.get_wavelength_grid(), "minWavelength", value)

    def set_maxwavelength(self, value):
        self.set_quantity(self.get_wavelength_grid(), "maxWavelength", value)

    def minwavelength(self):
        return self.get_quantity(self.get_wavelength_grid(), "minWavelength")

    ## This function returns a list of the wavelengths as Quantities, sorted from lowest to highest
    def get_wavelengths(self, input_path=None, as_grid=False):

        from .wavelengthgrid import WavelengthGrid
        from ..units.parsing import parse_unit as u

        # Wavelengths file
        if self.uses_wavelength_file:

            # Check that input path is specified
            if input_path is None: raise ValueError("Input path(s) should be specified")

            # Determine wavelengths file path
            path = self.wavelengthsfile(input_path)

            # Load the wavelength grid
            grid = WavelengthGrid.from_skirt_input(path)

            # Return the wavelengths as a list
            if as_grid: return grid
            else: return sorted(grid.wavelengths(unit="micron"))

        # Wavelengths defined in ski file
        else:

            if as_grid: return WavelengthGrid.from_wavelengths(self.wavelengths(), "micron", sort=True)
            else: return sorted([wavelength * u("micron") for wavelength in self.wavelengths()])

    ## This function returns the minimum wavelength in the grid, but also works when the wavelengths are defined in a file
    def get_min_wavelength(self, input_path):
        wavelengths = self.get_wavelengths(input_path)
        return wavelengths[0]

    ## This function returns the maximum wavelength in the grid, but also works when the wavelengths are defined in a file
    def get_max_wavelength(self, input_path):
        wavelengths = self.get_wavelengths(input_path)
        return wavelengths[-1]

    @property
    def min_wavelength(self):
        return self.minwavelength()

    def maxwavelength(self):
        return self.get_quantity(self.get_wavelength_grid(), "maxWavelength")

    @property
    def max_wavelength(self):
        return self.maxwavelength()

    ## This function sets the number of wavelength points
    def set_nwavelengths(self, value):

        # Get the wavelength grid
        grid = self.get_wavelength_grid()

        # Set the number of points
        grid.set("points", str(value))

    ## This function removes the current wavelength grid and returns the parent
    def remove_wavelength_grid(self):

        try:

            # Get the wavelength grid
            wavelength_grid = self.get_wavelength_grid()

            # Get the parent
            parent = wavelength_grid.getparent()

            # Remove the old wavelength grid
            parent.remove(wavelength_grid)

        except ValueError:

            try: parent = xml.get_unique_element_direct(self.get_simulation(), "wavelengthGrid")
            except ValueError:

                parent = self.get_simulation().makeelement("wavelengthGrid", {"type": "PanWavelengthGrid"})
                self.get_simulation().append(parent)

        return parent

    ## This function sets the wavelengths for an oligochromatic simulation
    def set_wavelengths(self, *args):

        from ..units.stringify import represent_quantity

        # Remove the old wavelength grid
        parent = self.remove_wavelength_grid()

        # Make the oligochromatic wavelength grid
        attrs = {"wavelengths": ", ".join(map(represent_quantity, args))}
        parent.append(parent.makeelement("OligoWavelengthGrid", attrs))

    ## This function sets the wavelength grid to a file
    def set_file_wavelength_grid(self, filename):

        parent = self.remove_wavelength_grid()

        # Make and add the new wavelength grid
        attrs = {"filename": filename}
        parent.append(parent.makeelement("FileWavelengthGrid", attrs))

    ## This function sets the wavelength grid to a NestedLogWavelengthGrid
    def set_nestedlog_wavelength_grid(self, min_lambda, max_lambda, points, min_lambda_sub, max_lambda_sub, points_sub, write):

        from ..units.stringify import represent_quantity

        parent = self.remove_wavelength_grid()

        # Make and add the new wavelength grid
        attrs = {"minWavelength": represent_quantity(min_lambda), "maxWavelength": represent_quantity(max_lambda),
                 "points": str(points), "minWavelengthSubGrid": represent_quantity(min_lambda_sub),
                 "maxWavelengthSubGrid": represent_quantity(max_lambda_sub), "pointsSubGrid": str(points_sub),
                 "writeWavelengths": str_from_bool(write, lower=True)}
        parent.append(parent.makeelement("NestedLogWavelengthGrid", attrs))

    ## This functions sets the wavelength grid to a LogWavelengthGrid
    def set_log_wavelength_grid(self, min_lambda, max_lambda, points, write):

        from ..units.stringify import represent_quantity

        parent = self.remove_wavelength_grid()

        # Make and add the new wavelength grid
        attrs = {"minWavelength": represent_quantity(min_lambda), "maxWavelength": represent_quantity(max_lambda),
                 "points": str(points), "writeWavelengths": str_from_bool(write, lower=True)}
        parent.append(parent.makeelement("LogWavelengthGrid", attrs))

    ## This function returns the geometry of the stellar component with the specified id
    def get_stellar_component_geometry(self, component_id):

        # Get the stellar component
        stellar_component = self.get_stellar_component(component_id)

        # Return the geometry element of the stellar component
        return xml.get_unique_element(stellar_component, "geometry")

    ## This function returns the geometry hierarchy of the stellar component with the specified ID
    def get_stellar_component_geometry_hierarchy_names(self, component_id):

        # Get the stellar component
        geometry = self.get_stellar_component_geometry(component_id)

        # Initialize list for the hierarchy
        hierarchy = []

        # Fill hierarchy
        self._fill_geometry_hierarchy(geometry, hierarchy)

        # Return the hierarchy
        return hierarchy

    ## This function returns the geometry hierarchy of the dust component with the specified ID
    def get_dust_component_geometry_hierarchy_names(self, component_id):

        # Get the dust component geometry
        geometry = self.get_dust_component_geometry(component_id)

        # Initialize list for the hierarchy
        hierarchy = []

        # Fill hiearchy
        self._fill_geometry_hierarchy(geometry, hierarchy)

        # Return the hierarchy
        return hierarchy

    ## This function ...
    def _fill_geometry_hierarchy(self, geometry, hierarchy):

        hierarchy.append(geometry.tag)

        # Decorator
        if geometry.tag.endswith("Decorator"):

            child = xml.get_unique_element(geometry, "geometry")
            return self._fill_geometry_hierarchy(child, hierarchy)

        # Not a decorator
        else: return

    ## This function returns the geometry of the dust component with the specified id
    def get_dust_component_geometry(self, component_id):

        # Get the dust component
        dust_component = self.get_dust_component(component_id)

        # Return the geometry element of the dust component
        return xml.get_unique_element(dust_component, "geometry")

    ## This function rotates the geometry of the specified stellar component
    def rotate_stellar_component(self, component_id, alpha, beta, gamma):

        # alpha: 0 to 360 degrees
        # beta: 0 to 180 degrees
        # gamma: 0 to 360 degrees

        # Get the geomery of the stellar component
        geometry = self.get_stellar_component_geometry(component_id)

        # Get the parent
        parent = geometry.getparent()

        # Remove the old geometry
        parent.remove(geometry)

        # Create the new rotated geometry
        attrs = {"euleralpha": str_from_angle(alpha), "eulerbeta": str_from_angle(beta), "eulergamma": str_from_angle(gamma)}
        new_geometry = parent.makeelement("RotateGeometryDecorator", attrs)

        attrs = {"type": "Geometry"}
        geometry_of_new_geometry = new_geometry.makeelement("geometry", attrs)
        new_geometry.append(geometry_of_new_geometry)

        # Add the original geometry that has to be rotated
        geometry_of_new_geometry.append(geometry)

        # Add the new geometry to the parent
        parent.append(new_geometry)

    ## This function rotates the geometry of the specified dust component
    def rotate_dust_component(self, component_id, alpha, beta, gamma):

        # alpha: 0 to 360 degrees
        # beta: 0 to 180 degrees
        # gamma: 0 to 360 degrees

        # Get the geomery of the dust component
        geometry = self.get_dust_component_geometry(component_id)

        # Get the parent
        parent = geometry.getparent()

        # Remove the old geometry
        parent.remove(geometry)

        # Create the new rotated geometry
        attrs = {"euleralpha": str_from_angle(alpha), "eulerbeta": str_from_angle(beta), "eulergamma": str_from_angle(gamma)}
        new_geometry = parent.makeelement("RotateGeometryDecorator", attrs)

        attrs = {"type": "Geometry"}
        geometry_of_new_geometry = new_geometry.makeelement("geometry", attrs)
        new_geometry.append(geometry_of_new_geometry)

        # Add the original geometry that has to be rotated
        geometry_of_new_geometry.append(geometry)

        # Add the new geometry to the parent
        parent.append(new_geometry)

    ## This function sets the geometry of the specified stellar component to a FITS file
    def set_stellar_component_fits_geometry(self, component_id, filename, pixelscale, position_angle, inclination, x_size, y_size, x_center, y_center, scale_height):

        from ..units.stringify import represent_quantity

        # Get parent
        parent = self.remove_stellar_component_geometry(component_id)

        # Create and add the new geometry
        attrs = {"filename": filename, "pixelScale": represent_quantity(pixelscale), "positionAngle": str_from_angle(position_angle),
                 "inclination": str_from_angle(inclination), "xelements": str(x_size), "yelements": str(y_size),
                 "xcenter": str(x_center), "ycenter": str(y_center), "axialScale": represent_quantity(scale_height)}
        new_geometry = parent.makeelement("ReadFitsGeometry", attrs)
        parent.append(new_geometry)

    ## This function sets the filename of a stellar component ReadFits geometry
    def set_stellar_component_fits_geometry_filename(self, component_id, filename):

        # Get the geometry
        geometry = self.get_stellar_component_geometry(component_id)

        # Set
        self.set_value(geometry, "filename", filename)

    ## This function sets the geometry of the specified dust component to a FITS file
    def set_dust_component_fits_geometry(self, component_id, filename, pixelscale, position_angle, inclination, x_size, y_size, x_center, y_center, scale_height):

        from ..units.stringify import represent_quantity

        # Get parent
        parent = self.remove_dust_component_geometry(component_id)

        # Create and add the new geometry
        attrs = {"filename": filename, "pixelScale": represent_quantity(pixelscale), "positionAngle": str_from_angle(position_angle),
                 "inclination": str_from_angle(inclination), "xelements": str(x_size), "yelements": str(y_size),
                 "xcenter": str(x_center), "ycenter": str(y_center), "axialScale": represent_quantity(scale_height)}
        new_geometry = parent.makeelement("ReadFitsGeometry", attrs)
        parent.append(new_geometry)

    ## This fucntion sets the filename of a dust component ReadFits geometry
    def set_dust_component_fits_geometry_filename(self, component_id, filename):

        # Get the geometry
        geometry = self.get_dust_component_geometry(component_id)

        # Set
        self.set_value(geometry, "filename", filename)

    ## This function sets the geometry of the specified stellar component to a ring geometry
    def set_stellar_component_ring_geometry(self, component_id, radius, width, height):

        from ..units.stringify import represent_quantity

        # Get parent
        parent = self.remove_stellar_component_geometry(component_id)

        # Create and add the new geometry
        attrs = {"radius": represent_quantity(radius), "width": represent_quantity(width), "height": represent_quantity(height)}
        new_geometry = parent.makeelement("RingGeometry", attrs)
        parent.append(new_geometry)

    ## This function sets the geometry of the specified dust component to a ring geometry
    def set_dust_component_ring_geometry(self, component_id, radius, width, height):

        from ..units.stringify import represent_quantity

        # Get parent
        parent = self.remove_dust_component_geometry(component_id)

        # Create and add the new geometry
        attrs = {"radius": represent_quantity(radius), "width": represent_quantity(width),
                 "height": represent_quantity(height)}
        new_geometry = parent.makeelement("RingGeometry", attrs)
        parent.append(new_geometry)

    ## This function returns the geometry of the specified stellar component
    def get_stellar_component_geometry_object(self, component_id):

        from ...modeling.basics.models import SersicModel3D, ExponentialDiskModel3D, DeprojectionModel3D, RingModel3D

        # Get the geometry object
        geometry = self.get_stellar_component_geometry(component_id)

        # RotateGeometryDecorator
        if geometry.tag == "RotateGeometryDecorator":

            # Sersic or exponential
            base_geometry = xml.get_unique_element(geometry, "geometry")

            if base_geometry.tag == "SersicGeometry":

                return SersicModel3D()

            elif base_geometry.tag == "ExpDiskGeometry":

                return ExponentialDiskModel3D()

            elif base_geometry.tag == "SpheroidalGeometryDecorator":

                base_base_geometry = xml.get_unique_element(base_geometry, "geometry")

                if base_base_geometry == "SersicGeometry":

                    return SersicModel3D()

                else: raise NotImplementedError("Rotated version of " + base_base_geometry.tag + " with spheroidal decorator is not supported")

        # Sersic model
        elif geometry.tag == "SersicGeometry":

            return SersicModel3D()

        # Exponential disk
        elif geometry.tag == "ExpDiskGeometry":

            return ExponentialDiskModel3D()

        # Spheroidal
        elif geometry.tag == "SpheroidalGeometryDecorator":

            base_geometry = xml.get_unique_element(geometry, "geometry")

            if base_geometry.tag == "SersicGeometry":

                return SersicModel3D()

            elif base_geometry.tag == "RotateGeometryDecorator":

                base_base_geometry = xml.get_unique_element(base_geometry, "geometry")

                if base_base_geometry == "SersicGeometry":

                    return SersicModel3D()

                else: raise NotImplementedError("Rotated version of " + base_base_geometry.tag + " with spheroidal decorator is not supported")

        # Deprojection
        elif geometry.tag == "ReadFitsGeometry":

            return DeprojectionModel3D()

        # Ring
        elif geometry.tag == "RingGeometry":

            return RingModel3D()

        # Other
        else: raise NotImplementedError("Geometry " + geometry.tag + " not supported")

    ## This function sets the geometry of the specified stellar component.
    def set_stellar_component_geometry(self, component_id, model):

        from astropy.coordinates import Angle
        from ...modeling.basics.models import SersicModel3D, ExponentialDiskModel3D, DeprojectionModel3D, RingModel3D

        # Rotation:
        #  alpha: 0 to 360 degrees
        #  beta: 0 to 180 degrees
        #  gamma: 0 to 360 degrees

        # the first rotation is by an angle α about the Z axis.
        # the second rotation is by an angle β about the new X' axis.
        # the third rotation is by an angle γ about the new Z'' axis.

        # Sersic model
        if isinstance(model, SersicModel3D):

            # Set the Sersic geometry (with flattening)
            self.set_stellar_component_sersic_geometry(component_id, model.index, model.effective_radius, y_flattening=model.y_flattening, z_flattening=model.z_flattening)

            # Determine the Euler angles
            alpha = model.azimuth
            beta = model.tilt
            gamma = Angle(0.0, "deg")

            # Check angles
            if alpha < Angle(0.0, "deg"): # alpha must be between 0 and 360 degrees
                alpha = Angle(360., "deg") + alpha
            if beta < Angle(0.0, "deg"): # beta must be between 0 and 180 degrees, if beta is negative, rotate over z axis with 180 degrees first
                if alpha <= Angle(180, "deg"):# beta must be between 0 and 180 degrees, if beta is negative, rotate over z axis with 180 degrees first
                    alpha += Angle(180, "deg")
                elif alpha > Angle(180, "deg"):
                    alpha = alpha - Angle(180, "deg")
                beta = - beta
            if gamma < Angle(0.0,"deg"): # gamma must be between 0 and 360 degrees
                gamma = Angle(360., "deg") + gamma
            self.rotate_stellar_component(component_id, alpha, beta, gamma)

        # Exponential Disk
        elif isinstance(model, ExponentialDiskModel3D):

            # Set the exponential disk geometry
            radial_scale = model.radial_scale
            axial_scale = model.axial_scale
            radial_truncation = model.radial_truncation
            axial_truncation = model.axial_truncation
            inner_radius = model.inner_radius
            self.set_stellar_component_expdisk_geometry(component_id, radial_scale, axial_scale, radial_truncation, axial_truncation, inner_radius)

            # Rotate the exponential disk geometry with the tilt angle
            alpha = Angle(0.0, "deg")
            beta = model.tilt
            #print("beta", beta)
            gamma = Angle(0.0, "deg")
            if beta < Angle(0.0, "deg"): # beta must be between 0 and 180 degrees, if beta is negative, rotate over z axis with 180 degrees first
                alpha = Angle(180, "deg")
                beta = - beta
            #print(alpha, beta, gamma)
            self.rotate_stellar_component(component_id, alpha, beta, gamma)

        # Deprojection model
        elif isinstance(model, DeprojectionModel3D):

            # Set the ReadFitsGeometry
            filename = model.filename
            scale = model.pixelscale
            pa = model.position_angle
            i = model.inclination
            nx = model.x_size
            ny = model.y_size
            xc = model.x_center
            yc = model.y_center
            hz = model.scale_height
            self.set_stellar_component_fits_geometry(component_id, filename, scale, pa, i, nx, ny, xc, yc, hz)

        # Ring model
        elif isinstance(model, RingModel3D):

            # Get the properties
            radius = model.radius
            width = model.width
            height = model.height

            # Set the geometry
            self.set_stellar_component_ring_geometry(component_id, radius, width, height)

        # Unsupported model
        else: raise ValueError("Models other than SersicModel3D, ExponentialDiskModel3D, RingModel3D, and DeprojectionModel3D are not supported yet. This model is of type " + str(type(model)))

    ## This function returns the geometry of the specified dust component
    def get_dust_component_geometry_object(self, component_id):

        from ...modeling.basics.models import SersicModel3D, ExponentialDiskModel3D, DeprojectionModel3D, RingModel3D
        pass

    ## This function sets the geometry of the specified dust component
    def set_dust_component_geometry(self, component_id, model):

        from astropy.coordinates import Angle
        from ...modeling.basics.models import SersicModel3D, ExponentialDiskModel3D, DeprojectionModel3D, RingModel3D

        # Rotation:
        #  alpha: 0 to 360 degrees
        #  beta: 0 to 180 degrees
        #  gamma: 0 to 360 degrees

        # Sersic model
        if isinstance(model, SersicModel3D):

            # Set the Sersic geometry (with flattening)
            self.set_dust_component_sersic_geometry(component_id, model.index, model.effective_radius, z_flattening=model.flattening)

            # Rotate the Sersic geometry with the tilt angle
            alpha = Angle(0.0, "deg")
            beta = model.tilt
            gamma = Angle(0.0, "deg")
            if beta < Angle(0.0, "deg"): # beta must be between 0 and 180 degrees, if beta is negative, rotate over z axis with 180 degrees first
                alpha = Angle(180, "deg")
                beta = - beta
            self.rotate_dust_component(component_id, alpha, beta, gamma)

        # Exponential Disk
        elif isinstance(model, ExponentialDiskModel3D):

            # Set the exponential disk geometry
            radial_scale = model.radial_scale
            axial_scale = model.axial_scale
            radial_truncation = model.radial_truncation
            axial_truncation = model.axial_truncation
            inner_radius = model.inner_radius
            self.set_dust_component_expdisk_geometry(component_id, radial_scale, axial_scale, radial_truncation, axial_truncation, inner_radius)

            # Rotate the exponential disk geometry with the tilt angle
            alpha = Angle(0.0, "deg")
            beta = model.tilt
            #print("beta", beta)
            gamma = Angle(0.0, "deg")
            if beta < Angle(0.0, "deg"): # beta must be between 0 and 180 degrees, if beta is negative, rotate over z axis with 180 degrees first
                alpha = Angle(180, "deg")
                beta = - beta
            #print(alpha, beta, gamma)
            self.rotate_dust_component(component_id, alpha, beta, gamma)

        # Deprojection model
        elif isinstance(model, DeprojectionModel3D):

            # Set the ReadFitsGeometry
            filename = model.filename
            scale = model.pixelscale
            pa = model.position_angle
            i = model.inclination
            nx = model.x_size
            ny = model.y_size
            xc = model.x_center
            yc = model.y_center
            hz = model.scale_height
            self.set_dust_component_fits_geometry(component_id, filename, scale, pa, i, nx, ny, xc, yc, hz)

        # Ring model
        elif isinstance(model, RingModel3D):

            # Get the properties
            radius = model.radius
            width = model.width
            height = model.height

            # Set the geometry
            self.set_dust_component_ring_geometry(component_id, radius, width, height)

        # Unsupported model
        else: raise ValueError("Models other than SersicModel3D, ExponentialDiskModel3D, RingModel3D, and DeprojectionModel3D are not supported yet")

    ## This function adds clumpiness to a stellar component geometry
    def add_stellar_component_clumpiness(self, component_id, fraction, count, radius, cutoff=False, kernel_type="uniform"):

        from ..units.stringify import represent_quantity

        # Get the geometry
        geometry = self.get_stellar_component_geometry(component_id)

        # Remove the old geometry from the tree
        parent = geometry.getparent()
        parent.remove(geometry)

        # Create decorator
        class_name = "ClumpyGeometryDecorator"

        # Set attributes dictionary
        attrs = dict()
        attrs["clumpFraction"] = repr(fraction) # min: 0, max: 1
        attrs["clumpCount"] = str(count) # min: 1
        attrs["clumpRadius"] = represent_quantity(radius) # min: 0 m
        attrs["cutoff"] = "true" if cutoff else "false" # default: false
        decorator = self.tree.getroot().makeelement(class_name, attrs)

        # Add smoothing kernel
        kernel_parent = decorator.makeelement("kernel", {"type": "SmoothingKernel"})
        if kernel_type == "uniform": kernel = kernel_parent.makeelement("UniformSmoothingKernel")
        elif kernel_type == "cubic_spline": kernel = kernel_parent.makeelement("CubicSplineSmoothingKernel")
        else: raise ValueError("Invalid kernel type")
        kernel_parent.append(kernel)
        decorator.append(kernel_parent)

        # Add the underlying geometry
        geometry_parent = decorator.makeelement("geometry", {"type": "Geometry"})
        geometry_parent.append(geometry)
        decorator.append(geometry_parent)

        # Set new stellar component geometry
        parent.append(decorator)

    ## This function adds clumpiness to a dust component geometry
    def add_dust_component_clumpiness(self, component_id, fraction, count, radius, cutoff=False, kernel_type="uniform"):

        from ..units.stringify import represent_quantity

        # Get the geometry
        geometry = self.get_dust_component_geometry(component_id)

        # Remove the old geometry from the tree
        parent = geometry.getparent()
        parent.remove(geometry)

        # Create decorator
        class_name = "ClumpyGeometryDecorator"

        # Set attributes dictionary
        attrs = dict()
        attrs["clumpFraction"] = repr(fraction)  # min: 0, max: 1
        attrs["clumpCount"] = str(count)  # min: 1
        attrs["clumpRadius"] = represent_quantity(radius)  # min: 0 m
        attrs["cutoff"] = "true" if cutoff else "false"  # default: false
        decorator = self.tree.getroot().makeelement(class_name, attrs)

        # Add smoothing kernel
        kernel_parent = decorator.makeelement("kernel", {"type": "SmoothingKernel"})
        if kernel_type == "uniform": kernel = kernel_parent.makeelement("UniformSmoothingKernel")
        elif kernel_type == "cubic_spline": kernel = kernel_parent.makeelement("CubicSplineSmoothingKernel")
        else: raise ValueError("Invalid kernel type")
        kernel_parent.append(kernel)
        decorator.append(kernel_parent)

        # Add the underlying geometry
        geometry_parent = decorator.makeelement("geometry", {"type": "Geometry"})
        geometry_parent.append(geometry)
        decorator.append(geometry_parent)

        # Set new dust component geometry
        parent.append(decorator)

    ## This function adds spiral structure to a stellar component
    def add_stellar_component_spiral_structure(self, component_id, radius, perturbation_weight, arms=1, pitch=0.1745329252, phase=0, index=1):

        class_name = "SpiralStructureGeometryDecorator"

        from ..units.stringify import represent_quantity

        # Get the geometry
        geometry = self.get_dust_component_geometry(component_id)

        # Remove the old geometry from the tree
        parent = geometry.getparent()
        parent.remove(geometry)

        # Set attributes dictionary
        attrs = dict()

        # Set the properties
        attrs["arms"] = str(arms)
        attrs["pitch"] = str_from_angle(pitch) # min: 0 rad, 1.570796327 rad, default: 0.1745329252 rad
        attrs["radius"] = represent_quantity(radius) # min: 0 m
        attrs["phase"] = str_from_angle(phase) # min: 0 rad, max: 6.283185307 rad, default: 0 rad
        attrs["perturbWeight"] = repr(perturbation_weight) # min: 0, max: 1
        attrs["index"] = str(index) # min: 0, max: 10, default: 1
        decorator = self.tree.getroot().makeelement(class_name, attrs)

        # Add the underlying geometry
        geometry_parent = decorator.makeelement("geometry", {"type": "AxGeometry"})
        geometry_parent.append(geometry)
        decorator.append(geometry_parent)

        # Set the new stellar component geometry
        parent.append(decorator)

    ## THis function adds spiral structure to a dust component
    def add_dust_component_spiral_structure(self, component_id, radius, perturbation_weight, arms=1, pitch=0.1745329252, phase=0, index=1):

        class_name = "SpiralStructureGeometryDecorator"

        from ..units.stringify import represent_quantity

        # Get the geometry
        geometry = self.get_dust_component_geometry(component_id)

        # Remove the old geometry from the tree
        parent = geometry.getparent()
        parent.remove(geometry)

        # Set attributes dictionary
        attrs = dict()

        # Set the properties
        attrs["arms"] = str(arms)
        attrs["pitch"] = str_from_angle(pitch)  # min: 0 rad, 1.570796327 rad, default: 0.1745329252 rad
        attrs["radius"] = represent_quantity(radius)  # min: 0 m
        attrs["phase"] = str_from_angle(phase)  # min: 0 rad, max: 6.283185307 rad, default: 0 rad
        attrs["perturbWeight"] = repr(perturbation_weight)  # min: 0, max: 1
        attrs["index"] = str(index)  # min: 0, max: 10, default: 1
        decorator = self.tree.getroot().makeelement(class_name, attrs)

        # Add the underlying geometry
        geometry_parent = decorator.makeelement("geometry", {"type": "AxGeometry"})
        geometry_parent.append(geometry)
        decorator.append(geometry_parent)

        # Set the new dust component geometry
        parent.append(decorator)

    ## This function removes the stellar component geometry, and returns the 'geometry' parent
    def remove_stellar_component_geometry(self, component_id):

        # Try getting the current geometry
        try:

            # Get the stellar component geometry
            geometry = self.get_stellar_component_geometry(component_id)

            # Get the parent
            parent = geometry.getparent()

            # Remove the old geometry
            parent.remove(geometry)

        # No geometry yet
        except ValueError:

            # Get the stellar component
            stellar_component = self.get_stellar_component(component_id)

            # Get the 'geometry' element
            try: parent = xml.get_unique_element_direct(stellar_component, "geometry")
            except ValueError:

                parent = stellar_component.makeelement("geometry", {"type": "Geometry"})
                stellar_component.append(parent)

        # Return the 'geometry' parent
        return parent

    ## This function removes the dust component geometry, and returns the '
    def remove_dust_component_geometry(self, component_id):

        # Try getting the current geometry
        try:

            # Get the dust component geometry
            geometry = self.get_dust_component_geometry(component_id)

            # Get the parent
            parent = geometry.getparent()

            # Remove the old geometry
            parent.remove(geometry)

        # No geometry yet
        except ValueError:

            # Get the dust component
            dust_component = self.get_dust_component(component_id)

            # Get the 'geometry' element
            try: parent = xml.get_unique_element_direct(dust_component, "geometry")
            except ValueError:

                parent = dust_component.makeelement("geometry", {"type": "Geometry"})
                dust_component.append(parent)

        # Return the 'geometry' parent
        return parent

    ## This function sets the geometry of the specified stellar component to a Sersic profile with an specific y and z flattening
    def set_stellar_component_sersic_geometry(self, component_id, index, radius, y_flattening=1, z_flattening=1):

        # Get the parent, remove current geometry
        parent = self.remove_stellar_component_geometry(component_id)

        # Create and add the new geometry
        attrs = {"yFlattening": str(y_flattening), "zFlattening": str(z_flattening)}
        new_geometry = parent.makeelement("TriaxialGeometryDecorator", attrs)

        attrs = {"type": "SpheGeometry"}
        geometry_of_new_geometry = new_geometry.makeelement("geometry", attrs)
        new_geometry.append(geometry_of_new_geometry)

        # Add sersic profile to the geometry
        attrs = {"index": str(index), "radius": str(radius)}
        sersic_geometry = geometry_of_new_geometry.makeelement("SersicGeometry", attrs)
        geometry_of_new_geometry.append(sersic_geometry)

        # Add the new geometry
        parent.append(new_geometry)

    ## This function sets the geometry of the specified dust component to a Sersic profile with a specific y and z flattening
    def set_dust_component_sersic_geometry(self, component_id, index, radius, y_flattening=1, z_flattening=1):

        # Get the parent, remove current geometry
        parent = self.remove_dust_component_geometry(component_id)

        # Create and add the new geometry
        attrs = {"yFlattening": str(y_flattening), "zFlattening": str(z_flattening)}
        new_geometry = parent.makeelement("TriaxialGeometryDecorator", attrs)

        attrs = {"type": "SpheGeometry"}
        geometry_of_new_geometry = new_geometry.makeelement("geometry", attrs)
        new_geometry.append(geometry_of_new_geometry)

        # Add sersic profile to the geometry
        attrs = {"index": str(index), "radius": str(radius)}
        sersic_geometry = geometry_of_new_geometry.makeelement("SersicGeometry", attrs)
        geometry_of_new_geometry.append(sersic_geometry)

        # Add the new geometry
        parent.append(new_geometry)

    ## This function sets the geometry of the specified stellar component to an exponential disk profile
    def set_stellar_component_expdisk_geometry(self, component_id, radial_scale, axial_scale, radial_truncation=0, axial_truncation=0, inner_radius=0):

        # Get the parent, remove current geometry
        parent = self.remove_stellar_component_geometry(component_id)

        # Create and add the new exponential disk geometry
        attrs = {"radialScale": str(radial_scale), "axialScale": str(axial_scale), "radialTrunc": str(radial_truncation), "axialTrunc": str(axial_truncation), "innerRadius": str(inner_radius)}
        new_geometry = parent.makeelement("ExpDiskGeometry", attrs)

        # Add the new geometry
        parent.append(new_geometry)

    ## This function sets the geometry of the specified dust component to an exponential disk profile
    def set_dust_component_expdisk_geometry(self, component_id, radial_scale, axial_scale, radial_truncation=0, axial_truncation=0, inner_radius=0):

        # Get the parent, remove current geometry
        parent = self.remove_dust_component_geometry(component_id)

        # Create and add the new exponential disk geometry
        attrs = {"radialScale": str(radial_scale), "axialScale": str(axial_scale), "radialTrunc": str(radial_truncation), "axialTrunc": str(axial_truncation), "innerRadius": str(inner_radius)}
        new_geometry = parent.makeelement("ExpDiskGeometry", attrs)

        # Add the new geometry
        parent.append(new_geometry)

    ## This function returns the SED template of the specified stellar component
    def get_stellar_component_sed(self, component_id):

        # Get the stellar component
        component = self.get_stellar_component(component_id)

        # Get the SED element
        return xml.get_unique_element(component, "sed")

    ## This function removes the current stellar component SED template, and returns the parent
    def remove_stellar_component_sed(self, component_id):

        # Try getting the current SED
        try:

            # Get the stellar component SED
            sed = self.get_stellar_component_sed(component_id)

            # Get the parent
            parent = sed.getparent()

            # Remove the old sED
            parent.remove(sed)

        # No geometry yet
        except ValueError:

            # Get the stellar component
            stellar_component = self.get_stellar_component(component_id)

            # Get the 'sed' element
            try: parent = xml.get_unique_element_direct(stellar_component, "sed")
            except ValueError:

                parent = stellar_component.makeelement("sed", {"type": "StellarSED"})
                stellar_component.append(parent)

        # Return the 'sed' parent
        return parent

    ## This function sets the SED template of the specified stellar component to a certain model with a specific age
    #  and metallicity (but not MAPPINGS SED)
    def set_stellar_component_sed(self, component_id, template, age, metallicity):

        # Get the parent, remove current SED
        parent = self.remove_stellar_component_sed(component_id)

        # The name of the template class in SKIRT
        template_class = template + "SED"

        # Create and add the new geometry
        attrs = {"age": str(age), "metallicity": str(metallicity)}
        parent.append(parent.makeelement(template_class, attrs))

    ## This function sets a MAPPINGS SED template for the stellar component with the specified id
    def set_stellar_component_mappingssed(self, component_id, metallicity, compactness, pressure, covering_factor):

        from ..units.stringify import represent_quantity

        # Get the parent, remove current SED
        parent = self.remove_stellar_component_sed(component_id)

        # Create and add the new geometry
        attrs = {"metallicity": str(metallicity), "compactness": str(compactness), "pressure": represent_quantity(pressure), "coveringFactor": str(covering_factor)}
        parent.append(parent.makeelement("MappingsSED", attrs))

    ## This function returns the dust emissivity
    def get_dust_emissivity(self):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Return the dust emissivity element
        return xml.get_unique_element(dust_system, "dustEmissivity")

    ## This property returns whether a dust emissivity object is present in the ski file
    @property
    def has_dust_emissivity(self):
        try:
            em = self.get_dust_emissivity()
            return True
        except ValueError: return False

    @property
    def transient_dust_emissivity(self):

        if self.has_dust_emissivity: return self.get_dust_emissivity().tag == "TransientDustEmissivity"
        else: return None

    @property
    def grey_body_dust_emissivity(self):

        if self.has_dust_emissivity: return self.get_dust_emissivity().tag == "GreyBodyDustEmissivity"
        else: return None

    ## This function removes the current dust emissivity and returns the parent
    def remove_dust_emissivity(self):

        try:
            # Get the dust emissivity
            emissivity = self.get_dust_emissivity()

            # Get the parent
            parent = emissivity.getparent()

            # Remove the old emissivity
            parent.remove(emissivity)

        except ValueError:

            dust_system = self.get_dust_system()
            parent = dust_system.makeelement("dustEmissivity", {"type": "DustEmissivity"})
            dust_system.append(parent)

        return parent

    ## This function sets a transient dust emissivity for the simulation
    def set_transient_dust_emissivity(self):

        parent = self.remove_dust_emissivity()

        # Create and add the new emissivity
        parent.append(parent.makeelement("TransientDustEmissivity", {}))

    ## This function sets a grey body dust emissivity for the simulation
    def set_grey_body_dust_emissivity(self):

        parent = self.remove_dust_emissivity()

        # Create and add the new emissivity
        parent.append(parent.makeelement("GreyBodyDustEmissivity", {}))

    ## This function returns the dust library
    def get_dust_lib(self):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Return the dust lib element
        return xml.get_unique_element(dust_system, "dustLib")

    ## This function removes the current dust library and returns the parent
    def remove_dust_lib(self):

        try:
            # Get the dust lib
            lib = self.get_dust_lib()

            # Get the parent
            parent = lib.getparent()

            # Remove the old DustLib element
            parent.remove(lib)

        except ValueError:

            dust_system = self.get_dust_system()
            parent = dust_system.makeelement("dustLib", {"type": "DustLib"})
            dust_system.append(parent)

        return parent

    ## This function sets the dust library to an AllCellsDustLib
    def set_allcells_dust_lib(self):

        parent = self.remove_dust_lib()

        # Create and add the new library
        parent.append(parent.makeelement("AllCellsDustLib", {}))

    ## This function sets the dust library to a 2D dust library
    def set_2d_dust_lib(self, temperature_points=25, wavelength_points=10):

        parent = self.remove_dust_lib()

        # Create and add the new library
        attrs = {"pointsTemperature": str(temperature_points), "pointsWavelength": str(wavelength_points)}
        parent.append(parent.makeelement("Dim2DustLib", attrs))

    ## This function sets the dust library to a 1D dust library
    def set_1d_dust_lib(self, points):

        parent = self.remove_dust_lib()

        # Create and add the new library
        attrs = {"entries": str(points)}
        parent.append(parent.makeelement("Dim1DustLib", attrs))

    ## This function returns the dust grid
    def get_dust_grid(self):

        # Get the dust system
        dust_system = self.get_dust_system()

        # Return the dust grid
        return xml.get_unique_element(dust_system, "dustGrid")

    @property
    def has_tree_dust_grid(self):
        # Get the dust grid
        grid = self.get_dust_grid()
        return grid.tag == "BinTreeDustGrid" or grid.tag == "OctTreeDustGrid"

    ## This function returns the dust grid as a DustGrid object
    def get_dust_grid_object(self):

        from .grids import BinaryTreeDustGrid, OctTreeDustGrid, CartesianDustGrid, FileTreeDustGrid, CylindricalGrid

        # Get the dust grid
        grid = self.get_dust_grid()

        # Cartesian grid
        if grid.tag == "CartesianDustGrid":

            write = self.get_boolean(grid, "writeGrid")

            min_x = self.get_quantity(grid, "minX")
            max_x = self.get_quantity(grid, "maxX")
            min_y = self.get_quantity(grid, "minY")
            max_y = self.get_quantity(grid, "maxY")
            min_z = self.get_quantity(grid, "minZ")
            max_z = self.get_quantity(grid, "maxZ")

            mesh_x = xml.get_unique_element(grid, "meshX")
            xbins = int(mesh_x.get("numBins"))
            xratio = int(mesh_x.get("ratio"))

            mesh_y = xml.get_unique_element(grid, "meshY")
            ybins = int(mesh_y.get("numBins"))
            yratio = int(mesh_y.get("ratio"))

            mesh_z = xml.get_unique_element(grid, "meshZ")
            zbins = int(mesh_z.get("numBins"))
            zratio = int(mesh_z.get("ratio"))

            # Create and return the grid
            return CartesianDustGrid(x_bins=xbins, y_bins=ybins, z_bins=zbins, mesh_type="symmetric_power", ratio=xratio,
                                     min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y, min_z=min_z, max_z=max_z, write=write)

        # Binary tree dust grid
        elif grid.tag == "BinTreeDustGrid":

            write = self.get_boolean(grid, "writeGrid")

            min_x = self.get_quantity(grid, "minX")
            max_x = self.get_quantity(grid, "maxX")
            min_y = self.get_quantity(grid, "minY")
            max_y = self.get_quantity(grid, "maxY")
            min_z = self.get_quantity(grid, "minZ")
            max_z = self.get_quantity(grid, "maxZ")

            min_level = int(grid.get("minLevel"))
            max_level = int(grid.get("maxLevel"))
            search_method =  grid.get("searchMethod")
            sample_count = int(grid.get("sampleCount"))
            maxoptdepth = float(grid.get("maxOpticalDepth"))
            maxmassfraction = float(grid.get("maxMassFraction"))
            maxdensdispfraction = float(grid.get("maxDensDispFraction"))
            directionmethod = grid.get("directionMethod")

            # Create and return the grid
            return BinaryTreeDustGrid(min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y, min_z=min_z, max_z=max_z,
                                      min_level=min_level, max_level=max_level, search_method=search_method, sample_count=sample_count,
                                      max_optical_depth=maxoptdepth, max_mass_fraction=maxmassfraction, max_dens_disp_fraction=maxdensdispfraction,
                                      direction_method=directionmethod, write=write)

        # Oct tree dust grid
        elif grid.tag == "OctTreeDustGrid":

            write = self.get_boolean(grid, "writeGrid")

            min_x = self.get_quantity(grid, "minX")
            max_x = self.get_quantity(grid, "maxX")
            min_y = self.get_quantity(grid, "minY")
            max_y = self.get_quantity(grid, "maxY")
            min_z = self.get_quantity(grid, "minZ")
            max_z = self.get_quantity(grid, "maxZ")

            min_level = int(grid.get("minLevel"))
            max_level = int(grid.get("maxLevel"))
            search_method = grid.get("searchMethod")
            sample_count = int(grid.get("sampleCount"))
            maxoptdepth = float(grid.get("maxOpticalDepth"))
            maxmassfraction = float(grid.get("maxMassFraction"))
            maxdensdispfraction = float(grid.get("maxDensDispFraction"))
            barycentric = self.get_boolean(grid, "barycentric")

            # Create and return the grid
            return OctTreeDustGrid(min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y, min_z=min_z, max_z=max_z,
                                   min_level=min_level, max_level=max_level, search_method=search_method, sample_count=sample_count,
                                   max_optical_depth=maxoptdepth, max_mass_fraction=maxmassfraction, max_dens_disp_fraction=maxdensdispfraction,
                                   barycentric=barycentric, write=write)

        # File tree dust grid
        elif grid.tag == "FileTreeDustGrid":

            write = self.get_boolean(grid, "writeGrid")

            filename = grid.get("filename")
            search_method = grid.get("searchMethod")

            # Create and return the grid
            return FileTreeDustGrid(filename=filename, search_method=search_method, write=write)

        # Cylindrical grid
        elif grid.tag == "Cylinder2DDustGrid": raise NotImplementedError("Cylindrical grid not supported yet in this function")

        # Invalid
        else: raise NotImplementedError("Other grid types not yet supported")

    ## This function sets the dust grid
    def set_dust_grid(self, grid):

        from .grids import BinaryTreeDustGrid, OctTreeDustGrid, CartesianDustGrid, CylindricalGrid, FileTreeDustGrid

        # Cartesian
        if isinstance(grid, CartesianDustGrid):

            # Set cartesian dust grid
            self.set_cartesian_dust_grid(grid.min_x, grid.max_x, grid.min_y, grid.max_y, grid.min_z, grid.max_z,
                                         grid.x_bins, grid.y_bins, grid.mesh_type, grid.ratio, grid.write)

        # Binary tree
        elif isinstance(grid, BinaryTreeDustGrid):

            # Set binary tree dust grid
            self.set_binary_tree_dust_grid(grid.min_x, grid.max_x, grid.min_y, grid.max_y, grid.min_z, grid.max_z,
                                           grid.write, grid.min_level, grid.max_level, grid.search_method,
                                           grid.sample_count, grid.max_optical_depth, grid.max_mass_fraction,
                                           grid.max_dens_disp_fraction, grid.direction_method, write_tree=grid.write_tree)

        # Octtree
        elif isinstance(grid, OctTreeDustGrid):

            # Set octtree dust grid
            self.set_octtree_dust_grid(grid.min_x, grid.max_x, grid.min_y, grid.max_y, grid.min_z, grid.max_z,
                                       grid.write, grid.min_level, grid.max_level, grid.search_method,
                                       grid.sample_count, grid.max_optical_depth, grid.max_mass_fraction,
                                       grid.max_dens_disp_fraction, grid.barycentric, write_tree=grid.write_tree)

        # Cylindrical
        elif isinstance(grid, CylindricalGrid):

            # Set cylindrical grid
            # self, max_r, min_z, max_z, nbins_r, nbins_z, fraction_r=None, fraction_z=None,
            # ratio_r=None, ratio_z=None, write_grid=False
            self.set_cylindrical_dust_grid(grid.max_r, grid.min_z, grid.max_z, grid.nbins_r, grid.nbins_z,
                                           fraction_r=grid.central_bin_fraction_r, fraction_z=grid.central_bin_fraction_z,
                                           ratio_r=grid.ratio_r, ratio_z=grid.ratio_z, write_grid=grid.write)

        # File tree grid
        elif isinstance(grid, FileTreeDustGrid):

            # Set file tree dust grid
            self.set_filetree_dust_grid(grid.filename, grid.search_method, grid.write)

        # Invalid
        else: raise ValueError("Invalid grid type")

    ## This function removes the dust grid and returns the parent
    def remove_dust_grid(self):

        try:
            # Get the dust grid
            grid = self.get_dust_grid()

            # Get the parent
            parent = grid.getparent()

            # Remove the old grid element
            parent.remove(grid)

        except ValueError:

            dust_system = self.get_dust_system()
            parent = dust_system.makeelement("dustGrid", {"type": "DustGrid"})
            dust_system.append(parent)

        return parent

    ## This function sets a cartesian dust grid for the dust system
    def set_cartesian_dust_grid(self, min_x, max_x, min_y, max_y, min_z, max_z, x_bins, y_bins, z_bins, mesh_type="linear", ratio=1., write_grid=True):

        from ..units.stringify import represent_quantity

        parent = self.remove_dust_grid()

        # Create and add the new grid
        attrs = {"minX": represent_quantity(min_x), "maxX": represent_quantity(max_x), "minY": represent_quantity(min_y),
                 "maxY": represent_quantity(max_y), "minZ": represent_quantity(min_z), "maxZ": represent_quantity(max_z),
                 "writeGrid": str_from_bool(write_grid, lower=True)}
        grid = parent.makeelement("CartesianDustGrid", attrs)
        parent.append(grid)

        # Create the X mesh
        attrs = {"type": "MoveableMesh"}
        x_mesh = grid.makeelement("meshX", attrs)
        grid.append(x_mesh)
        if mesh_type == "linear":
            attrs = {"numBins": str(x_bins)}
            x_mesh.append(x_mesh.makeelement("LinMesh", attrs))
        elif mesh_type == "power":
            attrs = {"numBins": str(x_bins), "ratio": str(ratio)}
            x_mesh.append(x_mesh.makeelement("PowMesh", attrs))
        elif mesh_type == "symmetric_power":
            attrs = {"numBins": str(x_bins), "ratio": str(ratio)}
            x_mesh.append(x_mesh.makeelement("SymPowMesh", attrs))
        else: raise ValueError("Unrecognized mesh type")

        # Create the Y mesh
        attrs = {"type": "MoveableMesh"}
        y_mesh = grid.makeelement("meshY", attrs)
        grid.append(y_mesh)
        if mesh_type == "linear":
            attrs = {"numBins": str(y_bins)}
            y_mesh.append(y_mesh.makeelement("LinMesh", attrs))
        elif mesh_type == "power":
            attrs = {"numBins": str(y_bins), "ratio": str(ratio)}
            y_mesh.append(y_mesh.makeelement("PowMesh", attrs))
        elif mesh_type == "symmetric_power":
            attrs = {"numBins": str(z_bins), "ratio": str(ratio)}
            y_mesh.append(y_mesh.makeelement("SymPowMesh", attrs))
        else: raise ValueError("Unrecognized mesh type")

        # Create the Z mesh
        attrs = {"type": "MovableMesh"}
        z_mesh = grid.makeelement("meshZ", attrs)
        grid.append(z_mesh)
        if mesh_type == "linear":
            attrs = {"numBins": str(z_bins)}
            z_mesh.append(z_mesh.makeelement("LinMesh", attrs))
        elif mesh_type == "power":
            attrs = {"numBins": str(z_bins), "ratio": str(ratio)}
            y_mesh.append(z_mesh.makeelement("PowMesh", attrs))
        elif mesh_type == "symmetric_power":
            attrs = {"numBins": str(z_bins), "ratio": str(ratio)}
            z_mesh.append(z_mesh.makeelement("SymPowMesh", attrs))
        else: raise ValueError("Unrecognized mesh type")

    ## This function sets a binary tree dust grid for the dust system
    def set_binary_tree_dust_grid(self, min_x, max_x, min_y, max_y, min_z, max_z, write_grid=True, min_level=2,
                                  max_level=10, search_method="Neighbor", sample_count=100, max_optical_depth=0,
                                  max_mass_fraction=1e-6, max_dens_disp_fraction=0, direction_method="Alternating",
                                  write_tree=False):

        parent = self.remove_dust_grid()

        # Create and add the new grid
        attrs = {"minX": str(min_x), "maxX": str(max_x), "minY": str(min_y), "maxY": str(max_y), "minZ": str(min_z),
                 "maxZ": str(max_z), "writeGrid": str_from_bool(write_grid, lower=True), "minLevel": str(min_level),
                 "maxLevel": str(max_level), "searchMethod": search_method, "sampleCount": str(sample_count),
                 "maxOpticalDepth": str(max_optical_depth), "maxMassFraction": str(max_mass_fraction),
                 "maxDensDispFraction": str(max_dens_disp_fraction), "directionMethod": direction_method} #"writeTree": str_from_bool(write_tree, lower=True)}

        from ..prep.smile import SKIRTSmileSchema
        smile = SKIRTSmileSchema()
        if smile.supports_file_tree_grids: attrs["writeTree"] = str_from_bool(write_tree, lower=True)

        # Create and add the grid
        parent.append(parent.makeelement("BinTreeDustGrid", attrs))

    ## This function sets the maximal optical depth
    def set_binary_tree_max_optical_depth(self, value):

        # Get the dust grid
        grid = self.get_dust_grid()

        if grid.tag != "BinTreeDustGrid": raise ValueError("The ski file does not specify a binary tree dust grid")

        # Set the optical depth
        self.set_value(grid, "maxOpticalDepth", str(value))

    ## This function sets the maximal mass fraction
    def set_binary_tree_max_mass_fraction(self, value):

        # Get the dust grid
        grid = self.get_dust_grid()

        if grid.tag != "BinTreeDustGrid": raise ValueError("The ski file does not specify a binary tree dust grid")

        # Set the max mass fraction
        self.set_value(grid, "maxMassFraction", str(value))

    ## This function sets an octtree dust grid for the dust system
    def set_octtree_dust_grid(self, min_x, max_x, min_y, max_y, min_z, max_z, write_grid=True, min_level=2,
                              max_level=6, search_method="Neighbor", sample_count=100, max_optical_depth=0,
                              max_mass_fraction=1e-6, max_dens_disp_fraction=0, barycentric=False, write_tree=False):

        parent = self.remove_dust_grid()

        # Create and add the new grid
        attrs = {"minX": str(min_x), "maxX": str(max_x), "minY": str(min_y), "maxY": str(max_y), "minZ": str(min_z),
                 "maxZ": str(max_z), "writeGrid": str_from_bool(write_grid, lower=True), "minLevel": str(min_level),
                 "maxLevel": str(max_level), "searchMethod": search_method, "sampleCount": sample_count,
                 "maxOpticalDepth": str(max_optical_depth), "maxMassFraction": str(max_mass_fraction),
                 "maxDensDispFraction": str(max_dens_disp_fraction), "barycentric": str_from_bool(barycentric, lower=True)}
                 #"writeTree": str_from_bool(write_tree, lower=True)}

        from ..prep.smile import SKIRTSmileSchema
        smile = SKIRTSmileSchema()
        if smile.supports_file_tree_grids: attrs["writeTree"] = str_from_bool(write_tree, lower=True)

        # Create and add the grid
        parent.append(parent.makeelement("OctTreeDustGrid", attrs))

    ## This function sets a cylindrical grid
    def set_cylindrical_dust_grid(self, max_r, min_z, max_z, nbins_r, nbins_z, fraction_r=None, fraction_z=None,
                                  ratio_r=None, ratio_z=None, write_grid=False):

        from ..units.stringify import represent_quantity

        parent = self.remove_dust_grid()

        # Set attrs
        attrs = {"writeGrid": str_from_bool(write_grid, lower=True), "maxR": represent_quantity(max_r), "minZ": represent_quantity(min_z),
                 "maxZ": represent_quantity(max_z)}

        # Create new grid
        grid = parent.makeelement("Cylinder2DDustGrid", attrs)

        # Mesh r
        mesh_r = grid.makeelement("meshR", {"type": "Mesh"})
        mesh_r_log = mesh_r.makeelement("LogMesh", {"numBins": str(nbins_r), "centralBinFraction": repr(fraction_r)})
        mesh_r.append(mesh_r_log)
        grid.append(mesh_r)

        # Mesh z
        mesh_z = grid.makeelement("meshZ", {"type": "MoveableMesh"})
        mesh_z_sympow = mesh_z.makeelement("SymPowMesh", {"numBins": str(nbins_z), "ratio": repr(ratio_z)})
        mesh_z.append(mesh_z_sympow)
        grid.append(mesh_z)

        # Add the grid
        parent.append(grid)

    ## This function sets a file tree dust grid
    def set_filetree_dust_grid(self, filename, search_method="Neighbor", write_grid=False):

        # Remove existing grid
        parent = self.remove_dust_grid()

        # Set attributes
        attrs = {"writeGrid": str_from_bool(write_grid, lower=True), "filename": filename, "searchMethod": search_method}

        # Create new grid
        grid = parent.makeelement("FileTreeDustGrid", attrs)

        # Add the grid
        parent.append(grid)

    ## Range of x in length units
    def get_dust_grid_x_range(self):
        from ..basics.range import QuantityRange
        grid = self.get_dust_grid()
        min_x = self.get_quantity(grid, "minX")
        max_x = self.get_quantity(grid, "maxX")
        return QuantityRange(min_x, max_x)

    ## Range of y in length units
    def get_dust_grid_y_range(self):
        from ..basics.range import QuantityRange
        grid = self.get_dust_grid()
        min_y = self.get_quantity(grid, "minY")
        max_y = self.get_quantity(grid, "maxY")
        return QuantityRange(min_y, max_y)

    ## Range of z in length units
    def get_dust_grid_z_range(self):
        from ..basics.range import QuantityRange
        grid = self.get_dust_grid()
        min_z = self.get_quantity(grid, "minZ")
        max_z = self.get_quantity(grid, "maxZ")
        return QuantityRange(min_z, max_z)

    ## This function returns the instrument system
    def get_instrument_system(self):

        try: return self.get_unique_base_element("instrumentSystem")
        except ValueError:

            # Get parent
            try: parent = self.get_unique_base_element_direct("instrumentSystem")
            except ValueError:
                parent = self.get_simulation().makeelement("instrumentSystem", {"type": "InstrumentSystem"})
                self.get_simulation().append(parent)

            # Set actual instrument system
            system = parent.makeelement("InstrumentSystem")
            parent.append(system)

            # Return the new instrument system
            return system

    ## This funcion removes the complete instrument system
    def remove_instrument_system(self):

        instrument_system = self.get_instrument_system()
        parent = instrument_system.getparent()
        parent.getparent().remove(parent)

    ## This function returns a list of the instruments in the ski file, or the 'instruments' element if as_list is False
    def get_instruments(self, as_list=True):

        # Get the instrument system
        instrument_system = self.get_instrument_system()

        try: instruments = xml.get_unique_element_direct(instrument_system, "instruments")
        except ValueError:
            instruments = instrument_system.makeelement("instruments", {"type": "Instrument"})
            instrument_system.append(instruments)

        # Return the instruments as a list
        if as_list: return instruments.getchildren()
        else: return instruments

    ## This function returns the distance, if possible
    def get_distance(self, return_none=False):
        from ..tools import sequences
        distances = []
        for name in self.get_instrument_names():
            distance = self.get_instrument_distance(name)
            distances.append(distance)
        if len(distances) == 0:
            if return_none: return None
            else: raise ValueError("Distance not defined in ski file")
        if not sequences.all_close(distances):
            if return_none: return None
            else: raise ValueError("Distances are not equal")
        else: return distances[0]

    ## This function returns the names of all the instruments in the ski file as a list
    def get_instrument_names(self):

        # Initialize a list to contain the names
        names = []

        # Get the list of instruments
        instruments = self.get_instruments()

        # Loop over the instruments
        for instrument in instruments:

            # Get the instrument name
            instrument_name = instrument.get("instrumentName")

            # Add the name to the list
            names.append(instrument_name)

        # Return the list of names
        return names

    ## This function removes the instrument with the specified name
    def remove_instrument(self, name):

        # Get the instrument with the specified name
        instrument = self.get_instrument(name)

        # Get element that holds the instrument class
        parent = instrument.getparent()

        # Remove the instrument
        parent.remove(instrument)

    ## This function removes all instruments
    def remove_all_instruments(self):

        for name in self.get_instrument_names():
            self.remove_instrument(name)

    ## This function adds an instrument
    def add_instrument(self, name, instrument):

        from ...modeling.basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument, FullInstrument, FullSEDInstrument, MultiFrameInstrument

        distance = instrument.distance
        inclination = instrument.inclination
        azimuth = instrument.azimuth
        position_angle = instrument.position_angle

        # SED instrument
        if isinstance(instrument, SEDInstrument):

            # Add the SED instrument to the ski file
            self.add_sed_instrument(name, distance, inclination, azimuth, position_angle)

        # Frame instrument
        elif isinstance(instrument, FrameInstrument):

            field_x = instrument.field_x
            field_y = instrument.field_y
            pixels_x = instrument.pixels_x
            pixels_y = instrument.pixels_y
            center_x = instrument.center_x
            center_y = instrument.center_y

            # Add the simple instrument to the ski file
            self.add_frame_instrument(name, distance, inclination, azimuth, position_angle, field_x, field_y, pixels_x, pixels_y, center_x, center_y)

        # Simple instrument
        elif isinstance(instrument, SimpleInstrument):

            field_x = instrument.field_x
            field_y = instrument.field_y
            pixels_x = instrument.pixels_x
            pixels_y = instrument.pixels_y
            center_x = instrument.center_x
            center_y = instrument.center_y

            # Add the simple instrument to the ski file
            self.add_simple_instrument(name, distance, inclination, azimuth, position_angle, field_x, field_y, pixels_x, pixels_y, center_x, center_y)

        # Full instruemnt
        elif isinstance(instrument, FullInstrument):

            field_x = instrument.field_x
            field_y = instrument.field_y
            pixels_x = instrument.pixels_x
            pixels_y = instrument.pixels_y
            center_x = instrument.center_x
            center_y = instrument.center_y
            scattering_levels = instrument.scattering_levels
            counts = instrument.counts

            # Add the full instrument to the ski file
            self.add_full_instrument(name, distance, inclination, azimuth, position_angle, field_x, field_y, pixels_x,
                                     pixels_y, center_x, center_y, scattering_levels, counts)

        # Full SED instrument
        elif isinstance(instrument, FullSEDInstrument):

            # Get the number of scattering levels
            scattering_levels = instrument.scattering_levels

            from ..units.parsing import parse_quantity

            # Set properties
            pixels_x = 1
            pixels_y = 1
            field_x = parse_quantity("100 kpc")
            field_y = parse_quantity("100 kpc")
            center_x = parse_quantity("0 pc")
            center_y = parse_quantity("0 pc")

            # Add the full instrument to the ski file
            self.add_full_instrument(name, distance, inclination, azimuth, position_angle, field_x, field_y, pixels_x,
                                     pixels_y, center_x, center_y, scattering_levels)

        # Multi frame instrument
        elif isinstance(instrument, MultiFrameInstrument):

            # Add the multi frame instrument

            from ..units.stringify import represent_quantity

            # Get the 'instruments' element
            instruments = self.get_instruments(as_list=False)

            distance = instrument.distance
            inclination = instrument.inclination
            azimuth = instrument.azimuth
            position_angle = instrument.position_angle

            # Make and add the new FullInstrument
            attrs = {"instrumentName": name, "distance": represent_quantity(distance),
                     "inclination": str_from_angle(inclination),
                     "azimuth": str_from_angle(azimuth), "positionAngle": str_from_angle(position_angle),
                     "writeTotal": str_from_bool(instrument.write_total, lower=True),
                     "writeStellarComps": str_from_bool(instrument.write_stellar_components, lower=True)}
            instr = instruments.makeelement("MultiFrameInstrument", attrs)

            # Children
            frames = instr.makeelement("frames", {"type": "InstrumentFrame"})

            # Loop over the frames
            for frame in instrument.frames:

                fr_attrs = {"pixelsX": str(frame.pixels_x), "pixelsY": str(frame.pixels_y), "fieldOfViewX": represent_quantity(frame.field_x), "fieldOfViewY": represent_quantity(frame.field_y)}

                fr = frames.makeelement("InstrumentFrame", fr_attrs)

                frames.append(fr)

            # Add the instrument with its frames
            instr.append(frames)
            instruments.append(instr)

        # Unrecognized instrument
        else: raise ValueError("Instruments other than SimpleInstrument, SEDInstrument, FullInstrument, FullSEDInstrument, and MultiFrameInstrument are not yet supported")

    ## This function adds a FrameInstrument to the instrument system
    def add_frame_instrument(self, name, distance, inclination, azimuth, position_angle, field_x, field_y,
                                  pixels_x, pixels_y, center_x, center_y):

        from ..units.stringify import represent_quantity

        # Get the 'instruments' element
        instruments = self.get_instruments(as_list=False)

        # Make and add the new FrameInstrument
        attrs = {"instrumentName": name, "distance": represent_quantity(distance), "inclination": str_from_angle(inclination),
                 "azimuth": str_from_angle(azimuth), "positionAngle": str_from_angle(position_angle),
                 "fieldOfViewX": represent_quantity(field_x), "fieldOfViewY": represent_quantity(field_y), "pixelsX": str(pixels_x),
                 "pixelsY": str(pixels_y), "centerX": represent_quantity(center_x), "centerY": represent_quantity(center_y)}
        instruments.append(instruments.makeelement("FrameInstrument", attrs))

    ## This function adds a FullInstrument to the instrument system
    def add_full_instrument(self, name, distance, inclination, azimuth, position_angle, field_x, field_y,
                            pixels_x, pixels_y, center_x, center_y, scattering_levels=0, counts=False):

        from ..units.stringify import represent_quantity

        # Get the 'instruments' element
        instruments = self.get_instruments(as_list=False)

        # Make and add the new FullInstrument
        attrs = {"instrumentName": name, "distance": represent_quantity(distance), "inclination": str_from_angle(inclination),
                 "azimuth": str_from_angle(azimuth), "positionAngle": str_from_angle(position_angle), "fieldOfViewX": str(field_x),
                 "fieldOfViewY": represent_quantity(field_y), "pixelsX": str(pixels_x), "pixelsY": str(pixels_y),
                 "centerX": represent_quantity(center_x), "centerY": represent_quantity(center_y), "scatteringLevels": str(scattering_levels)}
        if counts: attrs["counts"] = "true"
        instruments.append(instruments.makeelement("FullInstrument", attrs))

    ## This function adds a SimpleInstrument to the instrument system
    def add_simple_instrument(self, name, distance, inclination, azimuth, position_angle, field_x, field_y,
                              pixels_x, pixels_y, center_x, center_y):

        from ..units.stringify import represent_quantity

        # Get the 'instruments' element
        instruments = self.get_instruments(as_list=False)

        # Make and add the new SimpleInstrument
        attrs = {"instrumentName": name, "distance": represent_quantity(distance), "inclination": str_from_angle(inclination),
                 "azimuth": str_from_angle(azimuth), "positionAngle": str_from_angle(position_angle), "fieldOfViewX": represent_quantity(field_x),
                 "fieldOfViewY": represent_quantity(field_y), "pixelsX": str(pixels_x), "pixelsY": str(pixels_y),
                 "centerX": represent_quantity(center_x), "centerY": represent_quantity(center_y)}
        instruments.append(instruments.makeelement("SimpleInstrument", attrs))

    ## This function adds an SEDInstrument to the instrument system
    def add_sed_instrument(self, name, distance, inclination, azimuth, position_angle):

        from ..units.stringify import represent_quantity

        # Get the 'instruments' element
        instruments = self.get_instruments(as_list=False)

        # Make and add the new SEDInstrument
        attrs = {"instrumentName": name, "distance": represent_quantity(distance), "inclination": str_from_angle(inclination),
                 "azimuth": str_from_angle(azimuth), "positionAngle": str_from_angle(position_angle)}
        instruments.append(instruments.makeelement("SEDInstrument", attrs))

    ## This function returns the instrument object
    def get_instrument_object(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Import the instrument classes
        from ...modeling.basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument, FullInstrument, FullSEDInstrument, MultiFrameInstrument

        # Get the instrument
        instrument = self.get_instrument(name)

        # Frame instrument
        if instrument.tag == "FrameInstrument":

            # Get the instrument properties
            distance = self.get_quantity(instrument, "distance")
            inclination = self.get_angle(instrument, "inclination")
            azimuth = self.get_angle(instrument, "azimuth")
            pa = self.get_angle(instrument, "positionAngle")
            fieldx = self.get_quantity(instrument, "fieldOfViewX")
            fieldy = self.get_quantity(instrument, "fieldOfViewY")
            pixelsx = instrument.get("pixelsX")
            pixelsy = instrument.get("pixelsY")
            centerx = self.get_quantity(instrument, "centerX")
            centery = self.get_quantity(instrument, "centerY")

            # Create and return the instrument
            return FrameInstrument(field_x=fieldx, field_y=fieldy, pixels_x=pixelsx, pixels_y=pixelsy, center_x=centerx,
                                   center_y=centery, distance=distance, inclination=inclination, azimuth=azimuth,
                                   position_angle=pa)

        # SED instrument
        elif instrument.tag == "SEDInstrument":

            # Get the instrument properties
            distance = self.get_quantity(instrument, "distance")
            inclination = self.get_angle(instrument, "inclination")
            azimuth = self.get_angle(instrument, "azimuth")
            pa = self.get_angle(instrument, "positionAngle")

            # Create and return the instrument
            return SEDInstrument(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=pa)

        # Simple instrument
        elif instrument.tag == "SimpleInstrument":

            # Get the instrument properties
            distance = self.get_quantity(instrument, "distance")
            inclination = self.get_angle(instrument, "inclination")
            azimuth = self.get_angle(instrument, "azimuth")
            pa = self.get_angle(instrument, "positionAngle")
            fieldx = self.get_quantity(instrument, "fieldOfViewX")
            fieldy = self.get_quantity(instrument, "fieldOfViewY")
            pixelsx = instrument.get("pixelsX")
            pixelsy = instrument.get("pixelsY")
            centerx = self.get_quantity(instrument, "centerX")
            centery = self.get_quantity(instrument, "centerY")

            # Create and return the instrument
            return SimpleInstrument(field_x=fieldx, field_y=fieldy, pixels_x=pixelsx, pixels_y=pixelsy, center_x=centerx,
                                    center_y=centery, distance=distance, inclination=inclination, azimuth=azimuth,
                                    position_angle=pa)

        # Full instrument
        elif instrument.tag == "FullInstrument":

            # Get the instrument properties
            distance = self.get_quantity(instrument, "distance")
            inclination = self.get_angle(instrument, "inclination")
            azimuth = self.get_angle(instrument, "azimuth")
            pa = self.get_angle(instrument, "positionAngle")
            fieldx = self.get_quantity(instrument, "fieldOfViewX")
            fieldy = self.get_quantity(instrument, "fieldOfViewY")
            pixelsx = instrument.get("pixelsX")
            pixelsy = instrument.get("pixelsY")
            centerx = self.get_quantity(instrument, "centerX")
            centery = self.get_quantity(instrument, "centerY")
            scattlevels = int(instrument.get("scatteringLevels"))

            # Is a full SED instrument?
            if pixelsx == 1 and pixelsy == 1:
                return FullSEDInstrument(distance=distance, inclination=inclination, azimuth=azimuth,
                                         position_angle=pa, scattering_levels=scattlevels)

            # Create and return the instrument
            else:
                return FullInstrument(field_x=fieldx, field_y=fieldy, pixels_x=pixelsx, pixels_y=pixelsy, center_x=centerx,
                                      center_y=centery, distance=distance, inclination=inclination, azimuth=azimuth,
                                      position_angle=pa, scattering_levels=scattlevels)

        # MultiFrameInstrument
        elif instrument.tag == "MultiFrameInstrument":

            from ...modeling.basics.instruments import InstrumentFrame

            distance = self.get_quantity(instrument, "distance")
            inclination = self.get_angle(instrument, "inclination")
            azimuth = self.get_angle(instrument, "azimuth")
            pa = self.get_angle(instrument, "positionAngle")

            write_total = self.get_boolean(instrument, "writeTotal")
            write_stellar_components = self.get_boolean(instrument, "writeStellarComps")

            # Create the instrument
            instr = MultiFrameInstrument(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=pa,
                                         write_total=write_total, write_stellar_components=write_stellar_components)

            # Get the frames
            frames = self.get_child_with_name(instrument, "frames")

            # Loop over the frames
            for frame in frames.getchildren():

                # Get properties
                pixelsx = int(self.get_value(frame, "pixelsX"))
                pixelsy = int(self.get_value(frame, "pixelsY"))
                fieldx = self.get_quantity(frame, "fieldOfViewX")
                fieldy = self.get_quantity(frame, "fieldOfViewY")

                # Create the frame
                frm = InstrumentFrame(pixels_x=pixelsx, pixels_y=pixelsy, field_x=fieldx, field_y=fieldy)

                # Add the frame
                instr.add_frame(frm)

            # Return the instrument
            return instr

        # Unrecognized instrument
        else: raise ValueError("Unrecognized instrument: " + instrument.tag)

    ## Get the instrument distance
    def get_instrument_distance(self, name):

        # Get the instrument object
        instrument = self.get_instrument_object(name)
        return instrument.distance

    ## This function returns the instrument with the specified name
    def get_instrument(self, name):

        # Get the list of instruments
        instruments = self.get_instruments()

        # Loop over the instruments
        for instrument in instruments:

            # Get the instrument name
            instrument_name = instrument.get("instrumentName")

            # If the name matches, return
            if name == instrument_name: return instrument

        raise ValueError("No instrument with the name '" + name + "'")

    ## This function changes the name of the specified instrument
    def set_instrument_name(self, old_name, new_name):

        # Get the instrument with the specified name
        instrument = self.get_instrument(old_name)

        # Set the new name
        instrument.set("instrumentName", new_name)

    ## This function returns the distance of the specified instrument as an Astropy quantity
    def get_instrument_distance(self, name):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Return the distance
        return self.get_quantity(instrument, "distance")

    ## This function sets the distance of the specified instruments. The distance should be an Astropy quantity.
    def set_instrument_distance(self, name, value):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Set the distance
        self.set_quantity(instrument, "distance", value)

    ## This function returns the inclination of the specified instrument as an Astropy Angle.
    def get_instrument_inclination(self, name):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Return the inclination
        return self.get_quantity(instrument, "inclination")

    ## This function sets the inclination of the specified instrument. The inclination should be an Astropy Angle or quantity.
    def set_instrument_inclination(self, name, value):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Set the inclination
        self.set_quantity(instrument, "inclination", value)

    ## This function returns the azimuth angle of the specified instrument as an Astropy Angle.
    def get_instrument_azimuth(self, name):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Return the azimuth
        return self.get_quantity(instrument, "azimuth")

    ## This function sets the azimuth angle of the specified instrument. The angle should be an Astropy Angle or quantity.
    def set_instrument_azimuth(self, name, value):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Set the azimuth angle
        self.set_quantity(instrument, "azimuth", value)

    ## This function returns the position angle of the specified instrument as an Astropy Angle.
    def get_instrument_pa(self, name):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Return the position angle
        return self.get_quantity(instrument, "positionAngle")

    ## This function sets the position angle of the specified instrument. The angle should be an Astropy Angle or quantity.
    def set_instrument_pa(self, name, value):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Set the position angle
        self.set_quantity(instrument, "positionAngle", value)

    ## This function sets the orientation of the specified instrument. The angles should be Astropy Angle or Quantity instances.
    def set_instrument_orientation(self, name, inclination, position_angle, azimuth):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Set the inclination
        self.set_quantity(instrument, "inclination", inclination)
        self.set_quantity(instrument, "positionAngle", position_angle)
        self.set_quantity(instrument, "azimuth", azimuth)

    ## This function sets the orientation of the specified instrument to a face-on orientation.
    def set_instrument_orientation_faceon(self, name):

        from astropy.coordinates import Angle

        # XY plane
        inclination = Angle(0., "deg")
        position_angle = Angle(90., "deg")
        azimuth = Angle(0.0, "deg")

        # Set the angles
        self.set_instrument_orientation(name, inclination, position_angle, azimuth)

    ## This function sets the orientation of the specified instrument to an edge-on orientation
    def set_instrument_orientation_edgeon(self, name):

        from astropy.coordinates import Angle

        # XZ plane
        inclination = Angle(90., "deg")
        position_angle = Angle(0., "deg")
        azimuth = Angle(-90., "deg")

        # Set the angles
        self.set_instrument_orientation(name, inclination, position_angle, azimuth)

    ## This function returns the size of the specified instrument as a tuple (size_x, size_y)
    def get_instrument_size(self, name):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Return the size
        return int(instrument.get("pixelsX")), int(instrument.get("pixelsY"))

    ## This function sets the size of the specified instrument
    def set_instrument_size(self, name, x_size, y_size):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Set the size
        instrument.set("pixelsX", str(x_size))
        instrument.set("pixelsY", str(y_size))

    ## This function returns the field of view of the specified instrument as a tuple (field_x, field_y)
    def get_instrument_field(self, name):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Get the field of view
        return self.get_quantity(instrument, "fieldOfViewX"), self.get_quantity(instrument, "fieldOfViewY")

    ## This function sets the field of view of the specified instrument
    def set_instrument_field(self, name, x_field, y_field):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Set the field of view
        self.set_quantity(instrument, "fieldOfViewX", x_field)
        self.set_quantity(instrument, "fieldOfViewY", y_field)

    ## This function gets the center of the specified instrument
    def get_instrument_center(self, name):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        return self.get_quantity(instrument, "centerX"), self.get_quantity(instrument, "centerY")

    ## This function sets the center of the specified instrument
    def set_instrument_center(self, name, x_center, y_center):

        # Get the instrument with this name
        instrument = self.get_instrument(name)

        # Set the center
        self.set_quantity(instrument, "centerX", x_center)
        self.set_quantity(instrument, "centerY", y_center)

    # -----------------------------------------------------------------

    def create_element(self, tag, properties):

        """
        This function ...
        :param tag:
        :param properties:
        :return:
        """

        #from ..tools.stringify import stringify_not_list
        from ..tools.stringify import tostr

        direct_children = []

        children = OrderedDict()
        children_types = OrderedDict()
        attrs = OrderedDict()

        #print(properties)

        if "children" in properties:

            # Loop over the children
            for property_name in properties["children"]:

                #print("PROPERTY NAME:", property_name)
                #print("PARAMETERS:", properties["children"][property_name])

                element = self.create_element(property_name, properties["children"][property_name])

                # Add child
                direct_children.append(element)

        # Loop over the properties, create the children
        for property_name in properties:

            # Children are defined now
            if property_name == "children": continue

                # Loop over the child properties
                #for child_property in properties["children"]: # e.g. geometry, mix, etc.
                #for child_name in properties["children"]:
                #    child_element = self.create_element(child_name, properties["children"][child_name])
                #    children[child_name] = child_element
                #    if "type" in properties["children"][child_name]: children_types[child_name] = properties["children"][child_name]["type"]
                #    children[child_property] = children
                #continue

                # Loop over
                #for child_property_name in properties["children"]

                #print("PROPERTY NAME:", property_name)
                #print("PARAMETERS:", properties["children"][property_name])

                #element = self.create_element(property_name, properties["children"][property_name])

            #print(property_name)
            value = properties[property_name]

            if types.is_tuple(value):

                child = self.create_element(value[0], value[1])
                children[property_name] = [child]

            elif types.is_dictionary(value):

                # Create list of children
                children[property_name] = []

                # Create multiple children
                for key in value:
                    child = self.create_element(key, value[key])
                    if "type" in value[key]:
                        child_type = value[key]["type"]
                        if property_name not in children_types: children_types[property_name] = child_type
                        elif children_types[property_name] != child_type: raise ValueError("Child with type '" + child_type + "' but previous type was '" + children_types[property_name] + "'")
                    children[property_name].append(child)

            # Regular value (string, int, float, quantity)
            else:
                attrs[property_name] = tostr(value, scientific_int=False) #stringify_not_list(value) # can also be 'type'
                #print(property_name, value, type(value))

        # Make element
        # example of attrs:
        #attrs = {"instrumentName": str(index),
        #         "pixelsX": str(pixels[0]), "pixelsY": str(pixels[1]), "width": str(size[0]),
        #         "viewX": str(view[0]), "viewY": str(view[1]), "viewZ": str(view[2]),
        #         "crossX": str(cross[0]), "crossY": str(cross[1]), "crossZ": str(cross[2]),
        #         "upX": str(up[0]), "upY": str(up[1]), "upZ": str(up[2]), "focal": str(focal)}
        #parent.append(parent.makeelement("PerspectiveInstrument", attrs))
        #print(tag)
        #print(attrs)
        element = self.root.makeelement(tag, attrs)

        #print(children)
        #print(children_types)

        # Make children
        for property_name in children:

            if property_name in children_types: attrs = {"type": children_types[property_name]}
            else: attrs = {"type": ""}

            list_element = element.makeelement(property_name, attrs)

            # Add child elements
            for child in children[property_name]:

                # Add to parent
                list_element.append(child)

            # Add the list element to the base element
            element.append(list_element)

        # Add direct children
        for child in direct_children: element.append(child)

        # Return the new element
        return element

    # -----------------------------------------------------------------

    def add_stellar_component(self, properties, title=None):

        """
        This function ...
        :param properties:
        :param title:
        :return:
        """

        tag = "PanStellarComp" if self.panchromatic() else "OligoStellarComp"
        element = self.create_element(tag, properties)

        # Add the component
        components = self.get_stellar_components_object()
        components.append(element)

    # -----------------------------------------------------------------

    def add_dust_component(self, properties, title=None):

        """
        This function ...
        :param properties:
        :param title:
        :return:
        """

        tag = "DustComp"
        element = self.create_element(tag, properties)

        # Add the component
        components = self.get_dust_components_object()
        components.append(element)

    # -----------------------------------------------------------------

    ## This (experimental) function converts the ski file structure into a (nested) python dictionary
    def to_dict(self):
        return xml.recursive_dict(self.tree.getroot())

    ## This (experimental) function converts the ski file structure into json format
    def to_json(self):

        import json
        return json.dumps(self.to_dict())

    ## This function returns the xml tree element with the specified name that is at the base level of the simulation hierarchy
    def get_unique_base_element(self, name):
        return xml.get_unique_element(self.get_simulation(), "//"+name)

    ## This function ...
    def get_unique_base_element_direct(self, name):
        return xml.get_unique_element_direct(self.get_simulation(), "//"+name)

    # -----------------------------------------------------------------

    ## This function returns the boolean value of a certain parameter of the specified tree elemtn
    def get_boolean(self, element, name):
        from ..tools.parsing import boolean
        return boolean(self.get_value(element, name))

    # -----------------------------------------------------------------

    ## This functions returns the value of a certain parameter of the specified tree element as an Astropy Angle.
    def get_angle(self, element, name, default_unit=None):

        # Import
        from astropy.coordinates import Angle

        # Parse as quantity and then convert to Angle
        quantity = self.get_quantity(element, name, default_unit=default_unit)
        return Angle(quantity.value, quantity.unit)

    # -----------------------------------------------------------------

    ## This function returns the child of the passed element with a given name
    def get_child_with_name(self, element, child_name):
        for child in element.getchildren():
            if child.tag == child_name: return child
        return None

    def get_value_for_path(self, property_path, element=None):

        # Determine the name of the property
        property_name = property_path.split("/")[-1]

        # Start with the simulation object as the root for searching through the simulation hierarchy
        if element is None: last = self.get_simulation()

        # Or start with the passed element
        else: last = element

        # Find the property
        for link in property_path.split("/")[:-1]:
            last = self.get_child_with_name(last, link)

        # Get the current property value
        return self.get_value(last, property_name)

    def set_value_for_path(self, property_path, value, element=None):

        # Determine the name of the property
        property_name = property_path.split("/")[-1]

        # Start with the simulation object as the root for searching through the simulation hierarchy
        if element is None: last = self.get_simulation()

        # Or start with the passed element
        else: last = element

        # Find the property
        for link in property_path.split("/")[:-1]:
            last = self.get_child_with_name(last, link)

        # Set the property value
        self.set_value(last, property_name, value)

    # Upgrade the ski file to the latest SKIRT version
    def upgrade(self):

        from ..prep.upgradeskifile import _get_upgrade_definitions

        # Get the upgrade conditions and transforms
        upgrades = _get_upgrade_definitions()

        # Apply the upgrades if needed
        changed = False
        for condition, templates in upgrades:
            changed |= self.transformif(condition, templates)

    # -----------------------------------------------------------------

    ## FROM HERE ON, FROM OLD LABELED SKI FILE CLASS

    @property
    def labels(self):

        """
        This function returns all labels
        :return:
        """
        from .skifile import get_label, is_labeled

        labels = set()

        # Loop over all elements in the tree
        for element in self.tree.getiterator():

            # Loop over the settings of the element
            for setting_name, setting_value in element.items():

                # Get label, if labeled
                if is_labeled(setting_value): labels.add(get_label(setting_value))

        # Return the list of labels
        return list(labels)

    # -----------------------------------------------------------------

    def delabel_all(self):

        """
        This function ...
        :return:
        """

        # Loop over the labels, delabel
        for label in self.labels: self.delabel(label)

    # -----------------------------------------------------------------

    def delabel_all_except(self, *args):

        """
        This function ...
        :param args
        :return:
        """

        # Loop over the labels, delabel
        for label in self.labels:

            # Skip if label is passed as argument
            if label in args: continue

            # Otherwise, delabel
            self.delabel(label)

    # -----------------------------------------------------------------

    def delabel(self, label):

        """
        This function removes the label from a certain property
        :param label:
        :return:
        """

        # Loop over all elements in the tree
        for element in self.tree.getiterator():

            # Loop over the settings of the element
            for setting_name, setting_value in element.items():

                if setting_value.startswith("[") and setting_value.endswith("]"):

                    label_item = setting_value.split("[")[1].split(":")[0]
                    value_item = setting_value[1:-1].split(":")[1]

                    if label == label_item: element.set(setting_name, value_item)

    # -----------------------------------------------------------------

    def set_labeled_value(self, label, value):

        """
        This function ...
        :param label:
        :param value:
        :return:
        """

        from ..tools.stringify import stringify_not_list

        if label not in self.labels: raise ValueError("The label '" + label + "' is not present in the ski file")

        # Get the labeled elements
        elements = self.get_labeled_elements(label)

        # Labeled value
        labeled_value = "[" + label + ":" + stringify_not_list(value)[1] + "]"

        # Set the new value for each corresponding element
        for element, setting_name in elements: element.set(setting_name, labeled_value)

    # -----------------------------------------------------------------

    def set_labeled_values(self, values):

        """
        This function ...
        :param values: a dictionary, with the keys a subset of the labels in the ski file
        :return:
        """

        from ..tools.stringify import stringify_not_list

        existing_labels = self.labels

        # Loop over the labels
        for label in values:

            # Check for label existence
            if label not in existing_labels: raise ValueError("The label '" + label + "' is not present in the ski file")

            # Get the labeled elements
            elements = self.get_labeled_elements(label)

            # Set the new value for each corresponding element
            for element, setting_name in elements:

                # Convert the value into a string
                if (element.tag, setting_name) in fake_quantities: string = repr(values[label].to(fake_quantities[(element.tag, setting_name)]).value)
                else: string = stringify_not_list(values[label])[1]

                # Label the value string and set it in the ski file
                labeled_value = "[" + label + ":" + string + "]"
                element.set(setting_name, labeled_value)

    # -----------------------------------------------------------------

    def get_labeled_elements(self, label):

        """
        This function ...
        :return:
        """
        from .skifile import get_label, is_labeled

        # The list of elements
        elements = []

        # Loop over all elements in the tree
        for element in self.tree.getiterator():

            # Loop over the settings of the element
            for setting_name, setting_value in element.items():

                if is_labeled(setting_value):
                    label_item = get_label(setting_value)
                    if label_item == label: elements.append((element, setting_name))

        # Return the list of elements
        return elements

    # -----------------------------------------------------------------

    def get_labeled_value(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Get all elements with this value
        elements = self.get_labeled_elements(label)

        the_value = None

        # Loop over the elements
        for element, setting_name in elements:

            #print(setting_name)
            #print(element.get(setting_name))

            # Get the value
            value = self.get_quantity(element, setting_name)

            if the_value is None: the_value = value
            elif the_value != value: raise ValueError("The '" + label + "' property has different values throughout the ski file (" + str(value) + " and " + str(the_value) + ")")

        # Return the value
        return the_value

    # -----------------------------------------------------------------

    def get_labeled_values(self):

        """
        This function ...
        :return:
        """

        from .skifile import get_label, is_labeled

        values = dict()

        # Loop over all elements in the tree
        for element in self.tree.getiterator():

            # Loop over the settings of the element
            for setting_name, setting_value in element.items():

                if is_labeled(setting_value):

                    label = get_label(setting_value)
                    value = self.get_quantity(element, setting_name)

                    if label in values and values[label] != value:
                        warnings.warn("The '" + label + "' property has different values throughout the SkiFile (" + str(values[label]) + " and " + str(value) + ")")
                    else: values[label] = value

        # Return the dictionary of values
        return values

    # -----------------------------------------------------------------

    def get_element_for_path(self, property_path, element=None):

        """
        This function ...
        :param property_path:
        :param element:
        :return:
        """

        from ..tools.numbers import is_number
        from ..tools.strings import is_quoted, unquote

        # Start with the simulation object as the root for searching through the simulation hierarchy
        if element is None: last = self.get_simulation()

        # Or start with the passed element
        else: last = element

        # Find the property to be set
        for link in property_path.split("/"):

            # Bracket specification
            if link.startswith("[") and link.endswith("]"):

                specification = link.split("[")[1].split("]")[0]

                # All children
                if specification == ":":

                    elements = self.get_children(last)
                    if len(elements) == 0: raise ValueError("The path defines to existing elements")
                    elif len(elements) > 1: raise ValueError("The path defines multiple elements")
                    else: last = elements[0]

                # Child with certain index
                elif is_number(specification):

                    index = int(specification)
                    children = self.get_children(last)
                    child = children[index]
                    #elements = [child]
                    last = child

                # Child with certain title
                elif is_quoted(specification):

                    title = unquote(specification)
                    children = self.get_ids_children(last)
                    child = None
                    for key in children:
                        if key == title:
                            child = children[key]
                            break
                    if child is None: raise RuntimeError("Something went wrong")
                    #elements = [child]
                    last = child

                # Invalid specification
                else: raise ValueError("Invalid specification: '" + specification + "'")

            # REGULAR SPECIFICATION
            else: last = self.get_child_with_name(last, link)

        # Return the last element
        return last

    # -----------------------------------------------------------------

    def add_label_to_path(self, property_path, label, element=None, allow_change=True):

        """
        This function ...
        :param property_path:
        :param label:
        :param element:
        :param allow_change:
        :return:
        """
        from .skifile import get_label, is_labeled

        # Determine the name of the property
        property_name = property_path.split("/")[-1]
        property_element_path = "/".join(property_path.split("/")[:-1])

        # Start with the simulation object as the root for searching through the simulation hierarchy
        if element is None: last = self.get_simulation()

        # Or start with the passed element
        else: last = element

        # Start with the base element
        elements = self.get_elements(last, property_element_path)

        #print(elements)

        # Loop over the elements
        for el in elements:

            # Get the current property value
            value = el.get(property_name)

            # If labeled, strip the label
            if is_labeled(value):
                old_label = get_label(value)
                value = self.get_value(el, property_name)
                if not allow_change and old_label != label: raise ValueError("The already present label '" + old_label + "' does not correspond to the new label '" + label + "'")

            # Add the label
            value = "[" + label + ":" + value + "]"

            # Set the property with the label
            el.set(property_name, value)

    # -----------------------------------------------------------------

    def get_ids_children(self, element):

        """
        This function ...
        :param element:
        :return:
        """

        # Get all children
        children = element.getchildren()

        # Get the child elements
        child_elements = [child for child in children if child.tag is not etree.Comment]

        # Loop over the components (also includes the comments)
        number_of_components = 0
        i = 0
        ids = []
        while i < len(children):
            if children[i].tag is etree.Comment:
                ids.append(children[i].text.strip())
                i += 2  # skip the next component -> it is the component corresponding to this comment
            # No name -> add the index of this component as the ID
            else:
                ids.append(number_of_components)
                i += 1
            # Increment the number of components
            number_of_components += 1

        # Create resulting dictionary
        result = OrderedDict()
        for id, element in zip(ids, child_elements):
            result[id] = element

        # Return the dictionary
        return result

    # -----------------------------------------------------------------

    def get_children(self, element, include_comments=False):

        """
        This function ...
        :param element:
        :param include_comments:
        :return:
        """

        # Return the stellar components as a list
        if include_comments: return element.getchildren()
        else: return [child for child in element.getchildren() if child.tag is not etree.Comment]

    # -----------------------------------------------------------------

    def add_labels(self, paths, allow_change=True):

        """
        This function ...
        :param paths:
        :param allow_change:
        :return:
        """

        # Loop over the labels
        for label in paths:

            # Loop over the paths, add the labels
            label_paths = paths[label]
            for label_path in label_paths: self.add_label_to_path(label_path, label, allow_change=allow_change)

    # -----------------------------------------------------------------

    def fix_labels(self, base_path):

        """
        This function ...
        :param base_path:
        :return:
        """

        paths = self.get_labels_and_paths(merge=False)

        # Loop over all labels
        for label in paths:

            label_paths = paths[label]
            if len(label_paths) == 1: continue # We cannot fix anything, only one so no comparison

            reference_label_path = None

            for label_path in label_paths:
                if not label_path.startswith(base_path):
                    reference_label_path = label_path
                    break
            if reference_label_path is None:
                #raise ValueError("No comparison for the '" + label + "' labeled property")
                continue # We CANNOT FIX ANYTHING

            reference_property_name = reference_label_path.split("/")[-1]
            reference_base_path = "/".join(reference_label_path.split("/")[:-1])
            #print(reference_base_path)
            reference_element = self.get_element_for_path(reference_base_path)
            reference_value = self.get_value(reference_element, reference_property_name)

            # Fix
            for label_path in label_paths:

                if not label_path.startswith(base_path): continue

                property_name = label_path.split("/")[-1]
                base_path = "/".join(label_path.split("/")[:-1])
                element = self.get_element_for_path(base_path)
                self.set_value(element, property_name, reference_value)

    # -----------------------------------------------------------------

    def get_elements_for_label(self, label):

        """
        Thisf unction ...
        :param label:
        :return:
        """

        # Initialize a list for the elements
        elements = []

        # Loop over all elements in the tree
        for element in self.tree.getiterator():

            # Loop over the settings of the element
            for setting_name, setting_value in element.items():

                # This is a labeled value
                if setting_value.startswith("[") and setting_value.endswith("]"):

                    # Get the label
                    label_i = setting_value.split("[")[1].split(":")[0]

                    # Add the element if the label corresponds
                    if label_i == label: elements.append((element, setting_name))

        # Return the list of elements
        return elements

    # -----------------------------------------------------------------

    def get_paths_for_label(self, label, merge=True):

        """
        This function ...
        :param label:
        :param merge:
        :return:
        """

        from .skifile import get_element_path

        # Initialize a list for the different paths
        paths = []

        # Loop over the elements (with setting name) that contain the label
        for element, setting_name in self.get_elements_for_label(label):

            # Get element path
            path = get_element_path(element)

            # Make relative to simulation
            path = path.split(self.simulation_class + "/")[1]

            # Add the setting name
            path += "/" + setting_name

            # Add the path
            paths.append(path)

        # MERGE?
        if merge:

            from ..tools.numbers import find_numbers_in_string
            from ..tools.strings import split_cumulative
            from ..tools.sequences import same_contents

            # Check the paths, if some can be merged together
            merge_candidates = []
            for path in paths:
                if "[" not in path: continue
                if "]" not in path: continue
                numbers = find_numbers_in_string(path, left="[", right="]")
                if len(numbers) > 0: merge_candidates.append(path)

            # More than one candidate
            if len(merge_candidates) > 1:
                merged = []
                merged_indices = []  # merged indices
                for index in range(len(merge_candidates)):

                    # Skip already merged
                    if index in merged_indices: continue

                    path = merge_candidates[index]
                    parts = split_cumulative(path, "[", include_total=False)

                    merged_with_index = []

                    for other_index in range(len(merge_candidates)):

                        # Skip already merged
                        if other_index in merged_indices: continue

                        # Do not check yourself
                        if other_index == index: continue

                        other_path = merge_candidates[other_index]
                        for part in parts:
                            if part in other_path:
                                #merged_with_index.append(other_path)
                                merged_with_index.append(other_index)
                                merged_indices.append(other_index)
                                break

                    if len(merged_with_index) > 0:
                        #merged_with_index.append(path)
                        merged_with_index.append(index)
                        merged_indices.append(index)
                        merged.append(merged_with_index)

                # Check whether the merged paths span the entire list (e.g. , xxx[0]yyy, xxx[1]yyy, xxx[2]yyy until [n-1])
                for mergers in merged:
                    first = merge_candidates[mergers[0]]
                    #first = mergers[0]
                    common = first.split("/[")[0] # HERE WE ASSUME THERE IS ONLY ONE [I] !!!
                    nchildren = len(self.get_children_for_path(common))
                    #print(common, nchildren)
                    merger_indices = [int(merge_candidates[index].split("/[")[1].split("]")[0]) for index in mergers]
                    #print(nchildren, merger_indices)
                    # CHECK IF THEY'RE ALL THERE
                    needed_indices = range(nchildren)

                    # SUCCES: MERGE
                    if same_contents(needed_indices, merger_indices):

                        start = common + "/["
                        end = "]" + first.split("/[")[1].split("]")[1]
                        merged_path = start + ":" + end
                        #print(merged_path)

                        # REMOVE THE MERGER INDICES FROM THE PATHS
                        for merger_index in sorted(mergers, reverse=True):
                            del paths[merger_index]

                        # Add the merged path
                        paths.append(merged_path)

        # Return the paths
        return paths

    # -----------------------------------------------------------------

    def get_children_for_path(self, path, include_comments=False):

        """
        This function ...
        :param path:
        :param include_comments:
        :return:
        """

        element = self.get_element_for_path(path)
        return self.get_children(element, include_comments=include_comments)

    # -----------------------------------------------------------------

    def get_labels_and_paths(self, merge=True):

        """
        This function ...
        :param merge:
        :return:
        """

        # Initialize a dictionary to contain the paths
        paths = dict()

        # Loop over the labels
        for label in self.labels:

            # Get the paths
            label_paths = self.get_paths_for_label(label, merge=merge)

            # Set the paths
            paths[label] = label_paths

        # Return the dictionary
        return paths

    # -----------------------------------------------------------------

    def get_value(self, element, name):

        """
        This function returns the string value of a property of an element
        :param element:
        :param name:
        :return:
        """

        if name not in element.attrib: raise ValueError("A property '" + name + "' does not exist for this element")

        prop = element.get(name)
        if prop.startswith("[") and prop.endswith("]"): prop = prop[1:-1].split(":")[1].strip()

        return prop

    # -----------------------------------------------------------------

    def set_value(self, element, name, string):

        """
        This function sets the string value of a property of an element
        :param element:
        :param name:
        :param string:
        :return:
        """
        from .skifile import is_labeled

        # Get the property
        prop = element.get(name)

        # Not a new property
        if prop is not None:
            # Labeled value, add label to stringified quantity
            if is_labeled(prop):
                label = prop[1:-1].split(":")[0]
                string = "[" + label + ":" + string + "]"

        # Set the value in the tree element
        element.set(name, string)

    # -----------------------------------------------------------------

    def get_quantity(self, element, name, default_unit=None):

        """
        This function returns the value of a certain parameter of the specified tree element as an Astropy quantity. The
        default unit can be specified which is used when the unit is not described in the ski file.
        :param element:
        :param name:
        :param default_unit:
        :return:
        """

        # Import Astropy here to avoid import errors for this module for users without an Astropy installation
        from ..units.parsing import parse_quantity
        from ..units.parsing import parse_unit

        # Get string value
        prop = self.get_value(element, name)

        try:
            quantity = parse_quantity(prop)
            if quantity.unit == "" and default_unit is not None:
                quantity = quantity * parse_unit(default_unit)
            return quantity
        except ValueError:
            value = float(prop)
            if default_unit is not None: quantity = value * parse_unit(default_unit)
            else: quantity = value
            return quantity

    # -----------------------------------------------------------------

    def set_quantity(self, element, name, value, default_unit=None):

        """
        This function sets the value of a certain parameter of the specified tree element from an Astropy quantity.
        :param element:
        :param name:
        :param value:
        :param default_unit:
        :return:
        """

        # Import here to avoid import errors for this module for users without an Astropy installation
        from ..units.stringify import represent_quantity, represent_unit
        from ..units.parsing import parse_unit

        # Stringify the value
        if hasattr(value, "unit"): string = represent_quantity(value)
        elif default_unit is not None: string = repr(value) + " " + represent_unit(parse_unit(default_unit))
        else: string = repr(value)

        # Set the value
        self.set_value(element, name, string)

    # -----------------------------------------------------------------

    def get_elements(self, base, path):

        """
        This function ...
        :param base:
        :param path:
        :return:
        """

        elements = [base]

        # Find the property to be set
        for link in path.split("/"):

            # Make new elements
            new_elements = []
            for element in elements: new_elements.extend(self.get_next_elements(element, link))

            # Replace
            elements = new_elements

        # Return
        return elements

    # -----------------------------------------------------------------

    def get_next_elements(self, base, link):

        """
        This function ...
        :param base: 
        :param link: 
        :return: 
        """

        from ..tools.numbers import is_number
        from ..tools.strings import is_quoted, unquote

        # Bracket specification
        if link.startswith("[") and link.endswith("]"):

            specification = link.split("[")[1].split("]")[0]

            # All children
            if specification == ":": elements = self.get_children(base)

            # Child with certain index
            elif is_number(specification):

                index = int(specification)
                children = self.get_children(base)
                child = children[index]
                elements = [child]

            # Child with certain title
            elif is_quoted(specification):

                title = unquote(specification)
                children = self.get_ids_children(base)
                child = None
                for key in children:
                    if key == title:
                        child = children[key]
                        break
                if child is None: raise RuntimeError("Something went wrong")
                elements = [child]

            # Invalid specification
            else: raise ValueError("Invalid specification: '" + specification + "'")

        # Regular specification
        else:
            last = self.get_child_with_name(base, link)
            if last is None: raise ValueError("Could not find the element with the link '" + link + "'")
            elements = [last]

        # Return the elements
        return elements

# -----------------------------------------------------------------
