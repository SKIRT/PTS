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
from lxml import etree
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..tools import types
from .skifile7 import SkiFile7
from .skifile8 import SkiFile8
from ..tools.introspection import skirt_main_version, has_skirt

# -----------------------------------------------------------------

# Get the SKIRT version number
if not has_skirt(): version_number = 8
else: version_number = skirt_main_version()

# -----------------------------------------------------------------

# Set the appropriate ski file class, depending on the SKIRT version
if version_number == 8: SkiFile = SkiFile8
elif version_number == 7: SkiFile = SkiFile7
else: raise ValueError("Invalid SKIRT version: " + str(version_number))

# -----------------------------------------------------------------

def fix_ski_file(path):

    """
    This function ...
    :param path:
    :return:
    """

    from ...core.prep.smile import expected_types

    # Set replacement dictionary
    replacements = dict()
    for a, b in expected_types.items():
        from_string = a + ' type=""'
        to_string = a + ' type="' + b + '"'
        replacements[from_string] = to_string

    # Replace lines
    replace_ski_file_lines(path, replacements)

# -----------------------------------------------------------------

def replace_ski_file_lines(path, replacement_dict):

    """
    This function ...
    :param path:
    :param replacement_dict:
    :return:
    """

    new_lines = []

    which_system = None

    # Read the lines
    for line in fs.read_lines(path):

        if "<dustSystem" in line: which_system = "dust"
        elif "/dustSystem" in line: which_system = None

        if "<stellarSystem" in line: which_system = "stellar"
        elif "/stellarSystem" in line: which_system = None

        # Loop over the replacements
        for from_string in replacement_dict:
            to_string = replacement_dict[from_string]

            # Determine the new line
            if from_string in line: line = line.replace(from_string, to_string)

            if 'components type=""' in line:
                if which_system == "dust": line = line.replace('components type=""', 'components type="DustComp"')
                elif which_system == "stellar": line = line.replace('components type=""', 'components type="StellarComp"')
                else: raise RuntimeError("Something went wrong")

        # Add the line
        new_lines.append(line)

    # Remove the file
    fs.remove_file(path)

    # Write lines
    fs.write_lines(path, new_lines)

# -----------------------------------------------------------------

def get_element_path(element):

    """
    This function ...
    :param element:
    :return:
    """

    from ..tools import sequences

    name = element.tag

    parent = element.getparent()
    if parent is not None:

        # Get the children of the parent
        children = parent.getchildren()

        # The parent has multiple children
        nchildren = len(children)
        if nchildren > 1:

            # Get the child labels (names)
            labels = [child.tag for child in children if child.tag is not etree.Comment]
            nactual_children = len(labels)

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

            #print(ids)

            # Check the index of this element
            index = 0
            actual_index = 0
            while True:
                if children[actual_index].tag is etree.Comment:
                    actual_index += 1
                    continue
                if children[actual_index] == element: break
                index += 1
                actual_index += 1
                if actual_index == nchildren:
                    index = None
                    break
            if index is None: raise RuntimeError("Something went wrong")

            #print(labels)
            #print(index)

            # Get the label and the ID
            element_label = labels[index]
            element_id = ids[index]

            if types.is_integer_type(element_id):
                if element_id != index: raise RuntimeError("Something went wrong")
                identifier = str(index)
            elif types.is_string_type(element_id): identifier = "'" + element_id + "'"
            else: raise RuntimeError("Something went wrong")

            if sequences.is_unique(labels, element_label): path = get_element_path(parent) + "/" + name
            else: path = get_element_path(parent) + "/[" + identifier + "]"

        # No multiple children
        else: path = get_element_path(parent) + "/" + name

    # No parent
    else: path = name

    # Return the total path
    return path

# -----------------------------------------------------------------

def is_labeled(value):

    """
    This function ...
    :param value:
    :return:
    """

    return value.startswith("[") and value.endswith("]")

# -----------------------------------------------------------------

def get_label(value):

    """
    This function ...
    :param value:
    :return:
    """

    label = value.split("[")[1].split(":")[0]
    return label

#-----------------------------------------------------------------

def get_normalizations_from_file(ski_path):

    """
    This function ...
    :param ski_path:
    :return:
    """

    # Load ski file
    ski = SkiFile(ski_path)

    # Get normalizations
    return get_normalizations(ski)

# -----------------------------------------------------------------

def get_normalizations(ski):

    """
    This function ...
    :param ski:
    :return:
    """

    from ..filter.filter import Filter

    # Get stellar component IDs
    stellar_component_ids = ski.get_stellar_component_ids()

    # Initialize dictionaries for wavelength and filter
    wavelengths = defaultdict(dict)
    filters = defaultdict(dict)

    # Loop over the stellar components
    for component_id in stellar_component_ids:

        # Get normalization filter/wavelength and luminosity
        filter_or_wavelength = ski.get_stellar_component_normalization_wavelength_or_filter(component_id)
        luminosity = ski.get_stellar_component_luminosity(component_id, return_wavelength=False)

        # Add to dict
        if isinstance(filter_or_wavelength, Filter): filters[filter_or_wavelength][component_id] = luminosity
        else: wavelengths[filter_or_wavelength][component_id] = luminosity

    # Return normalizations
    return wavelengths, filters

# -----------------------------------------------------------------

def show_normalizations(ski, flux_unit=None):

    """
    This function ...
    :param ski:
    :param flux_unit:
    :return:
    """

    # Import tools
    from ..tools import formatting as fmt
    from ..tools.stringify import tostr

    # Get normalizations
    wavelengths, filters = get_normalizations(ski)

    nwavelengths = len(wavelengths)
    nfilters = len(filters)
    has_wavelengths = nwavelengths > 0
    has_filters = nfilters > 0

    # Get distance (returns None if couldn't define)
    distance = ski.get_distance(return_none=True)

    # Show wavelength normalizations
    if has_wavelengths:

        print("")
        print(fmt.underlined + fmt.blue + "Wavelength normalizations:" + fmt.reset)
        print("")

        # Loop over the wavelengths
        for wavelength in wavelengths:

            print(fmt.bold + fmt.green + tostr(wavelength) + fmt.reset + ":")
            print("")

            total_luminosity = 0.

            # Loop over the components
            for component_id in wavelengths[wavelength]:

                luminosity = wavelengths[wavelength][component_id]
                if flux_unit is not None: flux = luminosity.to(flux_unit, distance=distance, wavelength=wavelength)
                else: flux = None
                if flux is not None: print(" - " + fmt.bold + component_id + fmt.reset + ": " + tostr(luminosity) + " (" + tostr(flux) + ")")
                else: print(" - " + fmt.bold + component_id + fmt.reset + ": " + tostr(luminosity))
                total_luminosity += luminosity

            if flux_unit is not None: total_flux = total_luminosity.to(flux_unit, distance=distance, wavelength=wavelength)
            else: total_flux = None

            print("")
            if total_flux is not None: print(" - TOTAL: " + tostr(total_luminosity) + " (" + tostr(total_flux) + ")")
            else: print(" - TOTAL: " + tostr(total_luminosity))
            print("")

    # Show filter normalizations
    if has_filters:

        print("")
        print(fmt.underlined + fmt.blue + "Filter normalizations:" + fmt.reset)
        print("")

        # Loop over the filters
        for fltr in filters:

            print(fmt.bold + fmt.green + str(fltr) + fmt.reset + ":")
            print("")

            # Get wavelength
            wavelength = fltr.wavelength

            total_luminosity = 0.

            # Loop over the components
            for component_id in filters[fltr]:

                luminosity = filters[fltr][component_id]
                if flux_unit is not None: flux = luminosity.to(flux_unit, distance=distance, wavelength=wavelength)
                else: flux = None
                if flux is not None: print(" - " + fmt.bold + component_id + fmt.reset + ": " + tostr(luminosity) + " (" + tostr(flux) + ")")
                else: print(" - " + fmt.bold + component_id + fmt.reset + ": " + tostr(luminosity))
                total_luminosity += luminosity

            if flux_unit is not None: total_flux = total_luminosity.to(flux_unit, distance=distance, wavelength=wavelength)
            else: total_flux = None

            print("")
            if total_flux is not None: print(" - TOTAL: " + tostr(total_luminosity) + " (" + tostr(total_flux) + ")")
            else: print(" - TOTAL: " + tostr(total_luminosity))
            print("")

# -----------------------------------------------------------------
