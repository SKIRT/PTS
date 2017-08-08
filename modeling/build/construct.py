#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.construct Contains functions to construct ski files from model definitions.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.filter.filter import parse_filter

# -----------------------------------------------------------------

def add_stellar_component(ski, name, component):

    """
    This function ...
    :param ski: 
    :param name:
    :param component: 
    :return: 
    """

    # Debugging
    log.debug("Adding stellar component '" + name + "' to the ski file ...")

    # If an input map is required
    if "map_path" in component: filename = set_stellar_input_map(name, component)
    else: filename = None

    # Set geometry
    if "model" in component: set_stellar_component_model(ski, component)

    # Set deprojection
    elif "deprojection" in component: set_stellar_component_deprojection(ski, component)

    # From parameters
    if component.parameters is not None:

        # Check if this is a new component, add geometry, SED and normalization all at once
        if "geometry" in component.parameters: set_stellar_component_geometry_sed_and_normalization(ski, component)

        # Existing component, with MAPPINGS template
        elif "sfr" in component.parameters: set_stellar_component_mappings(ski, component)

        # Existing component, no MAPPINGS
        else: set_stellar_component(ski, component)

    # From properties
    if component.properties is not None:

        # Add component
        ski.add_stellar_component(component.properties, name)

    # Return the map filename
    return filename

# -----------------------------------------------------------------

def set_stellar_input_map(name, component):

    """
    This function ...
    :param name:
    :param component:
    :return: 
    """

    # Generate a filename for the map
    map_filename = "stars_" + name + ".fits"

    # Set the filename
    if "deprojection" in component: component.deprojection.filename = map_filename
    #elif "geometry" in component.parameters: component.properties["geometry"].filename = filename # HOW IT WAS
    elif component.parameters is not None and "geometry" in component.parameters: component.parameters["geometry"].filename = map_filename
    elif component.properties is not None: component.properties["children"]["geometry"]["children"]["ReadFitsGeometry"]["filename"] = map_filename
    else: raise RuntimeError("Stellar component based on an input map should either have a deprojection or geometry properties")

    # Add entry to the input maps dictionary
    #self.input_map_paths[filename] = component.map_path

    return map_filename

# -----------------------------------------------------------------

def set_stellar_component_model(ski, component):

    """
    This function ...
    :param ski: 
    :param component:
    :return: 
    """

    # Get title
    title = component.parameters.title

    # Set the geometry
    ski.set_stellar_component_geometry(title, component.model)

# -----------------------------------------------------------------

def set_stellar_component_deprojection(ski, component):

    """
    THis function ...
    :param ski:
    :param component:
    :return: 
    """

    # Get title
    title = component.parameters.title

    # Set the deprojection geometry
    ski.set_stellar_component_geometry(title, component.deprojection)

# -----------------------------------------------------------------

def set_stellar_component_geometry_sed_and_normalization(ski, component):

    """
    This function ...
    :param ski:
    :param component:
    :return: 
    """

    # Get title
    title = component.parameters.title

    # Get class names
    geometry_type = component.parameters.geometry
    sed_type = component.parameters.sed
    normalization_type = component.parameters.normalization

    # Get properties for each of the three classes
    geometry_properties = component.properties["geometry"]
    sed_properties = component.properties["sed"]
    normalization_properties = component.properties["normalization"]

    # Create stellar component
    ski.create_new_stellar_component(title, geometry_type, geometry_properties, sed_type, sed_properties,
                                          normalization_type, normalization_properties)

# -----------------------------------------------------------------

def set_stellar_component_mappings(ski, component):

    """
    THis function ...
    :param ski:
    :param component:
    :return: 
    """

    # Get title
    title = component.parameters.title

    # Get SED properties
    metallicity = component.parameters.metallicity
    compactness = component.parameters.compactness
    pressure = component.parameters.pressure
    covering_factor = component.parameters.covering_factor

    # Get normalization
    fltr = parse_filter(component.parameters.filter)
    luminosity = component.parameters.luminosity

    # Set SED
    ski.set_stellar_component_mappingssed(title, metallicity, compactness, pressure, covering_factor)  # SED

    # Set center wavelength of the filter as normalization wavelength (keeps label)
    ski.set_stellar_component_normalization_wavelength(title, fltr.center)

    # Set spectral luminosity at that wavelength (keeps label)
    ski.set_stellar_component_luminosity(title, luminosity)

    # Scale height doesn't need to be set as parameter, this is already in the deprojection model

# -----------------------------------------------------------------

def set_stellar_component(ski, component):

    """
    This function ...
    :return: 
    :param ski:
    :param component:
    """

    # Get title
    title = component.parameters.title

    # Get SED properties
    template = component.parameters.template
    age = component.parameters.age
    metallicity = component.parameters.metallicity

    # Get normalization
    fltr = parse_filter(component.parameters.filter)
    luminosity = component.parameters.luminosity

    # Set SED
    ski.set_stellar_component_sed(title, template, age, metallicity)

    # Set center wavelength of the filter as normalization wavelength (keeps label)
    ski.set_stellar_component_normalization_wavelength(title, fltr.center)

    # Set spectral luminosity at that wavelength (keeps label)
    ski.set_stellar_component_luminosity(title, luminosity)

    # Scale height doesn't need to be set as parameter, this is already in the deprojection model

# -----------------------------------------------------------------

def add_dust_component(ski, name, component):

    """
    This function ...
    :param ski: 
    :param name:
    :param component: 
    :return: 
    """

    # Debugging
    log.debug("Adding dust component '" + name + "' to the ski file ...")

    # If an input map is required
    if "map_path" in component: filename = set_dust_input_map(name, component)
    else: filename = None

    # Set geometry
    if "model" in component: set_dust_component_model(ski, component)

    # Set deprojection
    elif "deprojection" in component: set_dust_component_deprojection(ski, component)

    # From parameters
    if component.parameters is not None:

        # Check if this is a new dust component, add geometry, mix and normalization all at once
        if "geometry" in component.parameters: set_dust_component_geometry_mix_and_normalization(ski, component)

        # Existing component, THEMIS dust mix
        elif "hydrocarbon_pops" in component.parameters: set_dust_component_themis_mix(ski, component)

        # Existing component, not THEMIS dust mix
        else: raise NotImplementedError("Only THEMIS dust mixes are implemented at this moment")

        # TODO: implement 'Existing component, no THEMIX'
        #else: set_dust_component(ski, component)

    # From properties
    if component.properties is not None:

        # From unicode
        #if isinstance(name, unicode): name = name.encode('ascii','ignore')
        #print(name)
        # Create element
        #element = ski.create_element(element_name, component.properties)
        #print(element)

        ski.add_dust_component(component.properties, name)

    # Return the map filename
    return filename

# -----------------------------------------------------------------

def set_dust_input_map(name, component):

    """
    This function ...
    :param name:
    :param component:
    :return: 
    """

    # Generate a filename for the map
    map_filename = "dust_" + name + ".fits"

    # Set the filename
    if "deprojection" in component: component.deprojection.filename = map_filename
    #elif "geometry" in component.parameters: component.properties["geometry"].filename = map_filename # HOW IT WAS
    elif component.parameters is not None and "geometry" in component.parameters: component.parameters["geometry"].filename = map_filename
    elif component.properties is not None: component.properties["children"]["geometry"]["children"]["ReadFitsGeometry"]["filename"] = map_filename
    else: raise RuntimeError("Dust component based on an input map should either have a deprojection or geometry properties")

    # Add entry to the input maps dictionary
    #self.input_map_paths[filename] = component.map_path

    return map_filename

# -----------------------------------------------------------------

def set_dust_component_model(ski, component):

    """
    This function ...
    :param ski: 
    :param component: 
    :return: 
    """

    # Get title
    title = component.parameters.title

    # Set the geometry
    ski.set_dust_component_geometry(title, component.model)

# -----------------------------------------------------------------

def set_dust_component_deprojection(ski, component):

    """
    This function ...
    :param ski: 
    :param component: 
    :return: 
    """

    # Get title
    title = component.parameters.title

    # Set the deprojection geometry
    ski.set_dust_component_geometry(title, component.deprojection)

# -----------------------------------------------------------------

def set_dust_component_geometry_mix_and_normalization(ski, component):

    """
    This function ...
    :param ski: 
    :param component: 
    :return: 
    """

    # Get title
    title = component.parameters.title

    # Get class names
    geometry_type = component.parameters.geometry
    mix_type = component.parameters.sed
    normalization_type = component.parameters.normalization

    # Get properties for each of the three classes
    geometry_properties = component.properties["geometry"]
    mix_properties = component.properties["mix"]
    normalization_properties = component.properties["normalization"]

    # Create stellar component
    ski.create_new_dust_component(title, geometry_type, geometry_properties, mix_type, mix_properties, normalization_type, normalization_properties)

# -----------------------------------------------------------------

def set_dust_component_themis_mix(ski, component):

    """
    This function ...
    :param ski: 
    :param component: 
    :return: 
    """

    # Get title
    title = component.parameters.title

    # Get parameters
    mass = component.parameters.mass
    hydrocarbon_pops = component.parameters.hydrocarbon_pops
    enstatite_pops = component.parameters.enstatite_pops
    forsterite_pops = component.parameters.forsterite_pops

    # Set the dust mix
    ski.set_dust_component_themis_mix(title, hydrocarbon_pops, enstatite_pops, forsterite_pops)  # dust mix

    # Set the dust mass (keeps label)
    ski.set_dust_component_mass(title, mass)

# -----------------------------------------------------------------
