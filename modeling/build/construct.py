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
from ...core.tools.introspection import skirt_main_version, has_skirt
from ...core.tools.stringify import tostr
from ...core.filter.filter import Filter
from ...core.tools import types

# -----------------------------------------------------------------

# Check SKIRT version
if not has_skirt(): version_number = 8
else: version_number = skirt_main_version()

# Set flags
if version_number == 7:
    skirt7 = True
    skirt8 = False
elif version_number == 8:
    skirt7 = False
    skirt8 = True
else: raise RuntimeError("Invalid SKIRT version number")

# -----------------------------------------------------------------

def add_stellar_component(ski, name, component, title=None):

    """
    This function ...
    :param ski: 
    :param name:
    :param component:
    :param title:
    :return: 
    """

    # Debugging
    log.debug("Adding stellar component '" + name + "' to the ski file ...")

    # THIS HAS TO COME FIRST!!
    # If an input map is required
    if "map_path" in component: filename = set_stellar_input_map(name, component)
    else: filename = None

    # NEW COMPONENT OR ADJUST EXISTING
    if title is not None and not ski.has_stellar_component(title): add_new_stellar_component(ski, name, component, title=title)
    else: adjust_stellar_component(ski, name, component, title=title)

    # Return the input filename
    return filename

# -----------------------------------------------------------------

def add_new_stellar_component(ski, name, component, title=None):

    """
    This function ...
    :param ski:
    :param name:
    :param component:
    :param title:
    :return:
    """

    # Debugging
    log.debug("Adding new stellar component '" + name + "' to the ski file ...")

    # From properties
    if component.properties is not None:

        # Check title
        if title is None: log.warning("Title for the component '" + name + "' is not given")

        # Add component
        ski.add_stellar_component(component.properties, title=title)
        return

    # Initialize properties
    geometry = None
    geometry_type = None
    geometry_properties = None
    sed_type = None
    sed_properties = None
    normalization_type = None
    normalization_properties = None
    luminosities = [1]

    sed_template = None
    age = None
    metallicity = None
    compactness = None
    pressure = None
    covering_factor = None
    luminosity = None
    filter_or_wavelength = None

    # Set properties of the component
    if "model" in component: geometry = component.model
    elif "deprojection" in component: geometry = component.deprojection

    # Parameters are defined
    if component.parameters is not None:

        # Check if this is a new component (geometry not defined above): add geometry, SED and normalization all at once
        if "geometry" in component.parameters:

            # Get class names
            geometry_type = component.parameters.geometry
            sed_type = component.parameters.sed
            normalization_type = component.parameters.normalization

            # Get properties for each of the three classes
            geometry_properties = component.properties["geometry"]
            sed_properties = component.properties["sed"]
            normalization_properties = component.properties["normalization"]

        # Component with MAPPINGS template (geometry defined above)
        elif "sfr" in component.parameters: #set_stellar_component_mappings(ski, component)

            # Set template for MAPPINGS
            sed_template = "Mappings"

            # Get SED properties
            metallicity = component.parameters.metallicity
            compactness = component.parameters.compactness
            pressure = component.parameters.pressure
            covering_factor = component.parameters.covering_factor

            # Get normalization
            fltr = parse_filter(component.parameters.filter)
            luminosity = component.parameters.luminosity

            # Determine the normalization wavelength
            filter_or_wavelength = fltr.center

        # Existing component, no MAPPINGS
        else: # set_stellar_component(ski, component)

            # Get SED properties
            sed_template = component.parameters.template
            age = component.parameters.age
            metallicity = component.parameters.metallicity

            # Get normalization
            luminosity = component.parameters.luminosity

            # Determine the normalization wavelength
            if "wavelength" in component.parameters: wavelength = component.parameters.wavelength
            elif "filter" in component.parameters:
                fltr = parse_filter(component.parameters.filter)
                wavelength = fltr.wavelength
            else: raise ValueError("Neither wavelength nor filter is defined in the component parameters")

            # Set the normalization to the wavelength
            filter_or_wavelength = wavelength

    # Check whether title is defined
    if title is None: log.warning("Title for the component '" + name + "' is not defined")

    # Set normalization type
    if normalization_type is None:
        if filter_or_wavelength is None: raise ValueError("Cannot determine normalization type")
        if isinstance(filter_or_wavelength, Filter): normalization_type = "LuminosityStellarCompNormalization"
        elif types.is_length_quantity(filter_or_wavelength): normalization_type = "SpectralLuminosityStellarCompNormalization"
        else: normalization_type = "BolLuminosityStellarCompNormalization" #raise ValueError("Unrecognized filter of wavelength of type '" + str(type(filter_or_wavelength)))

    # Set stellar component properties
    properties = dict()
    properties["geometry"] = geometry
    properties["geometry_type"] = geometry_type
    properties["geometry_properties"] = geometry_properties
    properties["sed_type"] = sed_type
    properties["sed_properties"] = sed_properties
    properties["normalization_type"] = normalization_type
    properties["normalization_properties"] = normalization_properties
    properties["luminosities"] = luminosities
    properties["sed_template"] = sed_template
    properties["age"] = age
    properties["metallicity"] = metallicity
    properties["compactness"] = compactness
    properties["pressure"] = pressure
    properties["covering_factor"] = covering_factor
    properties["luminosity"] = luminosity
    properties["filter_or_wavelength"] = filter_or_wavelength

    # Show properties
    log.debug("")
    log.debug("Stellar component properties:")
    log.debug("")
    for label in properties:
        if label == "geometry":
            log.debug(" - geometry:")
            for parameter in properties[label]: log.debug("    * " + parameter + ": " + tostr(properties[label][parameter]))
        else: log.debug(" - " + label + ": " + tostr(properties[label]))
    log.debug("")

    # Create new component
    ski.create_new_stellar_component(title, **properties)

# -----------------------------------------------------------------

def adjust_stellar_component(ski, name, component, title=None):

    """
    This function ...
    :param ski:
    :param name:
    :param component:
    :param title:
    :return:
    """

    # Debugging
    log.debug("Adjusting existing stellar component '" + name + "' in the ski file ...")

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

        # Check if title is given
        if title is None: log.warning("Title for the component '" + name + "' is not specified")

        # Add component
        ski.add_stellar_component(component.properties, title=title)

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

    # Determine the normalization wavelength
    wavelength = fltr.center

    # Set SED
    ski.set_stellar_component_mappingssed(title, metallicity, compactness, pressure, covering_factor)  # SED

    # Set center wavelength of the filter as normalization wavelength (keeps label)
    ski.set_stellar_component_normalization_wavelength(title, wavelength)

    # Set spectral luminosity at that wavelength (keeps label)
    ski.set_stellar_component_luminosity(title, luminosity, filter_or_wavelength=wavelength)

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

    # Determine the normalization wavelength
    wavelength = fltr.center

    # Set SED
    ski.set_stellar_component_sed(title, template, age, metallicity)

    # Set center wavelength of the filter as normalization wavelength (keeps label)
    ski.set_stellar_component_normalization_wavelength(title, wavelength)

    # Set spectral luminosity at that wavelength (keeps label)
    ski.set_stellar_component_luminosity(title, luminosity, filter_or_wavelength=wavelength)

    # Scale height doesn't need to be set as parameter, this is already in the deprojection model

# -----------------------------------------------------------------

def add_dust_component(ski, name, component, title=None):

    """
    This function ...
    :param ski: 
    :param name:
    :param component:
    :param title:
    :return: 
    """

    # Debugging
    log.debug("Adding dust component '" + name + "' to the ski file ...")

    # THIS HAS TO COME FIRST!!
    # If an input map is required
    if "map_path" in component: filename = set_dust_input_map(name, component)
    else: filename = None

    # NEW COMPONENT OR ADJUST EXISTING
    if title is not None and not ski.has_dust_component(title): add_new_dust_component(ski, name, component, title=title)
    else: adjust_dust_component(ski, name, component, title=title)

    # Return the map filename
    return filename

# -----------------------------------------------------------------

def add_new_dust_component(ski, name, component, title=None):

    """
    This function ...
    :param ski:
    :param name:
    :param component:
    :param title:
    :return:
    """

    # Debugging
    log.debug("Adding new dust component '" + name + "' to the ski file ...")

    # From properties
    if component.properties is not None:

        # Check if title is given
        if title is None: log.warning("Title of the component '" + name + "' is not given")

        # Add component
        ski.add_dust_component(component.properties, title=title)
        return

    # Initialize properties
    geometry = None
    geometry_type = None
    geometry_properties = None
    mix_type = None
    mix_properties = None
    normalization_type = None
    normalization_properties = None

    mix = None
    mass = None

    # For THEMIS mix
    hydrocarbon_pops = None
    silicate_pops = None

    # For Zubko mix
    graphite_populations = None
    silicate_populations = None
    pah_populations = None

    # Set properties of the component
    if "model" in component: geometry = component.model
    elif "deprojection" in component: geometry = component.deprojection

    # Parameters are defined
    if component.parameters is not None:

        # Check title
        if title is not None and component.parameters.title != title: raise ValueError("The title of the component '" + title + "' doesn't match that defined in the component parameters")

        # Check if this is a new component (geometry not defined above): add geometry, mix and normalization
        if "geometry" in component.parameters:

            # Get class names
            geometry_type = component.parameters.geometry
            mix_type = component.parameters.sed
            normalization_type = component.parameters.normalization

            # Get properties for each of the three classes
            geometry_properties = component.properties["geometry"]
            mix_properties = component.properties["mix"]
            normalization_properties = component.properties["normalization"]

        # Existing component (geometry defined above), THEMIS dust mix
        elif "hydrocarbon_pops" in component.parameters: #set_dust_component_themis_mix(ski, component)

            # Set mix name
            mix = "themis"

            # Get parameters
            mass = component.parameters.mass
            hydrocarbon_pops = component.parameters.hydrocarbon_pops
            silicate_pops = component.parameters.silicate_pops

        # Existing component (geometry defined above), Zubko dust mix
        elif "graphite_populations" in component.parameters:

            # Set mix name
            mix = "zubko"

            # Get parameters
            mass = component.parameters.mass
            graphite_populations = component.parameters.graphite_populations
            silicate_populations = component.parameters.silicate_populations
            pah_populations = component.parameters.pah_populations

        # Existing component, not THEMIS dust mix
        else: raise NotImplementedError("Only THEMIS dust mixes are implemented at this moment")

    # Check whether the title is defined
    if title is None: log.warning("The title for the '" + name + "' dust component is not specified")

    # Set dust component properties
    properties = dict()
    properties["geometry"] = geometry
    properties["geometry_type"] = geometry_type
    properties["geometry_properties"] = geometry_properties
    properties["mix_type"] = mix_type
    properties["mix_properties"] = mix_properties
    properties["normalization_type"] = normalization_type
    properties["normalization_properties"] = normalization_properties
    properties["mix"] = mix
    properties["mass"] = mass
    properties["hydrocarbon_pops"] = hydrocarbon_pops
    properties["silicate_pops"] = silicate_pops
    properties["graphite_populations"] = graphite_populations
    properties["silicate_populations"] = silicate_populations
    properties["pah_populations"] = pah_populations

    # Show properties
    log.debug("")
    log.debug("Dust component properties:")
    log.debug("")
    for label in properties: log.debug(" - " + label + ": " + tostr(properties[label]))
    log.debug("")

    # Add the new component
    ski.create_new_dust_component(title, **properties)

# -----------------------------------------------------------------

def adjust_dust_component(ski, name, component, title=None):

    """
    Thisf unction ...
    :param ski:
    :param name:
    :param component:
    :param title:
    :return:
    """

    # Debugging
    log.debug("Adjusting existing dust component '" + name + "' in the ski file ...")

    # Set geometry
    if "model" in component: set_dust_component_model(ski, component)

    # Set deprojection
    elif "deprojection" in component: set_dust_component_deprojection(ski, component)

    # From parameters
    if component.parameters is not None:

        # Check title
        if title is not None and component.parameters.title != title: raise ValueError("The title of the component '" + title + "' doesn't match that defined in the component parameters")

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

        # Check if title is given
        if title is None: log.warning("Title for the component '" + name + "' is not given")
        ski.add_dust_component(component.properties, title=title)

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
    silicate_pops = component.parameters.silicate_pops

    # Set the dust mix
    ski.set_dust_component_themis_mix(title, hydrocarbon_pops, silicate_pops)  # dust mix

    # Set the dust mass (keeps label)
    ski.set_dust_component_mass(title, mass)

# -----------------------------------------------------------------
