#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.upgradeskifile Upgrade a ski file to the latest format version

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import subprocess
from lxml import etree

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ..tools.logging import log
from ..tools import introspection
from ..tools import filesystem as fs
from ..basics.map import Map
from ..tools import xml
from ..basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter
from ..basics.unit import parse_unit as u
from ..tools import parsing
from ..tools import formatting as fmt
from ..basics.configuration import print_mapping
from ..tools import stringify

# -----------------------------------------------------------------

# bulkmassdensity: g / cm3,kg / m3
# distance: AU,Mpc,cm,km,kpc,m,pc
# grainsize: A,cm,m,micron,mm,nm
# length: AU,Mpc,cm,km,kpc,m,pc
# mass: Msun,g,kg
# massvolumedensity: Msun / AU3,Msun / pc3,g / cm3,kg / m3
# monluminosity: Lsun / micron,W / m,W / micron
# opacity: m2 / kg
# posangle: deg,rad
# pressure: K / m3,Pa
# temperature: K
# wavelength: A,cm,m,micron,mm,nm

skirt_quantities_to_pts_quantities = dict()
skirt_quantities_to_pts_quantities["bulkmassdensity"] = "mass_density_quantity"
skirt_quantities_to_pts_quantities["distance"] = "length_quantity"
skirt_quantities_to_pts_quantities["grainsize"] = "length_quantity"
skirt_quantities_to_pts_quantities["length"] = "length_quantity"
skirt_quantities_to_pts_quantities["mass"] = "mass_quantity"
skirt_quantities_to_pts_quantities["massvolumedensity"] = "mass_density_quantity"
skirt_quantities_to_pts_quantities["monluminosity"] = "photometric_density_quantity"
skirt_quantities_to_pts_quantities["opacity"] = "quantity"
skirt_quantities_to_pts_quantities["posangle"] = "angle"
skirt_quantities_to_pts_quantities["pressure"] = "quantity"
skirt_quantities_to_pts_quantities["temperature"] = "temperature_quantity"
skirt_quantities_to_pts_quantities["wavelength"] = "length_quantity"

# -----------------------------------------------------------------

class SKIRTSmileSchema(object):

    """
    This function ...
    """

    def __init__(self):

        """
        This function ...
        """

        # Determine path to PTS temporary directory
        temp_path = introspection.pts_temp_dir

        # Create the command
        command = ["skirt", "-x"]

        # Run SKIRT
        if log.is_debug(): subprocess.call(command, cwd=temp_path)
        else: subprocess.call(command, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'), cwd=temp_path)

        # Load the smile scheme
        smile_path = fs.join(temp_path, "skirt.smile")

        # load the XML tree (remove blank text to avoid confusing the pretty printer when saving)
        self.tree = etree.parse(smile_path, parser=etree.XMLParser(remove_blank_text=True))

    # -----------------------------------------------------------------

    def get_types(self):

        """
        This function ...
        :return:
        """

        return self.tree.xpath("//types")[0].getchildren()

    # -----------------------------------------------------------------

    def get_concrete_types(self):

        """
        This function ...
        :return:
        """

        return [t for t in self.get_types() if "concrete" in t.attrib and t.attrib["concrete"] == "true"]

    # -----------------------------------------------------------------

    @lazyproperty
    def concrete_types(self):

        """
        This function ...
        :return:
        """

        types = dict()

        # Loop
        for t in self.get_concrete_types():

            name = t.attrib["name"]
            description = t.attrib["title"]

            # Add to dict
            types[name] = description

        # Return the dictionary
        return types

    # -----------------------------------------------------------------

    def is_derived_from(self, name, base_name):

        """
        This function ...
        :param name:
        :param base_name:
        :return:
        """

        t = self.find_type(name)

        # Get base
        if "base" not in t.attrib: return False
        base = t.get("base")
        if base == "": return False

        if base == base_name: return True
        else: return self.is_derived_from(base, base_name)

    # -----------------------------------------------------------------

    def get_concrete_types_with_base(self, base):

        """
        This function ...
        :param base:
        :return:
        """

        types = []

        # Loop
        for t in self.get_concrete_types():

            name = t.attrib["name"]

            # If this concrete type is the base itsself
            if name == base: types.append(t)

            # Derived from
            elif self.is_derived_from(name, base): types.append(t)

        # Return the types
        return types

    # -----------------------------------------------------------------

    def concrete_types_with_base(self, base):

        """
        This function ...
        :param base:
        :return:
        """

        types = dict()

        # Loop
        for t in self.get_concrete_types():

            name = t.attrib["name"]
            description = t.attrib["title"]

            # If this concrete type is the base itsself
            if name == base: #return {name: description}
                types[name] = description

            # Check if derived from specified base
            # If derived, add to dict
            elif self.is_derived_from(name, base): types[name] = description

        # Return the dictionary
        return types

    # -----------------------------------------------------------------

    def get_concrete_geometries(self):

        """
        This function ...
        :return:
        """

        return [geometry for geometry in self.get_concrete_types() if geometry.attrib["name"].endswith("Geometry")]

    # -----------------------------------------------------------------

    @lazyproperty
    def concrete_geometries(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        geometries = dict()

        # Loop
        for geometry in self.get_concrete_geometries():

            name = geometry.attrib["name"]
            description = geometry.attrib["title"]

            # Add to the dictionary
            geometries[name] = description

        # Return the dictionary
        return geometries

    # -----------------------------------------------------------------

    def get_concrete_stellar_normalizations(self):

        """
        This function ...
        :return:
        """

        return [normalization for normalization in self.get_concrete_types() if normalization.attrib["name"].endswith("StellarCompNormalization")]

    # -----------------------------------------------------------------

    @lazyproperty
    def concrete_stellar_normalizations(self):

        """
        This function ...
        :return:
        """

        # Initialize a dictionary
        normalizations = dict()

        # Loop
        for normalization in self.get_concrete_stellar_normalizations():

            name = normalization.attrib["name"]
            description = normalization.attrib["title"]

            # Add to the dictionary
            normalizations[name] = description

        # Return the dictionary
        return normalizations

    # -----------------------------------------------------------------

    def get_concrete_dust_normalizations(self):

        """
        This function ...
        :return:
        """

        return [normalization for normalization in self.get_concrete_types() if normalization.attrib["name"].endswith("DustCompNormalization")]

    # -----------------------------------------------------------------

    @lazyproperty
    def concrete_dust_normalizations(self):

        """
        This function ...
        :return:
        """

        # Initialize a dictionary
        normalizations = dict()

        # Loop
        for normalization in self.get_concrete_dust_normalizations():

            name = normalization.attrib["name"]
            description = normalization.attrib["title"]

            # Add to the dictionary
            normalizations[name] = description

        # Return the dictionary
        return normalizations

    # -----------------------------------------------------------------

    def get_concrete_stellar_seds(self):

        """
        This function ...
        :return:
        """

        return self.get_concrete_types_with_base("StellarSED")

    # -----------------------------------------------------------------

    @lazyproperty
    def concrete_stellar_seds(self):

        """
        This function ...
        :return:
        """

        seds = dict()

        # Loop
        for sed in self.get_concrete_stellar_seds():

            name = sed.attrib["name"]
            description = sed.attrib["title"]

            # Add to the dictionary
            seds[name] = description

        # Return the dictionary
        return seds

    # -----------------------------------------------------------------

    def get_concrete_dust_mixes(self):

        """
        This function ...
        :return:
        """

        return [mix for mix in self.get_concrete_types() if mix.attrib["name"].endswith("DustMix")]

    # -----------------------------------------------------------------

    @lazyproperty
    def concrete_dust_mixes(self):

        """
        This function ...
        :return:
        """

        # Initialize a dictionary
        mixes = dict()

        # Loop
        for mix in self.get_concrete_dust_mixes():

            name = mix.attrib["name"]
            description = mix.attrib["title"]

            # Add to the dictionary
            mixes[name] = description

        # Return the dictionary
        return mixes

    # -----------------------------------------------------------------

    def find_type(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Search for type
        for t in self.get_types():
            if t.get("name") == name: return t

        raise RuntimeError("Could not find '" + name + "' as a type in the smile schema")

    # -----------------------------------------------------------------

    def get_properties_for_type(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Find the type in the tree
        t = self.find_type(name)

        # Get list of properties
        properties = xml.get_list(t, "properties")

        # Look for derived properties
        if "base" in t.attrib and t.attrib["base"] != "":
            base_name = t.get("base")
            properties += self.get_properties_for_type(base_name)

        # Return the properties
        return properties

    # -----------------------------------------------------------------

    def properties_for_type(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        properties = dict()

        for property in self.get_properties_for_type(name):

            ptype = get_ptype(property.tag)
            pname = property.get("name")
            description = property.get("title").encode('utf-8').strip()
            default = property.get("default") if "default" in property.attrib else None

            # Enumeration
            if ptype == "ENUM":

                values = xml.get_list(property, "enumValues")
                enum = dict()
                for value in values:
                    value_name = value.get("name")
                    value_description = value.get("title")
                    enum[value_name] = value_description

                if len(enum) == 0: print(name, pname)

                # Add property
                properties[pname] = Map(ptype="string", description=description, default=default, choices=enum)

            # Enumeration list
            elif ptype == "ENUM_LIST":

                values = xml.get_list(property, "enumValues")
                enum = dict()
                for value in values:
                    value_name = value.get("name")
                    value_description = value.get("title")
                    enum[value_name] = value_description

                if len(enum) == 0: print(name, pname)

                # Add property
                properties[pname] = Map(ptype="string_list", description=description, choices=enum)

            # Simulation item
            elif ptype == "ITEM":

                # Get allowed concrete types
                base = property.get("base")
                choices = self.concrete_types_with_base(base)

                if len(choices) == 0: print(name, pname)

                # Add property
                properties[pname] = Map(ptype="string", description=description, default=default, choices=choices, item=True)

            # Simulation item list
            elif ptype == "ITEM_LIST":

                # Get allowed concrete types
                base = property.get("base")
                choices = self.concrete_types_with_base(base)

                if len(choices) == 0: print(name, pname)

                # Add property
                properties[pname] = Map(ptype="string_list", description=description, choices=choices, item=True)

            # Other
            else:

                quantity = property.get("quantity") if "quantity" in property.attrib else None
                if quantity == "":
                    ptype = "real"
                    quantity = None

                if quantity is not None: ptype = skirt_quantities_to_pts_quantities[quantity]

                # Parse
                parser = getattr(parsing, ptype)
                default = parser(default) if default is not None else None

                min_value = property.get("min") if "min" in property.attrib else None
                max_value = property.get("max") if "max" in property.attrib else None

                min_value = parser(min_value) if min_value is not None else None
                max_value = parser(max_value) if max_value is not None else None

                # Add property
                properties[pname] = Map(ptype=ptype, description=description, min=min_value, max=max_value, default=default)

        # Return the properties dictionary
        return properties

    # -----------------------------------------------------------------

    def definition_for_type(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get properties
        properties = self.properties_for_type(name)

        # Keep track of the properties that are actually simulation items
        simulation_items = []

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Show the properties
        for prop_name in properties:

            description = properties[prop_name].description
            ptype = properties[prop_name].ptype
            min_value = properties[prop_name].min
            max_value = properties[prop_name].max
            default = properties[prop_name].default
            choices = properties[prop_name].choices
            item = properties[prop_name].item

            # Add setting
            if ptype == "boolean": definition.add_flag(prop_name, description, default)
            elif default is None: definition.add_required(prop_name, ptype, description, choices=choices, min_value=min_value, max_value=max_value, dynamic_list=True)
            else: definition.add_optional(prop_name, ptype, description, default=default, choices=choices, min_value=min_value, max_value=max_value, dynamic_list=True)

            # Add to list if simulation item
            if item: simulation_items.append(prop_name)

        # Return the definition
        return definition, simulation_items

    # -----------------------------------------------------------------

    def prompt_parameters_for_type(self, name, merge=False):

        """
        This function ...
        :param name:
        :param merge:
        :return:
        """

        # Get the configuration definition
        definition, simulation_items = self.definition_for_type(name)

        # Prompt for parameters
        setter = InteractiveConfigurationSetter("configuration of " + name + " simulation item", add_cwd=False, add_logging=False)
        config = setter.run(definition, prompt_optional=True)

        children = dict()

        # Prompt for parameters of simulation item properties
        for prop_name in simulation_items:

            # List of simulation items
            if isinstance(config[prop_name], list):

                item_names = config[prop_name]

                for item_name in item_names:
                    child_parameters, child_children = self.prompt_parameters_for_type(item_name)
                    children[item_name] = child_parameters, child_children

            # Single simulation item
            else:

                item_name = config[prop_name]
                child_parameters, child_children = self.prompt_parameters_for_type(item_name)
                children[item_name] = child_parameters, child_children

        # Return the parameters of this item and of its children
        if merge: return merge_parameters(config, children)
        else: return config, children

    # -----------------------------------------------------------------

    def get_quantities(self):

        """
        This function ...
        :return:
        """

        return self.tree.xpath("//quantities")[0].getchildren()

    # -----------------------------------------------------------------

    @lazyproperty
    def quantities(self):

        """
        This function ...
        :return:
        """

        return [quantity.get("name") for quantity in self.get_quantities()]

    # -----------------------------------------------------------------

    def find_quantity(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Search for quantity
        for t in self.get_quantities():
            if t.get("name") == name: return t

        raise RuntimeError("Could not find '" + name + "' as a quantity in the smile schema")

    # -----------------------------------------------------------------

    def units_for_quantity(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        units = xml.get_list(self.find_quantity(name), "units")
        return [u(unit.get("name")) for unit in units]

    # -----------------------------------------------------------------

    def get_unit_systems(self):

        """
        This function ...
        :return:
        """

        return self.tree.xpath("//unitSystems")[0].getchildren()

    # -----------------------------------------------------------------

    @lazyproperty
    def unit_systems(self):

        """
        This function ...
        :return:
        """

        return [system.get("name") for system in self.get_unit_systems()]

    # -----------------------------------------------------------------

    def find_unit_system(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Search for unit system
        for t in self.get_unit_systems():
            if t.get("name") == name: return t

        raise RuntimeError("Could not find '" + name + "' as a unit system in the smile schema")

    # -----------------------------------------------------------------

    def units_for_unit_system(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        units = dict()

        defaults = xml.get_list(self.find_unit_system(name), "defaultUnits")

        for default in defaults:

            quantity_name = default.get("quantity")
            unit = u(default.get("unit"))

            # Add to units dictionary
            units[quantity_name] = unit

        return units

# -----------------------------------------------------------------

def get_ptype(tag):

    """
    This function ...
    :param tag:
    :return:
    """

    if tag == "DoubleProperty": return "real"
    elif tag == "IntProperty": return "integer"
    elif tag == "BoolProperty": return "boolean"
    elif tag == "StringProperty": return "string"
    elif tag == "DoubleListProperty": return "real_list"
    elif tag == "IntListProperty": return "integer_list"

    # SPECIAL
    elif tag == "EnumProperty": return "ENUM"
    elif tag == "EnumListProperty": return "ENUM_LIST"
    elif tag == "ItemProperty": return "ITEM"
    elif tag == "ItemListProperty": return "ITEM_LIST"

    # Not recognized
    else: raise RuntimeError("Did not recognize tag '" + tag + "' for property type")

# -----------------------------------------------------------------

def show_parameters(parameters, children, indent="", name=None):

    """
    This function ...
    :param parameters:
    :param children:
    :param indent:
    :param name:
    :return:
    """

    print("")
    if len(parameters) == 0:
        if name is not None: print(indent + fmt.red + fmt.underlined + name + ": no parameters" + fmt.reset)
        else: print(indent + fmt.red + fmt.underlined + "No parameters" + fmt.reset)
    else:
        if name is not None: print(indent + fmt.red + fmt.underlined + name + " parameters" + fmt.reset)
        else: print(indent + fmt.red + fmt.underlined + "Parameters" + fmt.reset)

    # Print mapping
    print_mapping(parameters, indent=indent)

    # Show simulation items
    if len(children) > 0:
        print(indent + fmt.blue + "Simulation items: " + fmt.reset + stringify.stringify(children.keys())[1].replace(",", ", "))
    print("")

    # Print children
    for name in children:

        # Get parameters of child (and children)
        parameters, child_children = children[name]

        # Show
        show_parameters(parameters, child_children, indent=indent+"    ", name=name)

# -----------------------------------------------------------------

def merge_parameters(parameters, children):

    """
    This function ...
    :param parameters:
    :param children:
    :return:
    """

    new = dict()

    children_names = children.keys()

    for name in parameters:

        if name.startswith("_"): continue

        value = parameters[name]

        if isinstance(value, list):

            new[name] = dict()

            for val in value:

                if val in children_names:

                    child_parameters, child_children = children[val]
                    new[name][val] = merge_parameters(child_parameters, child_children)

        elif value in children_names:

            child_parameters, child_children = children[value]
            new[name] = (value, merge_parameters(child_parameters, child_children))

        else: new[name] = parameters[name]

    # Return the new parameters dict
    return new

# -----------------------------------------------------------------
