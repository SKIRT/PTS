#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.composite Contains the SimplePropertyComposite class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
from collections import OrderedDict
import warnings
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ..tools import parsing
from ..basics.log import log
from ..tools import formatting as fmt
from ..tools.stringify import stringify, tostr
from .map import Map

# -----------------------------------------------------------------

class SimplePropertyComposite(object):

    """
    This class ...
    @DynamicAttrs
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self):

        """
        This function ...
        """

        # The path
        self._path = None

        # The descriptions
        self._descriptions = dict()

        # The parsing types
        self._ptypes = dict()

        # The choices
        self._choices = dict()

    # -----------------------------------------------------------------

    def add_property(self, name, ptype, description, default_value=None, choices=None, convert_default=False):

        """
        This function ...
        :param name:
        :param ptype:
        :param description:
        :param default_value:
        :param choices:
        :param convert_default:
        :return:
        """

        # Check
        if hasattr(self, name): raise ValueError("A property with the name '" + name + "' already exists")

        # Check
        if " " in name: raise ValueError("Name cannot contain spaces")

        # Set the ptype
        self._ptypes[name] = ptype

        # Set the description
        self._descriptions[name] = description

        # Set the choices
        self._choices[name] = choices

        # Convert default
        if default_value is not None and convert_default:
            parsing_function = getattr(parsing, ptype)
            default_value = parsing_function(default_value)

        # Set the attribute with the default value
        setattr(self, name, default_value)

    # -----------------------------------------------------------------

    def add_integer_property(self, name, description, default_value=None, choices=None):

        """
        This function ...
        :param name:
        :param description:
        :param default_value:
        :param choices:
        :return:
        """

        self.add_property(name, "integer", description, default_value=default_value, choices=choices)

    # -----------------------------------------------------------------

    def add_real_property(self, name, description, default_value=None, choices=None):

        """
        This function ...
        :param name:
        :param description:
        :param default_value:
        :param choices:
        :return:
        """

        self.add_property(name, "real", description, default_value=default_value, choices=choices)

    # -----------------------------------------------------------------

    def add_string_property(self, name, description, default_value=None, choices=None):

        """
        This function ...
        :param name:
        :param description:
        :param default_value:
        :param choices:
        :return:
        """

        self.add_property(name, "string", description, default_value=default_value, choices=choices)

    # -----------------------------------------------------------------

    def add_boolean_property(self, name, description, default_value=None, choices=None):

        """
        This function ...
        :param name:
        :param description:
        :param default_value:
        :param choices:
        :return:
        """

        self.add_property(name, "boolean", description, default_value=default_value, choices=choices)

    # -----------------------------------------------------------------

    def add_section(self, name, description, dynamic=False):

        """
        This function ...
        :param name:
        :param description:
        :param dynamic:
        :return:
        """

        # Set the description
        self._descriptions[name] = description

        # Set an attribute that is a nested SimplePropertyComposite (or a Map)
        if dynamic: self.__dict__[name] = Map() #setattr(self, name, Map())
        else: self.__dict__[name] = SimplePropertyComposite() #setattr(self, name, SimplePropertyComposite())

    # -----------------------------------------------------------------

    @property
    def sections(self):

        """
        This function ...
        :return:
        """

        sections_dict = dict()
        for name in self.section_names:
            sections_dict[name] = getattr(self, name)
        return sections_dict

    # -----------------------------------------------------------------

    def __setattr__(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        # Hidden variable
        if name.startswith("_"):
            self.__dict__[name] = value
            return

        # Is property setter
        if hasattr(self.__class__, name) and isinstance(getattr(self.__class__, name), property):

            # IS PROPERTY
            class_attribute = getattr(self.__class__, name)
            assert isinstance(class_attribute, property)
            # class_attribute.fget # this is the properties 'get' entry
            if class_attribute.fset is None: raise ValueError("The property '" + name + "' does not have a setter")

            # Set the property
            class_attribute.fset(self, value)
            return

        # Set 'simple' property
        if value is None: pass
        elif isinstance(value, SimplePropertyComposite): assert name in self._descriptions
        elif isinstance(value, dict):  # elif isinstance(value, Map):
            assert name in self._descriptions
            #print(value)
            #print(getattr(self, name))
            #print(name in self.__dict__)
            for key in value:
                keyvalue = value[key]
                #print(self.__dict__)
                self.__dict__[name].__setattr__(key, keyvalue)
            return
        else:

            # Check the type
            #print(value, type(value), hasattr(value, "__array__"))
            ptype, string = stringify(value)

            # None value
            if string == "None": value = None

            # Actual value
            else:

                # Check
                if name not in self._ptypes: raise AttributeError("A " + self.__class__.__name__ + " object has no attribute '" + name + "'")

                # Try converting the string back to the type it actually needs to be
                the_type = self._ptypes[name]
                parsing_function = getattr(parsing, the_type)
                try: value = parsing_function(string)
                except ValueError: raise ValueError("The value of '" + str(value) + "' for '" + name +  "' given is of the wrong type: '" + ptype + "', must be '" + the_type + "' (value is " + string + ")")

        # Set the attribute
        self.__dict__[name] = value

    # -----------------------------------------------------------------

    def set_property(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        # Check
        if name not in self.property_names: raise ValueError("Not a property of this class: '" + name + "'")

        # Set attribute
        setattr(self, name, value)

    # -----------------------------------------------------------------

    def set_properties(self, properties):

        """
        This function ...
        :param properties:
        :return:
        """

        # Loop over all the options defined in the 'options' dictionary
        for name in properties:

            # Check whether an option with this name exists in this class
            if hasattr(self, name):

                # Check if the option is composed of other options (a Map), or if it is just a simple variable
                #if isinstance(getattr(self, name), Map):
                if isinstance(getattr(self, name), SimplePropertyComposite):

                    assert isinstance(properties[name], dict)  # dict, or Map (is derived from dict)
                    getattr(self, name).set_properties(properties[name])

                # If it is a simple variable, just use setattr to set the attribute of this class
                else: setattr(self, name, properties[name])

            # If the option does not exist, ignore it but give a warning
            else: warnings.warn("The property '" + name + "' does not exist")

    # -----------------------------------------------------------------

    def prompt_properties(self, recursive=True, contains=None, not_contains=None, exact_name=None, exact_not_name=None):

        """
        This function ...
        :param recursive:
        :param contains:
        :param not_contains:
        :param exact_name:
        :param exact_not_name:
        :return:
        """

        from .configuration import prompt_variable

        # Loop over the properties
        for name in self.property_names:

            # Checks
            if contains is not None and contains not in name: continue
            if not_contains is not None and not_contains in name: continue
            if exact_name is not None and name != exact_name: continue
            if exact_not_name is not None and name == exact_not_name: continue

            description = self.get_description(name)
            ptype = self.get_ptype(name)
            default = self.get_value(name)
            choices = self.get_choices(name)

            # Ask for the new value
            value = prompt_variable(name, ptype, description, choices=choices, default=default, required=False)

            # Set the new value
            self.set_property(name, value)

        # Recursive: also loop over the settings
        if recursive:

            # Loop over the sections
            for name in self.section_names:

                # Prompt the settings of the section
                self.sections[name].prompt_properties(recursive=True, contains=contains, not_contains=not_contains, exact_name=exact_name, exact_not_name=exact_not_name)

    # -----------------------------------------------------------------

    @classmethod
    def from_dict(cls, dictionary):

        """
        This function ...
        :param dictionary
        :return:
        """

        # Create a new instance
        composite = cls()

        properties = OrderedDict() # to keep order of original dict (possibly an ordered dict as well)
        for key in dictionary:
            ptype, pvalue = stringify(dictionary[key])
            description = key
            name = key.replace(" ", "_").lower()
            properties[name] = dictionary[key] # set the value
            composite.add_property(name, ptype, description)

        # Set the properties from the dictionary
        composite.set_properties(properties)

        # Return the new property composite
        return composite

    # -----------------------------------------------------------------

    @property
    def all_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for name in vars(self):

            # Skip internal variables
            if name.startswith("_"): continue

            # Add the name
            names.append(name)

        # Return the names
        return names

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        This function ...
        :return:
        """

        for name in self.property_names: yield name

    # -----------------------------------------------------------------

    def get_value(self, name):

        """
        This funtion ...
        :param name:
        :return:
        """

        return getattr(self, name)

    # -----------------------------------------------------------------

    def __getitem__(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_value(name)

    # -----------------------------------------------------------------

    def get_ptype(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self._ptypes[name]

    # -----------------------------------------------------------------

    def get_description(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self._descriptions[name]

    # -----------------------------------------------------------------

    def get_choices(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return self._choices[name]

    # -----------------------------------------------------------------

    @property
    def property_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for name in vars(self):

            # Skip internal variables
            if name.startswith("_"): continue

            if name not in self._ptypes:
                if not (isinstance(getattr(self, name), SimplePropertyComposite) or isinstance(getattr(self, name), Map)): raise Exception("Property '" + name + "' doesn't have its type defined")
                else: continue

            # Add the name
            names.append(name)

        # Return the names
        return names

    # -----------------------------------------------------------------

    @property
    def section_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for name in vars(self):

            # Skip internal
            if name.startswith("_"): continue

            # Skip simple properties
            if name in self._ptypes: continue

            # Should be composite
            if not (isinstance(getattr(self, name), SimplePropertyComposite) or isinstance(getattr(self, name), Map)): raise Exception("Property '" + name + "' doesn't have its type defined")

            # Add the name
            names.append(name)

        # Return the names
        return names

    # -----------------------------------------------------------------

    @property
    def property_types(self):

        """
        This function ...
        :return:
        """

        types = []
        for name in self.property_names:

            # Get the type
            ptype = self._ptypes[name]
            types.append(ptype)

        # Return the types
        return types

    # -----------------------------------------------------------------

    def type_for_property(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self._ptypes[name]

    # -----------------------------------------------------------------

    @property
    def property_units(self):

        """
        This function ...
        :return:
        """

        units = []
        for name in self.property_names:

            # Get the value
            value = getattr(self, name)

            # Check whether there is a unit
            if hasattr(value, "unit"): unit = value.unit
            else: unit = None

            # Add the unit
            units.append(unit)

        # Return the list of units
        return units

    # -----------------------------------------------------------------

    def unit_for_property(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the value
        value = getattr(self, name)

        # Check whether there is a unit
        if hasattr(value, "unit"): unit = value.unit
        else: unit = None

        # Return the unit
        return unit

    # -----------------------------------------------------------------

    @property
    def property_descriptions(self):

        """
        This function ...
        :return:
        """

        descriptions = []
        for name in self.property_names:

            # Get the description
            description = self._descriptions[name]
            descriptions.append(description)

        # Return the descriptions
        return descriptions

    # -----------------------------------------------------------------

    def description_for_property(self, name):

        """
        This function ...
        """

        return self._descriptions[name]

    # -----------------------------------------------------------------

    def get_properties(self):

        """
        This function ...
        :return:
        """

        properties = dict()
        for name in vars(self):

            # Skip internal variables
            if name.startswith("_"): continue

            # Set property
            properties[name] = getattr(self, name)

        # Return the properties as a dictionary
        return properties

    # -----------------------------------------------------------------

    def to_dict(self):

        """
        This function ...
        :return:
        """

        # Create dict with the properties
        result = self.get_properties()

        # Add the sections (as dicts as well)
        for name in self.section_names: result[name] = self.sections[name].to_dict()

        # Return the resulting dict
        return result

    # -----------------------------------------------------------------

    def to_lines(self, line_prefix="", ignore=None, ignore_none=False):

        """
        This function ...
        :param line_prefix:
        :param ignore:
        :param ignore_none:
        :return:
        """

        lines = []

        # The simple properties
        for name in self.property_names:

            # Get the value
            value = getattr(self, name)

            # Check ignore
            if ignore is not None and name in ignore: continue
            if ignore_none and value is None: continue

            # Generate line
            dtype, value = stringify(value)
            line = " - " + fmt.bold + name + fmt.reset + ": " + value
            lines.append(line_prefix + line)

        # Sections
        for name in self.section_names:

            # Check ignore
            if ignore is not None and name in ignore: continue

            line = " - " + fmt.bold + name + fmt.reset + ":"
            lines.append(line_prefix + line)

            if isinstance(getattr(self, name), SimplePropertyComposite): section_lines = [line for line in getattr(self, name).to_lines(line_prefix=line_prefix+"    ")]
            elif isinstance(getattr(self, name), Map):

                section_lines = []
                section = getattr(self, name)
                for key in section:

                    line = line_prefix + "    " + " - " + fmt.bold + key + fmt.reset + ": " + tostr(section[key])
                    section_lines.append(line)

            else: raise ValueError("Unknown type for section: " + str(type(getattr(self, name))))

            lines += section_lines

        # Return the lines
        return lines

    # -----------------------------------------------------------------

    def to_string(self, line_prefix="", ignore=None, ignore_none=False):

        """
        This function ...
        :param line_prefix:
        :param ignore:
        :param ignore_none:
        :return:
        """

        # Return
        return "\n".join(self.to_lines(line_prefix=line_prefix, ignore=ignore, ignore_none=ignore_none))

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        return self.to_string()

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, remote=None):

        """
        This function ...
        :param path:
        :param remote:
        :return:
        """

        # FROM REMOTE FILE
        if remote is not None: return cls.from_remote_file(path, remote)

        # Inform the user
        log.debug("Loading " + cls.__name__ + " from '" + path + "' ...")

        # Initialize dictionary to contain the properties
        properties = dict()

        # Open the file
        with open(path, 'r') as f:

            # Read first line
            for line in f:
                first = line
                break

            #composite_type = first.split("Type: ")[1]
            #assert composite_type == cls.__name__

            # Load
            load_properties(properties, f)

        #print(properties)

        # Create the class instance
        composite = cls(**properties)

        # Set the path
        composite._path = path

        # Return
        return composite

    # -----------------------------------------------------------------

    @classmethod
    def from_remote_file(cls, path, remote):

        """
        This function ...
        :param path:
        :param remote:
        :return:
        """

        # Inform the user
        log.debug("Loading " + cls.__name__ + " from '" + path + "' on remote host '" + remote.host_id + "' ...")

        # Initialize dictionary to contain the properties
        properties = dict()

        # Create iterator for the lines
        lines_iterator = remote.read_lines(path)

        # Read first line (skip it)
        for line in lines_iterator:
            first = line
            break

        # Load the properties
        load_properties(properties, lines_iterator)

        # Create the class instance
        composite = cls(**properties)

        # Set the path: NO
        #composite._path = path

        # Return
        return composite

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check whether the path is valid
        if self._path is None: raise RuntimeError("The path is not defined")

        # Save
        self.saveto(self._path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Saving the " + self.__class__.__name__ + " to " + path + " ...")

        # Write the properties
        with open(path, 'w') as fh:

            # Print the type
            print("Type:", self.__class__.__name__, file=fh)

            # Write
            write_properties(self, fh)

        # Update the path
        self._path = path

    # -----------------------------------------------------------------

    def as_tuples(self):

        """
        This function ...
        :return:
        """

        tuples = []

        for name in self.property_names:
            value = getattr(self, name)
            tuples.append((name, value))

        for name in self.section_names:
            value = getattr(self, name)
            tuples.append((name, value.as_tuples()))

        # Retunr the tuples
        return tuples

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

# -----------------------------------------------------------------

def write_properties(properties, fh, indent=""):

    """
    This function ...
    :param properties:
    :param fh:
    :param indent:
    :return:
    """

    # Property composite
    if isinstance(properties, SimplePropertyComposite):

        # Loop over the properties
        for name in properties.property_names:

            dtype, value = stringify(getattr(properties, name))
            actual_dtype = properties._ptypes[name]
            print(indent + name + ":", value + " [" + actual_dtype + "]", file=fh)

        # Loop over the sections
        for name in properties.section_names:

            print(indent + name + ":", file=fh)
            write_properties(getattr(properties, name), fh, indent=indent+"  ")

    # Dict-like
    elif isinstance(properties, dict):

        # Loop over the keys
        for key in properties:

            value = properties[key]
            if isinstance(value, dict):

                print(indent + key + ":", file=fh)
                write_properties(properties[key], fh, indent=indent+"  ")

            else:

                dtype, string = stringify(value)
                print(indent + key + ":", string + " [" + dtype + "]", file=fh)

    # Invalid
    else: raise ValueError("Cannot write the properties of type '" + str(type(properties)))

# -----------------------------------------------------------------

def load_properties(properties, fh_or_iterator, indent=""):

    """
    Thisf unction ...
    :param properties:
    :param fh_or_iterator:
    :param indent:
    :return:
    """

    # Loop over the lines
    for line in fh_or_iterator:

        #if "Type:" in line: continue

        if line.endswith("\n"): line = line[:-1]

        # Empty line
        if not line: continue

        name, rest = line.split(":")
        rest = rest.strip()

        # End of section?
        if not name.startswith(indent): return

        # Not end of section: remove the indent
        name = name.lstrip(indent)

        if rest == "":

            #expect_section = name

            # Load properties
            properties[name] = dict()
            load_properties(properties[name], fh_or_iterator, indent=indent+"  ")

        else:

            value, dtype = rest.split(" [")
            dtype = dtype.split("]")[0]

            # Set the property value
            if dtype == "None" or value.strip() == "None": properties[name] = None
            else: properties[name] = getattr(parsing, dtype)(value)

# -----------------------------------------------------------------

