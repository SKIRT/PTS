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
import numpy as np

# Import the relevant PTS classes and modules
from ..tools import parsing
from ..basics.log import log
from ..tools import formatting as fmt
from ..tools.stringify import stringify, tostr, stringify_list_fancy
from ..tools.utils import defaultlen
from .map import Map

# -----------------------------------------------------------------

def load_composite(path):

    """
    This function ...
    :param path:
    :return:
    """

    from ..tools import introspection
    from ..tools import filesystem as fs

    first_line = fs.get_first_line(path)
    class_name = first_line.split("Type: ")[1].strip()

    # Get the PTS composite class
    cls = introspection.load_class(class_name)

    # Load the object
    return cls.from_file(path)

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

        # The descriptions: ordered to remember the order
        self._descriptions = OrderedDict()

        # The fixed flags
        self._fixed = dict()

        # The parsing types
        self._ptypes = dict()

        # The choices
        self._choices = dict()

    # -----------------------------------------------------------------

    def add_property(self, name, ptype, description, default_value=None, choices=None, convert_default=False, fixed=False):

        """
        This function ...
        :param name:
        :param ptype:
        :param description:
        :param default_value:
        :param choices:
        :param convert_default:
        :param fixed:
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

        # Set fixed flag to False
        self._fixed[name] = False  # to be able to set value

        # Convert default
        if default_value is not None and convert_default:
            parsing_function = getattr(parsing, ptype)
            default_value = parsing_function(default_value)

        # Set the attribute with the default value
        setattr(self, name, default_value)

        # Set the fixed flag
        self._fixed[name] = fixed

    # -----------------------------------------------------------------

    def add_fixed(self, name, description, value):

        """
        This function ...
        :param name:
        :param description:
        :param value:
        :return:
        """

        # Check
        if hasattr(self, name): raise ValueError("A property with the name '" + name + "' already exists")

        # Check
        if " " in name: raise ValueError("Name cannot contain spaces")

        # Get pytpe
        ptype, string = stringify(value)

        # Set the ptype
        self._ptypes[name] = ptype

        # Set the description
        self._descriptions[name] = description

        # Set fixed flag to False
        self._fixed[name] = False  # to be able to set value

        # Set the attribute with the default value
        setattr(self, name, value)

        # Set the fixed flag
        self._fixed[name] = True

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

    def add_section(self, name, description, dynamic=False, fixed=False):

        """
        This function ...
        :param name:
        :param description:
        :param dynamic:
        :param fixed:
        :return:
        """

        # Set the description
        self._descriptions[name] = description

        # Set fixed flag
        self._fixed[name] = fixed

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

        # Set the value
        self.set_value(name, value)

    # -----------------------------------------------------------------

    def set_value(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        #print(name, value)

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

        # None
        if value is None: pass

        # Property composite
        elif isinstance(value, SimplePropertyComposite): assert name in self._descriptions

        # Dict-like
        elif isinstance(value, dict):  # elif isinstance(value, Map):

            #print("NAME", name)
            #print("VALUE", value)

            # There is a subsection with this name: fill in the attributes for that subsection
            #if name in self._descriptions:

            assert name in self._descriptions

            if self.get_fixed(name): raise ValueError("This is a fixed section")

            #print(value)
            #print(getattr(self, name))
            #print(name in self.__dict__)

            # Add the section if not yet present or if None
            if name not in self.__dict__ or self.__dict__[name] is None:

                # Set dict
                self.__dict__[name] = value

            # Fill in sub-section properties
            elif isinstance(self.__dict__[name], SimplePropertyComposite):

                for key in value:
                    keyvalue = value[key]
                    # print(self.__dict__)
                    self.__dict__[name].__setattr__(key, keyvalue)
                return

            # Set already existing dictionary's items
            elif isinstance(self.__dict__[name], dict):

                for key in value: self.__dict__[name][key] = value[key]
                return

            # Just set the whole dictionary as an attribute
            #else: pass #self.__dict__[name] = value

        # Not dict-like
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
                #the_type = self._ptypes[name]
                #parsing_function = getattr(parsing, the_type)
                parsing_function = self.parser_for_property(name)
                the_type = self.type_for_property(name)
                try: value = parsing_function(string)
                except ValueError: raise ValueError("The value of '" + str(value) + "' for '" + name +  "' given is of the wrong type: '" + ptype + "', must be '" + the_type + "' (value is " + string + ")")

        # Check whether not fixed
        if self.get_fixed(name): raise ValueError("This is a fixed property or section")

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

    def show_properties(self, recursive=True, contains=None, not_contains=None, exact_name=None, exact_not_name=None,
                          startswith=None, endswith=None, label=None, bold=True):

        """
        This function ...
        :param recursive:
        :param contains:
        :param not_contains:
        :param exact_name:
        :param exact_not_name:
        :param startswith:
        :param endswith:
        :param label:
        :param bold:
        :return:
        """

        # Loop over the properties
        for name in self.property_names:

            #if suggestions is not None and name in suggestions: used_suggestions.append(name)

            # Checks
            if contains is not None and contains not in name: continue
            if not_contains is not None and not_contains in name: continue
            if exact_name is not None and name != exact_name: continue
            if exact_not_name is not None and name == exact_not_name: continue
            if startswith is not None and not name.startswith(startswith): continue
            if endswith is not None and not name.endswith(endswith): continue

            description = self.get_description(name)
            ptype = self.get_ptype(name)
            default = self.get_value(name)
            choices = self.get_choices(name)

            # Get suggestions
            #suggestns = suggestions[name] if suggestions is not None and name in suggestions else None

            # Add label to description
            if label is not None: description = description + " [" + label + "]"

            # Fixed variable: show value and description
            #if self.get_fixed(name):
            #    value = prompt_fixed(name, description, self.get_value(name))
            #    continue

            if self.get_fixed(name): suffix = " FIXED"
            else: suffix = ""

            # Set ptype string
            if defaultlen(default, 1) == 0: ptype_string = "[empty " + ptype + "]"
            else: ptype_string = " [" + ptype + "]"

            # Show
            if bold: print(" - " + fmt.bold + name + fmt.reset + " (" + description + "): " + tostr(default) + ptype_string + suffix)
            else: print(" - " + name + " (" + description + "): " + tostr(default) + ptype_string + suffix)

        # Recursive: also loop over the sections
        if recursive:

            # Loop over the sections
            for name in self.section_names:

                # Show the properties of the section
                self.sections[name].show_properties(recursive=True, contains=contains, not_contains=not_contains, exact_name=exact_name, exact_not_name=exact_not_name, bold=bold)

    # -----------------------------------------------------------------

    def prompt_properties(self, recursive=True, contains=None, not_contains=None, exact_name=None, exact_not_name=None,
                          startswith=None, endswith=None, required=True, label=None, suggestions=None,
                          add_suggestions=False, replace_string=None, types=None, only_replacements=False):

        """
        This function ...
        :param recursive:
        :param contains:
        :param not_contains:
        :param exact_name:
        :param exact_not_name:
        :param startswith:
        :param endswith:
        :param required:
        :param label:
        :param suggestions:
        :param add_suggestions:
        :param replace_string:
        :param types:
        :param only_replacements:
        :return:
        """

        from .configuration import prompt_variable, prompt_fixed

        has_changed = False
        used_suggestions = []

        # Loop over the properties
        for name in self.property_names:

            if suggestions is not None and name in suggestions: used_suggestions.append(name)

            # Checks
            if contains is not None and contains not in name: continue
            if not_contains is not None and not_contains in name: continue
            if exact_name is not None and name != exact_name: continue
            if exact_not_name is not None and name == exact_not_name: continue
            if startswith is not None and not name.startswith(startswith): continue
            if endswith is not None and not name.endswith(endswith): continue

            description = self.get_description(name)
            ptype = self.get_ptype(name)
            default = self.get_value(name)
            choices = self.get_choices(name)

            # Check type
            if types is not None and ptype not in types: continue

            # Get suggestions
            suggestns = suggestions[name] if suggestions is not None and name in suggestions else None

            # Add label to description
            if label is not None: description = description + " [" + label + "]"

            # Fixed variable: show value and description
            if self.get_fixed(name):
                value = prompt_fixed(name, description, self.get_value(name))
                continue

            # Replace?
            if ptype == "string" and replace_string is not None and default is not None and replace_string[0] in default: value = default.replace(replace_string[0], replace_string[1])

            # No replacement: skip
            elif only_replacements: continue

            # No replacement: prompt for value
            else:

                # Ask for the new value
                value = prompt_variable(name, ptype, description, choices=choices, default=default, required=required, suggestions=suggestns)
                if default is None and value == "": continue

            # Set the property
            if value != default:

                # Debugging
                log.debug("Changing the value of '" + name + "' to '" + tostr(value) + "' ...")
                log.debug("Original value: '" + tostr(default) + "'")

                # Set the new value
                self.set_property(name, value)

                # Set flag
                has_changed = True

        # Recursive: also loop over the sections
        if recursive:

            # Loop over the sections
            for name in self.section_names:

                # Prompt the settings of the section
                has_changed_section = self.sections[name].prompt_properties(recursive=True, contains=contains, not_contains=not_contains,
                                                                            exact_name=exact_name, exact_not_name=exact_not_name,
                                                                            startswith=startswith, endswith=endswith,
                                                                            required=required, label=label, suggestions=suggestions,
                                                                            add_suggestions=add_suggestions, replace_string=replace_string,
                                                                            types=types, only_replacements=only_replacements)
                if has_changed_section: has_changed = True

        # Add suggested
        if suggestions is not None and add_suggestions:
            for name in suggestions:

                if name in used_suggestions: continue
                values = suggestions[name]
                if len(values) > 1: raise ValueError("Multiple suggestions")
                value = values[0]
                self.add_fixed(name, "no description", value)
                has_changed = True

                # Show the suggested value
                description = "no description"
                value = prompt_fixed(name, description, value)

        # Return whether any property changed
        return has_changed

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

    def get_value(self, name, add_unit=True, unit=None):

        """
        This funtion ...
        :param name:
        :param add_unit:
        :param unit:
        :return:
        """

        value = getattr(self, name)

        # Return
        if value is None: return None
        elif self.has_unit(name):
            if add_unit: return value
            else:
                if unit is None: return value.value
                else: return value.to(unit).value
        else: return value

    # -----------------------------------------------------------------

    def __getitem__(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_value(name)

    # -----------------------------------------------------------------

    def __setitem__(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        return self.set_value(name, value)

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

    def get_fixed(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self._fixed[name]

    # -----------------------------------------------------------------

    def set_fixed(self, name, value=None):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        if value is not None:
            self._fixed[name] = False # to be able to set value
            self.set_property(name, value)
        self._fixed[name] = True

    # -----------------------------------------------------------------

    @property
    def property_names(self):

        """
        This function ...
        :return:
        """

        from ..tools.utils import lazyproperty

        names = []
        for name in vars(self):

            # Skip internal variables
            if name.startswith("_"): continue

            #print(type(getattr(self, name)))

            if hasattr(self.__class__, name):
                class_attribute = getattr(self.__class__, name)
                assert isinstance(class_attribute, lazyproperty)
                continue

            if name not in self._ptypes:
                if not (isinstance(getattr(self, name), SimplePropertyComposite) or isinstance(getattr(self, name), Map)): raise Exception("Property '" + name + "' doesn't have its type defined")
                else: continue

            # Add the name
            names.append(name)

        # Return the names
        return names

    # -----------------------------------------------------------------

    @property
    def ordered_property_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for name in self._descriptions:

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

        from ..tools.utils import lazyproperty

        names = []
        for name in vars(self):

            # Skip internal
            if name.startswith("_"): continue

            if hasattr(self.__class__, name):
                class_attribute = getattr(self.__class__, name)
                assert isinstance(class_attribute, lazyproperty)
                continue

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
    def ordered_section_names(self):

        """
        This function ...
        :return:
        """

        names = []

        for name in self._descriptions:

            # Skip simple properties
            if name in self._ptypes: continue
            if not (isinstance(getattr(self, name), SimplePropertyComposite) or isinstance(getattr(self, name), Map)): raise Exception("Section '" + name + "' is not a property composite or mapping")

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

    def parser_for_property(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        ptype = self.type_for_property(name)
        return getattr(parsing, ptype)

    # -----------------------------------------------------------------

    def get_type(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.type_for_property(name)

    # -----------------------------------------------------------------

    def get_parent_type(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        from .configuration import parent_type
        return parent_type(self.get_type(name))

    # -----------------------------------------------------------------

    def is_boolean(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_parent_type(name) == "boolean"

    # -----------------------------------------------------------------

    def is_integer(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_parent_type(name) == "integer"

    # -----------------------------------------------------------------

    def is_real(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_parent_type(name) == "real"

    # -----------------------------------------------------------------

    def is_string(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_parent_type(name) == "string"

    # -----------------------------------------------------------------

    def is_quantity(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_parent_type(name) == "quantity"

    # -----------------------------------------------------------------

    def is_unit(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_parent_type(name) == "unit"

    # -----------------------------------------------------------------

    def is_filter(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_parent_type(name) == "filter"

    # -----------------------------------------------------------------

    def is_list(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_parent_type(name) == "list"

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

    def has_unit(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        value = getattr(self, name)
        return hasattr(value, "unit")

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

    def differences(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        differences = OrderedDict()
        none_a = []
        none_b = []
        different = []

        # Loop over the parameters
        for name in self.property_names:

            # Get the values
            value_a = self.get_value(name)
            value_b = other.get_value(name)

            # Both None
            if value_a is None and value_b is None: continue

            # Value b is None
            if value_a is not None and value_b is None:
                none_b.append(name)
                continue

            # Value a is None
            if value_b is not None and value_a is None:
                none_a.append(name)
                continue

            # Get difference
            if self.is_integer(name) or self.is_real(name) or self.is_quantity(name): differences[name] = value_a - value_b
            elif value_a != value_b: different.append(name)
            else: continue

        # Return
        return differences, none_a, none_b, different

    # -----------------------------------------------------------------

    def isclose(self, other, ignore_none=False, ignore_other=False, rtol=1e-05, atol=1e-08):

        """
        This function ...
        :param other:
        :param ignore_none:
        :param ignore_other: check things other than numbers
        :param rtol:
        :param atol:
        :return:
        """

        # Get differences
        differences, none_a, none_b, different = self.differences(other)

        # Check none
        if not ignore_none and len(none_a) > 0: return False
        if not ignore_none and len(none_b) > 0: return False

        # Check other
        if not ignore_other and len(different) > 0: return False

        # Check numbers
        for name in differences:

            # Get unit for this property
            if self.has_unit(name): unit = self.unit_for_property(name)
            else: unit = None

            # Get the values
            value_a = self.get_value(name, unit=unit, add_unit=False)
            value_b = other.get_value(name, unit=unit, add_unit=False)

            #print(value_a, value_b, np.isclose(value_a, value_b, rtol=rtol, atol=atol))

            # Check whether close
            if not np.isclose(value_a, value_b, rtol=rtol, atol=atol): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def __eq__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Get differences
        differences, none_a, none_b, different = self.differences(other)

        if len(none_a) > 0: return False
        if len(none_b) > 0: return False
        if len(different) > 0: return False

        # Check numbers
        for name in differences:

            value_a = self.get_value(name)
            value_b = other.get_value(name)
            if value_a != value_b: return False

        # ALl checks passed
        return True

    # -----------------------------------------------------------------

    def to_lines(self, line_prefix="", ignore=None, ignore_none=False, bullet="-", bold=True):

        """
        This function ...
        :param line_prefix:
        :param ignore:
        :param ignore_none:
        :param bullet:
        :param bold:
        :return:
        """

        from ..tools import types

        lines = []

        # The simple properties
        for name in self.ordered_property_names:

            # Get the value
            value = getattr(self, name)

            # Check ignore
            if ignore is not None and name in ignore: continue
            if ignore_none and value is None: continue

            # Multiline sequences
            if types.is_sequence(value):

                if bold: line =  " " + bullet + " " + fmt.bold + name + fmt.reset + ": "
                else: line =  " " + bullet + " " + name + ": "

                lines.append(line_prefix + line)

                if len(value) > 0:
                    #lines.append("")
                    dtype, string = stringify_list_fancy(value)
                    for line in string.split("\n"): lines.append(line_prefix + "    " + line)
                    #lines.append("")

            # Multiline dictionaries
            elif types.is_dictionary(value):

                if bold: line = " " + bullet + " " + fmt.bold + name + fmt.reset + ": "
                else: line = " " + bullet + " " + name + ": "

                lines.append(line_prefix + line)

                if len(value) > 0:
                    for key in value:
                        line = line_prefix + "    * " + key + ": " + tostr(value[key])
                        lines.append(line)

            # Single-line
            else:

                # Generate line
                dtype, value = stringify(value)

                if bold: line = " " + bullet + " " + fmt.bold + name + fmt.reset + ": " + value
                else: line = " " + bullet + " " + name + ": " + value

                lines.append(line_prefix + line)

        # Sections
        for name in self.ordered_section_names:

            # Check ignore
            if ignore is not None and name in ignore: continue

            if bold: line = " - " + fmt.bold + name + fmt.reset + ":"
            else: line = " - " + name + ":"

            lines.append(line_prefix + line)

            if isinstance(getattr(self, name), SimplePropertyComposite): section_lines = [line for line in getattr(self, name).to_lines(line_prefix=line_prefix+"    ")]
            elif isinstance(getattr(self, name), Map):

                section_lines = []
                section = getattr(self, name)
                for key in section:

                    if bold: line = line_prefix + "    " + " " + bullet + " " + fmt.bold + key + fmt.reset + ": " + tostr(section[key])
                    else: line = line_prefix + "    " + " " + bullet + " " + key + ": " + tostr(section[key])

                    section_lines.append(line)

            else: raise ValueError("Unknown type for section: " + str(type(getattr(self, name))))

            lines += section_lines

        # Return the lines
        return lines

    # -----------------------------------------------------------------

    def to_string(self, line_prefix="", ignore=None, ignore_none=False, bullet="-", bold=True):

        """
        This function ...
        :param line_prefix:
        :param ignore:
        :param ignore_none:
        :param bullet:
        :param bold:
        :return:
        """

        # Return
        return "\n".join(self.to_lines(line_prefix=line_prefix, ignore=ignore, ignore_none=ignore_none, bullet=bullet, bold=bold))

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
        #if cls.__name__ == "DeprojectionModel3D": raise ValueError("Test")

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

        # Add directory path
        from ..tools import filesystem as fs
        dirpath = fs.directory_of(path)
        properties["dirpath"] = dirpath

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

            #print(name, properties[name])

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

