#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.table Contains the SmartTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import StringIO
from collections import OrderedDict

# Import astronomical modules
from astropy.table import Table, MaskedColumn

# Import the relevant PTS classes and modules
from ..units.unit import PhotometricUnit
from ..units.parsing import parse_unit as u
from ..tools import filesystem as fs
from ..tools import types
from ..tools import sequences

# -----------------------------------------------------------------

class SmartTable(Table):

    """
    This class ...
    """

    # default extension
    default_extension = "dat"

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Always used masked tables
        kwargs["masked"] = True

        # Call the constructor of the base class
        super(SmartTable, self).__init__(*args, **kwargs)

        # Column info
        self.column_info = []

        # Path
        self.path = None

        # Initialize 'density' meta object
        if "density" not in self.meta: self.meta["density"] = []
        if "brightness" not in self.meta: self.meta["brightness"] = []

        # The column descriptions
        self._descriptions = dict()

    # -----------------------------------------------------------------

    @classmethod
    def from_composites(cls, *composites, **kwargs):

        """
        This function ...
        :param composites:
        :param kwargs:
        :return:
        """

        # Import tostr function
        from ..tools.stringify import tostr

        # Check number
        if len(composites) == 0: raise ValueError("No input is provided")

        # Get names
        labels = kwargs.pop("labels", None)
        label = kwargs.pop("label", "-")

        # Get tostr kwargs
        tostr_kwargs = kwargs.pop("tostr_kwargs", {})

        #print(composites)

        # Get list of lists of property names
        property_names = [composite.property_names for composite in composites]
        #property_types = [composite.property_types for composite in composites]
        #property_units = [composite.property_units for composite in composites]
        #property_descriptions = [composite.property_descriptions for composite in composites]

        # Determine the column names, types and descriptions
        column_names = sequences.union(*property_names)
        #column_units = sequences.union(*property_units)
        #column_types = sequences.union(*property_types)
        #column_descriptions = sequences.union(*property_descriptions)

        # Add the label column name
        if labels is not None: sequences.prepend(column_names, label)
        else: labels = [None] * len(composites) # labels was None

        # Create the table
        table = cls()

        # Create lists to contain the column types and units
        column_types = []
        column_units = []

        # COLUMNS FOR WHICH THE VALUE HAS TO BE CONVERTED TO STRING
        to_string = []

        # Make the columns
        for name in column_names:

            if name == label:

                real_type = str
                column_unit = None
                column_description = label

            else:

                # Get the types, units and descriptions
                types = [composite.type_for_property(name) for composite in composites]
                units = [composite.unit_for_property(name) for composite in composites]
                descriptions = [composite.description_for_property(name) for composite in composites]

                # Determine column type, unit and description
                if sequences.all_equal_to(types, 'None') or sequences.all_none(types): column_type = 'None'
                else: column_type = sequences.get_all_equal_value(types, ignore_none=True, ignore='None')
                column_unit = sequences.get_first_not_none_value(units)
                column_description = sequences.get_first_not_none_value(descriptions)

                column_types.append(column_type)
                column_units.append(column_unit)

                if column_type.endswith("real"): real_type = float
                elif column_type.endswith("integer"): real_type = int
                elif column_type.endswith("string"): real_type = str
                elif column_type.endswith("boolean"): real_type = bool
                elif column_type.endswith("quantity"): real_type = float

                # UNIT HAVE TO BE CONVERTED TO STRINGS
                elif column_type.endswith("unit"):
                    to_string.append(name)
                    real_type = str

                # FILTERS HAVE TO BE CONVERTED TO STRINGS
                elif column_type.endswith("filter"):
                    to_string.append(name)
                    real_type = str

                # SPECIAL CASE: WAS NONE FOR EACH COMPOSITE
                elif column_type == "None": real_type = str

                # LISTS OF THINGS -> STRINGS
                elif column_type.endswith("_list"):
                    to_string.append(name)
                    real_type = str

                # NOT RECOGNIZED
                else: raise ValueError("Column type not recognized: " + str(column_type) + " (" + str(type(column_type)) + ")")

            # Add the column info
            table.add_column_info(name, real_type, column_unit, column_description)

        # Make the columns
        #for name, dtype, description in zip(column_names, column_types, column_units, column_descriptions):
        #    # Add column info
        #    table.add_column_info(name, simple_dtype, unit, description)

        # Add the rows
        for composite_label, composite in zip(labels, composites):

            values = []

            # Fill the row
            #for name, dtype, unit in zip(column_names, column_types, column_units):
            for name in column_names:

                if name == label: value = composite_label
                else:

                    # Get the value
                    if hasattr(composite, name):
                        value = getattr(composite, name)
                        if name in to_string: value = tostr(value, **tostr_kwargs)
                    else: value = None

                # Add the value
                values.append(value)

            # Add the row: unit conversion is done here
            table.add_row(values)

        # Return the table
        return table

    # -----------------------------------------------------------------

    def add_column_info(self, name, dtype, unit, description):

        """
        This function ...
        :param name:
        :param dtype:
        :param unit:
        :param description:
        :return:
        """

        if types.is_string_type(unit): unit = u(unit)
        self.column_info.append((name, dtype, unit, description))

    # -----------------------------------------------------------------

    def column_unit(self, column_name):

        """
        This function ...
        :param column_name:
        :return:
        """

        if self[column_name].unit is None: return None

        # Construct unit
        if column_name in self.meta["density"]: density = True
        else: density = False
        if column_name in self.meta["brightness"]: brightness = True
        else: brightness = False
        return u(self[column_name].unit, density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def column_unit_string(self, column_name):

        """
        This function ...
        :param column_name:
        :return:
        """

        unit = self.column_unit(column_name)
        if unit is None: return ""
        else: return str(unit)

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        #print("Performing setup of the table ...")

        # Setup has already been called
        if len(self.colnames) != 0: return

        # Create the table names and types lists
        for entry in self.column_info:

            name = entry[0]
            dtype = entry[1]
            unit = entry[2]
            description = entry[3]

            data = []

            # Add column
            col = MaskedColumn(data=data, name=name, dtype=dtype, unit=unit)
            self.add_column(col)

            # Set the description
            self._descriptions[name] = description

            # Set whether this column is a spectral density
            if isinstance(unit, PhotometricUnit) and unit.density:
                if "density" not in self.meta: self.meta["density"] = []
                self.meta["density"].append(name)

            # Set whether this column is a surface brightness
            if isinstance(unit, PhotometricUnit) and unit.brightness:
                if "brightness" not in self.meta: self.meta["brightness"] = []
                self.meta["brightness"].append(name)

    # -----------------------------------------------------------------

    def __getitem__(self, item):

        """
        This function ...
        :param item: 
        :return: 
        """

        # Run the setup if not yet performed
        if len(self.colnames) == 0: self.setup()

        # Call the implementation of the base class
        return super(SmartTable, self).__getitem__(item)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :return:
        """

        fill_values = [('--', '0')]
        #fill_values = [("", '0')]
        #fill_values = None

        # Check the path
        if not fs.is_file(path): raise IOError("The file '" + path + "' does not exist")

        # Open the table
        table = super(SmartTable, cls).read(path, format="ascii.ecsv", fill_values=fill_values)

        #print(table.meta)

        # Set masks
        #if "masks" in table.meta:
        #    for name in table.meta["masks"]:
        #        for index in range(len(table)):
        #            table[name].mask[index] = table.meta["masks"][name][index]
        #    # Remove masks from meta
        #    del table.meta["masks"]

        # Look for masks
        for colname in table.colnames:
            key = colname + " mask"
            if key not in table.meta: continue
            if len(table) != len(table.meta[key]): raise IOError("Length of the table does not correspond to the length of the masks")
            for index in range(len(table)):
                table[colname].mask[index] = table.meta[key][index]
            del table.meta[key]

        # Set the path
        table.path = path

        # Clear the column info so that we can rebuild it
        table.column_info = []

        # Set the column info
        # Loop over the columns
        for name in table.colnames:

            # Get the type
            dtype = table[name].dtype
            if np.issubdtype(dtype, np.string_): simple_dtype = str
            elif np.issubdtype(dtype, np.float): simple_dtype = float
            elif np.issubdtype(dtype, np.int): simple_dtype = int
            elif np.issubdtype(dtype, np.bool): simple_dtype = bool
            else: raise ValueError("Did not recognize the dtype of column '" + name + "'")

            # Get unit of the column
            unit = table.column_unit(name)

            # Add column info
            table.add_column_info(name, simple_dtype, unit, None)

            # Initialize "density" meta
            if "density" not in table.meta: table.meta["density"] = []

            # Initialize "brightness" meta
            if "brightness" not in table.meta: table.meta["brightness"] = []

        # Return the table
        return table

    # -----------------------------------------------------------------

    def _resize_string_columns(self, values):

        """
        This function ...
        :return:
        """

        new_sizes = dict()

        colnames = self.colnames
        for index in range(len(colnames)):

            colname = colnames[index]

            dtype_str = str(self[colname].dtype)

            if not dtype_str.startswith("|S"): continue

            current_string_length = int(dtype_str.split("S")[1])
            if values[index] is None: new_string_length = 0
            else: new_string_length = len(values[index])

            if new_string_length > current_string_length: new_sizes[colname] = new_string_length

        # Resize columns
        for colname in new_sizes:

            current_resized = self[colname].astype("S" + str(new_sizes[colname]))
            self.replace_column(colname, current_resized)

    # -----------------------------------------------------------------

    def _resize_string_column(self, colname, value):

        """
        This function ...
        :param colname:
        :param value:
        :return:
        """

        if value is None: return

        dtype_str = str(self[colname].dtype)

        if not dtype_str.startswith("|S"): raise ValueError("Column " + colname + " is not a column of strings")

        current_string_length = int(dtype_str.split("S")[1])

        new_string_length = len(value)

        if new_string_length > current_string_length:

            new_size = new_string_length

            # Replace the column by a resized one
            current_resized = self[colname].astype("S" + str(new_size))
            self.replace_column(colname, current_resized)

    # -----------------------------------------------------------------

    def _strip_units(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        scalar_values = []

        for i, value in enumerate(values):

            # If this value has a unit, we have to make sure it is converted into the proper column unit
            if hasattr(value, "unit"):

                #print(self.column_info)

                column_unit = self.column_info[i][2]
                assert column_unit is not None

                # Quantity with photometric unit
                if isinstance(value.unit, PhotometricUnit):

                    factor = value.unit.conversion_factor(column_unit)
                    scalar_value = value.value * factor

                # Quantity with regular Astropy Unit
                else: scalar_value = value.to(column_unit).value

                # Add the value
                scalar_values.append(scalar_value)

            # A scalar value (or string, int, ...)
            else: scalar_values.append(value)

        # Return the values without unit
        return scalar_values

    # -----------------------------------------------------------------

    def _convert_lists(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        converted_values = []

        for i, value in enumerate(values):

            if isinstance(value, list):

                column_type = self.column_info[i][1]

                if len(value) == 0: converted_value = None
                elif column_type == str: converted_value = ",".join(map(str, value))
                else: raise ValueError("Cannot have a list element in the row at the column that is not of string type")

            else: converted_value = value

            converted_values.append(converted_value)

        # Return the converted values
        return converted_values

    # -----------------------------------------------------------------

    def get_quantity(self, colname, index):

        """
        This function ...
        :param colname:
        :param index:
        :return:
        """

        value = self[colname][index]

        if self[colname].mask[index]: return None
        elif self[colname].unit is not None:
            quantity = value * self[colname].unit
        else: quantity = value

        # Return the quantity
        return quantity

    # -----------------------------------------------------------------

    def get_value(self, colname, index, add_unit=True):

        """
        This function ...
        :param colname:
        :param index:
        :param add_unit:
        :return:
        """

        value = self[colname][index]

        if self[colname].mask[index]: value = None
        elif self[colname].unit is not None and add_unit: value = value * self[colname].unit

        # Return the value
        return value

    # -----------------------------------------------------------------

    def is_masked_value(self, colname, index):

        """
        This function ...
        :param colname:
        :param index:
        :return:
        """

        return self[colname].mask[index]

    # -----------------------------------------------------------------

    def get_row(self, index, add_units=True):

        """
        This function ...
        :param index:
        :param add_units:
        :return:
        """

        row = OrderedDict()

        for name in self.colnames:

            # Get the value
            value = self.get_value(name, index, add_unit=add_units)

            # Add the value
            row[name] = value

        # Return the row
        return row

    # -----------------------------------------------------------------

    def add_row(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        # Setup if necessary
        if len(self.colnames) == 0: self.setup()

        #print(len(self.colnames))
        #print(len(values))

        # CHECK TYPES BEFORE RESIZE STRING COLUMNS?

        # Resize string columns for the new values
        self._resize_string_columns(values)

        # Strip units
        values = self._strip_units(values)

        # Convert lists to string
        values = self._convert_lists(values)

        # Create mask
        mask = [value is None for value in values]

        # Set masked values to have a default value (None will not work for Astropy)
        new_values = []
        for i in range(len(values)):

            if values[i] is not None: new_values.append(values[i])
            else:
                colname = self.colnames[i]
                if self.is_string_type(colname): new_values.append("")
                elif self.is_real_type(colname): new_values.append(0.)
                elif self.is_integer_type(colname): new_values.append(0)
                elif self.is_boolean_type(colname): new_values.append(False)
                else: raise ValueError("Unknown column type for '" + colname + "'")

        # Add the row
        super(SmartTable, self).add_row(new_values, mask=mask)

    # -----------------------------------------------------------------

    def column_type(self, column_name):

        """
        This function ...
        :param column_name: 
        :return: 
        """

        coltype = self[column_name].dtype.name

        if coltype.startswith("string"): return "string"
        elif coltype.startswith("float"): return "real"
        elif coltype.startswith("int"): return "integer"
        elif coltype.startswith("bool"): return "boolean"
        else: raise ValueError("Unknown column type: " + coltype)

    # -----------------------------------------------------------------

    def is_string_type(self, column_name):

        """
        This function ...
        :param column_name: 
        :return: 
        """

        return self.column_type(column_name) == "string"

    # -----------------------------------------------------------------

    def is_real_type(self, column_name):

        """
        This function ...
        :param column_name: 
        :return: 
        """

        return self.column_type(column_name) == "real"

    # -----------------------------------------------------------------

    def is_integer_type(self, column_name):

        """
        This function ...
        :param column_name: 
        :return: 
        """

        return self.column_type(column_name) == "integer"

    # -----------------------------------------------------------------

    def is_boolean_type(self, column_name):

        """
        This function ...
        :param column_name: 
        :return: 
        """

        return self.column_type(column_name) == "boolean"

    # -----------------------------------------------------------------

    def all_equal(self, column_name):

        """
        This function ...
        :param column_name:
        :return:
        """

        if len(self[column_name]) == 0: return True

        # Doesn't work with strings I think ...
        #return (self[column_name] == self[column_name][0]).all() # all equal to the first element

        first = self[column_name][0]
        for i in range(len(self[column_name])):
            if self[column_name][i] != first: return False
        return True

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        if self.path is None: raise RuntimeError("Path has not been set yet")

        # Save to the current path
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        #print(self.colnames)

        # Setup if necessary
        if len(self.colnames) == 0: self.setup()

        # Get masks
        masks = self.get_masks()

        # Replace masked values (not masked anymore)
        self.replace_masked_values()

        # Set masks in meta
        for name in masks: self.meta[name + " mask"] = masks[name]

        # If the file already exists, remove
        if fs.is_file(path): fs.remove_file(path)

        # Write the table in ECSV format
        #self.write(path, format="ascii.ecsv", overwrite=True)
        self.write(path, format="ascii.ecsv")

        # Set the path
        self.path = path

        # Set the masks back (because they were set to False by replace_masked_values, necessary to avoid writing out
        # '""' (empty string) for each masked value, which is unreadable by Astropy afterwards)
        self.set_masks(masks)

    # -----------------------------------------------------------------

    def get_masks(self):

        """
        This function ...
        :return: 
        """

        masks = dict()
        for name in self.colnames:
            masks[name] = [int(boolean) for boolean in self[name].mask] #list(self[name].mask)
        return masks

    # -----------------------------------------------------------------

    def set_masks(self, masks):

        """
        This function ...
        :param masks:
        :return: 
        """

        # Loop over the columns for which there is a mask in the 'masks' dictionary
        for colname in masks:

            # Loop over the rows, set mask elements
            for index in range(len(self)): self[colname].mask[index] = masks[colname][index]

    # -----------------------------------------------------------------

    def replace_masked_values(self):

        """
        This function ...
        :return: 
        """

        # Loop over the columns
        for colname in self.colnames:

            # Loop over the rows
            for index in range(len(self)):

                # If not masked, skip
                if not self[colname].mask[index]: continue

                # Set value
                if self.is_string_type(colname): value = ""
                elif self.is_real_type(colname): value = 0.
                elif self.is_integer_type(colname): value = 0
                elif self.is_boolean_type(colname): value = False
                else: raise ValueError("Unknown column type for '" + colname + "'")

                # Set value
                self[colname][index] = value

    # -----------------------------------------------------------------

    def set_masked_constants(self):

        """
        This function ...
        :return: 
        """

        # Loop over the columns
        for colname in self.colnames:

            # Loop over the rows
            for index in range(len(self)):

                # If not masked, skip
                if not self[colname].mask[index]: continue

                # Else, set masked value
                #if self.column_type(colname) == "string": self[colname][index] = "--"
                #elif self.column_type(colname) == "integer": self[colname][index] = 0
                #elif self.column_type(colname) == "real": self[colname][index] = 0
                #elif self.column_type(colname) == "boolean": self[colname][index] = False

                self[colname][index] = np.ma.masked
                self[colname].mask[index] = True

    # -----------------------------------------------------------------

    def to_html(self):

        """
        This function ...
        :return:
        """

        # Make string output
        output = StringIO.StringIO()

        # Write as HTML
        self.write(output, format='html')
        contents = output.getvalue()

        # Close object and discard memory buffer --
        # .getvalue() will now raise an exception.
        output.close()

        # Return the HTML string
        return contents

    # -----------------------------------------------------------------

    @property
    def column_names(self):

        """
        Thisn function ...
        :return:
        """

        return self.colnames

    # -----------------------------------------------------------------

    @property
    def units(self):

        """
        This function ...
        :return:
        """

        units = []
        for name in self.column_names: units.append(self.column_unit(name))
        return units

    # -----------------------------------------------------------------

    @property
    def unit_strings(self):

        """
        This function ...
        :return:
        """

        strings = []
        for name in self.column_names: strings.append(self.column_unit_string(name))
        return strings

    # -----------------------------------------------------------------

    @property
    def descriptions(self):

        """
        This function ...
        :return:
        """

        strings = []
        for name in self.column_names:
            if name in self._descriptions: description = self._descriptions[name]
            else: description = None
            strings.append(description)
        return strings

    # -----------------------------------------------------------------

    @property
    def ncolumns(self):

        """
        This function ...
        :return:
        """

        return len(self.column_names)

    # -----------------------------------------------------------------

    @property
    def nrows(self):

        """
        This function ...
        :return:
        """

        return len(self)

    # -----------------------------------------------------------------

    def as_tuples(self, add_units=True):

        """
        This function ...
        :param add_units:
        :return:
        """

        tuples = []
        for index in range(len(self)):
            values = self.get_row(index, add_units=add_units).values()
            tuples.append(tuple(values))
        return tuples

    # -----------------------------------------------------------------

    def as_lists(self, add_units=True):

        """
        This function ...
        :param add_units:
        :return:
        """

        lists = []
        for index in range(len(self)):
            values = self.get_row(index, add_units=add_units).values()
            lists.append(values)
        return lists

    # -----------------------------------------------------------------

    def print_latex(self):

        """
        This function ...
        :return: 
        """

        self.show_in_browser()

        header = " & ".join(self.colnames) + " \\\\"
        print(header)

        units = []
        for name in self.colnames:
            unit = self.column_unit(name)
            if unit is None:
                units.append("")
            else:
                units.append(str(unit))
        units_string = " & ".join(units) + " \\\\"
        print(units_string)

        for index in range(len(self)):

            row = []
            for name in self.colnames: row.append(str(self[name][index]))
            row_string = " & ".join(row) + " \\\\"

            print(row_string)

# -----------------------------------------------------------------
