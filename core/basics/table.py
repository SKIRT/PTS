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
import warnings
from collections import OrderedDict, defaultdict
import copy

# Import astronomical modules
from astropy.table import Table, MaskedColumn
from astropy.io.registry import IORegistryError

# Import the relevant PTS classes and modules
from ..units.unit import PhotometricUnit, get_common_unit
from ..units.parsing import parse_unit as u
from ..tools import filesystem as fs
from ..tools import types
from ..tools import sequences
from .containers import DefaultOrderedDict

# -----------------------------------------------------------------

string_column_none_default = "--"

# -----------------------------------------------------------------

def is_pts_data_format(path):

    """
    This function ...
    :param path:
    :return:
    """

    first = fs.get_first_line(path)
    return "PTS data format" in first

# -----------------------------------------------------------------

def has_same_values_for_column(tables, column_name, reference_column_name, reference_values):

    """
    This function ...
    :param tables:
    :param column_name:
    :param reference_column_name:
    :param reference_values:
    :return:
    """

    from ..tools.tables import find_index

    for value in reference_values:

        table_column_values = []

        for table in tables:

            index = find_index(table, value, column_name=reference_column_name)
            if index is None: continue
            table_column_values.append(table[column_name][index])

        #print(table_column_values)

        if not sequences.all_equal(table_column_values): return False

    return True

# -----------------------------------------------------------------

def merge_tables(*tables, **kwargs):

    """
    This function ...
    :param tables:
    :param kwargs:
    :return:
    """

    from ..tools.tables import find_index

    # Get flags
    differences = kwargs.pop("differences", False)
    rel_differences = kwargs.pop("rel_differences", False)
    percentual = kwargs.pop("percentual", False)

    # Columns to use
    columns = kwargs.pop("columns", None)
    not_columns = kwargs.pop("not_columns", None)

    # Get table labels
    labels = kwargs.pop("labels", None)
    if labels is not None and len(labels) != len(tables): raise ValueError("Number of labels must be equal to number of tables")

    # Only shared columns
    only_shared = kwargs.pop("only_shared", False)

    # Determine column name
    reference_column_name = kwargs.pop("column_name", None)
    if reference_column_name is None:

        first_column_names = [table.colnames[0] for table in tables]
        if not sequences.all_equal(first_column_names): raise ValueError("Tables have different first columns")
        reference_column_name = first_column_names[0]

    reference_column_dtypes = []
    reference_column_units = []
    for table in tables:
        reference_column_dtypes.append(table.get_column_dtype(reference_column_name))
        reference_column_units.append(table.get_column_unit(reference_column_name))
    # Check
    if not sequences.all_equal(reference_column_dtypes): raise ValueError("Different reference column dtypes")
    if not sequences.all_equal(reference_column_units): raise ValueError("Different reference column units")
    reference_column_dtype = reference_column_dtypes[0]
    reference_column_unit = reference_column_units[0]

    all_column_types = DefaultOrderedDict(list)
    all_column_units = defaultdict(list)
    all_column_tableids = defaultdict(list)

    for table_index, table in enumerate(tables):
        #for colname in table.colu
        for name, dtype, unit, description in table.column_info:
            if name == reference_column_name: continue
            all_column_types[name].append(dtype)
            all_column_units[name].append(unit)
            all_column_tableids[name].append(table_index)

    all_values = []
    for table in tables: all_values.extend(table[reference_column_name])
    # print(all_values)
    reference_values = sequences.unique_values(all_values, ignore_none=True)

    # Create a table
    merged = SmartTable()

    unique_column_names = dict()
    equal_column_names = []
    original_column_names = dict()
    new_column_names = dict()

    # Add the columns
    merged.add_column_info(reference_column_name, reference_column_dtype, reference_column_unit, None)
    for column_name in all_column_types:

        # Use column?
        if columns is not None and column_name not in columns: continue
        if not_columns is not None and column_name in not_columns: continue

        # Get types and units
        dtypes = all_column_types[column_name]
        units = all_column_units[column_name]
        tableids = all_column_tableids[column_name]

        if len(dtypes) == 1:

            if only_shared: continue
            assert len(units) == 1
            assert len(tableids) == 1
            merged.add_column_info(column_name, dtypes[0], units[0], None)
            unique_column_names[column_name] = tableids[0]

        else:

            # Check if the columns are equal
            if has_same_values_for_column(tables, column_name, reference_column_name, reference_values):

                dtype = tables[tableids[0]].get_column_dtype(column_name)
                unit = tables[tableids[0]].get_column_unit(column_name)
                merged.add_column_info(column_name, dtype, unit, None)
                equal_column_names.append(column_name)

            else:

                for dtype, unit, tableid in zip(dtypes, units, tableids):
                    if labels is not None: new_column_name = column_name + " " + labels[tableid]
                    else: new_column_name = column_name + " " + str(tableid)
                    merged.add_column_info(new_column_name, dtype, unit, None)
                    original_column_names[new_column_name] = (column_name, tableid)

                # Add differences?
                if differences:
                    #print(tableids)
                    if len(tableids) > 2: raise ValueError("Cannot add differences for more than 2 tables")
                    if not sequences.all_equal(dtypes): raise ValueError("Not the same types")
                    if not sequences.all_equal(units): raise ValueError("Not the same units")
                    dtype = dtypes[0]
                    unit = units[0]
                    difference_column_name = column_name + " difference"
                    merged.add_column_info(difference_column_name, dtype, unit, "difference")
                    new_column_names[difference_column_name] = tableids

                # Add relative differences?
                if rel_differences:
                    #print(tableids)
                    if len(tableids) > 2: raise ValueError("Cannot add relative differences for more than 2 tables")
                    if not sequences.all_equal(dtypes): raise ValueError("Not the same types")
                    if not sequences.all_equal(units): raise ValueError("Not the same units")
                    dtype = dtypes[0]
                    unit = units[0]
                    if percentual: rel_difference_column_name = column_name + " percentual difference"
                    else: rel_difference_column_name = column_name + " relative difference"
                    merged.add_column_info(rel_difference_column_name, dtype, unit, "relative difference")
                    new_column_names[rel_difference_column_name] = tableids

    # All columns are added
    merged._setup()
    column_names = merged.column_names

    # Loop over all reference values
    for value in reference_values:

        # Find index of this value for each table
        indices = []
        for table in tables:
            index = find_index(table, value, column_name=reference_column_name)
            indices.append(index)

        #print(value, indices)

        values = [value]

        # Loop over the column names
        for column_name in column_names:
            if column_name == reference_column_name: continue

            #tableids = all_column_tableids[column_name]

            if column_name in equal_column_names:
                #print(column_name, indices[0])
                #value = tables[0][column_name][indices[0]]
                values_tables = [tables[tableid][column_name][indices[tableid]] for tableid in range(len(tables)) if column_name in tables[tableid].colnames and indices[tableid] is not None]
                #print(column_name, values_tables)
                value = sequences.get_first_not_none_value(values_tables)
                #print(value)

            elif column_name in original_column_names:

                original_column_name, tableid = original_column_names[column_name]
                if indices[tableid] is not None:
                    #print(indices[tableid])
                    value = tables[tableid][original_column_name][indices[tableid]]
                else: value = None

            elif column_name in unique_column_names:

                tableid = unique_column_names[column_name]
                if indices[tableid] is not None:
                    #print(indices[tableid])
                    value = tables[tableid][column_name][indices[tableid]]
                else: value = None

            elif column_name in new_column_names:

                tableids = new_column_names[column_name]
                id_a = tableids[0]
                id_b = tableids[1]
                index_a = indices[id_a]
                index_b = indices[id_b]
                #print(index_a, index_b)

                if index_a is None or index_b is None: value = None
                else:

                    if column_name.endswith("relative difference"):

                        original_column_name = column_name.split(" relative difference")[0]
                        value_a = tables[id_a][original_column_name][index_a]
                        value_b = tables[id_b][original_column_name][index_b]
                        value = abs(value_a - value_b) / value_a

                    elif column_name.endswith("percentual difference"):

                        original_column_name = column_name.split(" percentual difference")[0]
                        value_a = tables[id_a][original_column_name][index_a]
                        value_b = tables[id_b][original_column_name][index_b]
                        value = abs(value_a - value_b) / value_a * 100.

                    elif column_name.endswith("difference"):

                        original_column_name = column_name.split(" difference")[0]
                        value_a = tables[id_a][original_column_name][index_a]
                        value_b = tables[id_b][original_column_name][index_b]
                        value = abs(value_a - value_b)

                    else: raise ValueError("Column name not recognized")

            else: raise RuntimeError("Something went wrong: " + column_name)

            #print(value)

            # Add the value
            values.append(value)

        # Add row
        merged.add_row(values)

    # Return the table
    return merged

# -----------------------------------------------------------------

delattr(Table, "__copy__")
delattr(Table, "__deepcopy__")

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

        #print("SMARTTABLE")
        #print(args)
        #print(kwargs)
        #print("")

        # Always used masked tables
        kwargs["masked"] = True

        #if "names" in kwargs: kwargs.pop("names")
        if "units" in kwargs: kwargs.pop("units")
        #if "dtypes" in kwargs: kwargs["dtype"] = kwargs.pop("dtypes")
        #print(kwargs)

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

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    @property
    def has_column_info(self):
        return len(self.column_info) > 0

    # -----------------------------------------------------------------

    @classmethod
    def from_columns(cls, *columns, **kwargs):

        """
        This function ...
        :param columns:
        :param kwargs:
        :return:
        """

        # Import tostr function
        from ..tools.stringify import tostr

        # Get number of columns
        ncolumns = len(columns)
        nrows = len(columns[0])

        # Get options
        names = kwargs.get("names", None) # not pop, so appears in constructor (of e.g. Relation) as well
        units = kwargs.get("units", None) # not pop, so appears in constructor (of e.g. Relation) as well
        dtypes = kwargs.pop("dtypes", None) ### not pop, so appears in constructor (of e.g. Relation) as well
        descriptions = kwargs.pop("descriptions", None)
        as_columns = kwargs.pop("as_columns", False)
        meta = kwargs.pop("meta", {})
        if meta is None: meta = {}
        #print("names", names)

        # Get tostr kwargs
        tostr_kwargs = kwargs.pop("tostr_kwargs", {})

        # Create the table
        #print(kwargs)
        table = cls(**kwargs) # WILL ALREADY CREATE THE COLUMN INFO

        has_info = table.has_column_info
        #print(type(table).__name__)
        #print(table.column_info)

        # Get the names if necessary
        if names is None:
            if has_info: names = [info[0] for info in table.column_info]
            elif hasattr(cls, "_column_info"): names = cls._column_info.keys()
            else: names = ["col" + str(index) for index in range(ncolumns)]

        #print("names:", names)

        # The column units
        #column_units = dict()

        # Initialize list for the column names for which we have to convert the values to strings
        to_string = []
        #to_split = []

        # Loop over the columns
        for index, name in enumerate(names):

            # Determine dtype
            if dtypes is not None and dtypes[index] is not None: dtype = dtypes[index]
            elif as_columns:
                coltype = columns[index].dtype.name
                column_type = dtype_name_to_column_type(coltype)
                dtype = column_type_to_builtin(column_type)
            else:
                #print(columns[index])
                ptype = get_common_property_type(columns[index])
                dtype, needs_tostring = property_type_to_builtin(ptype, return_tostring=True)
                if needs_tostring: to_string.append(name)

            # Determine the unit
            unit = None
            if units is not None:
                if types.is_sequence_or_tuple(units): unit = units[index]
                elif types.is_dictionary(units):
                    if name in units: unit = units[name]
                    else: unit = None
                else: raise ValueError("Invalid type for 'units': must be list or dictionary")

            # Determine unit from quantities
            if unit is None: unit = get_common_unit(columns[index])

            # Set column unit?
            #self.column_units[name] = unit

            # Determine the description
            if descriptions is not None and descriptions[index] is not None: description = descriptions[index]
            else: description = "no description"

            # Add the column
            #print(unit)
            if not has_info: table.add_column_info(name, dtype, unit, description)

        # Set None string
        if "none_string" not in tostr_kwargs: tostr_kwargs["none_string"] = string_column_none_default

        # Setup the table
        table._setup()
        #print(table.colnames)

        # Add data as columns to the table? -> tostring columns not possible
        if as_columns:

            if len(to_string) > 0: raise ValueError("Cannot add data as columns if some columns have complex types")
            #if len(to_split) > 0: raise ValueError("Cannot add data as columns if some column objects have to be splitted into builtin types")

            # Loop over the columns
            table.remove_all_columns()
            for j, name in enumerate(names):

                # Create column
                col = MaskedColumn(data=columns[j], name=name, dtype=table.get_column_dtype(name), unit=table.get_column_unit(name), copy=False)
                table.add_column(col)

        # Add data as rows
        else:

            # Add the rows, using the Table add_row implementation directly
            # because the add_row function may be prohibited in the actual class (because of lazy features)
            for i in range(nrows):

                # row = []
                # for j, name in enumerate(names):
                #     unit = table.get_column_unit(name)
                #     value = columns[j][i]
                #     if name in to_string: value = tostr(value, **tostr_kwargs)
                #     elif value is None: pass
                #     elif unit is not None and hasattr(value, "unit"): value = value.to(unit).value
                #     else: pass
                #     row.append(value)

                #row = [column[i] for column in columns]
                row = []
                for j, name in enumerate(names):
                    value = columns[j][i]
                    if name in to_string: value = tostr(value, **tostr_kwargs)
                    row.append(value)

                # Add the row
                #super(SmartTable, table).add_row(row)
                #print(row)
                #print([table.column_unit(colname) for colname in table.colnames])
                #print("ROW", row)
                SmartTable.add_row(table, row)

        # Set meta info
        for key in meta: table.meta[key] = meta[key]

        # Return the table
        return table

    # -----------------------------------------------------------------

    @classmethod
    def from_dictionary(cls, dictionary, key_label="Property", value_label="Value", tostr_kwargs=None,
                        key_description="property name", value_description="property value"):

        """
        This function ...
        :param dictionary:
        :param key_label:
        :param value_label:
        :param tostr_kwargs:
        :param key_description:
        :param value_description:
        :return:
        """

        # Import tostr function
        from ..tools.stringify import tostr

        if tostr_kwargs is None: tostr_kwargs = dict()

        # Create the table
        table = cls()

        # Add the column info
        table.add_column_info(key_label, str, None, key_description)
        table.add_column_info(value_label, str, None, value_description)
        table._setup()

        for label in dictionary:

            value = dictionary[label]
            value_string = tostr(value, **tostr_kwargs)

            values = [label, value_string]
            table.add_row(values)

        # Return the table
        return table

    # -----------------------------------------------------------------

    @classmethod
    def from_dictionaries(cls, *dictionaries, **kwargs):

        """
        This function ...
        :param dictionaries:
        :param kwargs:
        :return:
        """

        from ..tools.stringify import stringify

        # Get list of lists of property names
        property_names = [dictionary.keys() for dictionary in dictionaries]

        # Determine the column names, types and descriptions
        column_names = sequences.union(*property_names)

        # Sort column names?
        first = kwargs.pop("first", None)
        last = kwargs.pop("last", None)
        if first is not None or last is not None: column_names = sequences.sort_with_first_last(column_names, first=first, last=last)

        # Ignore all None?
        ignore_none = kwargs.pop("ignore_none", False)
        remove_columns = []

        # Get the column types, units and descriptions
        prop_types = dict()
        prop_units = dict()
        prop_descriptions = dict()
        for name in column_names:

            # Get the types, units and descriptions
            #types = [composite.type_for_property(name) for composite in composites]
            #units = [composite.unit_for_property(name) for composite in composites]
            #descriptions = [composite.description_for_property(name) for composite in composites]

            types = []
            units = []
            descriptions = []

            for dictionary in dictionaries:

                # Get value
                if name in dictionary:

                    value = dictionary[name]

                    # Check whether there is a unit
                    if hasattr(value, "unit"): unit = value.unit
                    else: unit = None

                    # Get the type
                    dtype, string = stringify(value)

                # Not in this dictionary
                else: unit = dtype = None

                # Add the type and
                types.append(dtype)
                units.append(unit)

            # Determine column type, unit and description
            if sequences.all_equal_to(types, 'None') or sequences.all_none(types):
                if ignore_none: remove_columns.append(name)
                column_type = 'None'
            else: column_type = sequences.get_all_equal_value(types, ignore_none=True, ignore='None')

            # Determine column unit
            column_unit = sequences.get_first_not_none_value(units)

            # Determine column description
            #column_description = sequences.get_first_not_none_value(descriptions)
            column_description = None

            # Set type, unit and description
            prop_types[name] = column_type
            prop_units[name] = column_unit
            prop_descriptions[name] = column_description

        # Remove columns
        if len(remove_columns) > 0: column_names = sequences.removed(column_names, remove_columns)

        # Create and return
        return cls.from_properties(column_names, prop_types, prop_units, prop_descriptions, dictionaries, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_objects(cls, objects):

        """
        This function creates a table from objects of the SAME TYPE (same attributes)
        :param objects:
        :return:
        """

        # Create list of attribute dictionaries
        dictionaries = [obj.__dict__ for obj in objects]

        # Create the table and return
        return cls.from_dictionaries(*dictionaries)

    # -----------------------------------------------------------------

    @classmethod
    def from_composite(cls, composite, key_label="Property", value_label="Value", tostr_kwargs=None,
                       key_description="property name", value_description="property value"):

        """
        This function ...
        :param composite:
        :param key_label:
        :param value_label:
        :param tostr_kwargs:
        :param key_description:
        :param value_description:
        :return:
        """

        # Import tostr function
        from ..tools.stringify import tostr

        if tostr_kwargs is None: tostr_kwargs = dict()

        # Create the table
        table = cls()

        #column_types = [str, str]
        #column_units = [None, None]
        #keys_values = composite.as_tuples()
        #column_names = [key_label, value_label]

        # Add the column info
        table.add_column_info(key_label, str, None, key_description)
        table.add_column_info(value_label, str, None, value_description)
        table._setup()

        #print(table.column_info)

        for key in composite:

            value = composite[key]
            value_string = tostr(value, **tostr_kwargs)

            values = [key, value_string]
            table.add_row(values)

        # Return the table
        return table

    # -----------------------------------------------------------------

    @classmethod
    def from_composites(cls, *composites, **kwargs):

        """
        This function ...
        :param composites:
        :param kwargs:
        :return:
        """

        # Check number
        if len(composites) == 0: raise ValueError("No input is provided")

        # Get list of lists of property names
        property_names = [composite.property_names for composite in composites]

        # Determine the column names, types and descriptions
        column_names = sequences.union(*property_names)

        # Sort column names?
        first = kwargs.pop("first", None)
        last = kwargs.pop("last", None)
        if first is not None or last is not None: column_names = sequences.sort_with_first_last(column_names, first=first, last=last)

        # Get the column types, units and descriptions
        prop_types = dict()
        prop_units = dict()
        prop_descriptions = dict()
        for name in column_names:

            # Get the types, units and descriptions
            types = [composite.type_for_property(name) for composite in composites]
            units = [composite.unit_for_property(name) for composite in composites]
            descriptions = [composite.description_for_property(name) for composite in composites]

            # Determine column type, unit and description
            if sequences.all_equal_to(types, 'None') or sequences.all_none(types): column_type = 'None'
            else: column_type = sequences.get_all_equal_value(types, ignore_none=True, ignore='None')

            # Determine column unit
            column_unit = sequences.get_first_not_none_value(units)

            # Determine column description
            column_description = sequences.get_first_not_none_value(descriptions)

            # Set type, unit and description
            prop_types[name] = column_type
            prop_units[name] = column_unit
            prop_descriptions[name] = column_description

        # Create and return
        kwargs["attr"] = True
        return cls.from_properties(column_names, prop_types, prop_units, prop_descriptions, composites, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_properties(cls, property_names, property_types, property_units, property_descriptions, objects, **kwargs):

        """
        This function ...
        :param property_names:
        :param objects:
        :param property_types:
        :param property_units:
        :param property_descriptions:
        :param kwargs:
        :return:
        """

        # Import tostr function
        from ..tools.stringify import tostr

        # Get names
        labels = kwargs.pop("labels", None)
        label = kwargs.pop("label", "-")

        # Get tostr kwargs
        tostr_kwargs = kwargs.pop("tostr_kwargs", {})

        # Get flag
        attr = kwargs.pop("attributes", False)

        # Add the label column name
        if labels is not None: column_names = [label] + property_names #sequences.prepend(column_names, label)
        else:
            column_names = property_names
            labels = [None] * len(objects) # labels was None

        # Create the table
        table = cls()

        # Create lists to contain the column types and units
        column_types = []
        column_units = []

        # COLUMNS FOR WHICH THE VALUE HAS TO BE CONVERTED TO STRING
        # OR SPLITTED INTO MULTIPLE TYPES
        to_string = []
        to_split = []
        actual_column_names = OrderedDict()

        # Make the columns
        for name in column_names:

            if name == label:

                real_type = str
                column_unit = None
                column_description = label

                # Add the column info
                table.add_column_info(name, real_type, column_unit, column_description)
                actual_column_names[name] = (name, None, None, None)

            else:

                # Get type, unit and description
                column_type = property_types[name]
                column_unit = property_units[name]
                column_description = property_descriptions[name]

                # Add column type and unit to lists
                column_types.append(column_type)
                column_units.append(column_unit)

                # Get builtin type(s) for the column
                if composed_of_multiple_builtins(column_type):

                    all_items, itemspec = property_type_to_builtins(column_type)
                    if all_items:
                        keys_lists = [getattr(obj, name).keys() for obj in objects if hasattr(obj, name)] if attr else [obj[name].keys() for obj in objects if name in obj]
                        #print(keys_lists)
                        item_keys = sequences.union(*keys_lists)
                        for key in item_keys:
                            item_colname = name + " " + key
                            #print(item_colname)
                            # Get type, unit
                            raise NotImplementedError("Not yet implemented")
                            #table.add_column_info(name, real_type, column_unit, column_description)
                            #actual_column_names.append(item_colname)
                    for itemname, itemtype, item_is_attr in itemspec:
                        #print(itemname, itemtype, item_is_attr)
                        item_colname = name + " " + itemname
                        #print(item_colname)
                        itemdescription = column_description + " " + itemname if column_description is not None else None
                        if item_is_attr: itemvalues = [getattr(getattr(obj, name), itemname) for obj in objects if hasattr(obj, name)] if attr else [getattr(obj[name], itemname) for obj in objects if name in obj]
                        else: itemvalues = [getattr(obj, name)[itemname] for obj in objects if hasattr(obj, name)] if attr else [obj[name][itemname] for obj in objects if name in obj]
                        # Determine item unit
                        itemunits = [value.unit for value in itemvalues if hasattr(value, "unit")]
                        itemunit = sequences.get_first_not_none_value(itemunits)
                        #print(itemname, itemvalues, itemunit)
                        table.add_column_info(item_colname, itemtype, itemunit, itemdescription)
                        actual_column_names[item_colname] = (name, itemname, itemunit, item_is_attr)
                else:

                    #print(column_type)
                    real_type, needs_tostring = property_type_to_builtin(column_type, return_tostring=True)
                    if needs_tostring: to_string.append(name)

                    # Add the column info
                    table.add_column_info(name, real_type, column_unit, column_description)
                    actual_column_names[name] = (name, None, None, None)

        # Set None string
        if "none_string" not in tostr_kwargs: tostr_kwargs["none_string"] = string_column_none_default

        # Add the rows
        for composite_label, obj in zip(labels, objects):

            values = []

            # Fill the row
            #for name, dtype, unit in zip(column_names, column_types, column_units):
            #for name in column_names:
            for name in actual_column_names:

                pname, itemname, itemunit, itemattr = actual_column_names[name]
                #print(name, pname, punit, pisattr)

                if name == label: value = composite_label

                else:

                    # Properties are attributes of the objects
                    if attr:

                        # Get the value
                        if hasattr(obj, pname):
                            value = getattr(obj, pname)
                            if pname in to_string: value = tostr(value, **tostr_kwargs)
                        else: value = None

                    # Properties are items of the objects
                    else:

                        if pname in obj:
                            value = obj[pname]
                            if pname in to_string: value = tostr(value, **tostr_kwargs)
                        else: value = None

                    if name != pname:
                        if value is not None:
                            if itemattr: value = getattr(value, itemname)
                            else: value = value[itemname]

                # Add the value
                #print(value)
                values.append(value)

            # Add the row: unit conversion is done here
            table.add_row(values)

        #print(table.colnames)
        #print(table)

        # Return the table
        return table

    # -----------------------------------------------------------------

    def remove_other_columns(self, names):
        remove_names = sequences.get_other(self.column_names, names)
        self.remove_columns(remove_names)

    # -----------------------------------------------------------------

    def remove_all_columns(self):
        self.remove_columns(self.column_names)

    # -----------------------------------------------------------------

    def remove_all_rows(self):
        self.remove_rows(range(self.nrows))

    # -----------------------------------------------------------------

    def get_column_index(self, column_name):
        for index in range(len(self.column_info)):
            if self.column_info[index][0] == column_name: return index
        return None

    # -----------------------------------------------------------------

    def get_column_dtype(self, column_name):
        index = self.get_column_index(column_name)
        if index is None: raise ValueError("No column '" + column_name + "'")
        return self.column_info[index][1]

    # -----------------------------------------------------------------

    def get_column_array_dtype(self, column_name):
        return self[column_name].dtype

    # -----------------------------------------------------------------

    def get_column_array_dtype_string(self, column_name):
        return str(self.get_column_array_dtype(column_name))

    # -----------------------------------------------------------------

    def is_string_column(self, column_name):
        return self.get_column_array_dtype_string(column_name).startswith("|S")

    # -----------------------------------------------------------------

    def get_string_column_size(self, column_name):
        return int(self.get_column_array_dtype_string(column_name).split("S")[1])

    # -----------------------------------------------------------------

    def get_column_unit(self, column_name):
        index = self.get_column_index(column_name)
        if index is None: raise ValueError("No column '" + column_name + "'")
        return self.column_info[index][2]

    # -----------------------------------------------------------------

    def get_column_description(self, column_name):
        index = self.get_column_index(column_name)
        if index is None: raise ValueError("No column '" + column_name + "'")
        return self.column_info[index][3]

    # -----------------------------------------------------------------

    def add_column_info(self, name, dtype, unit, description, index=None):

        """
        This function ...
        :param name:
        :param dtype:
        :param unit:
        :param description:
        :param index:
        :return:
        """

        #print("adding column info", name, dtype, unit, description, index)
        if types.is_string_type(unit): unit = u(unit)
        if index is not None: self.column_info.insert(index, (name, dtype, unit, description))
        else: self.column_info.append((name, dtype, unit, description))

    # -----------------------------------------------------------------

    @property
    def column_info_names(self):
        return [info[0] for info in self.column_info]

    # -----------------------------------------------------------------

    def add_all_column_info(self, info_dict):

        """
        This function ...
        :param info_dict:
        :return:
        """

        # Loop over the columns
        for column_name in info_dict:

            info = info_dict[column_name]
            self.add_column_info(column_name, info[0], info[1], info[2])

    # -----------------------------------------------------------------

    def has_column_unit(self, column_name):

        """
        This function ...
        :param column_name:
        :return:
        """

        return self[column_name].unit is not None

    # -----------------------------------------------------------------

    def column_unit(self, column_name):

        """
        This function ...
        :param column_name:
        :return:
        """

        if not self.has_column_unit(column_name): return None

        # Construct unit
        if column_name in self.meta["density"]: density = True
        else: density = False
        if column_name in self.meta["brightness"]: brightness = True
        else: brightness = False
        return u(self[column_name].unit, density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def is_photometric_column(self, column_name):
        return self.has_column_unit(column_name) and isinstance(self.column_unit(column_name), PhotometricUnit)

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

    def _setup(self):

        """
        This function ...
        :return:
        """

        #print("Performing setup of the table ...")

        # Setup has already been called
        if len(self.colnames) != 0: return

        #print(self.column_info)

        # Create the table names and types lists
        for entry in self.column_info:

            name = entry[0]
            dtype = entry[1]
            unit = entry[2]
            description = entry[3]

            data = []

            # Add column
            #print(name, dtype)
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
        if len(self.colnames) == 0: self._setup()

        # Call the implementation of the base class
        return super(SmartTable, self).__getitem__(item)

    # -----------------------------------------------------------------

    @classmethod
    def from_remote_file(cls, path, remote, format=None):

        """
        This function ...
        :param path:
        :param remote:
        :param format:
        :return:
        """

        # Check if file exists
        if not remote.is_file(path): raise IOError("The file '" + path + "' does not exist")

        # Get filename
        filename = fs.strip_extension(fs.name(path))

        # Get the lines
        lines = remote.get_lines(path)

        # Get the table from the lines and return
        return cls.from_lines(lines, format=format, table_name=filename)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, format=None, method="lines"):

        """
        This function ...
        :param path:
        :param format:
        :param method:
        :return:
        """

        # Check the path
        if not fs.is_file(path): raise IOError("The file '" + path + "' does not exist")

        # FITS
        if path.endswith(".fits"): return cls.from_fits_file(path)

        # Other: assume ascii
        else: return cls.from_ascii_file(path, format=format, method=method)

    # -----------------------------------------------------------------

    @classmethod
    def from_fits_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get filename
        filename = fs.strip_extension(fs.name(path))

        # Open table
        table = super(SmartTable, cls).read(path, format="fits")

        # Initialize
        initialize_table(table, table_name=filename)

        # Re-order the columns
        reorder_columns(table)

        # Set the path
        table.path = path

        # Return the table
        return table

    # -----------------------------------------------------------------

    @classmethod
    def from_ascii_file(cls, path, format=None, method="lines"):

        """
        This function ...
        :param path:
        :param format:
        :param method:
        :return:
        """

        # Get filename
        filename = fs.strip_extension(fs.name(path))

        # Guess the format
        if format is None:
            first_line = fs.get_first_line(path)
            if "ECSV" in first_line: format = "ecsv"
            elif "PTS data format" in first_line: format = "pts"

        # More? ecsv as well?
        allowed_numpy_pands_formats = ["pts", "csv"]

        # Read as lines and parse each line individually
        if method == "lines":

            # Read lines: NO, astropy doesn't like generators: 'Input "table" must be a string (filename or data) or an iterable')
            #lines = fs.read_lines(path)
            lines = fs.get_lines(path)

            # Create table from the lines
            table = cls.from_lines(lines, format=format, table_name=filename)

        # Read using Pandas
        elif method == "pandas":

            # Check
            if format not in allowed_numpy_pands_formats: raise IOError("Reading through Pandas is currently only supported for PTS style tables and CSV ASCII tables")

            # Read the header
            column_names, column_types, column_units, meta = parse_header_file(path, format)

            import pandas as pd
            df = pd.read_csv(path, sep=" ", comment="#", header=None)
            ncolumns = len(df.columns)
            cols = [df[index].values for index in range(ncolumns)]
            table = cls.from_columns(*cols, as_columns=True, names=column_names, dtypes=column_types, units=column_units, meta=meta)

        # Read using NumPy
        elif method == "numpy":

            # Check
            if format not in allowed_numpy_pands_formats: raise IOError("Reading through NumPy is currently only supported for PTS style tables and CSV ASCII tables")

            # Read the header
            column_names, column_types, column_units, meta = parse_header_file(path, format)

            columns = np.loadtxt(path, unpack=True)
            ncolumns = len(columns)
            table = cls.from_columns(*columns, as_columns=True, names=column_names, dtypes=column_types, units=column_units, meta=meta)

        # Invalid
        else: raise ValueError("Invalid option for 'method'")

        # Set the path
        table.path = path

        # Return the table
        return table

    # -----------------------------------------------------------------

    @classmethod
    def from_lines(cls, lines, format=None, always_check_types=False, table_name=None):

        """
        This function ...
        :param lines:
        :param format:
        :param always_check_types:
        :param table_name:
        :return:
        """

        from ..tools import strings

        # If format is undefined, read the format from the first line
        if format is None:
            #lines = list(lines) # in case it is an iterator (NOW WE ALWAYS PASS A LIST)
            first_line = lines[0]
            if "ECSV" in first_line: format = "ecsv"
            elif "PTS data format" in first_line: format = "pts"

        # fill_values = [('--', '0')]
        fill_values = [('--', '0')]

        # PTS format
        if format == "pts":

            # Initialize the data lines and header lines
            data_lines = []
            header = []

            # Loop over the lines of the file
            for line in lines:

                # Header line
                if line.startswith("#"): header.append(line[2:])
                else: data_lines.append(line)

            # Put last header line (colum names) as first data line (and remove it from the header)
            sequences.prepend(data_lines, "# " + header[-1])
            header = header[:-1]

            # Search for density and brightness, set meta info
            index = -1
            density = None
            brightness = None
            for index, line in enumerate(header):

                if line.startswith("PTS data format"): continue

                if line.startswith("density:"):

                    density_string = line.split("density: ")[1]
                    string = "[" + ",".join('"' + s + '"' for s in density_string.split(",")) + "]"
                    density = eval(string)

                elif line.startswith("brightness:"):

                    brightness_string = line.split("brightness: ")[1]
                    string = "[" + ",".join('"' + s + '"' for s in brightness_string.split(",")) + "]"
                    brightness = eval(string)

                else: break

            # Get the META info
            from ..tools import parsing
            meta = OrderedDict()
            old_header = header[:]
            header = []
            for line in old_header:
                if line.startswith("META"):
                    spec = line.split("META: ")[1]
                    key = spec.split(" [")[0]
                    ptype = spec.split("[")[1].split("]")[0]
                    if ptype == "None": meta[key] = None
                    else:
                        string = spec.split("] ")[1]
                        parsing_function = getattr(parsing, ptype)
                        value = parsing_function(string)
                        meta[key] = value
                else: header.append(line)
            #print(meta)

            # Contains type line
            if index == len(header) - 2:

                # Get types
                types_string = header[-2]
                column_type_strings = types_string.split()
                column_types = []
                for type_string in column_type_strings:
                    dtype = column_type_to_builtin(type_string)
                    column_types.append(dtype)

            # Doesn't contain type line
            elif index == len(header) - 1: column_types = column_type_strings = None

            # Invalid
            else: raise IOError("Something is wrong with the file")

            ndata_lines = len(data_lines)
            nrows = ndata_lines - 1

            # Call the constructor from Astropy, to read in plain ascii format
            if nrows == 0:
                if column_types is None: raise IOError("The table file doesn't contain the column types and neither does it have any rows")
                colnames_string = data_lines[0][2:]
                column_names = strings.split_except_within_double_quotes(colnames_string, add_quotes=False)
                table = cls(data=None, names=column_names, masked=True, dtype=column_types, copy=True)
            else:
                # Was just to check:
                #from ..tools import strings
                #colnames_string = data_lines[0][2:]
                #column_names = strings.split_except_within_double_quotes(colnames_string, add_quotes=False)
                #first_row_strings = strings.split_except_within_double_quotes(data_lines[1])
                #print(len(column_names), len(first_row_strings))
                #for line in data_lines: print(line)
                table = super(SmartTable, cls).read(data_lines, format="ascii.commented_header")

            # Set density and brightness meta
            if density is not None: table.meta["density"] = density
            if brightness is not None: table.meta["brightness"] = brightness

            # Set other metadata
            for key in meta: table.meta[key] = meta[key]

            # FIX BOOLEAN COLUMNS
            # and also check column types
            to_boolean = []
            if column_type_strings is not None:

                for k, column_name in enumerate(table.colnames):
                    column_type_string = column_type_strings[k]
                    actual_column_type = dtype_name_to_column_type(table[column_name].dtype.name)
                    #print(column_type_string, actual_column_type)
                    if actual_column_type == column_type_string: continue
                    elif column_type_string == "boolean" and actual_column_type == "string": to_boolean.append(column_name)
                    else:
                        #raise RuntimeError("Something went wrong: column '" + column_name + "' with type '" + column_type_string + "' is parsed as '" + actual_column_type + "' column")
                        warnings.warn("Column '" + column_name + "' with type '" + column_type_string + "' is parsed as '" + actual_column_type + "' column: fixing ...")
                        real_type = column_type_to_builtin(column_type_string)
                        table[column_name] = table[column_name].astype(real_type)

            #else:
            # NOW ALWAYS CHECK AGAIN?
            if column_type_strings is None or always_check_types:
                for column_name in table.colnames:
                    #values = list(table[column_name]) # THIS IS VERY SLOW FOR LARGE TABLES
                    values = table[column_name] # DOES THIS WORK AS WELL?
                    if len(values) == 0: continue
                    #value_strings = [str(value) for value in values] # actually not necessary
                    #print(column_name, values)
                    if not sequences.all_strings(values, ignore_instance=np.ma.core.MaskedConstant): continue
                    #print(column_name)
                    #for value in values: print(value, type(value))
                    if sequences.all_in(values, ["True", "False"], ignore_instance=np.ma.core.MaskedConstant):
                        if column_type_strings is not None: warnings.warn("Column type strings were specified but still I think that column '" + column_name + "' actually represents boolean values and not strings: loading as boolean column ...")
                        to_boolean.append(column_name)

            #print(index, header[index], len(header))

            # Set units
            unit_string = header[-1]
            #unit_strings = unit_string.split()
            unit_strings = strings.split_except_within_round_brackets_and_double_quotes(unit_string)
            #print(unit_string)
            #print(unit_strings)
            #print(table.colnames)
            assert len(unit_strings) == len(table.colnames)
            #print()
            #print(unit_strings)
            #print(table.colnames)
            #print(len(unit_strings), len(table.colnames))
            #print(unit_strings)
            for unit_string, colname in zip(unit_strings, table.colnames):
                if unit_string == '""': continue
                table[colname].unit = unit_string

            #print(to_boolean)
            # DO THIS AT THE END BECAUSE OTHERWISE UNITS ASSIGNED TO THE WRONG COLUMNS
            # Loop over the columns to convert to booleans
            for column_name in to_boolean:
                booleans = []
                masked = []
                for index in range(len(table)):

                    # Masked?
                    if table[column_name].mask[index]:
                        # boolean = None # doesn't work
                        # boolean = np.ma.core.MaskedConstant
                        boolean = False # just a default value
                        masked.append(index)
                    else:
                        value = table[column_name][index]
                        #print(value, list(value))
                        #print(value, type(value))
                        if types.is_string_type(value): boolean = eval(value)
                        else: boolean = bool(value)
                    booleans.append(boolean)
                # remove original column
                #print(table.colnames.index(column_name))
                #table.remove_column(column_name)
                #table[column_name] = booleans
                #print(table.colnames.index(column_name))
                table.replace_column(column_name, booleans)  # TO AVOID CHANGING THE ORDER OF COLUMNS
                #print(table.colnames.index(column_name))
                for index in masked: table[column_name].mask[index] = True

        # ECSV format (with masks and units in the meta info)
        elif format == "ecsv":

            # Read
            table = super(SmartTable, cls).read(lines, fill_values=fill_values, format="ascii.ecsv")

            # Set masks
            set_table_masks(table)

            # Fix
            fix_column_types(table)

            # Fix None strings
            fix_nones(table)

        # Write the table in the desired format (by Astropy)
        elif format == "csv": table = super(SmartTable, cls).read(lines, fill_values=fill_values, format="ascii.csv")

        # HTML
        elif format == "html": table = super(SmartTable, cls).read(lines, fill_values=fill_values, format="ascii.html")

        # Latex
        elif format == "latex": table = super(SmartTable, cls).read(lines, fill_values=fill_values, format="ascii.latex")

        # Unknown
        elif format is None:

            # Try guessing from Astropy
            try: table = super(SmartTable, cls).read(lines, fill_values=fill_values)
            except IORegistryError:
                # Try ASCII
                table = super(SmartTable, cls).read(lines, fill_values=fill_values, format="ascii")

        # All other
        else: table = super(SmartTable, cls).read(lines, fill_values=fill_values, format=format)

        # Initialize
        initialize_table(table, table_name=table_name)

        # Re-order the columns
        reorder_columns(table)

        # Return the table
        return table

    # -----------------------------------------------------------------

    def _resize_string_columns(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        # Initialize dictionary for the new sizes of the columns
        #new_sizes = dict()

        # Loop over the columns
        for index, colname in enumerate(self.column_names):

            # Skip non-string columns
            if not self.is_string_column(colname): continue

            # Resize if necessary
            #print(colname, values[index])
            self.resize_string_column_for_string(colname, values[index])

        # Resize columns
        #for colname in new_sizes: self.resize_string_column(colname, new_sizes[colname])

    # -----------------------------------------------------------------

    def resize_string_column_for_string(self, colname, string):

        """
        This function ...
        :param colname:
        :param string:
        :return:
        """

        # Get current column string length
        current_string_length = self.get_string_column_size(colname)

        # Get new string length
        if string is None: new_string_length = 0
        else: new_string_length = len(string)

        # Doesn't need resize?
        if new_string_length <= current_string_length: return

        # Resize
        self.resize_string_column(colname, new_string_length)

    # -----------------------------------------------------------------

    def resize_string_column(self, colname, size):

        """
        This function ...
        :param colname:
        :param size:
        :return:
        """

        # Create new version of the column
        resized = self[colname].astype("S" + str(size))

        # Replace the column
        self.replace_column(colname, resized)

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

    def _strip_units(self, values, conversion_info=None):

        """
        This function ...
        :param values:
        :param conversion_info:
        :return:
        """

        from ..units.unit import get_converted_value

        scalar_values = []

        # Loop over the column values
        for i, value in enumerate(values):

            # Get the column name
            colname = self.column_info[i][0]

            # If this value has a unit, we have to make sure it is converted into the proper column unit
            if hasattr(value, "unit"):

                column_unit = self.column_unit(colname)
                assert column_unit is not None

                # Get the conversion info for this column
                if conversion_info is not None and colname in conversion_info: conv_info = conversion_info[colname]
                else: conv_info = dict()

                # Get converted scalar value
                scalar_value = get_converted_value(value, column_unit, conversion_info=conv_info)

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

            #if isinstance(value, list):
            if types.is_sequence_or_tuple(value):
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

    def get_value(self, colname, index, add_unit=True, unit=None, conversion_info=None):

        """
        This function ...
        :param colname:
        :param index:
        :param add_unit:
        :param unit:
        :param conversion_info:
        :return:
        """

        # Masked? -> return None
        if self[colname].mask[index]: return None

        # Get raw, scalar value
        value = self[colname][index]

        # Add unit
        if self.has_column_unit(colname): value = value * self.column_unit(colname)

        # Convert?
        if unit is not None:

            if conversion_info is None: conversion_info = {}
            if not self.has_column_unit(colname): raise ValueError("Column '" + colname + "' does not have a unit")

            # Convert
            if self.is_photometric_column(colname): value = value.to(unit, **conversion_info)
            else: value = value.to(unit)

        # Return the value
        if add_unit: return value
        elif self.has_column_unit(colname): return value.value
        else: return value

    # -----------------------------------------------------------------

    def get_column_values(self, colname, add_unit=True, unit=None):

        """
        This function ...
        :param colname:
        :param add_unit:
        :param unit:
        :return:
        """

        values = []
        for index in range(self.nrows): values.append(self.get_value(colname, index, add_unit=add_unit, unit=unit))
        return values

    # -----------------------------------------------------------------

    def get_array(self, colname):
        return np.asarray(self[colname])

    # -----------------------------------------------------------------

    def get_unit(self, colname):
        return self[colname].unit

    # -----------------------------------------------------------------

    def get_column_array(self, colname, unit=None, masked=True):

        """
        This function ...
        :param colname:
        :param unit:
        :param masked:
        :return:
        """

        # Column unit
        if self.has_column_unit(colname):

            # Check if unit is defined
            if unit is None: raise ValueError("Unit has to be defined")

            # Get conversion factor
            #conversion_factor = self.get_column_unit(colname).conversion_factor(unit)
            from ..units.unit import get_conversion_factor
            conversion_factor = get_conversion_factor(self.get_column_unit(colname), unit, parse=False)

            # Return
            if masked: return self[colname].data * conversion_factor
            else: return np.asarray(self[colname].data) * conversion_factor

        # No unit
        elif unit is not None: raise ValueError("Unit of column is not defined")
        else:
            if masked: return self[colname].data # masked array
            else: return np.asarray(self[colname].data) # asarray to make regular array

    # -----------------------------------------------------------------

    def get_column_mask(self, colname, invert=False):

        """
        This function ...
        :param colname:
        :param invert:
        :return:
        """

        if invert: return np.logical_not(self[colname].mask)
        else: return self[colname].mask

    # -----------------------------------------------------------------

    def get_column_nmasked(self, colname):
        return np.sum(self.get_column_mask(colname))

    # -----------------------------------------------------------------

    def get_column_nnotmasked(self, colname):
        return np.sum(self.get_column_mask(colname, invert=True))

    # -----------------------------------------------------------------

    def get_column_sum(self, colname, unit=None, add_unit=True):

        """
        This function ...
        :param colname:
        :param unit:
        :param add_unit:
        :return:
        """

        # Get masked array
        array = self.get_column_array(colname, unit=unit, masked=True)

        # Get sum
        value = np.sum(array.compressed())

        # Return
        if add_unit and unit is not None: return value * unit
        else: return value

    # -----------------------------------------------------------------

    def get_values(self, colnames, index, add_unit=True, as_dict=False):

        """
        This function ...
        :param colnames:
        :param index:
        :param add_unit:
        :param as_dict:
        :return:
        """

        # Initialize dictionary
        values = OrderedDict()

        # Loop over the columns
        for name in colnames: values[name] = self.get_value(name, index, add_unit=add_unit)

        # Return the values
        if as_dict: return values
        else: return values.values()

    # -----------------------------------------------------------------

    def set_value(self, colname, index, value, conversion_info=None, return_previous=False):

        """
        This function ...
        :param colname:
        :param index:
        :param value:
        :param return_previous:
        :param conversion_info:
        :return:
        """

        # Get the current value
        if return_previous: previous = self.get_value(colname, index)
        else: previous = None

        # Value is None?
        if value is None: self[colname].mask[index] = True

        # Column with unit
        elif self.has_column_unit(colname):

            # Set the value in the correct unit
            if conversion_info is None: conversion_info = dict()
            self[colname][index] = value.to(self.column_unit(colname), **conversion_info).value

        # Column without unit: check that value has no unit
        elif hasattr(value, "unit"): raise ValueError("Value has unit but unit of column is undefined")

        # String column
        elif self.is_string_column(colname):
            self.resize_string_column_for_string(colname, value)
            self[colname][index] = value

        # Other column type
        else: self[colname][index] = value

        # Return the previous value
        return previous

    # -----------------------------------------------------------------

    def is_masked_value(self, colname, index):
        return self[colname].mask[index]

    # -----------------------------------------------------------------

    def get_row(self, index, add_units=True, as_list=False, unit=None, add_unit=True):

        """
        This function ...
        :param index:
        :param add_units:
        :param as_list:
        :param unit:
        :param add_unit:
        :return:
        """

        # Check
        if not add_unit and unit is None: raise ValueError("You cannot know in which unit the values are going to be")

        # Initialize
        row = OrderedDict()

        # Loop over the columns
        for name in self.colnames:

            # Get the value
            value = self.get_value(name, index, add_unit=add_units)
            if unit is not None: value = value.to(unit)
            if not add_unit and hasattr(value, "unit"): value = value.value

            # Add the value
            row[name] = value

        # Return the row
        if as_list: return row.values()
        else: return row

    # -----------------------------------------------------------------

    def add_row(self, values, conversion_info=None):

        """
        This function ...
        :param values:
        :param conversion_info:
        :return:
        """

        # Setup if necessary
        if len(self.colnames) == 0: self._setup()

        #print(values, self.column_info)

        # Resize string columns for the new values
        self._resize_string_columns(values)

        # Get cleaned values
        values, mask = self._prepare_row_values(values, conversion_info=conversion_info)

        # Add the row
        super(SmartTable, self).add_row(values, mask=mask)

    # -----------------------------------------------------------------

    def _prepare_row_values(self, values, conversion_info=None):

        """
        This function ...
        :param values:
        :param conversion_info:
        :return:
        """

        # Strip units
        values = self._strip_units(values, conversion_info=conversion_info)

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

        # Return the new values
        return new_values, mask

    # -----------------------------------------------------------------

    def add_col(self, name, values, dtype=None, unit=None, description=None, as_column=False):

        """
        This function ...
        :param name:
        :param values:
        :param dtype:
        :param unit:
        :param description:
        :param as_column: for LARGE columns (big arrays)
        :return:
        """

        # Determine dtype
        to_string = []
        if dtype is None:
            if as_column:
                coltype = values.dtype.name
                column_type = dtype_name_to_column_type(coltype)
                dtype = column_type_to_builtin(column_type)
            else:
                ptype = get_common_property_type(values)
                dtype, needs_tostring = property_type_to_builtin(ptype, return_tostring=True)
                if needs_tostring: to_string.append(name)

        # Determine unit
        if unit is None:
            if as_column: pass
            # Determine unit from quantities
            else: unit = get_common_unit(values)

        # Column info
        if description is None: description = "no description"
        self.add_column_info(name, dtype, unit, description)

        # Directly use values as column?
        if as_column:
            if len(to_string) > 0: raise ValueError("Cannot add data as columns if some columns have complex types")
            column = MaskedColumn(name=name, data=values, dtype=dtype, unit=unit, copy=False)

        # Prepare first
        else:

            # Prepare values
            values, mask = self._prepare_column_values(values, name)

            # Add
            column = MaskedColumn(name=name, data=values, mask=mask, unit=unit)

        # Add to table
        self.add_column(column)

    # -----------------------------------------------------------------

    def _prepare_column_values(self, values, index_or_name, conversion_info=None):

        """
        This function ...
        :param values:
        :param index_or_name:
        :param conversion_info:
        :return:
        """

        from ..units.unit import get_converted_value

        # Get the conversion info for this column
        if conversion_info is None: conversion_info = dict()

        # Get column name
        if types.is_string_type(index_or_name): colname, index = index_or_name, self.column_info_names.index(index_or_name)
        elif types.is_integer_type(index_or_name): colname, index = self.colnames[index_or_name], index_or_name
        else: raise ValueError("Invalid type for 'index_or_name'")

        # Initialize list for new values
        new_values = []
        mask = []
        for value in values:

            # If this value has a unit, we have to make sure it is converted into the proper column unit
            if hasattr(value, "unit"):

                column_unit = self.column_unit(colname)
                assert column_unit is not None

                # Get converted scalar value
                converted_value = get_converted_value(value, column_unit, conversion_info=conversion_info)

            # Value is a sequence (or tuple)
            elif types.is_sequence_or_tuple(value):

                column_type = self.column_info[index][1]
                if len(value) == 0: converted_value = None
                elif column_type == str: converted_value = ",".join(map(str, value))
                else: raise ValueError("Cannot have a list element in the row at the column that is not of string type")

            # Other kind of values
            else: converted_value = value

            if converted_value is None:

                if self.is_string_type(colname): new_values.append("")
                elif self.is_real_type(colname): new_values.append(0.)
                elif self.is_integer_type(colname): new_values.append(0)
                elif self.is_boolean_type(colname): new_values.append(False)
                else: raise ValueError("Unknown column type for '" + colname + "'")
                mask.append(True)

            else:
                new_values.append(converted_value)
                mask.append(False)

        # Return the new values and mask
        return new_values, mask

    # -----------------------------------------------------------------

    def add_row_from_dict(self, dictionary, conversion_info=None):

        """
        This function ...
        :param dictionary:
        :param conversion_info:
        :return:
        """

        # Initialize list for the values
        values = []

        # Loop over the column names
        for column_name in self.column_names:

            if column_name not in dictionary: value = None
            else: value = dictionary[column_name]

            # Add the value
            values.append(value)

        #print(values)

        # Add row
        self.add_row(values, conversion_info=conversion_info)

    # -----------------------------------------------------------------

    def column_type(self, column_name):

        """
        This function ...
        :param column_name: 
        :return: 
        """

        coltype = self[column_name].dtype.name
        return dtype_name_to_column_type(coltype)

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
        #print(first)
        for i in range(len(self[column_name])):
            #print(self[column_name][i])
            if self[column_name][i] != first: return False
        return True

    # -----------------------------------------------------------------

    @property
    def filename(self):

        """
        This function ...
        :return:
        """

        if self.path is None: return None
        return fs.strip_extension(fs.name(self.path))

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

    def saveto(self, path, format="pts"):

        """
        This function ...
        :param path:
        :param format:
        :return:
        """

        # Setup if necessary
        if len(self.colnames) == 0: self._setup()

        # If the file already exists, remove
        if fs.is_file(path): fs.remove_file(path)

        # PTS format
        if format == "pts": self.saveto_pts(path)

        # ECSV format (with masks and units in the meta info)
        elif format == "ecsv": self.saveto_ecsv(path)

        # Write the table in the desired format (by Astropy)
        elif format == "csv": self.write(path, format="ascii.csv")

        # HTML
        elif format == "html": self.write(path, format="ascii.html")

        # Latex
        elif format == "latex": self.write(path, format="ascii.latex")

        # All other
        else: self.write(path, format=format)

        # Set the path
        self.path = path

    # -----------------------------------------------------------------

    def saveto_pts(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Import tostr function
        from ..tools.stringify import tostr, stringify

        # Create string buffer
        import StringIO
        output = StringIO.StringIO()

        # Write to buffer, get the lines
        self.write(output, format="ascii.commented_header")
        data_lines = output.getvalue().split("\n")
        #print("datalines", len(data_lines))

        # Get masks
        # masks = self.get_masks()

        # Create header
        header = []
        header.append("PTS data format")
        # for name in masks: header.append(name + " mask: " + tostr(masks[name])) # WILL BE READ FROM THE QUOTE CHARACTERS in the data lines

        # Set density and brightness lists
        if "density" in self.meta and len(self.meta["density"]) > 0: header.append("density: " + tostr(self.meta["density"]))
        if "brightness" in self.meta and len(self.meta["brightness"]) > 0: header.append("brightness: " + tostr(self.meta["brightness"]))

        # Set type string line for the header
        type_string = ""
        for name in self.column_names:
            dtype = self.column_type(name)
            type_string += " " + str(dtype)
        header.append(type_string.strip())

        # Set unit string line for the header
        unit_string = ""
        for name in self.column_names:
            unit = self.column_unit(name)
            if unit is None or unit == "": unit_string += ' ""'
            else: unit_string += " " + tostr(unit)
        header.append(unit_string.strip())

        # SET META
        for key in self.meta:
            if key == "density": continue
            if key == "brightness": continue
            if not types.is_string_type(key): raise ValueError("Keys must be strings")
            value = self.meta[key]
            ptype, string = stringify(value)
            line = "META: " + key + " [" + ptype + "] " + string
            header.append(line)

        # Create lines
        lines = []
        for line in header: lines.append("# " + line)
        # lines.append("# " + data_lines[0]) # add line with the column names
        # for line in data_lines[1:]:
        for line in data_lines:
            if not line: continue  # empty line at the end
            lines.append(line)
        #print("lines", len(lines))

        # Write the lines
        fs.write_lines(path, lines)

    # -----------------------------------------------------------------

    def saveto_csv(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Check or add extension
        if fs.has_extension(path):
            if fs.get_extension(path) != "csv": warnings.warn("The extension is not 'csv'")
        else: path = path + ".csv"

        # Save
        self.saveto(path, format="csv")

    # -----------------------------------------------------------------

    def saveto_ecsv(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get masks
        masks = self.get_masks_int()

        # Set masks in meta
        for name in masks: self.meta[name + " mask"] = masks[name]

        # Replace masked values (not masked anymore)
        self.replace_masked_values()

        # Save
        self.write(path, format="ascii.ecsv")

        # Set the masks back (because they were set to False by replace_masked_values, necessary to avoid writing out
        # '""' (empty string) for each masked value, which is unreadable by Astropy afterwards)
        self.set_masks(masks)

        # Check or add extension
        #if fs.has_extension(path):
        #    if fs.get_extension(path) != "ecsv": warnings.warn("The extension is not 'ecsv'")
        #else: path = path + ".ecsv"

        # Save
        #self.saveto(path, format="ecsv")

    # -----------------------------------------------------------------

    def saveto_html(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Check or add extension
        if fs.has_extension(path):
            if fs.get_extension(path) != "html": warnings.warn("The extension if not 'html'")
        else: path = path + ".html"

        # Save
        self.saveto(path, format="html")

    # -----------------------------------------------------------------

    def saveto_latex(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Check or add extension
        if fs.has_extension(path):
            if fs.get_extension(path) != "text": warnings.warn("The extension is not 'tex'")
        else: path = path + ".tex"

        # Save
        self.saveto(path, format="latex")

    # -----------------------------------------------------------------

    def saveto_votable(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Check or add extension
        if fs.has_extension(path):
            if fs.get_extension(path) != "votable": warnings.warn("The extension is not 'xml'")
        else: path = path + ".xml"

        # Save
        self.saveto(path, format="votable")

    # -----------------------------------------------------------------

    def saveto_ascii(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Check or add extension
        if fs.has_extension(path):
            if fs.get_extension(path) != "ascii": warnings.warn("The extension is not 'ascii'")
        else: path = path + ".ascii"

        # Save
        self.saveto(path, format="ascii")

    # -----------------------------------------------------------------

    # def saveto_pts(self, path):
    #
    #     """
    #     This function ...
    #     :param path:
    #     :return:
    #     """
    #
    #     # Check or add extension
    #     if fs.has_extension(path):
    #         if fs.get_extension(path) != "dat": warnings.warn("The extension is not 'dat'")
    #     else: path = path + ".dat"
    #
    #     # Save
    #     self.saveto(path, format="pts")

    # -----------------------------------------------------------------

    def get_masks(self):

        """
        This function ...
        :return: 
        """

        masks = dict()
        for name in self.colnames:
            masks[name] = [boolean for boolean in self[name].mask]
        return masks

    # -----------------------------------------------------------------

    def get_masks_int(self):

        """
        This function ...
        :return:
        """

        masks = dict()
        for name in self.colnames:
            masks[name] = [int(boolean) for boolean in self[name].mask]
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

    def find_index(self, key, column_name=None, where=None):

        """
        This function ...
        :param key:
        :param column_name:
        :param where:
        :return:
        """

        from ..tools import tables
        return tables.find_index(self, key, column_name, where=where)

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
        return self.colnames

    # -----------------------------------------------------------------

    @property
    def units(self):
        units = []
        for name in self.column_names: units.append(self.column_unit(name))
        return units

    # -----------------------------------------------------------------

    @property
    def column_units(self):
        return self.units

    # -----------------------------------------------------------------

    @property
    def unit_strings(self):
        strings = []
        for name in self.column_names: strings.append(self.column_unit_string(name))
        return strings

    # -----------------------------------------------------------------

    @property
    def column_unit_strings(self):
        return self.unit_strings

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
        return len(self.column_names)

    # -----------------------------------------------------------------

    @property
    def nrows(self):
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

    def print_latex(self, round=False, ndecimal_places=3):

        """
        This function ...
        :param round:
        :param ndecimal_places:
        :return: 
        """

        from ..tools.stringify import tostr

        #self.show_in_browser()

        colnames_escaped = [name.replace("_", "\_") for name in self.colnames]

        header = " & ".join(colnames_escaped) + " \\\\"
        print(header)

        units = []
        has_units = False
        for name in self.colnames:
            unit = self.column_unit(name)
            if unit is None:
                units.append("")
            else:
                has_units = True
                units.append(str(unit))

        if has_units:
            units_string = " & ".join(units) + " \\\\"
            print(units_string)

        for index in range(len(self)):

            row = []
            for name in self.colnames:
                string = tostr(self[name][index], round=round, decimal_places=ndecimal_places)
                row.append(string.replace("_", "\_"))
            row_string = " & ".join(row) + " \\\\"

            print(row_string)

# -----------------------------------------------------------------

def set_table_masks(table):

    """
    This function ...
    :param table:
    :return:
    """

    # Look for masks
    for colname in table.colnames:
        key = colname + " mask"
        if key not in table.meta: continue
        if len(table) != len(table.meta[key]): raise IOError("Length of the table does not correspond to the length of the masks")
        for index in range(len(table)): table[colname].mask[index] = table.meta[key][index]
        del table.meta[key]

# -----------------------------------------------------------------

def fix_column_types(table):

    """
    This function ...
    :param table:
    :return:
    """

    to_convert = dict()

    # Loop over the columns
    for name in table.colnames:

        # Get the type
        dtype_name = table[name].dtype.name

        # OK?
        if is_allowed_dtype_name(dtype_name): continue

        # Get property type
        values = list(table[name])
        ptype = get_common_property_type(values)
        #print(name, ptype)

        # Get builtin type
        real_type = property_type_to_builtin(ptype)

        # Set to convert
        to_convert[name] = real_type

    #print(to_convert)

    # Do the conversion
    for name in to_convert: table[name] = table[name].astype(to_convert[name])

# -----------------------------------------------------------------

def fix_nones(table):

    """
    This function ...
    :param table:
    :return:
    """

    # Loop over the string columns
    for name in table.colnames:
        if dtype_to_column_type(table[name].dtype) != "string": continue

        values = table[name]
        indices = sequences.find_indices(values, "None")

        #print(indices)

        # Mask
        for index in indices:
            table[name][index] = string_column_none_default
            table[name].mask[index] = True

# -----------------------------------------------------------------

def initialize_table(table, table_name=None):

    """
    This function ...
    :param table:
    :param table_name:
    :return:
    """

    from .log import log

    if table_name is None: table_name = ""
    else: table_name = table_name + " "

    # Debugging
    log.debug("Initializing the " + table_name + "table ...")

    # Clear the column info so that we can rebuild it
    table.column_info = []

    # Set the column info
    # Loop over the columns
    for name in table.colnames:

        # Get the type
        simple_dtype = dtype_to_builtin(table[name].dtype)

        # Get unit of the column
        unit = table.column_unit(name)

        # Set this unit object as the actual column unit (so it can be a PhotometricUnit)
        table[name].unit = unit

        # Add column info
        table.add_column_info(name, simple_dtype, unit, None)

        # Initialize "density" meta
        if "density" not in table.meta: table.meta["density"] = []

        # Initialize "brightness" meta
        if "brightness" not in table.meta: table.meta["brightness"] = []

# -----------------------------------------------------------------

def reorder_columns(table):

    """
    This function ...
    :param table:
    :return:
    """

    if not hasattr(table, "_column_info"): return
    #print(table._column_info)

    # Get the column names in the intended order
    column_names_ordered = table._column_info.keys()

    # Check if the table has its columns in the right order
    if table.column_names == column_names_ordered: return

    # Check whether the table contains the same columns as it should
    #if not sequences.same_contents(table.column_names, column_names_ordered): raise ValueError("Table does not contain the correct columns")
    extra_columns = sequences.elements_not_in_other(table.column_names, column_names_ordered, check_existing=True)
    #print(extra_columns)

    # Give warning
    warnings.warn("The order of the columns does not correspond to the order of the columns as defined by the table subclass: re-ordering the columns ...")
    if len(extra_columns) > 0: warnings.warn("The following extra columns exist: '" + ",".join(extra_columns) + "'. Assumed to be the last columns in this order.")

    # Loop over the column names in the desired order
    for index, column_name in enumerate(column_names_ordered):

        # Delete column and readd at right index
        column = table[column_name]
        del table[column_name]
        table.add_column(column, index=index)

        # Remove and re-add the column info
        info = table._column_info[column_name]
        table.column_info.pop(index)
        table.add_column_info(column_name, info[0], info[1], info[2], index=index)

# -----------------------------------------------------------------

all_column_types = ["string", "real", "integer", "boolean"]
dtype_name_startswith_strings = ["string", "float", "int", "bool"]

# -----------------------------------------------------------------

def dtype_to_column_type(dtype):

    """
    This function ...
    :param dtype:
    :return:
    """

    if np.issubdtype(dtype, np.string_): return "string"
    elif np.issubdtype(dtype, np.float): return "real"
    elif np.issubdtype(dtype, np.int): return "integer"
    elif np.issubdtype(dtype, np.bool): return "boolean"
    else: raise ValueError("Did not recognize the dtype '" + str(dtype) + "'")

# -----------------------------------------------------------------

def dtype_to_builtin(dtype):

    """
    This function ...
    :param dtype:
    :return:
    """

    if np.issubdtype(dtype, np.string_): return str
    elif np.issubdtype(dtype, np.float): return float
    elif np.issubdtype(dtype, np.int): return int
    elif np.issubdtype(dtype, np.bool): return bool
    else: raise ValueError("Did not recognize the dtype '" + str(dtype) + "'")

# -----------------------------------------------------------------

def dtype_name_to_column_type(name):

    """
    This function ...
    :param name:
    :return:
    """

    if name.startswith("string"): return "string"
    elif name.startswith("float"): return "real"
    elif name.startswith("int"): return "integer"
    elif name.startswith("bool"): return "boolean"
    else: raise ValueError("Unknown column type: '" + name + "'")

# -----------------------------------------------------------------

def column_type_to_builtin(column_type):

    """
    This function ...
    :param column_type:
    :return:
    """

    if column_type == "string": return str
    elif column_type == "boolean": return bool
    elif column_type == "integer": return int
    elif column_type == "real": return float
    else: raise ValueError("Unrecognized column type string: '" + column_type + "'")

# -----------------------------------------------------------------

def composed_of_multiple_builtins(ptype):

    """
    Thisf unction ...
    :param ptype:
    :return:
    """

    if ptype.endswith("coordinate"): return True
    elif ptype.endswith("extent"): return True
    elif ptype.endswith("dictionary"): return True
    else: return False

# -----------------------------------------------------------------

def property_type_to_builtins(ptype):

    """
    This function ...
    :param ptype:
    :return:
    """

    if ptype == "skycoordinate": return False, [("ra", float, True), ("dec", float, True)]
    elif ptype == "pixelcoordinate": return False,  [("x", float, True), ("y", float, True)]
    elif ptype == "physicalcoordinate": return False, [("x", float, True), ("y", float, True)]
    elif ptype.endswith("extent"): return False, [("x", float, True), ("y", float, True)]
    elif ptype.endswith("dictionary"): return False, []
    else: raise ValueError("Not recognized as property type that consists of multiple builtins")

# -----------------------------------------------------------------

def property_type_to_builtin(ptype, return_tostring=False):

    """
    This function ...
    :param ptype:
    :param return_tostring:
    :return:
    """

    # Real numbers
    if ptype.endswith("real"): real_type, tostring = float, False

    # Integers
    elif ptype.endswith("integer"): real_type, tostring = int, False

    # Strings
    elif ptype.endswith("string"): real_type, tostring = str, False

    # Booleans
    elif ptype.endswith("boolean"): real_type, tostring = bool, False

    # Quantities
    elif ptype.endswith("quantity"): real_type, tostring = float, False

    # UNIT HAVE TO BE CONVERTED TO STRINGS
    elif ptype.endswith("unit"): real_type, tostring = str, True

    # FILTERS HAVE TO BE CONVERTED TO STRINGS
    elif ptype.endswith("filter"): real_type, tostring = str, True

    # SPECIAL CASE: WAS NONE FOR EACH COMPOSITE
    elif ptype == "None": real_type, tostring = str, False

    # LISTS OF THINGS -> STRINGS
    elif ptype.endswith("_list"): real_type, tostring = str, True

    # Dictionary
    #elif ptype.endswith("dictionary"): real_type, tostring = str, True

    # Coordinate
    #elif ptype.endswith("coordinate"): real_type, tostring = str, True

    # NOT RECOGNIZED
    else: raise ValueError("Column type not recognized: " + str(ptype) + " (" + str(type(ptype).__name__) + ")")

    # Return
    if return_tostring: return real_type, tostring
    else: return real_type

# -----------------------------------------------------------------

def is_allowed_dtype_name(name):

    """
    This function ...
    :param name:
    :return:
    """

    from ..tools import strings
    return strings.startswith_any(name, dtype_name_startswith_strings)

# -----------------------------------------------------------------

def get_common_property_type(values):

    """
    This function ...
    :param values:
    :return:
    """

    from ..tools.stringify import stringify

    types = []
    units = []

    # Loop over the values
    for value in values:

        if value is None: dtype = unit = None
        else:

            # Check whether there is a unit
            if hasattr(value, "unit"): unit = value.unit
            else: unit = None

            # Get the type
            dtype, string = stringify(value)

        # Add the type and
        types.append(dtype)
        units.append(unit)

    # Determine column type, unit and description
    if sequences.all_equal_to(types, 'None') or sequences.all_none(types): ptype = 'None'
    else: ptype = sequences.get_all_equal_value(types, ignore_none=True, ignore='None')

    # Return the common type
    return ptype

# -----------------------------------------------------------------

def parse_header_file(filepath, format):

    """
    This function ...
    :param filepath:
    :param format:
    :return:
    """

    # PTS table format
    if format == "pts": return parse_pts_header_file(filepath)

    # CSV format
    elif format == "csv":

        column_names, column_units = fs.get_column_names(filepath, return_units=True)
        return column_names, None, column_units, None

    # Invalid format
    else: raise ValueError("Invalid format")

# -----------------------------------------------------------------

def parse_pts_header_file(filepath):

    """
    This function ...
    :param header:
    :return:
    """

    header = fs.get_header_lines(filepath)
    return parse_pts_header(header)

# -----------------------------------------------------------------

def parse_pts_header(header):

    """
    This function ...
    :param lines:
    :return:
    """

    # Put last header line (colum names) as first data line (and remove it from the header)
    #sequences.prepend(data_lines, "# " + header[-1])
    #header = header[:-1]

    # Search for density and brightness, set meta info
    index = -1
    density = None
    brightness = None
    for index, line in enumerate(header):

        if line.startswith("PTS data format"): continue

        if line.startswith("density:"):

            density_string = line.split("density: ")[1]
            string = "[" + ",".join('"' + s + '"' for s in density_string.split(",")) + "]"
            density = eval(string)

        elif line.startswith("brightness:"):

            brightness_string = line.split("brightness: ")[1]
            string = "[" + ",".join('"' + s + '"' for s in brightness_string.split(",")) + "]"
            brightness = eval(string)

        else: break

    # Get the META info
    from ..tools import parsing
    meta = OrderedDict()
    old_header = header[:]
    header = []
    for line in old_header:
        if line.startswith("META"):
            spec = line.split("META: ")[1]
            key = spec.split(" [")[0]
            ptype = spec.split("[")[1].split("]")[0]
            if ptype == "None":
                meta[key] = None
            else:
                string = spec.split("] ")[1]
                parsing_function = getattr(parsing, ptype)
                value = parsing_function(string)
                meta[key] = value
        else: header.append(line)

    # print(meta)

    # Contains type line
    if index == len(header) - 3:

        # Get types
        types_string = header[-3]
        column_type_strings = types_string.split()
        column_types = []
        for type_string in column_type_strings:
            dtype = column_type_to_builtin(type_string)
            column_types.append(dtype)

    # Doesn't contain type line
    elif index == len(header) - 1: column_types = column_type_strings = None

    # Invalid
    else: raise IOError("Something is wrong with the file")

    # Get column names
    from ..tools import strings
    column_names = strings.split_except_within_double_quotes(header[-1], add_quotes=False)

    # Get column units
    column_units = dict()
    unit_string = header[-2]
    #unit_strings = unit_string.split()
    #print(unit_strings, column_names)
    unit_strings = strings.split_except_within_round_brackets_and_double_quotes(unit_string)
    #print(unit_strings, len(unit_strings))
    #print(column_names, len(column_names))
    assert len(unit_strings) == len(column_names)

    for unit_string, colname in zip(unit_strings, column_names):
        if unit_string == '""': continue
        #table[colname].unit = unit_string
        is_brightness = colname in brightness if brightness is not None else False
        is_density = colname in density if density is not None else False
        unit = u(unit_string, brightness=is_brightness, density=is_density)
        column_units[colname] = unit

    # Return
    return column_names, column_types, column_units, meta

# -----------------------------------------------------------------
