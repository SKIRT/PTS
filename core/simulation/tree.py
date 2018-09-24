#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.tree Contains the DustGridTree class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import tables
from . import textfile
from ..basics.range import QuantityRange
from ..tools.utils import lazyproperty
from ..basics.distribution import Distribution
from ..tools import filesystem as fs
from ..basics.table import SmartTable
from ..basics.log import log
from ..tools.progress import Bar
from ..tools import numbers

# -----------------------------------------------------------------

# column 1: node ID
# column 2: dust cell index
# column 3: minimum x coordinate of the node (pc)
# column 4: maximum x coordinate of the node (pc)
# column 5: minimum y coordinate of the node (pc)
# column 6: maximum y coordinate of the node (pc)
# column 7: minimum z coordinate of the node (pc)
# column 8: maximum z coordinate of the node (pc)
# column 9: ID of the father node
# column 10: ID of child node 0
# column 11: ID of child node 1
# column 12: ID of child node 2
# column 13: ID of child node 3
# column 14: ID of child node 4
# column 15: ID of child node 5
# column 16: ID of child node 6
# column 17: ID of child node 7

original_column_names = ["node ID", "dust cell index", "minimum x coordinate of the node", "maximum x coordinate of the node",
                         "minimum y coordinate of the node", "maximum y coordinate of the node",
                         "minimum z coordinate of the node", "maximum z coordinate of the node",
                         "ID of the father node", "ID of child node 0", "ID of child node 1", "ID of child node 2",
                         "ID of child node 3", "ID of child node 4", "ID of child node 5", "ID of child node 6", "ID of child node 7"]

column_names = ["ID",  "Cell index", "Min x", "Max x", "Min y", "Max y", "Min z", "Max z", "Parent ID", "Child 0 ID", "Child 1 ID", "Child 2 ID", "Child 3 ID", "Child 4 ID", "Child 5 ID", "Child 6 ID", "Child 7 ID"]

# -----------------------------------------------------------------

class TreeNode(object):

    """
    This function ...
    """

    def __init__(self, id, cell, x_range, y_range, z_range, parent, children):

        """
        THis function ...
        :param id: 
        :param cell: 
        :param x_range: 
        :param y_range: 
        :param z_range: 
        :param parent: 
        :param children: 
        """

        self.id = id
        self.cell = cell
        self.x_range = x_range
        self.y_range = y_range
        self.z_range = z_range
        self.parent = parent
        self.children = children

    # -----------------------------------------------------------------

    @property
    def nchildren(self):

        """
        This function ...
        :return: 
        """

        return len(self.children)

# -----------------------------------------------------------------

def get_nnodes(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    # Example:
    # Total number of nodes: 5787413

    # Read header of the file
    lines = fs.get_header_lines(filepath)

    # Look for number of nodes
    nnodes = None
    for line in lines:
        if "Total number of nodes" in line:
            nnodes = int(line.split(": ")[1])
            break
    if nnodes is None: raise IOError("Something went wrong")

    # Return
    return nnodes

# -----------------------------------------------------------------

def get_nleaves(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    # Example:
    # Total number of leaves: 2893707

    # Read header of the file
    lines = fs.get_header_lines(filepath)

    # Look for the number of leaves
    nleaves = None
    for line in lines:
        if "Total number of leaves" in line:
            nleaves = int(line.split(": ")[1])
            break
    if nleaves is None: raise IOError("Something went wrong")

    # Return
    return nleaves

# -----------------------------------------------------------------

class DustGridTree(object):

    """
    This function ...
    """

    def __init__(self):

        """
        This function ...
        """

        # The filepath
        self.path = None

        # The nodes
        self.nodes = []

    # -----------------------------------------------------------------

    @property
    def nchildren_per_node(self):
        for node in self.nodes:
            if node.nchildren != 0: return node.nchildren
        return None

    # -----------------------------------------------------------------

    @property
    def tree_type(self):
        if self.nchildren_per_node == 2: return "bintree"
        elif self.nchildren_per_node == 8: return "octtree"
        else: raise ValueError("Invalid number of children per node: " + str(self.nchildren_per_node))

    # -----------------------------------------------------------------

    @property
    def nnodes(self):
        return len(self.nodes)

    # -----------------------------------------------------------------

    @property
    def nleaves(self):
        nchildless = 0
        for node in self.nodes:
            if node.nchildren == 0: nchildless += 1
        return nchildless

    # -----------------------------------------------------------------

    @lazyproperty
    def root(self):
        return self.nodes[0]

    # -----------------------------------------------------------------

    @property
    def x_range(self):
        return self.root.x_range

    # -----------------------------------------------------------------

    @property
    def y_range(self):
        return self.root.y_range

    # -----------------------------------------------------------------

    @property
    def z_range(self):
        return self.root.z_range

    # -----------------------------------------------------------------

    @classmethod
    def from_table(cls, table, units, fix_column_names=True):

        """
        This function ...
        :param table:
        :param units:
        :param fix_column_names:
        :return:
        """

        # Get the unit of length
        length_unit = units[2]

        # Create the tree
        tree = cls()

        # Rename the columns
        if fix_column_names:
            for index in range(len(table.colnames)):

                column_name = "col" + str(index+1)
                original_column_name = original_column_names[index]
                new_column_name = column_names[index]

                # Original column names have been read in
                if original_column_name in table.colnames: table.rename_column(original_column_name, new_column_name)

                # Column names are col0, col1, etc.
                else: table.rename_column(column_name, new_column_name)

        # Debugging
        log.debug("Processing the tree nodes, this may take a while ...")

        # Show progress bar
        nrows = len(table)
        with Bar(label='', expected_size=nrows, every=1, add_datetime=True) as bar:

            # Loop over the rows
            for index in range(len(table)):

                # Show progress
                bar.show(float(index + 1))

                # Get ID and index
                id = table["ID"][index]
                m = table["Cell index"][index]
                if m == -1: m = None

                # Get extents
                min_x = table["Min x"][index]
                max_x = table["Max x"][index]
                min_y = table["Min y"][index]
                max_y = table["Max y"][index]
                min_z = table["Min z"][index]
                max_z = table["Max z"][index]

                # Create ranges
                x_range = QuantityRange(min_x, max_x, unit=length_unit)
                y_range = QuantityRange(min_y, max_y, unit=length_unit)
                z_range = QuantityRange(min_z, max_z, unit=length_unit)

                # Get parent
                parent = table["Parent ID"][index]
                if parent == -1: parent = None

                # Get children
                children = []
                for j in range(8):
                    child_id = table["Child " + str(j) + " ID"][index]
                    if child_id == -1: break
                    else: children.append(child_id)

                # Create the node
                node = TreeNode(id, m, x_range, y_range, z_range, parent, children)

                # Add the node
                tree.nodes.append(node)

        # Return the tree
        return tree

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, method="pandas"):

        """
        This function ...
        :param path:
        :param method:
        :return:
        """

        # Load the table
        if method == "astropy": table = tables.from_file(path, format="ascii")
        else: table = SmartTable.from_file(path, method=method, format="csv")

        # Load the descriptions and the units
        # descriptions, units = textfile.get_descriptions_and_units(path)
        units = textfile.get_units(path)
        #print(units)

        # Create
        tree = cls.from_table(table, units)

        # Set path
        tree.path = path

        # Return the tree
        return tree

    # -----------------------------------------------------------------

    @classmethod
    def from_remote_file(cls, path, remote):

        """
        This function ...
        :param path:
        :param remote:
        :return:
        """

        # Load the table
        table = tables.from_remote_file(path, remote, format="ascii")

        # Load the descriptions and the units
        # descriptions, units = textfile.get_descriptions_and_units(path)
        units = textfile.get_units(path)

        # Create
        tree = cls.from_table(table, units)

        # Return the tree
        return tree

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path: 
        :return: 
        """

        # Inform the user
        #log.info("Saving the dust grid tree to '" + path + "' ...")

        id_column = []
        index_column = []
        min_x_column = []
        max_x_column = []
        min_y_column = []
        max_y_column = []
        min_z_column = []
        max_z_column = []
        parent_id_column = []
        child0_column = []
        child1_column = []
        child2_column = []
        child3_column = []
        child4_column = []
        child5_column = []
        child6_column = []
        child7_column = []

        length_unit = None

        # Loop over the nodes
        for node in self.nodes:

            id_column.append(node.id)
            index_column.append(node.cell)

            if length_unit is None: length_unit = node.x_range.min.unit
            elif length_unit != node.x_range.min.unit: raise ValueError("Something went wrong")

            min_x_column.append(node.x_range.min.to(length_unit).value)
            max_x_column.append(node.x_range.max.to(length_unit).value)
            min_y_column.append(node.y_range.min.to(length_unit).value)
            max_y_column.append(node.y_range.max.to(length_unit).value)
            min_z_column.append(node.z_range.min.to(length_unit).value)
            max_z_column.append(node.z_range.max.to(length_unit).value)

            # Parent
            parent_id = node.parent if node.parent is not None else -1
            parent_id_column.append(parent_id)

            # Children
            if len(node.children) > 0:
                child0_column.append(node.children[0])
                child1_column.append(node.children[1])
            else:
                child0_column.append(-1)
                child1_column.append(-1)

            if len(node.children) > 2:

                child2_column.append(node.children[2])
                child3_column.append(node.children[3])

            else:
                child2_column.append(-1)
                child3_column.append(-1)

            if len(node.children) > 4:

                child4_column.append(node.children[4])
                child5_column.append(node.children[5])
                child6_column.append(node.children[6])
                child7_column.append(node.children[7])

            else:
                child4_column.append(-1)
                child5_column.append(-1)
                child6_column.append(-1)
                child7_column.append(-1)

        # Fill the columns
        data = [id_column, index_column, min_x_column, max_x_column, min_y_column, max_y_column, min_z_column, max_z_column, parent_id_column, child0_column, child1_column, child2_column, child3_column, child4_column, child5_column, child6_column, child7_column]

        # Create table
        table = tables.new(data, column_names)

        # Save the table
        tables.write(table, path)

        # Set the path
        self.path = path

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return: 
        """

        # Save the pathc
        self.saveto(self.path)

# -----------------------------------------------------------------

class DustGridTreeDistribution(Distribution):

    """
    This class ...
    """

    @property
    def min_level_index(self):
        for index in range(len(self)):
            if self.frequencies[index] != 0: return index

    # -----------------------------------------------------------------

    @property
    def max_level_index(self):
        for index in reversed(range(len(self))):
            if self.frequencies[index] != 0: return index

    # -----------------------------------------------------------------

    @property
    def min_level(self):
        return numbers.as_integer_check(self.values[self.min_level_index])

    # -----------------------------------------------------------------

    @property
    def max_level(self):
        return numbers.as_integer_check(self.values[self.max_level_index])

    # -----------------------------------------------------------------

    @property
    def ncells(self):
        total = 0
        for index in range(self.min_level_index, self.max_level_index+1):
            total += self.frequencies[index]
        return numbers.as_integer_check(total)

    # -----------------------------------------------------------------

    def get_ncells_below_level(self, level, including=False):
        total = 0
        stop = level + 1 if including else level
        for index in range(self.min_level_index, stop):
            total += self.frequencies[index]
        return numbers.as_integer_check(total)

    # -----------------------------------------------------------------

    def get_ncells_above_level(self, level, including=False):
        total = 0
        start = level if including else level + 1
        for index in range(start, self.max_level_index+1):
            total += self.frequencies[index]
        return numbers.as_integer_check(total)

# -----------------------------------------------------------------

def get_cell_coordinates(tree_filepath, read_method="numpy"):

    """
    This function ...
    :param tree_filepath:
    :param read_method:
    :return:
    """

    import numpy as np

    # Get the number of cells
    ncells = get_nleaves(tree_filepath)

    # Initialize arrays
    x_min = np.full((ncells,), np.nan)
    x_max = np.full((ncells,), np.nan)
    y_min = np.full((ncells,), np.nan)
    y_max = np.full((ncells,), np.nan)
    z_min = np.full((ncells,), np.nan)
    z_max = np.full((ncells,), np.nan)

    # Read columns
    indices, xmin, xmax, ymin, ymax, zmin, zmax = fs.get_columns(tree_filepath, indices=[1,2,3,4,5,6,7], method=read_method)

    # Loop over the tree nodes
    for i in range(len(indices)):

        # Not a leaf?
        index = indices[i]
        if index == - 1: continue
        #print(index)

        # Add
        x_min[index] = xmin[i]
        x_max[index] = xmax[i]
        y_min[index] = ymin[i]
        y_max[index] = ymax[i]
        z_min[index] = zmin[i]
        z_max[index] = zmax[i]

    # Return
    return x_min, x_max, y_min, y_max, z_min, z_max

# -----------------------------------------------------------------
