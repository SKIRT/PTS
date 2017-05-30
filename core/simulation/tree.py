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

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ..tools import tables
from . import textfile
from ..basics.range import QuantityRange

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

        """
        This function ...
        :return: 
        """

        for node in self.nodes:
            if node.nchildren != 0: return node.nchildren
        return None

    # -----------------------------------------------------------------

    @property
    def tree_type(self):

        """
        This function ...
        :return: 
        """

        if self.nchildren_per_node == 2: return "bintree"
        elif self.nchildren_per_node == 8: return "octtree"
        else: raise ValueError("Invalid number of children per node: " + str(self.nchildren_per_node))

    # -----------------------------------------------------------------

    @property
    def nnodes(self):

        """
        THis function ...
        :return: 
        """

        return len(self.nodes)

    # -----------------------------------------------------------------

    @property
    def nleaves(self):

        """
        This function ...
        :return: 
        """

        nchildless = 0
        for node in self.nodes:
            if node.nchildren == 0: nchildless += 1
        return nchildless

    # -----------------------------------------------------------------

    @lazyproperty
    def root(self):

        """
        This function ...
        :return: 
        """

        return self.nodes[0]

    # -----------------------------------------------------------------

    @property
    def x_range(self):

        """
        This function ...
        :return: 
        """

        return self.root.x_range

    # -----------------------------------------------------------------

    @property
    def y_range(self):

        """
        This function ...
        :return: 
        """

        return self.root.y_range

    # -----------------------------------------------------------------

    @property
    def z_range(self):

        """"
        This function ...
        """

        return self.root.z_range
    
    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path: 
        :return: 
        """

        # Load the table
        table = tables.from_file(path, format="ascii")

        # Load the descriptions and the units
        #descriptions, units = textfile.get_descriptions_and_units(path)
        units = textfile.get_units(path)

        # Get the unit of length
        length_unit = units[2]

        # Create tre
        tree = cls()
        tree.path = path

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
        table.rename_column("col1", "ID")
        table.rename_column("col2", "Cell index")
        table.rename_column("col3", "Min x")
        table.rename_column("col4", "Max x")
        table.rename_column("col5", "Min y")
        table.rename_column("col6", "Max y")
        table.rename_column("col7", "Min z")
        table.rename_column("col8", "Max z")
        table.rename_column("col9", "Parent ID")
        table.rename_column("col10", "Child 0 ID")
        table.rename_column("col11", "Child 1 ID")
        table.rename_column("col12", "Child 2 ID")
        table.rename_column("col13", "Child 3 ID")
        table.rename_column("col14", "Child 4 ID")
        table.rename_column("col15", "child 5 ID")
        table.rename_column("col16", "child 6 ID")
        table.rename_column("col17", "child 7 ID")

        # Loop over the rows
        for index in range(len(table)):

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

    def saveto(self, path):

        """
        This function ...
        :param path: 
        :return: 
        """

        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return: 
        """

        self.saveto(self.path)

# -----------------------------------------------------------------
