#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.layers Contains the Layers class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# -----------------------------------------------------------------

class Layers(OrderedDict):

    """
    This class ...
    """

    def __setitem__(self, key, value):

        """
        This function ...
        :param key:
        :param value:
        :return:
        """

        if isinstance(key, int):
            if key >= len(self): raise IndexError("layer index out of range")
            key = self.keys()[key]
        elif isinstance(key, basestring):
            if " " in key: raise KeyError("Cannot use a string key that contains spaces")
        else: raise KeyError("Cannot use a key that is not a string or an integer")

        super(Layers, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    # -----------------------------------------------------------------

    def __getitem__(self, key):

        """
        This function ...
        :param key:
        :return:
        """

        if isinstance(key, int):
            if key >= len(self): raise IndexError("layer index out of range")
            key = self.keys()[key]
        elif isinstance(key, basestring):
            if " " in key: raise KeyError("Cannot use a string key that contains spaces")
        else: raise KeyError("Cannot use a key that is not a string or an integer")

        return super(Layers, self).__getitem__(key)

    # -----------------------------------------------------------------

    def as_list(self, copy=False):

        """
        This function ...
        :param copy:
        :return:
        """

        layer_list = []
        for name in self:
            if copy: layer_list.append(self[name].copy())
            else: layer_list.append(self[name])
        return layer_list

    # -----------------------------------------------------------------

    def get_first(self):

        """
        This function ...
        :return:
        """

        return self[self.keys()[0]]

    # -----------------------------------------------------------------

    def selected(self, require_single=False, allow_none=True):

        """
        This function ...
        :param require_single:
        :param allow_none:
        :return:
        """

        # Initialize a list to contain the selected frames
        layers = []

        # Loop over all layers, add them to the list if selected
        for name, layer in self.items(): layers.append(layer)

        # Check options
        if require_single and len(layers) > 1: raise RuntimeError("More than one layer is selected")
        if not allow_none and len(layers) == 0: raise RuntimeError("No layer is selected")

        # Return the list of layers (or the single selected layer if requested)
        return layers[0] if require_single else layers

    # -----------------------------------------------------------------

    def get_selected(self, require_single=False, allow_none=True):

        """
        This function returns a list of the names of the currently selected layers
        :param require_single:
        :param allow_none:
        :return:
        """

        # Initialize an empty list
        names = []

        # Loop over all layers
        for name in self.keys():

            # If this layer is currently selected, add it to the list
            if self[name].selected: names.append(name)

        if require_single and len(names) > 1: raise RuntimeError('More than one layer is selected')
        if not allow_none and len(names) == 0: raise RuntimeError('There is no layer selected')

        # TODO: do not use this hack, we should adopt an ordering for the frames, regions and masks (orderedDict?)
        if 'primary' in names:

            names.remove('primary')
            names.insert(0, 'primary')

        # Return the name(s) of the currently selected layers
        return names[0] if require_single else names

    # -----------------------------------------------------------------

    def select_all(self):

        """
        This function selects all the layers
        :return:
        """

        # Select each layer
        for name in self.keys(): self[name].select()

    # -----------------------------------------------------------------

    def deselect_all(self):

        """
        This function deselects all the layers
        :return:
        """

        # Deselect each layer
        for name in self.keys(): self[name].deselect()

    # -----------------------------------------------------------------

    def get_state(self):

        """
        This function ...
        :return:
        """

        # Create an empty dictionary to contain the state of the layers
        state = dict()

        # Loop over all frames, regions and masks and record whether they are selected
        for layer_name in self: state[layer_name] = self[layer_name].selected

        # Return the state dictionary
        return state

    # -----------------------------------------------------------------

    def set_state(self, state):

        """
        This function ...
        :param state:
        :return:
        """

        # Deselect all layers
        self.deselect_all()

        # Loop over the entries in the state dictionary
        for layer_name, selected in state.items():

            # Check if a layer with this name exists and set the appropriate flag
            if layer_name in self: self[layer_name].selected = selected
            else: raise ValueError("Invalid state dictionary")

# -----------------------------------------------------------------
