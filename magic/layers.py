#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

class Layers(dict):

    """
    This class is a wrapper around the dict class, with the additional benefit of being able to access its values
    with the 'dot' notation. It is a quite genious way of dealing with a set of layers (image frames, masks or regions
    in this case), with high user-friendliness and easy programming interface.
    """

    # *****************************************************************

    # Use a trick to be able to access attributes of this class by using 'dot' notation ("object.attribute")
    def __getattr__(self, attr): return self.get(attr, None)
    __setattr__= dict.__setitem__   # Set an item of the dictionary
    __delattr__= dict.__delitem__   # Delete an item from the dictionary

    # *****************************************************************

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

        # Return the name(s) of the currently selected layers
        return names[0] if require_single else names

    # *****************************************************************

    def select_all(self):

        """
        This function selects all the layers
        :return:
        """

        # Select each layer
        for name in self.keys(): self[name].select()

    # *****************************************************************

    def deselect_all(self):

        """
        This function deselects all the layers
        :return:
        """

        # Deselect each layer
        for name in self.keys(): self[name].deselect()

    # *****************************************************************

    def get_state(self):

        # Create an empty dictionary to contain the state of the layers
        state = dict()

        # Loop over all frames, regions and masks and record whether they are selected
        for layer_name in self: state[layer_name] = self[layer_name].selected

        # Return the state dictionary
        return state

    # *****************************************************************

    def set_state(self, state):

        # Deselect all layers
        self.deselect_all()

        # Loop over the entries in the state dictionary
        for layer_name, selected in state.items():

            # Check if a layer with this name exists and set the appropriate flag
            if layer_name in self: self[layer_name].selected = selected
            else: raise ValueError("Invalid state dictionary")

# *******************************************************************
