#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.attenuation Contains the AttenuationDustMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from copy import copy

# Import the relevant PTS classes and modules
from ....magic.core.frame import Frame
from ....core.basics.log import log
from ....core.basics.configurable import Configurable
from ....magic.core.frame import AllZeroError
from ....magic.core.image import Image
# -----------------------------------------------------------------

def make_map(attenuation):

    """
    This function ...
    :param attenuation:
    :return: 
    """

    # Initialize the map maker
    maker = AttenuationDustMapsMaker()

    # Set input
    attenuations = {"attenuation": attenuation}

    # Run the map maker
    maker.run(attenuation=attenuations)

    # Run the map
    return maker.single_map

# -----------------------------------------------------------------

class AttenuationDustMapsMaker(Configurable):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(AttenuationDustMapsMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Input
        self.attenuation = None
        self.attenuation_origins = None
        self.attenuation_methods = None

        # The method name
        self.method_name = None

        # The maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Make the maps
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Call the setup function of the base class
        super(AttenuationDustMapsMaker, self).setup(**kwargs)

        # Get input
        self.attenuation = kwargs.pop("attenuation")

        # Get origins
        self.attenuation_origins = kwargs.pop("attenuation_origins", None)

        # Get methods
        self.attenuation_methods = kwargs.pop("attenuation_methods", None)

        # Get method name for this class
        self.method_name = kwargs.pop("method_name", None)
        if self.has_methods and self.method_name is None: raise ValueError("Method name has to be specified when methods are specified")

        # Get already created maps
        self.maps = kwargs.pop("maps", dict())

    # -----------------------------------------------------------------

    @property
    def has_origins(self):

        """
        This function ...
        :return:
        """

        return self.attenuation_origins is not None

    # -----------------------------------------------------------------

    @property
    def has_methods(self):

        """
        This function ...
        :return:
        """

        return self.attenuation_methods is not None

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the dust maps ...")

        # Loop over the colour maps
        for name in self.attenuation:

            # Set origins
            if self.has_origins:

                # Set the origins
                origins = copy(self.attenuation_origins[name])
                self.origins[name] = origins

            # Set methods
            if self.has_methods:

                # Set the methods
                methods = copy(self.attenuation_methods[name])
                methods.append(self.method_name)
                self.methods[name] = methods

            # Check whether a dust map is already present
            if name in self.maps:
                log.success("The " + name + " dust map is already created: not creating it again")
                continue

            # Get the map
            dust = self.attenuation[name]
            if isinstance(dust, Image): dust_frame = dust.primary
            elif isinstance(dust, Frame): dust_frame = dust
            else: raise ValueError("Something went wrong")

            # Normalized
            try: dust_frame.normalize()
            except AllZeroError:
                log.warning("The " + name + " dust map cannot be normalized: sum is zero")

            # Set as dust map
            self.maps[name] = dust_frame

    # -----------------------------------------------------------------

    @property
    def single_map(self):

        """
        This function ...
        :return:
        """

        if len(self.maps) != 1: raise ValueError("Not a single map")
        return self.maps[self.maps.keys()[0]]

# -----------------------------------------------------------------
