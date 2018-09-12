#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.representations.base Contains the RepresentationBuilderBase class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component import BuildComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ..dustgrid import DustGridBuilder
from ..representation import Representation

# -----------------------------------------------------------------

class RepresentationBuilderBase(BuildComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RepresentationBuilderBase, self).__init__(*args, **kwargs)

        # The model definition
        self.definition = None

        # The representation
        self.representation = None

        # The projections
        self.projections = dict()

        # The instruments
        self.instruments = dict()

        # The dust grid
        self.dust_grid = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(RepresentationBuilderBase, self).setup(**kwargs)

        # Create the model definition
        self.definition = self.get_model_definition(self.config.model_name)

        # Create the representation
        path = fs.create_directory_in(self.representations_path, self.config.name)
        self.representation = Representation(self.config.name, self.config.model_name, path)

    # -----------------------------------------------------------------

    @property
    def representation_name(self):
        return self.representation.name

    # -----------------------------------------------------------------

    @property
    def representation_path(self):
        return self.representation.path

    # -----------------------------------------------------------------

    @property
    def model_name(self):
        return self.representation.model_name

    # -----------------------------------------------------------------

    def build_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the dust grid ...")

        # Create the builder
        builder = DustGridBuilder()

        # Set output path
        builder.config.output = self.representation.grid_path

        # Set simulation path
        builder.config.simulation_path = self.representation.grid_out_path

        # Set whether quality has to be calculated
        builder.config.quality = self.config.check_dust_grid_quality

        # Run the builder
        builder.run(definition=self.definition, dust_grid=self.dust_grid)

# -----------------------------------------------------------------
