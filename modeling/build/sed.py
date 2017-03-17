#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.sed Contains the SEDModelBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import BuildComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ..component.sed import get_ski_template
from ...core.tools.serialization import write_dict

# -----------------------------------------------------------------

class SEDModelBuilder(BuildComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(SEDModelBuilder, self).__init__(config, interactive)

        # The path for this model
        self.model_path = None

        # The path for the stellar components
        self.model_stellar_path = None

        # The path for the dust components
        self.model_dust_path = None

        # The ski file template
        self.ski = None

        # The stellar properties
        self.stellar_properties = dict()

        # The dust properties
        self.dust_properties = dict()

        # The paths to the component directories
        self.stellar_paths = dict()
        self.dust_paths = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Load the ski file template
        self.load_ski()

        # Get the stellar components
        self.get_stellar_components()

        # Get the dust components
        self.get_dust_components()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.config.name

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SEDModelBuilder, self).setup(**kwargs)

        # Set the model path and create it
        self.model_path = fs.create_directory_in(self.models_path, self.model_name)

        # Set the path of the directory for the stellar components
        self.model_stellar_path = fs.create_directory_in(self.model_path, "stellar")

        # Set the path of the directory for the dust components
        self.model_dust_path = fs.create_directory_in(self.model_path, "dust")

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Loading the ski file
        self.ski = get_ski_template(self.config.path)

    # -----------------------------------------------------------------

    def get_stellar_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the stellar components ...")

        # Loop over the stellar components
        for component_id in self.ski.get_stellar_component_ids():

            # Get the properties
            properties = self.ski.get_stellar_component_properties(component_id)

            # Set the properties
            self.stellar_properties[component_id] = properties

    # -----------------------------------------------------------------

    def get_dust_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the dust components ...")

        # Loop over the dust components
        for component_id in self.ski.get_dust_component_ids():

            # Get the properties
            properties = self.ski.get_dust_component_properties(component_id)

            # Set the properties
            self.dust_properties[component_id] = properties

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write directories
        self.write_stellar_component_directories()

        # Write directories
        self.write_dust_component_directories()

        # Write the stellar properties
        self.write_stellar_properties()

        # Write the dust properties
        self.write_dust_properties()

        # Write the table
        self.write_table()

    # -----------------------------------------------------------------

    def write_stellar_component_directories(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the stellar component directories ...")

        # Loop over the components
        for name in self.stellar_properties:

            # Create a directory
            component_path = self.output_path_file(name)
            fs.create_directory(component_path)

            # Set the path
            self.stellar_paths[name] = component_path

    # -----------------------------------------------------------------

    def write_dust_component_directories(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust component directories ...")

        # Loop over the components
        for name in self.dust_properties:

            # Create a directory
            component_path = self.output_path_file(name)
            fs.create_directory(component_path)

            # Set the path
            self.dust_paths[name] = component_path

    # -----------------------------------------------------------------

    def write_stellar_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing stellar properties ...")

        # Loop over the components
        for name in self.stellar_properties:

            # Write the properties
            path = fs.join(self.stellar_paths[name], "properties.dat")
            write_dict(self.stellar_properties[name], path)

    # -----------------------------------------------------------------

    def write_dust_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing dust properties ...")

        # Loop over the components
        for name in self.dust_properties:

            # Write the properties
            path = fs.join(self.dust_paths[name], "properties.dat")
            write_dict(self.dust_properties[name], path)

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the model table ...")

        # Add the model
        table = self.models_table
        table.add_model(self.model_name, None, None, None, None, None)

        # Save the table
        table.saveto(self.models_table_path)

# -----------------------------------------------------------------
