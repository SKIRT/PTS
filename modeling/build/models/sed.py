#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.models.sed Contains the SEDModelBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ...component.sed import get_ski_template, get_ski_input_path
from ....core.tools.serialization import write_dict
from ....magic.core.frame import Frame
from ..suite import model_map_filename
from .base import ModelBuilderBase

# -----------------------------------------------------------------

class SEDModelBuilder(ModelBuilderBase):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SEDModelBuilder, self).__init__(*args, **kwargs)

        # The path for the other input files
        self.model_input_path = None

        # The ski file template
        self.ski = None

        # The stellar properties
        self.stellar_properties = dict()

        # The dust properties
        self.dust_properties = dict()

        # The paths to the component directories
        self.stellar_paths = dict()
        self.dust_paths = dict()

        # The stellar maps
        self.stellar_maps = dict()

        # The dust maps
        self.dust_maps = dict()

        # The original input map filenames
        self.original_map_filenames = []

        # Other input paths
        self.other_input = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the ski file template
        self.load_ski()

        # 3. Get the stellar components
        self.get_stellar_components()

        # 4. Get the dust components
        self.get_dust_components()

        # 5. Load input other than the input maps defining
        if self.ski.needs_input: self.load_other_input()

        # 6. Write
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

        # Set the path of the input directory
        self.model_input_path = fs.create_directory_in(self.model_path, "input")

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

            # Check the geometry of the stellar component
            geometries = properties["children"]["geometry"]["children"]
            if len(geometries) != 1: raise ValueError("Cannot locate the geometry for stellar component '" + component_id + "'")
            geometry_type = geometries.keys()[0]
            geometry_properties = geometries[geometry_type]
            if geometry_type == "ReadFitsGeometry": self.load_stellar_map(component_id, geometry_properties)

            # Set the properties
            self.stellar_properties[component_id] = properties

    # -----------------------------------------------------------------

    def load_stellar_map(self, name, parameters):

        """
        This function ...
        :param name:
        :param parameters:
        :return:
        """

        # Inform the user
        log.info("Loading the input map for the '" + name + "' stellar component ...")

        # Add to original filenames
        self.original_map_filenames.append(parameters["filename"])

        # Get the absolute file path
        path = fs.join(get_ski_input_path(self.config.path), parameters["filename"])

        # Debugging
        log.debug("Loading the map from '" + path + "' ...")

        # Load the map
        self.stellar_maps[name] = Frame.from_file(path)

        # Change the filename in the geometry parameters
        parameters["filename"] = model_map_filename

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

            # Check the geometry of the dust component
            geometries = properties["children"]["geometry"]["children"]
            if len(geometries) != 1: raise ValueError("Cannot locate the geometry for dust component '" + component_id + "'")
            geometry_type = geometries.keys()[0]
            geometry_properties = geometries[geometry_type]
            if geometry_type == "ReadFitsGeometry": self.load_dust_map(component_id, geometry_properties)

            # Set the properties
            self.dust_properties[component_id] = properties

    # -----------------------------------------------------------------

    def load_dust_map(self, name, parameters):

        """
        This function ...
        :param name:
        :param parameters:
        :return:
        """

        # Inform the user
        log.info("Loading the input map for the '" + name + "' dust component ...")

        # Add to original filenames
        self.original_map_filenames.append(parameters["filename"])

        # Get the absolute file path
        path = fs.join(get_ski_input_path(self.config.path), parameters["filename"])

        # Debugging
        log.debug("Loading the map from '" + path + "' ...")

        # Load the map
        self.dust_maps[name] = Frame.from_file(path)

        # Change the filename in the geometry parameters
        parameters["filename"] = model_map_filename

    # -----------------------------------------------------------------

    def load_other_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading other input files for the ski file (expect wavelength grid files) ...")

        # Get the wavelength grid filename (if any)
        wavelength_grid_filename = self.ski.wavelengthsfile()

        # Loop over the input filenames
        for filename in self.ski.input_files:

            # Skip wavelength grid file
            if wavelength_grid_filename is not None and filename == wavelength_grid_filename: continue

            # Skip maps
            if filename in self.original_map_filenames: continue

            # Add others
            self.other_input[filename] = fs.join(get_ski_input_path(self.config.path), filename)

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

        # Write stellar maps
        self.write_stellar_maps()

        # Write dust maps
        self.write_dust_maps()

        # Write other input
        self.write_other_input()

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
            component_path = fs.join(self.model_stellar_path, str(name)) # name (component ID can be an integer)
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
            component_path = fs.join(self.model_dust_path, str(name)) # name (component ID) can be an integer
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

    def write_stellar_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing stellar maps ...")

        # Loop over the components with a map
        for name in self.stellar_maps:

            # Determine path and save
            path = fs.join(self.stellar_paths[name], model_map_filename)
            self.stellar_maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing dust maps ...")

        # Loop over the components with a map
        for name in self.dust_maps:

            # Determine path and save
            path = fs.join(self.dust_paths[name], model_map_filename)
            self.dust_maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_other_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing other input files ...")

        # Loop over the files
        for filename in self.other_input:

            # Original filepath
            filepath = self.other_input[filename]

            # Copy
            fs.copy_file(filepath, self.model_input_path)

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
