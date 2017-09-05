#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.models.dust Contains the DustBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.basics.configuration import prompt_proceed, ConfigurationDefinition, InteractiveConfigurationSetter, prompt_string, prompt_yn, prompt_filepath, PassiveConfigurationSetter
from ....core.basics.log import log
from .general import GeneralBuilder
from ..suite import model_map_filename
from ....core.tools import filesystem as fs
from ....magic.core.frame import Frame
from ...component.galaxy import GalaxyModelingComponent
from ....core.tools.utils import lazyproperty

# -----------------------------------------------------------------

basic_dust_map_name = "dust_disk"
basic_dust_map_names = [basic_dust_map_name]

# -----------------------------------------------------------------

# Define titles for the fixed components
titles = dict()
titles["disk"] = "Dust disk"

# -----------------------------------------------------------------

component_name_for_map_name = dict()
component_name_for_map_name["dust_disk"] = "disk"

# -----------------------------------------------------------------

class DustBuilder(GeneralBuilder, GalaxyModelingComponent):
    
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
        #super(DustBuilder, self).__init__(*args, **kwargs)
        GeneralBuilder.__init__(self, no_config=True)
        GalaxyModelingComponent.__init__(self, *args, **kwargs)

        # The scaleheight of the old stellar population
        self.old_scaleheight = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Build dust
        if self.config.disk: self.build_dust_disk()

        # Additional
        if self.config.additional: self.build_additional()

        # Write
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
        #super(DustBuilder, self).setup()
        GeneralBuilder.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

        # Get the scaleheight of the old stars
        self.old_scaleheight = kwargs.pop("old_scaleheight")

    # -----------------------------------------------------------------

    def build_dust_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the dust disk component ...")

        # 1. Get parameters
        self.get_dust_disk_parameters()

        # 2. Load the map
        self.load_dust_disk_map()

        # 3. Create the deprojection model
        self.create_deprojection_disk()

    # -----------------------------------------------------------------

    #@lazyproperty
    #def old_scaleheight(self):

        #"""
        #This function ...
        #:return:
        #"""

        #
        #definition = self.get_model_definition(self.model_name)
        #return definition.old_stars_scaleheight
        # NO: don't load the definition because e.g. the property old_stars_scaleheight depends on the models table to be completed, and it isn't at this point

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_scaleheight(self):

        """
        This fucntion ...
        :return: 
        """

        return self.config.dust_scaleheight_ratio * self.old_scaleheight

    # -----------------------------------------------------------------

    def get_dust_disk_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the dust component ...")

        # Check defaults
        if self.config.default_dust_mass is None: raise ValueError("Default dust mass cannot be undefined")
        if self.config.default_hydrocarbon_pops is None: raise ValueError("Default number of hydrocarbon populations cannot be undefined")
        if self.config.default_enstatite_pops is None: raise ValueError("Default number of enstatite populations cannot be undefined")
        if self.config.default_forsterite_pops is None: raise ValueError("Default number of forsterite populations cannot be undefined")

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_optional("scale_height", "quantity", "scale height", default=self.dust_scaleheight)
        definition.add_optional("mass", "quantity", "dust mass", default=self.config.default_dust_mass)
        definition.add_optional("hydrocarbon_pops", "positive_integer", "number of hydrocarbon populations", default=self.config.default_hydrocarbon_pops)
        definition.add_optional("enstatite_pops", "positive_integer", "number of enstatite populations", default=self.config.default_enstatite_pops)
        definition.add_optional("forsterite_pops", "positive_integer", "number of forsterite populations", default=self.config.default_forsterite_pops)

        # Use default settings
        if self.config.use_defaults:

            setter = PassiveConfigurationSetter("dust disk", add_logging=False, add_cwd=False)
            config = setter.run(definition)

        # Prompt for the settings
        else:

            # Prompt for settings
            setter = InteractiveConfigurationSetter("dust disk", add_logging=False, add_cwd=False)
            config = setter.run(definition, prompt_optional=True)

        # Set the title
        config.title = titles["disk"]

        # Set the parameters
        self.parameters["disk"] = config

    # -----------------------------------------------------------------

    def load_dust_disk_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the dust disk map ...")

        # Ask whether a custom map should be used
        custom = prompt_yn("custom", "use a custom map for the dust disk (instead of one of those created in the modelling pipeline)", default=False)

        # Load a custom dust disk map
        if custom: self.load_custom_dust_disk_map()

        # Load a dust disk map from modeling pipeline
        else: self.load_modeling_dust_disk_map()

    # -----------------------------------------------------------------

    def load_custom_dust_disk_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading a custom dust disk map ...")

        # Get the path
        path = prompt_filepath("filepath", "dust disk map path")

        # Load the map
        self.maps["disk"] = Frame.from_file(path)

    # -----------------------------------------------------------------

    def load_modeling_dust_disk_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading a dust disk map from the modeling pipeline ...")

        # Get the map names
        names = self.static_maps_selection.dust_map_names

        # Ask for the dust map to use
        name = prompt_string("dust_map", "dust disk map to use for this model", choices=names)

        # Set the path
        filepath = self.static_maps_selection.dust_map_paths[name]

        # Set the map
        self.maps["disk"] = Frame.from_file(filepath)

    # -----------------------------------------------------------------

    def create_deprojection_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the deprojection model for the dust disk ...")

        # Create the deprojection model
        deprojection = self.create_deprojection_for_map(self.galaxy_properties, self.disk_position_angle, self.maps["disk"], model_map_filename, self.parameters["disk"].scale_height)

        # Set the deprojection model
        self.deprojections["disk"] = deprojection

    # -----------------------------------------------------------------

    def build_additional(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building additional dust components ...")

        # Proceed?
        while prompt_proceed():

            # Set parameters
            name = self.set_additional_parameters()

            # Set properties
            normalization_parameters = self.set_normalization_properties(name)

            # Set geometry properties
            geometry_parameters = self.set_geometry_properties(name)

            # Set mix properties
            mix_parameters = self.set_mix_properties(name)

            # Create the properties
            properties = {"normalization": normalization_parameters, "geometry": geometry_parameters, "mix": mix_parameters}

            # Load map if necessary
            if self.parameters[name].geometry == "ReadFitsGeometry": self.load_additional_map(name, geometry_parameters)

            # Add the properties
            self.properties[name] = properties

    # -----------------------------------------------------------------

    def set_additional_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring an additional dust component ...")

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_required("name", "string", "name for this dust component")
        definition.add_optional("title", "string", "short description for this component")
        definition.add_optional("geometry", "string", "SKIRT base geometry for the component", self.smile.concrete_geometries)
        definition.add_optional("normalization", "string", "normalization for the component", self.smile.concrete_dust_normalizations)
        definition.add_optional("mix", "string", "dust mix for the component", self.smile.concrete_dust_mixes)

        # Prompt for settings
        setter = InteractiveConfigurationSetter("additional dust component", add_cwd=False, add_logging=False)
        config = setter.run(definition, prompt_optional=True)

        # Check the name
        if config.name in self.parameters: raise ValueError("You cannot use this name for the dust component: already in use")

        # Set the parameters
        self.parameters[config.name] = config

        # Return the name of the new dust component
        return config.name

    # -----------------------------------------------------------------

    def set_normalization_properties(self, name):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the normalization of dust component '" + name + "' ...")

        # Get the selected type of normalization
        normalization_type = self.parameters[name].normalization

        # Get parameters for this simulation item
        parameters = self.smile.prompt_parameters_for_type(normalization_type, merge=True)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def set_geometry_properties(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.info("Configuring the geometry of dust component '" + name + "' ...")

        # Get the selected type of geometry
        geometry_type = self.parameters[name].geometry

        # Get parameters for this simulation item
        parameters = self.smile.prompt_parameters_for_type(geometry_type, merge=True)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def set_mix_properties(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.info("Configuring the dust mix of dust component '" + name + "' ...")

        # Get the selected type of dust mix
        mix_type = self.parameters[name].mix

        # Get parameters for this simulation item
        parameters = self.smile.prompt_parameters_for_type(mix_type, merge=True)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def load_additional_map(self, name, parameters):

        """
        This function ...
        :param name:
        :param parameters:
        :return:
        """

        # Inform the user
        log.info("Loading the input map for stellar component '" + name + "' ...")

        # Get the absolute file path to the map
        path = fs.absolute_path(parameters["filename"])

        # Load the map
        self.maps[name] = Frame.from_file(path)

        # Change the filename in the geometry parameters
        parameters["filename"] = model_map_filename

    # -----------------------------------------------------------------

    @property
    def dust_map_path(self):

        """
        This function ...
        :return: 
        """

        component_name = component_name_for_map_name[basic_dust_map_name]
        return self.map_paths[component_name]

# -----------------------------------------------------------------
