#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.models.galaxy Contains the GalaxyModelBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .dust import DustBuilder
from .stars import StarsBuilder
from ....core.basics.log import log
from ...component.galaxy import GalaxyModelingComponent
from .base import ModelBuilderBase
from ....core.tools.utils import lazyproperty
from ....core.basics.configuration import prompt_yn, prompt_automatic, prompt_variable
from .stars import old_component_name, young_component_name, ionizing_component_name
from .dust import disk_component_name
from .general import write_component

# -----------------------------------------------------------------

class GalaxyModelBuilder(ModelBuilderBase, GalaxyModelingComponent):
    
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
        #super(ModelBuilder, self).__init__(*args, **kwargs)
        ModelBuilderBase.__init__(self, no_config=True)
        GalaxyModelingComponent.__init__(self, *args, **kwargs)

        # The scaleheight of the old stars
        self.old_scaleheight = None

        # Model component paths
        self.bulge_path = None
        self.old_path = None
        self.young_path = None
        self.ionizing_path = None
        self.dust_path = None
        self.extra_stellar_paths = []
        self.extra_dust_paths = []

        # Map paths
        self.old_stars_map_path = None
        self.young_stars_map_path = None
        self.ionizing_stars_map_path = None
        self.dust_map_path = None
        self.extra_stellar_map_paths = []
        self.extra_dust_map_paths = []

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Adjust stellar components from previous model
        if self.from_previous: self.adjust_stars()

        # Adjust dust components from previous model
        if self.from_previous: self.adjust_dust()

        exit()

        # 2. Build stars
        self.build_stars()

        # 3. Build dust component
        self.build_dust()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        #super(ModelBuilder, self).setup(**kwargs)
        ModelBuilderBase.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

    # -----------------------------------------------------------------

    @property
    def from_previous(self):

        """
        This function ...
        :return:
        """

        return self.config.from_previous is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def previous(self):

        """
        This function ...
        :return:
        """

        if not self.from_previous: return None
        else: return self.suite.get_model_definition(self.config.from_previous)

    # -----------------------------------------------------------------

    @property
    def previous_name(self):

        """
        This function ...
        :return:
        """

        return self.previous.name

    # -----------------------------------------------------------------

    def adjust_stars(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting stellar components from previous model [" + self.previous_name + "] ...")

        # Loop over the stellar components
        for name in self.previous.stellar_component_names:

            # Include?
            if not prompt_yn("include_" + name, "include the '" + name + "' stellar component of the '" + self.previous_name + "' model", default=True): continue

            # Adjust?
            adjust = prompt_yn("adjust_" + name, "adjust the properties of the '" + name + "' stellar component of the '" + self.previous_name + "' model definition")

            # Adjust the component
            if adjust: self.adjust_stellar_component(name)

            # Don't adjust the component
            else: self.set_previous_stellar_component(name)

    # -----------------------------------------------------------------

    def adjust_stellar_component(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.info("Adjusting the '" + name + "' stellar component ...")

        # Load the component
        component = self.previous.load_stellar_component(name, add_map=False)

        # Adjust
        adjusted = adjust_component(component)

        print(adjusted)

        # Set reference to the previous model for this dust component
        if not adjusted: self.set_previous_stellar_component(name)
        else:

            write_component(component)

    # -----------------------------------------------------------------

    def set_previous_stellar_component(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.info("Setting reference to the previous model for the '" + name + "' stellar component ...")

        # Get the component path
        path = self.previous.get_stellar_component_path(name)

        # Set the appropriate attribute
        if name == old_component_name: self.old_path = path
        elif name == young_component_name: self.young_path = path
        elif name == ionizing_component_name: self.ionizing_path = path
        else: self.extra_stellar_paths.append(path)

    # -----------------------------------------------------------------

    def adjust_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting dust components from previous model [" + self.previous_name + "] ...")

        # Loop over the dust components
        for name in self.previous.dust_component_names:

            # Include?
            if not prompt_yn("include_" + name, "include the '" + name + "' dust component of the '" + self.previous_name + "' model", default=True): continue

            # Adjust?
            adjust = prompt_yn("adjust_" + name, "adjust the '" + name + "' dust component of the '" + self.previous_name + "' model definition")

            # Adjust the component
            if adjust: self.adjust_dust_component(name)

            # Don't adjust
            else: self.set_previous_dust_component(name)

    # -----------------------------------------------------------------

    def adjust_dust_component(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.info("Adjusting the '" + name + "' dust component ...")

        # Load the component
        component = self.previous.load_dust_component(name, add_map=False)

        # Adjust
        adjusted = adjust_component(component)

        print(adjusted)

        # Set reference to the previous model for this dust component
        if not adjusted: self.set_previous_dust_component(name)
        else:

            write_component(component)

    # -----------------------------------------------------------------

    def set_previous_dust_component(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the component path
        path = self.previous.get_dust_component_path(name)

        # Set the appropriate attribute
        if name == disk_component_name: self.dust_path = path
        else: self.extra_dust_paths.append(path)

    # -----------------------------------------------------------------

    @property
    def has_bulge(self):

        """
        This function ...
        :return:
        """

        return self.bulge_path is not None

    # -----------------------------------------------------------------

    @property
    def has_old(self):

        """
        This function ...
        :return:
        """

        return self.old_path is not None

    # -----------------------------------------------------------------

    @property
    def has_young(self):

        """
        This function ...
        :return:
        """

        return self.young_path is not None

    # -----------------------------------------------------------------

    @property
    def has_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_path is not None

    # -----------------------------------------------------------------

    def build_stars(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the stellar components ...")

        # Create configuration
        config = dict()
        config["name"] = self.model_name
        config["output"] = self.model_stellar_path
        config["default_sfr"] = self.config.sfr

        # Set options to build different components
        config["bulge"] = not self.has_bulge
        config["old"] = not self.has_old
        config["young"] = not self.has_young
        config["ionizing"] = not self.has_ionizing
        config["additional"] = self.config.additional

        # Create the builder
        builder = StarsBuilder(interactive=True, cwd=self.config.path, config=config, prompt_optional=True)

        # Run
        builder.run()

        # Set the scaleheight of the old stars
        self.old_scaleheight = builder.old_scaleheight

        # Set map paths
        self.old_stars_map_path = builder.old_stars_map_path
        self.young_stars_map_path = builder.young_stars_map_path
        self.ionizing_stars_map_path = builder.ionizing_stars_map_path

    # -----------------------------------------------------------------

    @property
    def has_dust(self):

        """
        This function ...
        :return:
        """

        return self.dust_path is not None

    # -----------------------------------------------------------------

    def build_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the dust components ...")

        # Create configuration
        config = dict()
        config["name"] = self.model_name
        config["output"] = self.model_dust_path
        config["default_dust_mass"] = self.config.dust_mass

        # Set options to build different components
        config["disk"] = not self.has_dust
        config["additional"] = self.config.additional

        # Create the builder
        builder = DustBuilder(interactive=True, cwd=self.config.path, config=config, prompt_optional=True)

        # Run
        builder.run(old_scaleheight=self.old_scaleheight)

        # Set map paths
        self.dust_map_path = builder.dust_map_path

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write models table
        self.write_table()

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the models table ...")

        # Add the model
        table = self.models_table
        table.add_model(self.model_name, self.config.description, self.old_stars_map_path, self.young_stars_map_path, self.ionizing_stars_map_path, self.dust_map_path)

        # Save the table
        table.saveto(self.models_table_path)

# -----------------------------------------------------------------

def adjust_component(component):

    """
    Thisj function ...
    :param component:
    :return:
    """

    adjusted = False

    # Parameters
    if "parameters" in component:

        # Adjust?
        adjust = prompt_yn("adjust_parameters", "adjust the parameters", default=False)
        if adjust: adjusted_parameters = adjust_parameters(component.parameters)
        else: adjusted_parameters = False

        # Update adjusted flag
        if adjusted_parameters: adjusted = True

    # Deprojection
    if "deprojection" in component:

        # Adjust?
        adjust = prompt_yn("adjust_deprojection", "adjust the deprojection", default=False)
        if adjust: adjusted_deprojection = adjust_deprojection(component.deprojection)
        else: adjusted_deprojection = False

        # Update adjusted flag
        if adjusted_deprojection: adjusted = True

    # Model
    if "model" in component:

        # Adjust?
        adjust = prompt_yn("adjust_model", "adjust the model", default=False)
        if adjust: adjusted_model = adjust_model(component.model)
        else: adjusted_model = False

        # Update adjusted flag
        if adjusted_model: adjusted = True

    # Properties
    if "properties" in component:

        # Adjust?
        adjust = prompt_yn("adjust_properties", "adjust the properties", default=False)
        if adjust: adjusted_properties = adjust_properties(component.properties)
        else: adjusted_properties = False

        # Update adjusted flag
        if adjusted_properties: adjusted = True

    # Return adjusted flag
    return adjusted

# -----------------------------------------------------------------

def adjust_parameters(parameters):

    """
    This function ...
    :param parameters:
    :return:
    """

    adjusted = False

    # Loop over the parameters
    for name in parameters:

        default = parameters[name]
        value = prompt_automatic(name, "adjusted '" + name + "' value", default)
        if value != default: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------

def adjust_deprojection(deprojection):

    """
    This function ...
    :param deprojection:
    :return:
    """

    adjusted = False

    # Loop over the properties
    for name in deprojection:

        default = deprojection[name]
        ptype = deprojection.get_ptype(name)
        value = prompt_variable(name, ptype, "adjusted '" + name + "' value", default=default)
        if value != default: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------

def adjust_model(model):

    """
    This function ...
    :param model:
    :return:
    """

    adjusted = False

    # Loop over the properties
    for name in model:

        default = model[name]
        ptype = model.get_ptype(name)
        value = prompt_variable(name, ptype, "adjusted '" + name + "' value", default=default)
        if value != default: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------

def adjust_properties(properties):

    """
    This function ...
    :param properties:
    :return:
    """

    adjusted = False

    adjust_geometry = prompt_yn("adjust_geometry", "adjust geometry")
    if adjust_geometry: adjusted_geometry = adjust_geometry_properties(properties["geometry"])
    else: adjusted_geometry = False

    # Set flag
    if adjusted_geometry: adjusted = True

    adjust_sed = prompt_yn("adjust_sed", "adjust SED properties")
    if adjust_sed: adjusted_sed = adjust_sed_properties(properties["sed"])
    else: adjusted_sed = False

    # Set flag
    if adjusted_sed: adjusted = True

    adjust_normalization = prompt_yn("adjust_normalization", "adjust normalization properties")
    if adjust_normalization: adjusted_normalization = adjust_normalization_properties(properties["normalization"])
    else: adjusted_normalization = False

    # Set flag
    if adjusted_normalization: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------

def adjust_geometry_properties(properties):

    """
    This function ...
    :param properties:
    :return:
    """

    adjusted = False

    # Loop over the properties
    for name in properties:

        default = properties[name]
        value = prompt_automatic(name, "adjusted '" + name + "' value", default)
        if value != default: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------

def adjust_sed_properties(properties):

    """
    This function ...
    :param properties:
    :return:
    """

    adjusted = False

    # Loop over the properties
    for name in properties:

        default = properties[name]
        value = prompt_automatic(name, "adjusted '" + name + "' value", default)
        if value != default: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------

def adjust_normalization_properties(properties):

    """
    This function ...
    :param properties:
    :return:
    """

    adjusted = False

    # Loop over the properties
    for name in properties:

        default = properties[name]
        value = prompt_automatic(name, "adjusted '" + name + "' value", default)
        if value != default: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------
