#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.dust Contains the DustBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configuration import prompt_proceed, ConfigurationDefinition, InteractiveConfigurationSetter
from ...core.tools.logging import log
from ...core.basics.unit import parse_unit as u
from .general import GeneralBuilder

# -----------------------------------------------------------------

class DustBuilder(GeneralBuilder):
    
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
        super(DustBuilder, self).__init__(config, interactive)

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Build dust
        if self.config.disk: self.build_dust_disk()

        # Additional
        if self.config.additional: self.build_additional()

        # Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DustBuilder, self).setup()

    # -----------------------------------------------------------------

    def build_dust_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the dust disk component ...")

        # Get parameters
        self.get_dust_disk_parameters()

        # Load the map
        self.load_dust_disk_map()

        # Create the deprojection model
        self.create_deprojection_disk()

    # -----------------------------------------------------------------

    def get_dust_disk_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the dust component ...")

        # scale_height = 260.5 * Unit("pc") # first models
        scale_height = 200. * u("pc")  # M51
        dust_mass = 1.5e7 * u("Msun")

        hydrocarbon_pops = 25
        enstatite_pops = 25
        forsterite_pops = 25

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_optional("scale_height", "quantity", "scale height", default=scale_height)
        definition.add_optional("mass", "quantity", "dust mass", default=dust_mass)
        definition.add_optional("hydrocarbon_pops", "positive_integer", "number of hydrocarbon populations", default=hydrocarbon_pops)
        definition.add_optional("enstatite_pops", "positive_integer", "number of enstatite populations", default=enstatite_pops)
        definition.add_optional("forsterite_pops", "positive_integer", "number of forsterite populations", default=forsterite_pops)

        # Prompt for settings
        setter = InteractiveConfigurationSetter("dust disk")
        config = setter.run(definition)

        # Set the parameters
        self.parameters["disk"] = config

        ## NEW: SET FIXED PARAMETERS
        #self.fixed["dust_scaleheight"] = scale_height
        ##

        # Set the parameters of the dust component
        #deprojection = self.deprojection.copy()
        #deprojection.filename = self.dust_map_filename
        #deprojection.scale_height = scale_height
        #self.deprojections["dust"] = deprojection

        # Adjust the ski file
        #self.ski_template.set_dust_component_geometry(0, deprojection)
        #self.ski_template.set_dust_component_themis_mix(0, hydrocarbon_pops, enstatite_pops, forsterite_pops)  # dust mix
        # self.ski_template.set_dust_component_mass(0, dust_mass) # dust mass
        #self.ski_template.set_labeled_value("dust_mass", dust_mass)  # keep label

    # -----------------------------------------------------------------

    def load_dust_disk_map(self):

        """
        This function ...
        :return:
        """

        # Set the map
        self.maps["disk"] = None

    # -----------------------------------------------------------------

    def create_deprojection_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the deprojection model for the dust disk ...")

        # Create the deprojection model
        deprojection = self.create_deprojection_for_map(self.maps["disk"])

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
        definition.add_required("name", "string", "name for this stellar component")
        definition.add_optional("description", "string", "description for the component")
        definition.add_optional("geometry", "string", "SKIRT base geometry for the component", self.smile.concrete_geometries)
        definition.add_optional("normalization", "string", "normalization for the component", self.smile.concrete_dust_normalizations)
        definition.add_optional("mix", "string", "dust mix for the component", self.smile.concrete_dust_mixes)

        # Prompt for settings
        setter = InteractiveConfigurationSetter("additional dust component", add_cwd=False, add_logging=False)
        config = setter.run(definition)

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

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
