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
from .component import BuildComponent
from ...core.basics.configuration import prompt_proceed, ConfigurationDefinition, InteractiveConfigurationSetter
from ...core.tools.logging import log
from ...core.basics.unit import parse_unit as u
from ...core.prep.smile import SKIRTSmileSchema

# -----------------------------------------------------------------

class DustBuilder(BuildComponent):
    
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

        # The parameters
        self.parameters = dict()

        # The components
        self.components = dict()

        # The SKIRT smile schema
        self.smile = None

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

        # Get parameters
        self.get_dust_disk_parameters()

    # -----------------------------------------------------------------

    def get_dust_disk_parameters(self):

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

    def build_additional(self):

        """
        This function ...
        :return:
        """

        # Proceed?
        while prompt_proceed():

            definition = ConfigurationDefinition()
            definition.add_required("name", "string", "name for the dust component")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
