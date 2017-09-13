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

        # Map paths
        self.old_stars_map_path = None
        self.young_stars_map_path = None
        self.ionizing_stars_map_path = None
        self.dust_map_path = None

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


    # -----------------------------------------------------------------

    def adjust_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting dust components from previous model [" + self.previous_name + "] ...")

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
