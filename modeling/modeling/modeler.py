#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeling.modeler Contains the Modeler class, which selects the appropriate modeler and runs it.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ..component.component import get_config_file_path, load_modeling_configuration
from .galaxy import GalaxyModeler
from .sed import SEDModeler
from .images import ImagesModeler
from ..core.steps import galaxy_modeling, sed_modeling, images_modeling

# -----------------------------------------------------------------

class Modeler(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(Modeler, self).__init__(*args, **kwargs)

        # The modeling path
        self.modeling_path = None

        # The modeling configuration
        self.modeling_config = None

        # The specific modeler
        self.modeler = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Perform the modeling
        self.model(**kwargs)

    # -----------------------------------------------------------------

    @property
    def galaxy_modeling(self):

        """
        This function ...
        :return:
        """

        return self.modeling_config.modeling_type == galaxy_modeling

    # -----------------------------------------------------------------

    @property
    def sed_modeling(self):

        """
        This function ...
        :return:
        """

        return self.modeling_config.modeling_type == sed_modeling

    # -----------------------------------------------------------------

    @property
    def images_modeling(self):

        """
        This function ...
        :return:
        """

        return self.modeling_config.modeling_type == images_modeling

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup fucntion of the base class
        super(Modeler, self).setup(**kwargs)

        # Set the path to the modeling directory
        self.modeling_path = self.config.path

        # Check for the presence of the configuration file
        if not fs.is_file(get_config_file_path(self.modeling_path)): raise ValueError("The current working directory (" + self.config.path + ") is not a radiative transfer modeling directory (the configuration file is missing)")
        else: self.modeling_config = load_modeling_configuration(self.modeling_path)

    # -----------------------------------------------------------------

    def model(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the modeling ...")

        # Debugging
        log.debug("Modeling type: " + self.modeling_config.modeling_type)

        # Galaxy modeling
        if self.galaxy_modeling: self.model_galaxy(**kwargs)

        # SED modeling
        elif self.sed_modeling: self.model_sed(**kwargs)

        # Images modeling
        elif self.images_modeling: self.model_images(**kwargs)

    # -----------------------------------------------------------------

    def model_galaxy(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Create galaxy modeler
        self.modeler = GalaxyModeler(self.config)

        # Run the modeler
        self.modeler.run(**kwargs)

    # -----------------------------------------------------------------

    def model_sed(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Create SED modeler
        self.modeler = SEDModeler(self.config)

        # Run the modeler
        self.modeler.run(**kwargs)

    # -----------------------------------------------------------------

    def model_images(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Create images modeler
        self.modeler = ImagesModeler(self.config)

        # Run the modeler
        self.modeler.run(**kwargs)

# -----------------------------------------------------------------
