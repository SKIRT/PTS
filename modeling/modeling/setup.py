#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.setup Contains the ModelingSetupTool class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from ...magic.tools.catalogs import get_ngc_name, get_hyperleda_name
from ...core.tools import filesystem as fs
from ..component.component import get_config_file_path
from ..component.sed import get_observed_sed_file_path, get_ski_template_path
from ...core.basics.configuration import Configuration, ConfigurationDefinition, InteractiveConfigurationSetter
from .galaxy import modeling_methods
from ...core.remote.host import find_host_ids
from ...core.data.sed import ObservedSED
from ...core.simulation.skifile import LabeledSkiFile

# -----------------------------------------------------------------

class ModelingSetupTool(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(ModelingSetupTool, self).__init__(config)

        # The path to the modeling directory
        self.modeling_path = None

        # The NGC name of the galaxy
        self.ngc_name = None

        # The HYPERLEDA name of the galaxy
        self.hyperleda_name = None

        # The configuration for the object depending on the specific type or modeling
        self.object_config = None

        # The modeling configuration
        self.modeling_config = None

        # The observed SED
        self.sed = None

        # The ski template
        self.ski = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the modeling directory
        self.create_directory()

        # 3. Get options
        if self.config.type == "galaxy": self.set_galaxy_options()
        elif self.config.type == "other": self.set_other_options()
        else: raise ValueError("Invalid option for 'type': " + self.config.type)

        # 4. load the SED
        if self.config.type == "other": self.load_sed()

        # 5. Load the ski template
        if self.config.type == "other": self.load_ski()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ModelingSetupTool, self).setup(**kwargs)

        # Set the path to the modeling directory
        self.modeling_path = fs.join(self.config.path, self.config.name)

        # Initialize the modeling configuration
        self.modeling_config = Configuration()

        # Get kwargs
        if "object_config" in kwargs: self.object_config = kwargs.pop("object_config")
        if "modeling_config" in kwargs: self.modeling_config = kwargs.pop("modeling_config")


    # -----------------------------------------------------------------

    def create_directory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the modeling directory ...")

        # Check whether a directory with this name is already present
        if fs.is_directory(self.modeling_path): raise ValueError("A directory with the name '" + self.config.name + "' already exists in the current working directory (" + self.modeling_path + ")")

        # Create the directory
        fs.create_directory(self.modeling_path)

    # -----------------------------------------------------------------

    def set_galaxy_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting options for the galaxy modeling ...")

        # Resolve the name of the galaxy
        self.resolve_name()

        # Prompt for galaxy settings
        self.prompt_galaxy()

        # Create configuration for galaxy modeling
        self.create_galaxy_config()

    # -----------------------------------------------------------------

    def resolve_name(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Resolving the galaxy name ...")

        # Get the NGC name of the galaxy
        self.ngc_name = get_ngc_name(self.config.name)

        # Inform the user
        log.info("Galaxy NGC ID is '" + self.ngc_name + "'")

        # Get the name in the HYPERLEDA catalog
        self.hyperleda_name = get_hyperleda_name(self.config.name)

        # Inform the user
        log.info("Galaxy HYPERLEDA ID is '" + self.hyperleda_name + "'")

    # -----------------------------------------------------------------

    def prompt_galaxy(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for options relevant for galaxy modeling ...")

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_required("host_ids", "string_list", "remote hosts to use for heavy computations (in order of preference)", choices=find_host_ids(schedulers=False))
        definition.add_required("method", "string", "method to use for the modeling", choices=modeling_methods)

        # Create configuration setter
        setter = InteractiveConfigurationSetter("galaxy modeling", "options for 3D modeling of a galaxy", add_cwd=False, add_logging=False)

        # Create the object config
        self.object_config = setter.run(definition, prompt_optional=False)

    # -----------------------------------------------------------------

    def create_galaxy_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the modeling configuration ...")

        # Set the configuration settings
        self.modeling_config.name = self.config.name
        self.modeling_config.modeling_type = self.config.type
        self.modeling_config.ngc_name = self.ngc_name
        self.modeling_config.hyperleda_name = self.hyperleda_name
        self.modeling_config.method = self.object_config.method
        self.modeling_config.host_ids = self.object_config.host_ids
        self.modeling_config.fitting_host_ids = self.config.fitting_host_ids

    # -----------------------------------------------------------------

    def set_other_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting options for modeling the SED of an object ...")

        # 1. Prompt for options
        self.prompt_other()

        # 2. Create the configuration
        self.create_other_config()

    # -----------------------------------------------------------------

    def prompt_other(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for options relevant for SED modeling ...")

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_required("ski", "file_path", "path/name of the template ski file")
        definition.add_flag("use_sed_file", "import an SED file produced with PTS (instead of manually entering the flux points)", False)

        # Create configuration setter
        setter = InteractiveConfigurationSetter("SED modeling", "options for SED modeling of an object", add_cwd=False, add_logging=False)

        # Create the object config
        self.object_config = setter.run(definition, prompt_optional=True)

    # -----------------------------------------------------------------

    def load_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SED ...")

        # Use SED file as input
        if self.object_config.use_sed_file:

            # Create definition
            definition = ConfigurationDefinition()
            definition.add_required("sed_path", "file_path", "path/name of the SED file")

            # Prompt for the SED path
            setter = InteractiveConfigurationSetter("SED", "name/path of the input SED", add_cwd=False, add_logging=False)
            config = setter.run(definition, prompt_optional=False)

            # Load the sed
            self.sed = ObservedSED.from_file(config.sed_path)

        # Manually enter flux points
        else:

            # Create definition
            definition = ConfigurationDefinition()
            definition.add_required("flux_points", "sed_entry_list", "flux points", dynamic_list=True)

            # Prompt for the flux points
            setter = InteractiveConfigurationSetter("Fluxes", "flux points for the SED", add_cwd=False, add_logging=False)
            config = setter.run(definition, prompt_optional=False)

            # Create new observed SED
            self.sed = ObservedSED(photometry_unit="Jy")

            # Add the flux points
            for fltr, flux, error in config.flux_points: self.sed.add_point(fltr, flux, error)

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski template ...")

        # Load the ski file
        self.ski = LabeledSkiFile(self.object_config.ski)

    # -----------------------------------------------------------------

    def create_other_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the modeling configuration ...")

        # Set the settings
        self.modeling_config.name = self.config.name
        self.modeling_config.modeling_type = self.config.type
        self.modeling_config.fitting_host_ids = self.config.fitting_host_ids

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the configuration
        self.write_config()

        # Write the SED
        if self.config.type == "other": self.write_sed()

        # Write the ski template
        if self.config.type == "other": self.write_ski()

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the modeling configuration ...")

        # Determine the path
        path = get_config_file_path(self.modeling_path)

        # Save the config
        self.modeling_config.saveto(path)

    # -----------------------------------------------------------------

    def write_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the input SED ...")

        # Determine the path
        path = get_observed_sed_file_path(self.modeling_path)

        # Save the SED
        self.sed.saveto(path)

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski template ...")

        # Determine the path
        path = get_ski_template_path(self.modeling_path)

        # Save the ski template
        self.ski.saveto(path)

# -----------------------------------------------------------------
