#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.composer Contains the ModelComposer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.configurable import InteractiveConfigurable
from .smile import SKIRTSmileSchema
from ..simulation.skifile import SkiFile
from ..tools.utils import lazyproperty
from ..basics.configuration import prompt_choice

# -----------------------------------------------------------------

oligo_type = "oligo"
pan_type = "pan"
oligo_or_pan = [oligo_type, pan_type]

# -----------------------------------------------------------------

commands = OrderedDict()

# -----------------------------------------------------------------

subcommands = OrderedDict()

# -----------------------------------------------------------------

class ModelComposer(InteractiveConfigurable):

    """
    This function ...
    """

    _commands = commands
    _subcommands = subcommands

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ModelComposer, self).__init__(*args, **kwargs)

        # The SKIRT smile schema creator
        self.smile = SKIRTSmileSchema()

        # The skifile
        self.ski = None

    # -----------------------------------------------------------------

    @property
    def do_commands(self):

        """
        This function ...
        :return:
        """

        return self.config.commands is not None and len(self.config.commands) > 0

    # -----------------------------------------------------------------

    @property
    def do_interactive(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_showing(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.show

    # -----------------------------------------------------------------

    @property
    def do_writing(self):

        """
        This function ...
        :return:
        """

        return self.config.write

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        Thisf unction ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Run commands
        if self.do_commands: self.run_commands()

        # 3. Interactive
        if self.do_interactive: self.interactive()

        # 6. Show
        if self.do_showing: self.show()

        # 7. Write
        if self.do_writing: self.write()

        # 8. Plotting
        #if self.do_plotting: self.plot()

        # 12. Write the history
        if self.has_commands: self.write_history()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ModelComposer, self).setup(**kwargs)

        # Get the ski file
        if kwargs.get("skifile", None) is not None: self.ski = kwargs.pop("skifile")
        elif kwargs.get("ski", None) is not None: self.ski = kwargs.pop("ski")
        elif self.config.from_file is not None: self.ski = SkiFile(self.config.skifile)
        else: self.create_ski_template()

    # -----------------------------------------------------------------

    @lazyproperty
    def name(self):

        """
        This function ...
        :return:
        """

        if self.config.name is not None: return self.config.name
        else: return self.ski.prefix

    # -----------------------------------------------------------------

    def create_ski_template(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating empty ski template ...")

        # Oligo or pan?
        if self.config.type is None: simulation_type = prompt_choice("type", "simulation type", "")
        else: simulation_type = self.config.type

        # Create template ski file
        if simulation_type == oligo_type: self.ski = self.smile.create_oligochromatic_template()
        elif simulation_type == pan_type: self.ski =self.smile.create_panchromatic_template()
        else: raise ValueError("Invalid value for 'simulation_type'")

    # -----------------------------------------------------------------

    @property
    def simulation_type(self):

        """
        This function ...
        :return:
        """

        return self.ski.simulation_type

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Skifile
        self.write_ski()

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

    # -----------------------------------------------------------------

    @property
    def history_filename(self):

        """
        This function ...
        :return:
        """

        return "composer"

# -----------------------------------------------------------------
