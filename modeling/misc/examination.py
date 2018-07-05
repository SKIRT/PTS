#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.misc.examination Contains the ModelExamination class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.configurable import InteractiveConfigurable, InvalidCommandError
from ...core.basics.configuration import ConfigurationDefinition
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

# Standard commands
_help_command_name = "help"
_history_command_name = "history"

# Other
_project_command_name = "project"
_parameters_command_name = "parameters"

# -----------------------------------------------------------------

# Define commands
commands = OrderedDict()

# Standard commands
commands[_help_command_name] = ("show_help", False, "show help", None)
commands[_history_command_name] = ("show_history_command", True, "show history of executed commands", None)

# Other
commands[_project_command_name] = ("project_command", True, "project the model from one or multiple orientations", None)
commands[_parameters_command_name] = ("show_parameters_command", True, "show model parameters", None)

# -----------------------------------------------------------------

# Subcommands
subcommands = OrderedDict()

# -----------------------------------------------------------------

all_name = "all"
intrinsic_name = "intrinsic"
derived_name = "derived"
stellar_name = "stellar"
bulge_name = "bulge"
disk_name = "disk"
old_name = "old"
young_name = "young"
sfr_name = "sfr"
unevolved_name = "unevolved"
dust_name = "dust"
parameter_categories = [all_name, intrinsic_name, derived_name, stellar_name, bulge_name, disk_name, old_name, young_name, sfr_name, unevolved_name, dust_name]

# -----------------------------------------------------------------

class ModelExamination(InteractiveConfigurable):

    """
    This class ...
    """

    _commands = commands
    _subcommands = subcommands
    _log_section = "MODEL"

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(ModelExamination, self).__init__(*args, **kwargs)

        # The model
        self.model = None

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

        return self.config.interactive

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Run commands
        if self.do_commands: self.run_commands()

        # 3. Interactive
        if self.do_interactive: self.interactive()

        # 4. Show
        self.show()

        # 5. Write the history
        if self.has_commands: self.write_history()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ModelExamination, self).setup(**kwargs)

        # Get the model
        if kwargs.get("model", None) is not None: self.model = kwargs.pop("model")
        else: self.load_model()

    # -----------------------------------------------------------------

    def load_model(self):

        """
        This function ...
        :return:
        """

        # Debuggging
        log.debug("Loading the model ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def project_definition(self):

        """
        Thisf unction ...
        :return: 
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def project_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def show_parameters_definition(self):

        """
        Thisf unction ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("category", "string", "which parameters to show", all_name, choices=parameter_categories)
        return definition

    # -----------------------------------------------------------------

    def show_parameters_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.show_parameters_definition, **kwargs)

        # All parameters
        if config.category == all_name: self.show_all_parameters()

        # Intrinsic
        elif config.category == intrinsic_name: self.show_intrinsic_parameters()

        # Derived
        elif config.category == derived_name: self.show_derived_parameters()

        # Stellar
        elif config.category == stellar_name: self.show_stellar_parameters()

        # Bulge
        elif config.category == bulge_name: self.show_derived_parameters_bulge()

        # Disk
        elif config.category == disk_name: self.show_derived_parameters_disk()

        # Old
        elif config.category == old_name: self.show_derived_parameters_old()

        # Young
        elif config.category == young_name: self.show_derived_parameters_young()

        # SFR
        elif config.category == sfr_name: self.show_derived_parameters_sfr()

        # Unevolved
        elif config.category == unevolved_name: self.show_derived_parameters_unevolved()

        # Dust
        elif config.category == dust_name: self.show_dust_parameters()

        # Invalid
        else: raise ValueError("Invalid category: '" + config.category + "'")

    # -----------------------------------------------------------------

    def show_all_parameters(self):

        """
        This function ...
        :return:
        """

        # Intrinsic
        self.show_intrinsic_parameters()

        # Derived
        self.show_derived_parameters()

    # -----------------------------------------------------------------

    def show_intrinsic_parameters(self):

        """
        This function ...
        :return:
        """

        # Free
        self.show_free_parameters()

        # Other
        self.show_other_parameters()

    # -----------------------------------------------------------------

    def show_stellar_parameters(self):

        """
        This function ...
        :return:
        """

        # Bulge
        self.show_derived_parameters_bulge()

        # Disk
        self.show_derived_parameters_disk()

        # Old
        self.show_derived_parameters_old()

        # Young
        self.show_derived_parameters_young()

        # SFR
        self.show_derived_parameters_sfr()

        # Unevolved
        self.show_derived_parameters_unevolved()

    # -----------------------------------------------------------------

    def show_dust_parameters(self):

        """
        This function ...
        :return:
        """

        # Dust
        self.show_derived_parameters_dust()

    # -----------------------------------------------------------------

    def show_derived_parameters(self):

        """
        This function ...
        :return:
        """

        # Stars
        self.show_stellar_parameters()

        # Dust
        self.show_dust_parameters()

    # -----------------------------------------------------------------

    @property
    def free_parameter_values(self):
        return self.model.free_parameter_values

    # -----------------------------------------------------------------

    def show_free_parameters(self):

        """
        This function ...
        :return:
        """

        # Show the free parameter values
        print(fmt.cyan + fmt.underlined + "Free parameter values:" + fmt.reset)
        print("")
        for label in self.free_parameter_values: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.free_parameter_values[label]))
        print("")

    # -----------------------------------------------------------------

    @property
    def other_parameter_values(self):
        return self.model.other_parameter_values

    # -----------------------------------------------------------------

    def show_other_parameters(self):

        """
        This function ...
        :return:
        """

        # Show the other parameter values
        print(fmt.cyan + fmt.underlined + "Other parameter values:" + fmt.reset)
        print("")
        for label in self.other_parameter_values: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.other_parameter_values[label]))
        print("")

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_total(self):
        return self.model.derived_parameter_values_total

    # -----------------------------------------------------------------

    def show_derived_parameters_total(self):

        """
        This function ...
        :return:
        """

        # Derived parameter values of total model
        print(fmt.cyan + fmt.underlined + "Derived parameter values of total model:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_total: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_total[label]))
        print("")

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_bulge(self):
        return self.model.derived_parameter_values_bulge

    # -----------------------------------------------------------------

    def show_derived_parameters_bulge(self):

        """
        This function ...
        :return:
        """

        # Derived parameter values of bulge
        print(fmt.cyan + fmt.underlined + "Derived parameter values of old bulge stellar component:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_bulge: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_bulge[label]))
        print("")

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_disk(self):
        return self.model.derived_parameter_values_disk

    # -----------------------------------------------------------------

    def show_derived_parameters_disk(self):

        """
        This function ...
        :return:
        """

        # Derived parameter values of disk
        print(fmt.cyan + fmt.underlined + "Derived parameter values of old disk stellar component:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_disk: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_disk[label]))
        print("")

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_old(self):
        return self.model.derived_parameter_values_old

    # -----------------------------------------------------------------

    def show_derived_parameters_old(self):

        """
        This function ...
        :return:
        """

        # Derived parameter values of old component
        print(fmt.cyan + fmt.underlined + "Derived parameter values of old stellar component:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_old: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_old[label]))
        print("")

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_young(self):
        return self.model.derived_parameter_values_young

    # -----------------------------------------------------------------

    def show_derived_parameters_young(self):

        """
        This function ...
        :return:
        """

        # Derived parameter values of young component
        print(fmt.cyan + fmt.underlined + "Derived parameter values of young stellar component:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_young: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_young[label]))
        print("")

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_sfr(self):
        return self.model.derived_parameter_values_sfr

    # -----------------------------------------------------------------

    def show_derived_parameters_sfr(self):

        """
        This function ...
        :return:
        """

        # Derived parameter values of SF component
        print(fmt.cyan + fmt.underlined + "Derived parameter values of SFR component:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_sfr: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_sfr[label]))
        print("")

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_unevolved(self):
        return self.model.derived_parameter_values_unevolved

    # -----------------------------------------------------------------

    def show_derived_parameters_unevolved(self):

        """
        This function ...
        :return:
        """

        # Derived parameter values of unevolved components
        print(fmt.cyan + fmt.underlined + "Derived parameter values of unevolved stars:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_unevolved: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_unevolved[label]))
        print("")

    # -----------------------------------------------------------------

    @property
    def derived_parameter_values_dust(self):
        return self.model.derived_parameter_values_dust

    # -----------------------------------------------------------------

    def show_derived_parameters_dust(self):

        """
        This function ...
        :return:
        """

        # Derived parameter values of dust component
        print(fmt.cyan + fmt.underlined + "Derived parameter values of dust component:" + fmt.reset)
        print("")
        for label in self.derived_parameter_values_dust: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.derived_parameter_values_dust[label]))
        print("")

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @property
    def history_filename(self):

        """
        This function ...
        :return:
        """

        return "model_examination"

# -----------------------------------------------------------------
