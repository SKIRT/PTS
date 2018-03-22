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
from ..basics.configuration import ConfigurationDefinition
from .smile import SKIRTSmileSchema, skirt_quantities_to_pts_quantities
from ..simulation.skifile import SkiFile
from ..tools.utils import lazyproperty
from ..basics.configuration import prompt_choice
from ..tools import formatting as fmt
from ..tools.stringify import stringify, stringify_not_list, stringify_string_fancy
from ..units.stringify import stringify_unit

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

        # Create ski file
        self.create_ski()

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
        if kwargs.get("skifile", None) is not None: self.load_ski(kwargs.pop("skifile"))
        elif kwargs.get("ski", None) is not None: self.load_ski(kwargs.pop("ski"))
        elif self.config.from_ski is not None: self.load_ski(SkiFile(self.config.from_ski))

    # -----------------------------------------------------------------

    def load_ski(self, ski):

        """
        This function ...
        :param ski:
        :return:
        """

        # Debugging
        log.debug("Loading components from ski file ...")



    # -----------------------------------------------------------------

    def make_stellar_component(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Making stellar component ...")

        # Get the configuration parameters interactively
        parameters, children = self.smile.prompt_parameters_for_type("PanStellarComp", merge=True)

        # Show
        print(parameters)

    # -----------------------------------------------------------------

    def make_dust_component(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Making dust component ...")

        # Get the configuration parameters interactively
        parameters = self.smile.prompt_parameters_for_type("DustComp", merge=True)

        # Show
        print(parameters)

    # -----------------------------------------------------------------

    def make_item(self):

        """
        This function ...
        :return:
        """

        # Create the SKIRT smile schema
        smile = SKIRTSmileSchema()

        # Get the configuration parameters interactively
        parameters, children = smile.prompt_parameters_for_type(name)

        # Show the parameters
        show_parameters(parameters, children)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_items_definition(self):

        """
        Thisf unction ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("match", "string", "only show types with names that contain this string")
        definition.add_flag("definitions", "format the properties for each item as a configuration definition")
        return definition

    # -----------------------------------------------------------------

    def show_items_command(self, command, **kwargs):

        """
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    def show_items(self):

        """
        This function ...
        :return:
        """

        print("")

        # Loop over the concrete types
        for name in self.smile.concrete_types:

            if config.match is not None and config.match not in name: continue

            print(fmt.green + fmt.underlined + name + fmt.reset)
            print("")

            description = self.smile.concrete_types[name]
            print(fmt.darkgray + stringify_string_fancy(description, lines_prefix=" ")[1] + fmt.reset)
            print("")

            if config.definitions:

                # Create definition
                definition, simulation_items = self.smile.definition_for_type(name)

                # Create output string
                output = StringIO.StringIO()

                # Write definition to string buffer
                write_definition(definition, output, indent=" ")

                # Print the definition
                contents = output.getvalue()
                line_description = None
                for line in contents.split("\n"):
                    if line.strip().startswith("#"):
                        line_description = line.split("# ")[1]
                        continue
                    else:
                        if line_description is not None: print(line + "   " + fmt.darkgray + line_description + fmt.reset)
                        else: print(line)
                        line_description = None

                # Close string buffer
                output.close()

            else:

                # Get properties
                properties = self.smile.properties_for_type(name)

                # Show the properties
                for prop_name in properties:

                    description = properties[prop_name].description
                    ptype = properties[prop_name].ptype
                    min_value = properties[prop_name].min
                    max_value = properties[prop_name].max
                    default = properties[prop_name].default
                    choices = properties[prop_name].choices
                    item = properties[prop_name].item

                    string = fmt.blue + prop_name + fmt.reset

                    string += " [type: " + ptype + "]"
                    if min_value is not None: string += " [min: " + stringify_not_list(min_value)[1] + "]"
                    if max_value is not None: string += " [max: " + stringify_not_list(max_value)[1] + "]"
                    if default is not None: string += " [default: " + stringify_not_list(default)[1] + "]"
                    if choices is not None: string += " [choices: " + stringify(choices.keys())[1] + "]"

                    string += " " + fmt.darkgray + stringify_string_fancy(description, lines_prefix="    ")[1] + fmt.reset

                    if item: string += fmt.red + " [SIMULATION ITEM] " + fmt.reset

                    print(" -  " + string)
                    print("")

            # print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def show_quantites_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("match", "string",
                                           "only show quantities with names that contain this string")
        return definition

    # -----------------------------------------------------------------

    def show_quantities(self):

        """
        This function ...
        :return:
        """

        print("")
        for quantity in self.smile.quantities:
            units = self.smile.units_for_quantity(quantity)
            print(fmt.green + fmt.underlined + quantity + fmt.reset + ": " + stringify(units)[1].replace(" ", "").replace(",", ", ") + " [PTS: " + skirt_quantities_to_pts_quantities[quantity] + "]")
            print("")

        # print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def show_units_definition(self):

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("match", "string",
                                           "only show unit systems with names that contain this string")
        return definition

    # -----------------------------------------------------------------

    def show_units(self):

        """
        This function ...
        :return:
        """

        print("")

        # Loop over
        for unit_system in self.smile.unit_systems:

            # Get default units
            units = self.smile.units_for_unit_system(unit_system)

            print(fmt.green + fmt.underlined + unit_system + fmt.reset)
            print("")

            for quantity_name in units:
                print(" - " + fmt.blue + quantity_name + fmt.reset + ": " + stringify_unit(units[quantity_name])[1])

            print("")

    # -----------------------------------------------------------------

    # def create_ski_template(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Debugging
    #     log.debug("Creating empty ski template ...")
    #
    #     # Oligo or pan?
    #     if self.config.type is None: simulation_type = prompt_choice("type", "simulation type", "")
    #     else: simulation_type = self.config.type
    #
    #     # Create template ski file
    #     if simulation_type == oligo_type: self.ski = self.smile.create_oligochromatic_template()
    #     elif simulation_type == pan_type: self.ski =self.smile.create_panchromatic_template()
    #     else: raise ValueError("Invalid value for 'simulation_type'")

    # -----------------------------------------------------------------

    @property
    def simulation_type(self):

        """
        This function ...
        :return:
        """

        return self.ski.simulation_type

    # -----------------------------------------------------------------

    def create_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating ski file ...")

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
