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
from ..basics.configuration import prompt_choice, prompt_variable
from ..tools import formatting as fmt
from ..tools.stringify import stringify, stringify_not_list, stringify_string_fancy
from ..units.stringify import stringify_unit
from ..tools import strings

# -----------------------------------------------------------------

oligo_type = "oligo"
pan_type = "pan"
oligo_or_pan = [oligo_type, pan_type]

# -----------------------------------------------------------------

_help_command_name = "help"
_history_command_name = "history"
_show_command_name = "show"
_add_command_name = "add"
_remove_command_name = "remove"
_set_command_name = "set"
_plot_command_name = "plot"

# -----------------------------------------------------------------

_stellar_command_name = "stellar"
_dust_command_name = "dust"
_instrument_command_name = "instrument"
_item_command_name = "item"
_wavelengths_command_name = "wavelengths"
_grid_command_name = "grid"
_npackages_command_name = "npackages"
_random_command_name = "random"
_units_command_name = "units"
#_emissivity_command_name = ""

# -----------------------------------------------------------------

# Define show commands
show_commands = OrderedDict()
show_commands.description = "show"
show_commands[_stellar_command_name] = ("show_stellar_command", True, "show a stellar component", "stellar")
show_commands[_dust_command_name] = ("show_dust_command", True, "show a dust component", "dust")
show_commands[_instrument_command_name] = ("show_instrument_command", True, "show an instrument", "instrument")

# -----------------------------------------------------------------

# Define add commands
add_commands = OrderedDict()
add_commands.description = "add"
add_commands[_stellar_command_name] = ("add_stellar_command", True, "add a stellar component", None)
add_commands[_dust_command_name] = ("add_dust_command", True, "add a dust component", None)
add_commands[_instrument_command_name] = ("add_instrument_command", True, "add an instrument", None)
add_commands[_item_command_name] = ("add_item_command", True, "add an item", None)

# -----------------------------------------------------------------

# Define remove commands
remove_commands = OrderedDict()
remove_commands.description = "remove"
remove_commands[_stellar_command_name] = ("remove_stellar_command", True, "remove a stellar component", "stellar")
remove_commands[_dust_command_name] = ("remove_dust_command", True, "remove a dust component", "dust")
remove_commands[_instrument_command_name] = ("remove_instrument_command", True, "remove an instrument", "instrument")

# -----------------------------------------------------------------

# Define set commands
set_commands = OrderedDict()
set_commands.description = "set"
set_commands[_wavelengths_command_name] = ("set_wavelengths_command", True, "set the wavelength grid", None)
set_commands[_grid_command_name] = ("set_grid_command", True, "set the dust grid", None)
set_commands[_npackages_command_name] = ("set_npackages_command", True, "set the number of photon packages", None)
set_commands[_random_command_name] = ("set_random_command", True, "set the random seed", None)
set_commands[_units_command_name] = ("set_units_command", True, "set the units system", None)

# -----------------------------------------------------------------

# Define plot commands
plot_commands = OrderedDict()
plot_commands.description = "plot"
plot_commands[_wavelengths_command_name] = ("plot_wavelengths_command", True, "plot the wavelength grid", None)

# -----------------------------------------------------------------

# Set commands
commands = OrderedDict()

# Standard commands
commands[_help_command_name] = ("show_help", False, "show help", None)
commands[_history_command_name] = ("show_history_command", True, "show history of executed commands", None)

# Showing stuff
commands[_show_command_name] = show_commands # (None, False, "show", None)

# Add & Remove
commands[_add_command_name] = add_commands # (None, False, "add", None)
commands[_remove_command_name] = remove_commands # (None, False, "remove", None)

# Set
commands[_set_command_name] = set_commands # (None, False, "set", None)

# Plot
commands[_plot_command_name] = plot_commands #(None, False, "plot", None)

# -----------------------------------------------------------------

class ModelComposer(InteractiveConfigurable):

    """
    This function ...
    """

    _commands = commands

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

        # The stellar components
        self.stellar = OrderedDict()

        # The dust components
        self.dust = OrderedDict()

        # Other items
        self.items = OrderedDict()

        # The instruments
        self.instruments = OrderedDict()

        # The dust grid and wavelength grid
        self.dust_grid = None
        #self.wavelength_grid = None

        # Other properties
        self.npackages = None
        self.selfabsorption = None
        self.transient_heating = None
        self.dustlib = None # type of dust lib

        # The ski file
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

    def _run(self, **kwargs):

        """
        Thisf unction ...
        :param kwargs:
        :return:
        """

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

    @property
    def stellar_component_names(self):

        """
        This function ...
        :return:
        """

        return self.stellar.keys()

    # -----------------------------------------------------------------

    @property
    def dust_component_names(self):

        """
        This function ...
        :return:
        """

        return self.dust.keys()

    # -----------------------------------------------------------------

    @property
    def instrument_names(self):

        """
        This function ...
        :return:
        """

        return self.instruments.keys()

    # -----------------------------------------------------------------

    def get_stellar_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("component_name", "string", "stellar component name", choices=self.stellar_component_names)

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_stellar_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True,
                              interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        if command_definition is not None: parse_command = splitted[index:]
        else: parse_command = splitted[index:index + 1]  # only component name

        # Get the definition
        definition = self.get_stellar_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get the configuration
        config = self.get_config_from_definition(name, definition, parse_command, interactive=interactive)

        # Get component name
        component_name = config.pop("component_name")

        # Return
        return splitted, component_name, config

    # -----------------------------------------------------------------

    def get_stellar_component_name_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :return:
        """

        # Parse the command
        splitted, component_name, config = self.parse_stellar_command(command, name=name, interactive=interactive)

        # Return the component name
        return component_name

    # -----------------------------------------------------------------

    def get_stellar_component_name_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, component_name, config = self.parse_stellar_command(command, command_definition, name=name, interactive=interactive)

        # Return the component name
        return component_name, config

    # -----------------------------------------------------------------

    def get_dust_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("component_name", "string", "dust component name", choices=self.stellar_component_names)

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_dust_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True,
                              interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        if command_definition is not None: parse_command = splitted[index:]
        else: parse_command = splitted[index:index + 1]  # only component name

        # Get the definition
        definition = self.get_dust_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get the configuration
        config = self.get_config_from_definition(name, definition, parse_command, interactive=interactive)

        # Get component name
        component_name = config.pop("component_name")

        # Return
        return splitted, component_name, config

    # -----------------------------------------------------------------

    def get_dust_component_name_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, component_name, config = self.parse_dust_command(command, name=name, interactive=interactive)

        # Return the component name
        return component_name

    # -----------------------------------------------------------------

    def get_dust_component_name_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        Thisf unction ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, component_name, config = self.parse_dust_command(command, command_definition, name=name, interactive=interactive)

        # Return the component name
        return component_name, config

    # -----------------------------------------------------------------

    def get_instrument_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("instrument_name", "string", "instrument name", choices=self.instrument_names)

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_instrument_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True,
                              interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        if command_definition is not None: parse_command = splitted[index:]
        else: parse_command = splitted[index:index + 1]  # only instrument name

        # Get the definition
        definition = self.get_instrument_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get the configuration
        config = self.get_config_from_definition(name, definition, parse_command, interactive=interactive)

        # Get instrument name
        instrument_name = config.pop("instrument_name")

        # Return
        return splitted, instrument_name, config

    # -----------------------------------------------------------------

    def get_instrument_name_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, instrument_name, config = self.parse_instrument_command(command, name=name, interactive=interactive)

        # Return the instrument name
        return instrument_name

    # -----------------------------------------------------------------

    def get_instrument_name_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, instrument_name, config = self.parse_instrument_command(command, command_definition, name=name, interactive=interactive)

        # Return the instrument name
        return instrument_name, config

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

        # Load settings
        self.load_settings(ski)

        # Load stellar components
        self.load_stellar_components(ski)

        # Load dust components
        self.load_dust_components(ski)

        # Load instruments
        self.load_instruments(ski)

        # Load wavelength grid
        self.load_wavelength_grid(ski)

        # Load dust grid
        self.load_dust_grid(ski)

    # -----------------------------------------------------------------

    def load_settings(self, ski):

        """
        This function ...
        :return:
        """

        # Set type
        self.config.type = ski.simulation_type

    # -----------------------------------------------------------------

    def load_stellar_components(self, ski):

        """
        This function ...
        :return:
        """

        # Loop over the stellar components
        for component_id in ski.get_stellar_component_ids():

            #properties = ski.get_stellar_component_properties(component_id)

            # Get normalization
            lum, wavelength_or_filter = ski.get_stellar_component_luminosity(component_id)

            # Print
            #if isinstance(wavelength_or_filter, Filter):
            #    print(" - filter: " + str(wavelength_or_filter))
            #    print(" - (neutral) spectral luminosity: " + represent_quantity(lum, scientific=True))
            #elif isinstance(wavelength_or_filter, Quantity):
            #    print(" - wavelength: " + represent_quantity(wavelength_or_filter, scientific=True))
            #    print(" - spectral luminosity: " + represent_quantity(lum, scientific=True))
            #elif wavelength_or_filter is None: print(" - bolometric luminosity: " + represent_quantity(lum, scientific=True))
            #else: raise ValueError("Unrecognized filter or wavelength: " + str(wavelength_or_filter))

            # print(" - geometry: " + ski.get_stellar_component_geometry(name).tag)

            #print(" - geometry: " + " > ".join(ski.get_stellar_component_geometry_hierarchy_names(name)))
            #print(" - sed: " + ski.get_stellar_component_sed(name).tag)



    # -----------------------------------------------------------------

    def load_dust_components(self, ski):

        """
        This function ...
        :return:
        """

        # Loop over the dust components
        for component_id in ski.get_dust_component_ids():

            mass = ski.get_dust_component_mass(component_id)

            #print(" - mass: " + represent_quantity(mass, scientific=True))
            #print(" - geometry: " + " > ".join(ski.get_dust_component_geometry_hierarchy_names(name)))
            #print(" - mix: " + ski.get_dust_component_mix(name).tag)
            #if ski.transient_dust_emissivity: print(" - emissivity: transient")
            #if ski.grey_body_dust_emissivity: print(" - emissivity: grey body")



    # -----------------------------------------------------------------

    def load_instruments(self, ski):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def load_wavelength_grid(self, ski):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def load_dust_grid(self, ski):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    @property
    def oligo(self):

        """
        This function ...
        :return:
        """

        return self.config.type == oligo_type

    # -----------------------------------------------------------------

    @property
    def pan(self):

        """
        This function ...
        :return:
        """

        return self.config.type == pan_type

    # -----------------------------------------------------------------

    def show_stellar_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get stellar component name
        component_name = self.get_stellar_component_name_from_command(command, **kwargs)

        # Show
        self.show_stellar_component(component_name)

    # -----------------------------------------------------------------

    def show_stellar_component(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Debugging
        log.debug("Showing stellar component '" + component_name + "' ...")

        # Show
        #show_stellar_component(ski, id)

    # -----------------------------------------------------------------

    def show_dust_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get dust component name
        component_name = self.get_dust_component_name_from_command(command, **kwargs)

        # Show
        self.show_dust_component(component_name)

    # -----------------------------------------------------------------

    def show_dust_component(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Debugging
        log.debug("Showing dust component '" + component_name + "' ...")

        # Show


    # -----------------------------------------------------------------

    def show_instrument_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get instrument name
        instrument_name = self.get_instrument_name_from_command(command, **kwargs)

        # Show
        self.show_instrument(instrument_name)

    # -----------------------------------------------------------------

    def show_instrument(self, instrument_name):

        """
        This function ...
        :param instrument_name:
        :return:
        """

        # Debugging
        log.debug("Showing instrument '" + instrument_name + "' ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def add_stellar_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Component name
        definition.add_required("name", "string", "component name")

        # Options
        #definition.add_required("geometry", "string", "SKIRT base geometry for the component", choices=self.smile.concrete_geometries)
        #definition.add_positional_optional("sed", "string", "SED template for the component", choices=self.smile.concrete_stellar_seds)
        #definition.add_positional_optional("normalization", "string", "normalization for the component", choices=self.smile.concrete_stellar_normalizations)
        #definition.add_optional("luminosities")

        definition.add_flag("by_parameters", "prompt all parameters", False)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def add_stellar_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Parse
        config = self.get_config_from_command(command, self.add_stellar_definition, **kwargs)

        # Prompt parameters
        if config.by_parameters: self.add_stellar_parameters(config.name)

        # Add the component
        else:
            #self.add_stellar_component(config.name, config.geometry, sed_type=config.sed, normalization_type=config.normalization)
            self.add_stellar_component(config.name)

    # -----------------------------------------------------------------

    def add_stellar_parameters(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Debugging
        log.debug("Adding stellar component '" + component_name + "' ...")

        # Get the configuration parameters interactively
        if self.pan: parameters = self.prompt_parameters("PanStellarComp")
        else: parameters = self.prompt_parameters("OligoStellarComp")

        # Add
        self.stellar[component_name] = parameters

    # -----------------------------------------------------------------

    #def add_stellar_component(self, component_name, geometry_type, sed_type=None, normalization_type=None):
    def add_stellar_component(self, component_name):

        """
        This function ...
        :param component_name:
        """

        # Get geometry
        geometry_type, geometry_parameters = self.prompt_geometry()

        # Panchromatic
        if self.pan:

            sed_type, sed_parameters = self.prompt_sed()

            normalization_type, normalization_parameters = self.prompt_normalization()

        # Oligo
        else:

            for wavelength in self.wavelengths:

                luminosity = prompt_variable("luminosity", "photometric_quantity", "luminosity for ...", required=True)


        # Prompt stellar parameters
        #properties = self.prompt_stellar_properties(component_name, geometry_type, sed_type=sed_type, normalization_type=normalization_type)

        # Add the component
        self.stellar[component_name] = properties

    # -----------------------------------------------------------------

    def prompt_geometry(self):

        """
        This function ...
        :return:
        """

        # Prompt geometry
        geometry_type = prompt_choice("geometry", "geometry for the component", self.smile.concrete_geometries)

        # Set geometry properties
        parameters = self.smile.prompt_parameters_for_type(geometry_type, merge=True)

        # Create geometry
        #geometry = Geometry(geometry_type, parameters)

        # Return
        #return geometry

        return geometry_type, parameters

    # -----------------------------------------------------------------

    def prompt_sed(self):

        """
        This function ...
        :return:
        """

        # Prompt
        sed_type = prompt_choice("sed", "SED template for the component", self.smile.concrete_stellar_seds)

        # Get parameters
        #sed_parameters = self.prompt_sed_properties(component_name, sed_type)
        parameters = self.smile.prompt_parameters_for_type(sed_type, merge=True)

        # Create SED
        #SEDTemplate()

        return sed_type, parameters

    # -----------------------------------------------------------------

    def prompt_normalization(self):

        """
        Thisf unction ...
        :return:
        """

        # Prompt
        normalization_type = prompt_choice("normalization", "normalization for the component", self.smile.concrete_stellar_normalizations)

        # Get parameters for this simulation item
        parameters = self.smile.prompt_parameters_for_type(normalization_type, merge=True)

        return normalization_type, parameters

    # -----------------------------------------------------------------

    # def prompt_stellar_properties(self, component_name, geometry_type, sed_type=None, normalization_type=None):
    #
    #     """
    #     This function ...
    #     :param component_name:
    #     :param geometry_type:
    #     :param sed_type:
    #     :param normalization_type:
    #     :return:
    #     """
    #
    #     # Set geometry properties
    #     geometry_parameters = self.prompt_geometry_properties(component_name, geometry_type)
    #
    #     print(geometry_parameters)
    #
    #     # Set SED properties
    #     sed_parameters = self.prompt_sed_properties(component_name, sed_type)
    #
    #     print(sed_parameters)
    #
    #     # Set properties
    #     normalization_parameters = self.prompt_normalization_properties(component_name, normalization_type)
    #
    #     print(normalization_parameters)
    #
    #     # Set the properties
    #     properties = {"sed": sed_parameters, "normalization": normalization_parameters, "geometry": geometry_parameters}
    #
    #     return properties

    # -----------------------------------------------------------------

    def prompt_geometry_properties(self, name, geometry_type):

        """
        This function ...
        :param name:
        :param geometry_type:
        :return:
        """

        # Inform the user
        log.info("Configuring the geometry of stellar component '" + name + "' ...")

        # Get parameters for this simulation item
        parameters = self.smile.prompt_parameters_for_type(geometry_type, merge=True)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def prompt_sed_properties(self, name, sed_type):

        """
        This function ...
        :param name:
        :param sed_type:
        :return:
        """

        # Inform the user
        log.info("Configuring the SED template of stellar component '" + name + "' ...")

        # Get parameters for this simulation item
        parameters = self.smile.prompt_parameters_for_type(sed_type, merge=True)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def prompt_normalization_properties(self, name, normalization_type):

        """
        This function ...
        :param name:
        :param normalization_type:
        :return:
        """

        # Inform the user
        log.info("Configuring the normalization of stellar component '" + name + "' ...")

        # Get parameters for this simulation item
        parameters = self.smile.prompt_parameters_for_type(normalization_type, merge=True)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    @lazyproperty
    def add_dust_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Component name
        definition.add_required("name", "string", "component name")

        #definition.add_optional("title", "string", "short description for this component")
        definition.add_optional("geometry", "string", "SKIRT base geometry for the component", self.smile.concrete_geometries)
        definition.add_optional("normalization", "string", "normalization for the component", self.smile.concrete_dust_normalizations)
        definition.add_optional("mix", "string", "dust mix for the component", self.smile.concrete_dust_mixes)

        definition.add_flag("by_parameters", "prompt all parameters", False)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def add_dust_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get
        config = self.get_config_from_command(command, self.add_dust_definition, **kwargs)

        # Prompt parameters
        if config.by_parameters: self.add_dust_parameters(config.name)

        # Add component
        else: self.add_dust_component(config.name)

    # -----------------------------------------------------------------

    def add_dust_parameters(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Debugging
        log.debug("Adding dust component '" + component_name + "' ...")

        # Get the configuration parameters interactively
        parameters = self.prompt_parameters("DustComp")

        # Add
        self.dust[component_name] = parameters

    # -----------------------------------------------------------------

    def add_dust_component(self, component_name, geometry_type, mix_type=None, normalization_type=None):

        """
        This function ...
        :param component_name:
        :param geometry_type:
        :param mix_type:
        :param normalization_type:
        :return:
        """

        # Set geometry properties
        geometry_parameters = self.set_geometry_properties(name)

        # Set mix properties
        mix_parameters = self.set_mix_properties(name)

        # Set properties
        normalization_parameters = self.set_normalization_properties(name)

        # Create the properties
        properties = {"normalization": normalization_parameters, "geometry": geometry_parameters, "mix": mix_parameters}

        # Add
        self.dust[component_name] = properties

    # -----------------------------------------------------------------

    @lazyproperty
    def add_item_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("name", "string", "name of the item")
        return definition

    # -----------------------------------------------------------------

    def add_item_command(self, command, **kwargs):

        """
        Thisf unction ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.add_item_definition, **kwargs)

        # Add item
        self.add_item(config.name)

    # -----------------------------------------------------------------

    def add_item(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Prompt
        parameters = self.prompt_parameters(name)

        # Show the parameters
        # show_parameters(parameters, children)

        # Add
        self.items[name] = parameters

    # -----------------------------------------------------------------

    def prompt_parameters(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the configuration parameters interactively
        #parameters, children = self.smile.prompt_parameters_for_type(name)
        parameters = self.smile.prompt_parameters_for_type(name, merge=True)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def add_instrument_command(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def add_instrument(self, instrument_name):

        """
        This function ...
        :param instrument_name:
        :return:
        """

    # -----------------------------------------------------------------

    def remove_stellar_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get stellar component name
        component_name = self.get_stellar_component_name_from_command(command, **kwargs)

    # -----------------------------------------------------------------

    def remove_stellar_component(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Debugging
        log.debug("Removing stellar component '" + component_name + "' ...")

    # -----------------------------------------------------------------

    def remove_dust_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get dust component name
        component_name = self.get_dust_component_name_from_command(command, **kwargs)

        # Remove
        self.remove_dust_component(component_name)

    # -----------------------------------------------------------------

    def remove_dust_component(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Debugging
        log.debug("Removing dust component '" + component_name + "' ...")

    # -----------------------------------------------------------------

    def remove_instrument_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get instrument name
        instrument_name = self.get_instrument_name_from_command(command, **kwargs)

        # Remove
        self.remove_instrument(instrument_name)

    # -----------------------------------------------------------------

    def remove_instrument(self, instrument_name):

        """
        This function ...
        :param instrument_name:
        :return:
        """

        # Debugging
        log.debug("Removing instrument '" + instrument_name + "' ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def set_wavelengths_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def set_wavelengths_command(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def set_wavelengths(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def set_grid_command(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def set_npackages_command(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def set_random_command(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def set_units_command(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_wavelengths_command(self):

        """
        This function ...
        :return:
        """

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

        # Create template
        self.create_template()

        # Add the instruments
        self.add_instruments()

        # Set the number of photon packages
        self.ski.setpackages(self.npackages)

        # Set the name of the wavelength grid file
        self.set_wavelength_grid()

        # Set the dust emissivityex
        self.set_dust_emissivity()

        # Set the dust grid
        self.set_dust_grid()

        # Set all-cells dust library
        #self.ski.set_allcells_dust_lib()

        # Set the dust selfabsorption
        self.set_selfabsorption()

    # -----------------------------------------------------------------

    def create_template(self):

        """
        Thisf unction ...
        :return:
        """

        # Oligo or pan?
        #if self.config.type is None: simulation_type = prompt_choice("type", "simulation type", "")
        #else: simulation_type = self.config.type

        # Create template ski file
        if simulation_type == oligo_type: self.ski = self.smile.create_oligochromatic_template()
        elif simulation_type == pan_type: self.ski =self.smile.create_panchromatic_template()
        else: raise ValueError("Invalid value for 'simulation_type'")

    # -----------------------------------------------------------------

    def add_instruments(self):

        """
        This function ...
        :return:
        """

        # Add the instruments
        for name in self.instruments: self.ski.add_instrument(name, self.instruments[name])

    # -----------------------------------------------------------------

    def set_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # For file:
        #self.ski.set_file_wavelength_grid("wavelengths.txt")

    # -----------------------------------------------------------------

    def set_dust_emissivity(self):

        """
        This function ...
        :return:
        """

        # Enable or disable
        if self.config.transient_heating: self.ski.set_transient_dust_emissivity()
        else: self.ski.set_grey_body_dust_emissivity()

    # -----------------------------------------------------------------

    def set_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust grid ...")

        # Check whether we can use the file tree dust grid
        if self.smile.supports_file_tree_grids and self.representation.has_dust_grid_tree:

            # Create file tree dust grid
            dust_grid = self.representation.create_file_tree_dust_grid(write=False)

            # Make sure it is only the file name, not a complete path
            dust_grid.filename = fs.name(dust_grid.filename)

        # Just take the real dust grid object
        else: dust_grid = self.representation.dust_grid

        # Set the lowest-resolution dust grid
        self.ski.set_dust_grid(dust_grid)

    # -----------------------------------------------------------------

    def set_selfabsorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust self-absorption ...")

        # Dust self-absorption
        if self.config.selfabsorption: self.ski.enable_selfabsorption()
        else: self.ski.disable_selfabsorption()

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    def show_summary(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        print("")
        print(fmt.yellow + "GENERAL:" + fmt.reset)
        print("")

        if self.ski.oligochromatic():
            print(" - oligochromatic simulation")
        elif self.ski.panchromatic():
            print(" - panchromatic simulation")

        print("")
        print(fmt.yellow + "WAVELENGTH GRID:" + fmt.reset)
        print("")

        print(" - grid type: " + self.ski.get_wavelength_grid().tag)
        if self.ski.oligochromatic():
            print(" - wavelengths: " + ",".join(self.ski.wavelengths()))
        else:
            print(" - number of wavelengths: " + str(self.ski.nwavelengths()))
            print(" - minimum wavelength: " + represent_quantity(self.ski.min_wavelength))
            print(" - maximum wavelength: " + represent_quantity(self.ski.max_wavelength))

        print("")
        print(fmt.yellow + "INSTRUMENTS:" + fmt.reset)
        print("")

        # Loop over the instruments
        for name in self.ski.get_instrument_names():

            # Get the instrument
            instrument = self.ski.get_instrument_object(name)

            print(fmt.red + fmt.underlined + name + fmt.reset)
            print("")

            properties = instrument.get_properties()
            for label in properties:
                print(" - " + label + ": " + stringify.stringify(properties[label])[1])
            print("")

        print(fmt.yellow + "STELLAR COMPONENTS:" + fmt.reset)
        print("")

        # Loop over the stellar components
        for name in self.ski.get_stellar_component_ids(): show_stellar_component(self.ski, name)

        print(fmt.yellow + "DUST COMPONENTS:" + fmt.reset)
        print("")

        # Loop over the dust components
        for name in self.ski.get_dust_component_ids(): show_dust_component(self.ski, name)

        print(fmt.yellow + "DUST GRID:" + fmt.reset)
        print("")

        print(" - grid type: " + self.ski.get_dust_grid().tag)

        grid = self.ski.get_dust_grid_object()
        properties = grid.get_properties()
        for label in properties:
            print(" - " + label + ": " + stringify.stringify(properties[label], scientific=True)[1])
        print("")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the stellar components
        self.write_stellar()

        # Write the dust components
        self.write_dust()

        # Write the instruments
        self.write_instruments()

        # Write the wavelength grid
        self.write_wavelength_grid()

        # Write the dust grid
        self.write_dust_grid()

        # Skifile
        self.write_ski()

    # -----------------------------------------------------------------

    def write_stellar(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_dust(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_instruments(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_wavelength_grid(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_dust_grid(self):

        """
        This function ...
        :return:
        """

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
