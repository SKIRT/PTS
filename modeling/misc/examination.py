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
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr
from ..build.models.stars import basic_stellar_component_names
from ..build.models.stars import bulge_component_name, old_component_name, young_component_name, ionizing_component_name
from ..build.models.dust import basic_dust_component_names, disk_component_name
from ...core.tools import strings
from ...core.basics.configuration import prompt_settings, parse_arguments
from ..projection.model import create_faceon_projection_from_earth_projection, create_edgeon_projection_from_earth_projection
from ..projection.model import ComponentProjections
from ..basics.projection import GalaxyProjection
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ..basics.models import DeprojectionModel3D
from ...core.tools import sequences
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.tools import plotting

# -----------------------------------------------------------------

earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"
projection_names = [earth_name, faceon_name, edgeon_name]

# -----------------------------------------------------------------

# Standard commands
_help_command_name = "help"
_history_command_name = "history"

# Show
_components_command_name = "components"

# Other
_project_command_name = "project"
_parameters_command_name = "parameters"

# -----------------------------------------------------------------

# Define commands
commands = OrderedDict()

# Standard commands
commands[_help_command_name] = ("show_help", False, "show help", None)
commands[_history_command_name] = ("show_history_command", True, "show history of executed commands", None)

# Show stuff
commands[_components_command_name] = ("show_components", False, "show model components", None)

# Other
commands[_project_command_name] = ("project_command", True, "project the model from one or multiple orientations", "component")
commands[_parameters_command_name] = ("show_parameters_command", True, "show model parameters", None)

# -----------------------------------------------------------------

# Subcommands
subcommands = OrderedDict()

# -----------------------------------------------------------------

all_name = "all"
intrinsic_name = "intrinsic"
free_name = "free"
other_name = "other"
derived_name = "derived"
stellar_name = "stellar"
bulge_name = "bulge"
disk_name = "disk"
old_name = "old"
young_name = "young"
sfr_name = "sfr"
unevolved_name = "unevolved"
dust_name = "dust"
parameter_categories = [all_name, intrinsic_name, free_name, other_name, derived_name, stellar_name, bulge_name, disk_name, old_name, young_name, sfr_name, unevolved_name, dust_name]

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

        if self.config.interactive is None: return not self.has_any and not self.do_commands
        else: return self.config.interactive

    # -----------------------------------------------------------------

    @property
    def has_showing(self):
        return self.config.show and self.config.show_components

    # -----------------------------------------------------------------

    @property
    def has_plotting(self):
        return self.config.plot and False

    # -----------------------------------------------------------------

    @property
    def has_writing(self):
        return self.config.write and False

    # -----------------------------------------------------------------

    @property
    def do_show(self):
        return self.has_showing

    # -----------------------------------------------------------------

    @property
    def do_plot(self):
        return self.has_plotting

    # -----------------------------------------------------------------

    @property
    def do_write(self):
        return self.has_writing

    # -----------------------------------------------------------------

    @property
    def has_any(self):

        """
        This function ...
        :return:
        """

        return self.has_showing or self.has_plotting or self.has_writing

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
        if self.do_show: self.show()

        # Plot
        if self.do_plot: self.plot()

        # Write
        if self.do_write: self.write()

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

    @property
    def definition(self):
        return self.model.definition

    # -----------------------------------------------------------------

    @property
    def stellar_component_names(self):

        """
        This function ...
        :return:
        """

        return self.definition.stellar_component_names

    # -----------------------------------------------------------------

    def is_stellar_component(self, component_name):
        return component_name in self.stellar_component_names

    # -----------------------------------------------------------------

    @property
    def dust_component_names(self):

        """
        This function ...
        :return:
        """

        return self.definition.dust_component_names

    # -----------------------------------------------------------------

    def is_dust_component(self, component_name):
        return component_name in self.dust_component_names

    # -----------------------------------------------------------------

    @lazyproperty
    def component_names(self):

        """
        This function ...
        :return:
        """

        if not sequences.all_different_contents(self.stellar_component_names, self.dust_component_names): raise ValueError("Stellar and dust component names cannot overlap")
        return self.stellar_component_names + self.dust_component_names

    # -----------------------------------------------------------------

    @memoize_method
    def has_deprojection(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        geometry = self.get_component_geometry(component_name)
        return isinstance(geometry, DeprojectionModel3D)

    # -----------------------------------------------------------------

    @memoize_method
    def get_component_scaleheight(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Get geometry
        geometry = self.get_component_geometry(component_name)

        # Deprojection?
        if self.has_deprojection(component_name): return geometry.scale_height

        # Has axial scale defined
        elif hasattr(geometry, "axial_scale"): return geometry.axial_scale

        # Has effective radius defined
        elif hasattr(geometry, "effective_radius"):
            radius = geometry.effective_radius
            if hasattr(geometry, "z_flattening"):
                flattening = geometry.z_flattening
                return radius * flattening
            else: return radius

        # Scale height cannot be determined
        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def old_component_scaleheight(self):
        return self.get_component_scaleheight(old_component_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_component_scaleheight(self):
        return self.get_component_scaleheight(young_component_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_component_scaleheight(self):
        return self.get_component_scaleheight(ionizing_component_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_component_scaleheight(self):
        return self.get_component_scaleheight(disk_component_name)

    # -----------------------------------------------------------------

    @memoize_method
    def get_component(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        if self.is_stellar_component(component_name): return self.definition.load_stellar_component(component_name)
        elif self.is_dust_component(component_name): return self.definition.load_dust_component(component_name)
        else: raise ValueError("Invalid component name: '" + component_name + "'")

    # -----------------------------------------------------------------

    @memoize_method
    def get_component_geometry(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        # Get the component
        component = self.get_component(component_name)

        # Return the geometry
        if "deprojection" in component:
            if "model" in component: raise ValueError("Something went wrong")
            return component.deprojection
        elif "model" in component: return component.model
        else: return None

    # -----------------------------------------------------------------

    @memoize_method
    def get_component_deprojection(self, component_name):
        if not self.has_deprojection(component_name): return None
        return self.get_component_geometry(component_name)

    # -----------------------------------------------------------------

    @memoize_method
    def get_component_map_path(self, component_name):
        if not self.has_deprojection(component_name): return None
        return self.get_component_geometry(component_name).filepath

    # -----------------------------------------------------------------

    @memoize_method
    def get_component_map(self, component_name):
        if not self.has_deprojection(component_name): return None
        return self.get_component_geometry(component_name).map

    # -----------------------------------------------------------------

    @memoize_method
    def get_component_map_wcs(self, component_name):
        if not self.has_deprojection(component_name): return None
        return CoordinateSystem.from_file(self.get_component_map_path(component_name))

    # -----------------------------------------------------------------

    @memoize_method
    def get_component_input(self, component_name):

        """
        This function ...
        :param component_name:
        :return:
        """

        model_map_basename = "map"
        model_map_filename = model_map_basename + ".fits"
        map_filename = "map.fits"

        # Add the map filepath
        paths = dict()
        if self.has_deprojection(component_name): paths[map_filename] = self.get_component_geometry(component_name).filepath
        return paths

    # -----------------------------------------------------------------

    def get_component_command_definition(self, command_definition=None, required=True, choices=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required:
        :param choices:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Add host setting
        if required: definition.add_required("component", "string", "component name", choices=choices)
        else: definition.add_positional_optional("component", "string", "component name", choices=choices)

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_component_command(self, command, command_definition=None, name=None, index=1, required=True, choices=None,
                           required_to_optional=True, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required:
        :param choices:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        parse_command = splitted[index:]

        # Get the definition
        definition = self.get_component_command_definition(command_definition, required=required, choices=choices, required_to_optional=required_to_optional)

        # Get settings interactively
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False, add_config_path=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get the component name
        component_name = config.pop("component")
        if component_name not in self.component_names: raise InvalidCommandError("Invalid component name: '" + component_name + "'", command)

        # Return
        return splitted, component_name, config

    # -----------------------------------------------------------------

    def get_component_name_from_command(self, command, name=None, index=1, required=True, choices=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param index:
        :param required:
        :param choices:
        :param interactive:
        :return:
        """

        # Parse
        splitted, component_name, config = self.parse_component_command(command, name=name, index=index, required=required, choices=choices, interactive=interactive)

        # Return the name
        return component_name

    # -----------------------------------------------------------------

    def get_component_name_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, component_name, config = self.parse_component_command(command, command_definition=command_definition, name=name, interactive=interactive)

        # Return the component name and config
        return component_name, config

    # -----------------------------------------------------------------

    def show_components(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        print("")
        print(fmt.green + fmt.underlined + "STELLAR COMPONENTS" + fmt.reset)

        # Stellar
        self.show_stellar_components()

        # print("")
        print(fmt.green + fmt.underlined + "DUST COMPONENTS" + fmt.reset)

        # Dust
        self.show_dust_components()

    # -----------------------------------------------------------------

    def show_stellar_components(self):

        """
        This function ...
        :return:
        """

        from ..build.models.galaxy import show_component

        print("")
        print("  " + fmt.magenta + "OLD STELLAR BULGE:" + fmt.reset)
        print("")

        # Show
        bulge_path = self.definition.get_stellar_component_path(bulge_component_name)
        show_component(bulge_path, line_prefix="    ")

        print("  " + fmt.magenta + "OLD STELLAR DISK: " + fmt.reset)
        print("")

        # Show
        old_path = self.definition.get_stellar_component_path(old_component_name)
        show_component(old_path, line_prefix="    ")

        print("  " + fmt.magenta + "YOUNG STELLAR DISK:" + fmt.reset)
        print("")

        # Show
        young_path = self.definition.get_stellar_component_path(young_component_name)
        show_component(young_path, line_prefix="    ")

        print("  " + fmt.magenta + "IONIZING STELLAR DISK: " + fmt.reset)
        print("")

        # Show
        ionizing_path = self.definition.get_stellar_component_path(ionizing_component_name)
        show_component(ionizing_path, line_prefix="    ")

        # Loop over the extra stellar components
        for name in self.definition.additional_stellar_names:

            # Show the name
            print("  " + fmt.magenta + name.upper() + ":" + fmt.reset)
            print("")

            path = self.definition.get_stellar_component_path(name)

            # Show the component
            show_component(path, line_prefix="    ")

    # -----------------------------------------------------------------

    def show_dust_components(self):

        """
        This function ...
        :return:
        """

        from ..build.models.galaxy import show_component

        print("  " + fmt.magenta + "DUST DISK:" + fmt.reset)
        print("")

        # Show component
        dust_path = self.definition.get_dust_component_path(disk_component_name)
        show_component(dust_path, line_prefix="    ")

        # Loop over the extra dust components
        for name in self.definition.additional_dust_names:

            # Show the name
            print("  " + fmt.magenta + name.upper() + ": " + fmt.reset)
            print("")

            path = self.definition.get_dust_component_path(name)

            # Show
            show_component(path, line_prefix="    ")

    # -----------------------------------------------------------------

    @lazyproperty
    def project_definition(self):

        """
        Thisf unction ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Projections
        definition.add_positional_optional("orientations", "string_list", "orientation(s) from which to make a projection", [earth_name], choices=projection_names)

        # Paths
        definition.add_optional("path", "directory_path", "path for projections and maps")
        definition.add_optional("output_path", "directory_path", "path for (temporary) SKIRT output")

        # WCS path or component
        definition.add_optional("wcs_path", "file_path", "path of file containing the earth coordinate system")
        definition.add_optional("wcs_component", "string", "name of the component from which to take the coordinate system (must be defined as deprojection model)")

        # Vertical extent of the total model
        definition.add_optional("old_scale_heights", "real", "number of times to take the old stellar scale height as the vertical radius of the model", 2.5)
        definition.add_optional("scale_height", "length_quantity", "scale height for the component")
        definition.add_optional("radial_factor", "real", "factor with which to multiply the radial extent of the projections", 1.5)

        # Flags
        definition.add_flag("open_path", "open the path when the projections are created")
        definition.add_flag("open_output_path", "open the output path when the projections are created")
        definition.add_flag("open_maps", "open the projected maps")
        definition.add_flag("plot", "plot the projected maps")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def project_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the component name and configuration
        component_name, config = self.get_component_name_and_config_from_command(command, self.project_definition, **kwargs)

        ## Get the coordinate system

        # WCS from file
        if config.wcs_path is not None:
            if config.wcs_component is not None: raise ValueError("Cannot specify both 'wcs_path' and 'wcs_component'")
            earth_wcs = CoordinateSystem.from_file(config.wcs_path)

        # WCS from other component
        elif config.wcs_component is not None: earth_wcs = self.get_component_map_wcs(config.wcs_component)

        # WCS not explicitly given
        else:
            if self.has_deprojection(component_name): earth_wcs = self.get_component_map_wcs(component_name)
            else: raise ValueError("Projection cannot be defined for component '" + component_name + "' without a coordinate system as input")

        ## Get the projections

        # Create earth projection
        azimuth = 0.0
        earth_projection = GalaxyProjection.from_wcs(earth_wcs, self.model.center, self.model.distance, self.model.inclination, azimuth, self.model.position_angle)

        # Create other projections if necessary
        projections = OrderedDict()
        if earth_name in config.orientations: projections[earth_name] = earth_projection
        if faceon_name in config.orientations:
            faceon_projection = create_faceon_projection_from_earth_projection(earth_projection, radial_factor=config.radial_factor)
            projections[faceon_name] = faceon_projection
        if edgeon_name in config.orientations:

            if config.scale_height is not None: scale_height = config.scale_height
            else: scale_height = config.old_scale_heights * self.old_component_scaleheight

            edgeon_projection = create_edgeon_projection_from_earth_projection(earth_projection, scale_height, radial_factor=config.radial_factor)
            projections[edgeon_name] = edgeon_projection

        # Projection
        self.project_component(component_name, projections, path=config.path, output_path=config.output_path,
                               earth_wcs=earth_wcs, open_path=config.open_path, open_output_path=config.open_output_path,
                               open_maps=config.open_maps, plot=config.plot)

    # -----------------------------------------------------------------

    def project_component(self, component_name, projections, path=None, output_path=None, earth_wcs=None, open_path=False,
                          open_output_path=False, open_maps=False, plot=False):

        """
        This function ...
        :param component_name:
        :param projections:
        :param path:
        :param output_path:
        :param earth_wcs: wcs to give to the earth map
        :param open_path:
        :param open_output_path:
        :param open_maps:
        :param plot:
        :return:
        """

        # Debugging
        log.debug("Projecting component '" + component_name + "' from orientations " + tostr(projections.keys()) + " ...")

        # Get the component model geometry
        geometry = self.get_component_geometry(component_name)

        # Create a temporary path if necessary
        if output_path is None: output_path = introspection.create_unique_temp_dir("projection__" + component_name)
        elif not fs.is_directory(output_path): raise IOError("Output directory '" + output_path + "' does not exist")

        # Check the projections
        has_earth = earth_name in projections
        has_faceon = faceon_name in projections
        has_edgeon = edgeon_name in projections
        if has_earth: earth_projection = projections[earth_name]
        else: earth_projection = None
        if has_faceon: faceon_projection = projections[faceon_name]
        else: faceon_projection = None
        if has_edgeon: edgeon_projection = projections[edgeon_name]
        else: edgeon_projection = None

        # Set the input filepaths
        input_filepaths = self.get_component_input(component_name)

        #print(earth_projection, faceon_projection, edgeon_projection)
        #print(earth_wcs)

        # Create the projections
        projs = ComponentProjections(component_name, geometry, path=output_path, earth=has_earth, faceon=has_faceon, edgeon=has_edgeon,
                                            projection=earth_projection, projection_faceon=faceon_projection,
                                            projection_edgeon=edgeon_projection, center=self.model.center,
                                            earth_wcs=earth_wcs, distance=self.model.distance, input_filepaths=input_filepaths)

        # Get the maps
        maps = OrderedDict()
        if has_earth: maps[earth_name] = projs.earth
        if has_faceon: maps[faceon_name] = projs.faceon
        if has_edgeon: maps[edgeon_name] = projs.edgeon

        # Write?
        if path is not None:
            earth_map_path = fs.join(path, "earth.fits")
            earth_projection_path = fs.join(path, "earth.proj")
            faceon_map_path = fs.join(path, "faceon.fits")
            faceon_projection_path = fs.join(path, "faceon.proj")
            edgeon_map_path = fs.join(path, "edgeon.fits")
            edgeon_projection_path = fs.join(path, "edgeon.proj")
            if has_earth:
                projections[earth_name].saveto(earth_projection_path)
                maps[earth_name].saveto(earth_map_path)
            if has_faceon:
                projections[faceon_name].saveto(faceon_projection_path)
                maps[faceon_name].saveto(faceon_map_path)
            if has_edgeon:
                projections[edgeon_name].saveto(edgeon_projection_path)
                maps[edgeon_name].saveto(edgeon_map_path)

        # Plot?
        if plot:
            earth_map_plot_path = fs.join(path, "earth.pdf") if path is not None else None
            faceon_map_plot_path = fs.join(path, "faceon.pdf") if path is not None else None
            edgeon_map_plot_path = fs.join(path, "edgeon.pdf") if path is not None else None
            if has_earth: plotting.plot_map(maps[earth_name], scale="auto", path=earth_map_plot_path)
            if has_faceon: plotting.plot_map(maps[faceon_name], scale="auto", path=faceon_map_plot_path)
            if has_edgeon: plotting.plot_map(maps[edgeon_name], scale="auto", path=edgeon_map_plot_path)

        # Open?
        if open_path:
            if path is None: log.warning("Cannot open the path because it was not defined")
            else: fs.open_directory(path)
        if open_output_path: fs.open_directory(output_path)
        if open_maps:
            if has_earth:
                if path is not None: fs.open_file(fs.join(path, "earth.fits")) # has been saved with WCS
                else: fs.open_file(projs.earth_map_path)
            if has_faceon:
                if path is not None: fs.open_file(fs.join(path, "faceon.fits"))
                else: fs.open_file(projs.faceon_map_path)
            if has_edgeon:
                if path is not None: fs.open_file(fs.join(path, "edgeon.fits"))
                else: fs.open_file(projs.edgeon_map_path)

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

        # Empty line
        print("")

        # All parameters
        if config.category == all_name: self.show_all_parameters()

        # Intrinsic
        elif config.category == intrinsic_name: self.show_intrinsic_parameters()

        # Free
        elif config.category == free_name: self.show_free_parameters()

        # Other
        elif config.category == other_name: self.show_other_parameters()

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

        # Inform the user
        log.info("Showing ...")

        # Components
        if self.config.show_components: self.show_components()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    @property
    def history_filename(self):

        """
        This function ...
        :return:
        """

        return "model_examination"

# -----------------------------------------------------------------
