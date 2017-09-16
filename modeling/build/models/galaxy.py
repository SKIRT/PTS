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
from ....core.basics.configuration import prompt_yn, prompt_automatic, prompt_variable
from .stars import bulge_component_name, old_component_name, young_component_name, ionizing_component_name
from .dust import disk_component_name
from .general import write_component
from ....core.tools import formatting as fmt
from ....core.tools import filesystem as fs
from ....core.tools.stringify import tostr

# -----------------------------------------------------------------

bulge_default_parameters = ["default_old_bulge_template", "default_old_bulge_metallicity", "default_old_bulge_age"]
old_default_parameters = ["default_old_disk_template", "default_old_disk_metallicity", "default_old_disk_age"]
young_default_parameters = ["default_young_template", "default_young_metallicity", "default_young_age", "young_scaleheight_ratio"]
ionizing_default_parameters = ["default_ionizing_metallicity", "ionizing_scaleheight_ratio", "default_sfr", "default_ionizing_compactness", "default_ionizing_pressure", "default_covering_factor", "fuv_ionizing_contribution"]
dust_default_parameters = ["default_hydrocarbon_pops", "default_enstatite_pops", "default_forsterite_pops", "default_dust_mass", "dust_scaleheight_ratio"]

# -----------------------------------------------------------------

general_stellar_default_parameters = ["scalelength_to_scaleheight", "use_defaults"]

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
        self.extra_stellar_paths = []
        self.extra_dust_paths = []

        # Map paths
        self.old_stars_map_path = None
        self.young_stars_map_path = None
        self.ionizing_stars_map_path = None
        self.dust_map_path = None
        self.extra_stellar_map_paths = []
        self.extra_dust_map_paths = []

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Adjust stellar components from previous model
        if self.from_previous: self.adjust_stars()

        # 3. Adjust dust components from previous model
        if self.from_previous: self.adjust_dust()

        # 4. Build stars
        self.build_stars()

        # 5. Build dust component
        self.build_dust()

        # 6. Write
        self.write()

        # 7. Show
        if self.config.show: self.show()

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
        for name in self.previous.stellar_component_names:

            # Include?
            if not prompt_yn("include_" + name, "include the '" + name + "' stellar component of the '" + self.previous_name + "' model", default=True): continue

            # Adjust?
            adjust = prompt_yn("adjust_" + name, "adjust the properties of the '" + name + "' stellar component of the '" + self.previous_name + "' model definition")

            # Adjust the component
            if adjust: self.adjust_stellar_component(name)

            # Don't adjust the component
            else: self.set_previous_stellar_component(name)

    # -----------------------------------------------------------------

    def adjust_stellar_component(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.info("Adjusting the '" + name + "' stellar component ...")

        # Load the component
        component = self.previous.load_stellar_component(name, add_map=False)

        # Adjust
        adjusted, map_path = adjust_component(component, return_map_path=True)

        # Set reference to the previous model for this dust component
        if not adjusted: self.set_previous_stellar_component(name)
        else:

            # Write the component
            component_path = write_component(self.model_stellar_path, name, component)

            # Set the appropriate path
            if is_bulge_name(name): self.bulge_path = component_path
            elif is_old_name(name):
                self.old_path = component_path
                self.old_stars_map_path = map_path
            elif is_young_name(name):
                self.young_path = component_path
                self.young_stars_map_path = map_path
            elif is_ionizing_name(name):
                self.ionizing_path = component_path
                self.ionizing_stars_map_path = map_path
            else:
                self.extra_stellar_paths.append(component_path)
                self.extra_stellar_paths.append(map_path)

    # -----------------------------------------------------------------

    def set_previous_stellar_component(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.info("Setting reference to the previous model for the '" + name + "' stellar component ...")

        # Get the component path
        path = self.previous.get_stellar_component_path(name)

        # Get the map path
        map_path = self.previous.get_stellar_component_map_path(name)

        # Set the appropriate attribute
        if is_bulge_name(name): self.bulge_path = path
        elif is_old_name(name):
            self.old_path = path
            self.old_stars_map_path = map_path
        elif is_young_name(name):
            self.young_path = path
            self.young_stars_map_path = map_path
        elif is_ionizing_name(name):
            self.ionizing_path = path
            self.ionizing_stars_map_path = map_path
        else:
            self.extra_stellar_paths.append(path)
            self.extra_stellar_map_paths.append(map_path)

    # -----------------------------------------------------------------

    def adjust_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting dust components from previous model [" + self.previous_name + "] ...")

        # Loop over the dust components
        for name in self.previous.dust_component_names:

            # Include?
            if not prompt_yn("include_" + name, "include the '" + name + "' dust component of the '" + self.previous_name + "' model", default=True): continue

            # Adjust?
            adjust = prompt_yn("adjust_" + name, "adjust the '" + name + "' dust component of the '" + self.previous_name + "' model definition")

            # Adjust the component
            if adjust: self.adjust_dust_component(name)

            # Don't adjust
            else: self.set_previous_dust_component(name)

    # -----------------------------------------------------------------

    def adjust_dust_component(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.info("Adjusting the '" + name + "' dust component ...")

        # Load the component
        component = self.previous.load_dust_component(name, add_map=False)

        # Adjust
        adjusted, map_path = adjust_component(component, return_map_path=True)

        # Set reference to the previous model for this dust component
        if not adjusted: self.set_previous_dust_component(name)
        else:

            # Write the component, get the path
            component_path = write_component(self.model_dust_path, name, component)

            # Set the appropriate path
            if is_dust_disk_name(name):
                self.dust_path = component_path
                self.dust_map_path = map_path
            else:
                self.extra_dust_paths.append(component_path)
                self.extra_dust_paths.append(map_path)

    # -----------------------------------------------------------------

    def set_previous_dust_component(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.info("Setting reference to the previous model for the '" + name + "' dust component ...")

        # Get the component path
        path = self.previous.get_dust_component_path(name)

        # Get the map path
        map_path = self.previous.get_dust_component_map_path(name)

        # Set the appropriate attribute
        if is_dust_disk_name(name):
            self.dust_path = path
            self.dust_map_path = map_path
        else:
            self.extra_dust_paths.append(path)
            self.extra_dust_map_paths.append(map_path)

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

    @property
    def has_all_basic_stellar(self):

        """
        This function ...
        :return:
        """

        return self.has_bulge and self.has_old and self.has_young and self.has_ionizing

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

        # Set default options for components thare are already present so they won't be prompted for anymore
        use_default = None
        if self.has_bulge:
            if use_default is None: use_default = bulge_default_parameters
            else: use_default.extend(bulge_default_parameters)
        if self.has_old:
            if use_default is None: use_default = old_default_parameters
            else: use_default.extend(old_default_parameters)
        if self.has_young:
            if use_default is None: use_default = young_default_parameters
            else: use_default.extend(young_default_parameters)
        if self.has_ionizing:
            if use_default is None: use_default = ionizing_default_parameters
            else: use_default.extend(ionizing_default_parameters)
        if self.has_all_basic_stellar:
            if use_default is None: use_default = general_stellar_default_parameters
            else: use_default.extend(general_stellar_default_parameters)

        # Create the builder
        builder = StarsBuilder(interactive=True, cwd=self.config.path, config=config, prompt_optional=True, use_default=use_default)

        # Run
        builder.run()

        # Set the scaleheight of the old stars
        self.old_scaleheight = builder.old_scaleheight

        # Set bulge paths
        if not self.has_bulge:
            self.bulge_path = builder.bulge_path

        # Set old stellar component paths
        if not self.has_old:
            self.old_path = builder.old_stars_path
            self.old_stars_map_path = builder.old_stars_map_path

        # Set young stellar component paths
        if not self.has_young:
            self.young_path = builder.young_stars_path
            self.young_stars_map_path = builder.young_stars_map_path

        # Set ionizing stellar component paths
        if not self.has_ionizing:
            self.ionizing_path = builder.ionizing_stars_path
            self.ionizing_stars_map_path = builder.ionizing_stars_map_path

        # Set additional stellar component paths
        for name in builder.additional_names:

            # Get the paths
            path = builder.paths[name]
            map_path = builder.map_paths[name] if name in builder.map_paths else None

            # Set the paths
            self.extra_stellar_paths.append(path)
            self.extra_stellar_map_paths.append(map_path)

    # -----------------------------------------------------------------

    @property
    def has_dust(self):

        """
        This function ...
        :return:
        """

        return self.dust_path is not None

    # -----------------------------------------------------------------

    @property
    def has_all_basic_dust(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_dust

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

        # Set default options for components thare are already present so they won't be prompted for anymore
        use_default = None
        if self.has_dust:
            if use_default is None: use_default = dust_default_parameters
            else: use_default.extend(dust_default_parameters)

        # Create the builder
        builder = DustBuilder(interactive=True, cwd=self.config.path, config=config, prompt_optional=True, use_default=use_default)

        # Run
        builder.run(old_scaleheight=self.old_scaleheight)

        # Set dust disk paths
        if not self.has_dust:
            self.dust_path = builder.dust_path
            self.dust_map_path = builder.dust_map_path

        # Set additional dust component paths
        for name in builder.additional_names:

            # Get the paths
            path = builder.paths[name]
            map_path = builder.map_paths[name] if name in builder.map_paths else None

            # Set the paths
            self.extra_dust_paths.append(path)
            self.extra_dust_map_paths.append(map_path)

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

    @property
    def nextra_stellar(self):

        """
        Thisf unction ...
        :return:
        """

        return len(self.extra_stellar_paths)

    # -----------------------------------------------------------------

    @property
    def has_extra_stellar(self):

        """
        This function ...
        :return:
        """

        return self.nextra_stellar > 0

    # -----------------------------------------------------------------

    @property
    def nextra_dust(self):

        """
        Thisf unction ...
        :return:
        """

        return len(self.extra_dust_paths)

    # -----------------------------------------------------------------

    @property
    def has_extra_dust(self):

        """
        This function ...
        :return:
        """

        return self.nextra_dust > 0

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the models table ...")

        # Set extra paths
        extra_stellar_paths = self.extra_stellar_paths if self.has_extra_stellar else None
        extra_dust_paths = self.extra_dust_paths if self.has_extra_dust else None

        # Add the model
        table = self.models_table
        # name, description, bulge_path, old_stars_path, young_stars_path, ionizing_stars_path, dust_path,
        # additional_stellar_paths=None, additional_dust_paths=None
        table.add_model(self.model_name, self.config.description, self.bulge_path, self.old_path, self.young_path,
                        self.ionizing_path, self.dust_path, extra_stellar_paths, extra_dust_paths)

        # Save the table
        table.saveto(self.models_table_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the model composition ...")

        # Show stellar
        self.show_stellar()

        # Show dust
        self.show_dust()

    # -----------------------------------------------------------------

    def show_stellar(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the model stellar components ...")

        print("")
        print(fmt.green + fmt.underlined + "STELLAR COMPONENTS" + fmt.reset)
        print("")

        # Show the components
        if self.has_bulge: self.show_bulge()
        if self.has_old: self.show_old()
        if self.has_young: self.show_young()
        if self.has_ionizing: self.show_ionizing()
        self.show_extra_stellar()

    # -----------------------------------------------------------------

    def show_bulge(self):

        """
        This function ...
        :return:
        """

        print("  " + fmt.magenta + "OLD STELLAR BULGE:" + fmt.reset)
        print("")

        # Show
        show_component(self.bulge_path, line_prefix="    ")
        print("")

    # -----------------------------------------------------------------

    def show_old(self):

        """
        This function ...
        :return:
        """

        print("  " + fmt.magenta + "OLD STELLAR DISK: " + fmt.reset)
        print("")

        # Show
        show_component(self.old_path, self.old_stars_map_path, line_prefix="    ")
        print("")

    # -----------------------------------------------------------------

    def show_young(self):

        """
        This function ...
        :return:
        """

        print("  " + fmt.magenta + "YOUNG STELLAR DISK:" + fmt.reset)
        print("")

        # Show
        show_component(self.young_path, self.young_stars_map_path, line_prefix="    ")
        print("")

    # -----------------------------------------------------------------

    def show_ionizing(self):

        """
        This function ...
        :return:
        """

        print("  " + fmt.magenta + "IONIZING STELLAR DISK: " + fmt.reset)
        print("")

        # Show
        show_component(self.ionizing_path, self.ionizing_stars_map_path, line_prefix="    ")
        print("")

    # -----------------------------------------------------------------

    @property
    def extra_stellar_names(self):

        """
        This function ...
        :return:
        """

        for path in self.extra_stellar_paths:
            name = fs.name(path)
            yield name

    # -----------------------------------------------------------------

    def index_for_extra_stellar(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return list(self.extra_stellar_names).index(name)

    # -----------------------------------------------------------------

    def paths_for_extra_stellar(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        index = self.index_for_extra_stellar(name)
        return self.extra_stellar_paths[index], self.extra_stellar_map_paths[index]

    # -----------------------------------------------------------------

    def show_extra_stellar(self):

        """
        This function ...
        :return:
        """

        # Loop over the extra stellar components
        for name in self.extra_stellar_paths:

            # Show the name
            print("  " + fmt.magenta + name.upper() + ":" + fmt.reset)
            print("")

            # Get paths
            path, map_path = self.paths_for_extra_stellar(name)

            # Show the component
            show_component(path, map_path, line_prefix="    ")
            print("")

    # -----------------------------------------------------------------

    def show_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the model dust components ...")

        print("")
        print(fmt.green + fmt.underlined + "DUST COMPONENTS" + fmt.reset)
        print("")

        # Show the components
        if self.has_dust: self.show_dust_disk()
        self.show_extra_dust()

    # -----------------------------------------------------------------

    def show_dust_disk(self):

        """
        This function ...
        :return:
        """

        print("  " + fmt.magenta + "DUST DISK:" + fmt.reset)
        print("")

        # Show component
        show_component(self.dust_path, self.dust_map_path, line_prefix="    ")
        print("")

    # -----------------------------------------------------------------

    @property
    def extra_dust_names(self):

        """
        This function ...
        :return:
        """

        for path in self.extra_dust_paths:
            name = fs.name(path)
            yield name

    # -----------------------------------------------------------------

    def index_for_extra_dust(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return list(self.extra_dust_names).index(name)

    # -----------------------------------------------------------------

    def paths_for_extra_dust(self, name):

        """
        Thisf function ...
        :param name:
        :return:
        """

        index = self.index_for_extra_dust(name)
        return self.extra_dust_paths[index], self.extra_dust_map_paths[index]

    # -----------------------------------------------------------------

    def show_extra_dust(self):

        """
        This function ...
        :return:
        """

        # Loop over the extra dust components
        for name in self.extra_dust_names:

            # Show the name
            print("  " + fmt.magenta + name.upper() + ": " + fmt.reset)
            print("")

            # Get paths
            path, map_path = self.paths_for_extra_dust(name)

            # Show
            show_component(path, map_path, line_prefix="    ")
            print("")

# -----------------------------------------------------------------

def adjust_component(component, return_map_path=False):

    """
    Thisj function ...
    :param component:
    :param return_map_path:
    :return:
    """

    adjusted = False
    map_path = None

    # Parameters
    if "parameters" in component:

        # Adjust?
        adjust = prompt_yn("adjust_parameters", "adjust the parameters", default=False)
        if adjust: adjusted_parameters = adjust_parameters(component.parameters)
        else: adjusted_parameters = False

        # Update adjusted flag
        if adjusted_parameters: adjusted = True

    # Deprojection
    if "deprojection" in component:

        # Adjust?
        adjust = prompt_yn("adjust_deprojection", "adjust the deprojection", default=False)
        if adjust: adjusted_deprojection, deprojection_map_path = adjust_deprojection(component.deprojection, return_map_path=True)
        else:
            adjusted_deprojection = False
            deprojection_map_path = None # TODO: THIS IS NOT GOOD!

        # Update adjusted flag
        if adjusted_deprojection: adjusted = True
        if deprojection_map_path is not None: map_path = deprojection_map_path

    # Model
    if "model" in component:

        # Adjust?
        adjust = prompt_yn("adjust_model", "adjust the model", default=False)
        if adjust: adjusted_model = adjust_model(component.model)
        else: adjusted_model = False

        # Update adjusted flag
        if adjusted_model: adjusted = True

    # Properties
    if "properties" in component:

        # Adjust?
        adjust = prompt_yn("adjust_properties", "adjust the properties", default=False)
        if adjust: adjusted_properties, properties_map_path = adjust_properties(component.properties, return_map_path=True)
        else:
            adjusted_properties = False
            properties_map_path = None # TODO: THIS IS NOT GOOD!

        # Update adjusted flag
        if adjusted_properties: adjusted = True
        if properties_map_path is not None: map_path = properties_map_path

    # Return adjusted flag
    if return_map_path: return adjusted, map_path
    else: return adjusted

# -----------------------------------------------------------------

def adjust_parameters(parameters):

    """
    This function ...
    :param parameters:
    :return:
    """

    adjusted = False

    # Loop over the parameters
    for name in parameters:

        default = parameters[name]
        value = prompt_automatic(name, "adjusted '" + name + "' value", default)
        if value != default: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------

def adjust_deprojection(deprojection, return_map_path=False):

    """
    This function ...
    :param deprojection:
    :param return_map_path:
    :return:
    """

    adjusted = False
    map_path = None

    # Loop over the properties
    for name in deprojection:

        default = deprojection[name]
        ptype = deprojection.get_ptype(name)
        value = prompt_variable(name, ptype, "adjusted '" + name + "' value", default=default)
        if value != default: adjusted = True

        # If this is a filename, assume it's for the map
        if name.lower() == "filename": map_path = value

    # Return the adjusted flag
    if return_map_path: return adjusted, map_path
    else: return adjusted

# -----------------------------------------------------------------

def adjust_model(model):

    """
    This function ...
    :param model:
    :return:
    """

    adjusted = False

    # Loop over the properties
    for name in model:

        default = model[name]
        ptype = model.get_ptype(name)
        value = prompt_variable(name, ptype, "adjusted '" + name + "' value", default=default)
        if value != default: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------

def adjust_properties(properties, return_map_path=False):

    """
    This function ...
    :param properties:
    :param return_map_path:
    :return:
    """

    adjusted = False
    map_path = None

    # Geometry
    adjust_geometry = prompt_yn("adjust_geometry", "adjust geometry")
    if adjust_geometry: adjusted_geometry, geometry_map_path = adjust_geometry_properties(properties["geometry"], return_map_path=True)
    else:
        adjusted_geometry = False
        geometry_map_path = None # TODO: THIS IS NOT RIGHT!

    # Set flag
    if adjusted_geometry: adjusted = True
    if geometry_map_path is not None: map_path = geometry_map_path

    # SED (stellar properties)
    if "sed" in properties:

        adjust_sed = prompt_yn("adjust_sed", "adjust SED properties")
        if adjust_sed: adjusted_sed = adjust_sed_properties(properties["sed"])
        else: adjusted_sed = False

        # Set flag
        if adjusted_sed: adjusted = True

    # Mix (dust properties)
    elif "mix" in properties:

        adjust_mix = prompt_yn("adjust_mix", "adjust mix properties")
        if adjust_mix: adjusted_mix = adjust_mix_properties(properties["mix"])
        else: adjusted_mix = False

        # Set flag
        if adjusted_mix: adjusted = True

    # Neither SED nor mix
    else: raise ValueError("Properties does not contain 'sed' or 'mix'")

    # Normalization
    adjust_normalization = prompt_yn("adjust_normalization", "adjust normalization properties")
    if adjust_normalization: adjusted_normalization = adjust_normalization_properties(properties["normalization"])
    else: adjusted_normalization = False

    # Set flag
    if adjusted_normalization: adjusted = True

    # Return the adjusted flag
    if return_map_path: return adjusted, map_path
    else: return adjusted

# -----------------------------------------------------------------

def adjust_geometry_properties(properties, return_map_path=False):

    """
    This function ...
    :param properties:
    :param return_map_path:
    :return:
    """

    adjusted = False
    map_path = None

    # Loop over the properties
    for name in properties:

        default = properties[name]
        value = prompt_automatic(name, "adjusted '" + name + "' value", default)
        if value != default: adjusted = True

        # Check for filename property, assume it's for a map
        if name.lower() == "filename": map_path = value

    # Return the adjusted flag
    if return_map_path: return adjusted, map_path
    else: return adjusted

# -----------------------------------------------------------------

def adjust_sed_properties(properties):

    """
    This function ...
    :param properties:
    :return:
    """

    adjusted = False

    # Loop over the properties
    for name in properties:

        default = properties[name]
        value = prompt_automatic(name, "adjusted '" + name + "' value", default)
        if value != default: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------

def adjust_mix_properties(properties):

    """
    This function ...
    :param properties:
    :return:
    """

    adjusted = False

    # Loop over the properties
    for name in properties:

        default = properties[name]
        value = prompt_automatic(name, "adjusted '" + name + "' value", default)
        if value != default: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------

def adjust_normalization_properties(properties):

    """
    This function ...
    :param properties:
    :return:
    """

    adjusted = False

    # Loop over the properties
    for name in properties:

        default = properties[name]
        value = prompt_automatic(name, "adjusted '" + name + "' value", default)
        if value != default: adjusted = True

    # Return the adjusted flag
    return adjusted

# -----------------------------------------------------------------

def is_bulge_name(name):

    """
    This function ...
    :param name:
    :return:
    """

    return name == bulge_component_name

# -----------------------------------------------------------------

def is_old_name(name):

    """
    This function ...
    :param name:
    :return:
    """

    return name == old_component_name

# -----------------------------------------------------------------

def is_young_name(name):

    """
    This function ...
    :param name:
    :return:
    """

    return name == young_component_name

# -----------------------------------------------------------------

def is_ionizing_name(name):

    """
    This function ...
    :param name:
    :return:
    """

    return name == ionizing_component_name

# -----------------------------------------------------------------

def is_dust_disk_name(name):

    """
    This function ...
    :param name:
    :return:
    """

    return name == disk_component_name

# -----------------------------------------------------------------

def show_component(path, map_path=None, line_prefix=""):

    """
    This function ...
    :param path:
    :param map_path:
    :param line_prefix:
    :return:
    """

    from ..suite import load_component

    # Load the component
    component = load_component(path, add_map=False)

    # Parameters
    if "parameters" in component:

        print(line_prefix + "PARAMETERS:")
        print("")

        for key in component.parameters:
            value = component.parameters[key]
            print(line_prefix + " - " + fmt.bold + key + fmt.reset + ": " + tostr(value))

        print("")

    # Deprojection
    if "deprojection" in component:

        print(line_prefix + "DEPROJECTION:")
        print("")

        for key in component.deprojection:
            value = component.deprojection[key]
            print(line_prefix + " - " + fmt.bold + key + fmt.reset + ": " + tostr(value))

        print("")

    # Model
    if "model" in component:

        print(line_prefix + "MODEL:")
        print("")

        for key in component.model:
            value = component.model[key]
            print(line_prefix + " - " + fmt.bold + key + fmt.reset + ": " + tostr(value))

        print("")

    # Properties
    if "properties" in component:

        # Geometry
        print(line_prefix + "GEOMETRY PROPERTIES:")
        print("")

        for key in component.properties["geometry"]:
            value = component.properties["geometry"][key]
            print(line_prefix + " - " + fmt.bold + key + fmt.reset + ": " + tostr(value))
        print("")

        # SED
        if "sed" in component.properties:

            print(line_prefix + "SED PROPERTIES:")
            print("")

            for key in component.properties["sed"]:
                value = component.properties["sed"][key]
                print(line_prefix + " - " + fmt.bold + key + fmt.reset + ": " + tostr(value))
            print("")

        # Mix
        elif "mix" in component.properties:

            print(line_prefix + "MIX PROPERTIES:")
            print("")

            for key in component.properties["mix"]:
                value = component.properties["mix"][key]
                print(line_prefix + " - " + fmt.bold + key + fmt.reset + ": " + tostr(value))
            print("")

        # Neither SED nor mix
        else: raise ValueError("Neither 'sed' or 'mix' is defined in the component properties")

        # Normalization
        print(line_prefix + "NORMALIZATION PROPERTIES:")
        print("")

        for key in component.properties["normalization"]:
            value = component.properties["normalization"][key]
            print(line_prefix + " - " + fmt.bold + key + fmt.reset + ": " + tostr(value))
        print("")

# -----------------------------------------------------------------
