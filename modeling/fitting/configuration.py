#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.configuration Contains the FittingConfigurer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.simulation.skifile import SkiFile
from ...core.basics.log import log
from ..config.parameters import definition as parameters_definition
from ...core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, Configuration
from ...core.basics.configuration import DictConfigurationSetter, combine_configs
from ..config.parameters import parsing_types_for_parameter_types, unit_parsing_type
from ..config.parameters import default_units, possible_parameter_types_descriptions
from .run import FittingRun
from ...core.basics.configuration import prompt_string
from ..build.component import get_representations_for_model, get_representation_path, get_representation
from ..build.representation import Representation
from ...core.tools.stringify import tostr
from ...evolve.solve.extremizer import genetic_definition
from ...core.tools.utils import lazyproperty
from ...core.basics.range import QuantityRange
from ..config.parameters import distance_name, ionizing_scaleheight_name, sfr_compactness_name, fuv_young_name
from ..config.parameters import old_scaleheight_name, position_angle_name, dust_mass_name, fuv_ionizing_name
from ..config.parameters import metallicity_name, young_scaleheight_name, sfr_covering_name, dust_scaleheight_name
from ..config.parameters import i1_old_name, sfr_pressure_name, inclination_name

# -----------------------------------------------------------------

template_ski_path = fs.join(introspection.pts_modeling_ski_templates_path(), "labeled_template.ski")

# -----------------------------------------------------------------

# Default magnitude ranges
magnitude_ranges = dict()
magnitude_ranges["fuv_young"] = 2
magnitude_ranges["dust_mass"] = 1
magnitude_ranges["fuv_ionizing"] = 3

# -----------------------------------------------------------------

# Default percentual deviation from the initial value
percentual_ranges = dict()
percentual_ranges["distance"] = 20
percentual_ranges["ionizing_scaleheight"] = 25
percentual_ranges["sfr_compactness"] = 25
percentual_ranges["old_scaleheight"] = 25
percentual_ranges["position_angle"] = 10
percentual_ranges["metallicity"] = 50
percentual_ranges["young_scaleheight"] = 25
percentual_ranges["sfr_covering"] = 25
percentual_ranges["dust_scaleheight"] = 25
percentual_ranges["i1_old"] = 10
percentual_ranges["sfr_pressure"] = 25
percentual_ranges["inclination"] = 15

# -----------------------------------------------------------------

default_magnitude_range = 2
default_percentual_range = 20

# -----------------------------------------------------------------

class FittingConfigurer(FittingComponent):
    
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
        super(FittingConfigurer, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The fitting run
        self.fitting_run = None

        # The initial model representation
        self.initial_representation = None

        # The default ranges
        self.default_ranges = dict()

        # The ski file template
        self.ski = None

        # The individual configurations
        self.parameters_config = None
        self.descriptions_config = None
        self.types_config = None
        self.units_config = None
        self.ndigits_config = None
        self.ranges_config = None
        self.filters_config = None
        self.genetic_config = None
        self.grid_config = None

        # Additional settings
        self.settings = None

        # The final fitting config
        self.fitting_config = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the necessary input
        self.load_input()

        # 3. Get the fitting parameters
        self.set_parameters()

        # 4. Get parameter descriptions
        self.set_descriptions()

        # 5. Get parameter types
        self.set_types()

        # 6. Get parameter units
        self.set_units()

        # 7. Get the physical parameter ranges
        self.set_ranges()

        # 8. Get the number of significant digits
        self.set_ndigits()

        # 9. Get the fitting filters
        self.set_filters()

        # 10. Get the settings for the genetic algorithm
        self.set_method()

        # 11. Create the fitting configuration
        self.create_config()

        # 12. Adjust the labels of the template ski file
        self.adjust_labels()

        # 13. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingConfigurer, self).setup()

        # Create the fitting run
        self.create_fitting_run()

        # Set the initial model representation
        self.set_representation()

        # Get the default ranges
        self.default_ranges = kwargs.pop("default_ranges", dict())

        # Get configs as dicts
        if "parameters_config" in kwargs: self.parameters_config = kwargs.pop("parameters_config")
        if "descriptions_config" in kwargs: self.descriptions_config = kwargs.pop("descriptions_config")
        if "types_config" in kwargs: self.types_config = kwargs.pop("types_config")
        if "units_config" in kwargs: self.units_config = kwargs.pop("units_config")
        if "ndigits_config" in kwargs: self.ndigits_config = kwargs.pop("ndigits_config")
        if "ranges_config" in kwargs: self.ranges_config = kwargs.pop("ranges_config")
        if "filters_config" in kwargs: self.filters_config = kwargs.pop("filters_config")
        if "genetic_config" in kwargs: self.genetic_config = kwargs.pop("genetic_config")
        if "grid_config" in kwargs: self.grid_config = kwargs.pop("grid_config")

        # Set settings dict
        if "settings" in kwargs: self.settings = kwargs.pop("settings")

    # -----------------------------------------------------------------

    @lazyproperty
    def previous_run(self):

        """
        This function ...
        :return:
        """

        if self.config.from_run is None: return None
        else: return self.fitting_runs.load(self.config.from_run)

    # -----------------------------------------------------------------

    @lazyproperty
    def previous_configuration(self):

        """
        This function ...
        :return:
        """

        return self.previous_run.fitting_configuration

    # -----------------------------------------------------------------

    @property
    def from_previous_run(self):

        """
        This function ...
        :return:
        """

        return self.config.from_run is not None

    # -----------------------------------------------------------------

    @property
    def parameters_from_previous_run(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.adapt_parameters and not self.config.except_parameters

    # -----------------------------------------------------------------

    @property
    def previous_parameters(self):

        """
        This function ...
        :return:
        """

        return self.previous_configuration.free_parameters

    # -----------------------------------------------------------------

    @property
    def previous_parameters_definition(self):

        """
        Thisnf unction ...
        :return:
        """

        return ConfigurationDefinition.from_defaults(self.previous_parameters)

    # -----------------------------------------------------------------

    @property
    def parameters_definition(self):

        """
        This function ...
        :return:
        """

        if self.from_previous_run and not self.config.except_parameters: return self.previous_parameters_definition
        else: return parameters_definition

    # -----------------------------------------------------------------

    @property
    def prompt_optional_parameters(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.except_parameters

    # -----------------------------------------------------------------

    @property
    def descriptions_from_previous_run(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.adapt_descriptions and not self.config.except_descriptions

    # -----------------------------------------------------------------

    @property
    def previous_descriptions(self):

        """
        Thisf unction ...
        :return:
        """

        return self.previous_configuration.descriptions

    # -----------------------------------------------------------------

    @property
    def previous_descriptions_definition(self):

        """
        This function ...
        :return:
        """

        return ConfigurationDefinition.from_defaults(self.previous_descriptions)

    # -----------------------------------------------------------------

    @property
    def descriptions_definition(self):

        """
        This function ...
        :return:
        """

        if self.from_previous_run and not self.config.except_descriptions:
            definition = self.previous_descriptions_definition
            #prompt_optional = True
        else:

            # Create configuration definition
            definition = ConfigurationDefinition()
            for name in self.parameter_labels: definition.add_required(name, "string", "description of the '" + name + "' parameter")

            # Set flag
            #prompt_optional = False

        # Return
        return definition

    # -----------------------------------------------------------------

    @property
    def prompt_optional_descriptions(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.except_descriptions

    # -----------------------------------------------------------------

    @property
    def types_from_previous_run(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.adapt_types and not self.config.except_types

    # -----------------------------------------------------------------

    @property
    def previous_types(self):

        """
        This function ...
        :return:
        """

        return self.previous_configuration.types

    # -----------------------------------------------------------------

    @property
    def previous_types_definition(self):

        """
        This function ...
        :return:
        """

        return ConfigurationDefinition.from_defaults(self.previous_types)

    # -----------------------------------------------------------------

    @property
    def types_definition(self):

        """
        Thisf unction ...
        :return:
        """

        if self.from_previous_run and not self.config.except_types:
            definition = self.previous_types_definition
            #prompt_optional = True
        else:

            # Create definition
            definition = ConfigurationDefinition()
            for name in self.parameter_labels: definition.add_required(name, "string", "type of the '" + name + "' parameter", choices=possible_parameter_types_descriptions)

            # Set flag
            #prompt_optional = False

        # Return
        return definition

    # -----------------------------------------------------------------

    @property
    def prompt_optional_types(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.except_types

    # -----------------------------------------------------------------

    @property
    def units_from_previous_run(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.adapt_units and not self.config.except_units

    # -----------------------------------------------------------------

    @property
    def previous_units(self):

        """
        This function ...
        :return:
        """

        return self.previous_configuration.units

    # -----------------------------------------------------------------

    @property
    def previous_units_definition(self):

        """
        This function ...
        :return:
        """

        return ConfigurationDefinition.from_defaults(self.previous_units)

    # -----------------------------------------------------------------

    @property
    def units_definition(self):

        """
        This function ...
        :return:
        """

        if self.from_previous_run and not self.config.except_units:
            definition = self.previous_units_definition
            #prompt_optional = True
        else:

            # Create definition
            definition = ConfigurationDefinition()
            for name in self.parameter_labels:

                # Get the type of quantity for this parameter
                parameter_type = self.types_config.types[name]

                # Don't ask for units for dimensionless quantities
                if parameter_type == "dimless": definition.add_fixed(name, name + " has no unit (dimensionless)", None)
                else: definition.add_optional(name, unit_parsing_type(parameter_type), "unit of the '" + name + "' parameter", default=default_units[parameter_type])

            # Set flag
            #prompt_optional = False

        # Return
        return definition

    # -----------------------------------------------------------------

    @property
    def prompt_optional_units(self):

        """
        Thisn function ...
        :return:
        """

        return self.from_previous_run and not self.config.except_units

    # -----------------------------------------------------------------

    @property
    def ranges_from_previous_run(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.adapt_ranges and not self.config.except_ranges

    # -----------------------------------------------------------------

    @property
    def previous_ranges(self):

        """
        This function ...
        :return:
        """

        ranges = Configuration()

        for label in self.parameter_labels:

            # Get range
            key = label + "_range"

            # Get value
            value = self.previous_configuration[key]

            # Set value
            ranges[key] = value

        # Return
        return ranges

    # -----------------------------------------------------------------

    @property
    def previous_ranges_definition(self):

        """
        This function ...
        :return:
        """

        return ConfigurationDefinition.from_defaults(self.previous_ranges)

    # -----------------------------------------------------------------

    @property
    def ranges_definition(self):

        """
        This function ...
        :return:
        """

        if self.from_previous_run and not self.config.except_ranges: definition = self.previous_ranges_definition
        else:

            # Create the configuration
            definition = ConfigurationDefinition(write_config=False)

            # Add the options for the ranges
            for label in self.parameter_labels:

                # Use free ranges
                if self.config.free_ranges:

                    # Get the unit
                    unit = self.units_config.units[label]
                    units_info_string = " (don't forget the units!) " if unit is not None else " (dimensionless)"

                    # Get the parsing type
                    parsing_type = parsing_types_for_parameter_types[self.types_config.types[label]]

                    # Get the description
                    description = "range of " + self.description_for_parameter(label) + units_info_string

                    # Get the default range
                    default_range = self.default_ranges[label] if label in self.default_ranges else None

                    # Make optional or required depending on whether default is given
                    if default_range is not None: definition.add_optional(label + "_range", parsing_type + "_range", description, default=default_range, convert_default=True)
                    else: definition.add_required(label + "_range", parsing_type + "_range", description)

                # Use magnitudes
                elif label in magnitude_ranges:

                    # Get default magnitude
                    default_magnitude = magnitude_ranges[label]

                    # Make optional setting
                    definition.add_optional(label + "_magnitude", "real", "magnitude range for " + self.description_for_parameter(label), default=default_magnitude)

                # Use percentage
                elif label in percentual_ranges:

                    # Get default
                    default_percentage = float(percentual_ranges[label]) / 100.

                    # Make optional setting
                    definition.add_optional(label + "_percentage", "percentage", "percentual deviation from initial value for " + self.description_for_parameter(label), default=default_percentage)

                # Use default percentage
                else:

                    # Get default
                    default_percentage = float(default_percentual_range) / 100.

                    # Make optional setting
                    definition.add_optional(label + "_percentage", "percentage", "percentual deviation from initial value for " + self.description_for_parameter(label), default=default_percentage)

        # Return
        return definition

    # -----------------------------------------------------------------

    @property
    def prompt_optional_ranges(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.except_ranges

    # -----------------------------------------------------------------

    @property
    def ndigits_from_previous_run(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.adapt_ndigits and not self.config.except_ndigits

    # -----------------------------------------------------------------

    @property
    def previous_ndigits(self):

        """
        This function ...
        :return:
        """

        return self.previous_configuration.ndigits

    # -----------------------------------------------------------------

    @property
    def previous_ndigits_definition(self):

        """
        This function ...
        :return:
        """

        return ConfigurationDefinition.from_defaults(self.previous_ndigits)

    # -----------------------------------------------------------------

    @property
    def ndigits_definition(self):

        """
        This function ...
        :return:
        """

        if self.from_previous_run and not self.config.except_ndigits: definition = self.previous_ndigits_definition
        else:

            # Create definition
            definition = ConfigurationDefinition()
            for name in self.parameter_labels:

                # Ask for the number of significant digits
                definition.add_optional(name, "positive_integer", "number of significant digits for the '" + name + "' parameter", default=self.config.default_ndigits)

        # Return
        return definition

    # -----------------------------------------------------------------

    @property
    def prompt_optional_ndigits(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def filters_from_previous_run(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.adapt_filters and not self.config.except_filters

    # -----------------------------------------------------------------

    @property
    def previous_filters(self):

        """
        This function ...
        :return:
        """

        return self.previous_configuration.filters

    # -----------------------------------------------------------------

    @property
    def previous_filters_definition(self):

        """
        This function ...
        :return:
        """

        return ConfigurationDefinition.from_defaults(self.previous_filters)

    # -----------------------------------------------------------------

    @property
    def filters_definition(self):

        """
        This function ...
        :return:
        """

        if self.from_previous_run and not self.config.except_filters:
            definition = self.previous_filters_definition
            #prompt_optional = True
        else:

            # Create the configuration
            definition = ConfigurationDefinition(write_config=False)

            # Choose from all the possible filter names
            definition.add_required("filters", "string_list", "the filters for which to use the observed flux as reference for the fitting procedure", choices=self.sed_filter_names)

            # Set flag
            #prompt_optional = False

        # Return
        return definition

    # -----------------------------------------------------------------

    @property
    def prompt_optional_filters(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.except_filters

    # -----------------------------------------------------------------

    @property
    def genetic_settings_from_previous_run(self):

        """
        This function ...
        :return:
        """

        return self.from_previous_run and not self.config.adapt_genetic and not self.config.except_genetic

    # -----------------------------------------------------------------

    @property
    def previous_genetic(self):

        """
        This function ...
        :return:
        """

        return self.previous_configuration.genetic

    # -----------------------------------------------------------------

    @property
    def previous_genetic_definition(self):

        """
        This function ...
        :return:
        """

        return ConfigurationDefinition.from_defaults(self.previous_genetic)

    # -----------------------------------------------------------------

    @property
    def genetic_definition(self):

        """
        This function ...
        :return:
        """

        if self.from_previous_run and not self.config.except_genetic: return self.previous_genetic_definition
        else: return genetic_definition

    # -----------------------------------------------------------------

    @property
    def prompt_optional_genetic(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def grid_settings_from_previous_run(self):

        """
        Thisf unction ...
        :return:
        """

        return self.previous_run and not self.config.adapt_grid and not self.config.except_grid

    # -----------------------------------------------------------------

    @property
    def previous_grid(self):

        """
        This function ...
        :return:
        """

        return self.previous_configuration.grid

    # -----------------------------------------------------------------

    @property
    def previous_grid_definition(self):

        """
        This function ...
        :return:
        """

        return ConfigurationDefinition.from_defaults(self.previous_grid)

    # -----------------------------------------------------------------

    @property
    def grid_definition(self):

        """
        This function ...
        :return:
        """

        if self.from_previous_run and not self.config.except_grid: return self.previous_grid_definition
        else:

            # Create
            definition = ConfigurationDefinition(write_config=False)

            # Loop over the fitting parameters
            for label in self.parameter_labels:

                # Linear or logarithmic
                definition.add_required(label + "_scale", "string", "scale for " + self.description_for_parameter(label), choices=["linear", "logarithmic"])

                # Number of grid points
                definition.add_positional_optional(label + "_npoints", "positive_integer", "number of grid points for " + self.description_for_parameter(label), default=10)

            # Return
            return definition

    # -----------------------------------------------------------------

    @property
    def prompt_optional_grid(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    def create_fitting_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the fitting run ...")

        # Create the run
        self.fitting_run = FittingRun(self.config.path, self.config.name, self.config.model_name)

        # Create the run directory
        fs.create_directory(self.fitting_run.path)

    # -----------------------------------------------------------------

    def set_representation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the initial model representation ...")

        # Dictionary for the options
        options = dict()

        highest_pixelscale = None
        highest_pixelscale_name = None

        # Loop over the different representations for the model
        for name in get_representations_for_model(self.config.path, self.model_name):

            # Get the pixelscale
            #pixelscale = get_pixelscale_for_representation(self.config.path, name).average

            # Get the representation
            representation = get_representation(self.config.path, name)
            pixelscale = representation.earth_pixelscale.average

            # Determine name and description
            option = name
            if highest_pixelscale is None or pixelscale > highest_pixelscale:
                highest_pixelscale = pixelscale
                highest_pixelscale_name = name
            #description = "representation '" + name + "' with: " + tostr(representation.properties)
            description = "representation with instrument resolution of " + tostr(representation.earth_pixelscale.average.to("arcsec"), round=True, ndigits=3) + " with a " + representation.dust_grid_type + " dust grid"
            if representation.has_dust_grid_tree_distribution: description += " with " + str(representation.ndust_cells) + " dust cells (min level = " + str(representation.dust_grid_min_level) + ", max level = " + str(representation.dust_grid_max_level) + ")"

            # Add the option
            options[option] = description

        # Get the answer
        name = prompt_string("representation", "initial representation to use for fitting the model", choices=options, default=highest_pixelscale_name)

        # Set the initial representation
        path = get_representation_path(self.config.path, name)
        self.initial_representation = Representation(name, self.model_name, path)

    # -----------------------------------------------------------------

    @property
    def run_name(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.name

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.model_name

    # -----------------------------------------------------------------

    @lazyproperty
    def model_definition(self):

        """
        This function ...
        :return:
        """

        from ..build.component import get_model_definition
        return get_model_definition(self.config.path, self.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_fuv_young(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.young_stars_luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_fuv_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.ionizing_stars_luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.dust_mass

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_distance(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.distance

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_old_scaleheight(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.old_stars_scaleheight

    # -----------------------------------------------------------------

    @property
    def initial_young_scaleheight(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.young_stars_scaleheight

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_ionizing_scaleheight(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.ionizing_stars_scaleheight

    # -----------------------------------------------------------------

    @property
    def initial_dust_scaleheight(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.dust_scaleheight

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_sfr_compactness(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.ionizing_stars_compactness

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_position_angle(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.position_angle

    # -----------------------------------------------------------------

    @property
    def initial_metallicity(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.metallicity

    # -----------------------------------------------------------------

    @property
    def initial_sfr_covering(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.ionizing_stars_covering_factor

    # -----------------------------------------------------------------

    @property
    def initial_i1_old(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.old_stars_luminosity

    # -----------------------------------------------------------------

    @property
    def initial_sfr_pressure(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.ionizing_stars_pressure

    # -----------------------------------------------------------------

    @property
    def initial_inclination(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.inclination

    # -----------------------------------------------------------------

    def get_initial_value(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        if label == "distance": return self.initial_distance
        elif label == "ionizing_scaleheight": return self.initial_ionizing_scaleheight
        elif label == "sfr_compactness": return self.initial_sfr_compactness
        elif label == "fuv_young": return self.initial_fuv_young
        elif label == "old_scaleheight": return self.initial_old_scaleheight
        elif label == "position_angle": return self.initial_position_angle
        elif label == "dust_mass": return self.initial_dust_mass
        elif label == "fuv_ionizing": return self.initial_fuv_ionizing
        elif label == "metallicity": return self.initial_metallicity
        elif label == "young_scaleheight": return self.initial_young_scaleheight
        elif label == "sfr_covering": return self.initial_sfr_covering
        elif label == "dust_scaleheight": return self.initial_dust_scaleheight
        elif label == "i1_old": return self.initial_i1_old
        elif label == "sfr_pressure": return self.initial_sfr_pressure
        elif label == "inclination": return self.initial_inclination
        else: raise ValueError("Unrecognized parameter label")

    # -----------------------------------------------------------------

    @property
    def run_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.path

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the input ...")

        # Load the template ski file
        self.load_template()

    # -----------------------------------------------------------------

    def load_template(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file template ...")

        # Load ski file template
        self.ski = SkiFile(template_ski_path)

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the free parameters ...")

        # Load or prompt for the parameters
        if self.config.parameters is not None: self.parameters_config = Configuration(free_parameters=self.config.parameters)
        elif isinstance(self.parameters_config, dict): pass
        elif self.parameters_from_previous_run: self.parameters_config = Configuration(free_parameters=self.previous_parameters)
        else: self.prompt_parameters()

    # -----------------------------------------------------------------

    def prompt_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the free parameters ...")

        # Create configuration setter
        setter = InteractiveConfigurationSetter("Free parameters", add_logging=False, add_cwd=False)

        # Create config
        self.parameters_config = setter.run(self.parameters_definition, prompt_optional=self.prompt_optional_parameters)

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return: 
        """

        return self.parameters_config.free_parameters

    # -----------------------------------------------------------------

    def set_descriptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter descriptions ...")

        # Load or prompt for the descriptions
        if self.config.descriptions is not None: self.descriptions_config = Configuration(descriptions=self.config.descriptions)
        elif isinstance(self.descriptions_config, dict): pass
        elif self.descriptions_from_previous_run: self.descriptions_config = Configuration(descriptions=self.previous_descriptions)
        else: self.prompt_descriptions()

    # -----------------------------------------------------------------

    def prompt_descriptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the parameter descriptions ...")

        # Create the configuration setter
        setter = InteractiveConfigurationSetter("Parameter descriptions", add_cwd=False, add_logging=False)

        # Get the config and set the descriptions configuration
        config = setter.run(self.descriptions_definition, prompt_optional=self.prompt_optional_descriptions)
        self.descriptions_config = Configuration(descriptions=config)

    # -----------------------------------------------------------------

    def description_for_parameter(self, label):

        """
        This function ...
        :param label: 
        :return: 
        """

        return self.descriptions_config.descriptions[label]

    # -----------------------------------------------------------------

    def set_types(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter types ...")

        # load or prompt for the tyopes
        if self.config.types is not None: self.types_config = Configuration(types=self.config.types)
        elif isinstance(self.types_config, dict): pass
        elif self.types_from_previous_run: self.types_config = Configuration(types=self.previous_types)
        else: self.prompt_types()

    # -----------------------------------------------------------------

    def prompt_types(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the parameter types ...")

        # Create configuration setter
        setter = InteractiveConfigurationSetter("Parameter types", add_cwd=False, add_logging=False)

        # Create the config and set the types configuration
        config = setter.run(self.types_definition, prompt_optional=self.prompt_optional_types)
        self.types_config = Configuration(types=config)

    # -----------------------------------------------------------------

    def set_units(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter units ...")

        # load or prompt for the units
        if self.config.units is not None: self.units_config = Configuration(units=self.config.units)
        elif isinstance(self.units_config, dict): pass
        elif self.units_from_previous_run: self.units_config = Configuration(units=self.previous_units)
        else: self.prompt_units()

    # -----------------------------------------------------------------

    def prompt_units(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the parameter units ...")

        # Create configuration setter
        setter = InteractiveConfigurationSetter("Parameter units", add_cwd=False, add_logging=False)

        # Create the config and set the units configuration
        config = setter.run(self.units_definition, prompt_optional=self.prompt_optional_units)
        self.units_config = Configuration(units=config)

    # -----------------------------------------------------------------

    def set_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter ranges ...")

        # Ranges are given
        if len(self.config.ranges) != 0:

            self.ranges_config = dict()
            for label in self.config.ranges:
                self.ranges_config[label + "_range"] = self.config.ranges[label]

        elif isinstance(self.ranges_config, dict): pass
        elif self.ranges_from_previous_run: self.ranges_config = self.previous_ranges
        else: self.prompt_ranges()

    # -----------------------------------------------------------------

    def prompt_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the parameter ranges ...")

        # Create configuration setter
        setter = InteractiveConfigurationSetter("free parameter ranges", add_logging=False, add_cwd=False)

        # Create config, get the range for each chosen free parameter
        self.ranges_config = setter.run(self.ranges_definition, prompt_optional=self.prompt_optional_ranges)

        # Convert magnitudes or percentages into actual ranges
        for label in self.parameter_labels:

            # Determine key
            magnitude_key = label + "_magnitude"
            percentage_key = label + "_percentage"
            range_key = label + "_range"

            # Check magnitude
            if magnitude_key in self.ranges_config:

                # Get the magnitude
                magnitude = self.ranges_config[magnitude_key]

                # Get the initial value for this parameter
                value = self.get_initial_value(label)

                # Construct the actual range
                parameter_range = QuantityRange.within_magnitude(value, magnitude)

                # Set the range
                self.ranges_config[range_key] = parameter_range

            # Check percentage
            elif percentage_key in self.ranges_config:

                # Get the percentage
                fraction = self.ranges_config[percentage_key]

                # Get the initial value for this parameter
                value = self.get_initial_value(label)

                # Construct the actual range
                rel_min = 1. - fraction
                rel_max = 1. + fraction
                parameter_range = QuantityRange.around(value, rel_min=rel_min, rel_max=rel_max)

                # Set the range
                self.ranges_config[range_key] = parameter_range

            # Not magnitude nor percentage
            else:

                if range_key not in self.ranges_config: raise ValueError("No range defined for parameter '" + label + "'")
                continue

    # -----------------------------------------------------------------

    def set_ndigits(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting the number of significant digits ...")

        # load or prompt for the ndigits
        if self.config.ndigits is not None: self.ndigits_config = Configuration(ndigits=self.config.ndigits)
        elif isinstance(self.ndigits_config, dict): pass
        elif self.ndigits_from_previous_run: self.ndigits_config = Configuration(ndigits=self.previous_ndigits)
        else: self.prompt_ndigits()

    # -----------------------------------------------------------------

    def prompt_ndigits(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Prompting for the number of significant digits ...")

        # Create configuration setter
        setter = InteractiveConfigurationSetter("Parameter ndigits", add_cwd=False, add_logging=False)

        # Create the config and set the ndigits configuration
        config = setter.run(self.ndigits_definition, prompt_optional=self.prompt_optional_ndigits)
        self.ndigits_config = Configuration(ndigits=config)

    # -----------------------------------------------------------------

    def set_filters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the fitting filters ...")

        # Load or prompt for the filters
        if self.config.filters is not None: self.filters_config = Configuration(filters=self.config.filters)
        elif isinstance(self.filters_config, dict): pass
        elif self.filters_from_previous_run: self.filters_config = Configuration(filters=self.previous_filters)
        else: self.prompt_filters()

    # -----------------------------------------------------------------

    def prompt_filters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the filters ...")

        # Create configuration setter
        setter = InteractiveConfigurationSetter("filters", add_logging=False, add_cwd=False)

        # Create config, get the filter choices
        self.filters_config = setter.run(self.filters_definition, prompt_optional=self.prompt_optional_filters)

    # -----------------------------------------------------------------

    @property
    def grid_fitting(self):

        """
        This function ...
        :return:
        """

        return self.config.fitting_method == "grid"

    # -----------------------------------------------------------------

    @property
    def genetic_fitting(self):

        """
        Thisn function ...
        :return:
        """

        return self.config.fitting_method == "genetic"

    # -----------------------------------------------------------------

    def set_method(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Setting options depending on the fitting method ...")

        # Set options depending on the fitting method
        if self.grid_fitting: self.set_grid()
        elif self.genetic_fitting: self.set_genetic()
        else: raise ValueError("Invalid fitting method: must be 'grid' or 'genetic'")

    # -----------------------------------------------------------------

    def set_genetic(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the genetic algorithm configuration ...")

        # Load or prompt for the settings
        if self.config.genetic is not None:
            setter = DictConfigurationSetter(self.config.genetic, "genetic")
            configuration = setter.run(self.genetic_definition)
            self.genetic_config = Configuration(genetic=configuration)
        elif isinstance(self.genetic_config, dict): pass
        elif self.genetic_settings_from_previous_run: self.genetic_config = Configuration(genetic=self.previous_genetic)
        else: self.prompt_genetic()

        # Set the grid configuration to None
        self.grid_config = Configuration(grid=None)

    # -----------------------------------------------------------------

    def prompt_genetic(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the settings of the genetic algorithm ...")

        # Create configuration setter
        setter = InteractiveConfigurationSetter("genetic", add_logging=False, add_cwd=False)

        # Create config, get the choices
        config = setter.run(self.genetic_definition, prompt_optional=self.prompt_optional_genetic)
        self.genetic_config = Configuration(genetic=config)

    # -----------------------------------------------------------------

    def set_grid(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting the grid fitting configuration ...")

        # Load or prompt for the settings
        if self.config.grid is not None:
            setter = DictConfigurationSetter(self.config.grid, "grid")
            configuration = setter.run(self.grid_definition)
            self.grid_config = Configuration(grid=configuration)
        elif isinstance(self.grid_config, dict): pass
        elif self.grid_settings_from_previous_run: self.grid_config = Configuration(grid=self.previous_grid)
        else: self.prompt_grid()

        # Set the genetic configuration to None
        self.genetic_config = Configuration(genetic=None)

    # -----------------------------------------------------------------

    def prompt_grid(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Prompting for the grid fitting settings ...")

        # Create configuration setter
        setter = InteractiveConfigurationSetter("grid", add_logging=False, add_cwd=False)

        # Create config, get the choices
        config = setter.run(self.grid_definition, prompt_optional=self.prompt_optional_grid)
        self.grid_config = Configuration(grid=config)

    # -----------------------------------------------------------------

    def set_settings(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting additional settings ...")

        # Prompt settings
        if self.settings is None: self.prompt_settings()

    # -----------------------------------------------------------------

    def prompt_settings(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for additional settings ...")

        # ...
        self.settings = dict()

    # -----------------------------------------------------------------

    def create_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the fitting run configuration ...")

        # Combine configs
        self.fitting_config = combine_configs(self.parameters_config, self.descriptions_config, self.types_config,
                                              self.units_config, self.ndigits_config, self.ranges_config,
                                              self.filters_config, self.genetic_config, self.grid_config)

        # NEW: Set the fitting method !!
        self.fitting_config.method = self.config.fitting_method

        # Set the name of the initial representation
        self.fitting_config.initial_representation = self.initial_representation.name

        # Set additional settings
        if self.settings is not None:
            for label in self.settings: self.fitting_config[label] = self.settings[label]

    # -----------------------------------------------------------------

    def adjust_labels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting labels of the free parameters in the template ski file ...")

        # Loop over the labels currently in the template ski file
        for label in self.ski.labels:

            # If the label is in the list of free parameter labels, skip (don't remove)
            if label in self.parameters_config.free_parameters: continue

            # Otherwise, delabel
            self.ski.delabel(label)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the fitting configuration
        self.write_config()

        # Write the ski file template
        self.write_ski()

        # Write the runs table
        self.write_table()

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the fitting configuration ...")

        # Write the configuration
        self.fitting_config.saveto(self.fitting_run.fitting_configuration_path)

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the template ski file ...")

        # Save the ski file template
        self.ski.saveto(self.fitting_run.template_ski_path)

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the runs table ...")

        # Open the runs table, add the new run and save
        table = self.runs_table
        table.add_run(self.fitting_run)
        table.save()

# -----------------------------------------------------------------

def set_definition_values(model_definition, values):

    """
    This function ...
    :param model_definition:
    :param values:
    :return:
    """

    # Loop over the parameters
    for label in values: set_definition_value(model_definition, label, values[label])

# -----------------------------------------------------------------

def set_definition_value(model_definition, label, value):

    """
    This function ...
    :param model_definition:
    :param label:
    :param value:
    :return:
    """

    if label == distance_name: model_definition.distance = value
    elif label == ionizing_scaleheight_name: model_definition.ionizing_stars_scaleheight = value
    elif label == sfr_compactness_name: model_definition.ionizing_stars_compactness = value
    elif label == fuv_young_name: model_definition.young_stars_luminosity = value
    elif label == old_scaleheight_name: model_definition.old_stars_scaleheight = value
    elif label == position_angle_name: model_definition.position_angle = value
    elif label == dust_mass_name: model_definition.dust_mass = value
    elif label == fuv_ionizing_name: model_definition.ionizing_stars_luminosity = value
    elif label == metallicity_name: model_definition.metallicity = value
    elif label == young_scaleheight_name: model_definition.young_stars_scaleheight = value
    elif label == sfr_covering_name: model_definition.ionizing_stars_covering_factor = value
    elif label == dust_scaleheight_name: model_definition.dust_scaleheight = value
    elif label == i1_old_name: model_definition.old_stars_luminosity = value
    elif label == sfr_pressure_name: model_definition.ionizing_stars_pressure = value
    elif label == inclination_name: model_definition.inclination = value
    else: raise ValueError("Unrecognized parameter label")

# -----------------------------------------------------------------

def get_definition_value(model_definition, label):

    """
    This function ...
    :param model_definition:
    :param label:
    :return:
    """

    if label == distance_name: return model_definition.distance
    elif label == ionizing_scaleheight_name: return model_definition.ionizing_stars_scaleheight
    elif label == sfr_compactness_name: return model_definition.ionizing_stars_compactness
    elif label == fuv_young_name: return model_definition.young_stars_luminosity
    elif label == old_scaleheight_name: return model_definition.old_stars_scaleheight
    elif label == position_angle_name: return model_definition.position_angle
    elif label == dust_mass_name: return model_definition.dust_mass
    elif label == fuv_ionizing_name: return model_definition.ionizing_stars_luminosity
    elif label == metallicity_name: return model_definition.metallicity
    elif label == young_scaleheight_name: return model_definition.young_stars_scaleheight
    elif label == sfr_covering_name: return model_definition.ionizing_stars_covering_factor
    elif label == dust_scaleheight_name: return model_definition.dust_scaleheight
    elif label == i1_old_name: return model_definition.old_stars_luminosity
    elif label == sfr_pressure_name: return model_definition.ionizing_stars_pressure
    elif label == inclination_name: return model_definition.inclination
    else: raise ValueError("Unrecognized parameter label")

# -----------------------------------------------------------------
