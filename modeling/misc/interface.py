#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.misc.interface Contains the ModelSimulationInterface class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractproperty, abstractmethod

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.basics.log import log
from .select import select_from_model_suite, select_from_fitting_context
from ...core.prep.dustgrids import create_one_dust_grid_for_galaxy_from_deprojection
from ..build.representations.galaxy import create_projections_from_dust_grid, create_projections_from_deprojections
from ...core.basics.configuration import prompt_yn
from ..build.wavelengthgrid import WavelengthGridBuilder

# -----------------------------------------------------------------

earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

class ModelSimulationInterface(GalaxyModelingComponent):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ModelSimulationInterface, self).__init__(*args, **kwargs)

        # The ski file
        self.ski = None

        # The model definition
        self.definition = None

        # The parameter values
        self.parameter_values = None

        # The dictionary of input map paths
        self.input_paths = None

        # The wavelength grid
        self.wavelength_grid = None

        # The dust grid
        self.dust_grid = None

        # The deprojections
        self.deprojections = dict()

        # The projections and instruments
        self.projections = dict()
        self.instruments = dict()

    # -----------------------------------------------------------------

    @abstractmethod
    def get_model(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def prompt_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the model from the model suite ...")

        # Select the model
        model_name, ski, definition, input_paths, parameter_values = select_from_model_suite(self.model_suite, adapt=self.config.adapt, name=self.config.model_name)

        # Set attributes
        self.ski = ski
        self.definition = definition
        self.parameter_values = parameter_values
        self.input_paths = input_paths

        # Set the 'distance' parameter value, since instruments are still not adapted from the default template.
        # Instruments will only be added to the ski file later, so the parameter value obtained from the 'distance'
        # label in 'select_from_model_suite' is still incorrect
        # Other instrument properties should have been fixed (with the 'fix_labels' function)
        self.parameter_values["distance"] = self.galaxy_distance

        # Return the model name
        return model_name

    # -----------------------------------------------------------------

    def prompt_fitting(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the model from the fitting context ...")

        # Select the model
        run_id, generation_name, simulation_name, fitting_run, chi_squared, ski, input_paths, parameter_values = select_from_fitting_context(self.fitting_context)

        # Load the model definition
        definition = fitting_run.model_definition

        # Set attributes
        self.ski = ski
        self.definition = definition
        self.parameter_values = parameter_values
        self.input_paths = input_paths

        # Return the fitting run ID, generation name, simulation name, and chi_squared
        return run_id, generation_name, simulation_name, chi_squared

    # -----------------------------------------------------------------

    @property
    def normalization_filters(self):

        """
        This function ...
        :return:
        """

        return self.definition.normalization_filters

    # -----------------------------------------------------------------

    @property
    def normalization_wavelengths(self):

        """
        This function ...
        :return:
        """

        return self.definition.normalization_wavelengths

    # -----------------------------------------------------------------

    def create_wavelength_grid(self, check_filters=True, output_path=None, plot=False):

        """
        This function ...
        :param check_filters:
        :param output_path:
        :param plot:
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Create wavelength grid builder
        builder = WavelengthGridBuilder()

        # Set options
        builder.config.add_emission_lines = self.config.wg.add_emission_lines
        builder.config.min_wavelength = self.config.wg.range.min
        builder.config.max_wavelength = self.config.wg.range.max
        builder.config.check_filters = self.observed_filter_wavelengths_no_iras_planck if check_filters else None
        builder.config.npoints = self.config.wg.npoints
        builder.config.adjust_to = self.observed_filter_wavelengths_no_iras_planck #self.observed_filter_wavelengths
        builder.config.fixed = self.normalization_wavelengths

        # Set output path
        if output_path is not None: builder.config.write = True
        builder.config.output = output_path
        #builder.config.write_grid = False # DO write the grid in table format

        # Plot?
        builder.config.plot = plot

        # Run the wavelength grid builder
        builder.run()

        # Get the wavelength grid
        self.wavelength_grid = builder.wavelength_grid

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        return len(self.wavelength_grid)

    # -----------------------------------------------------------------

    def create_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust grid ...")

        # Load the dust disk deprojection
        deprojection = self.definition.dust_deprojection

        #print(self.config.dg.grid_type)

        # Set minimum level
        if self.config.dg.grid_type == "bintree": min_level = self.config.dg.bintree_min_level
        elif self.config.dg.grid_type == "octtree": min_level = self.config.dg.octtree_min_level
        else: min_level = None

        # Set max ndivisions per pixel
        max_ndivisions_per_pixel = 1. / self.config.dg.scale  # default 1/0.5 = 2 divisions along each direction per pixel

        # Create the dust grid
        # grid_type, deprojection, distance, sky_ellipse, min_level, max_mass_fraction, max_ndivisions_per_pixel=2, nscaleheights=10.
        self.dust_grid = create_one_dust_grid_for_galaxy_from_deprojection(self.config.dg.grid_type, deprojection,
                                                                           self.galaxy_distance,
                                                                           self.truncation_ellipse,
                                                                           min_level, self.config.dg.max_mass_fraction,
                                                                           max_ndivisions_per_pixel,
                                                                           self.config.dg.scale_heights)

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.definition.name

    # -----------------------------------------------------------------

    def load_deprojections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the deprojections ...")

        # Stellar
        self.load_stellar_deprojections()

        # Dust
        self.load_dust_deprojections()

    # -----------------------------------------------------------------

    def load_stellar_deprojections(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Loading the stellar deprojections ...")

        # Loop over the stellar components
        for name in self.model_suite.get_stellar_component_names(self.model_name):

            # Load the deprojection of the component, if applicable
            title, deprojection = self.model_suite.load_stellar_component_deprojection(self.model_name, name, load_map=True)
            if deprojection is not None: self.deprojections[(name, title)] = deprojection

    # -----------------------------------------------------------------

    def load_dust_deprojections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the dust deprojections ...")

        # Loop over the dust components
        for name in self.model_suite.get_dust_component_names(self.model_name):

            # Load the deprojection of the component, if applicable
            title, deprojection = self.model_suite.load_dust_component_deprojection(self.model_name, name, load_map=True)
            if deprojection is not None: self.deprojections[(name, title)] = deprojection

    # -----------------------------------------------------------------

    def create_projection_systems(self, make_edgeon=True, make_faceon=True, use_grid=None, reference_name=None):

        """
        This function ...
        :param make_edgeon:
        :param make_faceon:
        :param use_grid:
        :param reference_name:
        :return:
        """

        # Debugging
        log.debug("Creating the projection systems ...")

        azimuth = 0.0

        # Use grid resolution?
        if reference_name is not None:
            if use_grid is None: use_grid = False
            elif use_grid: raise ValueError("Cannot use '" + reference_name + "' as reference for the projection systems when 'use_grid' is enabled")
        if use_grid is None: use_grid = prompt_yn("grid_resolution", "use the resolution of the dust grid for setting up the instruments?", default=False)

        # Use grid?
        if use_grid:

            earth, faceon, edgeon = create_projections_from_dust_grid(self.dust_grid, self.galaxy_distance,
                                                                      self.disk_inclination, azimuth,
                                                                      self.disk_position_angle)
            deprojection_name = "grid"

        # Use deprojections
        else:

            # Create projection systems
            earth, faceon, edgeon, deprojection_name = create_projections_from_deprojections(self.deprojections, self.galaxy_distance,
                                                                                               azimuth, scale_heights=self.config.old_scale_heights,
                                                                                               return_deprojection_name=True,
                                                                                               scale_heights_reference="old",
                                                                                               reference_deprojection_name=reference_name,
                                                                                               radial_factor=self.config.radial_factor,
                                                                                               from_projection=self.config.from_projection)

        # Set the projection systems
        self.projections[earth_name] = earth
        if make_faceon: self.projections[faceon_name] = faceon
        if make_edgeon: self.projections[edgeon_name] = edgeon

        # Return the deprojection name, or 'grid' if grid resolution was used to create the instruments
        return deprojection_name

    # -----------------------------------------------------------------

    @property
    def earth_projection(self):

        """
        This function ...
        :return:
        """

        return self.projections[earth_name]

    # -----------------------------------------------------------------

    @property
    def faceon_projection(self):

        """
        This function ...
        :return:
        """

        if not self.has_faceon_projection: return None
        return self.projections[faceon_name]

    # -----------------------------------------------------------------

    @property
    def edgeon_projection(self):

        """
        This function ...
        :return:
        """

        if not self.has_edgeon_projection: return None
        return self.projections[edgeon_name]

    # -----------------------------------------------------------------

    @property
    def has_faceon_projection(self):

        """
        This function ...
        :return:
        """

        return faceon_name in self.projections

    # -----------------------------------------------------------------

    @property
    def has_edgeon_projection(self):

        """
        Thisfunction ...
        :return:
        """

        return edgeon_name in self.projections

    # -----------------------------------------------------------------

    @abstractproperty
    def instrument_class(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def earth_instrument_properties(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Create an earth instrument
        self.instruments[earth_name] = self.instrument_class.from_projection(self.earth_projection, **self.earth_instrument_properties)

        # Create a faceon instrument
        if self.has_faceon_projection: self.instruments[faceon_name] = self.instrument_class.from_projection(self.faceon_projection)

        # Create an edgeon instrument
        if self.has_edgeon_projection: self.instruments[edgeon_name] = self.instrument_class.from_projection(self.edgeon_projection)

    # -----------------------------------------------------------------

    @property
    def earth_instrument(self):

        """
        This function ...
        :return:
        """

        return self.instruments[earth_name]

    # -----------------------------------------------------------------

    @property
    def faceon_instrument(self):

        """
        This function ...
        :return:
        """

        if not self.has_faceon_instrument: return None
        return self.instruments[faceon_name]

    # -----------------------------------------------------------------

    @property
    def edgeon_instrument(self):

        """
        Thsi function ...
        :return:
        """

        if not self.has_edgeon_instrument: return None
        return self.instruments[edgeon_name]

    # -----------------------------------------------------------------

    @property
    def has_faceon_instrument(self):

        """
        Thisfunction ...
        :return:
        """

        return faceon_name in self.instruments

    # -----------------------------------------------------------------

    @property
    def has_edgeon_instrument(self):

        """
        This function ...
        :return:
        """

        return edgeon_name in self.instruments

    # -----------------------------------------------------------------

    @property
    def largest_instrument_npixels(self):

        """
        This function ...
        :return:
        """

        largest = None

        for name in self.instruments:

            npixels = self.instruments[name].npixels
            if npixels > largest: largest = npixels

        return largest

# -----------------------------------------------------------------
