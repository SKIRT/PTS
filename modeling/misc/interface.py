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
from ...core.tools.logging import log
from ...core.basics.emissionlines import EmissionLines
from ...core.prep.wavelengthgrids import create_one_subgrid_wavelength_grid
from .select import select_from_model_suite, select_from_fitting_context
from ...core.prep.dustgrids import create_one_dust_grid_for_galaxy_from_deprojection
from ..build.representation import create_projections_from_dust_grid, create_projections_from_deprojections
from ...core.basics.configuration import prompt_yn

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
        model_name, ski, definition, input_paths, parameter_values = select_from_model_suite(self.model_suite)

        # Set attributes
        self.ski = ski
        self.definition = definition
        self.parameter_values = parameter_values
        self.input_paths = input_paths

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

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Fixed wavelengths (always in the grid: for normalization)
        fixed = [self.i1_filter.pivot, self.fuv_filter.pivot]

        # Create the emission lines instance
        if self.config.wg.add_emission_lines: emission_lines = EmissionLines()
        else: emission_lines = None

        # Create the grid
        # grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added
        grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added = \
            create_one_subgrid_wavelength_grid(self.config.wg.npoints, emission_lines, fixed,
                                               min_wavelength=self.config.wg.range.min,
                                               max_wavelength=self.config.wg.range.max)

        # Set the grid
        self.wavelength_grid = grid

        # Debugging
        log.debug("Created a wavelength grid with:")
        log.debug("")
        log.debug(" - number of points: " + str(len(grid)))
        log.debug(" - number of points in subgrids: " + str(subgrid_npoints))
        log.debug(" - number of emission points: " + str(emission_npoints))
        log.debug(" - number of fixed points: " + str(fixed_npoints))
        log.debug(" - filters for which extra sampling was performed: " + str(broad_resampled))
        log.debug(" - narrow band filters for which wavelength was added: " + str(narrow_added))
        log.debug("")
        if log.is_debug():
            print("")
            print(self.wavelength_grid)
            print("")

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

        # Set minimum level
        if self.config.dg.grid_type == "bintree": min_level = self.config.dg.bintree_max_level
        elif self.config.dg.grid_type == "octtree": min_level = self.config.dg.octtree_max_level
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
            title, deprojection = self.model_suite.load_stellar_component_deprojection(self.model_name, name)
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
        for name in self.model_suite.get_dust_component_names(self.config.path, self.model_name):

            # Load the deprojection of the component, if applicable
            title, deprojection = self.model_suite.load_dust_component_deprojection(self.model_name, name)
            if deprojection is not None: self.deprojections[(name, title)] = deprojection

    # -----------------------------------------------------------------

    def create_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the projection systems ...")

        azimuth = 0.0

        # Use grid?
        if prompt_yn("grid_resolution", "use the resolution of the dust grid for setting up the instruments?"):

            earth, faceon, edgeon = create_projections_from_dust_grid(self.dust_grid, self.galaxy_distance,
                                                                      self.galaxy_inclination, azimuth,
                                                                      self.disk_position_angle)

        # Use deprojections
        else:
            earth, faceon, edgeon = create_projections_from_deprojections(self.deprojections, self.galaxy_distance,
                                                                          azimuth)

        # Set the projection systems
        self.projections["earth"] = earth
        self.projections["faceon"] = faceon
        self.projections["edgeon"] = edgeon

    # -----------------------------------------------------------------

    @property
    def earth_projection(self):

        """
        This function ...
        :return:
        """

        return self.projections["earth"]

    # -----------------------------------------------------------------

    @property
    def faceon_projection(self):

        """
        This function ...
        :return:
        """

        return self.projections["faceon"]

    # -----------------------------------------------------------------

    @property
    def edgeon_projection(self):

        """
        This function ...
        :return:
        """

        return self.projections["edgeon"]

    # -----------------------------------------------------------------

    @abstractproperty
    def instrument_class(self):

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
        self.instruments["earth"] = self.instrument_class.from_projection(self.earth_projection)

        # Create a faceon instrument
        self.instruments["faceon"] = self.instrument_class.from_projection(self.faceon_projection)

        # Create an edgeon instrument
        self.instruments["edgeon"] = self.instrument_class.from_projection(self.edgeon_projection)

    # -----------------------------------------------------------------

    @property
    def earth_instrument(self):

        """
        This function ...
        :return:
        """

        return self.instruments["earth"]

    # -----------------------------------------------------------------

    @property
    def faceon_instrument(self):

        """
        This function ...
        :return:
        """

        return self.instruments["faceon"]

    # -----------------------------------------------------------------

    @property
    def edgeon_instrument(self):

        """
        Thsi function ...
        :return:
        """

        return self.instruments["edgeon"]

# -----------------------------------------------------------------
