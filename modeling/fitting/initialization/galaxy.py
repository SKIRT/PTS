#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.initialization.galaxy Contains the GalaxyFittingInitializer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ...component.galaxy import GalaxyModelingComponent
from .base import FittingInitializerBase
from ....core.prep.wavelengthgrids import WavelengthGridGenerator
from ....core.basics.emissionlines import get_id_strings, important_lines
from ....core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class GalaxyFittingInitializer(FittingInitializerBase, GalaxyModelingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructors of the base classes
        FittingInitializerBase.__init__(self, no_config=True)
        GalaxyModelingComponent.__init__(self, *args, **kwargs)

        # The wavelength grid generators
        self.basic_generator = None
        self.refined_generator = None
        self.highres_generator = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the ski file
        self.load_ski()

        # 3. Load the model representation
        self.load_representation()

        # 4. Set the stellar and dust components
        self.set_components()

        # 5. Create the wavelength grids
        self.create_wavelength_grids()

        # 6. Adjust the ski template
        self.adjust_ski()

        # 7. Calculate the weight factor to give to each band
        self.calculate_weights()

        # 8. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        FittingInitializerBase.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the template ski file ...")

        # Load the ski template
        self.ski = self.fitting_run.ski_template

    # -----------------------------------------------------------------

    def set_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the stellar and dust components ...")

        # Add the components to the ski file and to the input map paths dictionary
        self.suite.add_model_components(self.model_name, self.ski, self.input_map_paths)

    # -----------------------------------------------------------------

    @property
    def fitting_filters(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.fitting_filters

    # -----------------------------------------------------------------

    def create_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grids ...")

        # Basic wavelength grids
        self.create_basic_wavelength_grids()

        # Refined wavelength grids
        self.create_refined_wavelength_grids()

        # High-resolution wavelength grids
        self.create_highres_wavelength_grids()

    # -----------------------------------------------------------------

    @property
    def important_emission_lines(self):

        """
        This function ...
        :return:
        """

        return important_lines

    # -----------------------------------------------------------------

    @lazyproperty
    def all_emission_lines(self):

        """
        Thins function ...
        :return:
        """

        return get_id_strings()

    # -----------------------------------------------------------------

    def create_basic_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating basic wavelength grids ...")

        # Create generator
        self.basic_generator = WavelengthGridGenerator()

        # Set basic settings
        self.basic_generator.config.ngrids = self.config.ngrids_basic
        self.basic_generator.config.npoints_range = self.config.wg.npoints_range_basic
        self.basic_generator.config.range = self.config.wg.range

        # Set other
        self.basic_generator.config.add_emission_lines = self.config.wg.add_emission_lines
        self.basic_generator.config.emission_lines = self.important_emission_lines
        self.basic_generator.config.check_filters = self.observed_filter_wavelengths_no_iras_planck
        self.basic_generator.config.adjust_to = self.fitting_filters
        self.basic_generator.config.fixed = self.normalization_wavelengths

        # Set other flags
        self.basic_generator.config.show = False
        self.basic_generator.config.write = False
        self.basic_generator.config.plot = True

        # Generate the wavelength grids
        self.basic_generator.run(table=self.wg_table)

    # -----------------------------------------------------------------

    def create_refined_wavelength_grids(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Creating refined wavelength grids ...")

        # Create generator
        self.refined_generator = WavelengthGridGenerator()

        # Set basic settings
        self.refined_generator.config.ngrids = self.config.ngrids_refined
        self.refined_generator.config.npoints_range = self.config.wg.npoints_range_refined
        self.refined_generator.config.range = self.config.wg.range

        # Set other
        self.refined_generator.config.add_emission_lines = self.config.wg.add_emission_lines
        self.refined_generator.config.emission_lines = self.important_emission_lines
        self.refined_generator.config.check_filters = self.observed_filter_wavelengths_no_iras_planck
        self.refined_generator.config.filters = self.fitting_filters
        self.refined_generator.config.fixed = self.normalization_wavelengths

        # Set other flags
        self.refined_generator.config.show = False
        self.refined_generator.config.write = False
        self.refined_generator.config.plot = True

        # Generate the refined wavelength grids
        self.refined_generator.run(table=self.wg_table)

    # -----------------------------------------------------------------

    def create_highres_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating high-resolution wavelength grids ...")

        # Create generator
        self.highres_generator = WavelengthGridGenerator()

        # Set basic settings
        self.highres_generator.config.ngrids = self.config.ngrids_highres
        self.highres_generator.config.npoints_range = self.config.wg.npoints_range_highres
        self.highres_generator.config.range = self.config.wg.range

        # Set other
        self.highres_generator.config.add_emission_lines = self.config.wg.add_emission_lines
        self.highres_generator.config.emission_lines = self.all_emission_lines
        self.highres_generator.config.check_filters = self.observed_filter_wavelengths_no_iras_planck
        self.highres_generator.config.filters = self.fitting_filters
        self.highres_generator.config.fixed = self.normalization_wavelengths

        # Set other flags
        self.highres_generator.config.show = False
        self.highres_generator.config.write = False
        self.highres_generator.config.plot = True

        # Generate the high-resolution wavelength grids
        self.highres_generator.run(table=self.wg_table)

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add the instrument
        self.ski.add_instrument("earth", self.representation.sed_instrument)

        # Set the number of photon packages
        self.ski.setpackages(self.config.npackages)

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid("wavelengths.txt")

        # Set the dust emissivityex
        self.set_dust_emissivity()

        # Set the dust grid
        self.set_dust_grid()

        # Set all-cells dust library
        self.ski.set_allcells_dust_lib()

        # Set the dust selfabsorption
        self.set_selfabsorption()

        # Disable all writing options
        self.ski.disable_all_writing_options()

    # -----------------------------------------------------------------

    def set_dust_emissivity(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust emissivity ...")

        # Set dust emissivity if present
        if self.ski.has_dust_emissivity:

            # Enable or disable
            if self.config.transient_heating: self.ski.set_transient_dust_emissivity()
            else: self.ski.set_grey_body_dust_emissivity()

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

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 1. Write the ski file
        self.write_ski()

        # 2. Write the weights table
        self.write_weights()

        # 3. Write the wavelength grids
        self.write_wavelength_grids()

        # 4. Write the wavelength grid table
        self.write_wavelength_grid_table()

        # 5. Write the paths to the input maps
        self.write_input_map_paths()

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file to " + self.fitting_run.template_ski_path + " ...")

        # Save the ski template file
        self.ski.save()

    # -----------------------------------------------------------------

    def write_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grids ...")

        # Basic
        self.write_basic_wavelength_grids()

        # Refined
        self.write_refined_wavelength_grids()

        # Highres
        self.write_highres_wavelength_grids()

    # -----------------------------------------------------------------

    def write_basic_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the basic wavelength grids ...")

    # -----------------------------------------------------------------

    def write_refined_wavelength_grids(self):

        """
        Thins function ...
        :return:
        """

        # Inform the user
        log.info("Writing the refined wavelength grids ...")

    # -----------------------------------------------------------------

    def write_highres_wavelength_grids(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the high-resolution wavelength grids ...")

# -----------------------------------------------------------------
