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
from ....core.data.sun import Sun
from ....core.basics.log import log
from ...component.galaxy import GalaxyModelingComponent
from .base import FittingInitializerBase

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

        # Solar luminosity units
        self.sun_fuv = None
        self.sun_i1 = None

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

        # Solar properties
        sun = Sun()
        self.sun_fuv = sun.luminosity_for_filter_as_unit(self.fuv_filter) # Get the luminosity of the Sun in the FUV band
        self.sun_i1 = sun.luminosity_for_filter_as_unit(self.i1_filter)   # Get the luminosity of the Sun in the IRAC I1 band

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

    def create_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grids ...")

        # Fixed wavelengths (always in the grid)
        fixed = [self.i1_filter.pivot, self.fuv_filter.pivot]

        # Set options
        self.wg_generator.config.show = False
        self.wg_generator.config.write = False

        # Generate the wavelength grids
        self.wg_generator.run(npoints_range=self.config.wg.npoints_range, ngrids=self.config.wg.ngrids,
                              fixed=fixed, add_emission_lines=self.config.wg.add_emission_lines,
                              min_wavelength=self.config.wg.range.min, max_wavelength=self.config.wg.range.max,
                              filters=self.fitting_run.fitting_filters)

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

        # 4. Write the paths to the input maps
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
