#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.initialization.images Contains the ImagesFittingInitializer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ...component.images import ImagesModelingComponent
from ....core.simulation.skifile import SkiFile
from .base import FittingInitializerBase

# -----------------------------------------------------------------

class ImagesFittingInitializer(FittingInitializerBase, ImagesModelingComponent):
    
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
        ImagesModelingComponent.__init__(self, *args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the ski template
        self.load_ski_template()

        # 3. Load the model representation
        self.load_representation()

        # 3. Set the input map paths (in the dictionary and in the ski file template)
        self.set_input_map_paths()

        # 3. Create the wavelength grids
        self.create_wavelength_grids()

        # 5. Adjust the ski template
        self.adjust_ski()

        # 6. Calculate the weight factor to give to each band
        self.calculate_weights()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup functions of the base classes
        FittingInitializerBase.setup(self, **kwargs)
        ImagesModelingComponent.setup(self, **kwargs)

    # -----------------------------------------------------------------

    def load_ski_template(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file template ...")

        # Load the ski file
        self.ski = SkiFile(self.environment.ski_path)

    # -----------------------------------------------------------------

    def set_input_map_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the paths to the input maps ...")

        ## Stellar

        # Loop over the stellar components
        for name in self.suite.get_stellar_component_names(self.model_name):

            # Load the component
            component = self.suite.load_stellar_component(self.model_name, name)

            # Check whether requires a map
            if component.map_path is None: continue

            # Debugging
            log.debug("Setting the path to the '" + name + "' stellar component map ...")

            # Generate a filename for the map
            filename = "stars_" + name + ".fits"

            # Set the filename in the ski file template
            #self.ski.change_input_filename(old_filename, filename)

            # Set the filename in the ski file template
            self.ski.set_stellar_component_fits_geometry_filename(component.name, filename)

            # Add entry to the input maps dictionary
            self.input_map_paths[filename] = component.map_path

        ## Dust

        # Loop over the dust components
        for name in self.suite.get_dust_component_names(self.model_name):

            # Load the component
            component = self.suite.load_dust_component(self.model_name, name)

            # Check whether requires a map
            if component.map_path is None: continue

            # Debugging
            log.debug("Setting the path to the '" + name + "' dust component map ...")

            # Generate a filename for the map
            filename = "dust_" + name + ".fits"

            # Set the filename in the ski file template
            #self.ski.change_input_filename(old_filename, filename)

            # Set the filename in the ski file template
            self.ski.set_dust_component_fits_geometry_filename(component.name, filename)

            # Add entry to the input maps dictionary
            self.input_map_paths[filename] = component.map_path

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # 1. Set the instrument
        self.set_instrument()

        # 2. Set the number of photon packages
        self.ski.setpackages(self.config.npackages)

        # 3. Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid("wavelengths.txt")

        # 4. Set parameter values
        self.set_initial_parameter_values()

        # 5. Set dust emissivity
        self.set_dust_emissivity()

        # Set the dust grid
        self.set_dust_grid()

        # 6. Set all-cells dust library
        self.ski.set_allcells_dust_lib()

        # 7. Dust self-absorption
        self.set_selfabsorption()

        # 8. Disable all writing options
        self.ski.disable_all_writing_options()

    # -----------------------------------------------------------------

    def set_instrument(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the instrument ...")

        instrument_names = self.ski.get_instrument_names()
        if len(instrument_names) > 1: raise ValueError("The ski file can only have one instrument")
        old_name = instrument_names[0]
        self.ski.set_instrument_name(old_name, "earth")

    # -----------------------------------------------------------------

    def set_initial_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the initial parameter values ...")

        # Loop over the free parameter labels
        for label in self.fitting_run.free_parameter_labels:

            # Determine value as the mean of the given range
            value = self.fitting_run.free_parameter_ranges[label].mean

            # Set the value for the parameter
            self.ski.set_labeled_value(label, value)

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

        # Set
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

        # Write the wavelength grid table
        self.write_wavelength_grid_table()

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
        self.ski.saveto(self.fitting_run.template_ski_path)

# -----------------------------------------------------------------
