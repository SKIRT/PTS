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
from ....core.tools.logging import log
from ...component.galaxy import GalaxyModelingComponent
from ...build.component import get_stellar_component_names, get_dust_component_names, load_stellar_component, load_dust_component
from ....core.filter.filter import parse_filter
from .base import FittingInitializerBase

# -----------------------------------------------------------------

class GalaxyFittingInitializer(FittingInitializerBase, GalaxyModelingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        FittingInitializerBase.__init__(self, config, interactive)
        GalaxyModelingComponent.__init__(self, config, interactive)

        # The initial model representation
        self.representation = None

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

        self.ski = self.fitting_run.ski_template

    # -----------------------------------------------------------------

    def load_representation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the model representation ...")

        # Load the initial representation
        self.representation = self.fitting_run.initial_representation

    # -----------------------------------------------------------------

    def set_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the stellar and dust components ...")

        # 1. Set stellar components
        self.set_stellar_components()

        # 2. Set dust components
        self.set_dust_components()

    # -----------------------------------------------------------------

    def set_stellar_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the stellar components ...")

        # Loop over the stellar components
        for name in get_stellar_component_names(self.config.path, self.model_name):

            # Load the component
            component = load_stellar_component(self.config.path, self.model_name, name)

            # If an input map is required
            if "map_path" in component:

                # Generate a filename for the map
                filename = "stars_" + name + ".fits"

                # Set the filename
                if "deprojection" in component: component.deprojection.filename = filename
                elif "geometry" in component.parameters: component.properties["geometry"].filename = filename
                else: raise RuntimeError("Stellar component based on an input map should either have a deprojection or geometry properties")

                # Add entry to the input maps dictionary
                self.input_map_paths[filename] = component.map_path

            # Set geometry
            if "model" in component:

                # Get title
                title = component.parameters.title

                # Set the geometry
                self.ski.set_stellar_component_geometry(title, component.model)

            # Set deprojection
            elif "deprojection" in component:

                # Get title
                title = component.parameters.title

                # Set the deprojection geometry
                self.ski.set_stellar_component_geometry(title, component.deprojection)

            # Check if this is a new component, add geometry, SED and normalization all at once
            if "geometry" in component.parameters:

                # Get title
                title = component.parameters.title

                # Get class names
                geometry_type = component.parameters.geometry
                sed_type = component.parameters.sed
                normalization_type = component.parameters.normalization

                # Get properties for each of the three classes
                geometry_properties = component.properties["geometry"]
                sed_properties = component.properties["sed"]
                normalization_properties = component.properties["normalization"]

                # Create stellar component
                self.ski.create_new_stellar_component(title, geometry_type, geometry_properties, sed_type, sed_properties, normalization_type, normalization_properties)

            # Existing component, with MAPPINGS template
            elif "sfr" in component.parameters:

                # Get title
                title = component.parameters.title

                # Get SED properties
                metallicity = component.parameters.metallicity
                compactness = component.parameters.compactness
                pressure = component.parameters.pressure
                covering_factor = component.parameters.covering_factor

                # Get normalization
                fltr = parse_filter(component.parameters.filter)
                luminosity = component.parameters.luminosity

                # Set SED
                self.ski.set_stellar_component_mappingssed(title, metallicity, compactness, pressure, covering_factor)  # SED

                # Set center wavelength of the filter as normalization wavelength (keeps label)
                self.ski.set_stellar_component_normalization_wavelength(title, fltr.center)

                # Set spectral luminosity at that wavelength (keeps label)
                self.ski.set_stellar_component_luminosity(title, luminosity)

                # Scale height doesn't need to be set as parameter, this is already in the deprojection model

            # Existing component, no MAPPINGS
            else:

                # Get title
                title = component.parameters.title

                # Get SED properties
                template = component.parameters.template
                age = component.parameters.age
                metallicity = component.parameters.metallicity

                # Get normalization
                fltr = parse_filter(component.parameters.filter)
                luminosity = component.parameters.luminosity

                # Set SED
                self.ski.set_stellar_component_sed(title, template, age, metallicity)

                # Set center wavelength of the filter as normalization wavelength (keeps label)
                self.ski.set_stellar_component_normalization_wavelength(title, fltr.center)

                # Set spectral luminosity at that wavelength (keeps label)
                self.ski.set_stellar_component_luminosity(title, luminosity)

                # Scale height doesn't need to be set as parameter, this is already in the deprojection model

    # -----------------------------------------------------------------

    def set_dust_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust components ...")

        # Loop over the dust components
        for name in get_dust_component_names(self.config.path, self.model_name):

            # Load the component
            component = load_dust_component(self.config.path, self.model_name, name)

            # If an input map is required
            if "map_path" in component:

                # Generate a filename for the map
                filename = "dust_" + name + ".fits"

                # Set the filename
                if "deprojection" in component: component.deprojection.filename = filename
                elif "geometry" in component.parameters: component.properties["geometry"].filename = filename
                else: raise RuntimeError("Dust component based on an input map should either have a deprojection or geometry properties")

                # Add entry to the input maps dictionary
                self.input_map_paths[filename] = component.map_path

            # Set geometry
            if "model" in component:

                # Get title
                title = component.parameters.title

                # Set the geometry
                self.ski.set_dust_component_geometry(title, component.model)

            # Set deprojection
            elif "deprojection" in component:

                # Get title
                title = component.parameters.title

                # Set the deprojection geometry
                self.ski.set_dust_component_geometry(title, component.deprojection)

            # Check if this is a new dust component, add geometry, mix and normalization all at once
            if "geometry" in component.parameters:

                # Get title
                title = component.parameters.title

                # Get class names
                geometry_type = component.parameters.geometry
                mix_type = component.parameters.sed
                normalization_type = component.parameters.normalization

                # Get properties for each of the three classes
                geometry_properties = component.properties["geometry"]
                mix_properties = component.properties["mix"]
                normalization_properties = component.properties["normalization"]

                # Create stellar component
                self.ski.create_new_dust_component(title, geometry_type, geometry_properties, mix_type, mix_properties, normalization_type, normalization_properties)

            # Existing component, THEMIS dust mix
            elif "hydrocarbon_pops" in component.parameters:

                # Get title
                title = component.parameters.title

                # Get parameters
                mass = component.parameters.mass
                hydrocarbon_pops = component.parameters.hydrocarbon_pops
                enstatite_pops = component.parameters.enstatite_pops
                forsterite_pops = component.parameters.forsterite_pops

                # Set the dust mix
                self.ski.set_dust_component_themis_mix(title, hydrocarbon_pops, enstatite_pops, forsterite_pops)  # dust mix

                # Set the dust mass (keeps label)
                self.ski.set_dust_component_mass(title, mass)

            # Existing component, not THEMIS dust mix
            else: raise NotImplementedError("Only THEMIS dust mixes are implemented at this moment")

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
                              min_wavelength=self.config.wg.min_wavelength, max_wavelength=self.config.wg.max_wavelength,
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

        # Set the dust emissivity
        self.set_dust_emissivity()

        # Set the lowest-resolution dust grid
        self.ski.set_dust_grid(self.representation.dust_grid)

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
