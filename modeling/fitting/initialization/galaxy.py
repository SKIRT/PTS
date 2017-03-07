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

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from ..component import FittingComponent
from ....core.tools import tables
from ....core.tools import filesystem as fs
from ...basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument
from ....core.data.sun import Sun
from ....magic.tools import wavelengths
from ....core.tools.logging import log
from ....core.prep.wavelengthgrids import WavelengthGridGenerator
from ...component.galaxy import GalaxyModelingComponent
from ..tables import WeightsTable
from ...build.component import get_stellar_component_names, get_dust_component_names, load_stellar_component, load_dust_component
from ....core.filter.filter import parse_filter

# -----------------------------------------------------------------

class GalaxyFittingInitializer(FittingComponent, GalaxyModelingComponent):
    
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
        FittingComponent.__init__(self, config, interactive)
        GalaxyModelingComponent.__init__(self, config, interactive)

        # The fitting run
        self.fitting_run = None

        # The ski file
        self.ski = None

        # The table of weights for each band
        self.weights = None

        # Solar luminosity units
        self.sun_fuv = None
        self.sun_i1 = None

        # The wavelength grid generator
        self.wg_generator = None

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

        # 3. Set the stellar and dust components
        self.set_components()

        # 3. Create the wavelength grids
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
        FittingComponent.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

        # Load the fitting run
        self.fitting_run = self.load_fitting_run(self.config.name)

        # Solar properties
        sun = Sun()
        self.sun_fuv = sun.luminosity_for_filter_as_unit(self.fuv_filter) # Get the luminosity of the Sun in the FUV band
        self.sun_i1 = sun.luminosity_for_filter_as_unit(self.i1_filter)   # Get the luminosity of the Sun in the IRAC I1 band

        # Create a WavelengthGridGenerator
        self.wg_generator = WavelengthGridGenerator()

        # Create the table to contain the weights
        self.weights = WeightsTable()

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.model_name

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        self.ski = self.fitting_run.ski_template

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

                # Add to the dictionary of deprojections
                #self.deprojections[(name, title)] = component.deprojection

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

                # Add to the dictionary of deprojections
                #self.deprojections[(name, title)] = component.deprojection

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
        fixed = [self.i1_filter.pivotwavelength(), self.fuv_filter.pivotwavelength()]

        # Set options
        self.wg_generator.config.show = False
        self.wg_generator.config.write = False

        # Generate the wavelength grids
        self.wg_generator.run(npoints_range=self.config.wg.npoints_range, ngrids=self.config.wg.ngrids,
                              fixed=fixed, add_emission_lines=self.config.wg.add_emission_lines,
                              min_wavelength=self.config.wg.min_wavelength, max_wavelength=self.config.wg.max_wavelength,
                              filters=self.fitting_run.fitting_filters)

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Create an SED instrument
        self.instruments["SED"] = SEDInstrument.from_projection(self.earth_projection)

        # Create a frame instrument to generate datacube
        self.instruments["frame"] = FrameInstrument.from_projection(self.earth_projection)

        # Create a simple instrument (SED + frame)
        self.instruments["simple"] = SimpleInstrument.from_projection(self.earth_projection)

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
        self.ski.add_instrument("earth", self.instruments["SED"])

        # Set the number of photon packages
        self.ski.setpackages(self.config.npackages)

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid("wavelengths.txt")

        # Set the stellar and dust components
        #self.set_components()

        # Set the dust emissivity
        self.set_dust_emissivity()

        # Set the lowest-resolution dust grid
        self.ski.set_dust_grid(self.dg_generator.grids[0])

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

    def calculate_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the weight to give to each band ...")

        # Initialize lists to contain the filters of the different wavelength ranges
        uv_bands = []
        optical_bands = []
        nir_bands = []
        mir_bands = []
        fir_bands = []
        submm_bands = []

        # Set the number of groups
        number_of_groups = 6

        # Loop over the observed SED filters
        for fltr in self.fitting_run.fitting_filters:

            # Get the central wavelength
            wavelength = fltr.center

            # Get a string identifying which portion of the wavelength spectrum this wavelength belongs to
            spectrum = wavelengths.name_in_spectrum(wavelength)

            # Determine to which group
            if spectrum[0] == "UV": uv_bands.append(fltr)
            elif spectrum[0] == "Optical": optical_bands.append(fltr)
            elif spectrum[0] == "Optical/IR": optical_bands.append(fltr)
            elif spectrum[0] == "IR":
                if spectrum[1] == "NIR": nir_bands.append(fltr)
                elif spectrum[1] == "MIR": mir_bands.append(fltr)
                elif spectrum[1] == "FIR": fir_bands.append(fltr)
                else: raise RuntimeError("Unknown IR range")
            elif spectrum[0] == "Submm": submm_bands.append(fltr)
            else: raise RuntimeError("Unknown wavelength range")

        # Determine the weight for each group of filters
        number_of_data_points = len(self.fitting_run.fitting_filters)
        uv_weight = 1. / (len(uv_bands) * number_of_groups) * number_of_data_points
        optical_weight = 1. / (len(optical_bands) * number_of_groups) * number_of_data_points
        nir_weight = 1. / (len(nir_bands) * number_of_groups) * number_of_data_points
        mir_weight = 1. / (len(mir_bands) * number_of_groups) * number_of_data_points
        fir_weight = 1. / (len(fir_bands) * number_of_groups) * number_of_data_points
        submm_weight = 1. / (len(submm_bands) * number_of_groups) * number_of_data_points

        # Debugging
        log.debug("UV: number of bands = " + str(len(uv_bands)) + ", weight = " + str(uv_weight))
        log.debug("Optical: number of bands = " + str(len(optical_bands)) + ", weight = " + str(optical_weight))
        log.debug("NIR: number of bands = " + str(len(nir_bands)) + ", weight = " + str(nir_weight))
        log.debug("MIR: number of bands = " + str(len(mir_bands)) + ", weight = " + str(mir_weight))
        log.debug("FIR: number of bands = " + str(len(fir_bands)) + ", weight = " + str(fir_weight))
        log.debug("Submm: number of bands = " + str(len(submm_bands)) + ", weight = " + str(submm_weight))

        # Loop over the bands in each group and set the weight in the weights table
        for fltr in uv_bands: self.weights.add_row([fltr.instrument, fltr.band, uv_weight])
        for fltr in optical_bands: self.weights.add_row([fltr.instrument, fltr.band, optical_weight])
        for fltr in nir_bands: self.weights.add_row([fltr.instrument, fltr.band, nir_weight])
        for fltr in mir_bands: self.weights.add_row([fltr.instrument, fltr.band, mir_weight])
        for fltr in fir_bands: self.weights.add_row([fltr.instrument, fltr.band, fir_weight])
        for fltr in submm_bands: self.weights.add_row([fltr.instrument, fltr.band, submm_weight])

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 2. Write the ski file
        self.write_ski()

        # 4. Write the weights table
        self.write_weights()

        # 6. Write the wavelength grids
        self.write_wavelength_grids()

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

    def write_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with weights to " + self.fitting_run.weights_table_path + " ...")

        # Write the table with weights
        self.weights.saveto(self.fitting_run.weights_table_path)

    # -----------------------------------------------------------------

    def write_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grids ...")

        # Loop over the grids
        index = 0
        for grid in self.wg_generator.grids:

            # Determine the path to the grid
            path = fs.join(self.fitting_run.wavelength_grids_path, str(index) + ".txt")

            # Save the wavelength grid
            grid.to_skirt_input(path)

            # Increment the index
            index += 1

        # Write the wavelength grids table
        tables.write(self.wg_generator.table, self.fitting_run.wavelength_grids_table_path)

# -----------------------------------------------------------------
