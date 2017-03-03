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
from ....core.prep.dustgrids import DustGridGenerator
from ....core.basics.range import RealRange, QuantityRange
from ....core.basics.configuration import write_mapping
from ....core.basics.map import Map
from ...component.galaxy import GalaxyModelingComponent
from ..tables import WeightsTable
from ....core.basics.unit import parse_unit as u

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

        # The instruments
        self.instruments = dict()

        # The table of weights for each band
        self.weights = None

        # Solar luminosity units
        self.sun_fuv = None
        self.sun_i1 = None

        # The wavelength grid and dust grid generators
        self.wg_generator = None
        self.dg_generator = None

        # The fixed parameters map
        self.fixed = Map()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Load the ski file
        self.load_ski()

        # 2. Create the wavelength grids
        self.create_wavelength_grids()

        # 4. Create the instruments
        self.create_instruments()

        # 5. Create the dust grids
        self.create_dust_grids()

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

        # Create the DustGridGenerator
        self.dg_generator = DustGridGenerator()

        # Create the table to contain the weights
        self.weights = WeightsTable()

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        self.ski = self.fitting_run.ski_template

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

    def create_dust_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the grids ...")

        # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
        semimajor_angular = self.truncation_ellipse.semimajor  # semimajor axis length of the sky ellipse
        radius_physical = (semimajor_angular * self.galaxy_properties.distance).to("pc", equivalencies=dimensionless_angles())

        # Get the pixelscale in physical units
        distance = self.galaxy_properties.distance
        pixelscale_angular = self.reference_wcs.average_pixelscale.to("deg")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # BINTREE: (smallest_cell_pixels, min_level, max_mass_fraction)
        # Low-resolution: 10., 6, 1e-5
        # High-resolution: 0.5, 9, 0.5e-6

        # OCTTREE:
        # Low-resolution: 10., 2, 1e-5
        # High-resolution: 0.5, 3, 0.5e-6

        # Because we (currently) can't position the grid exactly as the 2D pixels (rotation etc.),
        # take half of the pixel size to avoid too much interpolation
        min_scale = self.config.dg.scale_range.min * pixelscale
        max_scale = self.config.dg.scale_range.max * pixelscale
        scale_range = QuantityRange(min_scale, max_scale, invert=True)

        # The range of the max mass fraction
        mass_fraction_range = RealRange(self.config.dg.mass_fraction_range.min, self.config.dg.mass_fraction_range.max, invert=True) # must be inverted

        # Set fixed grid properties
        self.dg_generator.grid_type = self.config.dg.grid_type # set grid type
        self.dg_generator.x_radius = radius_physical
        self.dg_generator.y_radius = radius_physical
        self.dg_generator.z_radius = 3. * u("kpc")

        # Set options
        self.dg_generator.show = False
        self.dg_generator.write = False

        # Generate the dust grids
        self.dg_generator.run(scale_range=scale_range, level_range=self.config.dg.level_range, mass_fraction_range=mass_fraction_range, ngrids=10)

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
        self.set_components()

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

    def set_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the stellar and dust components ...")

        # Set stellar components
        self.set_stellar_components()

        # Set dust components
        self.set_dust_components()

    # -----------------------------------------------------------------

    def set_stellar_components(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def set_dust_components(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def set_bulge_component(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the bulge component ...")

        # Set the parameters of the bulge
        self.ski.set_stellar_component_geometry("Evolved stellar bulge", self.bulge_model)
        self.ski.set_stellar_component_sed("Evolved stellar bulge", bulge_template, bulge_age, bulge_metallicity) # SED
        #self.ski.set_stellar_component_luminosity("Evolved stellar bulge", luminosity, self.i1) # normalization by band
        self.ski.set_stellar_component_luminosity("Evolved stellar bulge", luminosity, self.i1_filter.centerwavelength() * u("micron"))

    # -----------------------------------------------------------------

    def set_old_stellar_component(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the old stellar component ...")

        ## NEW: SET FIXED PARAMETERS
        self.fixed["metallicity"] = disk_metallicity
        self.fixed["old_scaleheight"] = scale_height
        self.fixed["i1_old"] = luminosity
        ##

        # Set the parameters of the evolved stellar component
        #deprojection = self.deprojection.copy()
        #deprojection.filename = self.old_stellar_map_filename
        #deprojection.scale_height = scale_height
        #self.deprojections["old stars"] = deprojection

        # Adjust the ski file
        self.ski.set_stellar_component_geometry("Evolved stellar disk", deprojection)
        self.ski.set_stellar_component_sed("Evolved stellar disk", disk_template, disk_age, disk_metallicity) # SED
        #self.ski.set_stellar_component_luminosity("Evolved stellar disk", luminosity, self.i1) # normalization by band
        self.ski.set_stellar_component_luminosity("Evolved stellar disk", luminosity, self.i1_filter.centerwavelength() * u("micron"))

    # -----------------------------------------------------------------

    def set_young_stellar_component(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the young stellar component ...")

        ## NEW: SET FIXED PARAMETERS
        self.fixed["young_scaleheight"] = scale_height
        ##

        # Set the parameters of the young stellar component
        #deprojection = self.deprojection.copy()
        #deprojection.filename = self.young_stellar_map_filename
        #deprojection.scale_height = scale_height
        #self.deprojections["young stars"] = deprojection

        # Adjust the ski file
        self.ski.set_stellar_component_geometry("Young stars", deprojection)
        self.ski.set_stellar_component_sed("Young stars", young_template, young_age, young_metallicity) # SED
        #self.ski.set_stellar_component_luminosity("Young stars", luminosity, self.fuv) # normalization by band
        #self.ski_template.set_stellar_component_luminosity("Young stars", luminosity, self.fuv_filter.centerwavelength() * u("micron"))

        # SET NORMALIZATION (IS FREE PARAMETER)
        self.ski.set_stellar_component_normalization_wavelength("Young stars", self.fuv_filter.centerwavelength() * u("micron"))
        self.ski.set_labeled_value("fuv_young", luminosity)

    # -----------------------------------------------------------------

    def set_ionizing_stellar_component(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the ionizing stellar component ...")

        ## NEW: SET FIXED PARAMETERS
        self.fixed["ionizing_scaleheight"] = scale_height
        self.fixed["sfr_compactness"] = ionizing_compactness
        self.fixed["sfr_covering"] = ionizing_covering_factor
        self.fixed["sfr_pressure"] = ionizing_pressure
        ##

        # Set the parameters of the ionizing stellar component
        #deprojection = self.deprojection.copy()
        #deprojection.filename = self.ionizing_stellar_map_filename
        #deprojection.scale_height = scale_height
        #self.deprojections["ionizing stars"] = deprojection

        # Adjust the ski file
        self.ski.set_stellar_component_geometry("Ionizing stars", deprojection)
        self.ski.set_stellar_component_mappingssed("Ionizing stars", ionizing_metallicity, ionizing_compactness, ionizing_pressure, ionizing_covering_factor) # SED
        #self.ski.set_stellar_component_luminosity("Ionizing stars", luminosity, self.fuv) # normalization by band
        #self.ski_template.set_stellar_component_luminosity("Ionizing stars", luminosity, self.fuv_filter.centerwavelength() * Unit("micron"))

        # SET NORMALIZATION (IS FREE PARAMETER)
        self.ski.set_stellar_component_normalization_wavelength("Ionizing stars", self.fuv_filter.centerwavelength() * u("micron"))
        self.ski.set_labeled_value("fuv_ionizing", luminosity) # keep label

    # -----------------------------------------------------------------

    def set_dust_component(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the dust component ...")

        ## NEW: SET FIXED PARAMETERS
        self.fixed["dust_scaleheight"] = scale_height
        ##

        # Set the parameters of the dust component
        #deprojection = self.deprojection.copy()
        #deprojection.filename = self.dust_map_filename
        #deprojection.scale_height = scale_height
        #self.deprojections["dust"] = deprojection

        # Adjust the ski file
        self.ski.set_dust_component_geometry(0, deprojection)
        self.ski.set_dust_component_themis_mix(0, hydrocarbon_pops, enstatite_pops, forsterite_pops) # dust mix
        #self.ski_template.set_dust_component_mass(0, dust_mass) # dust mass

        self.ski.set_labeled_value("dust_mass", dust_mass) # keep label

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

        # 1. Write the instruments
        self.write_instruments()

        # 2. Write the ski file
        self.write_ski()

        # 3. Write the fixed parameters
        self.write_fixed()

        # 4. Write the weights table
        self.write_weights()

        # 5. Write the deprojection models
        self.write_deprojection_models()

        # 6. Write the wavelength grids
        self.write_wavelength_grids()

        # 7. Write the dust grids
        self.write_dust_grids()

    # -----------------------------------------------------------------

    def write_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SED, frame and simple instruments ...")

        # Write the SED instrument
        self.instruments["SED"].saveto(self.sed_instrument_path)

        # Write the frame instrument
        self.instruments["frame"].saveto(self.frame_instrument_path)

        # Write the simple instrument
        self.instruments["simple"].saveto(self.simple_instrument_path)

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

    def write_fixed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the fixed parameter values to " + self.fitting_run.fixed_parameters_path + " ...")

        # Write the fixed parameters map
        with open(self.fitting_run.fixed_parameters_path, 'w') as f: write_mapping(f, self.fixed)

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

    def write_deprojection_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojection model for the other components ...")

        # Write the deprojection models
        for label in self.deprojections:

            # Save the deprojection model
            path = fs.join(self.fit_geometries_path, label + ".mod")
            self.deprojections[label].saveto(path)

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
            path = fs.join(self.fitting_run.fit_wavelength_grids_path, str(index) + ".txt")

            # Save the wavelength grid
            grid.to_skirt_input(path)

            # Increment the index
            index += 1

        # Write the wavelength grids table
        tables.write(self.wg_generator.table, self.fitting_run.wavelength_grids_table_path)

    # -----------------------------------------------------------------

    def write_dust_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grids ...")

        # Loop over the grids
        index = 0
        for grid in self.dg_generator.grids:

            # Determine the path to the grid
            path = fs.join(self.fitting_run.fit_dust_grids_path, str(index) + ".dg")

            # Save the dust grid
            grid.saveto(path)

            # Increment the index
            index += 1

        # Write the dust grids table
        tables.write(self.dg_generator.table, self.fitting_run.dust_grids_table_path)

# -----------------------------------------------------------------
