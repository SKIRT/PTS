#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.initialization Contains the FittingInitializer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import astronomical modules
from astropy.units import Unit, dimensionless_angles
from astropy import constants
from astropy.table import Table

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import tables
from ...core.tools import filesystem as fs
from ...core.basics.filter import Filter
from ..basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument
from ..core.sun import Sun
from ..core.mappings import Mappings
from ...magic.tools import wavelengths
from ...core.tools.logging import log
from .wavelengthgrids import WavelengthGridGenerator
from .dustgrids import DustGridGenerator
from ...core.basics.range import IntegerRange, RealRange, QuantityRange
from ..basics.models import DeprojectionModel3D

# -----------------------------------------------------------------

class FittingInitializer(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(FittingInitializer, self).__init__(config)

        # -- Attributes --

        # The deprojection model
        self.deprojection = None
        self.deprojections = dict()

        # The instruments
        self.instruments = dict()

        # The table of weights for each band
        self.weights = None

        # Solar luminosity units
        self.sun_fuv = None
        self.sun_i1 = None

        # The ski files for simulating the contributions of the various stellar components
        self.ski_contributions = dict()

        # The ski file for generating simulated images
        self.ski_images = None

        # The wavelength grid and dust grid generators
        self.wg_generator = None
        self.dg_generator = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Create the wavelength grids
        self.create_wavelength_grids()

        # 3. Create the deprojection model
        self.create_deprojection_model()

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

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingInitializer, self).setup()

        # Solar properties
        sun = Sun()
        self.sun_fuv = sun.luminosity_for_filter_as_unit(self.fuv_filter) # Get the luminosity of the Sun in the FUV band
        self.sun_i1 = sun.luminosity_for_filter_as_unit(self.i1_filter)   # Get the luminosity of the Sun in the IRAC I1 band

        # Create a WavelengthGridGenerator
        self.wg_generator = WavelengthGridGenerator()

        # Create the DustGridGenerator
        self.dg_generator = DustGridGenerator()

        # Create the table to contain the weights
        self.weights = Table(names=["Instrument", "Band", "Weight"], dtype=["S5", "S7", "float64"])

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

        # Generate the wavelength grids
        self.wg_generator.run(self.config.wg.npoints_range, self.config.wg.ngrids, fixed=fixed)

    # -----------------------------------------------------------------

    def create_deprojection_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the deprojection parameters ...")

        filename = None
        hz = None

        # Get the galaxy distance, the inclination and position angle
        distance = self.galaxy_properties.distance
        inclination = self.galaxy_properties.inclination
        pa = self.earth_projection.position_angle

        # Get the center pixel
        pixel_center = self.galaxy_properties.center.to_pixel(self.reference_wcs)
        xc = pixel_center.x
        yc = pixel_center.y

        # Get the pixelscale in physical units
        pixelscale_angular = self.reference_wcs.average_pixelscale * Unit("pix") # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # Get the number of x and y pixels
        x_size = self.reference_wcs.xsize
        y_size = self.reference_wcs.ysize

        # Create the deprojection model
        self.deprojection = DeprojectionModel3D(filename, pixelscale, pa, inclination, x_size, y_size, xc, yc, hz)

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
        major_angular = self.truncation_ellipse.major  # major axis length of the sky ellipse
        radius_physical = (major_angular * self.galaxy_properties.distance).to("pc", equivalencies=dimensionless_angles())

        # Get the pixelscale in physical units
        distance = self.galaxy_properties.distance
        pixelscale_angular = self.reference_wcs.average_pixelscale * Unit("pix")  # in deg
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
        self.dg_generator.z_radius = 3. * Unit("kpc")

        # Generate the dust grids
        self.dg_generator.run(scale_range, self.config.dg.level_range, mass_fraction_range, 10)

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # Remove the existing instruments
        self.ski_template.remove_all_instruments()

        # Add the instrument
        self.ski_template.add_instrument("earth", self.instruments["SED"])

        # Set the number of photon packages
        self.ski_template.setpackages(self.config.packages)

        # Set the name of the wavelength grid file
        self.ski_template.set_file_wavelength_grid("wavelengths.txt")

        # Set the stellar and dust components
        self.set_components()

        # Set transient dust emissivity
        self.ski_template.set_transient_dust_emissivity()

        # Set the lowest-resolution dust grid
        self.ski_template.set_dust_grid(self.dg_generator.grids[0])

        # Set all-cells dust library
        self.ski_template.set_allcells_dust_lib()

        # Dust self-absorption
        if self.config.selfabsorption: self.ski_template.enable_selfabsorption()
        else: self.ski_template.disable_selfabsorption()

        # Disable all writing options
        self.ski_template.disable_all_writing_options()

    # -----------------------------------------------------------------

    def set_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the stellar and dust components ...")

        # Set the evolved stellar bulge component
        self.set_bulge_component()

        # Set the evolved stellar disk component
        self.set_old_stellar_component()

        # Set the young stellar component
        self.set_young_stellar_component()

        # Set the ionizing stellar component
        self.set_ionizing_stellar_component()

        # The dust component
        self.set_dust_component()

    # -----------------------------------------------------------------

    def set_bulge_component(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the bulge component ...")

        # Like M31
        bulge_template = "BruzualCharlot"
        bulge_age = 12
        #bulge_metallicity = 0.02
        bulge_metallicity = 0.03

        # Get the flux density of the bulge
        fluxdensity = self.bulge2d_model.fluxdensity

        # Convert the flux density into a spectral luminosity
        luminosity = fluxdensity_to_luminosity(fluxdensity, self.i1_filter.pivotwavelength() * Unit("micron"), self.galaxy_properties.distance)

        # Get the spectral luminosity in solar units
        #luminosity = luminosity.to(self.sun_i1).value

        # Set the parameters of the bulge
        self.ski_template.set_stellar_component_geometry("Evolved stellar bulge", self.bulge_model)
        self.ski_template.set_stellar_component_sed("Evolved stellar bulge", bulge_template, bulge_age, bulge_metallicity) # SED
        #self.ski.set_stellar_component_luminosity("Evolved stellar bulge", luminosity, self.i1) # normalization by band
        self.ski_template.set_stellar_component_luminosity("Evolved stellar bulge", luminosity, self.i1_filter.centerwavelength() * Unit("micron"))

    # -----------------------------------------------------------------

    def set_old_stellar_component(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the old stellar component ...")

        # Like M31
        disk_template = "BruzualCharlot"
        disk_age = 8
        #disk_metallicity = 0.02
        disk_metallicity = 0.03

        # Get the scale height
        #scale_height = 521. * Unit("pc") # first models
        scale_height = self.disk2d_model.scalelength / 8.26 # De Geyter et al. 2014
        bulge_fluxdensity = self.bulge2d_model.fluxdensity

        # Get the 3.6 micron flux density with the bulge subtracted
        fluxdensity = self.observed_flux(self.i1_filter, unit="Jy") - bulge_fluxdensity

        # Convert the flux density into a spectral luminosity
        luminosity = fluxdensity_to_luminosity(fluxdensity, self.i1_filter.pivotwavelength() * Unit("micron"), self.galaxy_properties.distance)

        # Get the spectral luminosity in solar units
        #luminosity = luminosity.to(self.sun_i1).value

        # Set the parameters of the evolved stellar component
        deprojection = self.deprojection.copy()
        deprojection.filename = "old_stars.fits"
        deprojection.scale_height = scale_height
        self.deprojections["Old stars"] = deprojection

        # Adjust the ski file
        self.ski_template.set_stellar_component_geometry("Evolved stellar disk", deprojection)
        self.ski_template.set_stellar_component_sed("Evolved stellar disk", disk_template, disk_age, disk_metallicity) # SED
        #self.ski.set_stellar_component_luminosity("Evolved stellar disk", luminosity, self.i1) # normalization by band
        self.ski_template.set_stellar_component_luminosity("Evolved stellar disk", luminosity, self.i1_filter.centerwavelength() * Unit("micron"))

    # -----------------------------------------------------------------

    def set_young_stellar_component(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the young stellar component ...")

        # Like M31
        young_template = "BruzualCharlot"
        young_age = 0.1
        #young_metallicity = 0.02
        young_metallicity = 0.03

        # Get the scale height
        #scale_height = 150 * Unit("pc") # first models
        scale_height = 100. * Unit("pc") # M51

        # Get the FUV flux density
        fluxdensity = 2. * self.observed_flux(self.fuv_filter, unit="Jy")

        # Convert the flux density into a spectral luminosity
        luminosity = fluxdensity_to_luminosity(fluxdensity, self.fuv_filter.pivotwavelength() * Unit("micron"), self.galaxy_properties.distance)

        # Get the spectral luminosity in solar units
        #luminosity = luminosity.to(self.sun_fuv).value # for normalization by band

        # Set the parameters of the young stellar component
        deprojection = self.deprojection.copy()
        deprojection.filename = "young_stars.fits"
        deprojection.scale_height = scale_height
        self.deprojections["Young stars"] = deprojection

        # Adjust the ski file
        self.ski_template.set_stellar_component_geometry("Young stars", deprojection)
        self.ski_template.set_stellar_component_sed("Young stars", young_template, young_age, young_metallicity) # SED
        #self.ski.set_stellar_component_luminosity("Young stars", luminosity, self.fuv) # normalization by band
        #self.ski_template.set_stellar_component_luminosity("Young stars", luminosity, self.fuv_filter.centerwavelength() * Unit("micron"))

        # SET NORMALIZATION (IS FREE PARAMETER)
        self.ski_template.set_stellar_component_normalization_wavelength("Young stars", self.fuv_filter.centerwavelength() * Unit("micron"))
        self.ski_template.set_labeled_value("fuv_young", luminosity)

    # -----------------------------------------------------------------

    def set_ionizing_stellar_component(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the ionizing stellar component ...")

        # Like M51 and M31
        #ionizing_metallicity = 0.02
        ionizing_metallicity = 0.03 # XU KONG et al. 2000
        ionizing_compactness = 6
        ionizing_pressure = 1e12 * Unit("K/m3")
        ionizing_covering_factor = 0.2

        # Get the scale height
        #scale_height = 150 * Unit("pc") # first models
        scale_height = 100. * Unit("pc") # M51

        # Convert the SFR into a FUV luminosity
        sfr = 0.8 # The star formation rate # see Perez-Gonzalez 2006 (mentions Devereux et al 1995)
        mappings = Mappings(ionizing_metallicity, ionizing_compactness, ionizing_pressure, ionizing_covering_factor, sfr)
        luminosity = mappings.luminosity_for_filter(self.fuv_filter)
        #luminosity = luminosity.to(self.sun_fuv).value # for normalization by band

        # Set the parameters of the ionizing stellar component
        deprojection = self.deprojection.copy()
        deprojection.filename = "ionizing_stars.fits"
        deprojection.scale_height = scale_height
        self.deprojections["Ionizing stars"] = deprojection

        # Adjust the ski file
        self.ski_template.set_stellar_component_geometry("Ionizing stars", deprojection)
        self.ski_template.set_stellar_component_mappingssed("Ionizing stars", ionizing_metallicity, ionizing_compactness, ionizing_pressure, ionizing_covering_factor) # SED
        #self.ski.set_stellar_component_luminosity("Ionizing stars", luminosity, self.fuv) # normalization by band
        #self.ski_template.set_stellar_component_luminosity("Ionizing stars", luminosity, self.fuv_filter.centerwavelength() * Unit("micron"))

        # SET NORMALIZATION (IS FREE PARAMETER)
        self.ski_template.set_stellar_component_normalization_wavelength("Ionizing stars", self.fuv_filter.centerwavelength() * Unit("micron"))
        self.ski_template.set_labeled_value("fuv_ionizing", luminosity) # keep label

    # -----------------------------------------------------------------

    def set_dust_component(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the dust component ...")

        #scale_height = 260.5 * Unit("pc") # first models
        scale_height = 200. * Unit("pc") # M51
        dust_mass = 1.5e7 * Unit("Msun")

        hydrocarbon_pops = 25
        enstatite_pops = 25
        forsterite_pops = 25

        # Set the parameters of the dust component
        deprojection = self.deprojection.copy()
        deprojection.filename = "dust.fits"
        deprojection.scale_height = scale_height
        self.deprojections["Dust"] = deprojection

        # Adjust the ski file
        self.ski_template.set_dust_component_geometry(0, deprojection)
        self.ski_template.set_dust_component_themis_mix(0, hydrocarbon_pops, enstatite_pops, forsterite_pops) # dust mix
        #self.ski_template.set_dust_component_mass(0, dust_mass) # dust mass

        self.ski_template.set_labeled_value("dust_mass", dust_mass) # keep label

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
        for fltr in self.fitting_filters:

            # Get the central wavelength
            wavelength = fltr.centerwavelength() * Unit("micron")

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
        number_of_data_points = len(self.fitting_filters)
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

        # 3. Write the weights table
        self.write_weights()

        # 4. Write the deprojection models
        self.write_deprojection_models()

        # 5. Write the wavelength grids
        self.write_wavelength_grids()

        # Write the dust grids
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
        self.instruments["SED"].save(self.sed_instrument_path)

        # Write the frame instrument
        self.instruments["frame"].save(self.frame_instrument_path)

        # Write the simple instrument
        self.instruments["simple"].save(self.simple_instrument_path)

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file to " + self.template_ski_path + " ...")

        # Save the ski template file
        self.ski_template.save()

    # -----------------------------------------------------------------

    def write_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with weights to " + self.weights_table_path + " ...")

        # Write the table with weights
        tables.write(self.weights, self.weights_table_path)

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
            self.deprojections[label].save(path)

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
            path = fs.join(self.fit_wavelength_grids_path, str(index) + ".txt")

            # Save the wavelength grid
            grid.to_skirt_input(path)

            # Increment the index
            index += 1

        # Write the wavelength grids table
        tables.write(self.wg_generator.table, self.wavelength_grids_table_path)

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
            path = fs.join(self.fit_dust_grids_path, str(index) + ".dg")

            # Save the dust grid
            grid.save(path)

            # Increment the index
            index += 1

        # Write the dust grids table
        tables.write(self.dg_generator.table, self.dust_grids_table_path)

# -----------------------------------------------------------------

# The speed of light
speed_of_light = constants.c

# -----------------------------------------------------------------

def spectral_factor_hz_to_micron(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    wavelength_unit = "micron"
    frequency_unit = "Hz"

    # Convert string units to Unit objects
    if isinstance(wavelength_unit, basestring): wavelength_unit = Unit(wavelength_unit)
    if isinstance(frequency_unit, basestring): frequency_unit = Unit(frequency_unit)

    conversion_factor_unit = wavelength_unit / frequency_unit

    # Calculate the conversion factor
    factor = (wavelength ** 2 / speed_of_light).to(conversion_factor_unit).value
    return 1. / factor

# -----------------------------------------------------------------

def fluxdensity_to_luminosity(fluxdensity, wavelength, distance):

    """
    This function ...
    :param fluxdensity:
    :param wavelength:
    :param distance:
    :return:
    """

    luminosity = (fluxdensity * 4. * math.pi * distance ** 2.).to("W/Hz")

    # 3 ways:
    #luminosity_ = luminosity.to("W/micron", equivalencies=spectral_density(wavelength)) # does not work
    luminosity_ = (speed_of_light * luminosity / wavelength**2).to("W/micron")
    luminosity = luminosity.to("W/Hz").value * spectral_factor_hz_to_micron(wavelength) * Unit("W/micron")
    #print(luminosity_, luminosity) # is OK!

    return luminosity

# -----------------------------------------------------------------
