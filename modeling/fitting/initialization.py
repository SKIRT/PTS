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

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import inspection, tables
from ...core.tools import filesystem as fs
from ...core.simulation.skifile import SkiFile
from ...core.basics.filter import Filter
from ..basics.models import SersicModel, DeprojectionModel
from ...magic.basics.coordinatesystem import CoordinateSystem
from ..decomposition.decomposition import load_parameters
from ...magic.basics.skyregion import SkyRegion
from ..basics.instruments import SEDInstrument, FrameInstrument
from ..core.sun import Sun
from ..core.mappings import Mappings
from ...magic.tools import wavelengths
from ...core.tools.logging import log
from ..basics.projection import GalaxyProjection
from ...core.simulation.execute import SkirtExec
from ...core.simulation.arguments import SkirtArguments
from ..core.sed import ObservedSED
from .wavelengthgrids import WavelengthGridGenerator
from .dustgrids import DustGridGenerator
from ...core.basics.range import IntegerRange, RealRange, QuantityRange

# -----------------------------------------------------------------

template_ski_path = fs.join(inspection.pts_dat_dir("modeling"), "ski", "template.ski")

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

        # The ski file
        self.ski = None

        # The structural parameters
        self.parameters = None

        # The projection system
        self.projection = None

        # The truncation ellipse
        self.ellipse = None

        # The geometric bulge model
        self.bulge = None

        # The deprojection model
        self.deprojection = None
        self.deprojections = dict()

        # The instrument
        self.instrument = None

        # The table of weights for each band
        self.weights = None

        # The observed SED
        self.observed_sed = None

        # Filters
        self.i1 = None
        self.fuv = None

        # Solar luminosity units
        self.sun_fuv = None
        self.sun_i1 = None

        # Coordinate system
        self.reference_wcs = None

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

        # 2. Load the template ski file
        self.load_template()

        # 3. Load the structural parameters of the galaxy
        self.load_parameters()

        # 4. Load the projection system
        self.load_projection()

        # 5. Load the truncation ellipse
        self.load_truncation_ellipse()

        # 6. Load the observed SED
        self.load_observed_sed()

        # 7. Create the wavelength grid
        self.create_wavelength_grids()

        # 8. Create the bulge model
        self.create_bulge_model()

        # 9. Create the deprojection model
        self.create_deprojection_model()

        # 10. Create the instrument
        self.create_instrument()

        # 11. Create the dust grids
        self.create_dust_grids()

        # 12. Adjust the ski file
        self.adjust_ski()

        # 13. Adjust the ski files for simulating the contributions of the various stellar components
        self.adjust_ski_contributions()

        # 14. Adjust the ski file for generating simulated images
        self.adjust_ski_images()

        # 15. Generate the grids
        self.write_input()
        self.generate_grids()

        # 16. Calculate the weight factor to give to each band
        self.calculate_weights()

        # 17. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingInitializer, self).setup()

        # Create filters
        self.i1 = Filter.from_string("I1")
        self.fuv = Filter.from_string("FUV")

        # Solar properties
        sun = Sun()
        self.sun_fuv = sun.luminosity_for_filter_as_unit(self.fuv) # Get the luminosity of the Sun in the FUV band
        self.sun_i1 = sun.luminosity_for_filter_as_unit(self.i1)   # Get the luminosity of the Sun in the IRAC I1 band

        # Reference coordinate system
        reference_path = fs.join(self.truncation_path, self.reference_image + ".fits")
        self.reference_wcs = CoordinateSystem.from_file(reference_path)

        # Create a WavelengthGridGenerator
        self.wg_generator = WavelengthGridGenerator()

        # Create the DustGridGenerator
        self.dg_generator = DustGridGenerator()

    # -----------------------------------------------------------------

    def load_template(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file template ...")

        # Open the template ski file
        self.ski = SkiFile(template_ski_path)

    # -----------------------------------------------------------------

    def load_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the decomposition parameters ...")

        # Determine the path to the parameters file
        path = fs.join(self.components_path, "parameters.dat")

        # Load the parameters
        self.parameters = load_parameters(path)

    # -----------------------------------------------------------------

    def load_projection(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the projection system ...")

        # Determine the path to the projection file
        path = fs.join(self.components_path, "earth.proj")

        # Load the projection system
        self.projection = GalaxyProjection.from_file(path)

    # -----------------------------------------------------------------

    def load_truncation_ellipse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ellipse region used for truncating the observed images ...")

        # Determine the path
        path = fs.join(self.truncation_path, "ellipse.reg")

        # Get the ellipse
        region = SkyRegion.from_file(path)
        self.ellipse = region[0]

    # -----------------------------------------------------------------

    def load_observed_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed SED ...")

        # Load the SED
        self.observed_sed = ObservedSED.from_file(self.observed_sed_path)

    # -----------------------------------------------------------------

    def create_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grids ...")

        # Create the range of npoints for the wavelength grids
        npoints_range = IntegerRange(150, 500)

        # Generate the wavelength grids
        self.wg_generator.run(npoints_range, 10)

    # -----------------------------------------------------------------

    def create_bulge_model(self):

        """
        :return:
        """

        # Inform the user
        log.info("Creating the bulge model ...")

        # Create a Sersic model for the bulge
        self.bulge = SersicModel.from_galfit(self.parameters.bulge, self.parameters.inclination, self.parameters.disk.PA)

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
        distance = self.parameters.distance
        inclination = self.parameters.inclination
        pa = self.parameters.disk.PA

        # Get the center pixel
        pixel_center = self.parameters.center.to_pixel(self.reference_wcs)
        xc = pixel_center.x
        yc = pixel_center.y

        # Get the pixelscale in physical units
        pixelscale_angular = self.reference_wcs.xy_average_pixelscale * Unit("pix") # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # Get the number of x and y pixels
        x_size = self.reference_wcs.xsize
        y_size = self.reference_wcs.ysize

        # Create the deprojection model
        self.deprojection = DeprojectionModel(filename, pixelscale, pa, inclination, x_size, y_size, xc, yc, hz)

    # -----------------------------------------------------------------

    def create_instrument(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instrument ...")

        # Create an SED instrument
        self.instrument = SEDInstrument.from_projection(self.projection)

    # -----------------------------------------------------------------

    def create_dust_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the grids ...")

        # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
        major_angular = self.ellipse.major  # major axis length of the sky ellipse
        radius_physical = (major_angular * self.parameters.distance).to("pc", equivalencies=dimensionless_angles())

        # Get the pixelscale in physical units
        distance = self.parameters.distance
        pixelscale_angular = self.reference_wcs.xy_average_pixelscale * Unit("pix")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # BINTREE: (smallest_cell_pixels, min_level, max_mass_fraction)
        # Low-resolution: 10., 6, 1e-5
        # High-resolution: 0.5, 9, 0.5e-6

        # OCTTREE:
        # Low-resolution: 10., 2, 1e-5
        # High-resolution: 0.5, 3, 0.5e-6

        # Because we (currently) can't position the grid exactly as the 2D pixels (rotation etc.),
        # take half of the pixel size to avoid too much interpolation
        min_scale = 0.5 * pixelscale
        max_scale = 10. * pixelscale
        scale_range = QuantityRange(min_scale, max_scale, invert=True)

        # The range of the maximum depth level of the tree
        level_range = IntegerRange(6, 9)

        # The range of the max mass fraction
        mass_fraction_range = RealRange(0.5e-6, 1e-5, invert=True)

        # Set fixed grid properties
        self.dg_generator.grid_type = "bintree"
        self.dg_generator.x_radius = radius_physical
        self.dg_generator.y_radius = radius_physical
        self.dg_generator.z_radius = 3. * Unit("kpc")

        # Generate the dust grids
        self.dg_generator.run(scale_range, level_range, mass_fraction_range, 10, grid_type="bintree")

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
        self.ski.add_instrument("earth", self.instrument)

        # Set the number of photon packages
        self.ski.setpackages(self.config.packages)

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid("wavelengths_lowres.txt")

        # Set the stellar and dust components
        self.set_components()

        # Set transient dust emissivity
        self.ski.set_transient_dust_emissivity()

        # Set the dust grid
        self.ski.set_dust_grid(self.lowres_dust_grid)

        # Set all-cells dust library
        self.ski.set_allcells_dust_lib()

        # Dust self-absorption
        if self.config.selfabsorption: self.ski.enable_selfabsorption()
        else: self.ski.disable_selfabsorption()

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
        fluxdensity = self.parameters.bulge.fluxdensity # In Jy

        # Convert the flux density into a spectral luminosity
        luminosity = fluxdensity_to_luminosity(fluxdensity, self.i1.pivotwavelength() * Unit("micron"), self.parameters.distance)

        # Get the spectral luminosity in solar units
        #luminosity = luminosity.to(self.sun_i1).value

        # Set the parameters of the bulge
        self.ski.set_stellar_component_geometry("Evolved stellar bulge", self.bulge)
        self.ski.set_stellar_component_sed("Evolved stellar bulge", bulge_template, bulge_age, bulge_metallicity) # SED
        #self.ski.set_stellar_component_luminosity("Evolved stellar bulge", luminosity, self.i1) # normalization by band
        self.ski.set_stellar_component_luminosity("Evolved stellar bulge", luminosity, self.i1.centerwavelength() * Unit("micron"))

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
        scale_height = self.parameters.disk.hr / 8.26 # De Geyter et al. 2014

        # Get the 3.6 micron flux density with the bulge subtracted
        i1_index = tables.find_index(self.observed_sed.table, "I1", "Band")
        fluxdensity = self.observed_sed.table["Flux"][i1_index] * Unit("Jy") - self.parameters.bulge.fluxdensity

        # Convert the flux density into a spectral luminosity
        luminosity = fluxdensity_to_luminosity(fluxdensity, self.i1.pivotwavelength() * Unit("micron"), self.parameters.distance)

        # Get the spectral luminosity in solar units
        #luminosity = luminosity.to(self.sun_i1).value

        # Set the parameters of the evolved stellar component
        deprojection = self.deprojection.copy()
        deprojection.filename = "old_stars.fits"
        deprojection.scale_height = scale_height
        self.deprojections["Old stars"] = deprojection

        # Adjust the ski file
        self.ski.set_stellar_component_geometry("Evolved stellar disk", deprojection)
        self.ski.set_stellar_component_sed("Evolved stellar disk", disk_template, disk_age, disk_metallicity) # SED
        #self.ski.set_stellar_component_luminosity("Evolved stellar disk", luminosity, self.i1) # normalization by band
        self.ski.set_stellar_component_luminosity("Evolved stellar disk", luminosity, self.i1.centerwavelength() * Unit("micron"))

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
        fuv_index = tables.find_index(self.observed_sed.table, "FUV", "Band")
        fluxdensity = 2. * self.observed_sed.table["Flux"][fuv_index] * Unit("Jy")

        # Convert the flux density into a spectral luminosity
        luminosity = fluxdensity_to_luminosity(fluxdensity, self.fuv.pivotwavelength() * Unit("micron"), self.parameters.distance)

        # Get the spectral luminosity in solar units
        #luminosity = luminosity.to(self.sun_fuv).value

        # Set the parameters of the young stellar component
        deprojection = self.deprojection.copy()
        deprojection.filename = "young_stars.fits"
        deprojection.scale_height = scale_height
        self.deprojections["Young stars"] = deprojection

        # Adjust the ski file
        self.ski.set_stellar_component_geometry("Young stars", deprojection)
        self.ski.set_stellar_component_sed("Young stars", young_template, young_age, young_metallicity) # SED
        #self.ski.set_stellar_component_luminosity("Young stars", luminosity, self.fuv) # normalization by band
        self.ski.set_stellar_component_luminosity("Young stars", luminosity, self.fuv.centerwavelength() * Unit("micron"))

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
        luminosity = mappings.luminosity_for_filter(self.fuv)
        #luminosity = luminosity.to(self.sun_fuv).value

        # Set the parameters of the ionizing stellar component
        deprojection = self.deprojection.copy()
        deprojection.filename = "ionizing_stars.fits"
        deprojection.scale_height = scale_height
        self.deprojections["Ionizing stars"] = deprojection

        # Adjust the ski file
        self.ski.set_stellar_component_geometry("Ionizing stars", deprojection)
        self.ski.set_stellar_component_mappingssed("Ionizing stars", ionizing_metallicity, ionizing_compactness, ionizing_pressure, ionizing_covering_factor) # SED
        #self.ski.set_stellar_component_luminosity("Ionizing stars", luminosity, self.fuv) # normalization by band
        self.ski.set_stellar_component_luminosity("Ionizing stars", luminosity, self.fuv.centerwavelength() * Unit("micron"))

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
        self.ski.set_dust_component_geometry(0, deprojection)
        self.ski.set_dust_component_themis_mix(0, hydrocarbon_pops, enstatite_pops, forsterite_pops) # dust mix
        self.ski.set_dust_component_mass(0, dust_mass) # dust mass

    # -----------------------------------------------------------------

    def adjust_ski_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting ski files for simulating the contribution of the various stellar components ...")

        # Loop over the different contributions, create seperate ski file instance
        contributions = ["old", "young", "ionizing"]
        component_names = {"old": ["Evolved stellar bulge", "Evolved stellar disk"],
                           "young": "Young stars",
                           "ionizing": "Ionizing stars"}
        for contribution in contributions:

            # Create a copy of the ski file instance
            ski = self.ski.copy()

            # Remove other stellar components
            ski.remove_stellar_components_except(component_names[contribution])

            # Add the ski file to the dictionary
            self.ski_contributions[contribution] = ski

    # -----------------------------------------------------------------

    def adjust_ski_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting ski files for generating simulated images ...")

        # Create a copy of the ski file instance
        self.ski_images = self.ski.copy()

        # Remove all instruments
        self.ski_images.remove_all_instruments()

        # Create frame instrument to generate datacube
        frame_instrument = FrameInstrument.from_projection(self.projection)

        # Add the frame instrument
        self.ski_images.add_instrument("earth", frame_instrument)

        # Add the SED instrument
        self.ski_images.add_instrument("earth", self.instrument)

    # -----------------------------------------------------------------

    def generate_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the dust grid data ...")

        # Generate the low-resolution dust grid data
        self.generate_low_res_grid()

        # Generate the high-resolution dust grid data
        self.generate_high_res_grid()

    # -----------------------------------------------------------------

    def generate_low_res_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the low-resolution grid data ...")

        # Generate the grid
        optical_depth = self.generate_grid(self.lowres_dust_grid, self.fit_grid_lowres_path)

        # Debugging
        log.debug("For the low-resolution dust grid, 90% of the cells have an optical depth smaller than " + str(optical_depth))

    # -----------------------------------------------------------------

    def generate_high_res_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining ideal optical depth criterion ...")

        # Calculate the optical depth (at 90% percentile) for the current grid parameters
        optical_depth = self.generate_grid(self.highres_dust_grid, self.fit_grid_highres_path)

        # Adapt the maximal optical depth criterion
        self.highres_dust_grid.max_optical_depth = optical_depth

        # Inform the user
        log.info("Generating the high-resolution grid data ...")

        # Rerun the simulation
        optical_depth = self.generate_grid(self.highres_dust_grid, self.fit_grid_highres_path)

        # Debugging
        log.debug("For the high-resolution grid, 90% of the cells have an optical depth smaller than " + str(optical_depth))

    # -----------------------------------------------------------------

    def generate_grid(self, grid, output_path):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Running a simulation just to generate the dust grid data files ...")

        # Create a copy of the ski file
        ski = self.ski.copy()

        # Set the dust grid
        ski.set_dust_grid(grid)

        # Convert to oligochromatic simulation
        ski.to_oligochromatic([1. * Unit("micron")])

        # Remove the instrument system
        ski.remove_instrument_system()

        # Set the number of photon packages to zero
        ski.setpackages(0)

        # Disable all writing options, except the one for writing the dust grid and cell properties
        ski.disable_all_writing_options()
        ski.set_write_grid()
        ski.set_write_cell_properties()

        # Write the ski file
        ski_path = fs.join(output_path, self.galaxy_name + ".ski")
        ski.saveto(ski_path)

        # Create the local SKIRT execution context
        skirt = SkirtExec()

        # Create the SKIRT arguments object
        arguments = SkirtArguments()
        arguments.ski_pattern = ski_path
        arguments.input_path = self.fit_in_path
        arguments.output_path = output_path

        # Run SKIRT to generate the dust grid data files
        skirt.run(arguments)

        # Determine the path to the cell properties file
        cellprops_path = fs.join(output_path, self.galaxy_name + "_ds_cellprops.dat")

        # Get the optical depth for which 90% of the cells have a smaller value
        optical_depth = None
        for line in reversed(open(cellprops_path).readlines()):
            if "of the cells have optical depth smaller than" in line:
                optical_depth = float(line.split("than: ")[1])
                break

        # Return the optical depth
        return optical_depth

    # -----------------------------------------------------------------

    def calculate_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the weight to give to each band ...")

        # Create the table to contain the weights
        self.weights = tables.new([[], [], []], names=["Instrument", "Band", "Weight"], dtypes=["S5", "S7", "float64"])

        # Initialize lists to contain the filters of the different wavelength ranges
        uv_bands = []
        optical_bands = []
        nir_bands = []
        mir_bands = []
        fir_bands = []
        submm_bands = []

        # Set the number of groups
        number_of_groups = 6

        # Loop over the entries in the observed fluxes table
        for i in range(len(self.observed_sed.table)):

            instrument = self.observed_sed.table["Instrument"][i]
            band = self.observed_sed.table["Band"][i]

            # Construct filter
            fltr = Filter.from_instrument_and_band(instrument, band)

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
        number_of_data_points = len(self.observed_sed.table)
        uv_weight = 1. / (len(uv_bands) * number_of_groups) * number_of_data_points
        optical_weight = 1. / (len(optical_bands) * number_of_groups) * number_of_data_points
        nir_weight = 1. / (len(nir_bands) * number_of_groups) * number_of_data_points
        mir_weight = 1. / (len(mir_bands) * number_of_groups) * number_of_data_points
        fir_weight = 1. / (len(fir_bands) * number_of_groups) * number_of_data_points
        submm_weight = 1. / (len(submm_bands) * number_of_groups) * number_of_data_points

        #print("UV", len(uv_bands), uv_weight)
        #print("Optical", len(optical_bands), optical_weight)
        #print("NIR", len(nir_bands), nir_weight)
        #print("MIR", len(mir_bands), mir_weight)
        #print("FIR", len(fir_bands), fir_weight)
        #print("Submm", len(submm_bands), submm_weight)

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

        # Write the ski file
        self.write_ski_file()

        # Write the ski files for simulating the contributions of the various stellar components
        self.write_ski_files_contributions()

        # Write the ski file for generating simulated images
        self.write_ski_file_images()

        # Write the weights table
        self.write_weights()

        # Write the geometries
        self.write_geometries()

        # Write the wavelength grids
        self.write_wavelength_grids()

        # Write the dust grids
        self.write_dust_grids()

    # -----------------------------------------------------------------



    # -----------------------------------------------------------------

    def write_ski_file(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file to " + self.fit_ski_path + " ...")

        # Save the ski file to the specified location
        self.ski.saveto(self.fit_ski_path)

    # -----------------------------------------------------------------

    def write_ski_files_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski files for simulating the contribution of the various stellar components ...")

        # Loop over the ski files
        for contribution in self.ski_contributions:

            # Determine the path to the ski file
            ski_path = fs.join(self.fit_best_contribution_paths[contribution], self.galaxy_name + ".ski")

            # Write the ski file
            self.ski_contributions[contribution].saveto(ski_path)

    # -----------------------------------------------------------------

    def write_ski_file_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file for creating simulated images ...")

        # Determine the path to the ski file
        ski_path = fs.join(self.fit_best_images_path, self.galaxy_name + ".ski")

        # Write the ski file
        self.ski_images.saveto(ski_path)

    # -----------------------------------------------------------------

    def write_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with weights to " + self.weights_table_path + " ...")

        # Write the table with weights
        tables.write(self.weights, self.weights_table_path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def write_geometries(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the Sersic model for the bulge and the deprojection model for the other components ...")

        # Write the bulge model
        bulge_path = fs.join(self.fit_geometries_path, "bulge.mod")
        self.bulge.save(bulge_path)

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
        table_path = fs.join(self.fit_wavelength_grids_path, "grids.dat")
        tables.write(self.wg_generator.table, table_path)

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
        table_path = fs.join(self.fit_dust_grids_path, "grids.dat")
        tables.write(self.dg_generator.table, table_path)

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
