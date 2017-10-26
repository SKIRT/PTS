#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractproperty
import numpy as np
from collections import defaultdict

# Import astronomical modules
from astropy.io.fits import Header
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from pts.core.test.implementation import TestImplementation
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.core.basics.log import log
from pts.modeling.basics.models import load_3d_model
from pts.modeling.basics.properties import GalaxyProperties
from pts.magic.core.frame import Frame
from pts.core.filter.broad import BroadBandFilter
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.core.simulation.skifile import SkiFile
from pts.modeling.basics.models import DeprojectionModel3D
from pts.core.plot.wavelengthgrid import WavelengthGridPlotter
from pts.core.plot.transmission import TransmissionPlotter
from pts.modeling.basics.models import load_2d_model
from pts.modeling.modeling.galaxy import fitting_filter_names
from pts.core.units.quantity import PhotometricQuantity
from pts.core.prep.wavelengthgrids import WavelengthGridGenerator
from pts.core.prep.dustgrids import DustGridGenerator
from pts.core.units.parsing import parse_quantity
from pts.magic.region.list import SkyRegionList
from pts.core.filter.filter import parse_filter
from pts.core.units.parsing import parse_unit as u
from pts.core.tools import sequences, parsing, stringify
from pts.modeling.config.parameters import parsing_types_for_parameter_types
from pts.modeling.config.parameters import default_units
from pts.core.remote.moderator import PlatformModerator
from pts.core.tools.stringify import tostr
from pts.evolve.analyse.database import get_scores_named_individuals, load_database
from pts.evolve.analyse.statistics import get_best_score_for_generation, load_statistics
from pts.do.commandline import Command
from pts.evolve.optimize.tables import ElitismTable
from pts.core.prep.templates import get_pan_template
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

# Determine the path to the dropbox path and the path of the directory with the data for M81
m81_data_path = fs.join(introspection.get_dropbox_tests_pts_path_for_subproject("modeling"), "M81")

# -----------------------------------------------------------------

models_path = fs.join(m81_data_path, "models")
disk2d_path = fs.join(models_path, "disk.mod")
bulge2d_path = fs.join(models_path, "bulge.mod")
disk2d_model = load_2d_model(disk2d_path)
bulge2d_model = load_2d_model(bulge2d_path)

# -----------------------------------------------------------------

instrument_name = "earth"

# -----------------------------------------------------------------

titles = dict()
titles["bulge"] = "Evolved stellar bulge"
titles["old"] = "Evolved stellar disk"
titles["young"] = "Young stars"
titles["ionizing"] = "Ionizing stars"
titles["dust"] = "Dust disk"

# -----------------------------------------------------------------

# Bulge
bulge_template = "BruzualCharlot"
bulge_age = 10
bulge_metallicity = 0.03

# Get the flux density of the bulge
bulge_fluxdensity = bulge2d_model.fluxdensity

# -----------------------------------------------------------------

# Old stellar disk
disk_template = "BruzualCharlot"
disk_age = 8
# disk_metallicity = 0.02
disk_metallicity = 0.03

# Get the scale height
old_scale_height = disk2d_model.scalelength / 8.26  # De Geyter et al. 2014

# Get the 3.6 micron flux density with the bulge subtracted
total_i1_fluxdensity = PhotometricQuantity(10.6552814592, "Jy")
old_fluxdensity = total_i1_fluxdensity - bulge_fluxdensity

# -----------------------------------------------------------------

# Young stellar disk
young_template = "BruzualCharlot"
young_age = 0.1
# young_metallicity = 0.02
young_metallicity = 0.03

# Get the scale height
# scale_height = 150 * Unit("pc") # first models
#young_scale_height = 100. * u("pc")  # M51
young_scale_height = 0.5 * old_scale_height

# -----------------------------------------------------------------

# Ionizing stellar disk
ionizing_metallicity = 0.03  # XU KONG et al. 2000
ionizing_compactness = 6.
ionizing_pressure = 1e12 * u("K/m3")
ionizing_covering_factor = 0.2

# Get the scale height
# scale_height = 150 * Unit("pc") # first models
#ionizing_scale_height = 100. * u("pc")  # M51
ionizing_scale_height = 0.25 * old_scale_height

# Convert the SFR into a FUV luminosity
sfr = 0.8  # The star formation rate # see Perez-Gonzalez 2006 (mentions Devereux et al 1995)

# -----------------------------------------------------------------

# fuv_young:6.0068695608165e+36 W/micron
# fuv_ionizing:2.4590756925069244e+33 W/micron]
# 15450820.890962543 Msun

fuv_young = PhotometricQuantity(1e36, "W/micron")
fuv_ionizing = PhotometricQuantity(2.5e33, "W/micron")
dust_mass = parse_quantity("1.5e7 Msun")

# scale_height = 260.5 * Unit("pc") # first models
#dust_scale_height = 200. * u("pc")  # M51
dust_scale_height = 0.25 * old_scale_height

#dust_mass = 1.5e7 * u("Msun")

hydrocarbon_pops = 25
silicate_pops = 25

# -----------------------------------------------------------------

old_filename = "old_stars.fits"
young_filename = "young_stars.fits"
ionizing_filename = "ionizing_stars.fits"
dust_filename = "dust.fits"

# -----------------------------------------------------------------

# Define the possible free parameters
possible_free_parameters = ["dust_mass", "fuv_young", "fuv_ionizing", "distance"]
default_free_parameters = ["dust_mass", "fuv_young", "fuv_ionizing"]

# Define the free parameter types
free_parameter_types = dict()

# Parameter types
free_parameter_types["dust_mass"] = "mass"
free_parameter_types["fuv_young"] = "spectral luminosity density"
free_parameter_types["fuv_ionizing"] = "spectral luminosity density"
free_parameter_types["distance"] = "length"

# Define free parameter units
free_parameter_units = dict()
for label in free_parameter_types:
    parameter_type = free_parameter_types[label]
    free_parameter_units[label] = default_units[parameter_type]

# Define the number of digits
parameter_ndigits = dict()
for label in free_parameter_types:
    parameter_ndigits[label] = 3

# Absolute ski file parameters
free_parameters_absolute_paths = dict()

# Stellar component parameters
free_parameters_relative_stellar_component_paths = dict()
free_parameters_relative_stellar_component_paths["fuv_young"] = ("normalization/SpectralLuminosityStellarCompNormalization/luminosity", titles["young"])
free_parameters_relative_stellar_component_paths["fuv_ionizing"] = ("normalization/SpectralLuminosityStellarCompNormalization/luminosity", titles["ionizing"])

# Dust component parameters
free_parameters_relative_dust_component_paths = dict()
free_parameters_relative_dust_component_paths["dust_mass"] = ("normalization/DustMassDustCompNormalization/dustMass", titles["dust"])

# Instrument parameters
free_parameters_relative_instruments_paths = dict()
free_parameters_relative_instruments_paths["distance"] = ("distance", instrument_name)

# Free parameter descriptions
free_parameter_descriptions = dict()
free_parameter_descriptions["dust_mass"] = "total dust mass"
free_parameter_descriptions["fuv_young"] = "FUV spectral luminosity of the young stellar component"
free_parameter_descriptions["fuv_ionizing"] = "FUV spectral luminosity of the ionizing stellar component"
free_parameter_descriptions["distance"] = "galaxy distance"

# -----------------------------------------------------------------

seds_path = fs.join(m81_data_path, "seds")
dustpedia_sed_path = fs.join(seds_path, "DustPedia.dat")

# -----------------------------------------------------------------

reference_wavelength_grid_filename = "wavelengths.txt"

# -----------------------------------------------------------------

class M81TestBase(TestImplementation):

    """
    This class runs the test on M81, but by only adjusting the normalizations (not by creating a model),
    and fitting to a mock observed SED
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor
        :param kwargs:
        """

        # Call the constructor of the base class
        super(M81TestBase, self).__init__(*args, **kwargs)

        # The platform moderator
        self.moderator = None

        # The galaxy properties
        self.properties = None

        # Bulge and disk model
        self.bulge = None
        self.disk = None

        # The input maps
        self.maps = dict()

        # The deprojections
        self.deprojections = dict()

        # The instrument
        self.instrument = None

        # Path to the ski file for the reference simulation
        self.reference_path = None
        self.reference_wcs_path = None
        self.reference_ski_path = None
        self.simulation_input_path = None
        self.simulation_output_path = None
        self.simulation_extract_path = None
        self.simulation_plot_path = None
        self.simulation_misc_path = None
        self.wavelength_grid_path = None

        # The plot path
        self.plot_path = None

        # The reference ski file
        self.ski = None

        # The wcs for the reference simulation
        self.wcs = None

        # The wavelength grid
        self.wavelength_grid = None

        # The dust grid
        self.dust_grid = None

        # The simulation launcher
        self.launcher = None

        # The flux calculator
        self.flux_calculator = None

        # The modeler
        self.modeler = None

        # The real parameter values
        self.real_parameter_values = dict()

        # The best parameter values
        self.best_parameter_values = dict()

        # The chi squared of the best model
        self.best_chi_squared = None

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_filter(self):

        """
        This function ...
        :return:
        """

        return BroadBandFilter("GALEX FUV")

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_filter(self):

        """
        This function ...
        :return:
        """

        return BroadBandFilter("IRAC I1")

    # -----------------------------------------------------------------

    @property
    def galaxy_name(self):

        """
        This function ...
        :return:
        """

        return self.properties.name

    # -----------------------------------------------------------------

    @property
    def galaxy_center(self):

        """
        This function ...
        :return:
        """

        return self.properties.center

    # -----------------------------------------------------------------

    @property
    def galaxy_position_angle(self):

        """
        This function ...
        :return:
        """

        return self.properties.position_angle

    # -----------------------------------------------------------------

    @property
    def galaxy_inclination(self):

        """
        This function ...
        :return:
        """

        return self.properties.inclination

    # -----------------------------------------------------------------

    @property
    def galaxy_distance(self):

        """
        This function ...
        :return:
        """

        return self.properties.distance

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(M81TestBase, self).setup(**kwargs)

        # Reference base path
        self.reference_path = fs.create_directory_in(self.path, "ref")

        # Reference wcs path
        self.reference_wcs_path = fs.join(self.reference_path, "wcs.txt")

        # Reference ski path
        self.reference_ski_path = fs.join(self.reference_path, "M81.ski")

        # Determine the simulation input and output path
        self.simulation_input_path = fs.create_directory_in(self.reference_path, "in")
        self.simulation_output_path = fs.create_directory_in(self.reference_path, "out")
        self.simulation_extract_path = fs.create_directory_in(self.reference_path, "extr")
        self.simulation_plot_path = fs.create_directory_in(self.reference_path, "plot")
        self.simulation_misc_path = fs.create_directory_in(self.reference_path, "misc")

        # Determine the path to the wavelength grid file
        self.wavelength_grid_path = fs.join(self.simulation_input_path, reference_wavelength_grid_filename)

        # Create the plot path
        self.plot_path = fs.create_directory_in(self.path, "plot")

        # Set the execution platforms
        self.set_platforms()

    # -----------------------------------------------------------------

    def set_platforms(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining execution platforms ...")

        # Setup the platform moderator
        self.moderator = PlatformModerator()

        # Set platform for the reference simulation
        if self.config.host_ids is None: self.moderator.add_local("reference")
        else: self.moderator.add_single("reference", self.config.host_ids)

        # Set platform(s) for fitting (simulations)
        if self.config.host_ids is None: self.moderator.add_local("fitting")
        else: self.moderator.add_ensemble("fitting", self.config.host_ids)

        # Run the moderator
        self.moderator.run()

    # -----------------------------------------------------------------

    def load_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the galaxy properties ...")

        # Determine the path
        path = fs.join(m81_data_path, "properties.dat")

        # Load
        self.properties = GalaxyProperties.from_file(path)

    # -----------------------------------------------------------------

    def load_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the components ...")

        # Determine paths
        path = fs.join(m81_data_path, "components")
        bulge_path = fs.join(path, "bulge.mod")
        disk_path = fs.join(path, "disk.mod")

        # Load bulge model
        self.bulge = load_3d_model(bulge_path)

        # Load disk model
        self.disk = load_3d_model(disk_path)

    # -----------------------------------------------------------------

    def load_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the input maps ...")

        # Determine path to maps directory
        maps_path = fs.join(m81_data_path, "maps")

        # Determine the path to the header file
        header_path = fs.join(maps_path, "header.txt")
        header = Header.fromtextfile(header_path)
        wcs = CoordinateSystem(header=header)

        # Old stars
        old_map_path = fs.join(maps_path, old_filename)
        old_map = Frame.from_file(old_map_path)
        old_map.wcs = wcs
        self.maps["old"] = old_map

        # young stars
        young_map_path = fs.join(maps_path, young_filename)
        young_map = Frame.from_file(young_map_path)
        young_map.wcs = wcs
        self.maps["young"] = young_map

        # Ionizing stars
        ionizing_map_path = fs.join(maps_path, ionizing_filename)
        ionizing_map = Frame.from_file(ionizing_map_path)
        ionizing_map.wcs = wcs
        self.maps["ionizing"] = ionizing_map

        # Dust
        dust_map_path = fs.join(maps_path, dust_filename)
        dust_map = Frame.from_file(dust_map_path)
        dust_map.wcs = wcs
        self.maps["dust"] = dust_map

    # -----------------------------------------------------------------

    def create_deprojections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the deprojections ...")

        # Set deprojection of old stars
        wcs = self.maps["old"].wcs
        if wcs is None: raise IOError("The map of old stars has no WCS information")
        self.deprojections["old"] = DeprojectionModel3D.from_wcs(wcs, self.galaxy_center, self.galaxy_distance, self.galaxy_position_angle, self.galaxy_inclination, old_filename, old_scale_height)

        # Set deprojection of young stars
        wcs = self.maps["young"].wcs
        if wcs is None: raise IOError("The map of young stars has no WCS information")
        self.deprojections["young"] = DeprojectionModel3D.from_wcs(wcs, self.galaxy_center, self.galaxy_distance, self.galaxy_position_angle, self.galaxy_inclination, young_filename, young_scale_height)

        # Set deprojection of ionizing stars
        wcs = self.maps["ionizing"].wcs
        if wcs is None: raise IOError("The map of ionizing stars has no WCS information")
        self.deprojections["ionizing"] = DeprojectionModel3D.from_wcs(wcs, self.galaxy_center, self.galaxy_distance, self.galaxy_position_angle, self.galaxy_inclination, ionizing_filename, ionizing_scale_height)

        # Set deprojection of dust map
        wcs = self.maps["dust"].wcs
        if wcs is None: raise IOError("The map of old stars has no WCS information")
        self.deprojections["dust"] = DeprojectionModel3D.from_wcs(wcs, self.galaxy_center, self.galaxy_distance, self.galaxy_position_angle, self.galaxy_inclination, dust_filename, dust_scale_height)

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Create the wavelength generator
        generator = WavelengthGridGenerator()

        # Set input
        input_dict = dict()
        input_dict["ngrids"] = 1
        input_dict["npoints"] = self.config.nwavelengths
        input_dict["fixed"] = [self.i1_filter.pivot, self.fuv_filter.pivot]
        input_dict["add_emission_lines"] = True
        input_dict["lines"] = ["Halpha"] # only the H-alpha line is of importance
        input_dict["min_wavelength"] = self.config.wavelength_range.min
        input_dict["max_wavelength"] = self.config.wavelength_range.max
        input_dict["filters"] = [parse_filter(string) for string in fitting_filter_names]

        # Run the generator
        generator.run(**input_dict)

        # Set the wavelength grid
        self.wavelength_grid = generator.single_grid

    # -----------------------------------------------------------------

    def create_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust grid ...")

        # Create the dust grid generator
        generator = DustGridGenerator()

        # Determine truncation ellipse
        disk_ellipse_path = fs.join(m81_data_path, "components", "disk.reg")
        disk_ellipse = SkyRegionList.from_file(disk_ellipse_path)[0]
        truncation_ellipse = self.config.physical_domain_disk_ellipse_factor * disk_ellipse

        # Determine the radius of the galaxy
        semimajor_angular = truncation_ellipse.semimajor  # semimajor axis length of the sky ellipse
        radius_physical = (semimajor_angular * self.galaxy_distance).to("pc", equivalencies=dimensionless_angles())

        # Set properties
        generator.grid_type = "bintree"  # set grid type
        generator.x_radius = radius_physical
        generator.y_radius = radius_physical
        generator.z_radius = 2. * u("kpc")

        # Set input
        input_dict = dict()
        input_dict["ngrids"] = 1
        input_dict["scale"] = self.config.dust_grid_relative_scale * self.deprojections["dust"].pixelscale # in pc
        input_dict["level"] = self.config.dust_grid_min_level
        input_dict["mass_fraction"] = self.config.dust_grid_max_mass_fraction

        # Generate the grid
        generator.run(**input_dict)

        # Set the dust grid
        self.dust_grid = generator.single_grid

    # -----------------------------------------------------------------

    def create_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the ski file ...")

        # Load the ski file template
        self.ski = get_pan_template()

        # Set components
        self.set_components()

        # Add the instrument
        self.ski.add_instrument(instrument_name, self.instrument)

        # Set the number of photon packages
        self.ski.setpackages(self.config.npackages)

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid(reference_wavelength_grid_filename)

        # Set the dust emissivity
        if self.config.transient_heating: self.ski.set_transient_dust_emissivity()
        else: self.ski.set_grey_body_dust_emissivity()

        # Set the lowest-resolution dust grid
        self.ski.set_dust_grid(self.dust_grid)

        # Set all-cells dust library
        self.ski.set_allcells_dust_lib()

        # Set the dust selfabsorption
        if self.config.selfabsorption: self.ski.enable_selfabsorption()
        else: self.ski.disable_selfabsorption()

        # Disable all writing options
        self.ski.disable_all_writing_options()

    # -----------------------------------------------------------------

    def set_components(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Setting the ski components ...")

        # Set bulge
        self.set_bulge()

        # Set old
        self.set_old()

        # Set young
        self.set_young()

        # Set ionizing
        self.set_ionizing()

        # Set dust
        self.set_dust()

    # -----------------------------------------------------------------

    def set_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the old stellar bulge component ...")

        # Get the title for this component
        title = titles["bulge"]

        # Create the new component
        self.ski.create_new_stellar_component(title)

        # Set the geometry
        self.ski.set_stellar_component_geometry(title, self.bulge)

        # Set the SED
        # component_id, template, age, metallicity
        self.ski.set_stellar_component_sed(title, bulge_template, bulge_age, bulge_metallicity)

        # Convert the flux density into a spectral luminosity
        luminosity = bulge_fluxdensity.to("W/micron", fltr=self.i1_filter, distance=self.galaxy_distance)

        # Set the normalization
        # luminosity, filter_or_wavelength=None
        self.ski.set_stellar_component_luminosity(title, luminosity, self.i1_filter.pivot)

    # -----------------------------------------------------------------

    def set_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the old stellar disk component ...")

        # Get the title
        title = titles["old"]

        # Create the new component
        self.ski.create_new_stellar_component(title)

        # Set the geometry
        self.ski.set_stellar_component_geometry(title, self.deprojections["old"])

        # Set the SED
        self.ski.set_stellar_component_sed(title, disk_template, disk_age, disk_metallicity)

        # Convert the flux density into a spectral luminosity
        luminosity = old_fluxdensity.to("W/micron", fltr=self.i1_filter, distance=self.galaxy_distance)

        # Set the normalization
        self.ski.set_stellar_component_luminosity(title, luminosity, self.i1_filter.pivot)

    # -----------------------------------------------------------------

    def set_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the young stellar disk component ...")

        # Get the title
        title = titles["young"]

        # Create the new component
        self.ski.create_new_stellar_component(title)

        # Set the geometry
        self.ski.set_stellar_component_geometry(title, self.deprojections["young"])

        # Set the SED
        self.ski.set_stellar_component_sed(title, young_template, young_age, young_metallicity)

        # Set the normalization
        self.ski.set_stellar_component_luminosity(title, fuv_young, self.fuv_filter.pivot)

    # -----------------------------------------------------------------

    def set_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the ionizing stellar disk component ...")

        # Get the title
        title = titles["ionizing"]

        # Create the new component
        self.ski.create_new_stellar_component(title)

        # Set the geometry
        self.ski.set_stellar_component_geometry(title, self.deprojections["ionizing"])

        # Set the SED
        # metallicity, compactness, pressure, covering_factor
        self.ski.set_stellar_component_mappingssed(title, ionizing_metallicity, ionizing_compactness, ionizing_pressure, ionizing_covering_factor)

        # Set the normalization
        self.ski.set_stellar_component_luminosity(title, fuv_ionizing, self.fuv_filter.pivot)

    # -----------------------------------------------------------------

    def set_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust component ...")

        # Get the title
        title = titles["dust"]

        # Create the new component
        self.ski.create_new_dust_component(title)

        # Set the geometry
        self.ski.set_dust_component_geometry(title, self.deprojections["dust"])

        # Set the mix
        self.ski.set_dust_component_themis_mix(title, hydrocarbon_pops, silicate_pops, write_mix=False, write_mean_mix=False, write_size=False)

        # Set the normalization
        self.ski.set_dust_component_normalization(title, dust_mass)

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        # Save
        self.ski.saveto(self.reference_ski_path)

    # -----------------------------------------------------------------

    def write_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the simulation input ...")

        # Write wavelength grid
        self.wavelength_grid.to_skirt_input(self.wavelength_grid_path)

        # Write maps
        self.maps["old"].saveto(fs.join(self.simulation_input_path, old_filename))
        self.maps["young"].saveto(fs.join(self.simulation_input_path, young_filename))
        self.maps["ionizing"].saveto(fs.join(self.simulation_input_path, ionizing_filename))
        self.maps["dust"].saveto(fs.join(self.simulation_input_path, dust_filename))

    # -----------------------------------------------------------------

    def plot_wavelengths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelengths ...")

        # Create the plotter
        plotter = WavelengthGridPlotter()

        # Add the wavelength grid
        plotter.add_wavelength_grid(self.wavelength_grid, "reference simulation")

        # Determine the plot path
        path = fs.join(self.reference_path, "wavelengths.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_filters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the filters ...")

        # Create the plotter
        plotter = TransmissionPlotter()

        # Set the filters
        plotter.config.filters = fitting_filter_names

        # Add the wavelengths of the wavelength grid
        for wavelength in self.wavelength_grid.wavelengths(): plotter.add_wavelength(wavelength)

        # Determine the plot path
        path = fs.join(self.reference_path, "filters.pdf")

        # Run the plotter
        plotter.run(output=path)

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return: 
        """

        return self.config.free_parameters

    # -----------------------------------------------------------------

    def get_real_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the real parameter values ...")

        # Store the different values encountered in the ski file
        values = defaultdict(list)

        # Add the labels for the free parameters
        # Loop over the free parameters
        for parameter_name in self.config.free_parameters:

            parameter_type = free_parameter_types[parameter_name]
            parsing_type = parsing_types_for_parameter_types[parameter_type]

            # Get the parsing function for this parameter
            parser = getattr(parsing, parsing_type)

            # Search in the absolute parameters
            if parameter_name in free_parameters_absolute_paths:

                # Determine the path to the property
                path = free_parameters_absolute_paths

                # Get the current value
                value = parser(self.ski.get_value_for_path(path))

                # Set the value
                values[parameter_name].append(value)

            # Search in the stellar components
            if parameter_name in free_parameters_relative_stellar_component_paths:

                # Determine the relative path to the property and the stellar component name
                path, component_name = free_parameters_relative_stellar_component_paths[parameter_name]

                if component_name is not None:

                    # Get the stellar component
                    stellar_component = self.ski.get_stellar_component(component_name)

                    # Get the current value
                    value = parser(self.ski.get_value_for_path(path, stellar_component))

                    # Set the value
                    values[parameter_name].append(value)

                else:

                    # Loop over the stellar components
                    for component_id in self.ski.get_stellar_component_ids():

                        # Get the stellar component
                        stellar_component = self.ski.get_stellar_component(component_id)

                        # Get the current value
                        value = parser(self.ski.get_value_for_path(path, stellar_component))

                        # Set the value
                        values[parameter_name].append(value)

            # Search in the dust components
            if parameter_name in free_parameters_relative_dust_component_paths:

                # Determine the relative path to the property and the dust component name
                path, component_name = free_parameters_relative_dust_component_paths[parameter_name]

                if component_name is not None:

                    # Get the dust component
                    dust_component = self.ski.get_dust_component(component_name)

                    # Get the current value
                    value = parser(self.ski.get_value_for_path(path, dust_component))

                    # Set the value
                    values[parameter_name].append(value)

                else:

                    # Loop over the dust components
                    for component_id in self.ski.get_dust_component_ids():

                        # Get the dust component
                        dust_component = self.ski.get_dust_component(component_id)

                        # Get the current value
                        value = parser(self.ski.get_value_for_path(path, dust_component))

                        # Set the value
                        values[parameter_name].append(value)

            # Search in instruments
            if parameter_name in free_parameters_relative_instruments_paths:

                # Determine the relative path to the property and the instrument name
                path, instrument_name = free_parameters_relative_instruments_paths[parameter_name]

                if instrument_name is not None:

                    # Get the instrument
                    instrument = self.ski.get_instrument(instrument_name)

                    # Get the current value
                    value = parser(self.ski.get_value_for_path(path, instrument))

                    # Set the value
                    values[parameter_name].append(value)

                else:

                    # Loop over the instruments
                    for instrument_name in self.ski.get_instrument_names():

                        # Get the instruemnt
                        instrument = self.ski.get_instrument(instrument_name)

                        # Get the current value
                        value = parser(self.ski.get_value_for_path(path, instrument))

                        # Set the value
                        values[parameter_name].append(value)

        # Check whether we have only one value for each parameter
        for parameter_name in self.config.free_parameters:

            # Check if any
            if len(values[parameter_name]) == 0: raise ValueError("No parameter values for '" + parameter_name + "' were found in the ski file")

            # Check if all equal
            if not sequences.all_equal(values[parameter_name]): raise ValueError("Parameter values for '" + parameter_name + "' are not equal throughout the ski file")

            # Set the unique real parameter value
            self.real_parameter_values[parameter_name] = values[parameter_name][0]

        # Debugging
        log.debug("The real parameter values are: ")
        log.debug("")
        for parameter_name in self.real_parameter_values: log.debug(" - " + parameter_name + ": " + tostr(self.real_parameter_values[parameter_name], scientific=True, fancy=True, ndigits=parameter_ndigits[parameter_name]))
        log.debug("")

    # -----------------------------------------------------------------

    def get_best_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the best parameter values ...")

        # Get the best parameter values
        self.best_parameter_values, self.best_chi_squared = self.modeler.modeler.fitter.fitting_run.best_parameter_values_and_chi_squared

        # Debugging
        log.debug("The best parameter values are:")
        log.debug("")
        for parameter_name in self.best_parameter_values: log.debug(" - " + parameter_name + ": " + tostr(self.best_parameter_values[parameter_name], scientific=True, fancy=True, ndigits=parameter_ndigits[parameter_name]))
        log.debug("")

        # Debugging
        log.debug("The best chi squared value is " + str(self.best_chi_squared))

    # -----------------------------------------------------------------

    def get_initial_generation_name(self):

        """
        This function ...
        :return: 
        """

        return self.modeler.modeler.explorer.get_initial_generation_name()

    # -----------------------------------------------------------------

    def get_generation_name(self, generation_index):

        """
        This function ...
        :param generation_index: 
        :return: 
        """

        return self.modeler.modeler.explorer.get_genetic_generation_name(generation_index)

    # -----------------------------------------------------------------

    @property
    def generation_names(self):

        """
        This function ...
        :return: 
        """

        names = [self.get_initial_generation_name()]
        for index in range(self.config.ngenerations):
            names.append(self.get_generation_name(index))
        return names

    # -----------------------------------------------------------------

    @property
    def fitting_run_name(self):

        """
        This function ...
        :return: 
        """

        return self.modeler.modeler.fitting_run_name

    # -----------------------------------------------------------------

    @property
    def modeling_name(self):

        """
        This function ...
        :return: 
        """

        return self.galaxy_name

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This function ...
        :return: 
        """

        return fs.join(self.path, self.modeling_name)

    # -----------------------------------------------------------------

    @abstractproperty
    def modeling_environment(self):

        """
        This function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    @property
    def fitting_run(self):

        """
        This function ...
        :return: 
        """

        return self.modeler.modeler.fitter.fitting_run

    # -----------------------------------------------------------------

    @property
    def generations_table(self):

        """
        This function ...
        :return: 
        """

        return self.fitting_run.generations_table

    # -----------------------------------------------------------------

    def individuals_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return self.fitting_run.individuals_table_for_generation(generation_name)

    # -----------------------------------------------------------------

    def parameters_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return self.fitting_run.parameters_table_for_generation(generation_name)

    # -----------------------------------------------------------------

    def chi_squared_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return self.fitting_run.chi_squared_table_for_generation(generation_name)

    # -----------------------------------------------------------------

    def elitism_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return fs.join(self.fitting_run.generations_path, generation_name, "elitism.dat")

    # -----------------------------------------------------------------

    def elitism_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return ElitismTable.from_file(self.elitism_table_path_for_generation(generation_name))

    # -----------------------------------------------------------------

    def get_parameter_values_for_simulation(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name: 
        :param simulation_name: 
        :return: 
        """

        # Get the parameters table
        parameters_table = self.parameters_table_for_generation(generation_name)

        # Return the parameters
        return parameters_table.parameter_values_for_simulation(simulation_name)

    # -----------------------------------------------------------------

    @property
    def database_path(self):

        """
        This function ...
        :return: 
        """

        return self.modeler.modeler.explorer.database_path

    # -----------------------------------------------------------------

    @property
    def database(self):

        """
        THis function ...
        :return: 
        """

        return load_database(self.database_path)

    # -----------------------------------------------------------------

    @property
    def statistics_path(self):

        """
        This function ...
        :return: 
        """

        return self.modeler.modeler.explorer.statistics_path

    # -----------------------------------------------------------------

    @property
    def statistics(self):

        """
        This function ...
        :return: 
        """

        return load_statistics(self.statistics_path)

    # -----------------------------------------------------------------

    @property
    def last_genetic_generation_index(self):

        """
        This function ...
        :return: 
        """

        return self.fitting_run.last_genetic_generation_index

    # -----------------------------------------------------------------

    @property
    def last_genetic_generation_name(self):

        """
        This function ...
        :return: 
        """

        return self.fitting_run.last_genetic_generation_name

    # -----------------------------------------------------------------

    def get_best_score_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        # Get the scores table
        scores_table = self.chi_squared_table_for_generation(generation_name)

        # Return the best (lowest) chi squared
        return scores_table.best_chi_squared

    # -----------------------------------------------------------------

    def get_best_parameter_values_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        # Get the scores table
        scores_table = self.chi_squared_table_for_generation(generation_name)

        # Get the best simulation name
        simulation_name = scores_table.best_simulation_name

        # Get the parameter values for this simulation
        return self.get_parameter_values_for_simulation(generation_name, simulation_name)

    # -----------------------------------------------------------------

    def get_best_parameter_values_last_generation(self):

        """
        This function ...
        :return: 
        """

        return self.get_best_parameter_values_for_generation(self.last_genetic_generation_name)

    # -----------------------------------------------------------------

    def get_best_parameter_values_all_generations(self):

        """
        This function ...
        :return: 
        """

        best_chi_squared = None
        best_parameters = None

        # Loop over the generation names
        for generation_name in self.generation_names:

            # Load the chi squared table
            scores_table = self.chi_squared_table_for_generation(generation_name)

            # Get the name of the simulation with the lowest chi squared value
            simulation_name, chi_squared = scores_table.best_simulation_name_and_chi_squared

            # If chi squared is better or first
            if best_chi_squared is None or chi_squared < best_chi_squared:

                # Get parameters values
                values = self.get_parameter_values_for_simulation(generation_name, simulation_name)

                # Set best chi squared
                best_chi_squared = chi_squared
                best_parameters = values

        # Return the best parameter values
        return best_parameters

    # -----------------------------------------------------------------

    def check_best(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the best individual ...")

        # Get the values from the tables
        values_last_generation = self.get_best_parameter_values_last_generation()
        values_all_generations = self.get_best_parameter_values_all_generations()

        print("")
        for label in self.real_parameter_values:

            # Get the values form the tables
            value_last_generation = values_last_generation[label]
            value_all_generations = values_all_generations[label]

            # Get the values and calculate the difference
            value = self.best_parameter_values[label]
            real_value = self.real_parameter_values[label]
            absolute_difference = abs(value - real_value)
            relative_difference = absolute_difference / real_value

            # Get the

            print(label + ":")
            print(" - Best value in last generation: " + tostr(value_last_generation))
            print(" - Best value across all generations: " + tostr(value_all_generations))
            print(" - Real value: " + tostr(real_value))
            print(" - Best fitted value: " + tostr(value))
            print(" - Absolute difference: " + tostr(absolute_difference))
            print(" - Relative difference: " + tostr(relative_difference) + " (" + tostr(relative_difference * 100) + "%)")
            print("")

    # -----------------------------------------------------------------

    def check_database(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the database ...")

        # TODO: use GenerationPlatform.check_database()

    # -----------------------------------------------------------------

    def check_statistics(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the statistics table ...")

        # TODO: use GenerationPlatform.check_database()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the scores
        self.plot_scores()

        # Plot the heatmap
        self.plot_heatmap()

    # -----------------------------------------------------------------

    def plot_scores(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting the scores ...")

        # Settings
        settings = dict()
        settings["fitness"] = False
        settings["output"] = self.plot_path

        # Input
        input_dict = dict()
        input_dict["database"] = self.database

        # Plot
        command = Command("plot_scores", "plot the scores", settings, input_dict)
        plotter = self.run_command(command)

    # -----------------------------------------------------------------

    def plot_heatmap(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting heatmap ...")

        # Settings
        settings = dict()
        settings["fitness"] = False
        settings["output"] = self.plot_path

        # Input
        input_dict = dict()
        input_dict["database"] = self.database

        # Plot
        command = Command("plot_heat_map", "plot heatmap", settings, input_dict)
        plotter = self.run_command(command)

# -----------------------------------------------------------------
