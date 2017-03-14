#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import inspect
import numpy as np

# Import astronomical modules
from astropy.utils import lazyproperty
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.unit import parse_unit as u
from pts.core.test.implementation import TestImplementation
from pts.core.tools import introspection
from pts.core.simulation.skifile import SkiFile
from pts.core.tools.logging import log
from pts.do.commandline import Command
from pts.magic.misc.kernels import AnianoKernels
from pts.modeling.modeling.galaxy import fitting_filter_names
from pts.modeling.basics.properties import GalaxyProperties
from pts.modeling.basics.models import load_3d_model, load_2d_model, DeprojectionModel3D
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.modeling.basics.instruments import FullInstrument
from pts.core.prep.wavelengthgrids import WavelengthGridGenerator
from pts.core.prep.dustgrids import DustGridGenerator
from pts.core.basics.quantity import PhotometricQuantity, parse_quantity
from pts.core.filter.broad import BroadBandFilter
from pts.magic.core.frame import Frame
from pts.magic.region.list import SkyRegionList
from pts.core.filter.filter import parse_filter

# -----------------------------------------------------------------

# Determine path of this directory
this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "fitting the parameters of a model of M81"

# -----------------------------------------------------------------

# Determine the path to the template ski file for panchromatic simulations
dat_ski_path = fs.join(introspection.pts_dat_dir("modeling"), "ski")
pan_ski_path = fs.join(dat_ski_path, "pan.ski")

# Determine the path to the dropbox path and the path of the directory with the data for M81
dropbox_path = introspection.get_dropbox_path()
m81_data_path = fs.join(dropbox_path, "Data", "Tests", "PTS", "modeling", "M81")

# -----------------------------------------------------------------

models_path = fs.join(m81_data_path, "models")
disk2d_path = fs.join(models_path, "disk.mod")
bulge2d_path = fs.join(models_path, "bulge.mod")
disk2d_model = load_2d_model(disk2d_path)
bulge2d_model = load_2d_model(bulge2d_path)

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

fuv_young = PhotometricQuantity(6e36, "W/micron")
fuv_ionizing = PhotometricQuantity(2.5e33, "W/micron")
dust_mass = parse_quantity("1.5e7 Msun")

# scale_height = 260.5 * Unit("pc") # first models
#dust_scale_height = 200. * u("pc")  # M51
dust_scale_height = 0.25 * old_scale_height

#dust_mass = 1.5e7 * u("Msun")

hydrocarbon_pops = 25
enstatite_pops = 25
forsterite_pops = 25

# -----------------------------------------------------------------

old_filename = "old_stars.fits"
young_filename = "young_stars.fits"
ionizing_filename = "ionizing_stars.fits"
dust_filename = "dust.fits"

# -----------------------------------------------------------------

class M81Test(TestImplementation):

    """
    This class runs the test on M81
    """

    def __init__(self, config=None):

        """
        This function ...
        :param config:
        """

        # Call the constructor of the base class
        super(M81Test, self).__init__(config)

        # The galaxy properties
        self.properties = None

        # The disk and bulge models
        self.disk = None
        self.bulge = None

        # The maps
        self.maps = dict()

        # The deprojections
        self.deprojections = dict()

        # The instrument
        self.instrument = None

        # Path to the ski file for the reference simulation
        self.reference_path = None
        self.ski_path = None
        self.simulation_input_path = None
        self.simulation_output_path = None
        self.wavelength_grid_path = None

        # The reference ski file
        self.ski = None

        # The wcs for the reference simulation
        self.wcs = None

        # The simulation launcher
        self.launcher = None

        # The flux calculator
        self.flux_calculator = None

        # The image maker
        self.image_maker = None

        # The wavelength grid
        self.wavelength_grid = None

        # The dust grid
        self.dust_grid = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the properties
        self.load_properties()

        # 3. Load the components
        self.load_components()

        # Load the input maps
        self.load_maps()

        # 4. Load the wcs
        self.create_wcs()

        # Create the instrument
        self.create_instrument()

        # 5. Create wavelength grid
        self.create_wavelength_grid()

        # Create dust grid
        self.create_dust_grid()

        # Create deprojections
        self.create_deprojections()

        # 6. Create the ski file
        self.create_ski()

        # Write
        self.write()

        # 7. Launch reference simulation
        self.launch_reference()

        # 8. Make observed SED
        self.make_sed()

        # 9. Make observed images
        self.make_images()

        # 10. Setup modelling
        self.setup_modelling()

        # 11. Model
        self.model()

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
        super(M81Test, self).setup(**kwargs)

        # Reference base path
        self.reference_path = fs.create_directory_in(self.path, "ref")

        # Reference ski path
        self.ski_path = fs.join(self.reference_path, "M81.ski")

        # Determine the simulation input and output path
        self.simulation_input_path = fs.create_directory_in(self.reference_path, "in")
        self.simulation_output_path = fs.create_directory_in(self.reference_path, "out")

        # Determine the path to the wavelength grid file
        self.wavelength_grid_path = fs.join(self.path, "wavelengths.txt")

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
        THis function ...
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

        # Old stars
        old_map_path = fs.join(maps_path, old_filename)
        old_map = Frame.from_file(old_map_path)
        self.maps["old"] = old_map

        # young stars
        young_map_path = fs.join(maps_path, young_filename)
        young_map = Frame.from_file(young_map_path)
        self.maps["young"] = young_map

        # Ionizing stars
        ionizing_map_path = fs.join(maps_path, ionizing_filename)
        ionizing_map = Frame.from_file(ionizing_map_path)
        self.maps["ionizing"] = ionizing_map

        # Dust
        dust_map_path = fs.join(maps_path, dust_filename)
        dust_map = Frame.from_file(dust_map_path)
        self.maps["dust"] = dust_map

    # -----------------------------------------------------------------

    def create_wcs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the WCS ...")

        min_pixelscale = None

        # Determine the path to the headers directory
        headers_path = fs.join(m81_data_path, "headers")

        for path, name in fs.files_in_path(headers_path, extension="txt", returns=["path", "name"]):

            # Get WCS
            wcs = CoordinateSystem.from_header_file(path)

            # Adjust the pixelscale
            if min_pixelscale is None:
                min_pixelscale = wcs.pixelscale
                self.wcs = wcs
            elif min_pixelscale > wcs.pixelscale:
                min_pixelscale = wcs.pixelscale
                self.wcs = wcs

    # -----------------------------------------------------------------

    def create_wcs_not_working(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the WCS ...")

        # Determine the path to the headers directory
        headers_path = fs.join(m81_data_path, "headers")

        # Loop over the files in the directory
        ra_range = None
        dec_range = None
        min_pixelscale = None
        for path, name in fs.files_in_path(headers_path, extension="txt", returns=["path", "name"]):

            # Get the filter
            #fltr = parse_filter(name)

            # Get WCS
            wcs = CoordinateSystem.from_header_file(path)

            # Adjust RA range
            if ra_range is None: ra_range = wcs.ra_range
            else: ra_range.adjust(wcs.ra_range)

            # Adjust DEC range
            if dec_range is None: dec_range = wcs.dec_range
            else: dec_range.adjust(wcs.dec_range)

            # Adjust the pixelscale
            if min_pixelscale is None: min_pixelscale = wcs.pixelscale
            elif min_pixelscale > wcs.pixelscale: min_pixelscale = wcs.pixelscale

        # Create coordinate system
        # size, center_pixel, center_sky, pixelscale
        self.wcs = CoordinateSystem.from_ranges(ra_range, dec_range, min_pixelscale)

    # -----------------------------------------------------------------

    def create_instrument(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instrument ...")

        # Create a full instrument
        # center, distance, inclination, azimuth, position_angle
        azimuth = 0.0
        self.instrument = FullInstrument.from_wcs(self.wcs, self.galaxy_center, self.galaxy_distance, self.galaxy_inclination, azimuth, self.galaxy_position_angle)

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
        input_dict["npoints"] = 150
        input_dict["fixed"] = [self.i1_filter.pivot, self.fuv_filter.pivot]
        input_dict["add_emission_lines"] = True
        input_dict["min_wavelength"] = parse_quantity("0.1 micron")
        input_dict["max_wavelength"] = parse_quantity("1000 micron")
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
        disk_ellipse = SkyRegionList.from_file(disk_ellipse_path)
        truncation_ellipse = 0.82 * disk_ellipse

        # Determine the radius of the galaxy
        semimajor_angular = truncation_ellipse.semimajor  # semimajor axis length of the sky ellipse
        radius_physical = (semimajor_angular * self.galaxy_distance).to("pc", equivalencies=dimensionless_angles())

        # Set properties
        generator.grid_type = "bintree"  # set grid type
        generator.x_radius = radius_physical
        generator.y_radius = radius_physical
        generator.z_radius = 3. * u("kpc")

        # Set input
        input_dict = dict()
        input_dict["ngrids"] = 1
        input_dict["scale"] = 0.5 * self.deprojections["dust"].pixelscale # in pc
        input_dict["level"] = 6
        input_dict["mass_fraction"] = 0.5e-6

        # Generate the grid
        generator.run(**input_dict)

        # Set the dust grid
        self.dust_grid = generator.single_grid

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
        # galaxy_center, distance, pa, inclination, filepath, scale_height
        self.deprojections["old"] = DeprojectionModel3D.from_wcs(wcs, self.galaxy_center, self.galaxy_distance, self.galaxy_position_angle, self.galaxy_inclination, old_filename, old_scale_height)

        # Set deprojection of young stars
        wcs = self.maps["young"].wcs
        self.deprojections["young"] = DeprojectionModel3D.from_wcs(wcs, self.galaxy_center, self.galaxy_distance, self.galaxy_position_angle, self.galaxy_inclination, young_filename, young_scale_height)

        # Set deprojection of ionizing stars
        wcs = self.maps["ionizing"].wcs
        self.deprojections["ionizing"] = DeprojectionModel3D.from_wcs(wcs, self.galaxy_center, self.galaxy_distance, self.galaxy_position_angle, self.galaxy_inclination, ionizing_filename, ionizing_scale_height)

        # Set deprojection of dust map
        wcs = self.maps["dust"].wcs
        self.deprojections["dust"] = DeprojectionModel3D.from_wcs(wcs, self.galaxy_center, self.galaxy_distance, self.galaxy_position_angle, self.galaxy_inclination, dust_filename, dust_scale_height)

    # -----------------------------------------------------------------

    def create_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the ski file ...")

        # Load the ski file template
        self.ski = SkiFile(pan_ski_path)

        # Set components
        self.set_components()

        # Add the instrument
        self.ski.add_instrument("earth", self.instrument)

        # Set the number of photon packages
        self.ski.setpackages(self.config.npackages)

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid(fs.name(self.wavelength_grid_path))

        # Set the dust emissivity
        #self.set_dust_emissivity()
        # Set dust emissivity if present
        self.ski.set_transient_dust_emissivity()
        #if self.ski.has_dust_emissivity:
            # Enable or disable
            #if self.config.transient_heating:
            #    self.ski.set_transient_dust_emissivity()
            #else:
            #    self.ski.set_grey_body_dust_emissivity()

        # Set the lowest-resolution dust grid
        self.ski.set_dust_grid(self.dust_grid)

        # Set all-cells dust library
        self.ski.set_allcells_dust_lib()

        # Set the dust selfabsorption
        #self.set_selfabsorption()
        # Dust self-absorption
        #if self.config.selfabsorption:
        #    self.ski.enable_selfabsorption()
        #else:
        #    self.ski.disable_selfabsorption()
        self.ski.enable_selfabsorption()

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
        title = titles["young"]

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
        # component_id, hydrocarbon_pops=25, enstatite_pops=25, forsterite_pops=25, write_mix=True, write_mean_mix=True, write_size=True
        self.ski.set_dust_component_themis_mix(title, hydrocarbon_pops, enstatite_pops, forsterite_pops, write_mix=False, write_mean_mix=False, write_size=False)

        # Set the normalization
        self.ski.set_dust_component_mass(dust_mass)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write the ski file
        self.write_ski()

        # Write the input
        self.write_input()

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Save
        self.ski.saveto(self.ski_path)

    # -----------------------------------------------------------------

    def write_input(self):

        """
        This function ...
        :return:
        """

        # Write wavelength grid
        self.wavelength_grid.saveto(self.wavelength_grid_path)

        # Write maps
        self.maps["old"].saveto(fs.join(self.simulation_input_path, old_filename))
        self.maps["young"].saveto(fs.join(self.simulation_input_path, young_filename))
        self.maps["ionizing"].saveto(fs.join(self.simulation_input_path, ionizing_filename))
        self.maps["dust"].saveto(fs.join(self.simulation_input_path, dust_filename))

    # -----------------------------------------------------------------

    def launch_reference(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the reference simulation ...")

        # Settings
        settings_launch = dict()
        settings_launch["ski"] = self.ski_path
        settings_launch["input"] = self.simulation_input_path
        settings_launch["output"] = self.simulation_output_path
        settings_launch["create_output"] = True

        # Input
        input_launch = dict()

        # Launch command
        launch = Command("launch_simulation", "launch the reference simulation", settings_launch, input_launch, cwd=".")
        self.launcher = self.run_command(launch)

    # -----------------------------------------------------------------

    def make_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the observed mock SED ...")

        # Settings
        settings_sed = dict()
        settings_sed["spectral_convolution"] = False

        # Input
        input_sed = dict()
        input_sed["simulation_output_path"] = self.simulation_output_path
        input_sed["output_path"] = "."

        # Launch command
        calculate = Command("observed_fluxes", "create the mock SED", settings_sed, input_sed, cwd=".")
        self.flux_calculator = self.run_command(calculate)

        # Determine the path to the mock SED
        #mock_sed_path = "spiral_earth_fluxes.dat"

    # -----------------------------------------------------------------

    def make_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the observed mock images ...")

        # Create Aniano kernels object
        aniano = AnianoKernels()

        # Set the paths to the kernel for each image
        kernel_paths = dict()
        for filter_name in fitting_filter_names: kernel_paths[filter_name] = aniano.get_psf_path(
            parse_filter(filter_name))

        # Settings
        settings_images = dict()
        settings_images["spectral_convolution"] = False
        # No output path is specified, so images won't be written out

        # Input
        input_images = dict()
        input_images["simulation_output_path"] = self.simulation_output_path
        input_images["output_path"] = "."
        input_images["filter_names"] = fitting_filter_names
        input_images["instrument_names"] = ["earth"]
        # input_images["wcs_path"] =
        input_images["wcs"] = self.wcs
        input_images["kernel_paths"] = kernel_paths
        input_images["unit"] = "Jy/pix"
        # input_images["host_id"] = "nancy"

        # Launch the command
        make = Command("observed_images", "create the mock images", settings_images, input_images, cwd=".")
        self.image_maker = self.run_command(make)

        # MAKE DATA:

        # Create directory for the images
        ref_path = fs.create_directory_in(self.path, "ref")
        images_path = fs.create_directory_in(ref_path, "images")

        # Determine the name of the datacube
        datacube_names = self.image_maker.images.keys()
        if len(datacube_names) > 1: raise RuntimeError("Unexpected number of datacubes")
        datacube_name = datacube_names[0]

        # Loop over the images
        for filter_name in self.image_maker.images[datacube_name]:

            # Get the image
            image = self.image_maker.images[datacube_name][filter_name]

            #

            # Save the image
            image_path = fs.join(images_path, filter_name + ".fits")
            image.saveto(image_path)

    # -----------------------------------------------------------------

    def setup_modelling(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting up the modelling ...")

        # Settings
        settings_setup = dict()
        settings_setup["type"] = "galaxy"
        settings_setup["name"] = "Galaxy"
        settings_setup["fitting_host_ids"] = None

        # Create input dict for setup
        input_setup = dict()
        input_setup["ngc_name"] = fake_name
        input_setup["hyperleda_name"] = fake_name

        # Construct the command
        # stp = Command("setup", "setup the modeling", settings_setup, input_setup, cwd=".")

        # Add the command
        # commands.append(stp)

    # -----------------------------------------------------------------

    def model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the modelling ...")

        # Settings
        settings_model = dict()
        settings_model["ngenerations"] = 4
        settings_model["nsimulations"] = 20
        settings_model["fitting_settings"] = {"spectral_convolution": False}

        # Input
        input_model = dict()

        # Set galaxy properties
        input_model["properties"] = self.properties

        # Set SEDs
        input_model["seds"] = dict()

        # Set images dictionary
        images = dict()
        for filter_name in fitting_filter_names:
            images[filter_name] = fs.join("../ref/images", filter_name + ".fits")
        input_model["images"] = images

        # Construct the command
        # command = Command("model", "perform the modelling", settings_model, input_model, "./Galaxy")

        # Add the command
        # commands.append(command)

# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    :param temp_path:
    :return:
    """

    pass

# -----------------------------------------------------------------
