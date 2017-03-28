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

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.tools import network
from pts.core.tools.logging import log
from pts.do.commandline import Command
from pts.core.test.implementation import TestImplementation
from pts.core.tools import introspection
from pts.core.basics.range import RealRange, QuantityRange
from pts.magic.core.frame import Frame
from pts.core.basics.quantity import parse_quantity, parse_angle
from pts.core.simulation.skifile import LabeledSkiFile
from pts.modeling.basics.instruments import InstrumentFrame, MultiFrameInstrument
from pts.core.filter.filter import parse_filter
from pts.modeling.fitting.fitskirt import FskiFile
from pts.core.simulation.grids import CylindricalGrid

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "Fitting the galaxy NGC4013 using a flattened Sersic profile for the central bulge and a double exponential for the stellar disk and dust disk"

# -----------------------------------------------------------------

input_url = "http://www.skirt.ugent.be/downloads/tutorial_NGC4013.tar.gz"
all_url = "http://www.skirt.ugent.be/downloads/tutorial_NGC4013_complete.tar.gz"

# -----------------------------------------------------------------

ski_path = fs.join(this_dir_path, "ngc4013.ski")
fski_path = fs.join(this_dir_path, "ngc4013.fski")

# -----------------------------------------------------------------

# Determine the path to the dropbox path and the path of the directory with the data for M81
ngc4013_data_path = fs.join(introspection.get_dropbox_tests_pts_path_for_subproject("modeling"), "NGC4013")

# Determine path to 'norm' directory
norm_path = fs.join(ngc4013_data_path, "norm")

# -----------------------------------------------------------------

# Free parameters
free_parameter_labels = ["inclination", "stellar_length", "stellar_height", "flattening", "sersic_index", "bulge_radius", "dust_length", "dust_height", "dust_mass"]

# Free parameters with descriptions
free_parameters = dict()
free_parameters["inclination"] = "inclination of the galaxy"
free_parameters["stellar_length"] = "scale length of the stellar disk"
free_parameters["stellar_height"] = "scale height of the stellar disk"
free_parameters["flattening"] = "flattening of the stellar bulge"
free_parameters["sersic_index"] = "sersic index of the stellar bulge"
free_parameters["bulge_radius"] = "scale length of the stellar bulge"
free_parameters["dust_length"] = "scale length of the dust disk"
free_parameters["dust_height"] = "scale height of the dust disk"
free_parameters["dust_mass"] = "total dust mass"

# Fitting filters
fitting_filter_names = ["SDSS u", "SDSS g"]

# Types of parameters
parameter_skirt_types = dict()
parameter_skirt_types["inclination"] = "posangle"
parameter_skirt_types["stellar_length"] = "length"
parameter_skirt_types["stellar_height"] = "length"
parameter_skirt_types["flattening"] = "dimless"
parameter_skirt_types["sersic_index"] = "length"
parameter_skirt_types["bulge_radius"] = "length"
parameter_skirt_types["dust_length"] = "length"
parameter_skirt_types["dust_height"] = "length"
parameter_skirt_types["dust_mass"] = "mass"

# Ranges of parameters
parameter_ranges = dict()
parameter_ranges["inclination"] = QuantityRange("88 deg", "92 deg")
parameter_ranges["stellar_length"] = QuantityRange("500 pc", "8000 pc")
parameter_ranges["stellar_height"] = QuantityRange("100 pc", "1000 pc")
parameter_ranges["flattening"] = RealRange(0.01, 1)
parameter_ranges["sersic_index"] = RealRange(0.51, 6.99)
parameter_ranges["bulge_radius"] = QuantityRange("200 pc", "5000 pc")
parameter_ranges["dust_length"] = QuantityRange("1000 pc", "12000 pc")
parameter_ranges["dust_height"] = QuantityRange("50 pc", "800 pc")
parameter_ranges["dust_mass"] = QuantityRange("1e6 Msun", "1e8 Msun")

# Initial guesses
initial_guesses = dict()
initial_guesses["inclination"] = parse_quantity("90 deg")
initial_guesses["stellar_length"] = parse_quantity("4400 pc")
initial_guesses["stellar_height"] = parse_quantity("500 pc")
initial_guesses["flattening"] = 0.5
initial_guesses["sersic_index"] = 2.5
initial_guesses["bulge_radius"] = parse_quantity("2500 pc")
initial_guesses["dust_length"] = parse_quantity("6600 pc")
initial_guesses["dust_height"] = parse_quantity("250 pc")
initial_guesses["dust_mass"] = parse_quantity("4e7 Msun")

# -----------------------------------------------------------------

instrument_name = "earth"

# -----------------------------------------------------------------

min_luminosity = 1e5
max_luminosity = 1e7

# -----------------------------------------------------------------

class NGC4013Test(TestImplementation):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(NGC4013Test, self).__init__(config, interactive)

        # The path with the test data
        self.data_path = None

        # The path for the fitskirt run
        self.reference_path = None
        self.reference_output_path = None
        self.reference_ski_path = None
        self.reference_fski_path = None

        # The ski file
        self.ski = None

        # The fski file
        self.fski = None

        # The images
        self.images = dict()

        # The original image paths
        self.image_paths = dict()

        # The instrument
        self.instrument = None

        # The reference wcs
        self.wcs = None

        # The wavelengths
        self.wavelengths = None

        # The dust grid
        self.dust_grid = None

        # The FitSKIRT launcher
        self.launcher = None

        # The modeler
        self.modeler = None

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

        # 3. Load the fski file
        self.load_fski()

        # 2. Load the images
        self.load_images()

        # 3. Create instrument
        self.create_instrument()

        # 7. Create the wavelength grid
        self.create_wavelength_grid()

        # 8. Create the dust grid
        self.create_dust_grid()

        # 9. Create the ski file
        self.adjust_ski()

        # Create the fski file
        self.adjust_fski()

        # Write
        self.write()

        # 3. Launch with FitSKIRT
        self.launch_fitskirt()

        # Setup the modeling
        self.setup_modelling()

        # Model
        self.model()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # CAll the setup function of the base class
        super(NGC4013Test, self).setup(**kwargs)

        # Check the data
        if fs.is_directory(norm_path): self.data_path = norm_path
        else:

            # Create the data directory and get the data
            self.data_path = fs.create_directory_in(self.path, "data")
            self.get_data()

        # Create the reference directory and subdirectories
        self.reference_path = fs.create_directory_in(self.path, "ref")
        self.reference_output_path = fs.create_directory_in(self.reference_path, "out")
        self.reference_ski_path = fs.join(self.reference_path, "NGC4013.ski")
        self.reference_fski_path = fs.join(self.reference_path, "NGC4013.fski")

    # -----------------------------------------------------------------

    def get_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Downloading the input ...")

        # Download the input
        network.download_and_decompress_directory(all_url, self.data_path, progress_bar=True, into_root=True)

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file ...")

        # Load ski
        self.ski = LabeledSkiFile(ski_path)

    # -----------------------------------------------------------------

    def load_fski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the fski file ...")

        # Load the fski file
        self.fski = FskiFile(fski_path)

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # The common wcs of the images
        wcs = None

        # Loop over the images in the data directory
        for path, filename in fs.files_in_path(self.data_path, extension="fits", returns=["path", "name"]):

            # Load the frame
            frame = Frame.from_file(path)

            # Set filter
            previous_filter = frame.filter
            frame.filter = parse_filter(filename.split("_norm")[0])
            if previous_filter != frame.filter: frame.save()

            # Check filter
            if str(frame.filter) not in [str(fltr) for fltr in self.config.fitting_filters]: continue

            # Determine name
            name = str(frame.filter)

            # Check wcs
            if wcs is None: wcs = frame.wcs
            elif wcs == frame.wcs: pass
            else: raise IOError("The coordinate system of image '" + filename + "' does not match that of other images")

            # Debugging
            log.debug("Adding frame '" + filename + "' ...")

            # Add to dictionary
            self.images[name] = frame

            # Set original path
            self.image_paths[name] = path

        # Set the wcs
        self.wcs = wcs

    # -----------------------------------------------------------------

    def create_instrument(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instrument ...")

        # Set general properties
        distance = parse_quantity("18.6 Mpc")
        inclination = initial_guesses["inclination"]
        azimuth = parse_angle("0 deg")
        position_angle = parse_angle("0 deg")

        # Create multiframe instrument
        self.instrument = MultiFrameInstrument(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=position_angle)

        # Loop over the filters, create frames
        for fltr in self.config.fitting_filters:

            # Determine xsize and ysize
            xsize = self.wcs.xsize
            ysize = self.wcs.ysize

            # Determine pixelscale
            pixelscale_physical_x = (self.wcs.pixelscale.x * distance).to("pc", equivalencies=dimensionless_angles())
            pixelscale_physical_y = (self.wcs.pixelscale.y * distance).to("pc", equivalencies=dimensionless_angles())

            # Determine field of view
            field_x = pixelscale_physical_x * xsize
            field_y = pixelscale_physical_y * ysize

            # Create frame
            frame = InstrumentFrame(pixels_x=xsize, pixels_y=ysize, field_x=field_x, field_y=field_y)

            # Add the frame
            self.instrument.add_frame(frame)

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Set the wavelengths
        self.wavelengths = [fltr.pivot for fltr in self.config.fitting_filters]

    # -----------------------------------------------------------------

    def create_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust grid ...")

        # Set boundaries
        max_r = parse_quantity("25000 pc")
        min_z = parse_quantity("-4000 pc")
        max_z = parse_quantity("4000 pc")

        # Create the grid
        self.dust_grid = CylindricalGrid(max_r=max_r, min_z=min_z, max_z=max_z, type_r="logarithmic", type_z="symmetric_power",
                                         nbins_r=250, nbins_z=250, central_bin_fraction_r=0.0004, ratio_z=50)

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file ...")

        # Set the number of photon packages
        self.ski.setpackages(self.config.npackages)

        # Set wavelength grid
        self.ski.set_wavelengths(*self.wavelengths)

        # Set dust grid
        self.ski.set_dust_grid(self.dust_grid)

        # Set instrument
        self.ski.remove_all_instruments()
        self.ski.add_instrument(instrument_name, self.instrument)

        # Add label
        element = self.ski.get_instrument(instrument_name)
        path = "inclination"
        if "inclination" in self.config.free_parameters: self.ski.add_label_to_path(path, "inclination", element=element)

    # -----------------------------------------------------------------

    def adjust_fski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the fski file ...")

        # Set ski name
        self.fski.set_ski_name(fs.name(self.reference_ski_path))

        # Remove reference images
        self.fski.remove_all_reference_images()

        # Set reference images
        for name in self.image_paths:

            # Get filename
            path = self.image_paths[name]
            filename = fs.name(path)

            # Set FWHM
            if "u" in filename: kernel_fwhm = 1.7
            elif "g" in filename: kernel_fwhm = 1.6
            else: raise NotImplementedError("Not implemented")

            # Ranges (for both stellar components)
            luminosity_range = RealRange(min_luminosity, max_luminosity)
            luminosity_ranges = [luminosity_range, luminosity_range]

            # Add reference image
            self.fski.add_reference_image(filename, luminosity_ranges, kernel_fwhm, self.config.kernel_type, self.config.kernel_dimension)

        # Set genetic algorithm properties
        self.fski.set_population_size(self.config.nmodels)
        self.fski.set_ngenerations(self.config.ngenerations)
        self.fski.set_mutation_rate(self.config.mutation_rate)
        self.fski.set_crossover_rate(self.config.crossover_rate)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write ski
        self.write_ski()

        # Write fski
        self.write_fski()

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        # Write
        self.ski.saveto(self.reference_ski_path)

    # -----------------------------------------------------------------

    def write_fski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the fski file ...")

        # Write the fski file
        self.fski.saveto(self.reference_fski_path)

    # -----------------------------------------------------------------

    def launch_fitskirt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching FitSKIRT ...")

        # Create configuration
        config = dict()

        # The path to the ski file and fski file
        config["ski"] = self.reference_ski_path
        config["fski"] = self.reference_fski_path

        config["input"] = self.data_path
        config["output"] = self.reference_output_path

        # Input
        input_dict = dict()

        # Construct the command
        command = Command("fitskirt", "run the reference fitting with FitSKIRT", config, input_dict, cwd=".")

        # Run the command
        self.launcher = self.run_command(command)

    # -----------------------------------------------------------------

    def setup_modelling(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting up the modeling ...")

        # Settings
        settings_setup = dict()
        settings_setup["type"] = "images"
        settings_setup["name"] = "NGC4013"
        settings_setup["fitting_host_ids"] = None

        # Create object config
        object_config = dict()
        ski_path = fs.join(this_dir_path, "NGC4013.ski")
        object_config["ski"] = ski_path

        object_config["images"] = self.image_paths.values()

        # Create input dict for setup
        input_setup = dict()
        input_setup["object_config"] = object_config

        # Construct the command
        setup_command = Command("setup", "setup the modelling", settings_setup, input_setup, ".")

        # Run the command
        tool = self.run_command(setup_command)

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
        settings_model["ngenerations"] = self.config.ngenerations

        # -----------------------------------------------------------------

        # Create input dict for model
        input_model = dict()
        #input_model["parameters_config"] = Configuration(free_parameters=free_parameter_names)
        #input_model["descriptions_config"] = Configuration(descriptions=descriptions)
        #input_model["types_config"] = Configuration(types=types)
        #input_model["units_config"] = Configuration(units=units)
        #input_model["ranges_config"] = Configuration(luminosity_range=luminosity_range, dustmass_range=dustmass_range, grainsize_range=grainsize_range, fsil_range=fsil_range)
        #input_model["filters_config"] = Configuration(filters=filter_names)

        # Fitting initializer config
        #input_model["initialize_config"] = Configuration(npackages=1e4)

        # Construct the command
        model_command = Command("model", "perform the modelling", settings_model, input_model, cwd="./NGC4013")

        # Run the command
        self.modeler = self.run_command(model_command)

# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    :param temp_path:
    """

    return

# -----------------------------------------------------------------
