#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.prepare.preparer Contains the ImagePreparer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.basics.configurable import InteractiveConfigurable
from ..core.list import NamedImageList
from ...core.tools import filesystem as fs
from ...core.basics.table import SmartTable
from ...core.tools.utils import lazyproperty
from ...core.basics.configuration import ConfigurationDefinition
from ..config.find_sources import definition as find_sources_definition
from ..config.extract import definition as extract_sources_definition

# -----------------------------------------------------------------

_help_command_name = "help"
_history_command_name = "history"
_status_command_name = "status"

_initialize_command_name = "initialize"
_find_command_name = "find"
_extract_command_name = "extract"
_correct_command_name = "correct"
_convolve_command_name = "convolve"
_rebin_command_name = "rebin"
_subtract_command_name = "subtract"
_errors_command_name = "errors"
_units_command_name = "units"

# -----------------------------------------------------------------

commands = OrderedDict()

# Standard commands
commands[_help_command_name] = ("show_help", False, "show help", None)
commands[_history_command_name] = ("show_history_command", True, "show history of executed commands", None)
commands[_status_command_name] = ("show_status_command", True, "show preparation status", None)

# Preparation of image commands
commands[_initialize_command_name] = ("initialize_image", True, "initialize an image", "image")
commands[_find_command_name] = ("find_sources_image", True, "find sources in an image", "image")
commands[_extract_command_name] = ("extract_sources_image", True, "extract sources in an image", "image")
commands[_correct_command_name] = ("correct_extinction_image", True, "correct an image for galactic extinction", "image")
commands[_convolve_command_name] = ("convolve_image", True, "convolve an image", "image")
commands[_rebin_command_name] = ("rebin_image", True, "rebin an image", "image")
commands[_subtract_command_name] = ("subtract_background_image", True, "subtract the background from an image", "image")
commands[_errors_command_name] = ("create_errors_image", True, "create errormap for an image", "image")
commands[_units_command_name] = ("convert_units_image", True, "convert the units of an image", "image")

# -----------------------------------------------------------------

subcommands = OrderedDict()

# -----------------------------------------------------------------

class PreparationStatusTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Image name"] = (str, None, "name of the image")
    _column_info["Status"] = (str, None, "preparation status")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PreparationStatusTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_row(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        raise RuntimeError("Cannot add rows to a preparation status table")

    # -----------------------------------------------------------------

    def is_initialized(self, name):

        """
        This function ...
        :param name:
        :return:
        """

    # -----------------------------------------------------------------

    def is_source_found(self, name):

        """
        This function ...
        :param name:
        :return:
        """

    # -----------------------------------------------------------------

    def is_source_extracted(self, name):

        """
        This function ...
        :param name:
        :return:
        """

    # -----------------------------------------------------------------

    def is_extinction_corrected(self, name):

        """
        This function ...
        :param name:
        :return:
        """

    # -----------------------------------------------------------------

    def is_convolved(self, name):

        """
        This function ...
        :param name:
        :return:
        """

    # -----------------------------------------------------------------

    def is_rebinned(self, name):

        """
        This function ...
        :param name:
        :return:
        """

    # -----------------------------------------------------------------

    def is_background_subtracted(self, name):

        """
        This function ...
        :param name:
        :return:
        """

    # -----------------------------------------------------------------

    def is_errors_created(self, name):

        """
        This function ...
        :param name:
        :return:
        """

    # -----------------------------------------------------------------

    def is_unit_converted(self, name):

        """
        This function ...
        :param name:
        :return:
        """

# -----------------------------------------------------------------

_fwhm_extra_name = "FWHM"
_size_extra_name = "size"
_xsize_extra_name = "xsize"
_ysize_extra_name = "ysize"

# -----------------------------------------------------------------

# Define extra columns
extra_columns = OrderedDict()
extra_columns[_fwhm_extra_name] = "screen session name"
extra_columns[_size_extra_name] = "simulation output disk size"
extra_columns[_xsize_extra_name] = "number of x pixels"
extra_columns[_ysize_extra_name] = "number of y pixels"

# Define extra column names
extra_column_names = dict()
extra_column_names[_fwhm_extra_name] = "FWHM of the PSF"
extra_column_names[_size_extra_name] = "Filesize"
extra_column_names[_xsize_extra_name] = "Nx"
extra_column_names[_ysize_extra_name] = "Ny"

# Define extra column units
extra_column_units = dict()
extra_column_units[_fwhm_extra_name] = "arcsec"
extra_column_units[_size_extra_name] = "GB"

# -----------------------------------------------------------------

class ImagePreparer(InteractiveConfigurable):
    
    """
    This class ...
    """

    _commands = commands
    _subcommands = subcommands

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):
        
        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ImagePreparer, self).__init__(*args, **kwargs)

        # The preparation status
        self.status = None

        # The images
        self.images = NamedImageList()

    # -----------------------------------------------------------------

    @property
    def do_commands(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_interactive(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_initialization(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_finding(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_extraction(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_extinction(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_convolution(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_rebinning(self):

        """
        This fucntion ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_subtraction(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_errors(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_conversion(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_showing(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_writing(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_plotting(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Run commands
        if self.do_commands: self.run_commands()

        # 2. Interactive mode
        if self.do_interactive: self.interactive()

        # 3. Initialize the images
        if self.do_initialization: self.initialize()

        # 4. Find sources
        if self.do_finding: self.find_sources()

        # 2. Extract sources
        if self.do_extraction: self.extract_sources()

        # 3. Correct for galactic extinction
        if self.do_extinction: self.correct_for_extinction()

        # Convolve
        if self.do_convolution: self.convolve()

        # Rebin
        if self.do_rebinning: self.rebin()

        # 4. Subtract background
        if self.do_subtraction: self.subtract_background()

        # 5. Calculate the error maps
        if self.do_errors: self.create_errormaps()

        # 6. If requested, convert the unit
        if self.do_conversion: self.convert_units()

        # Show
        if self.do_showing: self.show()

        # 8. Writing
        if self.do_writing: self.write()
        
        # Plot
        if self.do_plotting: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ImagePreparer, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        This function ...
        :return:
        """

        return self.images.names

    # -----------------------------------------------------------------

    def get_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.images[name]

    # -----------------------------------------------------------------

    def get_filter(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_image(name).filter

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return:
        """

        return [self.get_filter(name) for name in self.names]

    # -----------------------------------------------------------------

    def is_initialized(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.status.is_initialized(name)

    # -----------------------------------------------------------------

    def is_source_found(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.status.is_source_found(name)

    # -----------------------------------------------------------------

    def is_source_extracted(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.status.is_source_extracted(name)

    # -----------------------------------------------------------------

    def is_extinction_corrected(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.status.is_extinction_corrected(name)

    # -----------------------------------------------------------------

    def is_convolved(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.status.is_convolved(name)

    # -----------------------------------------------------------------

    def is_rebinned(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.status.is_rebinned(name)

    # -----------------------------------------------------------------

    def is_background_subtracted(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.status.is_background_subtracted(name)

    # -----------------------------------------------------------------

    def is_errors_created(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.status.is_errors_created(name)

    # -----------------------------------------------------------------

    def is_unit_converted(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.status.is_unit_converted(name)

    # -----------------------------------------------------------------

    @property
    def not_initialized_names(self):

        """
        This function ...
        :return:
        """

        for name in self.names:
            if self.is_initialized(name): continue
            yield name

    @property
    def not_

    # -----------------------------------------------------------------

    @lazyproperty
    def show_status_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("extra_columns", "string_list", "extra columns to show", choices=extra_columns)
        definition.add_optional("path", "string", "save the status information as a table at this path")
        definition.add_flag("refresh", "refresh the status info", False)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def show_status_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.show_status_definition, **kwargs)

        # Show
        self.show_status(extra=config.extra_columns, path=config.path, refresh=config.refresh)

    # -----------------------------------------------------------------

    def show_status(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Showing the status of the preparation ...")

        # Get settings
        extra = kwargs.pop("extra", None)
        path = kwargs.pop("path", None)
        refresh = kwargs.pop("refresh", False)

        # Refresh if requested
        if refresh: self.reset_status()

        # ...

    # -----------------------------------------------------------------

    def initialize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the images ...")

        # Loop over all images
        for name in self.not_initialized_names:

            # Initialize
            self.initialize_image(name)

    # -----------------------------------------------------------------

    @lazyproperty
    def initialize_image_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def initialize_image_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get image name
        image_name, config = self.get_image_name_and_config_from_command(command, self.initialize_image_definition, **kwargs)

        # Check whether the simulation is not initialized yet
        if self.is_initialized(image_name): raise RuntimeError("Image '" + image_name + "' is already initialized")

        # Relaunch the simulation
        self.initialize_image(image_name)

    # -----------------------------------------------------------------

    def initialize_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Debugging
        log.debug("Initializing image '" + name + "' ...")

        # Get the image path
        #image_path = self.paths[prep_name]

        # Determine error path
        #error_path = self.error_paths[prep_name] if prep_name in self.error_paths else None

        # Determine the output path for this image
        #output_path = self.get_prep_path(prep_name)

        # Set initialized path
        initialized_path = fs.join(output_path, "initialized.fits")
        self.initialized_paths[prep_name] = initialized_path

        # Check whether this image already has an initialized image
        if fs.is_file(initialized_path):

            # Already present
            log.success("Initialized '" + prep_name + "' is already present")

            # Check wether the sources directory is present
            # -> NO, we are going to loop over self.paths again anyway in the create_directories and get_sources functions ...

            # Check whether in dataset
            if prep_name not in self.set.names:

                # Give a warning
                log.warning("Initialized '" + prep_name + "' was not yet in the dataset: adding it now ...")
                # Add to the dataset
                # Add entry to the dataset
                self.set.add_path(prep_name, initialized_path)
                self.set.save() # Save

            # Cache
            self.cache_image(prep_name, image_path, error_path)

            # Now skip the rest
            continue

        # Debugging
        log.debug("Initializing image '" + image_path + "' ...")

        # Set the path to the region of bad pixels
        bad_region_path = fs.join(self.data_path, "bad", prep_name + ".reg")
        if not fs.is_file(bad_region_path): bad_region_path = None

        # Get the filter
        fltr = parse_filter(prep_name)

        # Set the FWHM if the instrument has a fixed PSF
        if has_variable_fwhm(fltr): fwhm = None
        else: fwhm = get_fwhm(fltr)

        # Debugging
        log.debug("Loading image " + image_path + " as " + prep_name + " ...")

        # Import the image
        importer = ImageImporter()
        importer.run(image_path, bad_region_path, fwhm=fwhm, find_error_frame=False) # don't look for error frames

        # Get the imported image
        image = importer.image

        # Set the image name
        image.name = prep_name

        # Remove all frames except for the primary frame
        image.remove_frames_except("primary")

        # If a poisson error map was found, add it to the image
        if prep_name in self.error_paths:

            # Debugging
            log.debug("Adding the poisson error frame to the " + prep_name + " image ...")

            # Add the error frame
            error_map = Frame.from_file(self.error_paths[prep_name])
            image.add_frame(error_map, "errors")

        # Save the image
        image.saveto(initialized_path)

        # Success
        log.success("Initialized the '" + prep_name + "' image")

        # NEW: ADD TO THE DATASET NOW AND SAVE IT IMMEDIATELY

        # Inform the user
        log.info("Adding the '" + prep_name + "' image to the dataset ...")

        # Add entry to the dataset
        self.set.add_path(prep_name, initialized_path)

        # Cache the original image
        self.cache_image(prep_name, image_path, error_path)

    # -----------------------------------------------------------------

    def find_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding sources in the images ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def find_sources_image_definition(self):

        """
        This function ...
        :return:
        """

        return find_sources_definition.copy(pos_optional=False)

    # -----------------------------------------------------------------

    def find_sources_image_command(self, command, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get image name
        image_name, config = self.get_image_name_and_config_from_command(command, self.find_sources_image_definition, **kwargs)

        # Check
        if self.is_source_found(image_name): raise RuntimeError("Sources are already found for image '" + image_name + "'")

        # Find sources
        self.find_sources_image(image_name)

    # -----------------------------------------------------------------

    def find_sources_image(self, name):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Finding sources in the '" + name + "' image ...")

    # -----------------------------------------------------------------

    def extract_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting sources in the images ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def extract_sources_image_definition(self):

        """
        This property ...
        :return:
        """

        return extract_sources_definition.copy(pos_optional=False)

    # -----------------------------------------------------------------

    def extract_sources_image_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get image name
        image_name, config = self.get_image_name_and_config_from_command(command, self.extract_sources_image_definition, **kwargs)

        # Check
        if self.is_source_extracted(image_name): raise RuntimeError("Sources are already extracted for image '" + image_name + "'")

        # Extract sources
        self.extract_sources_image(image_name)

    # -----------------------------------------------------------------

    def extract_sources_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Debugging
        log.debug("Extracting sources in the '" + name + "' image ...")

    # -----------------------------------------------------------------

    def correct_for_extinction(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting images for galactic extinction ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def correct_extinction_image_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def correct_extinction_image_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get image name
        image_name, config = self.get_image_name_and_config_from_command(command, self.correct_extinction_image_definition, **kwargs)

        # Check
        if self.is_extinction_corrected(image_name): raise RuntimeError("Image '" + image_name + "' is already corrected for galactic extinction")

        # Correct
        self.correct_extinction_image(image_name)

    # -----------------------------------------------------------------

    def correct_extinction_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Debugging
        log.debug("Correcting image '" + name + "' for galactic extinction ...")

    # -----------------------------------------------------------------

    def convolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Convolving the images ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def convolve_image_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def convolve_image_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get image name
        image_name, config = self.get_image_name_and_config_from_command(command, self.convolve_image_definition, **kwargs)

        # Check
        if self.is_convolved(image_name): raise RuntimeError("Image '" + image_name + "' is already convolved")

        # Convolve
        self.convolve_image(image_name)

    # -----------------------------------------------------------------

    def convolve_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Debugging
        log.debug("Convolving the '" + name + "' image ...")

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning the images ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def rebin_image_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def rebin_image_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    def rebin_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Debugging
        log.debug("Rebinning the '" + name + "' image ...")

    # -----------------------------------------------------------------

    def subtract_background(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the background from the images ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def subtract_background_image_definition(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def subtract_background_image_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    def subtract_background_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """


    # -----------------------------------------------------------------

    def create_errormaps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating errormaps for the images ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def create_errors_image_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def create_errors_image_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    def create_errors_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Debugging
        log.debug("Creating errormap for the '" + name + "' image ...")

    # -----------------------------------------------------------------

    def convert_units(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Converting the units of the images ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def convert_units_image_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("unit", "unit", "target unit")
        return definition

    # -----------------------------------------------------------------

    def convert_units_image_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

    # -----------------------------------------------------------------

    def convert_units_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Debugging
        log.debug("Converting the units of the '" + name + "' image ...")

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
