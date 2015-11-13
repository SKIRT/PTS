#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import inspect
import subprocess
from datetime import datetime

# Import Astromagic modules
from ..core.frames import Frame
from ..tools import configuration

# Import astronomical modules
from astropy.table import Table
from astropy import log
import astropy.logger

# *****************************************************************

class SExtractor(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        ## Configuration

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "sextractor.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config = configuration.open(config, default_config)

        ## Logging

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

        ## Attributes

        # Set the frame to None initially
        self.frame = None

        # Set the path to the dat/SExtractor directory
        self.dat_path = os.path.join(directory, "dat", "sextractor")

        # Set the path to the temporary directory to None initially
        self.directory = None

        # Set the segmentation map to None initially
        self.segments = None

        # Set the catalog to None initially
        self.catalog = None

    # *****************************************************************

    def run(self, frame, galaxyextractor=None, starextractor=None):

        """
        This function ...
        :return:
        """

        # Make a local reference to the frame
        self.frame = frame

        # Create a temporary directory for the input and output of SExtractor
        self.create_directory()

        # Copy the input files into the temporary directory
        self.copy_input()

        # Launch SExtractor
        self.launch_sextractor()

        # Read the SExtractor output
        self.read_output()

        # Remove the temporary directory, if requested
        if self.config.remove_temp: self.remove_directory()

    # *****************************************************************

    def create_directory(self):

        """
        This function ...
        :return:
        """

        # Generate a timestamp identifying this particular run for the scaling test
        timestamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

        # Set the path to the temporary directory
        self.directory = os.path.join(os.getcwd(), "sextractor_" + timestamp)

        # Create the directory
        os.mkdir(self.directory)

    # *****************************************************************

    def copy_input(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the SExtractor input file and open it
        input_file_path = os.path.join(self.dat_path, self.config.input_file)
        input_file = open(input_file_path, 'r')

        # Create an empty file for the
        self.new_input_file_path = os.path.join(self.directory, self.config.input_file)
        new_input_file = open(self.new_input_file_path, 'w')

        for line in input_file:

            # Look for the name of the parameters file
            if "PARAMETERS_NAME" in line or "FILTER_NAME" in line or "STARNNW_NAME" in line:

                name = line.split()[1]
                new_line = line.replace(name, os.path.join(self.dat_path, name))

            # Just copy all other lines
            else: new_line = line

            # Write the line to the new file
            new_input_file.write(new_line)

        # Close both the original input file and the temporary copy
        input_file.close()
        new_input_file.close()

        # Save the input frame to the temporary directory
        image_path = os.path.join(self.directory, "input.fits")
        self.frame.save(image_path)

    # *****************************************************************

    def launch_sextractor(self):

        """
        This function ...
        :param self:
        :param file_in:
        :param m0:
        :param GAIN:
        :param pix2sec:
        :param fwhm:
        :param se_file:
        :return:
        """

        # Switch to the temporary directory
        os.chdir(self.directory)

        # Check how to call SExtractor on this system
        if subprocess.call("which sex >/dev/null", shell=True) == 0:

            command = "sex input.fits"

        elif subprocess.call("which sextractor >/dev/null", shell=True) == 0:

            command = "sextractor input.fits"

        # If none of the above commands worked, raise an error
        else: raise RuntimeError("SExtractor was not found on this system.")

        # Add command-line arguments for SExtractor
        command += " -MAG_ZEROPOINT " + str(self.config.zero_point)
        command += " -GAIN " + str(self.config.gain)
        command += " -PIXEL_SCALE " + str(self.config.pixelscale)
        command += " -SEEING_FWHM " + str(self.config.fwhm)
        command += " -c " + self.config.input_file
        command += " -VERBOSE_TYPE=QUIET"

        # Inform the user
        log.debug("SExtractor command: " + command)

        # Inform the user
        log.info("Running SExtractor...")

        # Launch the SExtractor command as a seperate process
        subprocess.call(command, shell=True)

    # *****************************************************************

    def read_output(self):

        """
        This function ...
        :return:
        """

        # Read in the segmentation map
        segments_path = os.path.join(self.directory, "segm.fits")
        self.segments = Frame.from_file(segments_path)

        # Read in the catalog file
        catalog_path = os.path.join(self.directory, "field.cat")
        self.catalog = Table.read(catalog_path, format="ascii")

    # *****************************************************************

    def remove_directory(self):

        """
        This function ...
        :return:
        """

        pass

    # *****************************************************************