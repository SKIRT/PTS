#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.analysis.sextractor Contains the SExtractor class

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import subprocess

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ..core.frame import Frame
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.basics.log import log

# -----------------------------------------------------------------

class SExtractor(object):
    
    """
    This class ...
    """
    
    def __init__(self):
        
        """
        The constructor ...
        """

        # Set the path to the dat/sextractor directory
        self.data_path = os.path.join(introspection.pts_dat_dir("magic"), "sextractor")

    # -----------------------------------------------------------------

    def run(self, frame, input_file_name="default.sex", zero_point=20.0, gain=4.0, pixelscale=1.0, fwhm=1.0, keep=False):

        """
        This function ...
        :param frame:
        :param input_file_name: sextractor input file
        :param zero_point: zero point in [mag arcsec^-2]
        :param gain: gain in [e-/ADU]
        :param pixelscale: pixelscale in [arcsec/pix]
        :param fwhm: FWHM of the PSF
        :param keep: keep the temporary directory
        :return:
        """

        # Create a temporary directory for the input and output of SExtractor
        temp_path = fs.create_temporary_directory("sextractor")

        # Copy the input files into the temporary directory
        self.copy_input(temp_path, frame, input_file_name)

        # Launch SExtractor
        self.launch(temp_path, input_file_name, zero_point, gain, pixelscale, fwhm)

        # Read the SExtractor output
        segments, catalog = self.read_output(temp_path)

        # Remove the temporary directory, if requested
        if not keep: fs.remove_directory(temp_path)

        # Return the segmentation map and catalog
        return segments, catalog

    # -----------------------------------------------------------------

    def copy_input(self, directory, frame, input_file_name):

        """
        This function ...
        :param directory:
        :param frame:
        :param input_file_name:
        :return:
        """

        # Determine the path to the SExtractor input file and open it
        input_file_path = os.path.join(self.data_path, input_file_name)
        input_file = open(input_file_path, 'r')

        # Create an empty file for the
        new_input_file_path = os.path.join(directory, input_file)
        new_input_file = open(new_input_file_path, 'w')

        for line in input_file:

            # Look for the name of the parameters file
            if "PARAMETERS_NAME" in line or "FILTER_NAME" in line or "STARNNW_NAME" in line:

                name = line.split()[1]
                new_line = line.replace(name, os.path.join(self.data_path, name))

            # Just copy all other lines
            else: new_line = line

            # Write the line to the new file
            new_input_file.write(new_line)

        # Close both the original input file and the temporary copy
        input_file.close()
        new_input_file.close()

        # Save the input frame to the temporary directory
        image_path = os.path.join(directory, "input.fits")
        frame.saveto(image_path)

    # -----------------------------------------------------------------

    def launch(self, directory, input_file_name, zero_point, gain, pixelscale, fwhm):

        """
        This function ...
        :param directory:
        :param input_file_name:
        :param zero_point:
        :param gain:
        :param pixelscale:
        :param fwhm:
        :return:
        """

        # Switch to the temporary directory
        os.chdir(directory)

        # Check how to call SExtractor on this system
        if subprocess.call("which sex >/dev/null", shell=True) == 0: command = "sex input.fits"
        elif subprocess.call("which sextractor >/dev/null", shell=True) == 0: command = "sextractor input.fits"
        else: raise RuntimeError("SExtractor was not found on this system.") # If none of the above commands worked, raise an error

        # Add command-line arguments for SExtractor
        command += " -MAG_ZEROPOINT " + str(zero_point)
        command += " -GAIN " + str(gain)
        command += " -PIXEL_SCALE " + str(pixelscale)
        command += " -SEEING_FWHM " + str(fwhm)
        command += " -c " + input_file_name
        command += " -VERBOSE_TYPE=QUIET"

        # Inform the user
        log.debug("SExtractor command: " + command)

        # Inform the user
        log.info("Running SExtractor...")

        # Launch the SExtractor command as a seperate process
        subprocess.call(command, shell=True)

    # -----------------------------------------------------------------

    def read_output(self, directory):

        """
        This function ...
        :param directory:
        :return:
        """

        # Read in the segmentation map
        segments_path = os.path.join(directory, "segm.fits")
        segments = Frame.from_file(segments_path)

        # Read in the catalog file
        catalog_path = os.path.join(directory, "field.cat")
        catalog = Table.read(catalog_path, format="ascii")

        # Return the segments and catalog
        return segments, catalog

# -----------------------------------------------------------------
