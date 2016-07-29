#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.component Contains the ModelingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from abc import ABCMeta

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ..core.sed import ObservedSED
from ...core.basics.filter import Filter
from ...magic.core.dataset import DataSet
from ...magic.core.frame import Frame
from ...magic.basics.skyregion import SkyRegion
from ...core.tools import parsing
from ...core.basics.map import Map
from ...magic.basics.skygeometry import SkyCoordinate

# -----------------------------------------------------------------

class ModelingComponent(Configurable):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(ModelingComponent, self).__init__(config)

        # Attributes
        self.galaxy_name = None

        # Modeling directories
        self.data_path = None
        self.prep_path = None
        self.truncation_path = None
        self.phot_path = None
        self.maps_path = None
        self.components_path = None
        self.fit_path = None
        self.analysis_path = None
        self.reports_path = None
        self.visualisation_path = None
        self.plot_path = None
        self.log_path = None
        self.show_path = None

        # PTS directories
        self.kernels_path = None

        # Reference image
        self.reference_image = "Pacs red"

        # The path to the observed SEDs
        self.observed_sed_path = None
        self.observed_sed_dustpedia_path = None

        # The path to the free parameter file and the fitting filters file
        self.free_parameters_path = None
        self.fitting_filters_path = None

        # The path to the maps
        self.old_stellar_map_path = None
        self.young_stellar_map_path = None
        self.ionizing_stellar_map_path = None
        self.dust_map_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ModelingComponent, self).setup()

        # -- Attributes --

        # Get the name of the galaxy (the name of the base directory)
        self.galaxy_name = fs.name(self.config.path)

        # Get the full paths to the necessary subdirectories
        self.data_path = fs.join(self.config.path, "data")
        self.prep_path = fs.join(self.config.path, "prep")
        self.truncation_path = fs.join(self.config.path, "truncated")
        self.phot_path = fs.join(self.config.path, "phot")
        self.maps_path = fs.join(self.config.path, "maps")
        self.components_path = fs.join(self.config.path, "components")
        self.fit_path = fs.join(self.config.path, "fit")
        self.analysis_path = fs.join(self.config.path, "analysis")
        self.reports_path = fs.join(self.config.path, "reports")
        self.visualisation_path = fs.join(self.config.path, "visualisation")
        self.plot_path = fs.join(self.config.path, "plot")
        self.log_path = fs.join(self.config.path, "log")
        self.show_path = fs.join(self.config.path, "show")

        # Determine the path to the kernels user directory
        self.kernels_path = fs.join(introspection.pts_user_dir, "kernels")

        # Check whether the 'data' directory exists, otherwise exit with an error
        if fs.is_directory(self.data_path):

            # Create the prep path if it does not exist yet
            fs.create_directories(self.prep_path, self.truncation_path, self.maps_path, self.phot_path,
                                  self.maps_path, self.components_path, self.fit_path, self.analysis_path,
                                  self.reports_path, self.visualisation_path, self.plot_path, self.log_path,
                                  self.show_path)

        # Exit with an error
        else: raise ValueError("The current working directory is not a radiative transfer modeling directory (the data directory is missing)")

        # Set the path to the observed SED
        self.observed_sed_path = fs.join(self.phot_path, "fluxes.dat")

        # Set the path to the DustPedia observed SED
        self.observed_sed_dustpedia_path = fs.join(self.data_path, "fluxes.dat")

        # Set the path to the free parameter file and the fitting filters file
        self.free_parameters_path = fs.join(self.fit_path, "free_parameters.txt")
        self.fitting_filters_path = fs.join(self.fit_path, "fitting_filters.txt")

        # Set the paths to the input maps
        self.old_stellar_map_path = fs.join(self.maps_path, "old_stars.fits")
        self.young_stellar_map_path = fs.join(self.maps_path, "young_stars.fits")
        self.ionizing_stellar_map_path = fs.join(self.maps_path, "ionizing_stars.fits")
        self.dust_map_path = fs.join(self.maps_path, "dust.fits")

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed(self):

        """
        This function ...
        :return:
        """

        return ObservedSED.from_file(self.observed_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed.filters()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_names(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.observed_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filters(self):

        """
        This function ...
        :return:
        """

        return map(Filter.from_string, self.fitting_filter_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filter_names(self):

        """
        This function ...
        :return:
        """

        return list(np.loadtxt(self.fitting_filters_path, dtype=str))

    # -----------------------------------------------------------------

    @lazyproperty
    def free_parameter_labels(self):

        """
        This function ...
        :return:
        """

        return get_free_parameter_labels(self.free_parameters_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dataset(self):

        """
        This function ...
        :return:
        """

        # Initialize the dataset
        dataset = DataSet()

        # Loop over all FITS files found in the 'truncated' directory
        for path, name in fs.files_in_path(self.truncation_path, extension="fits", returns=["path", "name"]):

            # Ignore the bulge, disk and model images
            if name == "bulge" or name == "disk" or name == "model": continue

            # Ignore the H alpha image
            if "Halpha" in name: continue

            # Add the image path to the dataset
            dataset.add_path(name, path)

        # Return the dataset
        return dataset

    # -----------------------------------------------------------------

    @lazyproperty
    def halpha_frame(self):

        """
        This function ...
        :return:
        """

        # Determine the path
        path = fs.join(self.truncation_path, "Mosaic Halpha.fits")

        # Load and return the frame
        return Frame.from_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_frame(self):

        """
        This function ...
        :return:
        """

        # Determine the path
        path = fs.join(self.truncation_path, "disk.fits")

        # Open the frame and return it
        return Frame.from_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_frame(self):

        """
        This function ...
        :return:
        """

        # Determine the path
        path = fs.join(self.truncation_path, "bulge.fits")

        # Open the frame and return it
        return Frame.from_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_ellipse(self):

        """
        This function ...
        :return:
        """

        # Load the ellipse
        path = fs.join(self.truncation_path, "ellipse.reg")
        region = SkyRegion.from_file(path)
        ellipse = region[0]

        # Return the (sky) ellipse
        return ellipse

    # -----------------------------------------------------------------

    def truncation_mask(self, frame):

        """
        This function ...
        :param frame: frame, image, datacube ...
        :return:
        """

        # Convert sky ellipse to pixel ellipse and then to mask
        return self.truncation_ellipse.to_pixel(frame.wcs).to_mask(frame.xsize, frame.ysize)

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_parameters(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the parameters file
        path = fs.join(self.components_path, "parameters.dat")

        # Load the galaxy parameters
        parameters = load_galaxy_parameters(path)

        # Return the parameter map
        return parameters

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stars_map(self):

        """
        This function ...
        :return:
        """

        # Open and return the frame of the old stellar distribution
        return Frame.from_file(self.old_stellar_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stars_map(self):

        """
        This function ...
        :return:
        """

        # Open and return the frame of the young stellar distribution
        return Frame.from_file(self.young_stellar_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_stars_map(self):

        """
        This function ...
        :return:
        """

        # Open and return the frame of the ionizing stellar distribution
        return Frame.from_file(self.ionizing_stellar_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map(self):

        """
        This function ...
        :return:
        """

        # Open and return the frame of the dust distribution
        return Frame.from_file(self.dust_map_path)

# -----------------------------------------------------------------

def get_free_parameter_labels(table_path):

    """
    This function ...
    :param table_path:
    :return:
    """

    return dict(np.genfromtxt(table_path, delimiter=" | ", dtype=str)) if fs.is_file(table_path) else None

# -----------------------------------------------------------------

def load_galaxy_parameters(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Create parameters structure
    parameters = Map()
    parameters.bulge = Map()
    parameters.disk = Map()

    ra = None
    dec = None

    # Read the parameter file
    with open(path, 'r') as parameter_file:

        # Loop over all lines in the file
        for line in parameter_file:

            splitted = line.split(": ")

            # Bulge parameters
            if splitted[0] == "Bulge":

                if splitted[1] == "Relative contribution": parameters.bulge.f = float(splitted[2])
                elif splitted[1] == "IRAC 3.6um flux density": parameters.bulge.fluxdensity = parsing.quantity(splitted[2])
                elif splitted[1] == "Axial ratio": parameters.bulge.q = float(splitted[2])
                elif splitted[1] == "Position angle": parameters.bulge.PA = parsing.angle(splitted[2])
                elif splitted[1] == "Effective radius": parameters.bulge.Re = parsing.quantity(splitted[2])
                elif splitted[1] == "Sersic index": parameters.bulge.n = float(splitted[2])

            # Disk parameters
            elif splitted[0] == "Disk":

                if splitted[1] == "Relative contribution": parameters.disk.f = float(splitted[2])
                elif splitted[1] == "IRAC 3.6um flux density": parameters.disk.fluxdensity = parsing.quantity(splitted[2])
                elif splitted[1] == "Axial ratio": parameters.disk.q = float(splitted[2])
                elif splitted[1] == "Position angle": parameters.disk.PA = parsing.angle(splitted[2])
                elif splitted[1] == "Central surface brightness": parameters.disk.mu0 = parsing.quantity(splitted[2])
                elif splitted[1] == "Exponential scale length": parameters.disk.hr = parsing.quantity(splitted[2])

            # Other parameters
            elif len(splitted) == 2:

                if splitted[0] == "Name": parameters.galaxy_name = splitted[1]
                elif splitted[0] == "Center RA": ra = parsing.quantity(splitted[1])
                elif splitted[0] == "Center DEC": dec = parsing.quantity(splitted[1])
                elif splitted[0] == "Major axis length": parameters.major = parsing.quantity(splitted[1])
                elif splitted[0] == "Ellipticity": parameters.ellipticity = float(splitted[1])
                elif splitted[0] == "Position angle": parameters.position_angle = parsing.angle(splitted[1])
                elif splitted[0] == "Distance": parameters.distance = parsing.quantity(splitted[1])
                elif splitted[0] == "Distance error": parameters.distance_error = parsing.quantity(splitted[1])
                elif splitted[0] == "Inclination": parameters.inclination = parsing.angle(splitted[1])
                elif splitted[0] == "IRAC 3.6um flux density": parameters.i1_fluxdensity = parsing.quantity(splitted[1])
                elif splitted[0] == "IRAC 3.6um flux density error": parameters.i1_error = parsing.quantity(splitted[1])
                elif splitted[0] == "IRAC 4.5um flux density": parameters.i2_fluxdensity = parsing.quantity(splitted[1])
                elif splitted[0] == "IRAC 4.5um flux density error": parameters.i2_error = parsing.quantity(splitted[1])

    # Add the center coordinate
    parameters.center = SkyCoordinate(ra=ra, dec=dec)

    # Return the parameters
    return parameters

# -----------------------------------------------------------------
