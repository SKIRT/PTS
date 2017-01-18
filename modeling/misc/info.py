#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.misc.info Contains the InfoShower class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from ..component.component import ModelingComponent
from ...core.tools import formatting as fmt
from ...core.tools import filesystem as fs
from ...magic.core.frame import Frame
from ..component.galaxy import load_preparation_statistics

# -----------------------------------------------------------------

class InfoShower(ModelingComponent):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(InfoShower, self).__init__(config)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ..
        :param kwargs:
        :return:
        """

        # Call the setup function
        self.setup(**kwargs)

        # Show general info
        self.show_general()

        # Show depending on which type of object
        if self.modeling_type == "galaxy": self.show_galaxy()
        elif self.modeling_type == "other": self.show_other()

        # Show fitting info
        if self.fitting_configuration is not None: self.show_fitting()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(InfoShower, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def show_general(self):

        """
        This function ...
        :return:
        """

        print("")
        print("Object name: " + self.modeling_configuration.name)
        print("Modeling type: " + self.modeling_configuration.modeling_type)

    # -----------------------------------------------------------------

    def show_galaxy(self):

        """
        This function ...
        :return:
        """

        # Get preparation statistics
        statistics = load_preparation_statistics(self.config.path)

        # Maps paths
        maps_path = fs.join(self.config.path, "maps")
        old_stars_path = fs.join(maps_path, "old_stars.fits")
        young_stars_path = fs.join(maps_path, "young_stars.fits")
        ionizing_stars_path = fs.join(maps_path, "ionizing_stars.fits")
        dust_path = fs.join(maps_path, "dust.fits")

        # Open the old stars map
        old_stars = Frame.from_file(old_stars_path)

        # Get convolution and rebinning filter
        convolution_filter = statistics.convolution_filter
        rebinning_filter = statistics.rebinning_filter

        # Print info
        print("Galaxy NGC name: " + self.modeling_configuration.ngc_name)
        print("Modeling method: " + self.modeling_configuration.modeling_method)
        print("Model pixelscale: " + str(old_stars.average_pixelscale) + " (" + str(rebinning_filter) + ")")
        print("Model resolution (FWHM): " + str(old_stars.fwhm) + " (" + str(convolution_filter) + ")")

    # -----------------------------------------------------------------

    def show_other(self):

        """
        This function ...
        :return:
        """

        print("")

    # -----------------------------------------------------------------

    def show_fitting(self):

        """
        This function ...
        :return:
        """

        print("Reference fluxes (for the SED fitting): ")
        print("")
        for filter_name in self.fitting_configuration.filters: print(" - " + filter_name)
        print("")

# -----------------------------------------------------------------

