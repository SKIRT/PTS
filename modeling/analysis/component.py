#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.component Contains the AnalysisComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem

# -----------------------------------------------------------------

class AnalysisComponent(ModelingComponent):
    
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
        super(AnalysisComponent, self).__init__(config)

        # -- Attributes --

        # The path to the analysis/out directory
        self.analysis_out_path = None

        # The path to the analysis/attenuation directory
        self.analysis_attenuation_path = None

        # The path to the analysis/colours directory
        self.analysis_colours_path = None

        # The path to the analysis/residuals directory
        self.analysis_residuals_path = None

        # The path to the analysis/heating directory
        self.analysis_heating_path = None

        # The path to the ski file and the wavelength grid file
        self.analysis_ski_path = None
        self.analysis_wavelengths_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisComponent, self).setup()

        # Set the path to the analysis/out path
        self.analysis_out_path = filesystem.join(self.analysis_path, "out")

        # Set the path to the analysis/attenuation path
        self.analysis_attenuation_path = filesystem.join(self.analysis_path, "attenuation")

        # Set the path to the analysis/colours path
        self.analysis_colours_path = filesystem.join(self.analysis_path, "colours")

        # Set the path to the analysis/residuals path
        self.analysis_residuals_path = filesystem.join(self.analysis_path, "residuals")

        # Set the path to the analysis/heating path
        self.analysis_heating_path = filesystem.join(self.analysis_path, "heating")

        # Create the analysis/out and fit/out directories
        filesystem.create_directories([self.analysis_out_path, self.analysis_attenuation_path,
                                       self.analysis_colours_path, self.analysis_residuals_path,
                                       self.analysis_heating_path])

        # Set the path to the ski file and wavelength grid file
        self.analysis_ski_path = filesystem.join(self.analysis_path, self.galaxy_name + ".ski")
        self.analysis_wavelengths_path = filesystem.join(self.analysis_path, "wavelengths.txt")

# -----------------------------------------------------------------
