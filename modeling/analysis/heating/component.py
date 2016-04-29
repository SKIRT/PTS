#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating.component Contains the DustHeatingAnalysisComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component import AnalysisComponent
from ....core.tools import filesystem as fs

# -----------------------------------------------------------------

class DustHeatingAnalysisComponent(AnalysisComponent):
    
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
        super(DustHeatingAnalysisComponent, self).__init__(config)

        # -- Attributes --

        # The different contributing components
        self.contributions = ["total", "old", "young", "ionizing"]
        self.component_names = {"old": ["Evolved stellar bulge", "Evolved stellar disk"],
                                "young": "Young stars",
                                "ionizing": "Ionizing stars"}

        # The paths to the analysis/heating/total, analysis/heating/old, analysis/heating/young and analysis/heating/ionizing directory
        self.simulation_paths = dict()
        self.output_paths = dict()
        self.ski_paths = dict()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DustHeatingAnalysisComponent, self).setup()

        # Set the paths to the different simulation directories and corresponding output directories
        for contribution in self.contributions:

            # Set the simulation path
            simulation_path = fs.join(self.analysis_heating_path, contribution)

            # Create the simulation directory if it is not present
            if not fs.is_directory(simulation_path): fs.create_directory(simulation_path)

            # Set the path to the output directory
            output_path = fs.join(simulation_path, "out")

            # Create the output directory if it is not present
            if not fs.is_directory(simulation_path): fs.create_directory(output_path)

            # Add the paths to the appropriate dictionaries
            self.simulation_paths[contribution] = simulation_path
            self.output_paths[contribution] = output_path

            # Determine the path to the ski file for this contribution
            ski_path = fs.join(simulation_path, self.galaxy_name + ".ski")

            # Set the ski file path
            self.ski_paths[contribution] = ski_path

# -----------------------------------------------------------------
