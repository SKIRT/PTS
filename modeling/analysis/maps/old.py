#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.old Contains the OldMapsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import MapsAnalysisComponent
from ....core.basics.log import log
from ....magic.maps.oldstars.disk import DiskOldStellarMapMaker
from ....magic.core.list import FrameList
from ....core.tools import filesystem as fs
from ....magic.core.frame import Frame

# -----------------------------------------------------------------

class OldMapsAnalyser(MapsAnalysisComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(OldMapsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # 2. Make the maps
        self.make_maps()

        # 3. Analyse
        self.analyse_maps()

        # 4. Write
        self.write()

        # 5. Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(OldMapsAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making old stellar maps ...")

        # Set the method name
        method_name = "disk"

        # Get current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Create the maker
        maker = DiskOldStellarMapMaker()

        # Get the I1 frame
        i1 = self.get_frame_for_filter(self.i1_filter)
        frames = FrameList(i1)

        # Get the bulge frame
        bulge = self.bulge_frame
        bulges = FrameList(i1=bulge)

        # Run
        maker.run(frames=frames, bulges=bulges, method_name=method_name, maps=current)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def analyse_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the maps ...")

        # Analyse the residuals
        self.analyse_residuals()

    # -----------------------------------------------------------------

    def analyse_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the map residuals ...")

        # Get the method and name of the old stellar map used for the model, based on observation
        method_name, map_name = self.analysis_run.old_map_method_and_name

        # Create directory for the analysis of this map
        path = self.get_path_for_map(map_name, method=method_name, add_extension=False)
        if not fs.is_directory(path): fs.create_directory(path)

        # Determine the path to the residual map
        residuals_path = fs.join(path, "residuals.fits")

        # Determine the path to the distribution
        distribution_path = fs.join(path, "distribution.dat")

        # Not yet created
        if not fs.is_file(residuals_path):

            # Create the residuals frame
            residuals = self.create_residuals_for_map(self.analysis_run.model_old_map, map_name, method_name)

            # Write the residuals map
            residuals.saveto(residuals_path)

        # Load from file
        else: residuals = Frame.from_file(residuals_path)

        # Not yet created
        if not fs.is_file(distribution_path):

            distribution = self.create_residuals_distribution(residuals, nbins=20)

            # Save the distribution
            distribution.saveto(distribution_path)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return:
        """

        return self.old_path

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the colour maps
        self.write_maps()

        # Write origins
        self.write_origins()

        # Write the methods
        self.write_methods()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the maps
        self.plot_maps()

# -----------------------------------------------------------------
