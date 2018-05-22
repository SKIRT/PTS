#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsAnalysisComponent
from ....core.basics.log import log
from ....magic.maps.ionizingstars.ionizing import IonizingStellarMapsMaker
from ....core.tools import filesystem as fs
from ....magic.core.frame import Frame

# -----------------------------------------------------------------

class IonizingMapsAnalyser(MapsAnalysisComponent):

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
        super(IonizingMapsAnalyser, self).__init__(*args, **kwargs)

        # The Halpha map
        self.halpha = None

        # The maps of hot dust
        self.hots = None

        # Origins
        self.hots_origins = None

        # Methods
        self.hots_methods = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # 2. Load
        self.load_frames()

        # 3. Make the maps
        self.make_maps()

        # 4. Analyse the maps
        self.analyse_maps()

        # 5. Writing
        self.write()

        # 6. Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(IonizingMapsAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary data ...")

        # Load the MIPS 24 micron image (and convert to solar units == > NO?)
        self.load_hot()

        # Load the H alpha image (and convert to solar units == > NO?)
        self.load_halpha()

    # -----------------------------------------------------------------

    def load_hot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the maps of hot dust ...")

        # Get hot dust maps
        self.hots = self.get_hot_dust_maps()

        # Get hot dust origins
        self.hots_origins = self.get_hot_dust_origins()

        # Methods
        self.hots_methods = self.get_hot_dust_methods()

    # -----------------------------------------------------------------

    def load_halpha(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the H-alpha image and converting to solar units ...")

        # Get the H-alpha image (THE OBSERVED ONE)
        self.halpha = self.dataset.get_frame_for_filter("Halpha") #self.get_frame_for_filter(parse_filter("Halpha"))

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Making ionizing stellar maps ...")

        # Get the current maps
        if self.config.remake: current = dict()
        else: current = self.current_maps

        # Create
        maker = IonizingStellarMapsMaker()

        # Run
        maker.run(halpha=self.halpha, hots=self.hots, hots_origins=self.hots_origins, hots_methods=self.hots_methods, maps=current)

        # Set the maps
        self.maps = maker.maps

        # Set the origins
        self.origins = maker.origins

        # Set the methods
        self.methods = maker.methods

    # -----------------------------------------------------------------

    def analyse_maps(self):

        """
        Thisj function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the maps ...")

        # Analyse residuals
        self.analyse_residuals()

        # Analyse correlations
        self.analyse_correlations()

    # -----------------------------------------------------------------

    def analyse_residuals(self):

        """
        Thisj function ....
        :return:
        """

        # Inform the user
        log.info("Analysing the map residuals ...")

        # Get the method and name of the ionizing stellar map used for the model, based on observation
        method_name, map_name = self.analysis_run.ionizing_map_method_and_name

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
            residuals = self.create_residuals_for_map(self.analysis_run.model_ionizing_map, map_name, method_name)

            # Write the residuals map
            residuals.saveto(residuals_path)

        # Load from file
        else: residuals = Frame.from_file(residuals_path)

        # Not yet created
        if not fs.is_file(distribution_path):

            # Create distribution
            distribution = self.create_residuals_distribution(residuals, nbins=20)

            # Save the distribution
            distribution.saveto(distribution_path)

    # -----------------------------------------------------------------

    def analyse_correlations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the map correlations ...")

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_path

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
