#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.young Contains the YoungMapsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsAnalysisComponent
from ....core.basics.log import log
from ....magic.maps.youngstars.young import YoungStellarMapsMaker
from ....core.tools import filesystem as fs
from ....magic.core.frame import Frame

# -----------------------------------------------------------------

class YoungMapsAnalyser(MapsAnalysisComponent):

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
        super(YoungMapsAnalyser, self).__init__(*args, **kwargs)

        # The input FUV and FUV error maps
        self.fuv = None
        self.fuv_errors = None

        # The map of the old stellar disk
        self.old = None

        # The maps of FUV attenuation
        self.fuv_attenuations = None

        # The origins
        self.old_origin = None
        self.fuv_attenuations_origins = None

        # Methods
        self.old_method = None
        self.fuv_attenuations_methods = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # 2. Load the necessary input maps
        self.load_input()

        # 3. Make the maps
        self.make_maps()

        # 4. Analyse the maps
        self.analyse_maps()

        # 5. Write
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
        super(YoungMapsAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary input ...")

        # Load the GALEX FUV image and error map
        self.load_fuv()

        # Load FUV attenuation map
        self.load_fuv_attenuation_maps()

        # Load old stellar map
        self.load_old_stellar_map()

    # -----------------------------------------------------------------

    def load_fuv(self):

        """
        This function ...
        :return:
        """

        # Get FUV frame and error map
        self.fuv = self.dataset.get_frame("GALEX FUV") # in original MJy/sr units
        self.fuv_errors = self.dataset.get_errormap("GALEX FUV") # in original MJy/sr units

    # -----------------------------------------------------------------

    def load_fuv_attenuation_maps(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Loading the maps of the FUV attenuation ...")

        # Get the FUV attenuation maps
        #self.fuv_attenuations, self.fuv_attenuations_origins = self.get_fuv_attenuation_maps_and_origins(flatten=True)
        self.fuv_attenuations, self.fuv_attenuations_origins, self.fuv_attenuations_methods = self.get_fuv_attenuation_maps_origins_and_methods(flatten=True)

    # -----------------------------------------------------------------

    def load_old_stellar_map(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Loading the map of old stars ...")

        # Get the map
        self.old = self.get_old_stellar_disk_map(self.i1_filter)

        # Set the old origin
        self.old_origin = self.i1_filter

        # Set the old method
        self.old_method = "disk" #self.get_old_stellar_disk_methods()

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making young stellar maps ...")

        # Get the current maps
        if self.config.remake: current = dict()
        else: current = self.current_maps

        # Create the map maker
        maker = YoungStellarMapsMaker()

        # Set the factors
        factors = self.config.factor_range.linear(self.config.factor_nvalues, as_list=True)

        # Run the map maker
        maker.run(fuv=self.fuv, fuv_errors=self.fuv_errors, old=self.old, fuv_attenuations=self.fuv_attenuations,
                  factors=factors, old_origin=self.old_origin, fuv_attenuations_origins=self.fuv_attenuations_origins,
                  old_method=self.old_method, fuv_attenuations_methods=self.fuv_attenuations_methods, maps=current)

        # Set the maps
        self.maps = maker.maps

        # Set the origins
        self.origins = maker.origins

        # Set the methods
        self.methods = maker.methods

    # -----------------------------------------------------------------

    def analyse_maps(self):

        """
        Thsij function ...
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
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the map residuals ...")

        # Get the method and name of the young stellar map used for the model, based on observation
        method_name, map_name = self.analysis_run.young_map_method_and_name

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
            residuals = self.create_residuals_for_map(self.analysis_run.model_young_map, map_name, method_name)

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

        return self.young_path

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
