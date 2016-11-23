#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeler Contains the GalaxyModeler class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.basics.configurable import Configurable
from ..core.tools.logging import log
from ..core.tools import filesystem as fs
from .setup import ModelingSetupTool
from .data.properties import PropertyFetcher
from .data.images import ImageFetcher
from .data.seds import SEDFetcher
from .preparation.initialization import PreparationInitializer
from .preparation.preparer import DataPreparer
from .decomposition.decomposition import GalaxyDecomposer
from .truncation.truncation import Truncator
from .photometry.photometry import PhotoMeter
from .maps.stars.old import OldStellarMapMaker
from .maps.stars.young import YoungStellarMapMaker
from .maps.stars.ionizing import IonizingStellarMapMaker
from .maps.dust.dust import DustMapMaker
from .fitting.configuration import FittingConfigurer
from .fitting.initialization import FittingInitializer
from .fitting.explorer import ParameterExplorer

# -----------------------------------------------------------------

class GalaxyModeler(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(GalaxyModeler, self).__init__(config)

        # The path to the modeling directory
        self.modeling_path = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Get the data
        self.get_data()

        # Data preparation
        self.prepare_data()

        # Decomposition
        self.decompose()

        # Truncation
        self.truncate()

        # Do the photometry
        self.photometry()

        # Make the maps
        self.make_maps()

        # Do the fitting
        self.fit()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GalaxyModeler, self).setup()

        # Set the path to the modeling directory
        self.modeling_path = fs.join(self.config.path, self.config.galaxy_name)

        # Create the modeling setup tool
        tool = ModelingSetupTool()

        # Run the setup tool
        tool.run()

    # -----------------------------------------------------------------

    def get_data(self):

        """
        This function ...
        :return:
        """

        self.get_properties()

        self.get_images()

        self.get_seds()

    # -----------------------------------------------------------------

    def get_properties(self):

        """
        This function ...
        :return:
        """

        fetcher = PropertyFetcher()

        fetcher.run()

    # -----------------------------------------------------------------

    def get_images(self):

        """
        This function ...
        :return:
        """

        fetcher = ImageFetcher()

        fetcher.run()

    # -----------------------------------------------------------------

    def get_seds(self):

        """
        This function ...
        :return:
        """

        fetcher = SEDFetcher()

        fetcher.run()

    # -----------------------------------------------------------------

    def prepare_data(self):

        """
        This function ...
        :return:
        """

        self.initialize_preparation()

        self.prepare()

    # -----------------------------------------------------------------

    def initialize_preparation(self):

        """
        This function ...
        :return:
        """

        initializer = PreparationInitializer()

        initializer.run()

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        preparer = DataPreparer()

        preparer.run()

    # -----------------------------------------------------------------

    def decompose(self):

        """
        This function ...
        :return:
        """

        decomposer = GalaxyDecomposer()

        decomposer.run()

    # -----------------------------------------------------------------

    def truncate(self):

        """
        This function ...
        :return:
        """

        truncator = Truncator()

        truncator.run()

    # -----------------------------------------------------------------

    def photometry(self):

        """
        This function ...
        :return:
        """

        photometer = PhotoMeter()

        photometer.run()

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        self.make_old_stellar_map()

        self.make_young_stellar_map()

        self.make_ionizing_stellar_map()

        self.make_dust_map()

    # -----------------------------------------------------------------

    def make_old_stellar_map(self):

        """
        This function ...
        :return:
        """

        maker = OldStellarMapMaker()

        maker.run()

    # -----------------------------------------------------------------

    def make_young_stellar_map(self):

        """
        This function ...
        :return:
        """

        maker = YoungStellarMapMaker()

        maker.run()

    # -----------------------------------------------------------------

    def make_ionizing_stellar_map(self):

        """
        This function ...
        :return:
        """

        maker = IonizingStellarMapMaker()

        maker.run()

    # -----------------------------------------------------------------

    def make_dust_map(self):

        """
        This function ...
        :return:
        """

        maker = DustMapMaker()

        maker.run()

    # -----------------------------------------------------------------

    def fit(self):

        """
        This function ...
        :return:
        """

        self.configure_fit()

        self.initialize_fit()

        self.explore()

    # -----------------------------------------------------------------

    def configure_fit(self):

        """
        This function ...
        :return:
        """

        configurer = FittingConfigurer()

        configurer.run()

    # -----------------------------------------------------------------

    def initialize_fit(self):

        """
        This function ...
        :return:
        """

        initializer = FittingInitializer()

        initializer.run()

    # -----------------------------------------------------------------

    def explore(self):

        """
        This function ...
        :return:
        """

        explorer = ParameterExplorer()

        explorer.run()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
