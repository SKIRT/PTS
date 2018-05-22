#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.maps Contains the MapsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
#from .component import MapsAnalysisComponent
from ..component import AnalysisComponent
from ....core.basics.log import log
from .colours import ColourMapsAnalyser
from .ssfr import SSFRMapsAnalyser
from .tir import TIRMapsAnalyser
from .attenuation import AttenuationMapsAnalyser
from .old import OldMapsAnalyser
from .dust import DustMapsAnalyser
from .young import YoungMapsAnalyser
from .ionizing import IonizingMapsAnalyser

# -----------------------------------------------------------------

class MapsAnalyser(AnalysisComponent):

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
        super(MapsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def needs_colour_maps(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def needs_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def needs_tir_maps(self):

        """
        Thisf ucntion ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def needs_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def needs_old_stellar_maps(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def needs_dust_maps(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def needs_young_stellar_maps(self):

        """
        Thisn function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def needs_ionizing_stellar_maps(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # 2. Analyse colour maps
        if self.needs_colour_maps: self.analyse_colour_maps()

        # 3. Analyse sSFR maps
        if self.needs_ssfr_maps: self.analyse_ssfr_maps()

        # 4. Analyse TIR maps
        if self.needs_tir_maps: self.analyse_tir_maps()

        # 5. Analyse the attenuation map(s)
        if self.needs_attenuation_maps: self.analyse_attenuation_maps()

        # 6. Analyse the map of the old stellar disk
        if self.needs_old_stellar_maps: self.analyse_old_stellar_maps()

        # 7. Analyse the dust map
        if self.needs_dust_maps: self.analyse_dust_maps()

        # 8. Analyse the map of the young stellar population
        if self.needs_young_stellar_maps: self.analyse_young_stellar_maps()

        # 9. Analyse the map of the ionizing stellar population
        if self.needs_ionizing_stellar_maps: self.analyse_ionizing_stellar_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(MapsAnalyser, self).setup(**kwargs)

        # Load the analysis run
        #self.load_run()

    # -----------------------------------------------------------------

    def analyse_colour_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing colour maps ...")

        # Create the analyser
        analyser = ColourMapsAnalyser()

        # Set modeling path
        analyser.config.path = self.config.path

        # Set run name
        analyser.config.run = self.config.run

        # Run the analyser
        analyser.run()

    # -----------------------------------------------------------------

    def analyse_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing sSFR maps ...")

        # Create the analyser
        analyser = SSFRMapsAnalyser()

        # Set modeling path
        analyser.config.path = self.config.path

        # Set run name
        analyser.config.run = self.config.run

        # Run the analyser
        analyser.run()

    # -----------------------------------------------------------------

    def analyse_tir_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing TIR maps ...")

        # Create the analyser
        analyser = TIRMapsAnalyser()

        # Set modeling path
        analyser.config.path = self.config.path

        # Set run name
        analyser.config.run = self.config.run

        # Run the analyser
        analyser.run()

    # -----------------------------------------------------------------

    def analyse_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing attenuation maps ...")

        # Create the analyser
        analyser = AttenuationMapsAnalyser()

        # Set modeling path
        analyser.config.path = self.config.path

        # Set run name
        analyser.config.run = self.config.run

        # Run the analyser
        analyser.run()

    # -----------------------------------------------------------------

    def analyse_old_stellar_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing old stellar maps ...")

        # Create the analyser
        analyser = OldMapsAnalyser()

        # Set modeling path
        analyser.config.path = self.config.path

        # Set run name
        analyser.config.run = self.config.run

        # Run the analyser
        analyser.run()

    # -----------------------------------------------------------------

    def analyse_dust_maps(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Analysing dust maps ...")

        # Create the analyser
        analyser = DustMapsAnalyser()

        # Set modeling path
        analyser.config.path = self.config.path

        # Set run name
        analyser.config.run = self.config.run

        # Run the analyser
        analyser.run()

    # -----------------------------------------------------------------

    def analyse_young_stellar_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing young stellar maps ...")

        # Create the analyser
        analyser = YoungMapsAnalyser()

        # Set modeling path
        analyser.config.path = self.config.path

        # Set run name
        analyser.config.run = self.config.run

        # Run the analyser
        analyser.run()

    # -----------------------------------------------------------------

    def analyse_ionizing_stellar_maps(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Analysing ionizing stellar maps ...")

        # Create the analyser
        analyser = IonizingMapsAnalyser()

        # Set modeling path
        analyser.config.path = self.config.path

        # Set run name
        analyser.config.run = self.config.run

        # Run the analyser
        analyser.run()

# -----------------------------------------------------------------
