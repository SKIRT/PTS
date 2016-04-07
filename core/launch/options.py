#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.options Contains the AnalysisOptions class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class AnalysisOptions(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # From RemoteSimulation class:

        # Options for analysis of the simulation output
        self.extract_progress = False
        self.extract_timeline = False
        self.extract_memory = False
        self.plot_progress = False
        self.plot_timeline = False
        self.plot_memory = False
        self.plot_seds = False
        self.plot_grids = False
        self.make_rgb = False
        self.make_wave = False
        self.calculate_observed_fluxes = False
        self.make_observed_images = False

        # Extraction, plotting and 'misc' path
        self.extraction_path = None
        self.plotting_path = None
        self.misc_path = None

        # The path to a file containing an SED for which the points have to be plotted against the simulated curve (when plot_seds is enabled)
        self.reference_sed = None

        # The filters for which to recreate the observations
        self.observation_filters = None

        # Removal options
        self.remove_remote_input = True
        self.remove_remote_output = True
        self.remove_remote_simulation_directory = True

        # Local removal (after analysis)
        self.remove_local_output = False

        # Properties relevant for simulations part of a scaling test
        self.scaling_run_name = None
        self.scaling_data_file = None
        self.scaling_plot_path = None

        # Properties relevant for simulations part of radiative transfer modeling
        self.modeling_path = None

# -----------------------------------------------------------------
