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

# Import standard modules
import warnings

# Import the relevant PTS classes and modules
from ..basics.map import Map

# -----------------------------------------------------------------

class Options(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        pass

    # -----------------------------------------------------------------

    def set_options(self, options):

        """
        This function allows setting multiple options at once from a dictionary
        :param options:
        :return:
        """

        # Loop over all the options defined in the 'options' dictionary
        for option in options:

            # Check whether an option with this name exists in this class
            if hasattr(self, option):

                # Check if the option is composed of other options (a Map), or if it is just a simple variable
                if isinstance(getattr(self, option), Map):
                    getattr(self, option).set_items(options[option])

                # If it is a simple variable, just use setattr to set the attribute of this class
                else:
                    setattr(self, option, options[option])

            # If the option does not exist, ignore it but give a warning
            else: warnings.warn("The option " + option + " does not exist")

# -----------------------------------------------------------------

class SchedulingOptions(object):

    """
    This function ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(SchedulingOptions, self).__init__()

        # Scheduling options
        self.nodes = None
        self.ppn = None
        self.mail = None
        self.full_node = None
        self.walltime = None
        self.local_jobscript_path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_dict(cls, dictionary):

        """
        This function ...
        :param dictionary:
        :return:
        """

        # Create a new SchedulingOptions instance
        options = cls()

        if "nodes" in dictionary: options.nodes = dictionary["nodes"]
        if "ppn" in dictionary: options.ppn = dictionary["ppn"]
        if "mail" in dictionary: options.mail = dictionary["mail"]
        if "full_node" in dictionary: options.full_node = dictionary["full_node"]
        if "walltime" in dictionary: options.walltime = dictionary["walltime"]
        if "local_jobscript_path" in dictionary: options.local_jobscript_path = dictionary["local jobscript path"]

        # Return the scheduling options object
        return options

# -----------------------------------------------------------------

class AnalysisOptions(Options):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(AnalysisOptions, self).__init__()

        # Options for extracting data from the simulation's log files
        self.extraction = Map()
        self.extraction.path = None
        self.extraction.progress = False
        self.extraction.timeline = False
        self.extraction.memory = False

        # Options for plotting simulation output
        self.plotting = Map()
        self.plotting.path = None
        self.plotting.format = "pdf" # can be 'pdf', 'png' or other formats supported by MatplotLib
        self.plotting.progress = False
        self.plotting.timeline = False
        self.plotting.memory = False
        self.plotting.seds = False
        self.plotting.grids = False
        self.plotting.reference_sed = None # The path to a file containing an SED for which the points have to be
                                           # plotted against the simulated curve (when plotting seds is enabled)

        # Options for creating data of various formats
        self.misc = Map()
        self.misc.path = None
        self.misc.rgb = False
        self.misc.wave = False
        self.misc.fluxes = False
        self.misc.images = False
        self.misc.observation_filters = None # The filters for which to recreate the observations
        self.misc.make_images_remote = None  # Perform the calculation of the observed images on a remote machine (this is a memory and CPU intensive step)

        # Properties that are relevant for simulations launched as part of a batch (e.g. from an automatic launching procedure)
        self.timing_table_path = None
        self.memory_table_path = None

        # Properties relevant for simulations part of a scaling test
        self.scaling_run_name = None
        self.scaling_data_file = None
        self.scaling_plot_path = None

        # Properties relevant for simulations part of radiative transfer modeling
        self.modeling_path = None

    # -----------------------------------------------------------------

    @property
    def any_extraction(self):

        """
        This function ...
        :return:
        """

        return self.extraction.progress or self.extraction.timeline or self.extraction.memory

    # -----------------------------------------------------------------

    @property
    def any_plotting(self):

        """
        This function ...
        :return:
        """

        return self.plotting.seds or self.plotting.grids or self.plotting.progress or self.plotting.timeline or self.plotting.memory

    # -----------------------------------------------------------------

    @property
    def any_misc(self):

        """
        This function ...
        :return:
        """

        return self.misc.rgb or self.misc.wave or self.misc.fluxes or self.misc.images

# -----------------------------------------------------------------
