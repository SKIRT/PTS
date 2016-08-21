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
import copy
import warnings

# Import the relevant PTS classes and modules
from ..basics.map import Map
from ..tools.logging import log

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
                if isinstance(getattr(self, option), Map): getattr(self, option).set_items(options[option])

                # If it is a simple variable, just use setattr to set the attribute of this class
                else: setattr(self, option, options[option])

            # If the option does not exist, ignore it but give a warning
            else: warnings.warn("The option " + option + " does not exist")

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

# -----------------------------------------------------------------

class LoggingOptions(Options):

    """
    This function ...
    """

    def __init__(self, brief=False, verbose=False, memory=False, allocation=False, allocation_limit=1e-5):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(LoggingOptions, self).__init__()

        # Set options
        self.brief = brief     # Brief console logging
        self.verbose = verbose # Verbose logging
        self.memory = memory   # State the amount of used memory with each log message
        self.allocation = allocation    # Write log messages with the amount of (de)allocated memory
        self.allocation_limit = allocation_limit  # The lower limit for the amount of (de)allocated memory to be logged

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
        self.misc.observation_instruments = None # The instrument for which to recreate the observations
        self.misc.make_images_remote = None  # Perform the calculation of the observed images on a remote machine (this is a memory and CPU intensive step)
        self.misc.images_wcs = None  # the path to the FITS file from which the WCS should be set as the WCS of the simulated images
        self.misc.images_unit = None # the unit to which the simulated images should be converted (if None, the original unit is kept)
        self.misc.images_kernels = None # the paths to the FITS file of convolution kernel used for convolving the observed images (a dictionary where the keys are the filter names)

        # Properties that are relevant for simulations launched as part of a batch (e.g. from an automatic launching procedure)
        self.timing_table_path = None
        self.memory_table_path = None

        # Properties relevant for simulations part of a scaling test
        self.scaling_path = None
        self.scaling_run_name = None

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

    def check(self, logging_options):

        """
        This function ...
        :param logging_options:
        :return:
        """

        # Inform the user
        log.info("Checking the analysis options ...")

        # PLOTTING

        # If any plotting setting has been enabled, check whether the plotting path has been set
        if self.plotting.progress or self.plotting.memory or self.plotting.timeline or self.plotting.seds or self.plotting.grids:
            if self.plotting.path is None: raise ValueError("The plotting path has not been set")

        # If progress plotting has been enabled, enabled progress extraction
        if self.plotting.progress and not self.extraction.progress:
            log.warning("Progress plotting is enabled so progress extraction will also be enabled")
            self.extraction.progress = True

        # If memory plotting has been enabled, enable memory extraction
        if self.plotting.memory and not self.extraction.memory:
            log.warning("Memory plotting is enabled so memory extraction will also be enabled")
            self.extraction.memory = True

        # If timeline plotting has been enabled, enable timeline extraction
        if self.plotting.timeline and not self.extraction.timeline:
            log.warning("Timeline plotting is enabled so timeline extraction will also be enabled")
            self.extraction.timeline = True

        # EXTRACTION

        # If any extraction setting has been enabled, check whether the extraction path has been set
        if self.extraction.progress or self.extraction.memory or self.extraction.timeline:
            if self.extraction.path is None: raise ValueError("The extraction path has not been set")

        # If memory extraction has been enabled, enable memory logging
        if self.extraction.memory and not logging_options.memory:
            log.warning("Memory extraction is enabled so memory logging will also be enabled")
            logging_options.memory = True

# -----------------------------------------------------------------
