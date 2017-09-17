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

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.composite import SimplePropertyComposite

# -----------------------------------------------------------------

class Options(SimplePropertyComposite):

    """
    This class ...
    """

    def set_options(self, options):

        """
        This function allows setting multiple options at once from a dictionary
        :param options:
        :return:
        """

        # Call the set_properties function
        self.set_properties(options)

# -----------------------------------------------------------------

class LoggingOptions(Options):

    """
    This function ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(LoggingOptions, self).__init__()

        # Add properties
        self.add_property("brief", "boolean", "brief console logging", False)
        self.add_property("verbose", "boolean", "verbose logging", False)
        self.add_property("memory", "boolean", " state the amount of used memory with each log message", False)
        self.add_property("allocation", "boolean", "write log messages with the amount of (de)allocated memory", False)
        self.add_property("allocation_limit", "real", "lower limit for the amount of (de)allocated memory to be logged", 1e-5)

        # Set values
        self.set_properties(kwargs)

# -----------------------------------------------------------------

class SchedulingOptions(Options):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(SchedulingOptions, self).__init__()

        # Scheduling options
        self.add_property("nodes", "positive_integer", "number of nodes", None)
        self.add_property("ppn", "positive_integer", "number of processors per node", None)
        self.add_property("mail", "boolean", "send mails", None)
        self.add_property("full_node", "boolean", "use full nodes", None)
        self.add_property("walltime", "real", "expected walltime", None)
        self.add_property("local_jobscript_path", "string", None)

        # Set values
        self.set_properties(kwargs)

# -----------------------------------------------------------------

progress_name = "progress"
timeline_name = "timeline"
memory_name = "memory"
seds_name = "seds"
grids_name = "grids"
rgb_name = "rgb"
wave_name = "wave"
fluxes_name = "fluxes"
images_name = "images"

# -----------------------------------------------------------------

class AnalysisOptions(Options):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(AnalysisOptions, self).__init__()

        # Extraction
        self.add_section("extraction", "options for extractin data from the simulation's log files")
        self.extraction.add_property("path", "string", "extraction directory", None)
        self.extraction.add_property("progress", "boolean", "extract information about the progress in the different simulation phases", False)
        self.extraction.add_property("timeline", "boolean", "extract timeline information for the different simulation phases on the different processes", False)
        self.extraction.add_property("memory", "boolean", "extract information about the memory usage during the simulation", False)

        # Plotting
        self.add_section("plotting", "options for plotting simulation output")
        self.plotting.add_property("path", "string", "plotting directory", None)
        self.plotting.add_property("format", "string", "image format for the plots", "pdf", choices=["pdf", "png"])
        self.plotting.add_property("progress", "boolean", "make plots of the progress of the simulation phases as a function of time", False)
        self.plotting.add_property("timeline", "boolean", "plot the timeline for the different processes", False)
        self.plotting.add_property("memory", "boolean", "plot the memory consumption as a function of time", False)
        self.plotting.add_property("seds", "boolean", "make plots of the simulated SEDs", False)
        self.plotting.add_property("grids", "boolean", "make plots of the dust grid", False)
        self.plotting.add_property("reference_seds", "filepath_list", "path to a reference SED file against which the simulated SKIRT SEDs should be plotted", None)

        # Misc
        ## General
        self.add_section("misc", "settings for creating data of various types from the simulation output")
        self.misc.add_property("path", "string", "misc output directory", None)

        ## RGB images
        self.misc.add_property("rgb", "boolean", "make RGB images from the simulated datacube(s)", False)

        ## Wave movies
        self.misc.add_property("wave", "boolean", "make a wavelength movie through the simulated datacube(s)", False)

        ## Observed fluxes and images
        self.misc.add_property("observation_filters", "string_list", "the names of the filters for which to recreate the observations", None)
        self.misc.add_property("observation_instruments", "string_list", "the names of the instruments for which to recreate the observations", None)

        ## Observed fluxes
        self.misc.add_property("fluxes", "boolean", "calculate observed fluxes from the SKIRT output SEDs", False)
        self.misc.add_property("flux_errors", "string_string_dictionary", "errorbars for the different flux points of the mock observed SED")

        ## Observed images
        self.misc.add_property("images", "boolean", "make observed images form the simulated datacube(s)", False)
        self.misc.add_property("wcs_instrument", "string", "instrument for which to take the images_wcs as the WCS (if not a dictionary)")
        self.misc.add_property("images_wcs", "filepath_or_string_filepath_dictionary", "path to the FITS/txt file for which the WCS should be set as the WCS of the recreated observed images (single path or path for each datacube as a dictionary)", None)
        self.misc.add_property("images_unit", "string", "the unit to which the recreated observed images should be converted", None)

        ## CONVOLUTION
        self.misc.add_property("images_kernels", "string_string_dictionary", "paths to the FITS file of convolution kernel used for convolving the observed images (a dictionary where the keys are the filter names)", None)
        self.misc.add_property("images_psfs_auto", "boolean", "automatically determine the appropriate PSF kernel for each image", False)

        # Remote thresholds
        self.misc.add_property("make_images_remote", "string", "Perform the calculation of the observed images on a remote machine (this is a memory and CPU intensive step)", None)
        self.misc.add_property("images_remote_threshold", "data_quantity", "file size threshold for working with remote datacubes", "5 GB", convert_default=True) # 0.5 GB is a good size
        self.misc.add_property("rebin_remote_threshold", "data_quantity", "data size threshold for remote rebinning", "0.5 GB", convert_default=True)
        self.misc.add_property("convolve_remote_threshold", "data_quantity", "data size threshold for remote convolution", "1. GB", convert_default=True)

        ## REBINNING
        self.misc.add_property("rebin_instrument", "string", "instrument for which the rebin_wcs or rebin_dataset should be used")
        self.misc.add_property("rebin_wcs", "filepath_or_string_filepath_dictionary", "paths to the FITS/txt files of which the WCS should be used as the target for rebinning")
        self.misc.add_property("rebin_dataset", "file_path", "path to the dataset of which the coordinate systems are to be used as rebinning references")

        # OTHER FLAGS
        self.misc.add_property("spectral_convolution", "boolean", "use spectral convolution to calculate observed fluxes and create observed images", True)
        self.misc.add_property("group_images", "boolean", "group the images per instrument", False)


        # Properties that are relevant for simulations launched as part of a batch (e.g. from an automatic launching procedure)
        self.add_property("timing_table_path", "file_path", "path of the timing table", None)
        self.add_property("memory_table_path", "file_path", "path of the memory table", None)

        # Properties relevant for simulations part of a scaling test
        self.add_property("scaling_path", "string", "scaling directory path", None)
        self.add_property("scaling_run_name", "string", "name of scaling run", None)

        # Properties relevant for simulations part of radiative transfer modeling
        self.add_property("modeling_path", "string", "modeling directory path", None)

        # Set options
        self.set_options(kwargs)

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

    def check(self, logging_options=None, output_path=None, retrieve_types=None):

        """
        This function ...
        :param logging_options:
        :param output_path:
        :param retrieve_types:
        :return:
        """

        # Inform the user
        log.info("Checking the analysis options ...")

        # MISC
        self.check_misc(output_path=output_path, retrieve_types=retrieve_types)

        # PLOTTING
        self.check_plotting(output_path=output_path, retrieve_types=retrieve_types)

        # EXTRACTION
        self.check_extraction(output_path=output_path, retrieve_types=retrieve_types)

        # Check logging options
        self.check_logging(logging_options)

    # -----------------------------------------------------------------

    def check_misc(self, output_path=None, retrieve_types=None):

        """
        This function ...
        :param output_path:
        :param retrieve_types:
        :return:
        """

        # Debugging
        log.debug("Checking miscellaneous settings ...")

        from ..simulation.output import output_types as ot

        # If any misc setting has been enabled, check whether the misc path has been set
        if self.any_misc and self.misc.path is None:
            if output_path is None: raise ValueError("The misc output path has not been set")
            else:
                log.warning("Misc output will be written to " + output_path)
                self.misc.path = output_path

        # Adjust retrieve types
        if retrieve_types is not None:

            if self.misc.rgb and not (ot.images in retrieve_types or ot.total_images in retrieve_types):
                log.warning("Creating RGB images is enabled so total datacube retrieval will also be enabled")
                retrieve_types.append(ot.total_images)

            if self.misc.wave and not (ot.images in retrieve_types or ot.total_images in retrieve_types):
                log.warning("Creating wave movies is enabled so total datacube retrieval will also be enabled")
                retrieve_types.append(ot.total_images)

            if self.misc.fluxes and ot.seds not in retrieve_types:
                log.warning("Calculating observed fluxes is enabled so SED retrieval will also be enabled")
                retrieve_types.append(ot.seds)

            if self.misc.images and not (ot.images in retrieve_types or ot.total_images in retrieve_types):
                log.warning("Creating observed images is enabled so total datacube retrieval will also be enabled")
                retrieve_types.append(ot.images)

    # -----------------------------------------------------------------

    def check_plotting(self, output_path=None, retrieve_types=None):

        """
        This function ...
        :param output_path:
        :param retrieve_types:
        :return:
        """

        # Debugging
        log.debug("Checking plotting settings ...")

        from ..simulation.output import output_types as ot

        # If any plotting setting has been enabled, check whether the plotting path has been set
        if self.any_plotting and self.plotting.path is None:
            if output_path is None: raise ValueError("The plotting path has not been set")
            else:
                log.warning("Plots will be saved to " + output_path)
                self.plotting.path = output_path

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

        # Adjust retrieve types
        if retrieve_types is not None:

            # If SED plotting has been enabled, enable SED retrieval
            if self.plotting.seds and ot.seds not in retrieve_types:
                log.warning("SED plotting is enabled so SED file retrieval will also be enabled")
                retrieve_types.append(ot.seds)

            # If grid plotting has been enabled:
            if self.plotting.grids and ot.grid not in retrieve_types:
                log.warning("Grid plotting is enabled so grid data retrieval will also be enabled")
                retrieve_types.append(ot.grid)

    # -----------------------------------------------------------------

    def check_extraction(self, output_path=None, retrieve_types=None):

        """
        This function ...
        :param output_path:
        :param retrieve_types:
        :return:
        """

        # Debugging
        log.debug("Checking extraction settings ...")

        from ..simulation.output import output_types as ot

        # If any extraction setting has been enabled, check whether the extraction path has been set
        if self.any_extraction and self.extraction.path is None:
            if output_path is None: raise ValueError("The extraction path has not been set")
            else:
                log.warning("Extraction data will be placed in " + output_path)
                self.extraction.path = output_path

        # Adjust retrieve types
        if retrieve_types is not None:

            # If progress plotting has been enabled, enable log file retrieval
            if self.extraction.progress and ot.logfiles not in retrieve_types:
                log.warning("Progress extraction is enabled so log file retrieval will also be enabled")
                retrieve_types.append(ot.logfiles)

            # If memory extraction has been enabled, enable log file retrieval
            if self.extraction.memory and ot.logfiles not in retrieve_types:
                log.warning("Memory extraction is enabled so log file retrieval will also be enabled")
                retrieve_types.append(ot.logfiles)

            # If timeline extraction has been enabled, enable log file retrieval
            if self.extraction.timeline and ot.logfiles not in retrieve_types:
                log.warning("Timeline extraction is enabled so log file retrieval will also be enabled")
                retrieve_types.append(ot.logfiles)

    # -----------------------------------------------------------------

    def check_logging(self, logging_options):

        """
        This function ...
        :param logging_options:
        :return:
        """

        # Debugging
        log.debug("Checking logging options ...")

        # Check the logging options, and adapt if necessary
        if logging_options is not None:

            # If memory extraction has been enabled, enable memory logging
            if self.extraction.memory and not logging_options.memory:
                log.warning("Memory extraction is enabled so memory logging will also be enabled")
                logging_options.memory = True

        # Logging options are not passed
        elif self.extraction.memory: log.warning("Memory extraction is enabled but the logging options could not be verified")

# -----------------------------------------------------------------
