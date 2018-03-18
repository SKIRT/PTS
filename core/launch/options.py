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
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.composite import SimplePropertyComposite
from ..basics.plot import plotting_libraries, mpl, plotting_formats, pdf

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

    def __init__(self, **kwargs):

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
        self.add_property("walltime", "real", "expected walltime", None) # in seconds
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
animations_name = "animations"
fluxes_name = "fluxes"
fluxes_from_images_name = "fluxes_from_images"
images_name = "images"

# -----------------------------------------------------------------

extraction_names = [progress_name, timeline_name, memory_name]
plotting_names = [progress_name, timeline_name, memory_name, seds_name, grids_name]
misc_names = [rgb_name, animations_name, fluxes_name, fluxes_from_images_name, images_name]

# -----------------------------------------------------------------

def get_analysis_property_names():

    """
    This function ...
    :return:
    """

    return AnalysisOptions().property_names

# -----------------------------------------------------------------

def get_analysis_property_names_and_descriptions():

    """
    This function ...
    :return:
    """

    descriptions = OrderedDict()
    opts = AnalysisOptions()
    for name in opts.property_names:
        description = opts.description_for_property(name)
        descriptions[name] = description
    return descriptions

# -----------------------------------------------------------------

def get_analysis_section_names():

    """
    Thisf unction ...
    :return:
    """

    return AnalysisOptions().section_names

# -----------------------------------------------------------------

def get_analysis_section_names_and_descriptions():

    """
    This function ...
    :return:
    """

    descriptions = OrderedDict()
    opts = AnalysisOptions()
    for name in opts.section_names:
        description = opts.description_for_property(name)
        descriptions[name] = description
    return descriptions

# -----------------------------------------------------------------

def get_analysis_property_names_for_section(section_name):

    """
    This function ...
    :param section_name:
    :return:
    """

    return AnalysisOptions()[section_name].property_names

# -----------------------------------------------------------------

def get_analysis_property_names_and_descriptions_for_section(section_name):

    """
    This function ...
    :param section_name:
    :return:
    """

    descriptions = OrderedDict()
    opts = AnalysisOptions()
    section = opts[section_name]
    for name in section.property_names:
        description = section.description_for_property(name)
        descriptions[name] = description
    return descriptions

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
        self.add_section("extraction", "options for extracting data from the simulation's log files")
        self.extraction.add_property("path", "string", "extraction directory", None)
        self.extraction.add_property("progress", "boolean", "extract information about the progress in the different simulation phases", False)
        self.extraction.add_property("timeline", "boolean", "extract timeline information for the different simulation phases on the different processes", False)
        self.extraction.add_property("memory", "boolean", "extract information about the memory usage during the simulation", False)

        # Plotting
        self.add_section("plotting", "options for plotting simulation output")
        self.plotting.add_property("path", "string", "plotting directory", None)
        self.plotting.add_property("format", "string", "image format for the plots", pdf, choices=plotting_formats)
        self.plotting.add_property("progress", "boolean", "make plots of the progress of the simulation phases as a function of time", False)
        self.plotting.add_property("timeline", "boolean", "plot the timeline for the different processes", False)
        self.plotting.add_property("memory", "boolean", "plot the memory consumption as a function of time", False)
        self.plotting.add_property("seds", "boolean", "make plots of the simulated SEDs", False)
        self.plotting.add_property("grids", "boolean", "make plots of the dust grid", False)
        self.plotting.add_property("reference_seds", "filepath_list", "paths to reference SED files against which the simulated SKIRT SEDs should be plotted", None)
        self.plotting.add_property("ignore_filters", "filter_list", "filters to ignore for the plotting", [])
        self.plotting.add_property("library", "string", "plotting library", default_value=mpl, choices=plotting_libraries)

        # Misc
        ## General
        self.add_section("misc", "settings for creating data of various types from the simulation output")
        self.misc.add_property("path", "string", "misc output directory", None)

        ## RGB images
        self.misc.add_property("rgb", "boolean", "make RGB images from the simulated datacube(s)", False)

        ## Animation
        self.misc.add_property("animations", "boolean", "make an animation of the images of the simulated datacube(s)", False)
        self.misc.add_property("write_animation_frames", "boolean", "write the frames of the created animations separately", False)

        ## Observed fluxes and images
        self.misc.add_property("observation_filters", "string_list", "the names of the filters for which to recreate the observations", None)
        self.misc.add_property("observation_instruments", "string_list", "the names of the instruments for which to recreate the observations", None)

        ## Observed fluxes
        self.misc.add_property("fluxes", "boolean", "calculate observed fluxes from the SKIRT output SEDs", False)
        self.misc.add_property("flux_errors", "string_string_dictionary", "errorbars for the different flux points of the mock observed SED")
        self.misc.add_property("fluxes_from_images", "boolean", "calculate observed fluxes from the SKIRT output datacubes", False)
        self.misc.add_property("fluxes_from_images_instrument", "string", "calculate observed fluxes from images of this instrument", None)
        self.misc.add_property("fluxes_from_images_wcs", "file_path", "file containing the coordinate system of the instrument's coordinate system (necessary when masks are specified)")
        self.misc.add_property("fluxes_from_images_errors", "string_string_dictionary", "errorbars for the different flux points of the mock observed SED")
        self.misc.add_property("fluxes_from_images_masks", "string_filepath_dictionary", "filepaths of the FITS files that contain the image masks (with coordinate system!), keys must be filter names")
        self.misc.add_property("fluxes_from_images_mask_from_nans", "boolean", "load masks as the NaN pixels in the specified plane (or primary frame if not specified)")
        self.misc.add_property("fluxes_from_images_mask_plane", "string", "name of the image mask plane")
        self.misc.add_property("write_fluxes_images", "boolean", "write out the images created to calculate the fluxes", False)
        self.misc.add_property("plot_fluxes", "boolean", "plot the fluxes", False)
        self.misc.add_property("plot_fluxes_from_images", "boolean", "plot the fluxes from images", False)
        self.misc.add_property("plot_fluxes_reference_seds", "string_filepath_dictionary", "reference SEDs for plotting fluxes")
        self.misc.add_property("plot_fluxes_from_images_reference_seds", "string_filepath_dictionary", "reference SEDs for plotting fluxes")
        self.misc.add_property("plot_fluxes_images", "boolean", "plot the images created to calculate the fluxes", False)

        ## Observed fluxes from images REMOTE
        self.misc.add_property("fluxes_from_images_remote", "string", "perform the creation of the observed images for fluxes on a remote host", None)
        self.misc.add_property("fluxes_from_images_remote_spectral_convolution", "boolean", "perform spectral convolution with the datacubes for observed fluxes from images remotely", False)
        self.misc.add_property("fluxes_from_images_remote_threshold", "data_quantity", "file size threshold for working with remote datacubes for fluxes from images", "2 GB", convert_default=True)
        self.misc.add_property("fluxes_from_images_remote_npixels_threshold", "positive_integer", "threshold for working remotely of the number of pixels (nx * ny) of the datacubes", 1e4)
        self.misc.add_property("fluxes_from_images_rebin_remote_threshold", "data_quantity", "data size threshold for remote rebinning of images for fluxes", "0.5 GB", convert_default=True)

        ## Observed images
        self.misc.add_property("images", "boolean", "make observed images form the simulated datacube(s)", False)
        self.misc.add_property("wcs_instrument", "string", "instrument for which to take the images_wcs as the WCS (if not a dictionary)")
        self.misc.add_property("images_wcs", "filepath_or_string_filepath_dictionary", "path to the FITS/txt file for which the WCS should be set as the WCS of the recreated observed images (single path or path for each datacube as a dictionary)", None)
        self.misc.add_property("images_unit", "string", "the unit to which the recreated observed images should be converted", None)
        self.misc.add_property("write_intermediate_images", "boolean", "write intermediate results from the observed image making procedure", False)
        self.misc.add_property("write_convolution_kernels", "boolean", "write the convolution kernels used in the observed image making procedure", False)
        self.misc.add_property("no_images_filters", "filter_list", "don't create observed images for these filters", [])
        self.misc.add_property("plot_images", "boolean", "plot the mock observed images", False)

        ## CONVOLUTION
        self.misc.add_property("images_kernels", "string_string_dictionary", "paths to the FITS file of convolution kernel used for convolving the observed images (a dictionary where the keys are the filter names)", None)
        self.misc.add_property("images_psfs_auto", "boolean", "automatically determine the appropriate PSF kernel for each image", False)
        self.misc.add_property("fwhms_dataset", "file_path", "path to the dataset that should be used to get the target FWHM for the different images", None)

        # Remote thresholds
        self.misc.add_property("make_images_remote", "string", "perform the calculation of the observed images on a remote host (this is a memory and CPU intensive step)", None)
        self.misc.add_property("remote_spectral_convolution", "boolean", "perform spectral convolution with the datacube remotely", False)
        self.misc.add_property("images_remote_threshold", "data_quantity", "file size threshold for working with remote datacubes", "2 GB", convert_default=True)
        self.misc.add_property("images_remote_npixels_threshold", "positive_integer", "threshold for working remotely of the number of pixels (nx * ny) of the datacubes", 1e4)
        self.misc.add_property("rebin_remote_threshold", "data_quantity", "data size threshold for remote rebinning", "0.5 GB", convert_default=True)
        self.misc.add_property("convolve_remote_threshold", "data_quantity", "data size threshold for remote convolution", "1. GB", convert_default=True)

        ## REBINNING
        self.misc.add_property("rebin_instrument", "string", "instrument for which the rebin_wcs or rebin_dataset should be used")
        self.misc.add_property("rebin_wcs", "filepath_or_string_filepath_dictionary", "paths to the FITS/txt files of which the WCS should be used as the target for rebinning")
        self.misc.add_property("rebin_dataset", "file_path", "path to the dataset of which the coordinate systems are to be used as rebinning references")

        ## OTHER
        self.misc.add_property("images_nprocesses_local", "positive_integer", "number of parallel processes to use when creating the observed images (local execution)", 2)
        self.misc.add_property("images_nprocesses_remote", "positive_integer", "number of parallel processes to use when creating the observed images (remote execution)", 8)
        self.misc.add_property("group_images", "boolean", "group the images per instrument", False)
        self.misc.add_property("images_spectral_convolution", "boolean", "use spectral convolution to create observed images", True)
        self.misc.add_property("fluxes_spectral_convolution", "boolean", "use spectral convolution to calculate observed fluxes", True)
        self.misc.add_property("fluxes_from_images_spectral_convolution", "boolean", "use spectral convolution to calculate observed fluxes from the images", False)
        self.misc.add_property("no_images_spectral_convolution_filters", "filter_list", "don't spectrally convolve to create the observed images for these filters", [])
        self.misc.add_property("no_fluxes_spectral_convolution_filters", "filter_list", "don't spectrally convolve to calculate the observed fluxes for these filters", [])
        self.misc.add_property("no_fluxes_from_images_spectral_convolution_filters", "filter_list", "don't spectrally convolve to calculate the observed fluxes from images for these filters", [])

        # Properties that are relevant for simulations launched as part of a batch (e.g. from an automatic launching procedure)
        self.add_property("timing_table_path", "file_path", "path of the timing table", None)
        self.add_property("memory_table_path", "file_path", "path of the memory table", None)

        # Properties relevant for simulations part of a scaling test
        self.add_property("scaling_path", "string", "scaling directory path", None)
        self.add_property("scaling_run_name", "string", "name of scaling run", None)

        # Properties relevant for simulations part of radiative transfer modeling
        self.add_property("modeling_path", "string", "modeling directory path", None)

        # The paths to the extra simulation analysers
        self.add_property("analyser_paths", "string_list", "paths to the extra simulation analysers", [])

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

        return self.misc.rgb or self.misc.animations or self.misc.fluxes or self.misc.fluxes_from_images or self.misc.images

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
                log.warning("Making RGB images is enabled so total datacube retrieval will also be enabled")
                retrieve_types.append(ot.total_images)

            # if self.misc.wave and not (ot.images in retrieve_types or ot.total_images in retrieve_types):
            #     log.warning("Creating wave movies is enabled so total datacube retrieval will also be enabled")
            #     retrieve_types.append(ot.total_images)

            if self.misc.animations and not (ot.images in retrieve_types or ot.total_images in retrieve_types):
                log.warning("Making datacube animations is enabled so total datacube retrieval will also be enabled")
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
