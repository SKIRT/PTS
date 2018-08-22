#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.analyser Contains the BasicAnalyser class, used for analysing simulation output.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..extract.progress import ProgressExtractor, NoProgressData
from ..extract.timeline import TimeLineExtractor, NoTimingData
from ..extract.memory import MemoryExtractor, NoMemoryData
from ..plot.progress import ProgressPlotter
from ..plot.timeline import TimeLinePlotter
from ..plot.memory import MemoryPlotter
from ..plot.grids import plotgrids
from ..misc.rgb import RGBImageMaker
from ..misc.animations import DataCubeAnimationsMaker
from ..misc.fluxes import ObservedFluxCalculator
from ..misc.images import ObservedImageMaker
from ..basics.log import log
from ..tools import filesystem as fs
from ..plot.simulationseds import SimulationSEDPlotter
from ..tools import types
from .options import progress_name, timeline_name, memory_name, seds_name, grids_name, rgb_name, animations_name, fluxes_name, fluxes_from_images_name, images_name
from ..extract.progress import ProgressTable
from ..extract.timeline import TimeLineTable
from ..extract.memory import MemoryUsageTable
from ..tools import sequences
from ..tools.utils import lazyproperty
from ...magic.basics.coordinatesystem import CoordinateSystem
from ..filter.filter import parse_filter
from ...magic.core.mask import Mask
from ..simulation.remote import get_simulation_for_host

# -----------------------------------------------------------------

extraction = "extraction"
plotting = "plotting"
misc = "misc"
steps = [extraction, plotting, misc]

# -----------------------------------------------------------------

progress_filename = "progress.dat"
timeline_filename = "timeline.dat"
memory_filename = "memory.dat"

# -----------------------------------------------------------------

class BasicAnalyser(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(BasicAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

        # The analysis options
        self.extraction_options = None
        self.plotting_options = None
        self.misc_options = None

        # The tables with extracted information
        self.progress = None
        self.timeline = None
        self.memory = None

        # The RGB and animations makers
        self.rgb_maker = None
        self.animations_maker = None

        # The flux calculators and image maker
        self.flux_calculator = None
        self.image_flux_calculator = None
        self.image_maker = None

        # The observed fluxes and images
        self._mock_seds = None
        self._mock_seds_from_images = None
        self._mock_images = None

    # -----------------------------------------------------------------

    @property
    def extraction(self):

        """
        This function ...
        :return:
        """

        return self.config.extract and self.simulation.analysis.any_extraction

    # -----------------------------------------------------------------

    @property
    def plotting(self):

        """
        This function ...
        :return:
        """

        return self.config.plot and self.simulation.analysis.any_plotting

    # -----------------------------------------------------------------

    @property
    def miscellaneous(self):

        """
        This function ...
        :return:
        """

        return self.config.misc and self.simulation.analysis.any_misc

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

        # 2. Extract information from the simulation's log files
        if self.extraction: self.extract()

        # 3. Make plots based on the simulation output
        if self.plotting: self.plot()

        # 4. Miscellaneous output
        if self.miscellaneous: self.misc()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BasicAnalyser, self).setup()

        # Make a local reference to the simulation object
        if "simulation" in kwargs: self.simulation = kwargs.pop("simulation")
        elif self.config.remote is not None and self.config.id is not None: self.load_simulation()
        else: raise ValueError("No simulation is specified")

        # Update the analysis options for added features
        self.simulation.update_analysis_options()

        # Also make references to the simulation's analysis options for extraction, plotting and misc (for shorter notation)
        self.extraction_options = self.simulation.analysis.extraction
        self.plotting_options = self.simulation.analysis.plotting
        self.misc_options = self.simulation.analysis.misc

        # Make paths if necessary
        if self.extraction and not fs.is_directory(self.extraction_path): fs.create_directory(self.extraction_path, recursive=True)
        if self.plotting and not fs.is_directory(self.plotting_path): fs.create_directory(self.plotting_path, recursive=True)
        if self.miscellaneous and not fs.is_directory(self.misc_path): fs.create_directory(self.misc_path, recursive=True)

    # -----------------------------------------------------------------

    @property
    def extraction_path(self):

        """
        This function ...
        :return:
        """

        return self.extraction_options.path

    # -----------------------------------------------------------------

    @property
    def plotting_path(self):

        """
        This function ...
        :return:
        """

        return self.plotting_options.path

    # -----------------------------------------------------------------

    @property
    def misc_path(self):

        """
        This function ...
        :return:
        """

        return self.misc_options.path

    # -----------------------------------------------------------------

    def load_simulation(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the simulation ...")

        # Load simulation
        self.simulation = get_simulation_for_host(self.config.remote, self.config.id)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the simulation to None
        self.simulation = None

        # Set the options to None
        self.extraction_options = None
        self.plotting_options = None
        self.misc_options = None

        # Clear the extractors
        self.progress = None
        self.timeline = None
        self.memory = None

        # Clear the flux calculator and image maker
        self.flux_calculator = None
        self.image_maker = None

    # -----------------------------------------------------------------

    @property
    def extracted_progress(self):

        """
        This fucntion ...
        :return:
        """

        return progress_name in self.simulation.analysed_extraction

    # -----------------------------------------------------------------

    @property
    def extracted_timeline(self):

        """
        This function ...
        :return:
        """

        return timeline_name in self.simulation.analysed_extraction

    # -----------------------------------------------------------------

    @property
    def extracted_memory(self):

        """
        This function ...
        :return:
        """

        return memory_name in self.simulation.analysed_extraction

    # -----------------------------------------------------------------

    @property
    def progress_extraction(self):

        """
        This function ...
        :return:
        """

        return self.extraction_options.progress

    # -----------------------------------------------------------------

    @property
    def timeline_extraction(self):

        """
        Thisn function ...
        :return:
        """

        return self.extraction_options.timeline

    # -----------------------------------------------------------------

    @property
    def memory_extraction(self):

        """
        This function ...
        :return:
        """

        return self.extraction_options.memory

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting ...")

        # Show options
        self.show_extraction_options()

        # Progress
        if self.progress_extraction: self.extract_or_load_progress()

        # Timeline
        if self.timeline_extraction: self.extract_or_load_timeline()

        # Memory
        if self.memory_extraction: self.extract_or_load_memory()

    # -----------------------------------------------------------------

    def show_extraction_options(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Extraction options:")
        if log.is_debug: print(str(self.extraction_options))

    # -----------------------------------------------------------------

    def extract_or_load_progress(self):

        """
        This function ...
        :return:
        """

        # Extract
        if not self.extracted_progress: self.extract_progress()

        # Load
        else: self.load_progress()

    # -----------------------------------------------------------------

    def extract_or_load_timeline(self):

        """
        This function ...
        :return:
        """

        # Extract
        if not self.extracted_timeline: self.extract_timeline()

        # Load
        else: self.load_timeline()

    # -----------------------------------------------------------------

    def extract_or_load_memory(self):

        """
        This function ...
        :return:
        """

        # Extract
        if not self.extracted_memory: self.extract_memory()

        # Load
        else: self.load_memory()

    # -----------------------------------------------------------------

    def load_progress(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the progress table ...")

        # Determine the path to the progress file
        path = fs.join(self.extraction_options.path, progress_filename)

        # Load the table
        self.progress = ProgressTable.from_file(path)

    # -----------------------------------------------------------------

    def load_timeline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the timeline table ...")

        # Determine the path to the timeline file
        path = fs.join(self.extraction_options.path, timeline_filename)

        # Load the table
        self.timeline = TimeLineTable.from_file(path)

    # -----------------------------------------------------------------

    def load_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the memory table ...")

        # Determine the path to the memory file
        path = fs.join(self.extraction_options.path, memory_filename)

        # Load the table
        self.memory = MemoryUsageTable.from_file(path)

    # -----------------------------------------------------------------

    @property
    def plotted_seds(self):

        """
        This function ...
        :return:
        """

        return seds_name in self.simulation.analysed_plotting

    # -----------------------------------------------------------------

    @property
    def plotted_grids(self):

        """
        Thins function ...
        :return:
        """

        return grids_name in self.simulation.analysed_plotting

    # -----------------------------------------------------------------

    @property
    def plotted_progress(self):

        """
        This function ...
        :return:
        """

        return progress_name in self.simulation.analysed_plotting

    # -----------------------------------------------------------------

    @property
    def plotted_timeline(self):

        """
        This function ...
        :return:
        """

        return timeline_name in self.simulation.analysed_plotting

    # -----------------------------------------------------------------

    @property
    def plotted_memory(self):

        """
        This function ...
        :return:
        """

        return memory_name in self.simulation.analysed_plotting

    # -----------------------------------------------------------------

    @property
    def seds_plotting(self):

        """
        Thi function ...
        :return:
        """

        return self.plotting_options.seds

    # -----------------------------------------------------------------

    @property
    def needs_seds_plotting(self):

        """
        This function ...
        :return:
        """

        return self.seds_plotting and not self.plotted_seds

    # -----------------------------------------------------------------

    @property
    def grids_plotting(self):

        """
        This function ...
        :return:
        """

        return self.plotting_options.grids

    # -----------------------------------------------------------------

    @property
    def needs_grids_plotting(self):

        """
        This function ...
        :return:
        """

        return self.grids_plotting and not self.plotted_grids

    # -----------------------------------------------------------------

    @property
    def progress_plotting(self):

        """
        This function ...
        :return:
        """

        return self.plotting_options.progress

    # -----------------------------------------------------------------

    @property
    def needs_progress_plotting(self):

        """
        Thisn function ...
        :return:
        """

        return self.progress_plotting and not self.plotted_progress

    # -----------------------------------------------------------------

    @property
    def timeline_plotting(self):

        """
        Thisfunction ...
        :return:
        """

        return self.plotting_options.timeline

    # -----------------------------------------------------------------

    @property
    def needs_timeline_plotting(self):

        """
        This function ...
        :return:
        """

        return self.timeline_plotting and not self.plotted_timeline

    # -----------------------------------------------------------------

    @property
    def memory_plotting(self):

        """
        This function ...
        :return:
        """

        return self.plotting_options.memory

    # -----------------------------------------------------------------

    @property
    def needs_memory_plotting(self):

        """
        This function ...
        :return:
        """

        return self.memory_plotting and not self.plotted_memory

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Show
        self.show_plot_options()

        # If requested, plot the SED's
        if self.needs_seds_plotting: self.plot_seds()

        # If requested, make plots of the dust grid
        if self.needs_grids_plotting: self.plot_grids()

        # If requested, plot the simulation progress as a function of time
        if self.needs_progress_plotting: self.plot_progress()

        # If requested, plot a timeline of the different simulation phases
        if self.needs_timeline_plotting: self.plot_timeline()

        # If requested, plot the memory usage as a function of time
        if self.needs_memory_plotting: self.plot_memory()

    # -----------------------------------------------------------------

    def show_plot_options(self):

        """
        This ufnction ...
        :return:
        """

        # Debugging
        log.debug("Plotting options:")
        if log.is_debug: print(str(self.plotting_options))

    # -----------------------------------------------------------------

    @property
    def has_rgb(self):

        """
        This function ...
        :return:
        """

        return rgb_name in self.simulation.analysed_misc

    # -----------------------------------------------------------------

    @property
    def has_animations(self):

        """
        This fucntion ...
        :return:
        """

        return animations_name in self.simulation.analysed_misc

    # -----------------------------------------------------------------

    @property
    def has_fluxes(self):

        """
        This function ...
        :return:
        """

        return fluxes_name in self.simulation.analysed_misc

    # -----------------------------------------------------------------

    @property
    def has_fluxes_from_images(self):

        """
        This function ...
        :return:
        """

        return fluxes_from_images_name in self.simulation.analysed_misc

    # -----------------------------------------------------------------

    @property
    def has_images(self):

        """
        This function ...
        :return:
        """

        return images_name in self.simulation.analysed_misc

    # -----------------------------------------------------------------

    @property
    def rgb(self):

        """
        This function ...
        :return:
        """

        return self.misc_options.rgb

    # -----------------------------------------------------------------

    @property
    def animations(self):

        """
        This function ...
        :return:
        """

        return self.misc_options.animations

    # -----------------------------------------------------------------

    @property
    def observed_fluxes(self):

        """
        This function ...
        :return:
        """

        return self.misc_options.fluxes

    # -----------------------------------------------------------------

    @property
    def observed_fluxes_from_images(self):

        """
        This function ...
        :return:
        """

        return self.misc_options.fluxes_from_images

    # -----------------------------------------------------------------

    @property
    def observed_images(self):

        """
        This function ...
        :return:
        """

        return self.misc_options.images

    # -----------------------------------------------------------------

    def misc(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing miscellaneous analysis ...")

        # Show
        self.show_misc_options()

        # If requested, make RGB images of the output FITS files
        if self.rgb and not self.has_rgb: self.make_rgb()

        # If requested, make datacube animations from the output datacubes
        if self.animations and not self.has_animations: self.make_animations()

        # Observed fluxes
        if self.observed_fluxes: self.get_observed_fluxes()

        # Observed fluxes from images
        if self.observed_fluxes_from_images: self.get_observed_fluxes_from_images()

        # Observed images
        if self.observed_images: self.get_observed_images()

    # -----------------------------------------------------------------

    def show_misc_options(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Miscellaneous options:")
        if log.is_debug: print(str(self.misc_options))

    # -----------------------------------------------------------------

    def extract_progress(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the progress information ...")

        # Create a ProgressExtractor instance
        extractor = ProgressExtractor()

        # Determine the path to the progress file
        path = fs.join(self.extraction_options.path, progress_filename)

        # Run the progress extractor, get the progress table
        try:

            # Try extraction progress information
            extractor.run(self.simulation, path)
            self.progress = extractor.table

            # Done
            self.simulation.analysed_extraction.append(progress_name)
            self.simulation.save()

        # No progress data
        except NoProgressData:

            # Ignore this, but unset extract progress flag
            if self.config.ignore_missing_data:

                # Give warning
                log.warning("Missing progress data: disabling progress extraction option ...")

                # Unset progress extraction and plotting
                self.simulation.analysis.extraction.progress = False
                self.simulation.analysis.plotting.progress = False
                self.simulation.save()

            # Give error
            else: raise RuntimeError("Could not extract progress information")

    # -----------------------------------------------------------------

    def extract_timeline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the timeline information ...")

        # Create a TimeLineExtractor instance
        extractor = TimeLineExtractor()

        # Determine the path to the timeline file
        path = fs.join(self.extraction_options.path, timeline_filename)

        # Run the timeline extractor, get the timeline table
        try:

            # Try extracting timeline information
            extractor.run(simulation=self.simulation, output_path=path)
            self.timeline = extractor.table

            # Done
            self.simulation.analysed_extraction.append(timeline_name)
            self.simulation.save()

        # No timing data
        except NoTimingData:

            # Ignore this, but unset extract timeline flag
            if self.config.ignore_missing_data:

                # Give warning
                log.warning("Missing timing data: disabling timeline extraction option ...")

                # Unset timeline extraction and plotting
                self.simulation.analysis.extraction.timeline = False
                self.simulation.analysis.plotting.timeline = False
                self.simulation.save()

            # Give error
            else: raise RuntimeError("Could not extract timeline information")

    # -----------------------------------------------------------------

    def extract_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the memory information ...")

        # Create a MemoryExtractor instance
        extractor = MemoryExtractor()

        # Determine the path to the memory file
        path = fs.join(self.extraction_options.path, memory_filename)

        # Run the memory extractor, get the memory usage table
        try:

            # Try extracting memory information
            extractor.run(self.simulation, path)
            self.memory = extractor.table

            # Done
            self.simulation.analysed_extraction.append(memory_name)
            self.simulation.save()

        # No memory data
        except NoMemoryData:

            # Ignore this, but unset extract memory flag
            if self.config.ignore_missing_data:

                # Give warning
                log.warning("Missing memory data: disabling memory extraction option ...")

                # Unset memory extraction and plotting
                self.simulation.analysis.extraction.memory = False
                self.simulation.analysis.plotting.memory = False
                self.simulation.save()

            # Give error
            else: raise RuntimeError("Could not extract memory information")

    # -----------------------------------------------------------------

    def plot_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting SEDs ...")

        # Create a simulation SED plotter instance
        plotter = SimulationSEDPlotter()

        # Set the output directory path
        plotter.config.output = self.plotting_options.path

        # Set plotting options
        plotter.config.format = self.plotting_options.format
        plotter.config.library = self.plotting_options.library

        # Set ignore filters
        plotter.config.ignore_filters = self.plotting_options.ignore_filters

        # Set reference SED paths
        plotter.config.reference_seds = self.plotting_options.reference_seds

        # Run the plotter
        plotter.run(simulation=self.simulation)

        # Done
        self.simulation.analysed_plotting.append(seds_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    def plot_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting grids ...")

        # Plot the dust grid for the simulation
        plotgrids(self.simulation, output_path=self.plotting_options.path, silent=(not log.is_debug))

        # Done
        self.simulation.analysed_plotting.append(grids_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    def plot_progress(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the progress information ...")

        # Create a ProgressPlotter object
        plotter = ProgressPlotter()

        # Run the progress plotter
        plotter.run(self.progress, self.plotting_options.path)

        # Done
        self.simulation.analysed_plotting.append(progress_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    def plot_timeline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the timeline ...")

        # Create a TimeLinePlotter object
        plotter = TimeLinePlotter()

        # Set the output path
        plotter.config.output = self.plotting_options.path

        # Run the timeline plotter
        plotter.run(timeline=self.timeline)

        # Done
        self.simulation.analysed_plotting.append(timeline_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    def plot_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the memory information ...")

        # Create a MemoryPlotter object
        plotter = MemoryPlotter()

        # Run the memory plotter
        plotter.run(self.memory, self.plotting_options.path)

        # Done
        self.simulation.analysed_plotting.append(memory_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    @lazyproperty
    def rgb_output_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.misc_options.path, "rgb", recursive=True)

    # -----------------------------------------------------------------

    def make_rgb(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making RGB images ...")

        # Create the RGBImageMaker
        self.rgb_maker = RGBImageMaker()

        # Run the maker
        self.rgb_maker.run(simulation=self.simulation, output_path=self.rgb_output_path)

        # Done
        self.simulation.analysed_misc.append(rgb_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    @lazyproperty
    def animations_output_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.misc_options.path, "animations", recursive=True)

    # -----------------------------------------------------------------

    def make_animations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making datacube animations ...")

        # Create the animations maker
        self.animations_maker = DataCubeAnimationsMaker()

        # Set options
        self.animations_maker.config.write_frames = self.misc_options.write_animation_frames

        # Run the maker
        self.animations_maker.run(simulation=self.simulation, output_path=self.animations_output_path)

        # Done
        self.simulation.analysed_misc.append(animations_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    @property
    def filters_for_fluxes(self):

        """
        Thisf unction ...
        :return:
        """

        return self.misc_options.observation_filters

    # -----------------------------------------------------------------

    @lazyproperty
    def fluxes_output_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.misc_options.path, "fluxes", recursive=True)

    # -----------------------------------------------------------------

    # @lazyproperty
    # def fluxes_reference_sed(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     from ..data.sed import ObservedSED
    #     if self.misc_options.plot_fluxes_reference_sed is None: return None
    #     else: return ObservedSED.from_file(self.misc_options.plot_fluxes_reference_sed)

    # -----------------------------------------------------------------

    @lazyproperty
    def fluxes_reference_seds(self):

        """
        This function ...
        :return:
        """

        from ..data.sed import ObservedSED
        if self.misc_options.plot_fluxes_reference_seds is None: return None
        else:
            seds = OrderedDict()
            for label in self.misc_options.plot_fluxes_reference_seds: seds[label] = ObservedSED.from_file(self.misc_options.plot_fluxes_reference_seds[label])
            return seds

    # -----------------------------------------------------------------

    def get_observed_fluxes(self):

        """
        This function ...
        :return:
        """

        # Create
        if not self.has_fluxes: self.calculate_observed_fluxes()

        # Load
        else: self.load_observed_fluxes()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fluxes_input(self):

        """
        This function ...
        :return:
        """

        # Create input dictionary
        input_dict = dict()

        # Set input
        input_dict["simulation"] = self.simulation
        input_dict["output_path"] = self.fluxes_output_path
        input_dict["filter_names"] = self.filters_for_fluxes
        input_dict["instrument_names"] = self.misc_options.observation_instruments
        input_dict["errors"] = self.misc_options.flux_errors
        input_dict["no_spectral_convolution_filters"] = self.misc_options.no_fluxes_spectral_convolution_filters
        input_dict["reference_seds"] = self.fluxes_reference_seds

        # Return the input
        return input_dict

    # -----------------------------------------------------------------

    def calculate_observed_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the observed fluxes ...")

        # Create a ObservedFluxCalculator object
        self.flux_calculator = ObservedFluxCalculator()

        # Set spectral convolution flag
        self.flux_calculator.config.spectral_convolution = self.misc_options.fluxes_spectral_convolution

        # Set plot flag
        self.flux_calculator.config.plot = True
        self.flux_calculator.config.plot_seds = self.misc_options.plot_fluxes
        self.flux_calculator.config.plot_images = False

        # EXTRA OPTIONS
        self.flux_calculator.config.check_wavelengths = self.config.check_wavelengths
        self.flux_calculator.config.ignore_bad = self.config.ignore_bad
        self.flux_calculator.config.skip_ignored_bad_convolution = self.config.skip_ignored_bad_convolution
        self.flux_calculator.config.skip_ignored_bad_closest = self.config.skip_ignored_bad_closest

        # DEPLOYMENT
        self.flux_calculator.config.deploy_pts = self.config.deploy_pts
        self.flux_calculator.config.update_dependencies = self.config.update_dependencies

        # Run
        self.flux_calculator.run(**self.observed_fluxes_input)

        # Done
        self.simulation.analysed_misc.append(fluxes_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    def load_observed_fluxes(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the observed fluxes ...")

        from ..data.sed import ObservedSED

        # Get output path
        output_path = self.fluxes_output_path

        # Initialize dictionary
        self._mock_seds = OrderedDict()

        # Loop over the instruments for which observed fluxes were calculated
        for instr_name in self.misc_options.observation_instruments:

            # Determine the path to the output flux table
            path = fs.join(output_path, instr_name + "_fluxes.dat")

            # Debugging
            log.debug("Loading the mock SED '" + instr_name + "' from '" + path + "' ...")

            # Load the sed
            sed = ObservedSED.from_file(path)

            # Set the sed
            self._mock_seds[instr_name] = sed

    # -----------------------------------------------------------------

    @property
    def mock_seds(self):

        """
        This function ...
        :return:
        """

        if self._mock_seds is not None: return self._mock_seds
        elif self.flux_calculator is not None: return self.flux_calculator.mock_seds
        else: return None #raise ValueError("No mock SEDs")

    # -----------------------------------------------------------------

    @lazyproperty
    def fluxes_from_images_instruments(self):

        """
        This fucntion ...
        :return:
        """

        if self.misc_options.fluxes_from_images_instrument is not None: return [self.misc_options.fluxes_from_images_instrument]
        else: return self.misc_options.observation_instruments

    # -----------------------------------------------------------------

    @lazyproperty
    def single_fluxes_from_images_instrument(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_fluxes_from_images_instrument: raise ValueError("Not a single instrument")
        else: return self.fluxes_from_images_instruments[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def nfluxes_from_images_instruments(self):

        """
        This function ...
        :return:
        """

        return len(self.fluxes_from_images_instruments)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_single_fluxes_from_images_instrument(self):

        """
        This function ...
        :return:
        """

        return self.nfluxes_from_images_instruments == 1

    # -----------------------------------------------------------------

    @lazyproperty
    def fluxes_from_images_coordinate_systems(self):

        """
        This function ...
        :return:
        """

        # No WCS is specified
        if self.misc_options.fluxes_from_images_wcs is None:

            # Masks are given: WCS has to be specified
            if self.has_fluxes_from_images_masks: raise ValueError("When masks are specified for the images, coordinate system of the simulated datacube must also be specified")
            else: return None

        # Wcs is specified, but multiple instruments: fail
        elif not self.has_single_fluxes_from_images_instrument: raise ValueError("Cannot specifify coordinate system when multiple instruments should be used for creating image fluxes")

        # Return the coordinate system
        else:

            # Load the coordinate system
            wcs = CoordinateSystem.from_file(self.misc_options.fluxes_from_images_wcs)

            # Create dictionary for the sole instruments
            coordinate_systems = dict()
            coordinate_systems[self.single_fluxes_from_images_instrument] = wcs

            # Return the dictionary
            return coordinate_systems

    # -----------------------------------------------------------------

    @property
    def filters_for_fluxes_from_images(self):

        """
        This function ...
        :return:
        """

        if self.has_fluxes_from_images_masks: return self.fluxes_from_images_masks[self.single_fluxes_from_images_instrument].keys()
        else: return self.misc_options.observation_filters

    # -----------------------------------------------------------------

    @property
    def has_fluxes_from_images_masks(self):

        """
        This function ...
        :return:
        """

        return self.misc_options.fluxes_from_images_masks is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def fluxes_from_images_masks(self):

        """
        This function ...
        :return:
        """

        # No mask paths specified
        if self.misc_options.fluxes_from_images_masks is None: return None

        # Check whether only one instrument is specified
        if not self.has_single_fluxes_from_images_instrument: raise ValueError("When masks are specified, only one instrument can be used")

        # Initialize a dictionary for the masks per filter
        masks_instrument = dict()

        # Loop over the filter names
        for filter_name in self.misc_options.fluxes_from_images_masks:

            # Debugging
            log.debug("Loading mask for '" + filter_name + "' filter ...")

            # Get filter
            fltr = parse_filter(filter_name)

            # Load the mask from frame
            if self.misc_options.fluxes_from_images_mask_from_nans:

                # Determine path
                frame_path = self.misc_options.fluxes_from_images_masks[filter_name]

                # Load mask
                mask = Mask.nans_from_file(frame_path, plane=self.misc_options.fluxes_from_images_mask_plane)

            # Load mask directly
            else:

                # Determine path
                mask_path = self.misc_options.fluxes_from_images_masks[filter_name]

                # Load mask
                mask = Mask.from_file(mask_path, plane=self.misc_options.fluxes_from_images_mask_plane)

            # Set the mask
            masks_instrument[fltr] = mask

        # Create dictionary with the masks for the sole instrument
        masks = dict()
        masks[self.single_fluxes_from_images_instrument] = masks_instrument

        # Return the dictionary of masks
        return masks

    # -----------------------------------------------------------------

    @lazyproperty
    def fluxes_from_images_output_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.misc_options.path, "image fluxes", recursive=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def fluxes_from_images_reference_seds(self):

        """
        This function ...
        :return:
        """

        from ..data.sed import ObservedSED
        if self.misc_options.plot_fluxes_from_images_reference_seds is None: return None
        else:
            seds = OrderedDict()
            for label in self.misc_options.plot_fluxes_from_images_reference_seds: seds[label] = ObservedSED.from_file(self.misc_options.plot_fluxes_from_images_reference_seds[label])
            return seds

    # -----------------------------------------------------------------

    def get_observed_fluxes_from_images(self):

        """
        This function ...
        :return:
        """

        # Create
        if not self.has_fluxes_from_images: self.calculate_observed_fluxes_from_images()

        # Load
        else: self.load_observed_fluxes_from_images()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_fluxes_from_images_input(self):

        """
        This function ...
        :return:
        """

        # Set input
        input_dict = dict()

        # Set input
        input_dict["simulation"] = self.simulation
        input_dict["output_path"] = self.fluxes_from_images_output_path
        input_dict["filter_names"] = self.filters_for_fluxes_from_images
        input_dict["instrument_names"] = self.fluxes_from_images_instruments
        input_dict["errors"] = self.misc_options.fluxes_from_images_errors
        input_dict["no_spectral_convolution_filters"] = self.misc_options.no_fluxes_from_images_spectral_convolution_filters
        input_dict["coordinate_systems"] = self.fluxes_from_images_coordinate_systems
        input_dict["masks"] = self.fluxes_from_images_masks
        input_dict["reference_seds"] = self.fluxes_from_images_reference_seds

        # Remote
        if not self.config.local:
            input_dict["host_id"] = self.misc_options.fluxes_from_images_remote
            input_dict["remote_images_spectral_convolution"] = self.misc_options.fluxes_from_images_remote_spectral_convolution
            input_dict["remote_threshold"] = self.misc_options.fluxes_from_images_remote_threshold
            input_dict["remote_npixels_threshold"] = self.misc_options.fluxes_from_images_remote_npixels_threshold
            input_dict["remote_rebin_threshold"] = self.misc_options.fluxes_from_images_rebin_remote_threshold

        # Return the input
        return input_dict

    # -----------------------------------------------------------------

    def calculate_observed_fluxes_from_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the observed fluxes from the output images ...")

        # Create
        self.image_flux_calculator = ObservedFluxCalculator()

        # Set spectral convolution flag
        self.image_flux_calculator.config.spectral_convolution = self.misc_options.fluxes_from_images_spectral_convolution

        # Set plot SEDs flag
        self.image_flux_calculator.config.plot = True
        self.image_flux_calculator.config.plot_seds = self.misc_options.plot_fluxes_from_images

        # Set from images flag
        self.image_flux_calculator.config.from_images = True
        self.image_flux_calculator.config.write_images = self.misc_options.write_fluxes_images

        # Plot
        self.image_flux_calculator.config.plot_images = self.misc_options.plot_fluxes_images

        # EXTRA OPTIONS
        self.image_flux_calculator.config.check_wavelengths = self.config.check_wavelengths
        self.image_flux_calculator.config.ignore_bad = self.config.ignore_bad
        self.image_flux_calculator.config.skip_ignored_bad_convolution = self.config.skip_ignored_bad_convolution
        self.image_flux_calculator.config.skip_ignored_bad_closest = self.config.skip_ignored_bad_closest

        # DEPLOYMENT
        self.image_flux_calculator.config.deploy_pts = self.config.deploy_pts
        self.image_flux_calculator.config.update_dependencies = self.config.update_dependencies
        self.image_flux_calculator.config.pubkey_password = self.config.pubkey_password

        # Run
        self.image_flux_calculator.run(**self.observed_fluxes_from_images_input)

        # Done
        self.simulation.analysed_misc.append(fluxes_from_images_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    def load_observed_fluxes_from_images(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading observed fluxes from output images ...")

        from ..data.sed import ObservedSED

        # Get output path
        output_path = self.fluxes_from_images_output_path

        # Initialize dictionary
        self._mock_seds_from_images = OrderedDict()

        # Loop over the instruments for which observed fluxes were calculated
        for instr_name in self.fluxes_from_images_instruments:

            # Determine the path to the output flux table
            path = fs.join(output_path, instr_name + "_fluxes.dat")

            # Debugging
            log.debug("Loading the mock SED (from images) '" + instr_name + "' from '" + path + "' ...")

            # Load the sed
            sed = ObservedSED.from_file(path)

            # Set the sed
            self._mock_seds_from_images[instr_name] = sed

    # -----------------------------------------------------------------

    @property
    def mock_seds_from_images(self):

        """
        This function ...
        :return:
        """

        if self._mock_seds_from_images is not None: return self._mock_seds_from_images
        elif self.image_flux_calculator is not None: return self.image_flux_calculator.mock_seds
        else: return None #raise ValueError("No mock SEDs from images")

    # -----------------------------------------------------------------

    @property
    def filters_for_images(self):

        """
        This function ...
        :return:
        """

        return sequences.elements_not_in_other(self.misc_options.observation_filters, self.misc_options.no_images_filters, check_existing=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def images_output_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.misc_options.path, "images", recursive=True)

    # -----------------------------------------------------------------

    def get_observed_images(self):

        """
        Thisn function ...
        :return:
        """

        # Create
        if not self.has_images: self.make_observed_images()

        # Load
        else: self.load_observed_images()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_images_input(self):

        """
        This function ...
        :return:
        """

        # Set input
        input_dict = dict()

        # General things
        input_dict["simulation"] = self.simulation
        input_dict["output_path"] = self.images_output_path

        # Filters and instruments
        input_dict["filter_names"] = self.filters_for_images
        input_dict["instrument_names"] = self.misc_options.observation_instruments

        # Coordinate system of the datacubes
        if types.is_dictionary(self.misc_options.images_wcs): input_dict["wcs_paths"] = self.misc_options.images_wcs
        elif types.is_string_type(self.misc_options.images_wcs): input_dict["wcs_path"] = self.misc_options.images_wcs
        elif self.misc_options.images_wcs is not None: raise ValueError("Invalid value for 'images_wcs' misc option: " + str(self.misc_options.images_wcs))
        input_dict["wcs_instrument"] = self.misc_options.wcs_instrument

        # Unit conversion
        input_dict["unit"] = self.misc_options.images_unit

        # Convolution
        input_dict["auto_psfs"] = self.misc_options.images_psfs_auto
        input_dict["kernel_paths"] = self.misc_options.images_kernels
        if self.misc_options.fwhms_dataset is not None: input_dict["fwhms_dataset"] = self.misc_options.fwhms_dataset

        # Rebinning
        if types.is_dictionary(self.misc_options.rebin_wcs): input_dict["rebin_wcs_paths"] = self.misc_options.rebin_wcs
        elif types.is_string_type(self.misc_options.rebin_wcs): input_dict["rebin_wcs_path"] = self.misc_options.rebin_wcs
        elif self.misc_options.rebin_wcs is not None: raise ValueError("Invalid value for 'rebin_wcs' misc option: " + str(self.misc_options.rebin_wcs))
        input_dict["rebin_dataset"] = self.misc_options.rebin_dataset  # path
        input_dict["rebin_instrument"] = self.misc_options.rebin_instrument

        # Remote
        if not self.config.local:
            input_dict["host_id"] = self.misc_options.make_images_remote
            input_dict["remote_spectral_convolution"] = self.misc_options.remote_spectral_convolution
            input_dict["remote_threshold"] = self.misc_options.images_remote_threshold
            input_dict["remote_npixels_threshold"] = self.misc_options.images_remote_npixels_threshold
            input_dict["remote_rebin_threshold"] = self.misc_options.rebin_remote_threshold
            input_dict["remote_convolve_threshold"] = self.misc_options.convolve_remote_threshold

        # NO SPECTRAL CONVOLUTION FOR CERTAIN IMAGES?
        input_dict["no_spectral_convolution_filters"] = self.misc_options.no_images_spectral_convolution_filters

        # Return the input
        return input_dict

    # -----------------------------------------------------------------

    def make_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the observed images (this may take a while) ...")

        # Create and run an ObservedImageMaker object
        self.image_maker = ObservedImageMaker()

        # Set spectral convolution flag
        self.image_maker.config.spectral_convolution = self.misc_options.images_spectral_convolution

        # Set group flag
        self.image_maker.config.group = self.misc_options.group_images

        # Set number of processes
        self.image_maker.config.nprocesses_local = self.misc_options.images_nprocesses_local
        self.image_maker.config.nprocesses_remote = self.misc_options.images_nprocesses_remote

        # Set other
        self.image_maker.config.write_intermediate = self.misc_options.write_intermediate_images
        self.image_maker.config.write_kernels = self.misc_options.write_convolution_kernels

        # Plot
        self.image_maker.config.plot_images = self.misc_options.plot_images

        # Run
        self.image_maker.run(**self.observed_images_input)

        # Done
        self.simulation.analysed_misc.append(images_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        pass

        # Debugging
        #log.debug("Loading observed images ...")

        #raise NotImplementedError("Not yet implemented")

# -----------------------------------------------------------------
