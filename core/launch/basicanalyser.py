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

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..extract.progress import ProgressExtractor
from ..extract.timeline import TimeLineExtractor
from ..extract.memory import MemoryExtractor
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

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Extract information from the simulation's log files
        self.extract()

        # 3. Make plots based on the simulation output
        self.plot()

        # 4. Miscellaneous output
        self.misc()

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
        self.simulation = kwargs.pop("simulation")

        # Update the analysis options for added features
        self.simulation.update_analysis_options()

        # Also make references to the simulation's analysis options for extraction, plotting and misc (for shorter notation)
        self.extraction_options = self.simulation.analysis.extraction
        self.plotting_options = self.simulation.analysis.plotting
        self.misc_options = self.simulation.analysis.misc

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

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting ...")

        # Debugging
        log.debug("Extraction options:")
        if log.is_debug(): print(str(self.extraction_options))

        # Extract the progress information
        if self.extraction_options.progress and not self.extracted_progress: self.extract_progress()
        elif self.extraction_options.progress: self.load_progress()

        # Extract the timeline information
        if self.extraction_options.timeline and not self.extracted_timeline: self.extract_timeline()
        elif self.extraction_options.timeline: self.load_timeline()

        # Extract the memory information
        if self.extraction_options.memory and not self.extracted_memory: self.extract_memory()
        elif self.extraction_options.memory: self.load_memory()

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

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Debugging
        log.debug("Plotting options:")
        if log.is_debug(): print(str(self.plotting_options))

        # If requested, plot the SED's
        if self.plotting_options.seds and not self.plotted_seds: self.plot_seds()

        # If requested, make plots of the dust grid
        if self.plotting_options.grids and not self.plotted_grids: self.plot_grids()

        # If requested, plot the simulation progress as a function of time
        if self.plotting_options.progress and not self.plotted_progress: self.plot_progress()

        # If requested, plot a timeline of the different simulation phases
        if self.plotting_options.timeline and not self.plotted_timeline: self.plot_timeline()

        # If requested, plot the memory usage as a function of time
        if self.plotting_options.memory and not self.plotted_memory: self.plot_memory()

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

    def misc(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing miscellaneous analysis ...")

        # Debugging
        log.debug("Miscellaneous options:")
        if log.is_debug(): print(str(self.misc_options))

        # If requested, make RGB images of the output FITS files
        if self.misc_options.rgb and not self.has_rgb: self.make_rgb()

        # If requested, make datacube animations from the output datacubes
        if self.misc_options.animations and not self.has_animations: self.make_animations()

        # If requested, calculate observed fluxes from the output SEDs
        if self.misc_options.fluxes and not self.has_fluxes: self.calculate_observed_fluxes()

        # If requested, calculate observed fluxes from the output datacubes
        if self.misc_options.fluxes_from_images and not self.has_fluxes_from_images: self.calculate_observed_fluxes_from_images()

        # If requested, create observed imgaes from the output datacubes
        if self.misc_options.images and not self.has_images: self.make_observed_images()

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
        extractor.run(self.simulation, path)
        self.progress = extractor.table

        # Done
        self.simulation.analysed_extraction.append(progress_name)
        self.simulation.save()

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
        extractor.run(self.simulation, path)
        self.timeline = extractor.table

        # Done
        self.simulation.analysed_extraction.append(timeline_name)
        self.simulation.save()

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
        extractor.run(self.simulation, path)
        self.memory = extractor.table

        # Done
        self.simulation.analysed_extraction.append(memory_name)
        self.simulation.save()

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
        plotter.run(self.simulation)

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
        plotgrids(self.simulation, output_path=self.plotting_options.path, silent=(not log.is_debug()))

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

    @lazyproperty
    def fluxes_reference_sed(self):

        """
        This function ...
        :return:
        """

        from ..data.sed import ObservedSED
        if self.misc_options.plot_fluxes_reference_sed is None: return None
        else: return ObservedSED.from_file(self.misc_options.plot_fluxes_reference_sed)

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
        self.flux_calculator.config.plot = self.misc_options.plot_fluxes

        # Run
        self.flux_calculator.run(simulation=self.simulation, output_path=self.fluxes_output_path,
                                 filter_names=self.filters_for_fluxes,
                                 instrument_names=self.misc_options.observation_instruments,
                                 errors=self.misc_options.flux_errors,
                                 no_spectral_convolution_filters=self.misc_options.no_fluxes_spectral_convolution_filters,
                                 reference_sed=self.fluxes_reference_sed)

        # Done
        self.simulation.analysed_misc.append(fluxes_name)
        self.simulation.save()

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
    def fluxes_from_images_reference_sed(self):

        """
        This function ...
        :return:
        """

        from ..data.sed import ObservedSED
        if self.misc_options.plot_fluxes_from_images_reference_sed is None: return None
        else: return ObservedSED.from_file(self.misc_options.plot_fluxes_from_images_reference_sed)

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

        # Set plot flag
        self.image_flux_calculator.config.plot = self.misc_options.plot_fluxes_from_images

        # Set from images flag
        self.image_flux_calculator.config.from_images = True
        self.image_flux_calculator.config.write_images = self.misc_options.write_fluxes_images

        #print("instruments:", self.fluxes_from_images_instruments)
        #print("masks:" , self.fluxes_from_images_masks)
        #print("wcs:", self.fluxes_from_images_coordinate_systems)

        # Run
        self.image_flux_calculator.run(simulation=self.simulation, output_path=self.fluxes_from_images_output_path,
                                       filter_names=self.filters_for_fluxes_from_images,
                                       instrument_names=self.fluxes_from_images_instruments,
                                       errors=self.misc_options.fluxes_from_images_errors,
                                       no_spectral_convolution_filters=self.misc_options.no_fluxes_from_images_spectral_convolution_filters,
                                       coordinate_systems=self.fluxes_from_images_coordinate_systems,
                                       masks=self.fluxes_from_images_masks, reference_sed=self.fluxes_from_images_reference_sed)

        # Done
        self.simulation.analysed_misc.append(fluxes_from_images_name)
        self.simulation.save()

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

    def make_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the observed images (this may take a while) ...")

        # Create and run an ObservedImageMaker object
        self.image_maker = ObservedImageMaker()

        # Set spectral convolution flags
        #self.image_maker.config.spectral_convolution = self.misc_options.spectral_convolution
        self.image_maker.config.spectral_convolution = self.misc_options.images_spectral_convolution

        # Set group flag
        self.image_maker.config.group = self.misc_options.group_images

        # Set number of processes
        self.image_maker.config.nprocesses_local = self.misc_options.images_nprocesses_local
        self.image_maker.config.nprocesses_remote = self.misc_options.images_nprocesses_remote

        # Set other
        self.image_maker.config.write_intermediate = self.misc_options.write_intermediate_images
        self.image_maker.config.write_kernels = self.misc_options.write_convolution_kernels

        # Set input
        input_dict = dict()

        # General things
        input_dict["simulation"] = self.simulation
        input_dict["output_path"] = self.images_output_path

        # Filters and instruments
        input_dict["filter_names"] = self.filters_for_images
        input_dict["instrument_names"] = self.misc_options.observation_instruments

        # Coordinate system of the datacubes
        #input_dict["wcs"] =
        if types.is_dictionary(self.misc_options.images_wcs): input_dict["wcs_paths"] = self.misc_options.images_wcs
        elif types.is_string_type(self.misc_options.images_wcs): input_dict["wcs_path"] = self.misc_options.images_wcs
        elif self.misc_options.images_wcs is not None: raise ValueError("Invalid value for 'images_wcs' misc option: " + str(self.misc_options.images_wcs))
        input_dict["wcs_instrument"] = self.misc_options.wcs_instrument

        # Unit conversion
        input_dict["unit"] = self.misc_options.images_unit

        # Convolution
        input_dict["auto_psfs"] = self.misc_options.images_psfs_auto
        input_dict["kernel_paths"] = self.misc_options.images_kernels
        #input_dict["psf_paths"] =
        if self.misc_options.fwhms_dataset is not None: input_dict["fwhms_dataset"] = self.misc_options.fwhms_dataset

        # Rebinning
        if types.is_dictionary(self.misc_options.rebin_wcs): input_dict["rebin_wcs_paths"] = self.misc_options.rebin_wcs
        elif types.is_string_type(self.misc_options.rebin_wcs): input_dict["rebin_wcs_path"] = self.misc_options.rebin_wcs
        elif self.misc_options.rebin_wcs is not None: raise ValueError("Invalid value for 'rebin_wcs' misc option: " + str(self.misc_options.rebin_wcs))

        #input_dict["rebin_wcs"] =
        input_dict["rebin_dataset"] = self.misc_options.rebin_dataset # path
        input_dict["rebin_instrument"] = self.misc_options.rebin_instrument

        # Remote
        input_dict["host_id"] = self.misc_options.make_images_remote
        input_dict["remote_spectral_convolution"] = self.misc_options.remote_spectral_convolution
        input_dict["remote_threshold"] = self.misc_options.images_remote_threshold
        input_dict["remote_rebin_threshold"] = self.misc_options.rebin_remote_threshold
        input_dict["remote_convolve_threshold"] = self.misc_options.convolve_remote_threshold

        # NO SPECTRAL CONVOLUTION FOR CERTAIN IMAGES?
        input_dict["no_spectral_convolution_filters"] = self.misc_options.no_images_spectral_convolution_filters

        # Run
        self.image_maker.run(**input_dict)

        # Done
        self.simulation.analysed_misc.append(images_name)
        self.simulation.save()

# -----------------------------------------------------------------
