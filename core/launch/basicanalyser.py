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
from ..plot.seds import plotseds
from ..plot.grids import plotgrids
from ..plot.rgbimages import makergbimages
from ..plot.wavemovie import makewavemovie
from ..misc.fluxes import ObservedFluxCalculator
from ..misc.images import ObservedImageMaker
from ..basics.log import log
from ..tools import filesystem as fs
from ..plot.sed import SEDPlotter
from ..data.sed import SED, ObservedSED
from ..tools import types
from .options import progress_name, timeline_name, memory_name, seds_name, grids_name, rgb_name, wave_name, fluxes_name, images_name

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

        # The flux calculator and image maker
        self.flux_calculator = None
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
        if self.extraction_options.progress: self.extract_progress()

        # Extract the timeline information
        if self.extraction_options.timeline: self.extract_timeline()

        # Extract the memory information
        if self.extraction_options.memory: self.extract_memory()

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
    def has_wave(self):

        """
        Thins function ...
        :return:
        """

        return wave_name in self.simulation.analysed_misc

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

        # If requested, make wave movies from the ouput FITS files
        if self.misc_options.wave and not self.has_wave: self.make_wave()

        # If requested, calculate observed fluxes from the output SEDs
        if self.misc_options.fluxes and not self.has_fluxes: self.calculate_observed_fluxes()

        # If requested, create observed imgaes from the output FITS files
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
        path = fs.join(self.extraction_options.path, "progress.dat")

        # Run the progress extractor
        self.progress = extractor.run(self.simulation, path)

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
        path = fs.join(self.extraction_options.path, "timeline.dat")

        # Run the timeline extractor
        self.timeline = extractor.run(self.simulation, path)

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
        path = fs.join(self.extraction_options.path, "memory.dat")

        # Run the memory extractor
        self.memory = extractor.run(self.simulation, path)

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

        # If the simulated SED must be plotted against a set of reference flux points
        if self.plotting_options.reference_seds is not None:

            # Inform the user
            log.info("Plotting the SED with reference fluxes ...")

            # Create a new SEDPlotter instance
            plotter = SEDPlotter()

            # Get the simulation prefix
            prefix = self.simulation.prefix()

            # Loop over the simulated SED files and add the SEDs to the SEDPlotter
            for sed_path in self.simulation.seddatpaths():

                # Determine the name of the corresponding instrument
                instr_name = instrument_name(sed_path, prefix)

                # Load the SED
                sed = SED.from_skirt(sed_path)

                # Add the simulated SED to the plotter
                plotter.add_sed(sed, instr_name)

            # Add the reference SEDs
            for reference_sed_path in self.plotting_options.reference_seds:

                # Determine name
                reference_sed_name = fs.strip_extension(fs.name(reference_sed_path))

                # Add the reference SED
                reference_sed = ObservedSED.from_file(reference_sed_path)
                plotter.add_sed(reference_sed, reference_sed_name)

            # Determine the path to the plot file
            path = fs.join(self.plotting_options.path, "sed." + self.plotting_options.format)
            plotter.run(title=self.simulation.name, output=path)

            # Get the axis limits
            min_wavelength = plotter.min_wavelength
            max_wavelength = plotter.max_wavelength
            min_flux = plotter.min_flux
            max_flux = plotter.max_flux

            # Clear the SED plotter
            plotter.clear()

            # Check which SED files are produced by a FullInstrument (these files also contain the full SED of the various contributions)
            for sed_path in self.simulation.seddatpaths():

                # Determine the name of the corresponding instrument
                instr_name = instrument_name(sed_path, prefix)

                # Check how many columns the SED file contains
                ncols = number_of_columns(sed_path)

                # Check the type of the Instrument / SED
                if ncols == 2: continue # SEDInstrument

                for contribution in ["total", "direct", "scattered", "dust", "dustscattered", "transparent"]:

                    # Load the SED contribution
                    sed = SED.from_skirt(sed_path, contribution=contribution)

                    # Add the SED to the plotter
                    plotter.add_sed(sed, contribution, residuals=(contribution == "total"))

                # Add the reference SEDs
                for reference_sed_path in self.plotting_options.reference_seds:

                    # Determine name
                    reference_sed_name = fs.strip_extension(fs.name(reference_sed_path))

                    # Add the reference SED
                    reference_sed = ObservedSED.from_file(reference_sed_path)
                    plotter.add_sed(reference_sed, reference_sed_name)

                # Determine the path to the plot file
                path = fs.join(self.plotting_options.path, "sed_" + instr_name + "." + self.plotting_options.format)

                # Plot
                plotter.run(output=path, min_wavelength=min_wavelength, max_wavelength=max_wavelength, min_flux=min_flux, max_flux=max_flux)

                # Clear the SED plotter
                plotter.clear()

        # Use the simple plotseds function
        else: plotseds(self.simulation, output_path=self.plotting_options.path, format=self.plotting_options.format)

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

    def make_rgb(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making RGB images ...")

        # Make RGB images from the output images
        makergbimages(self.simulation, output_path=self.misc_options.path)

        # Done
        self.simulation.analysed_misc.append(rgb_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    def make_wave(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making wave movies ...")

        # Make wave movies from the output images
        #makewavemovie(self.simulation, output_path=self.misc_options.path)

        # Done
        self.simulation.analysed_misc.append(wave_name)
        self.simulation.save()

    # -----------------------------------------------------------------

    def calculate_observed_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the observed fluxes ...")

        # Create and run a ObservedFluxCalculator object
        self.flux_calculator = ObservedFluxCalculator()

        # Set spectral convolution flag
        self.flux_calculator.config.spectral_convolution = self.misc_options.spectral_convolution

        # Run
        self.flux_calculator.run(simulation=self.simulation, output_path=self.misc_options.path,
                                 filter_names=self.misc_options.observation_filters,
                                 instrument_names=self.misc_options.observation_instruments,
                                 errors=self.misc_options.flux_errors)

        # Done
        self.simulation.analysed_misc.append(fluxes_name)
        self.simulation.save()

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
        self.image_maker.config.spectral_convolution = self.misc_options.spectral_convolution

        # Set group flag
        self.image_maker.config.group = self.misc_options.group_images

        # Set input
        input_dict = dict()

        # General things
        input_dict["simulation"] = self.simulation
        input_dict["output_path"] = self.misc_options.path

        # Filters and instruments
        input_dict["filter_names"] = self.misc_options.observation_filters
        input_dict["instrument_names"] = self.misc_options.observation_instruments

        # Coordinate system of the datacubes
        #input_dict["wcs"] =
        if types.is_dictionary(self.misc_options.images_wcs): input_dict["wcs_paths"] = self.misc_options.images_wcs
        elif types.is_string_type(self.misc_options.images_wcs): input_dict["wcs_path"] = self.misc_options.images_wcs
        else: raise ValueError("Invalid value for 'images_wcs' misc option: " + str(self.misc_options.images_wcs))
        input_dict["wcs_instrument"] = self.misc_options.wcs_instrument

        # Unit conversion
        input_dict["unit"] = self.misc_options.images_unit

        # Convolution
        input_dict["auto_psfs"] = self.misc_options.images_psfs_auto
        input_dict["kernel_paths"] = self.misc_options.images_kernels
        #input_dict["psf_paths"] =

        # Rebinning
        if types.is_dictionary(self.misc_options.rebin_wcs): input_dict["rebin_wcs_paths"] = self.misc_options.rebin_wcs
        elif types.is_string_type(self.misc_options.rebin_wcs): input_dict["rebin_wcs_path"] = self.misc_options.rebin_wcs
        else: raise ValueError("Invalid value for 'rebin_wcs' misc option: " + str(self.misc_options.rebin_wcs))
        #input_dict["rebin_wcs"] =
        input_dict["rebin_dataset"] = self.misc_options.rebin_dataset # path
        input_dict["rebin_instrument"] = self.misc_options.rebin_instrument

        # Remote
        input_dict["host_id"] = self.misc_options.make_images_remote
        input_dict["remote_threshold"] = self.misc_options.images_remote_threshold
        input_dict["remote_rebin_threshold"] = self.misc_options.rebin_remote_threshold
        input_dict["remote_convolve_threshold"] = self.misc_options.convolve_remote_threshold

        # Run
        self.image_maker.run(**input_dict)

        # Done
        self.simulation.analysed_misc.append(images_name)
        self.simulation.save()

# -----------------------------------------------------------------

def instrument_name(sed_path, prefix):

    """
    This function ...
    :param sed_path:
    :param prefix:
    :return:
    """

    return fs.name(sed_path).split("_sed.dat")[0].split(prefix + "_")[1]

# -----------------------------------------------------------------

def number_of_columns(sed_path):

    """
    This function ...
    :param sed_path:
    :return:
    """

    with open(sed_path, 'r') as f:

        ncols = 0
        for line in f:

            if "# column" not in line: break
            else: ncols = int(line.split("column ")[1].split(": ")[0])

    return ncols

# -----------------------------------------------------------------
