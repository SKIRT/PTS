#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.fitting Contains the FittingPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import defaultdict, OrderedDict

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ..fitting.component import FittingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.basics.distribution import Distribution
from ...core.launch.timing import TimingTable
from ...core.data.transmission import TransmissionCurve
from ...core.plot.transmission import TransmissionPlotter
from ...core.plot.wavelengthgrid import WavelengthGridPlotter
from ...core.plot.grids import plotgrids
from ...core.simulation.simulation import SkirtSimulation
from ...core.simulation.logfile import LogFile
from ..core.emissionlines import EmissionLines
from ...core.data.sed import load_example_mappings_sed, load_example_bruzualcharlot_sed, load_example_zubko_sed
from ...core.data.sed import SED, ObservedSED
from ...magic.plot.imagegrid import ResidualImageGridPlotter
from ...magic.core.frame import Frame
from ...core.plot.sed import SEDPlotter
from ...magic.region.list import SkyRegionList
from ..misc.geometryplotter import GeometryPlotter
from ..basics.models import load_3d_model
from ...core.plot.distribution import DistributionPlotter

# -----------------------------------------------------------------

class FittingPlotter(PlottingComponent, FittingComponent):
    
    """
    This class...
    """

    # The load functions
    load_functions = dict()

    # The plot functions
    plot_functions = dict()

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base classes
        PlottingComponent.__init__(self, config)
        FittingComponent.__init__(self, config)

        # -- Attributes --

        # Paths

        self.plot_fitting_wavelength_grids_path = None
        self.plot_fitting_dust_grids_path = None
        self.plot_fitting_distributions_path = None
        self.plot_fitting_prior_distributions_path = None
        self.plot_fitting_posterior_distributions_path = None

        # The wavelength grids
        self.wavelength_grids = []

        # The transmission curves
        self.transmission_curves = dict()

        # The distribution of dust cells per level
        self.dust_cell_trees = []

        # The runtimes
        self.runtimes = None

        # The simulated SEDs of all models, as lists for each generation
        self.seds = dict()

        # The SEDs of the different stellar contributions (total, old, young, ionizing)
        self.sed_contributions = OrderedDict()

        # The observed SED
        self.observed_sed = None

        # The simulated images
        self.simulated_images = dict()

        # The observed imags
        self.observed_images = dict()

        # The geometries
        self.geometries = dict()

        # The prior and posterior probability distributions
        self.prior_distributions = dict()     # indexed on generation, then on the free parameter labels
        self.posterior_distributions = dict() # indexed on the free parameter labels

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingPlotter, self).setup(**kwargs)

        # Create directories
        self.create_directories()

    # -----------------------------------------------------------------

    def create_directories(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the necessary directories ...")

        # Directory for plotting the wavelength grids
        self.plot_fitting_wavelength_grids_path = fs.create_directory_in(self.plot_fitting_path, "wavelength grids")

        # Directory for plotting the dust grids
        self.plot_fitting_dust_grids_path = fs.create_directory_in(self.plot_fitting_path, "dust grids")

        # Loop over the finished generations
        for generation_name in self.finished_generations:

            # Determine the path
            path = fs.join(self.plot_fitting_path, generation_name)
            if fs.is_directory(path): continue # plotting has already been done for this finished generation

            # Create the directory and set the path
            #fs.create_directory(path)
            #self.plot_fitting_generation_paths[generation_name] = path

        # Directory for plotting probability distributions
        self.plot_fitting_distributions_path = fs.create_directory_in(self.plot_fitting_path, "distributions")

        # Directories for prior and posterior distributions
        self.plot_fitting_prior_distributions_path = fs.create_directory_in(self.plot_fitting_distributions_path, "prior")
        self.plot_fitting_posterior_distributions_path = fs.create_directory_in(self.plot_fitting_distributions_path, "posterior")

    # -----------------------------------------------------------------

    def load_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the wavelength grids used for the fitting ...")

        # Loop over the files found in the fit/wavelength grids directory
        for path in fs.files_in_path(self.fit_wavelength_grids_path, extension="txt", sort=int, not_contains="grids"):

            # Load the wavelength grid
            grid = WavelengthGrid.from_skirt_input(path)

            # Add the grid
            self.wavelength_grids.append(grid)

    # -----------------------------------------------------------------

    def load_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the distributions ...")

        # Load prior distributions
        self.load_prior_distributions()

        # Load the posterior distributions
        self.load_posterior_distributions()

    # -----------------------------------------------------------------

    def load_prior_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the prior probability distributions for the different free parameters for each generation ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Initialize a dictionary
            distributions_dict = dict()

            # Load the parameters table
            parameters_table = self.parameters_table_for_generation(generation_name)

            # Loop over the different parameters
            for label in self.free_parameter_labels:

                # Test if the values are not all equal
                if np.min(parameters_table[label]) == np.max(parameters_table[label]):

                    # Give a warning and skip this parameter (we can't make a distribution of values that are all the same)
                    log.warning("The prior values for the '" + label + "' parameter of the '" + generation_name + "' generation are all identical")
                    continue

                # Create a distribution from the list of parameter values
                distribution = Distribution.from_values(parameters_table[label])

                # Add the distribution to the dictionary
                distributions_dict[label] = distribution

            # Add the dictionary for this generation
            self.prior_distributions[generation_name] = distributions_dict

    # -----------------------------------------------------------------

    def load_posterior_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the probability distributions for the different fit parameters ...")

        # Loop over the different fit parameters
        for parameter_name in self.free_parameter_labels:

            # Check if the posterior distribution has been calculated
            if not self.has_distribution(parameter_name): continue

            # Load and set the distribution
            self.posterior_distributions[parameter_name] = self.get_parameter_distribution(parameter_name)

    # -----------------------------------------------------------------

    def load_transmission_curves(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the transmission curves for the observation filters ...")

        # Load the observed SED
        sed = ObservedSED.from_file(self.observed_sed_path)
        for fltr in sed.filters():

            # Create the transmission curve
            transmission = TransmissionCurve.from_filter(fltr)

            # Normalize the transmission curve
            transmission.normalize(value=1.0, method="max")

            # Add the transmission curve to the dictionary
            self.transmission_curves[str(fltr)] = transmission

    # -----------------------------------------------------------------

    def load_dust_cell_trees(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the dust cell tree information ...")

        # Loop over the directories inside the fit/dust grids directory
        for path in fs.directories_in_path(self.fit_dust_grids_path, sort=int):

            # Determine the path to the log file of the dust grid generating simulation
            log_file_path = fs.join(path, self.galaxy_name + "_log.txt")

            # If the log file does not exist, skip
            if not fs.is_file(log_file_path): continue

            # Open the log file
            log_file = LogFile(log_file_path)

            # Get the distribution of cells per level of the tree
            self.dust_cell_trees.append(log_file.tree_leaf_distribution)

    # -----------------------------------------------------------------

    def load_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the runtimes ...")

        # Determine the path to the timing table
        timing_table_path = fs.join(self.fit_path, "timing.dat")

        try:
            # Load the timing table
            timing_table = TimingTable.from_file(timing_table_path)
        except ValueError: return # ValueError: Column Timestamp failed to convert

        # Keep a list of all the runtimes recorded for a certain remote host
        self.runtimes = defaultdict(lambda: defaultdict(list))

        # Loop over the entries in the timing table
        for i in range(len(timing_table)):

            # Get the ID of the host and the cluster name for this particular simulation
            host_id = timing_table["Host id"][i]
            cluster_name = timing_table["Cluster name"][i]

            # Get the parallelization properties for this particular simulation
            cores = timing_table["Cores"][i]
            threads_per_core = timing_table["Threads per core"][i]
            processes = timing_table["Processes"][i]

            # Get the number of photon packages (per wavelength) used for this simulation
            packages = timing_table["Packages"][i]

            # Get the total runtime
            runtime = timing_table["Total runtime"][i]

            parallelization = (cores, threads_per_core, processes)

            # Add the runtime to the list of runtimes
            self.runtimes[host_id][packages, parallelization].append(runtime)

    # -----------------------------------------------------------------

    def load_seds(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SEDs ...")

        # Load the observed SEDs
        self.load_observed_sed()

        # Load the simulated SEDs
        self.load_model_seds()

        # Load the contributions to the SEDs
        self.load_sed_contributions()

    # -----------------------------------------------------------------

    def load_observed_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed SED ...")

        # Load the observed SED
        self.observed_sed = ObservedSED.from_file(self.observed_sed_path)

    # -----------------------------------------------------------------

    def load_model_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SEDs of all the fit models ...")

        # Loop over all finished generations
        for generation_name in self.finished_generations:

            # Initialize dictionary to contain SEDs for this generation
            #seds_generation = dict()
            seds_generation = []

            # Loop over all simulations in this generation
            for simulation_name in self.get_simulations_in_generation(generation_name):

                # Determine the path to the 'plot' directory for this simulation
                #plot_path = fs.join(self.fit_generations_path, generation_name, simulation_name, "plot")

                out_path = fs.join(self.fit_generations_path, generation_name, simulation_name, "out")

                # Determine the path to the SED plot
                #sed_path = fs.join(plot_path, self.galaxy_name + "_earth_sed.dat")

                # Determine the path to the SED data file
                sed_path = fs.join(out_path, self.galaxy_name + "_earth_sed.dat")

                # Check whether the SED file is present
                if not fs.is_file(sed_path):

                    # Give warning and continue
                    log.warning("The SED file for simulation " + simulation_name + " of generation " + generation_name + " is missing")
                    continue

                # Load the SED
                sed = SED.from_skirt(sed_path)

                # Add the SED
                seds_generation.append(sed)

            # Add the sed to the list of SEDs
            self.seds[generation_name] = seds_generation

    # -----------------------------------------------------------------

    def load_sed_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SEDs of the various stellar contributions for the best models of each generation ...")

        # Loop over the directories in the fit_best directory
        for path, generation_name in fs.directories_in_path(self.fit_best_path, returns=["path", "name"]):

            # Initialize ...
            seds_generation = dict()

            # Loop over the contributions that have been simulated
            for contribution in fs.directories_in_path(path, returns="name"):

                # Determine the path to the output directory
                out_path = fs.join(path, contribution, "out")

                # Determine the path to the SED file
                sed_path = fs.join(out_path, self.galaxy_name + "_earth_sed.dat")

                # Load the SED
                sed = SED.from_skirt(sed_path)

                # Add the SED
                seds_generation[contribution] = sed

            # Add the SEDS
            self.sed_contributions[generation_name] = seds_generation

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # Load simulated images
        self.load_simulated_images()

        # Load observed images
        self.load_observed_images()

    # -----------------------------------------------------------------

    def load_simulated_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated images ...")

        # Determine the path to the fit/best directory
        fit_best_path = fs.join(self.fit_path, "best")

        # Determine the path to the fit/best/images directory
        out_path = fs.join(fit_best_path, "images")

        # Loop over all FITS files found in the fit/best/images directory
        for path, name in fs.files_in_path(out_path, extension="fits", returns=["path", "name"], contains="__"):

            # Debugging
            log.debug("Loading the '" + name + "' image ...")

            # Get the filter name
            #filter_name = name.split("__")[1]

            # Open the image
            frame = Frame.from_file(path)

            # Get the filter name
            filter_name = str(frame.filter)

            # Add the image frame to the dictionary
            self.simulated_images[filter_name] = frame

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed images ...")

        # Loop over all FITS files found in the 'truncated' directory
        for path, name in fs.files_in_path(self.truncation_path, extension="fits", returns=["path", "name"]):

            # Ignore the bulge, disk and model images
            if name == "bulge" or name == "disk" or name == "model": continue

            # Ignore the H alpha image
            if "Halpha" in name: continue

            # Check whether a simulated image exists for this band
            if name not in self.simulated_images:
                log.warning("The simulated version of the " + name + " image could not be found, skipping " + name + " data ...")
                continue

            # Debugging
            log.debug("Loading the '" + name + "' image ...")

            # The filter name is the image name
            filter_name = name

            # Open the image
            frame = Frame.from_file(path)

            # Add the image frame to the dictionary
            self.observed_images[filter_name] = frame

    # -----------------------------------------------------------------

    def load_geometries(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the geometries ...")

        # Load all geometries in the directory
        for path, name in fs.files_in_path(self.fit_geometries_path, extension="mod", returns=["path", "name"]):

            # Load the geometry
            model = load_3d_model(path)

            # Add the geometry
            self.geometries[name] = model

    # -----------------------------------------------------------------

    def plot_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the progress of the best chi squared value over the course of generations ...")



    # -----------------------------------------------------------------

    def plot_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the probability distributions ...")

        # Plot prior distributions
        self.plot_prior_distributions()

        # Plot posterior distributions
        self.plot_posterior_distributions()

    # -----------------------------------------------------------------

    def plot_prior_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the posterior parameter distributions ...")

        # Create a distribution plotter
        plotter = DistributionPlotter()

        # Loop over the generations
        for generation_name in self.prior_distributions:

            # Create directory for this generation
            generation_path = fs.create_directory_in(self.plot_fitting_prior_distributions_path, generation_name)

            # Get the parameter ranges for this generation
            ranges = self.parameter_ranges_for_generation(generation_name)

            # Directories
            linear_path = fs.create_directory_in(generation_path, "linear")
            linear_smooth_path = fs.create_directory_in(generation_path, "linear smooth")
            log_path = fs.create_directory_in(generation_path, "log")
            log_smooth_path = fs.create_directory_in(generation_path, "log smooth")

            # Loop over the distributions (for the different parameters)
            for label in self.prior_distributions[generation_name]:

                # Get the limits
                limits = [ranges[label].min, ranges[label].max]

                # Show the limits
                #print(self.prior_distributions[generation_name][label].min_value)
                #print(self.prior_distributions[generation_name][label].max_value)

                # Define the title
                title = "Prior probability distribution of " + self.parameter_descriptions[label]

                # Add the distribution
                plotter.add_distribution(self.prior_distributions[generation_name][label], label)
                plotter.title = title

                # Plot linear
                path = fs.join(linear_path, label + ".pdf")
                plotter.run(output_path=path)

                plotter.clear(clear_distributions=False)

                # Plot linear smooth
                path = fs.join(linear_smooth_path, label + ".pdf")
                plotter.run(output_path=path, add_smooth=True)

                plotter.clear(clear_distributions=False)

                # Plot log
                path = fs.join(log_path, label + ".pdf")
                plotter.run(output_path=path, logscale=True)

                plotter.clear(clear_distributions=False)

                # Plot log smooth
                path = fs.join(log_smooth_path, label + ".pdf")
                plotter.run(output_path=path, add_smooth=True, logscale=True)

                # Clear
                plotter.clear()

    # -----------------------------------------------------------------

    def plot_posterior_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the posterior parameter distributions ...")

        # Loop over the parameter labels
        for label in self.posterior_distributions:

            # Path
            linear_path = fs.join(self.plot_fitting_posterior_distributions_path, label + "_linear.pdf")
            log_path = fs.join(self.plot_fitting_posterior_distributions_path, label + "_log.pdf")
            cumulative_path = fs.join(self.plot_fitting_posterior_distributions_path, label + "_cumulative.pdf")

            # Get the limits
            x_limits = [self.free_parameter_ranges[label].min, self.free_parameter_ranges[label].max]

            # Plot the distributions
            self.posterior_distributions[label].plot_smooth(x_limits=x_limits, title="Posterior probability distribution of " + self.parameter_descriptions[label], path=linear_path)
            self.posterior_distributions[label].plot_smooth(x_limits=x_limits, title="Posterior probability distribution of " + self.parameter_descriptions[label] + " (in log scale)", path=log_path)
            self.posterior_distributions[label].plot_cumulative_smooth(x_limits=x_limits, title="Cumulative distribution of " + self.parameter_descriptions[label], path=cumulative_path)

    # -----------------------------------------------------------------

    def plot_wavelengths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelength grids ...")

        # Loop over the different wavelength grids
        index = 0
        for grid in self.wavelength_grids:

            # Plot the low-resolution and high-resolution wavelength grids with the observation filters
            path = fs.join(self.plot_fitting_wavelength_grids_path, "wavelengths_" + str(index) + "_filters.pdf")
            self.plot_wavelengths_filters(grid, path)

            # Plot the wavelengths with the SEDs of stars and dust
            path = fs.join(self.plot_fitting_wavelength_grids_path, "wavelengths_" + str(index) + "_seds.pdf")
            self.plot_wavelengths_seds(grid, path)

            # Increment the index
            index += 1

    # -----------------------------------------------------------------

    def plot_wavelengths_filters(self, grid, filename):

        """
        This function ...
        :param grid:
        :param filename:
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelengths with")

        # Create the transmission plotter
        plotter = TransmissionPlotter()

        plotter.title = "Wavelengths used for fitting"
        plotter.transparent = True

        # Add the transmission curves
        for label in self.transmission_curves: plotter.add_transmission_curve(self.transmission_curves[label], label)

        # Add the wavelength points
        for wavelength in grid.wavelengths(): plotter.add_wavelength(wavelength)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, filename)

        # Run the plotter
        plotter.run(path, min_wavelength=grid.min_wavelength, max_wavelength=grid.max_wavelength, min_transmission=0.0, max_transmission=1.05)

    # -----------------------------------------------------------------

    def plot_wavelengths_seds(self, grid, filename):

        """
        This function ...
        :param grid:
        :param filename:
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelengths with SEDs of stars and dust ...")

        # Create the wavelength grid plotter
        plotter = WavelengthGridPlotter()

        # Set title
        plotter.title = "Wavelengths used for fitting"
        plotter.transparent = True

        # Add the wavelength grid
        plotter.add_wavelength_grid(grid, "fitting simulations")

        # Add MAPPINGS SFR SED
        mappings_sed = load_example_mappings_sed()
        plotter.add_sed(mappings_sed, "MAPPINGS")

        # Add Bruzual-Charlot stellar SED
        bc_sed = load_example_bruzualcharlot_sed()
        plotter.add_sed(bc_sed, "Bruzual-Charlot")

        # Add Zubko dust emission SED
        zubko_sed = load_example_zubko_sed()
        plotter.add_sed(zubko_sed, "Zubko")

        # Add emission lines
        emission_lines = EmissionLines()
        for line in emission_lines: plotter.add_emission_line(line)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, filename)

        # Run the plotter
        plotter.run(path, min_wavelength=grid.min_wavelength, max_wavelength=grid.max_wavelength)

    # -----------------------------------------------------------------

    def plot_dust_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the dust grids ...")

        # Loop over the different dust grids
        index = 0
        for path in fs.directories_in_path(self.fit_dust_grids_path, sort=int):

            ski_path = fs.join(path, self.galaxy_name + ".ski")
            simulation = SkirtSimulation(ski_path=ski_path, outpath=path)

            # Plot the grid
            plotgrids(simulation, output_path=self.plot_fitting_dust_grids_path, silent=True, prefix=str(index))

            # Increment the index
            index += 1

    # -----------------------------------------------------------------

    def plot_dust_cell_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the dust cell distribution ...")

        # Loop over the distributions
        index = 0
        for distribution in self.dust_cell_trees:

            title = "Dust cells in each tree level"
            path = fs.join(self.plot_fitting_dust_grids_path, "cells_tree_" + str(index) + ".pdf")

            # Plot
            distribution.plot(title=title, path=path)

            # Increment the index
            index += 1

    # -----------------------------------------------------------------

    def plot_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the runtimes ...")

        bins = 15

        # Loop over the different remote hosts
        for host_id in self.runtimes:

            # Loop over the different configurations (packages, parallelization)
            for packages, parallelization in self.runtimes[host_id]:

                runtimes = self.runtimes[host_id][packages, parallelization]

                distribution = Distribution.from_values(runtimes, bins)

                path = fs.join(self.plot_fitting_path, "runtimes_" + host_id + "_" + str(parallelization) + "_" + str(packages) + ".pdf")
                title = host_id + " " + str(parallelization) + " " + str(packages)

                # Plot the distribution
                distribution.plot(title=title, path=path)

    # -----------------------------------------------------------------

    def plot_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SEDs ...")

        # Plot model SEDs
        self.plot_model_seds()

        # Plot contributions to the SED
        self.plot_sed_contributions()

    # -----------------------------------------------------------------

    def plot_model_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SEDs of all the models ...")

        # Loop over the generations
        for generation_name in self.seds:

            # Create the SEDPlotter object
            plotter = SEDPlotter()

            # Add all model SEDs (these have to be plotted in gray)
            counter = 0
            for sed in self.seds[generation_name]:

                # Add the model SEDs
                plotter.add_modeled_sed(sed, str(counter), ghost=True)
                counter += 1

            # Add the 'best' model total SED
            if generation_name in self.sed_contributions: plotter.add_modeled_sed(self.sed_contributions[generation_name]["total"], "best") # this SED has to be plotted in black

            # Add the observed SED to the plotter
            plotter.add_observed_sed(self.observed_sed, "observation")

            # Determine the path to the SED plot file
            path = fs.join(self.plot_fitting_path, "model_seds_" + generation_name + ".pdf")

            # Run the plotter
            plotter.run(title=self.galaxy_name, output=path)

    # -----------------------------------------------------------------

    def plot_sed_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting SEDs of various contributions ...")

        # Loop over the generations
        for generation_name in self.sed_contributions:

            # Create the SEDPlotter object
            plotter = SEDPlotter(self.galaxy_name)

            # Loop over the contributions
            for contribution in self.sed_contributions[generation_name]:

                # Add the simulated SED to the plotter
                plotter.add_modeled_sed(self.sed_contributions[generation_name][contribution], contribution, residuals=(contribution=="total"))

            # Add the observed SED to the plotter
            plotter.add_observed_sed(self.observed_sed, "observation")

            # Determine the path to the SED plot file
            path = fs.join(self.plot_fitting_path, "sed_contributions_" + generation_name + ".pdf")

            # Run the plotter
            plotter.run(path)

    # -----------------------------------------------------------------

    def plot_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a grid with the observed, simulated and residual images ...")

        # Get the ellipse
        path = fs.join(self.truncation_path, "ellipse.reg")
        region = SkyRegion.from_file(path)
        ellipse = region[0]

        # Create the image grid plotter
        plotter = ResidualImageGridPlotter(title="Image residuals")

        # Create list of filter names sorted by increasing wavelength
        sorted_filter_names = sorted(self.observed_images.keys(), key=lambda key: self.observed_images[key].filter.pivotwavelength())

        # Loop over the filter names, add a row to the image grid plotter for each filter
        for filter_name in sorted_filter_names:

            observed = self.observed_images[filter_name]
            simulated = self.simulated_images[filter_name]

            plotter.add_row(observed, simulated, filter_name)

        # Set the bounding box for the plotter
        plotter.set_bounding_box(ellipse.bounding_box)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, "images.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_geometries(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the geometries of the model ...")

        # Create the geometry plotter
        plotter = GeometryPlotter()

        # Add the geometries
        for label in self.geometries: plotter.add_geometry(self.geometries[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, "geometries.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def create_animation(self):

        """
        This function ...
        :return:
        """

        return

        ## LOAD IMAGES:

        # Inform the user
        log.info("Loading the SED plot files ...")

        # Find all PNG files within the the fit/plot directory
        for path in fs.files_in_path(self.fit_plot_path, extension="png", recursive=True):

            # Load the image (as a NumPy array)
            image = imageio.imread(path)

            # Add the image to the list of frames
            self.frames.append(image)

        ## CREATE ANIMATION:

        # Create the animated GIF instance
        self.animation = Animation(self.frames)

        ## WRITE:

        # Inform the user
        log.info("Writing the GIF animation ...")

        # Determine the path to the animation file
        path = fs.join(self.plot_fitting_path, "fitting.gif")

        # Save the animation as a GIF file
        self.animation.save(path)

    # -----------------------------------------------------------------

    # Set the load functions
    load_functions["wavelength grids"] = load_wavelength_grids
    load_functions["distributions"] = load_distributions
    load_functions["transmission"] = load_transmission_curves
    load_functions["tree"] = load_dust_cell_trees
    load_functions["runtimes"] = load_runtimes
    load_functions["seds"] = load_seds
    load_functions["images"] = load_images
    load_functions["geometries"] = load_geometries

    # -----------------------------------------------------------------

    # Set the plot functions
    plot_functions["chi-squared"] = plot_chi_squared
    plot_functions["distributions"] = plot_distributions
    plot_functions["wavelength grids"] = plot_wavelengths
    plot_functions["dust grids"] = plot_dust_grids
    plot_functions["cell distribution"] = plot_dust_cell_distribution
    #if self.plot_feature("runtimes") and self.runtimes is not None: self.plot_runtimes()
    #if self.plot_feature("sed") and len(self.seds) > 0: self.plot_model_seds()
    #if self.plot_feature("sed") and len(self.sed_contributions) > 0: self.plot_sed_contributions()
    #if self.plot_feature("images") and len(self.simulated_images) > 0: self.plot_images()
    #if self.plot_feature("geometries") and len(self.geometries) > 0: self.plot_geometries()
    plot_functions["runtimes"] = plot_runtimes
    plot_functions["seds"] = plot_seds
    plot_functions["images"] = plot_images
    plot_functions["geometries"] = plot_geometries
    plot_functions["animation"] = create_animation

# -----------------------------------------------------------------
