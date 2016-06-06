#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.fitting Contains the FittingPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict, OrderedDict

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.basics.distribution import Distribution
from ...core.launch.timing import TimingTable
from ...core.tools import tables
from ...core.basics.filter import Filter
from ..core.transmission import TransmissionCurve
from ...core.plot.transmission import TransmissionPlotter
from ...core.plot.wavelengthgrid import WavelengthGridPlotter
from ...core.plot.grids import plotgrids
from ...core.simulation.simulation import SkirtSimulation
from ...core.simulation.logfile import LogFile
from ...core.simulation.skifile import SkiFile
from ..core.emissionlines import EmissionLines
from ..core.sed import load_example_mappings_sed, load_example_bruzualcharlot_sed, load_example_zubko_sed
from ..core.sed import SED, ObservedSED
from ...magic.plot.imagegrid import ResidualImageGridPlotter
from ...magic.core.frame import Frame
from ...core.plot.sed import SEDPlotter
from ...magic.basics.skyregion import SkyRegion

# -----------------------------------------------------------------

class FittingPlotter(PlottingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(FittingPlotter, self).__init__(config)

        # -- Attributes --

        # The ski file
        self.ski = None

        # The wavelength grid
        self.wavelength_grid = None

        # The transmission curves
        self.transmission_curves = dict()

        # The distribution of dust cells per level
        self.cell_distribution = None

        # The runtimes
        self.runtimes = None

        # The SEDs of the different stellar contributions (total, old, young, ionizing)
        self.seds = OrderedDict()

        # The observed SED
        self.observed_sed = None

        # The simulated images
        self.simulated_images = dict()

        # The observed imags
        self.observed_images = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the ski file
        self.load_ski_file()

        # 3. Load the wavelength grid
        self.load_wavelength_grid()

        # 4. Load the transmission curves
        self.load_transmission_curves()

        # 5. Load the dust cell tree data
        self.load_dust_cell_tree()

        # 6. Load the runtimes
        self.load_runtimes()

        # 7. Load the SEDs for the various contribution
        self.load_seds()

        # 8. Load the simulated images
        self.load_images()

        # 9. Plot
        self.plot()

    # -----------------------------------------------------------------

    def load_ski_file(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file ...")

        # Determine the path to the initialized ski file
        path = fs.join(self.fit_path, self.galaxy_name + ".ski")

        # Load the ski file
        self.ski = SkiFile(path)

    # -----------------------------------------------------------------

    def load_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the wavelength grid used for the fitting ...")

        # Determine the path to the wavelength grid file
        path = fs.join(self.fit_path, "in", "wavelengths.txt")

        # Load the wavelength grid
        self.wavelength_grid = WavelengthGrid.from_skirt_input(path)

    # -----------------------------------------------------------------

    def load_transmission_curves(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the transmission curves for the observation filters ...")

        # Determine the path to the fluxes table
        fluxes_path = fs.join(self.phot_path, "fluxes.dat")

        # Load the fluxes table
        fluxes = tables.from_file(fluxes_path, format="ascii.ecsv")

        for i in range(len(fluxes)):

            instrument = fluxes["Instrument"][i]
            band = fluxes["Band"][i]

            # Constructor filter
            fltr = Filter.from_instrument_and_band(instrument, band)

            # Create the transmission curve
            transmission = TransmissionCurve.from_filter(fltr)

            # Normalize the transmission curve
            transmission.normalize(value=1.0, method="max")

            # Add the transmission curve to the dictionary
            self.transmission_curves[str(fltr)] = transmission

    # -----------------------------------------------------------------

    def load_dust_cell_tree(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the dust cell tree information ...")

        # Determine the path to the log file of the dust grid generating simulation
        fit_grid_path = fs.join(self.fit_path, "grid")
        log_file_path = fs.join(fit_grid_path, self.galaxy_name + "_log.txt")

        # Open the log file
        log_file = LogFile(log_file_path)

        # Get the distribution of cells per level of the tree
        self.cell_distribution = log_file.tree_leaf_distribution

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
            timing_table = TimingTable.read(timing_table_path)
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
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SEDs ...")

        # Determine the path to the fit/best directory
        fit_best_path = fs.join(self.fit_path, "best")

        # Determine the path to the fit/best/images directory
        images_path = fs.join(fit_best_path, "images")

        # Determine the path to the SED file from the 'images' simulation
        sed_path = fs.join(images_path, self.galaxy_name + "_earth_sed.dat")

        # Load the SED if the file exists
        if fs.is_file(sed_path):

            # Load the SED
            sed = SED.from_skirt(sed_path)

            # Add the total SED to the dictionary of SEDs
            self.seds["total"] = sed

        # Add the SEDs of the simulations with the individual stellar populations
        contributions = ["old", "young", "ionizing"]
        for contribution in contributions:

            # Determine the output path for this simulation
            out_path = fs.join(fit_best_path, contribution)

            # Determine the path to the SED file
            sed_path = fs.join(out_path, self.galaxy_name + "_earth_sed.dat")

            # Skip if the file does not exist
            if not fs.is_file(sed_path): continue

            # Load the SED
            sed = SED.from_skirt(sed_path)

            # Add the SED to the dictionary of SEDs
            self.seds[contribution] = sed

        # Determine the path to the observed SED
        observed_sed_path = fs.join(self.phot_path, "fluxes.dat")

        # Load the observed SED
        self.observed_sed = ObservedSED.from_file(observed_sed_path)

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
            filter_name = name.split("__")[1]

            # Open the image
            frame = Frame.from_file(path)

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
                log.warning(
                    "The simulated version of the " + name + " image could not be found, skipping " + name + " data ...")
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

    def plot(self):

        """
        This function ...
        :return:
        """

        # Plot the model components
        #self.plot_components()

        # Plot the wavelength grid used for the fitting
        self.plot_wavelengths()

        # Plot the dust grid
        self.plot_dust_grid()

        # Plot the distribution of dust cells for the different tree levels
        self.plot_dust_cell_distribution()

        # Plot the distributions of the runtimes on different remote systems
        if self.runtimes is not None: self.plot_runtimes()

        # Plot the SEDs
        if len(self.seds) > 0: self.plot_seds()

        # Plot the images
        if len(self.simulated_images) > 0: self.plot_images()

    # -----------------------------------------------------------------

    def plot_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the model components ...")


        for component_id in self.ski.get_stellar_component_ids():

            props = self.ski.get_stellar_component_properties(component_id)
            print(component_id, props)

        for component_id in self.ski.get_dust_component_ids():

            props = self.ski.get_dust_component_properties(component_id)
            print(component_id, props)

        exit()

        import matplotlib.pyplot as plt
        from matplotlib.patches import Ellipse


        #ells = [Ellipse(xy=rnd.rand(2) * 10, width=rnd.rand(), height=rnd.rand(), angle=rnd.rand() * 360)
        #        for i in range(NUM)]

        # Create the figure
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        for e in ells:
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(rnd.rand())
            e.set_facecolor(rnd.rand(3))

        #ax.set_xlim(0, 10)
        #ax.set_ylim(0, 10)

        plt.savefig()

    # -----------------------------------------------------------------

    def plot_wavelengths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelength grid ...")

        # Plot the wavelength grid with the observation filters
        self.plot_wavelengths_filters()

        # Plot the wavelengths with the SEDs of stars and dust
        self.plot_wavelengths_seds()

    # -----------------------------------------------------------------

    def plot_wavelengths_filters(self):

        """
        This function ...
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
        for wavelength in self.wavelength_grid.wavelengths(): plotter.add_wavelength(wavelength)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, "wavelengths_filters.pdf")

        # Run the plotter
        plotter.run(path, min_wavelength=self.wavelength_grid.min_wavelength,
                    max_wavelength=self.wavelength_grid.max_wavelength, min_transmission=0.0, max_transmission=1.05)

    # -----------------------------------------------------------------

    def plot_wavelengths_seds(self):

        """
        This function ...
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
        plotter.add_wavelength_grid(self.wavelength_grid, "fitting simulations")

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
        path = fs.join(self.plot_fitting_path, "wavelengths_seds.pdf")

        # Run the plotter
        plotter.run(path, min_wavelength=self.wavelength_grid.min_wavelength, max_wavelength=self.wavelength_grid.max_wavelength)

    # -----------------------------------------------------------------

    def plot_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the dust grid ...")

        # Create a SkirtSimulation instance for the grid generating simulation
        fit_grid_path = fs.join(self.fit_path, "grid")
        ski_path = fs.join(fit_grid_path, self.galaxy_name + ".ski")
        simulation = SkirtSimulation(ski_path=ski_path, outpath=fit_grid_path)

        # Plot the grid
        plotgrids(simulation, output_path=self.plot_fitting_path, silent=True)

    # -----------------------------------------------------------------

    def plot_dust_cell_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the dust cell distribution ...")

        # Set plot title and path
        title = "Dust cells in each tree level"
        path = fs.join(self.plot_fitting_path, "cells_tree.pdf")

        # Make the plot
        self.cell_distribution.plot(title=title, path=path)

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
        log.info("Plotting SEDs of various contributions ...")

        # Create the SEDPlotter object
        plotter = SEDPlotter(self.galaxy_name)

        # Loop over the simulated SEDs of the various stellar contributions
        for label in self.seds:

            # Add the simulated SED to the plotter
            plotter.add_modeled_sed(self.seds[label], label)

        # Add the observed SED to the plotter
        plotter.add_observed_sed(self.observed_sed, "observation")

        # Determine the path to the SED plot file
        path = fs.join(self.plot_fitting_path, "sed_contributions.pdf")

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
        sorted_filter_names = sorted(self.observed_images.keys(), key=lambda key: self.observed[key].filter.pivotwavelength())

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
