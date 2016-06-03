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
from collections import defaultdict

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
from ...core.plot.grids import plotgrids
from ...core.simulation.simulation import SkirtSimulation
from ...core.simulation.logfile import LogFile
from ...core.simulation.skifile import SkiFile

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

        # 2. Load the wavelength grid
        self.load_wavelength_grid()

        # 3. Load the transmission curves
        self.load_transmission_curves()

        # 4. Load the dust cell tree data
        self.load_dust_cell_tree()

        # 5. Load the runtimes
        self.load_runtimes()

        # Plot
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

    def plot(self):

        """
        This function ...
        :return:
        """

        # Plot the model components
        self.plot_components()

        # Plot the wavelength grid used for the fitting
        self.plot_wavelengths()

        # Plot the dust grid
        self.plot_dust_grid()

        # Plot the distribution of dust cells for the different tree levels
        self.plot_dust_cell_distribution()

        # Plot the distributions of the runtimes on different remote systems
        if self.runtimes is not None: self.plot_runtimes()

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

        # Create the transmission plotter
        plotter = TransmissionPlotter()

        plotter.title = "Wavelengths used for fitting"
        plotter.transparent = True

        # Add the transmission curves
        for label in self.transmission_curves: plotter.add_transmission_curve(self.transmission_curves[label], label)

        # Add the wavelength points
        for wavelength in self.wavelength_grid.wavelengths(): plotter.add_wavelength(wavelength)

        # Determine the path to the plot file
        path = fs.join(self.plot_fitting_path, "wavelengths.pdf")

        # Run the plotter
        plotter.run(path, min_wavelength=self.wavelength_grid.min_wavelength, max_wavelength=self.wavelength_grid.max_wavelength, min_transmission=0.0, max_transmission=1.05)

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
