#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.analysis Contains the AnalysisPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict, OrderedDict

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ..analysis.component import AnalysisComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.basics.distribution import Distribution
from ...core.launch.timing import TimingTable
from ..core.transmission import TransmissionCurve
from ...core.plot.transmission import TransmissionPlotter
from ...core.plot.grids import plotgrids
from ...core.simulation.simulation import SkirtSimulation
from ...core.simulation.logfile import LogFile
from ..core.sed import SED, ObservedSED
from ...core.plot.sed import SEDPlotter

# -----------------------------------------------------------------

class AnalysisPlotter(PlottingComponent, AnalysisComponent):
    
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

        # Call the constructors of the base classes
        PlottingComponent.__init__(self, config)
        AnalysisComponent.__init__(self, config)

        # -- Attributes --

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

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the wavelength grid
        self.load_wavelength_grid()

        # 3. Load the transmission curves
        self.load_transmission_curves()

        # 4. Load the dust cell tree data
        self.load_dust_cell_tree()

        # 5. Load the runtimes
        self.load_runtimes()

        # 6. Load the SEDs
        self.load_seds()

        # 7. Plot
        self.plot()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing ...")

        # Set default values for all attributes
        self.wavelength_grid = None
        self.transmission_curves = dict()
        self.cell_distribution = None
        self.runtimes = None
        self.seds = OrderedDict()
        self.observed_sed = None

    # -----------------------------------------------------------------

    def load_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the wavelength grid ...")

        # Determine the path to the wavelength grid file
        path = fs.join(self.analysis_path, "in", "wavelengths.txt")

        # Load the wavelength grid
        if fs.is_file(path): self.wavelength_grid = WavelengthGrid.from_skirt_input(path)

    # -----------------------------------------------------------------

    def load_transmission_curves(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the transmission curves ...")

        # Load the observed SED
        sed = ObservedSED.from_file(self.observed_sed_path)

        # Loop over all filters for the points in the SED
        for fltr in sed.filters():

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

        # Determine the path to the log file of the best model simulation
        analysis_out_path = fs.join(self.analysis_path, "out")
        log_file_path = fs.join(analysis_out_path, self.galaxy_name + "_log.txt")

        # Break when the log file is not (yet) present
        if not fs.is_file(log_file_path): return

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
        timing_table_path = fs.join(self.analysis_path, "timing.dat")

        # Break if the timing table is not (yet) present
        if not fs.is_file(timing_table_path): return

        # Load the timing table
        timing_table = TimingTable.read(timing_table_path)

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

            # Add the runtime
            self.runtimes[host_id][packages, parallelization].append(runtime)

    # -----------------------------------------------------------------

    def load_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SEDs ...")

        # Determine the path to the analysis/heating directory
        analysis_heating_path = fs.join(self.analysis_path, "heating")

        # Load the SED of the simulation with the total stellar population
        total_out_path = fs.join(analysis_heating_path, "total", "out")
        total_sed_path = fs.join(total_out_path, self.galaxy_name + "_earth_sed.dat")
        if fs.is_file(total_sed_path):
            # Add the SED to the dictionary of SEDs
            total_sed = SED.from_skirt(total_sed_path)
            self.seds["total"] = total_sed

        # Add the SEDs of the simulations with the individual stellar populations
        for contribution in ["old", "young", "ionizing"]:

            # Determine the output path for the corresponding simulation
            out_path = fs.join(analysis_heating_path, contribution, "out")

            # Determine the path to the SED file
            sed_path = fs.join(out_path, self.galaxy_name + "_earth_sed.dat")

            # Skip if the file is not (yet) present
            if not fs.is_file(sed_path): continue

            # Load the SED
            sed = SED.from_skirt(sed_path)

            # Add the SED to the dictionary of SEDs
            self.seds[contribution] = sed

        # Break if no simulated SEDs were found
        if len(self.seds) == 0: return

        # Load the observed SED
        self.observed_sed = ObservedSED.from_file(self.observed_sed_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the wavelength grid
        if self.wavelength_grid is not None: self.plot_wavelengths()

        # Plot the dust grid
        self.plot_dust_grid()

        # Plot the distribution of dust cells for the different tree levels
        if self.cell_distribution is not None: self.plot_dust_cell_distribution()

        # Plot the distributions of the runtimes on different remote systems
        if self.runtimes is not None: self.plot_runtimes()

        # Plot the SED of the various stellar contributions
        if self.observed_sed is not None: self.plot_sed_contributions()

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

        plotter.title = "Wavelengths used for analysis"
        plotter.transparent = True

        # Add the transmission curves
        for label in self.transmission_curves: plotter.add_transmission_curve(self.transmission_curves[label], label)

        # Add the wavelength points
        for wavelength in self.wavelength_grid.wavelengths(): plotter.add_wavelength(wavelength)

        # Determine the path to the plot file
        path = fs.join(self.plot_analysis_path, "wavelengths.pdf")

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

        # Create a SkirtSimulation instance for the best model simulation
        analysis_out_path = fs.join(self.analysis_path, "out")
        ski_path = fs.join(self.analysis_path, self.galaxy_name + ".ski")
        if not fs.is_directory(analysis_out_path) or fs.is_empty(analysis_out_path): return # if the output path is non-existent or empty, break
        simulation = SkirtSimulation(ski_path=ski_path, outpath=analysis_out_path)

        # Plot the grid
        plotgrids(simulation, output_path=self.plot_analysis_path, silent=True)

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
        path = fs.join(self.plot_analysis_path, "cells_tree.pdf")

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

                # Get the runtimes
                runtimes = self.runtimes[host_id][packages, parallelization]

                distribution = Distribution.from_values(runtimes, bins)

                path = fs.join(self.plot_analysis_path, "runtimes_" + host_id + "_" + str(parallelization) + "_" + str(packages) + ".pdf")
                title = host_id + " " + str(parallelization) + " " + str(packages)

                # Plot the distribution
                distribution.plot(title=title, path=path)

    # -----------------------------------------------------------------

    def plot_sed_contributions(self):

        """
        This function ...
        :param self:
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
        path = fs.join(self.plot_analysis_path, "sed_contributions.pdf")

        # Run the plotter
        plotter.run(path)

# -----------------------------------------------------------------
