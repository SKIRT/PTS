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

        # The wavelength grid
        self.wavelength_grid = None

        # The transmission curves
        self.transmission_curves = dict()

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

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    def load_wavelength_grid(self):

        """
        This function ...
        :return:
        """

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

    def plot(self):

        """
        This function ...
        :return:
        """

        # Plot
        self.plot_wavelengths()

        # Plot
        self.plot_dust_grid()

        # Plot
        self.plot_runtimes()

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
        plotgrids(simulation, output_path=self.plot_fitting_path)

    # -----------------------------------------------------------------

    def plot_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the runtimes ...")

        # Determine the path to the timing table
        timing_table_path = fs.join(self.fit_path, "timing.dat")

        # Load the timing table
        timing_table = TimingTable.read(timing_table_path)

        # Keep a list of all the runtimes recorded for a certain remote host
        runtimes_for_hosts = defaultdict(lambda: defaultdict(list))

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

            runtimes_for_hosts[host_id][packages, parallelization].append(runtime)

        # -----------------------------------------------------------------

        bins = 25

        # Loop over the different remote hosts
        for host_id in runtimes_for_hosts:

            # Loop over the different configurations (packages, parallelization)
            for packages, parallelization in runtimes_for_hosts[host_id]:

                runtimes = runtimes_for_hosts[host_id][packages, parallelization]

                distribution = Distribution.from_values(runtimes, 15)

                path = fs.join(self.plot_fitting_path, "runtimes_" + host_id + "_" + str(parallelization) + "_" + str(packages) + ".pdf")
                title = host_id + " " + str(parallelization) + " " + str(packages)

                # Plot the distribution
                distribution.plot(title=title, path=path)

# -----------------------------------------------------------------
