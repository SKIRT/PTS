#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.fitting.chisquared Contains the ChiSquaredPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from .component import FittingPlottingComponent
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ....core.basics.distribution import Distribution
from ....core.launch.timing import TimingTable

# -----------------------------------------------------------------

class RuntimesPlotter(FittingPlottingComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param interactive:
        """

        # Call the constructor of the base class
        super(RuntimesPlotter, self).__init__(*args, **kwargs)

        # The runtimes
        self.runtimes = None

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

                distribution = Distribution.from_values("Runtime", runtimes, nbins=bins)

                path = fs.join(self.plot_fitting_path,
                               "runtimes_" + host_id + "_" + str(parallelization) + "_" + str(packages) + ".pdf")
                title = host_id + " " + str(parallelization) + " " + str(packages)

                # Plot the distribution
                distribution.plot(title=title, path=path)

# -----------------------------------------------------------------
