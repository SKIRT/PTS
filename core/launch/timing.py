#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.timing Contains the TimingTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.table import SmartTable
from ..simulation.simulation import RemoteSimulation, SkirtSimulation

# -----------------------------------------------------------------

class TimingTable(SmartTable):

    """
    This class ...
    """

    column_info = [("Simulation name", str, None, "Name of the simulation"),
                   ("Timestamp", str, None, "Timestamp"),
                   ("Host id", str, None, "Remote host ID"),
                   ("Cluster name", str, None, "Remote cluster name"),
                   ("Cores", int, None, "number of cores"),
                   ("Threads per core", int, None, "number of threads per core"),
                   ("Processes", int, None, "number of processes"),
                   ("Wavelengths", int, None, "number of wavelengths"),
                   ("Packages", int, None, "number of photon packages per wavelength"),
                   ("Dust cells", int, None, "number of dust cells"),
                   ("Grid type", str, None, "type of grid"),
                   ("Min level", int, None, "minimum division level for the tree"),
                   ("Max level", int, None, "maximum division level for the tree"),
                   ("Search method", str, None, "search method (TopDown, Neighbor or Bookkeeping)"),
                   ("Sample count", int, None, "sample count"),
                   ("Max optical depth", float, None, "maximum optical depth"),
                   ("Max mass fraction", float, None, "maximum mass fraction"),
                   ("Max density dispersion", float, None, "maximum density dispersion"),
                   ("Self-absorption", bool, None, "self-absorption enabled"),
                   ("Transient heating", bool, None, "transient (non-LTE) heating enabled"),
                   ("Data-parallel", bool, None, "data parallelization enabled"),
                   ("Total runtime", float, "s", "total simulation time"),
                   ("Setup time", float, "s", "time spent in simulation setup"),
                   ("Stellar emission time", float, "s", "time spent shooting stellar photon packages"),
                   ("Spectra calculation time", float, "s", "time spent in calculation of dust emission spectra"),
                   ("Dust emission time", float, "s", "time spent shooting dust emission photon packages"),
                   ("Writing time", float, "s", "time spent writing to disk"),
                   ("Waiting time", float, "s", "time spent waiting for other processes"),
                   ("Communication time", float, "s", "time spent in inter-process communication"),
                   ("Intermediate time", float, "s", "time spent in between other phases")]

    # -----------------------------------------------------------------

    def add_entry(self, name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths,
                  packages, ncells, grid_type, min_level, max_level, search_method, sample_count, max_optical_depth,
                  max_mass_fraction, max_density_dispersion, selfabsorption, transient_heating, data_parallel,
                  total_runtime, setup_time, stellar_time, spectra_time, dust_time, writing_time, waiting_time,
                  communication_time, intermediate_time):

        """
        This function ...
        :param name:
        :param timestamp:
        :param host_id:
        :param cluster_name:
        :param cores:
        :param threads_per_core:
        :param processes:
        :param wavelengths:
        :param packages:
        :param ncells:
        :param grid_type:
        :param min_level:
        :param max_level:
        :param search_method:
        :param sample_count:
        :param max_optical_depth:
        :param max_mass_fraction:
        :param max_density_dispersion:
        :param selfabsorption:
        :param transient_heating:
        :param data_parallel:
        :param total_runtime:
        :param setup_time:
        :param stellar_time:
        :param spectra_time:
        :param dust_time:
        :param writing_time:
        :param waiting_time:
        :param communication_time:
        :param intermediate_time:
        :return:
        """

        # Set the values
        values = [name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths, packages,
                  ncells, grid_type, min_level, max_level, search_method, sample_count, max_optical_depth,
                  max_mass_fraction, max_density_dispersion, selfabsorption, transient_heating, data_parallel,
                  total_runtime, setup_time, stellar_time, spectra_time, dust_time, writing_time, waiting_time,
                  communication_time, intermediate_time]

        # Add a row to the table
        self.add_row(values)

    # -----------------------------------------------------------------

    def add_from_simulation(self, simulation, ski, log_file, timeline, parameters=None):

        """
        This function ...
        :param simulation:
        :param ski:
        :param log_file:
        :param timeline:
        :param parameters:
        :return:
        """

        # Remote simulation
        if isinstance(simulation, RemoteSimulation):

            # Get the simulation name
            simulation_name = simulation.name

            # Time of submitting
            submitted_at = simulation.submitted_at

            # Get the name of the host on which the simulation was run
            host_id = simulation.host_id
            cluster_name = simulation.cluster_name

            # Get the parallelization object from the simulation
            parallelization = simulation.parallelization

            # Get the paralleliation properties
            cores = parallelization.cores
            hyperthreads = parallelization.threads_per_core
            processes = parallelization.processes

        # Basic simulation object
        elif isinstance(simulation, SkirtSimulation):

            simulation_name = simulation.prefix()

            # Time of submitting
            submitted_at = None

            # Host etc.
            host_id = log_file.host
            cluster_name = None
            #cores = None
            #hyperthreads = None

            # Parallelization
            processes = log_file.processes
            threads = log_file.threads

            # We don't know how many threads actually ran per core, guess 1 so we can put a number on the number of cores
            cores = processes * threads
            hyperthreads = 1

        # Invalid argument
        else: raise ValueError("Invalid argument for 'simulation'")

        # Get the total runtime (in seconds)
        total_runtime = log_file.total_runtime

        # Get the number of wavelengths
        wavelengths = log_file.wavelengths

        if ski is not None:

            # Get the number of photon packages
            packages = ski.packages()

            # Get the number of dust cells
            ncells = log_file.dust_cells

            # Get the dust grid type
            grid_type = ski.gridtype()

            # If the grid is a tree grid, get additional properties
            if ski.treegrid():

                min_level = ski.tree_min_level()
                max_level = ski.tree_max_level()
                search_method = ski.tree_search_method()
                sample_count = ski.tree_sample_count()
                max_optical_depth = ski.tree_max_optical_depth()
                max_mass_fraction = ski.tree_max_mass_fraction()
                max_dens_disp = ski.tree_max_dens_disp()

            # Else, set all properties to None
            else: min_level = max_level = search_method = sample_count = max_optical_depth = max_mass_fraction = max_dens_disp = None

            # Check whether dust self-absorption was enabled for the simulation
            selfabsorption = ski.dustselfabsorption()

            # Check whether transient heating was enabled for the simulation
            transient_heating = ski.transientheating()

        elif parameters is not None:

            packages = parameters.npackages
            ncells = parameters.ncells

            grid_type = None

            min_level = None
            max_level = None
            search_method = None
            sample_count = None
            max_optical_depth = None
            max_mass_fraction = None
            max_dens_disp = None

            selfabsorption = parameters.selfabsorption
            transient_heating = parameters.transient_heating

        else: raise ValueError("Ski file or parameters map must be specified")

        # Check whether data parallelization was enabled for the simulation
        data_parallel = log_file.data_parallel

        # Get the different contributions to the simulation's runtime
        setup_time = timeline.setup
        stellar_time = timeline.stellar
        spectra_time = timeline.spectra
        dust_time = timeline.dust
        writing_time = timeline.writing
        waiting_time = timeline.waiting
        communication_time = timeline.communication
        intermediate_time = timeline.other

        # Add an entry to the timing table
        self.add_entry(simulation_name, submitted_at, host_id, cluster_name, cores,
                               hyperthreads, processes, wavelengths, packages, ncells, grid_type, min_level, max_level,
                               search_method, sample_count, max_optical_depth, max_mass_fraction, max_dens_disp,
                               selfabsorption, transient_heating, data_parallel, total_runtime, setup_time,
                               stellar_time, spectra_time, dust_time, writing_time, waiting_time, communication_time,
                               intermediate_time)

# -----------------------------------------------------------------
