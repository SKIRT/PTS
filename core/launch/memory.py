#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.memory Contains the MemoryTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.table import SmartTable
from ..tools import tables, time
from ..simulation.simulation import SkirtSimulation, RemoteSimulation
from ..tools.logging import log
from ..units.parsing import parse_unit as u

# -----------------------------------------------------------------

class MemoryTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MemoryTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Simulation name", str, None, "name of the simulation")
        self.add_column_info("Timestamp", str, None, "timestamp")
        self.add_column_info("Host id", str, None, "remote host ID")
        self.add_column_info("Cluster name", str, None, "remote cluster name")
        self.add_column_info("Cores", int, None, "number of cores")
        self.add_column_info("Threads per core", int, None, "number of threads per core")
        self.add_column_info("Processes", int, None, "number of processes")
        self.add_column_info("Wavelengths", int, None, "number of wavelengths")
        self.add_column_info("Dust cells", int, None, "number of dust cells")
        self.add_column_info("Grid type", str, None, "type of grid")
        self.add_column_info("Min level", int, None, "minimum division level for the tree"),
        self.add_column_info("Max level", int, None, "maximum division level for the tree"),
        self.add_column_info("Search method", str, None, "search method (TopDown, Neighbor or Bookkeeping)")
        self.add_column_info("Sample count", int, None, "sample count")
        self.add_column_info("Max optical depth", float, None, "maximum optical depth")
        self.add_column_info("Max mass fraction", float, None, "maximum mass fraction")
        self.add_column_info("Max density dispersion", float, None, "maximum density dispersion")
        self.add_column_info("Self-absorption", bool, None, "self-absorption enabled")
        self.add_column_info("Transient heating", bool, None, "transient (non-LTE) heating enabled")
        self.add_column_info("Data-parallel", bool, None, "data parallelization enabled")
        self.add_column_info("Number of pixels", int, None, "total number of spatial pixels for all instruments")
        self.add_column_info("Total peak memory", float, u("Gbyte"), "peak memory usage during total simulation")
        self.add_column_info("Setup peak memory", float, u("Gbyte"), "peak memory usage during setup")
        self.add_column_info("Stellar emission peak memory", float, u("Gbyte"), "peak memory usage during stellar emission")
        self.add_column_info("Spectra calculation peak memory", float, u("Gbyte"), "peak memory usage during spectra calculation")
        self.add_column_info("Dust emission peak memory", float, u("Gbyte"), "peak memory usage during dust emission")
        self.add_column_info("Writing peak memory", float, u("Gbyte"), "peak memory usage during writing")

    # -----------------------------------------------------------------

    def add_entry(self, name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths,
                  ncells, grid_type, min_level, max_level, search_method, sample_count, max_optical_depth,
                  max_mass_fraction, max_density_dispersion, selfabsorption, transient_heating, data_parallel,
                  npixels, total_peak_memory, setup_peak_memory, stellar_peak_memory, spectra_peak_memory, dust_peak_memory,
                  writing_peak_memory):

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
        :param npixels
        :param total_peak_memory:
        :param setup_peak_memory:
        :param stellar_peak_memory:
        :param spectra_peak_memory:
        :param dust_peak_memory:
        :param writing_peak_memory:
        :return:
        """

        # Set the values
        values = [name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths,
                  ncells, grid_type, min_level, max_level, search_method, sample_count, max_optical_depth,
                  max_mass_fraction, max_density_dispersion, selfabsorption, transient_heating, data_parallel,
                  npixels, total_peak_memory, setup_peak_memory, stellar_peak_memory, spectra_peak_memory, dust_peak_memory,
                  writing_peak_memory]

        # Add a row to the table
        self.add_row(values)

    # -----------------------------------------------------------------

    def add_from_simulation(self, simulation, ski, log_file, parameters=None):

        """
        This function ...
        :param simulation:
        :param ski:
        :param log_file:
        :param parameters:
        :return:
        """

        # Get the simulation name
        simulation_name = simulation.name

        # Remote simulation
        if isinstance(simulation, RemoteSimulation):

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

            # Time of submitting
            submitted_at = None

            # Host etc.
            host_id = log_file.host
            cluster_name = None
            # cores = None
            # hyperthreads = None

            # Parallelization
            processes = log_file.processes
            threads = log_file.threads

            # We don't know how many threads actually ran per core, guess 1 so we can put a number on the number of cores
            cores = processes * threads
            hyperthreads = 1

        # Invalid argument
        else: raise ValueError("Invalid argument for 'simulation'")

        # Check whether the name is unique
        if simulation_name in self["Simulation name"]:
            log.warning("A simulation with the name '" + simulation_name + "' is already present in this memory table")
            simulation_name = time.unique_name(simulation_name)
            log.warning("Generating the unique name '" + simulation_name + "' for this simulation")

        # Get the peak memory usage
        peak_memory_usage = None
        try: peak_memory_usage = log_file.peak_memory
        except RuntimeError:
            for log_file in simulation.logfiles():
                try:
                    peak_memory_usage = log_file.peak_memory
                    break
                except RuntimeError:
                    pass
        if peak_memory_usage is None: raise RuntimeError("All log files were aborted")

        # Get the number of wavelengths
        wavelengths = log_file.wavelengths

        # Get the number of dust cells
        ncells = log_file.dust_cells

        # Check whether data parallelization was enabled for the simulation
        data_parallel = log_file.data_parallel

        if ski is not None:

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

            # Determine the total number of pixels from all the instruments defined in the ski file
            npixels = ski.nspatialpixels()

        elif parameters is not None:

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
            npixels = None

        else: raise ValueError("Ski file or parameters map must be specified")

        # Get additional memory info
        setup_peak_memory = log_file.setup_peak_memory
        stellar_peak_memory = log_file.stellar_peak_memory
        spectra_peak_memory = log_file.spectra_peak_memory
        dust_peak_memory = log_file.dust_peak_memory
        writing_peak_memory = log_file.writing_peak_memory

        # Add an entry to the memory table
        self.add_entry(simulation_name, submitted_at, host_id, cluster_name, cores,
                       hyperthreads, processes, wavelengths, ncells, grid_type, min_level, max_level,
                       search_method, sample_count, max_optical_depth, max_mass_fraction, max_dens_disp,
                       selfabsorption, transient_heating, data_parallel, npixels, peak_memory_usage, setup_peak_memory,
                       stellar_peak_memory, spectra_peak_memory, dust_peak_memory, writing_peak_memory)

        # Return the unique simulation name
        return simulation_name

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return list(self["Simulation name"])

    # -----------------------------------------------------------------

    def different_ski_parameters(self):

        """
        This function ...
        :return:
        """

        parameters = []
        ski_parameters = ["Wavelengths", "Dust cells", "Grid type", "Min level", "Max level", "Search method",
                          "Sample count", "Max optical depth", "Max mass fraction", "Max density dispersion",
                          "Self-absorption", "Transient heating", "Number of pixels"]
        for parameter in ski_parameters:
            if not self.all_equal(parameter): parameters.append(str(parameter))

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def index_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Find index of the simulation
        index = tables.find_index(self, simulation_name)
        return index

    # -----------------------------------------------------------------

    def ski_parameters_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Find index of the simulation
        index = self.index_for_simulation(simulation_name)

        # Initialize dictionary
        parameters = dict()

        # Set the parameter values
        for parameter in self.different_ski_parameters(): parameters[str(parameter)] = self[parameter][index] # dtype('S21') to str

        # Return the parameter values
        return parameters

    # -----------------------------------------------------------------

    def indices_for_parameters(self, parameters):

        """
        This function ...
        :return:
        """

        indices = []

        # Loop over the rows
        for index in range(len(self)):

            for label in parameters:

                if self[label][index] != parameters[label]: break

            # Break is not encountered: all parameters match for this row
            else: indices.append(index)

        return indices

# -----------------------------------------------------------------
