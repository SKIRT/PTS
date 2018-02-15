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

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..basics.table import SmartTable
from ..tools import tables, time
from ..simulation.simulation import SkirtSimulation, RemoteSimulation
from ..basics.log import log
from ..tools import sequences

# -----------------------------------------------------------------

all_ski_parameters = ["Wavelengths", "Dust cells", "Grid type", "Min level", "Max level", "Search method",
                      "Sample count", "Max optical depth", "Max mass fraction", "Max density dispersion",
                      "Self-absorption", "Transient heating", "Number of pixels"]

# -----------------------------------------------------------------

class MemoryTable(SmartTable):

    """
    This class ...
    """

    _column_info = OrderedDict()
    _column_info["Simulation name"] = (str, None, "name of the simulation")
    _column_info["Timestamp"] = (str, None, "timestamp")
    _column_info["Host id"] = (str, None, "remote host ID")
    _column_info["Cluster name"] = (str, None, "remote cluster name")
    _column_info["Cores"] = (int, None, "number of cores")
    _column_info["Threads per core"] = (int, None, "number of threads per core")
    _column_info["Processes"] = (int, None, "number of processes")
    _column_info["Wavelengths"] = (int, None, "number of wavelengths")
    _column_info["Dust cells"] = (int, None, "number of dust cells")
    _column_info["Grid type"] = (str, None, "type of grid")
    _column_info["Min level"] = (int, None, "minimum division level for the tree")
    _column_info["Max level"] = (int, None, "maximum division level for the tree")
    _column_info["Search method"] = (str, None, "search method (TopDown, Neighbor or Bookkeeping)")
    _column_info["Sample count"] = (int, None, "sample count")
    _column_info["Max optical depth"] = (float, None, "maximum optical depth")
    _column_info["Max mass fraction"] = (float, None, "maximum mass fraction")
    _column_info["Max density dispersion"] = (float, None, "maximum density dispersion")
    _column_info["Self-absorption"] = (bool, None, "self-absorption enabled")
    _column_info["Transient heating"] = (bool, None, "transient (non-LTE) heating enabled")
    _column_info["Data-parallel"] = (bool, None, "data parallelization enabled")
    _column_info["Number of pixels"] = (int, None, "total number of spatial pixels for all instruments")
    _column_info["Total peak memory"] = (float, "Gbyte", "peak memory usage during total simulation")
    _column_info["Setup peak memory"] = (float, "Gbyte", "peak memory usage during setup")
    _column_info["Stellar emission peak memory"] = (float, "Gbyte", "peak memory usage during stellar emission")
    _column_info["Spectra calculation peak memory"] = (float, "Gbyte", "peak memory usage during spectra calculation")
    _column_info["Dust emission peak memory"] = (float, "Gbyte", "peak memory usage during dust emission")
    _column_info["Writing peak memory"] = (float, "Gbyte", "peak memory usage during writing")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MemoryTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_entry(self, name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths,
                  ncells, grid_type, min_level, max_level, search_method, sample_count, max_optical_depth,
                  max_mass_fraction, max_density_dispersion, selfabsorption, transient_heating, data_parallel,
                  npixels, total_peak_memory, setup_peak_memory, stellar_peak_memory, spectra_peak_memory, dust_peak_memory,
                  writing_peak_memory, original_simulation_name=None):

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
        :param original_simulation_name:
        :return:
        """

        # Set the values
        values = [name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths,
                  ncells, grid_type, min_level, max_level, search_method, sample_count, max_optical_depth,
                  max_mass_fraction, max_density_dispersion, selfabsorption, transient_heating, data_parallel,
                  npixels, total_peak_memory, setup_peak_memory, stellar_peak_memory, spectra_peak_memory, dust_peak_memory,
                  writing_peak_memory]

        # Name was changed because already present
        if original_simulation_name is not None:
            index = self.index_for_simulation(original_simulation_name)
            previous_values = self.get_row(index, add_units=False, as_list=True)
            if sequences.equal_sequences(previous_values[1:], values[1:]):
                log.warning("Entry for simulation '" + original_simulation_name + "' is already present in the table: ignoring ...")
                return False
            else: pass  # new values, must be different simulation

        # Add a row to the table
        self.add_row(values)

        # Added
        return True

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
            original_simulation_name = simulation_name
            simulation_name = time.unique_name(simulation_name)
            log.warning("Generating the unique name '" + simulation_name + "' for this simulation ...")
        else: original_simulation_name = None

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
            if ski.treegrid_notfile():

                min_level = ski.tree_min_level()
                max_level = ski.tree_max_level()
                search_method = ski.tree_search_method()
                sample_count = ski.tree_sample_count()
                max_optical_depth = ski.tree_max_optical_depth()
                max_mass_fraction = ski.tree_max_mass_fraction()
                max_dens_disp = ski.tree_max_dens_disp()

            # File tree grid
            elif ski.filetreegrid():

                min_level = max_level = None
                search_method = ski.tree_search_method()
                sample_count = max_optical_depth = max_mass_fraction = max_dens_disp = None

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
        added = self.add_entry(simulation_name, submitted_at, host_id, cluster_name, cores,
                       hyperthreads, processes, wavelengths, ncells, grid_type, min_level, max_level,
                       search_method, sample_count, max_optical_depth, max_mass_fraction, max_dens_disp,
                       selfabsorption, transient_heating, data_parallel, npixels, peak_memory_usage, setup_peak_memory,
                       stellar_peak_memory, spectra_peak_memory, dust_peak_memory, writing_peak_memory, original_simulation_name=original_simulation_name)

        # Return the unique simulation name
        if not added:
            if original_simulation_name is None: raise RuntimeError("Something went wrong")
            return original_simulation_name
        else: return simulation_name

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
        for parameter in all_ski_parameters:
            if not self.all_equal(parameter): parameters.append(str(parameter))
        return parameters

    # -----------------------------------------------------------------

    def equal_ski_parameters(self):

        """
        This function ...
        :return:
        """

        parameters = []
        for parameter in all_ski_parameters:
            if self.all_equal(parameter): parameters.append(str(parameter))
        return parameters

    # -----------------------------------------------------------------

    def has_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return simulation_name in self["Simulation name"]

    # -----------------------------------------------------------------

    def remove_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        self.remove_row(index)

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

    def total_memory_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Total peak memory", index)

    # -----------------------------------------------------------------

    def setup_memory_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Setup peak memory", index)

    # -----------------------------------------------------------------

    def stellar_memory_for_simulation(self, simulation_name):

        """
        This fnuction ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Stellar emission peak memory", index)

    # -----------------------------------------------------------------

    def spectra_memory_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Spectra calculation peak memory", index)

    # -----------------------------------------------------------------

    def dust_memory_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Dust emission peak memory", index)

    # -----------------------------------------------------------------

    def writing_memory_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Writing peak memory", index)

    # -----------------------------------------------------------------

    def ski_parameters_for_simulation(self, simulation_name, which="different"):

        """
        This function ...
        :param simulation_name:
        :param which:
        :return:
        """

        # Find index of the simulation
        index = self.index_for_simulation(simulation_name)

        # Initialize dictionary
        parameters = dict()

        # Get the parameter names
        if which == "all": parameter_names = all_ski_parameters
        elif which == "different": parameter_names = self.different_ski_parameters()
        elif which == "equal": parameter_names = self.equal_ski_parameters()
        else: raise ValueError("Invalid value for 'which'")

        # Set the parameter values
        for parameter in parameter_names:
            parameters[str(parameter)] = self[parameter][index] if not self[parameter].mask[index] else None

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
