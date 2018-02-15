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

# Import standard modules
from collections import defaultdict, OrderedDict

# Import the relevant PTS classes and modules
from ..basics.table import SmartTable
from ..tools import tables, time
from ..simulation.simulation import RemoteSimulation, SkirtSimulation
from ..basics.log import log
from ..simulation.discover import matching_npackages
from ..tools import sequences

# -----------------------------------------------------------------

all_ski_parameters = ["Wavelengths", "Packages", "Dust cells", "Grid type", "Min level", "Max level",
                          "Search method", "Sample count", "Max optical depth", "Max mass fraction",
                          "Max density dispersion", "Self-absorption", "Transient heating"]

# -----------------------------------------------------------------

class TimingTable(SmartTable):

    """
    This class ...
    """

    _column_info = OrderedDict()
    _column_info["Simulation name"] = (str, None, "Name of the simulation")
    _column_info["Timestamp"] = (str, None, "Timestamp")
    _column_info["Host id"] = (str, None, "Remote host ID")
    _column_info["Cluster name"] = (str, None, "Remote cluster name")
    _column_info["Cores"] = (int, None, "number of cores")
    _column_info["Threads per core"] = (int, None, "number of threads per core")
    _column_info["Processes"] = (int, None, "number of processes")
    _column_info["Wavelengths"] = (int, None, "number of wavelengths")
    _column_info["Packages"] = (int, None, "number of photon packages per wavelength")
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
    _column_info["Total runtime"] = (float, "s", "total simulation time")
    _column_info["Setup time"] = (float, "s", "time spent in simulation setup")
    _column_info["Stellar emission time"] = (float, "s", "time spent shooting stellar photon packages")
    _column_info["Spectra calculation time"] = (float, "s", "time spent in calculation of dust emission spectra")
    _column_info["Dust emission time"] = (float, "s", "time spent shooting dust emission photon packages")
    _column_info["Writing time"] = (float, "s", "time spent writing to disk")
    _column_info["Waiting time"] = (float, "s", "time spent waiting for other processes")
    _column_info["Communication time"] = (float, "s", "time spent in inter-process communication")
    _column_info["Dust densities communication time"] = (float, "s", "time spent in communication of dust densities")
    _column_info["Stellar absorption communication time"] = (float, "s", "time spent in communication of stellar absorption luminosities")
    _column_info["Dust absorption communication time"] = (float, "s", "time spent in communication of dust absorption luminosities")
    _column_info["Emission spectra communication time"] = (float, "s", "time spent in communication of emission spectra")
    _column_info["Instruments communication time"] = (float, "s", "time spent in communication of instrument data")
    _column_info["Intermediate time"] = (float, "s", "time spent in between other phases")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(TimingTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_entry(self, name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths,
                  packages, ncells, grid_type, min_level, max_level, search_method, sample_count, max_optical_depth,
                  max_mass_fraction, max_density_dispersion, selfabsorption, transient_heating, data_parallel,
                  total_runtime, setup_time, stellar_time, spectra_time, dust_time, writing_time, waiting_time,
                  communication_time, densities_communication_time, stellar_absorption_communication_time,
                  dust_absorption_communication_time, emission_communication_time, instruments_communication_time,
                  intermediate_time, original_simulation_name=None):

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
        :param densities_communication_time:
        :param stellar_absorption_communication_time:
        :param dust_absorption_communication_time:
        :param emission_communication_time:
        :param instruments_communication_time:
        :param intermediate_time:
        :param original_simulation_name:
        :return:
        """

        # Set the values
        values = [name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths, packages,
                  ncells, grid_type, min_level, max_level, search_method, sample_count, max_optical_depth,
                  max_mass_fraction, max_density_dispersion, selfabsorption, transient_heating, data_parallel,
                  total_runtime, setup_time, stellar_time, spectra_time, dust_time, writing_time, waiting_time,
                  communication_time, densities_communication_time, stellar_absorption_communication_time,
                  dust_absorption_communication_time, emission_communication_time, instruments_communication_time,
                  intermediate_time]

        # Name was changed because already present
        if original_simulation_name is not None:
            index = self.index_for_simulation(original_simulation_name)
            previous_values = self.get_row(index, add_units=False, as_list=True)
            if sequences.equal_sequences(previous_values[1:], values[1:]):
                log.warning("Entry for simulation '" + original_simulation_name + "' is already present in the table: ignoring ...")
                return False
            else: pass # new values, must be different simulation

        # Add a row to the table
        self.add_row(values)

        # Added
        return True

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

        # Check whether the name is unique
        if simulation_name in self["Simulation name"]:
            log.warning("A simulation with the name '" + simulation_name + "' is already present in this timing table")
            original_simulation_name = simulation_name
            simulation_name = time.unique_name(simulation_name)
            log.warning("Generating the unique name '" + simulation_name + "' for this simulation ...")
        else: original_simulation_name = None

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
            if ski.treegrid_notfile():

                min_level = ski.tree_min_level()
                max_level = ski.tree_max_level()
                search_method = ski.tree_search_method()
                sample_count = ski.tree_sample_count()
                max_optical_depth = ski.tree_max_optical_depth()
                max_mass_fraction = ski.tree_max_mass_fraction()
                max_dens_disp = ski.tree_max_dens_disp()

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
        densities_communication_time = timeline.communication_densities
        stellar_absorption_communication_time = timeline.communication_stellar_absorption
        dust_absorption_communication_time = timeline.communication_dust_absorption
        emission_communication_time = timeline.communication_emission
        instruments_communication_time = timeline.communication_instruments
        intermediate_time = timeline.other

        # Add an entry to the timing table
        added = self.add_entry(simulation_name, submitted_at, host_id, cluster_name, cores,
                       hyperthreads, processes, wavelengths, packages, ncells, grid_type, min_level, max_level,
                       search_method, sample_count, max_optical_depth, max_mass_fraction, max_dens_disp,
                       selfabsorption, transient_heating, data_parallel, total_runtime, setup_time,
                       stellar_time, spectra_time, dust_time, writing_time, waiting_time, communication_time,
                       densities_communication_time, stellar_absorption_communication_time,
                       dust_absorption_communication_time, emission_communication_time,
                       instruments_communication_time, intermediate_time, original_simulation_name=original_simulation_name)

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
            if not self.all_equal(parameter): parameters.append(str(parameter)) # dtype('S21') to str
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

    def total_time_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Total runtime", index)

    # -----------------------------------------------------------------

    def setup_time_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Setup time", index)

    # -----------------------------------------------------------------

    def stellar_emission_time_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Stellar emission time", index)

    # -----------------------------------------------------------------

    def spectra_calculation_time_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Spectra calculation time", index)

    # -----------------------------------------------------------------

    def dust_emission_time_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Dust emission time", index)

    # -----------------------------------------------------------------

    def writing_time_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Writing time", index)

    # -----------------------------------------------------------------

    def waiting_time_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Waiting time", index)

    # -----------------------------------------------------------------

    def communication_time_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Communication time", index)

    # -----------------------------------------------------------------

    def intermediate_time_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Intermediate time", index)

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

    def simulation_names_for_parameters(self, parameters):

        """
        This function ...
        :param parameters:
        :return:
        """

        return self["Simulation name"][self.indices_for_parameters(parameters)]

    # -----------------------------------------------------------------

    def even_npackages(self):

        """
        This function ...
        :return:
        """

        different_npackages = defaultdict(list)

        # Loop over the rows
        for i in range(len(self)):

            # Get the number of photon packages
            npackages = self["Packages"][i]

            # Get the number of processes
            nprocesses = self["Processes"][i]

            # Loop over the previous npackages
            for np in different_npackages:

                if matching_npackages(np, npackages, nprocesses):
                    different_npackages[np].append(i)
                    break

            # If break is not encountered, add as unique value
            else: different_npackages[npackages].append(i)

        for npackages in different_npackages:

            for index in different_npackages[npackages]:

                self["Packages"][index] = npackages

# -----------------------------------------------------------------
