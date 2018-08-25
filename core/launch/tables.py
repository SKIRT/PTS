#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.tables Contains various table classes related to launching simulations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..remote.host import load_host
from ..tools import sequences
from ..tools.utils import lazyproperty, memoize_method, memoize_method_reset
from ..basics.table import SmartTable, property_type_to_builtin
from ..tools import tables
from ..simulation.remote import is_queued_status, is_finished_status, is_retrieved_status, is_running_status
from ..simulation.remote import is_analysed_status, is_aborted_status, is_cancelled_status, is_crashed_status, is_unknown_status
from ..simulation.simulation import SkirtSimulation, RemoteSimulation
from .options import AnalysisOptions
from ..simulation.definition import SingleSimulationDefinition
from ..basics.configuration import parent_type

# -----------------------------------------------------------------

class SimulationStatusTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Simulation name"] = (str, None, "name of the simulation")
    _column_info["Status"] = (str, None, "status of the simulation")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SimulationStatusTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_row(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        raise RuntimeError("Cannot add rows to a simulation status table")

    # -----------------------------------------------------------------

    def remove_row(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        raise RuntimeError("Cannot remove rows from a simulation status table")

    # -----------------------------------------------------------------

    @classmethod
    def from_previous(cls, table, new_statuses=None, remove_simulations=None, new_simulations=None):

        """
        This function ...
        :param table:
        :param new_statuses:
        :param remove_simulations:
        :param new_simulations:
        :return:
        """

        # Create lists
        simulation_names = []
        status_list = []

        # Loop over the existing simulations
        for simulation_name in table.simulation_names:

            # Remove simulation?
            if remove_simulations is not None and simulation_name in remove_simulations: continue

            # Get the correct status
            if new_statuses is not None and simulation_name in new_statuses:
                status = new_statuses[simulation_name]
            else:
                status = table.get_status(simulation_name)

            # Add to columns
            simulation_names.append(simulation_name)
            status_list.append(status)

        # Loop over the new simulations
        for simulation_name in new_simulations:

            # Get the status
            if isinstance(new_simulations, dict):
                status = new_simulations[simulation_name]
            elif new_statuses is not None and simulation_name in new_statuses:
                status = new_statuses[simulation_name]
            else:
                status = None

            # Add to columns
            simulation_names.append(simulation_name)
            status_list.append(status)

        # Create new status table (because status table class is full with lazyproperties and memoized methods)
        return cls.from_columns(simulation_names, status_list)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return list(self["Simulation name"])

    # -----------------------------------------------------------------

    def has_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return simulation_name in self.simulation_names

    # -----------------------------------------------------------------

    @memoize_method
    def index_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.simulation_names.index(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method_reset
    def get_status(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        return self.get_value("Status", index)

    # -----------------------------------------------------------------

    def set_status(self, simulation_name, simulation_status):

        """
        This function ...
        :param simulation_name:
        :param simulation_status:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        self.set_value("Status", index, simulation_status)

    # -----------------------------------------------------------------

    @lazyproperty
    def nsimulations(self):

        """
        This function ...
        :return:
        """

        return len(self.simulation_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def finished_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for simulation_name in self.simulation_names:
            if self.is_finished(simulation_name): names.append(simulation_name)
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def nfinished(self):

        """
        This function ...
        :return:
        """

        return len(self.finished_names)

    # -----------------------------------------------------------------

    @property
    def relative_nfinished(self):

        """
        This function ...
        :return:
        """

        return float(self.nfinished) / float(self.nsimulations)

    # -----------------------------------------------------------------

    @property
    def percentage_nfinished(self):

        """
        This function ...
        :return:
        """

        return self.relative_nfinished * 100.

    # -----------------------------------------------------------------

    @lazyproperty
    def retrieved_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for simulation_name in self.simulation_names:
            if self.is_retrieved(simulation_name): names.append(simulation_name)
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def nretrieved(self):

        """
        This function ...
        :return:
        """

        return len(self.retrieved_names)

    # -----------------------------------------------------------------

    @property
    def relative_nretrieved(self):

        """
        This function ...
        :return:
        """

        return float(self.nretrieved) / float(self.nsimulations)

    # -----------------------------------------------------------------

    @property
    def percentage_nretrieved(self):

        """
        This function ...
        :return:
        """

        return self.relative_nretrieved * 100.

    # -----------------------------------------------------------------

    @lazyproperty
    def analysed_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for simulation_name in self.simulation_names:
            if self.is_analysed(simulation_name): names.append(simulation_name)
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def nanalysed(self):

        """
        This function ...
        :return:
        """

        return len(self.analysed_names)

    # -----------------------------------------------------------------

    @property
    def relative_nanalysed(self):

        """
        This function ...
        :return:
        """

        return float(self.nanalysed) / float(self.nsimulations)

    # -----------------------------------------------------------------

    @property
    def percentage_nanalysed(self):

        """
        This function ...
        :return:
        """

        return self.relative_nanalysed * 100.

    # -----------------------------------------------------------------

    @lazyproperty
    def running_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for simulation_name in self.simulation_names:
            if self.is_running(simulation_name): names.append(simulation_name)
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def nrunning(self):
        return len(self.running_names)

    # -----------------------------------------------------------------

    @property
    def relative_nrunning(self):
        return float(self.nrunning) / float(self.nsimulations)

    # -----------------------------------------------------------------

    @property
    def percentage_nrunning(self):
        return self.relative_nrunning * 100.

    # -----------------------------------------------------------------

    @memoize_method_reset
    def is_queued(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        status = self.get_status(simulation_name)
        return is_queued_status(status)

    # -----------------------------------------------------------------

    @memoize_method_reset
    def is_running(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        status = self.get_status(simulation_name)
        return is_running_status(status)

    # -----------------------------------------------------------------

    def is_queued_or_running(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.is_queued(simulation_name) or self.is_running(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method_reset
    def is_finished(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        status = self.get_status(simulation_name)
        return is_finished_status(status)

    # -----------------------------------------------------------------

    def is_running_or_finished(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.is_running(simulation_name) or self.is_finished(simulation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def aborted_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for simulation_name in self.simulation_names:
            if self.is_aborted(simulation_name): names.append(simulation_name)
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def cancelled_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for simulation_name in self.simulation_names:
            if self.is_cancelled(simulation_name): names.append(simulation_name)
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def crashed_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for simulation_name in self.simulation_names:
            if self.is_crashed(simulation_name): names.append(simulation_name)
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def failed_names(self):

        """
        This function ...
        :return:
        """

        return self.aborted_names + self.cancelled_names + self.crashed_names

    # -----------------------------------------------------------------

    @memoize_method_reset
    def is_retrieved(self, simulation_name):
        status = self.get_status(simulation_name)
        return is_retrieved_status(status)

    # -----------------------------------------------------------------

    @memoize_method_reset
    def is_analysed(self, simulation_name):
        status = self.get_status(simulation_name)
        return is_analysed_status(status)

    # -----------------------------------------------------------------

    @memoize_method_reset
    def is_aborted(self, simulation_name):
        status = self.get_status(simulation_name)
        return is_aborted_status(status)

    # -----------------------------------------------------------------

    @memoize_method_reset
    def is_cancelled(self, simulation_name):
        status = self.get_status(simulation_name)
        return is_cancelled_status(status)

    # -----------------------------------------------------------------

    @memoize_method_reset
    def is_crashed(self, simulation_name):
        status = self.get_status(simulation_name)
        return is_crashed_status(status)

    # -----------------------------------------------------------------

    def is_running_or_finished_or_aborted_or_crashed(self, simulation_name):
        return self.is_running(simulation_name) or self.is_finished(simulation_name) or self.is_aborted(simulation_name) or self.is_crashed(simulation_name)

    # -----------------------------------------------------------------

    def is_failed(self, simulation_name):
        return self.is_aborted(simulation_name) or self.is_cancelled(simulation_name) or self.is_crashed(simulation_name)

    # -----------------------------------------------------------------

    def has_unknown_status(self, simulation_name):
        status = self.get_status(simulation_name)
        return is_unknown_status(status)

    # -----------------------------------------------------------------

    def reset_for_simulation(self, simulation_name, status):

        """
        This function ...
        :param simulation_name:
        :param status:
        :return:
        """

        # Set the status
        self.set_status(simulation_name, status)

        # Reset properties
        self.reset_properties()

        # Reset methods
        #self.reset_methods(simulation_name)

    # -----------------------------------------------------------------

    def reset_methods(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # TODO: detect memoized methods automatically

        raise NotImplementedError("This function is not working yet")

        #print(dir(self.get_status))
        #call_function = self.get_status.__dir__["__call__"]
        #print(dir(call_function))
        #print(vars(self.get_status))
        #call_function = vars(self.get_status)["__call__"]
        #call_function = self.get_status.__getattr__("__call__")
        #print(call_function)

        # Reset memoized methods
        self.get_status._reset_for_args(simulation_name)
        self.is_queued._reset_for_args(simulation_name)
        self.is_running._reset_for_args(simulation_name)
        self.is_finished._reset_for_args(simulation_name)
        self.is_retrieved._reset_for_args(simulation_name)
        self.is_analysed._reset_for_args(simulation_name)
        self.is_aborted._reset_for_args(simulation_name)
        self.is_cancelled._reset_for_args(simulation_name)
        self.is_crashed._reset_for_args(simulation_name)

    # -----------------------------------------------------------------

    def reset_properties(self):

        """
        This function ...
        :return:
        """

        # TODO: detect lazyproperties automatically

        del self.simulation_names
        del self.finished_names
        del self.nfinished
        del self.retrieved_names
        del self.nretrieved
        del self.analysed_names
        del self.nanalysed

# -----------------------------------------------------------------

class SimulationAssignmentTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Simulation name"] = (str, None, "name of the simulation")
    _column_info["Host ID"] = (str, None, "remote host ID")
    _column_info["Cluster"] = (str, None, "cluster name")
    _column_info["ID"] = (int, None, "simulation ID")
    _column_info["Success"] = (bool, None, "succesfully launched")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SimulationAssignmentTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        This function ...
        :return:
        """

        return self.get_column_values("Simulation name")

    # -----------------------------------------------------------------

    def has_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return simulation_name in self.names

    # -----------------------------------------------------------------

    def has_host_id(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_host_id_for_simulation(simulation_name) is not None

    # -----------------------------------------------------------------

    def set_host_id(self, simulation_name, host_id):

        """
        This function ...
        :param simulation_name:
        :param host_id:
        :return:
        """

        index = self.get_index_for_simulation(simulation_name)
        self.set_value("Host ID", index, host_id)

    # -----------------------------------------------------------------

    def has_id(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation_id_for_simulation(simulation_name) is not None

    # -----------------------------------------------------------------

    def set_id(self, simulation_name, id):

        """
        This function ...
        :param simulation_name:
        :param id:
        :return:
        """

        index = self.get_index_for_simulation(simulation_name)
        self.set_value("ID", index, id)

    # -----------------------------------------------------------------

    @property
    def ids(self):

        """
        This function ...
        :return:
        """

        return self.get_column_values("ID")

    # -----------------------------------------------------------------

    @property
    def host_ids(self):

        """
        This function ...
        :return:
        """

        return self.get_column_values("Host ID")

    # -----------------------------------------------------------------

    @property
    def unique_host_ids(self):

        """
        This function ...
        :return:
        """

        return sequences.unique_values(self.host_ids, ignore_none=True)

    # -----------------------------------------------------------------

    @property
    def hosts(self):

        """
        This function ...
        :return:
        """

        hosts = []
        for index in range(self.nsimulations):
            host = self.get_host_for_index(index)
            hosts.append(host)
        return hosts

    # -----------------------------------------------------------------

    @property
    def unique_hosts(self):

        """
        This function ...
        :return:
        """

        #print(self.hosts)
        return sequences.unique_values(self.hosts, ignore=[None])

    # -----------------------------------------------------------------

    def get_index_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.find_index(simulation_name, "Simulation name")

    # -----------------------------------------------------------------

    def get_simulation_name_for_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_value("Simulation name", index)

    # -----------------------------------------------------------------

    def get_index_for_simulation_id(self, simulation_id):

        """
        This function ...
        :param simulation_id:
        :return:
        """

        return tables.find_index(self, simulation_id, "ID")

    # -----------------------------------------------------------------

    def set_local_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get simulation index
        index = self.get_index_for_simulation(simulation_name)

        # Set values
        self.set_value("ID", index, None)
        self.set_value("Host ID", index, None)
        self.set_value("Cluster", index, None)

    # -----------------------------------------------------------------

    def set_host_id_for_simulation(self, simulation_name, host_id):

        """
        This function ...
        :param simulation_name:
        :param host_id:
        :return:
        """

        # Get simulation index
        index = self.get_index_for_simulation(simulation_name)

        # Set value
        self.set_value("Host ID", index, host_id)

    # -----------------------------------------------------------------

    def set_cluster_name_for_simulation(self, simulation_name, cluster_name):

        """
        This function ...
        :param simulation_name:
        :param cluster_name:
        :return:
        """

        # Get simulation index
        index = self.get_index_for_simulation(simulation_name)

        # Set value
        self.set_value("Cluster", index, cluster_name)

    # -----------------------------------------------------------------

    def set_host_for_simulation(self, simulation_name, host_id, cluster_name=None):

        """
        This function ...
        :param simulation_name:
        :param host_id:
        :param cluster_name:
        :return:
        """

        # Get index
        index = self.get_index_for_simulation(simulation_name)

        # Set values
        self.set_value("Host ID", index, host_id)
        self.set_value("Cluster", index, cluster_name)

    # -----------------------------------------------------------------

    def set_id_for_simulation(self, simulation_name, id):

        """
        This function ...
        :param simulation_name:
        :param id:
        :return:
        """

        # Get index
        index = self.get_index_for_simulation(simulation_name)

        # Set value
        self.set_value("ID", index, id)

    # -----------------------------------------------------------------

    def set_id_and_host_for_simulation(self, simulation_name, id, host_id, cluster_name=None):

        """
        This function ...
        :param simulation_name:
        :param id:
        :param host_id:
        :param cluster_name:
        :return:
        """

        from .batch import MissingSimulation

        # Get the index
        index = self.get_index_for_simulation(simulation_name)
        if index is None: raise MissingSimulation(simulation_name, "the assignment table")

        # Set values
        self.set_value("ID", index, id)
        self.set_value("Host ID", index, host_id)
        self.set_value("Cluster", index, cluster_name)

    # -----------------------------------------------------------------

    def update_simulation(self, simulation, success=None):

        """
        This function ...
        :param simulation:
        :param success:
        :return:
        """

        # Set
        self.set_id_and_host_for_simulation(simulation.name, simulation.id, simulation.host_id, simulation.cluster_name)
        if success is not None: self.set_success_for_simulation(simulation.name, success)

    # -----------------------------------------------------------------

    def add_or_update_simulation(self, simulation, success=None):

        """
        This function ...
        :param simulation:
        :param success:
        :return:
        """

        if self.has_simulation(simulation.name): self.update_simulation(simulation, success=success)
        else: self.add_simulation_object(simulation, success=success)

    # -----------------------------------------------------------------

    def get_host_id_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        from .batch import MissingSimulation

        index = self.get_index_for_simulation(simulation_name)
        if index is None: raise MissingSimulation(simulation_name, missing_from="the assignment table")
        return self.get_host_id_for_index(index)

    # -----------------------------------------------------------------

    def get_host_id_for_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_value("Host ID", index)

    # -----------------------------------------------------------------

    def get_cluster_name_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.get_index_for_simulation(simulation_name)
        return self.get_cluster_name_for_index(index)

    # -----------------------------------------------------------------

    def get_cluster_name_for_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_value("Cluster", index)

    # -----------------------------------------------------------------

    def get_host_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.get_index_for_simulation(simulation_name)
        return self.get_host_for_index(index)
    
    # -----------------------------------------------------------------

    def get_host_for_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """
        
        host_id = self.get_value("Host ID", index)
        cluster_name = self.get_value("Cluster", index)

        # Return None or host
        if host_id is None: return None
        else: return load_host(host_id, clustername=cluster_name)

    # -----------------------------------------------------------------

    def get_simulation_id_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        from .batch import MissingSimulation
        index = self.get_index_for_simulation(simulation_name)
        if index is None: return MissingSimulation(simulation_name, missing_from="assignment table")
        return self.get_value("ID", index)

    # -----------------------------------------------------------------

    def get_simulation_name_for_id(self, simulation_id):

        """
        This function ...
        :param simulation_id:
        :return:
        """

        index = self.get_index_for_simulation_id(simulation_id)
        return self.get_value("Simulation name", index)

    # -----------------------------------------------------------------

    def is_launched(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.get_index_for_simulation(simulation_name)
        return self.get_value("Success", index)

    # -----------------------------------------------------------------

    def is_local(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.is_masked_value("Host ID", index)

    # -----------------------------------------------------------------

    def is_local_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.get_index_for_simulation(simulation_name)
        return self.is_local(index)

    # -----------------------------------------------------------------

    def is_remote(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return not self.is_masked_value("Host ID", index)

    # -----------------------------------------------------------------

    def is_remote_host(self, index, host_id):

        """
        Thisf unction ....
        :param index:
        :param host_id:
        :return:
        """

        return self.get_value("Host ID", index) == host_id

    # -----------------------------------------------------------------

    def is_remote_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.get_index_for_simulation(simulation_name)
        return self.is_remote(index)

    # -----------------------------------------------------------------

    def set_success_for_simulation(self, simulation_name, value=True):

        """
        This function ...
        :param simulation_name:
        :param value:
        :return:
        """

        index = self.get_index_for_simulation(simulation_name)
        self.set_success(index, value=value)
        return index

    # -----------------------------------------------------------------

    def set_success(self, index, value=True):

        """
        Thisf unction ...
        :param index:
        :param value:
        :return:
        """

        self.set_value("Success", index, value)

    # -----------------------------------------------------------------

    @property
    def nsimulations(self):

        """
        This function ...
        :return:
        """

        return len(self)

    # -----------------------------------------------------------------

    def get_name(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_value("Simulation name", index)

    # -----------------------------------------------------------------

    def get_id(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_value("ID", index)

    # -----------------------------------------------------------------

    def get_host_id(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_value("Host ID", index)

    # -----------------------------------------------------------------

    @property
    def local_simulations(self):

        """
        Thisf unction ...
        :return:
        """

        # Keep list of simulation names
        names = []

        # Loop over the indices
        for index in range(self.nsimulations):
            if not self.is_local(index): continue
            name = self.get_name(index)
            names.append(name)

        # Return the list of names
        return names

    # -----------------------------------------------------------------

    @property
    def nlocal_simulations(self):

        """
        This function ...
        :return:
        """

        return len(self.local_simulations)

    # -----------------------------------------------------------------

    @property
    def remote_simulations(self):

        """
        This function ...
        :return:
        """

        # Keep list of simulation names
        names = []

        # Loop over the indices
        for index in range(self.nsimulations):
            if not self.is_remote(index): continue
            name = self.get_name(index)
            names.append(name)

        # Return the list of names
        return names

    # -----------------------------------------------------------------

    @property
    def nremote_simulations(self):

        """
        This function ...
        :return:
        """

        return len(self.remote_simulations)

    # -----------------------------------------------------------------

    def simulations_for_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Keep list of simulation names
        names = []

        # Loop over the indices
        for index in range(self.nsimulations):
            if not self.is_remote(index): continue
            host_id_simulation = self.get_host_id(index)
            if host_id_simulation != host_id: continue
            name = self.get_name(index)
            names.append(name)

        # Return the list of names
        return names

    # -----------------------------------------------------------------

    def nsimulations_for_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return len(self.simulations_for_remote(host_id))

    # -----------------------------------------------------------------

    def names_and_ids_for_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        names = []
        ids = []

        for index in range(self.nsimulations):

            if not self.is_remote(index): continue
            if not self.is_remote_host(index, host_id): continue

            name = self.get_name(index)
            id = self.get_id(index)
            names.append(name)
            ids.append(id)

        # Return the names and ids
        return names, ids

    # -----------------------------------------------------------------

    def ids_for_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        ids = []

        for index in range(self.nsimulations):

            if not self.is_remote(index): continue
            if not self.is_remote_host(index, host_id): continue
            id = self.get_id(index)
            ids.append(id)

        # Return the simulation IDs
        return ids

    # -----------------------------------------------------------------

    def add_local_simulation(self, name, success):

        """
        This function ...
        :param name:
        :param success:
        :return:
        """

        if name in self.names: raise ValueError("Already a simulation with the name '" + name + "'")
        values = [name, None, None, None, success]
        self.add_row(values)

    # -----------------------------------------------------------------

    def add_remote_simulation(self, name, host_id, cluster_name=None, simulation_id=None, success=True):

        """
        This function ...
        :param name:
        :param host_id:
        :param cluster_name:
        :param simulation_id:
        :param success:
        :return:
        """

        if name in self.names: raise ValueError("Already a simulation with the name '" + name + "'")
        values = [name, host_id, cluster_name, simulation_id, success]
        self.add_row(values)

    # -----------------------------------------------------------------

    def add_simulation(self, name, host_id=None, cluster_name=None, simulation_id=None, success=True):

        """
        This function ...
        :param name:
        :param host_id:
        :param cluster_name:
        :param simulation_id:
        :param success:
        :return:
        """

        if host_id is None:
            if cluster_name is not None: raise ValueError("Cannot pass cluster name if host ID is not defined")
            if simulation_id is not None: raise ValueError("Cannot pass simulation ID if host ID is not defined")
            self.add_local_simulation(name, success)
        else: self.add_remote_simulation(name, host_id, cluster_name=cluster_name, simulation_id=simulation_id, success=success)

    # -----------------------------------------------------------------

    def add_simulation_object(self, simulation, success=True):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Local simulation
        if isinstance(simulation, SkirtSimulation): self.add_local_simulation(simulation.name, success=success)

        # Remote simulation
        elif isinstance(simulation, RemoteSimulation):

            # Get simulation properties
            simulation_id = simulation.id
            host_id = simulation.host_id
            cluster_name = simulation.cluster_name

            # Add
            self.add_simulation(simulation.name, host_id=host_id, cluster_name=cluster_name, simulation_id=simulation_id, success=success)

        # Invalid
        else: raise ValueError("Invalid simulation object of type '" + str(type(simulation)) + "'")

# -----------------------------------------------------------------

class QueuedSimulationsTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Simulation name"] = (str, None, "name of the simulation")
    _column_info["Skifile path"] = (str, None, "path to the skifile")
    _column_info["Input path(s)"] = (str, None, "input filepaths / directory path")
    _column_info["Output path"] = (str, None, "output directory path")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(QueuedSimulationsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

        # Get list of lists of analysis property names
        analysis = AnalysisOptions()

        # Set property names
        property_names = analysis.ordered_property_names
        property_types = [analysis.get_ptype(name) for name in property_names]
        property_units = [analysis.unit_for_property(name) for name in property_names]
        property_descriptions = [analysis.description_for_property(name) for name in property_names]

        # Add properties of sections
        section_names = analysis.ordered_section_names
        #print(section_names)
        for section_name in section_names:
            section = analysis[section_name]

            for prop_name in section.ordered_property_names:
                #print(prop_name)

                ptype = section.get_ptype(prop_name)
                unit = section.unit_for_property(prop_name)
                description = section.description_for_property(prop_name)

                property_names.append(section_name + "/" + prop_name)
                property_types.append(ptype)
                property_units.append(unit)
                property_descriptions.append(description)

        # Add column info for analysis properties
        for name, ptype, unit, description in zip(property_names, property_types, property_units, property_descriptions):

            base_ptype = parent_type(ptype)
            dtype, _ = property_type_to_builtin(base_ptype)
            self.add_column_info(name, dtype, unit, description)

    # -----------------------------------------------------------------

    @classmethod
    def from_queue(cls, queue):

        """
        This function ...
        :param queue:
        :return:
        """

        # Create
        table = cls()

        # Setup
        table._setup()

        # Add the simulations
        for definition, name, analysis_options in queue: table.add_simulation(definition, name=name, analysis_options=analysis_options)

        # Return the table
        return table

    # -----------------------------------------------------------------

    def add_simulation(self, definition, name=None, analysis_options=None):

        """
        Thisf unction ...
        :param definition:
        :param name:
        :param analysis_options:
        :return:
        """

        # Check name
        if name is None: name = definition.name
        if name is None: raise ValueError("Name of the simulation is not defined in the definition")

        # Check type of definition
        if not isinstance(definition, SingleSimulationDefinition): raise ValueError("Definition should be single simulation definition")

        # Construct row
        input_string = str(definition.input_path) if definition.input_path is not None else None
        values = [name, definition.ski_path, input_string, definition.output_path]

        # Add analysis options
        for column_name, _, _, _ in self.column_info[4:]:

            # No analysis info?
            if analysis_options is None: value = None

            # Subsection property
            if "/" in column_name:
                section_name, prop_name = column_name.split("/")
                value = analysis_options[section_name][prop_name]
            # General property
            else: value = analysis_options[column_name]

            # Add the value
            values.append(value)

        # Add row
        self.add_row(values)

# -----------------------------------------------------------------
