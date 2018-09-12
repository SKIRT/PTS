#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.remote Contains the SKIRTRemote class, used for launching, checking and retrieving
#  remote SKIRT simulations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..remote.remote import Remote, TimeOutReached
from .jobscript import MultiJobScript, SKIRTJobScript
from ..tools import time, introspection
from ..tools import filesystem as fs
from .simulation import RemoteSimulation
from ..basics.log import log
from ..launch.options import SchedulingOptions
from ..simulation.parallelization import Parallelization
from ..simulation.arguments import SkirtArguments
from ..basics.handle import ExecutionHandle
from .status import LogSimulationStatus, get_status_from_last_log_line
from .input import SimulationInput
from ..tools import types
from ..tools import numbers
from .output import get_output_type, get_parent_type
from .data import SimulationData
from .screen import ScreenScript
from ..tools.stringify import tostr

# -----------------------------------------------------------------

queued_name = "queued"
running_name = "running"
finished_name = "finished"
retrieved_name = "retrieved"
analysed_name = "analysed"
aborted_name = "aborted"
cancelled_name = "cancelled"
crashed_name = "crashed"
unknown_name = "unknown"
invalid_name = "invalid"

# -----------------------------------------------------------------

# started means that we can expect log output: also crashed, aborted etc.
def is_started(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return is_finished_status(stat) or is_running_status(stat) or is_aborted_status(stat) or is_crashed_status(stat)

# -----------------------------------------------------------------

def is_analysed_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    # completely analysed
    return stat == analysed_name

# -----------------------------------------------------------------

def is_analysing_or_analysed_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    # still analysing or completely analysed
    return stat.startswith(analysed_name)

# -----------------------------------------------------------------

def is_analysing_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return is_analysing_or_analysed_status(stat) and not is_analysed_status(stat)

# -----------------------------------------------------------------

def is_retrieved_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return stat == retrieved_name or is_analysing_or_analysed_status(stat)

# -----------------------------------------------------------------

def is_finished_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return stat == finished_name or is_retrieved_status(stat)

# -----------------------------------------------------------------

def is_queued_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return stat == queued_name

# -----------------------------------------------------------------

def is_running_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return running_name in stat

# -----------------------------------------------------------------

def is_unknown_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return stat == unknown_name

# -----------------------------------------------------------------

def is_invalid_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return invalid_name in stat

# -----------------------------------------------------------------

def is_invalid_or_unknown_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return is_invalid_status(stat) or is_unknown_status(stat)

# -----------------------------------------------------------------

def is_aborted_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return stat == aborted_name

# -----------------------------------------------------------------

def is_cancelled_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return stat == cancelled_name

# -----------------------------------------------------------------

def is_crashed_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return crashed_name in stat

# -----------------------------------------------------------------

def get_simulation_path_for_host(host_id, simulation_id):

    """
    This function ...
    :param host_id:
    :param simulation_id:
    :return:
    """

    return fs.join(introspection.skirt_run_dir, host_id, str(simulation_id) + ".sim")

# -----------------------------------------------------------------

def has_simulation_for_host(host_id, simulation_id):

    """
    This function ...
    :param host_id:
    :param simulation_id:
    :return:
    """

    return fs.is_file(get_simulation_path_for_host(host_id, simulation_id))

# -----------------------------------------------------------------

def get_simulation_for_host(host_id, simulation_id):

    """
    This function ...
    :param host_id:
    :param simulation_id:
    :return:
    """

    # Determine the simulation path
    simulation_path = get_simulation_path_for_host(host_id, simulation_id)
    if not fs.is_file(simulation_path): raise ValueError("Simulation file does not exist")

    # load the simulation and return
    return RemoteSimulation.from_file(simulation_path)

# -----------------------------------------------------------------

def get_simulation_paths_for_host(host_id, as_dict=False, names=None):

    """
    This function ...
    :param host_id:
    :param as_dict:
    :param names:
    :return:
    """

    # Names are specified or we have to return as dict
    if names is not None or as_dict:

        # Get simulation objects, already filter on host
        simulations = get_simulations_for_host(host_id, names=names)
        nsimulations = len(simulations)

        # As dictionary
        if as_dict:

            dictionary = OrderedDict((simulation.name, simulation.path) for simulation in simulations)
            if len(dictionary) < nsimulations: raise ValueError("Something went wrong: simulation names not unique")
            return dictionary

        # Return the list of paths
        else: return [simulation.path for simulation in simulations]

    # No filtering on names is necessary, nor putting into dictionary: just return the list of paths
    else: return [get_simulation_path_for_host(host_id, simulation_id) for simulation_id in introspection.simulation_ids_for_host(host_id)]

# -----------------------------------------------------------------

def get_simulations_for_host(host_id, as_dict=False, names=None):

    """
    This function ...
    :param host_id:
    :param as_dict:
    :param names:
    :return:
    """

    # Get simulations for this host
    simulations = [get_simulation_for_host(host_id, simulation_id) for simulation_id in introspection.simulation_ids_for_host(host_id)]

    # Filter
    if names is not None: simulations = [simulation for simulation in simulations if simulation.name in names]

    # Get the number of simulations
    nsimulations = len(simulations)

    # Return
    if as_dict:
        dictionary = OrderedDict((simulation.name, simulation) for simulation in simulations)
        if len(dictionary) < nsimulations: raise ValueError("Something went wrong: simulation names not unique")
        return dictionary
    else: return simulations

# -----------------------------------------------------------------

def get_simulation_names_for_host(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    return [simulation.name for simulation in get_simulations_for_host(host_id)]

# -----------------------------------------------------------------

def get_simulation_ids(host_id, names, as_dict=False):

    """
    This function ...
    :param host_id:
    :param names:
    :param as_dict:
    :return:
    """

    # Get the simulations for this host
    simulations = get_simulations_for_host(host_id, names=names)
    nnames = len(names)

    # Return the IDs
    if as_dict:
        dictionary = OrderedDict((simulation.name, simulation.id) for simulation in simulations)
        if len(dictionary) < nnames: raise ValueError("Something went wrong: simulation names not unique")
        return dictionary
    else: return [simulation.id for simulation in simulations]

# -----------------------------------------------------------------

def get_simulation_id(host_id, name):

    """
    This function ...
    :param host_id:
    :param name:
    :return:
    """

    # Loop over the simulation IDs
    for simulation_id in introspection.simulation_ids_for_host(host_id):

        # Load simulation
        simulation = get_simulation_for_host(host_id, simulation_id)

        # Check name
        if simulation.name == name: return simulation_id

    # None found
    raise ValueError("Simulation with name '" + name + "' does not exist for remote host '" + host_id + "'")

# -----------------------------------------------------------------

def needs_retrieval(filename, types):

    """
    This function ...
    :param filename:
    :param types: retrieve types
    :return:
    """

    # Get the output type
    output_type = get_output_type(filename)

    # Check whether in specified retrieval types
    if output_type in types: return True

    # Get parent types
    parent_type = get_parent_type(output_type)

    # Check if (different from type itself) parent type is in the specified retrieval types
    if parent_type != output_type and parent_type in types: return True

    # Otherwise, return False
    return False

# -----------------------------------------------------------------

def get_status_simulations(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    from ..basics.configuration import prompt_proceed

    # Initialize a list to contain the statuses
    entries = []

    # Detemrine directory path
    local_skirt_host_run_dir = fs.join(introspection.skirt_run_dir, host_id)

    # Search for simulation files in the local SKIRT run/host_id directory
    for path in fs.files_in_path(local_skirt_host_run_dir, extension="sim", sort=int):

        # Open the simulation file
        try: simulation = RemoteSimulation.from_file(path)
        except EOFError:
            log.error("The simulation file '" + path + "' is not readable: perhaps serialization has failed")
            remove = prompt_proceed("remove the simulation file?")
            if remove: fs.remove_file(path)
            log.error("Skipping this simulation ...")
            continue

        # Check whether the handle is defined
        if simulation.handle is None:

            # Warning to get attention
            log.warning("Simulation '" + simulation.name + "' [" + str(simulation.id) + "] on remote host '" + host_id + "' doesn't appear to have an execution handle. Assuming it is still running in attached mode through another terminal.")
            entries.append((path, running_name))

        # Handle is defined
        else:

            # GET STATUS

            # The name of the ski file (the simulation prefix)
            #ski_name = simulation.prefix()

            # Determine the path to the remote log file
            #remote_log_file_path = simulation.remote_log_file_path

            # Check whether the simulation has already been analysed
            if simulation.analysed: simulation_status = analysed_name

            # Partly analysed
            elif simulation.analysed_any:

                analysed = []
                if simulation.analysed_all_extraction: analysed.append("extraction")
                if simulation.analysed_all_plotting: analysed.append("plotting")
                if simulation.analysed_all_misc: analysed.append("misc")
                if simulation.analysed_batch: analysed.append("batch")
                if simulation.analysed_scaling: analysed.append("scaling")
                if simulation.analysed_all_extra: analysed.append("extra")

                if len(analysed) > 0: simulation_status = "analysed: " + ", ".join(analysed)
                else: simulation_status = "analysed: started"

            # Check whether the simulation has already been retrieved
            elif simulation.retrieved: simulation_status = retrieved_name

            # Get the simulation status from the remote log file if not yet retrieved
            else: simulation_status = unknown_name

            # Add the simulation properties to the list
            entries.append((path, simulation_status))

    # Return the list of simulation properties
    return entries

# -----------------------------------------------------------------

def get_retrieved_simulations(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    # Initialize a list to contain the simulations that have been retrieved
    simulations = []

    # Loop over the different entries of the status list
    for path, simulation_status in get_status_simulations(host_id):

        # Skip already retrieved and analysed simulations
        #if simulation_status == "analysed": continue
        if simulation_status != retrieved_name: continue

        # If a simulation has been retrieved earlier, but is not yet analysed, also add it to the list of retrieved
        # simulations (again) so that its results can be analysed
        #elif simulation_status == "retrieved":

        # Open the simulation file
        simulation = RemoteSimulation.from_file(path)

        # Add the retrieved simulation to the list
        simulations.append(simulation)

        # Finished simulations
        #elif simulation_status == "finished":

    # Return the list of retrieved simulations
    return simulations

# -----------------------------------------------------------------

def get_status_tasks(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    from ..basics.task import Task

    # Initialize a list to contain the statuses
    entries = []

    # Determine directory path
    local_pts_host_run_dir = fs.join(introspection.pts_run_dir, host_id)

    # Search for task files in the local PTS run/host_id directory
    for path in fs.files_in_path(local_pts_host_run_dir, extension="task", sort=int):

        # Open the task file
        task = Task.from_file(path)

        # Check whether the task has already been analysed
        if task.analysed: task_status = analysed_name

        # Check whether the task has already been retrieved
        elif task.retrieved: task_status = retrieved_name

        # Unknown
        else: task_status = unknown_name

        # Add the task properties to the list
        entries.append((path, task_status))

    # Return the list of task properties
    return entries

# -----------------------------------------------------------------

def get_retrieved_tasks(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    from ..basics.task import Task

    # Initialize a list to contain the tasks that have been retrieved
    tasks = []

    # Loop over the different entries of the status list
    for path, task_status in get_status_tasks(host_id):

        # Check status
        if task_status != retrieved_name: continue

        # If a task has been retrieved earlier, but is not yet analysed, also add it to the list of retrieved
        # tasks (again) so that its results can be analysed
        #if task_status == "retrieved":

        # Open the task file
        task = Task.from_file(path)

        # Add the task to the list
        tasks.append(task)

    # Return the list of retrieved tasks
    return tasks

# -----------------------------------------------------------------

class SKIRTRemote(Remote):

    """
    This class ...
    """

    def __init__(self, host_id=None, cluster_name=None):

        """
        The constructor ...
        :param host_id:
        :param cluster_name:
        :return:
        """

        # Call the constructor of the base class
        super(SKIRTRemote, self).__init__()

        # -- Attributes --

        # Variables storing paths to the remote SKIRT installation location
        self.skirt_path = None
        self.skirt_dir = None
        self.skirt_repo_dir = None
        self.skirt_run_dir = None
        self.local_skirt_host_run_dir = None

        # Initialize an empty list for the simulation queue
        self.queue = OrderedDict()

        # Initialize a dictionary for the scheduling options
        self.scheduling_options = dict()

        # If host ID is given, setup
        if host_id is not None:
            if not self.setup(host_id, cluster_name=cluster_name): log.warning("The connection could not be made. Run setup().")

    # -----------------------------------------------------------------

    @classmethod
    def from_remote(cls, remote):

        """
        This function ...
        :param remote:
        :return:
        """

        # Create
        skirt = cls()

        # Give warning
        log.warning("When creating a SKIRTRemote instance from a regular Remote instance, the original Remote instance can not be used anymore")

        # Set attributes
        skirt.ssh = remote.ssh
        skirt.host = remote.host
        skirt.vpn = remote.vpn
        skirt.connected = remote.connected
        skirt.commands = remote.commands

        # Reset attributes of original remote
        remote.ssh = None
        remote.connected = False

        # Locate SKIRT, if we are already connected
        if skirt.connected:
            success = skirt.locate_skirt()
            if not success: raise RuntimeError("Could not locate SKIRT on the remote host")
        else: log.warning("Not yet connected to the remote host; use setup() to initiate the connection")

        # Return the instance
        return skirt

    # -----------------------------------------------------------------

    def setup(self, host_id, cluster_name=None):

        """
        This function ...
        :param host_id:
        :param cluster_name:
        :return:
        """

        # Call the setup function of the base class
        success = super(SKIRTRemote, self).setup(host_id, cluster_name)
        if not success: return False

        # Locate SKIRT
        return self.locate_skirt()

    # -----------------------------------------------------------------

    def locate_skirt(self):

        """
        This function ...
        :return:
        """

        # Initialize
        success = True

        # Obtain some information about the SKIRT installation on the remote machine
        self.skirt_path = self.find_executable("skirt")

        # Check whether a SKIRT installation is found on the remote host
        if self.skirt_path is None: raise RuntimeError("SKIRT is not installed on the remote host")

        # We want absolute paths
        if self.skirt_path.startswith("~"):

            # Change the SKIRT path to the full, absolute path
            self.skirt_path = fs.join(self.home_directory, self.skirt_path[2:])

        # Determine the path to the SKIRT directory
        self.skirt_dir = self.skirt_path.split("/release")[0]

        # Determine the path to the SKIRT repository directory ('git')
        self.skirt_repo_dir = fs.join(self.skirt_dir, "git")

        # Determine the path to the SKIRT run directory
        self.skirt_run_dir = fs.join(self.skirt_dir, "run")

        # Determine the path to the local SKIRT run directory
        self.local_skirt_host_run_dir = fs.join(introspection.skirt_run_dir, self.host.id)

        # Create the local SKIRT run directory for this host if it doesn't already exist
        if not fs.is_directory(self.local_skirt_host_run_dir): fs.create_directory(self.local_skirt_host_run_dir, recursive=True)

        # Give a warning if the remote SKIRT version is different from the local SKIRT version
        local_version = introspection.skirt_version().split("built on")[0]
        remote_version = self.skirt_version.split("built on")[0]
        if remote_version != local_version:
            log.warning("Remote SKIRT version (" + remote_version + ") is different from local SKIRT version (" + local_version + ")")

        # Return succes
        return success

    # -----------------------------------------------------------------

    @property
    def use_hyperthreading_skirt(self):

        """
        This function ...
        :return:
        """

        return False # don't use hyperthreading for SKIRT at the moment as we assume it doesn't give much performance gain

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.queue.keys()

    # -----------------------------------------------------------------

    @property
    def nsimulations(self):

        """
        This function ...
        :return:
        """

        return len(self.simulation_names)

    # -----------------------------------------------------------------

    def add_to_queue(self, definition, logging_options, parallelization, name=None, scheduling_options=None,
                     remote_input_path=None, analysis_options=None, emulate=False, has_remote_input=False,
                     simulation_id=None, save_simulation=True, clear_existing=False, dry=False, simulation=None):

        """
        This function ...
        :param definition:
        :param logging_options:
        :param parallelization:
        :param name: a name given to the simulation
        :param scheduling_options:
        :param remote_input_path:
        :param analysis_options:
        :param emulate:
        :param has_remote_input:
        :param simulation_id:
        :param save_simulation:
        :param clear_existing:
        :param dry:
        :param simulation:
        :return:
        """

        # Inform the user
        log.info("Adding simulation to the queue ...")

        # If a name is given for the simulation, check whether it doesn't contain spaces
        if name is not None and " " in name: raise ValueError("The simulation name cannot contain spaces")

        # Get local ski path, local input and output path
        local_ski_path = definition.ski_path
        local_input_path = definition.input_path
        local_output_path = definition.output_path

        # Create the remote simulation directory
        remote_simulation_path = self.create_simulation_directory(definition, clear_existing=clear_existing, dry=dry, remote_input_path=remote_input_path)

        # Set the name if none is given
        if name is None: name = fs.name(remote_simulation_path) # = remote_simulation_name

        # Make preparations for this simulation, create the SkirtArguments object
        arguments = self.prepare(definition, logging_options, parallelization, remote_simulation_path, remote_input_path, emulate=emulate, has_remote_input=has_remote_input, dry=dry)

        # Add the SkirtArguments object to the queue
        if name in self.simulation_names: raise ValueError("Simulation name '" + name + "' is not unique")
        self.queue[name] = (arguments, local_ski_path, local_input_path, local_output_path)

        # If scheduling options are defined, add them to the dictionary
        if scheduling_options is not None: self.scheduling_options[name] = scheduling_options

        # Generate a new simulation ID based on the ID's currently in use
        if simulation_id is None: simulation_id = self._new_simulation_id()
        # Check that the simulation ID is an integer
        if not types.is_integer_type(simulation_id): raise ValueError("Simulation ID should be an integer")

        # Create a simulation object (or adapt one)
        if simulation is None: simulation = self.create_simulation_object(arguments, name, simulation_id, remote_simulation_path, definition.ski_path, local_input_path, local_output_path)
        else: self.adapt_simulation_object(simulation, arguments, name, simulation_id, remote_simulation_path, definition.ski_path, local_input_path, local_output_path)

        # Set the parallelization properties to the simulation object
        processes = arguments.parallel.processes
        threads = arguments.parallel.threads
        threads_per_core = self.threads_per_core if self.use_hyperthreading_skirt else 1
        simulation.parallelization = Parallelization.from_processes_and_threads(processes, threads, threads_per_core)

        # If analysis options are defined, adjust the correponding settings of the simulation object
        if analysis_options is not None: simulation.set_analysis_options(analysis_options)

        # Save the simulation object
        if save_simulation:
            if not dry: simulation.save()
            else: log.warning("[DRY] Not saving the simulation object ...")

        # Return the simulation object
        return simulation

    # -----------------------------------------------------------------

    def start_queue(self, queue_name=None, local_script_path=None, screen_output_path=None, schedule_method="separate",
                    group_walltime=None, jobscripts_path=None, attached=False, dry=False):

        """
        This function ...
        :param queue_name:
        :param local_script_path:
        :param screen_output_path:
        :param schedule_method:
        :param group_walltime:
        :param jobscripts_path:
        :param attached:
        :param dry:
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Inform the user
        log.info("Starting the queued simulations remotely ...")

        # Check whether attached mode is not requested for a scheduling remote
        if self.scheduler and attached: raise ValueError("Attached mode is not possible for a remote with scheduling system")

        # If the remote host uses a scheduling system, schedule all the simulations in the queue
        if self.scheduler: handles = self.start_queue_jobs(schedule_method, group_walltime=group_walltime,
                                                           jobscripts_path=jobscripts_path, queue_name=queue_name, dry=dry)

        # Else, initiate a screen session in which the simulations are executed
        else: handles = self.start_queue_screen(queue_name, local_script_path, screen_output_path, attached=attached, dry=dry)

        # Clear the queue
        self.clear_queue()

        # Return the execution handles
        return handles

    # -----------------------------------------------------------------

    def start_queue_jobs(self, method="separate", group_walltime=None, jobscripts_path=None, queue_name=None, dry=False):

        """
        This function ...
        :param method: 'separate', 'worker', 'group', 'pts'
        :param group_walltime:
        :param jobscripts_path:
        :param queue_name:
        :param dry:
        :return:
        """

        # Inform the user
        log.info("Starting the queue by scheduling the simulations in one or multiple jobs ...")

        # Determine the preferred walltime per job
        # For group methods
        preferred_walltime = group_walltime if group_walltime is not None else self.host.preferred_walltime

        # Determine queue name
        if queue_name is None: queue_name = time.unique_name("simulations")

        # Schedule simulations as separate jobs
        if method == "separate": handles = self._start_queue_jobs_separate(jobscripts_path=jobscripts_path, dry=dry)

        # Use the worker framework
        elif method == "worker": handles = self._start_queue_jobs_worker(dry=dry)

        # Group simulations in jobs
        elif method == "group": handles = self._start_queue_jobs_groups(preferred_walltime, jobscripts_path, dry=dry)

        # Dynamic method
        elif method == "pts": handles = self._start_queue_jobs_pts(preferred_walltime, queue_name, dry=dry)

        # Invalid method
        else: raise ValueError("Invalid method: '" + method + "'")

        # Return execution handles
        return handles

    # -----------------------------------------------------------------

    def _start_queue_jobs_separate(self, jobscripts_path=None, dry=False):

        """
        This function ...
        :param jobscripts_path:
        :param dry:
        :return:
        """

        # Debugging
        log.debug("Starting simulations in the queue as separate jobs ...")

        # Initialize a list to contain the execution handles
        handles = dict()

        # Loop over the items in the queue
        for name in self.simulation_names:

            # Get entry in the queue
            arguments, local_ski_path, local_input_path, local_output_path = self.queue[name]

            # Check whether scheduling options are defined for this simulation
            scheduling_options = self.scheduling_options[name] if name in self.scheduling_options else None

            # Submit the simulation to the remote scheduling system
            job_id = self.schedule(arguments, name, scheduling_options, local_ski_path=local_ski_path, jobscript_dir_path=jobscripts_path, dry=dry)

            # Set the execution handle
            if job_id is not None: handle = ExecutionHandle.job(job_id, self.host_id)
            else: handle = ExecutionHandle.postponed(self.host_id)
            handles[name] = handle

        # Return the handles
        return handles

    # -----------------------------------------------------------------

    def _start_queue_jobs_worker(self, dry=False):

        """
        This function ...
        :param dry:
        :return:
        """

        # Debugging
        log.debug("Starting simulations in the queue within the worker framework ...")

        # Loop over the items in the queue
        #for arguments, name in self.queue:

        #handle = ExecutionHandle.group_job(job_id, self.host_id)

        # Return the handle(s)
        #return handle

        # EXAMPLE:
        ## !/bin/bash -l
        ## PBS -l nodes=1:ppn=1
        ## PBS -l walltime=00:15:00
        #cd $PBS_O_WORKDIR
        #INPUT_FILE = "input_${PBS_ARRAYID}.dat"
        #OUTPUT_FILE = "output_${PBS_ARRAYID}.dat"
        #my_prog - input ${INPUT_FILE} - output ${OUTPUT_FILE}

        # SUBMISSION:
        # wsub -t 1-100 -batch test_set.pbs
        # OUTPUT:
        # total number of work items: 100
        # 123456.master15.delcatty.gent.vsc

        raise NotImplementedError("")

    # -----------------------------------------------------------------

    def _start_queue_jobs_pts(self, preferred_walltime, queue_name, dry=False):

        """
        This function ...
        :param preferred_walltime:
        :param dry:
        :return:
        """

        # Debugging
        log.debug("Starting simulations in the queue using dynamic group assignment with PTS ...")

        if dry: raise NotImplementedError("Not yet implemented")

        # Initialize a list to contain the execution handles
        #handles = dict()

        # Start a remote python session
        session = self.start_python_session(assume_pts=True)

        # Initializing the simulation queue
        session.import_package("SimulationQueue", "pts.core.remote.queue")
        session.send_line("queue = SimulationQueue('" + queue_name + "')")
        session.send_line("queue.createtable()")
        #session.send_line("queue.close()")

        # Add the simulations
        with_statement = "with queue.transaction():"
        for_statement = "for record in records:"
        line = "db.insert(label, eaglesim, 0, record['galaxyid'], skitemplate)"
        lines = [line]

        # Execute the statements
        session.with_statement_and_loop(with_statement, for_statement, lines)

        # Close the queue
        session.send_line("queue.close()")

        # For each simulation:

        # update the records to indicate that a run has been scheduled and submit the job;
        # do this within a single transaction context to ensure that the scheduled job sees the updated records
        #print("Submitting job to queue", config.queue, "for run-ids", runids)
        #with db.transaction():
        #    db.updatestatus(runids, 'scheduled')
        #    db.updatefield(runids, 'queue', config.queue)
        #    subprocess.call(("bsub",), stdin=open(jobscriptname))

        #return handles

        # The execution handle is the same for all simulations (= the simulation queue database)
        handle = ExecutionHandle.sql(queue_name, self.host_id)

        # Return the execution handle
        return handle

    # -----------------------------------------------------------------

    def _start_queue_jobs_groups(self, preferred_walltime, jobscripts_path, dry=False):

        """
        This function ...
        :param preferred_walltime:
        :param jobscripts_path:
        :param dry:
        :return:
        """

        # Debugging
        log.debug("Starting simulations in the queue using grouping into jobs ...")

        if dry: raise NotImplementedError("Not implemented")

        # Initialize a list to contain the execution handles
        handles = dict()

        # Keep track of the walltime
        current_walltime = 0.

        threads_for_job = None
        processes_for_job = None
        scheduling_options_for_job = SchedulingOptions()
        scheduling_options_for_job.walltime = preferred_walltime * 3600.  # in seconds
        simulations_for_job = []

        # Loop over the items in the queue
        #for arguments, name in self.queue:
        for name in self.simulation_names:

            # Get entry in the queue
            arguments, local_ski_path, local_input_path, local_output_path = self.queue[name]

            # Get the estimated walltime
            estimated_walltime = self.scheduling_options[name].walltime

            # If this simulation doesn't fit in the current job anymore
            if current_walltime + estimated_walltime > preferred_walltime:

                # Schedule
                job_id = self.schedule_multisim(simulations_for_job, scheduling_options_for_job, jobscripts_path)
                handles[name] = ExecutionHandle.group_job(job_id, self.host_id)

                # Reset the current walltime
                current_walltime = 0.0

                # Reset threads and processes for job
                threads_for_job = None
                processes_for_job = None

                # Reset the scheduling options
                scheduling_options_for_job = SchedulingOptions()
                scheduling_options_for_job.walltime = preferred_walltime * 3600.  # in seconds

                # Reset the list of simulations for job
                simulations_for_job = []

            # If this simulation still fits in the current job
            else:

                # Check the scheduling options and parallelization
                threads = arguments.parallel.threads
                processes = arguments.parallel.processes

                if threads_for_job is None: threads_for_job = threads
                elif threads_for_job != threads: raise ValueError("Number of threads must be equal for all simulations in job")

                if processes_for_job is None: processes_for_job = processes
                elif processes_for_job != processes: raise ValueError("Number of processes must be equal for all simulations in job")

                # Add this simulation to the list for the current job
                simulations_for_job.append(arguments)

                # Update the current walltime
                current_walltime += estimated_walltime

        # Return the execution handles
        return handles

    # -----------------------------------------------------------------

    def start_queue_screen(self, screen_name, local_script_path=None, screen_output_path=None, attached=False, dry=False):

        """
        This function ...
        :param screen_name:
        :param local_script_path:
        :param screen_output_path:
        :param attached:
        :param dry:
        :return:
        """

        # Inform the user
        log.info("Starting the queue by initiating a screen session in which all simulations are executed ...")

        # Create a unique screen name indicating we are running SKIRT simulations if none is given
        if screen_name is None: screen_name = time.unique_name("SKIRT")

        # If no screen output path is set, create a directory
        if screen_output_path is None: screen_output_path = self.create_directory_in(self.pts_temp_path, screen_name, recursive=True)

        # Create screen script
        # name, host_id, output_path
        screen = ScreenScript(screen_name, self.host_id, screen_output_path, skirt_path=self.skirt_path, mpirun_path=self.mpirun_path, remote=self)

        # Loop over the items in the queue, add a line for each simulation
        #for arguments, name in self.queue: screen.add_simulation(name, arguments)
        for name in self.simulation_names:
            arguments, local_ski_path, local_input_path, local_output_path = self.queue[name]
            screen.add_simulation(name, arguments)

        # Determine local script path
        if local_script_path is None: local_script_path = fs.join(introspection.pts_temp_dir, screen_name + ".sh")

        # Save the screen script
        screen.saveto(local_script_path)

        # Start a screen session, UNLESS DRY MODE IS ENABLED, IN WHICH CASE THE USER HAS TO UPLOAD AND RUN THE SCRIPT HIM/HERSELF
        if not dry: remote_script_path = self.start_screen(screen_name, local_script_path, self.skirt_run_dir, screen_output_path, attached=attached)
        else: remote_script_path = None

        # Create the execution handle
        if attached: handle = ExecutionHandle.tty(self.tty, self.host_id)
        else: handle = ExecutionHandle.screen(screen_name, self.host_id, remote_screen_output_path=screen_output_path, remote_screen_script_path=remote_script_path)

        # Return the execution handle(s)
        return handle

    # -----------------------------------------------------------------

    def clear_queue(self):

        """
        This function ...
        :return:
        """

        # Empty the queue and clear the scheduling options dictionary
        #self.queue = []
        self.queue = OrderedDict()
        self.scheduling_options = dict()

    # -----------------------------------------------------------------

    def run(self, definition, logging_options, parallelization, name=None, scheduling_options=None,
            analysis_options=None, local_script_path=None, screen_output_path=None, attached=False,
            show_progress=False, remote_input_path=None, has_remote_input=False, debug_output=False,
            retrieve_types=None, remove_remote_input=True, remove_remote_output=True, remove_remote_simulation_directory=True):

        """
        This function ...
        :param definition:
        :param logging_options:
        :param parallelization:
        :param name:
        :param scheduling_options:
        :param analysis_options:
        :param local_script_path:
        :param screen_output_path:
        :param attached:
        :param show_progress:
        :param remote_input_path:
        :param has_remote_input:
        :param debug_output:
        :param retrieve_types:
        :param remove_remote_input:
        :param remove_remote_output:
        :param remove_remote_simulation_directory:
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Raise an error if there are other simulations currently waiting in the queue
        if len(self.queue) > 0: raise RuntimeError("The simulation queue is not empty")

        # Add the simulation arguments to the queue
        simulation = self.add_to_queue(definition, logging_options, parallelization, name, scheduling_options, analysis_options=analysis_options, remote_input_path=remote_input_path, has_remote_input=has_remote_input)

        # Check whether attached mode is not requested for a scheduling remote
        if self.scheduler and attached: raise ValueError("Attached mode is not possible for a remote with scheduling system")

        # Progress bar: ask attached but set detached so that the start_queue function returns immediately
        if show_progress and not attached: raise ValueError("Cannot show progress when 'attached' is False")
        if show_progress: attached = False

        # Start the queue, get execution handle(s)
        handles = self.start_queue(name, local_script_path, screen_output_path, attached=attached)

        # Set remote simulation properties
        simulation.retrieve_types = retrieve_types

        # Set remove remote files settings
        simulation.remove_remote_input = remove_remote_input
        simulation.remove_remote_output = remove_remote_output
        simulation.remove_remote_simulation_directory = remove_remote_simulation_directory

        # Set the execution handle for the simulation
        simulation.handle = handles if isinstance(handles, ExecutionHandle) else handles[0]
        simulation.save()

        # Show progress bar with progress
        if show_progress:

            #out_path = arguments.output_path if arguments.output_path is not None else fs.cwd()
            #prefix = arguments.prefix
            #log_path = fs.join(out_path, prefix + "_log.txt")
            status = LogSimulationStatus(simulation.remote_log_file_path, remote=self, debug_output=debug_output)

            # Get the execution handle for the simulation
            handle = handles

            # Show the simulation progress
            with log.no_debugging(): success = status.show_progress(handle)

            # Check whether not crashed
            if not success: raise RuntimeError("The simulation crashed")

        # Return the simulation object
        return simulation

    # -----------------------------------------------------------------

    def create_simulation_directory(self, definition, clear_existing=False, dry=False, remote_input_path=None):

        """
        This function ...
        :param definition:
        :param clear_existing:
        :param dry:
        :param remote_input_path:
        :return:
        """

        # Determine simulation directory name
        if definition.name is not None: remote_simulation_name = definition.name
        else:

            # Create a unique name for the simulation directory
            skifile_name = fs.name(definition.ski_path).split(".ski")[0]
            remote_simulation_name = time.unique_name(skifile_name, separator="__")

        # If an output path is defined in the remote host configuration file, use it for the simulation directory (input and output)
        if self.host.output_path is not None:

            # Get the output directory path
            output_path = self.absolute_path(self.host.output_path)

            # Get the simulation directory path
            remote_simulation_path = fs.join(output_path, remote_simulation_name)

        # Determine the full path of the simulation directory on the remote system
        else: remote_simulation_path = fs.join(self.skirt_run_dir, remote_simulation_name)

        # Check whether the simulation directory will contain input
        if remote_input_path is not None: contains_input = self.is_subdirectory(remote_input_path, remote_simulation_path)
        else: contains_input = False
        if self.is_directory(remote_simulation_path) and not contains_input:
            if clear_existing:
                if not dry: self.remove_directory(remote_simulation_path)
                else: log.warning("[DRY] Not removing directory '" + remote_simulation_path + "' ...")
            else: raise IOError("Simulation directory already exists")

        # Create the remote simulation directory
        if not self.is_directory(remote_simulation_path): # because remove can be skipped if 'contains_input'
            if not dry: self.create_directory(remote_simulation_path, recursive=True)
            else: log.warning("[DRY] Not creating directory '" + remote_simulation_path + "' ...")

        # Return the path to the remote simulation directory
        return remote_simulation_path

    # -----------------------------------------------------------------

    def prepare(self, definition, logging_options, parallelization, remote_simulation_path, remote_input_path=None,
                emulate=False, has_remote_input=False, dry=False):

        """
        This function ...
        :param definition:
        :param logging_options:
        :param parallelization:
        :param remote_simulation_path:
        :param remote_input_path:
        :param emulate:
        :param has_remote_input:
        :param dry:
        :return:
        """

        # Create the SkirtArguments object
        arguments = SkirtArguments(logging_options=logging_options, parallelization=parallelization, emulate=emulate)

        # If an output path is defined in the remote host configuration file, use it for the simulation output
        if self.host.output_path is not None:

            # Create the output directory
            output_path = self.absolute_path(self.host.output_path)
            if not self.is_directory(output_path):
                if not dry: self.create_directory(output_path)
                else: log.warning("[DRY] Not creating directory '" + output_path + "' ...")

            # Get the name of the remote simulation directory and use use that name for the output directory
            remote_simulation_name = fs.name(remote_simulation_path)
            remote_output_path = fs.join(output_path, remote_simulation_name)

            # Expand the alias to the user's home directory
            #remote_output_path = self.absolute_path(remote_output_path)

            # If the remote output path is the same as the remote simulation path, use a folder called 'out' inside
            # the simulation directory instead for the output
            if remote_output_path == remote_simulation_path: remote_output_path = fs.join(remote_output_path, "out")

        # If an output path is not specified by the user, place a directory called 'out' next to the simulation's 'in' directory
        else: remote_output_path = fs.join(remote_simulation_path, "out")

        # Prepare the input (upload it)
        remote_input_path = self.prepare_input(definition.input_path, remote_simulation_path, remote_input_path, has_remote_input=has_remote_input, dry=dry)

        # Create the remote output directory
        if not self.is_directory(remote_output_path):
            if not dry: self.create_directory(remote_output_path)
            else: log.warning("[DRY] Not creating directory '" + remote_output_path + "' ...")

        # Set the remote ski file path
        local_ski_path = definition.ski_path
        ski_name = fs.name(local_ski_path)
        remote_ski_path = fs.join(remote_simulation_path, ski_name)

        # Set the base simulation options such as ski path, input path and output path (remote)
        arguments.ski_pattern = remote_ski_path
        arguments.input_path = remote_input_path
        arguments.output_path = remote_output_path

        # Copy the input directory and the ski file to the remote host
        if not dry: self.upload(local_ski_path, remote_simulation_path)
        else: log.warning("[DRY] Not uploading '" + local_ski_path + "' to '" + remote_simulation_path + "' ...")

        # Return the SKIRT arguments instance
        return arguments

    # -----------------------------------------------------------------

    def prepare_input(self, input_path, remote_simulation_path, remote_input_path=None, has_remote_input=False, dry=False):

        """
        This function ...
        :param input_path:
        :param remote_simulation_path:
        :param remote_input_path:
        :param has_remote_input:
        :param dry:
        :return:
        """

        # The simulation does not require input
        if input_path is None: remote_input_path = None

        # The simulation input is defined in terms of a single local directory
        elif isinstance(input_path, basestring):

            # Check has_remote_input flag
            if has_remote_input: raise ValueError("Cannot enable 'has_remote_input' flag when input is a directory, only when input is defined in terms of a list or dictionary of filepaths (or SimulationInput object)")

            # A remote input path is not specified, this means that we have yet to copy the input
            if remote_input_path is None:

                # Determine the full path to the input directory on the remote system
                remote_input_path = fs.join(remote_simulation_path, "in")

                # Copy the input directory to the remote host
                if not dry: self.upload(input_path, remote_input_path, show_output=True)
                else: log.warning("[DRY] Not uploading '" + input_path + "' to '" + remote_input_path + "' ...")

            # The specified remote input directory (for re-usage of already uploaded input) does not exist
            else:
                if has_remote_input: raise ValueError("Cannot specify 'has_remote_input' and 'remote_input_path' simultaneously")
                if not self.is_directory(remote_input_path): raise RuntimeError("The remote input directory does not exist: '" + remote_input_path + "'")

        # If the simulation input is defined as a list of seperate file paths
        elif isinstance(input_path, list):

            local_input_file_paths = input_path # the list of file paths

            # A remote input path is not specified, this means that we have yet to copy the input files to a new
            # remote directory
            if remote_input_path is None:

                # Determine the full path to the remote input directory on the remote system
                remote_input_path = fs.join(remote_simulation_path, "in")

                # Create the remote directory
                if not dry: self.create_directory(remote_input_path)
                else: log.warning("[DRY] Not creating directory '" + remote_input_path + "' ...")

                # Check which files are actually local and which are already remote
                if has_remote_input:
                    nfiles = len(local_input_file_paths)
                    remote_input_file_paths = [filepath for filepath in local_input_file_paths if self.is_file(filepath)]
                    local_input_file_paths = [filepath for filepath in local_input_file_paths if fs.is_file(filepath)]
                    if len(remote_input_file_paths) + len(local_input_file_paths) != nfiles: raise ValueError("Some input files were not found")
                else: remote_input_file_paths = None

                # Upload the local input files to the new remote directory
                if not dry: self.upload(local_input_file_paths, remote_input_path)
                else: log.warning("[DRY] Not uploading '" + tostr(local_input_file_paths) + "' to '" + remote_input_path + "' ...")

                # Copy the already remote files to the new remote directory
                if remote_input_file_paths is not None:
                    if dry: log.warning("[DRY] Not copying files '" + tostr(remote_input_file_paths) + "' to '" + remote_input_path + "' ...")
                    else: self.copy_files(remote_input_file_paths, remote_input_path)

            # The specified remote directory (for re-usage of already uploaded input) does not exist
            else:
                if has_remote_input: raise ValueError("Cannot specify 'has_remote_input' and 'remote_input_path' simultaneously")
                if not self.is_directory(remote_input_path): raise RuntimeError("The remote input directory does not exist: '" + remote_input_path + "'")

        # If we have a SimulationInput instance
        elif isinstance(input_path, SimulationInput):

            # Check
            if has_remote_input: raise ValueError("Currently, has_remote_input is not possible with SimulationInput objects (due to internal checking in the latter class)")

            # If a remote input path is not specified
            if remote_input_path is None:

                # Determine the full path to the remote input directory on the remote system
                remote_input_path = fs.join(remote_simulation_path, "in")

                # Create the remote directory
                if not dry: self.create_directory(remote_input_path)
                else: log.warning("[DRY] Not creating directory '" + remote_input_path + "' ...")

                # Upload the local input files to the new remote directory
                for name, path in input_path:

                    # Debugging
                    log.debug("Uploading the '" + path + "' file to '" + remote_input_path + "' under the name '" + name + "'")

                    # Upload the file, giving it the desired name
                    if not dry: self.upload(path, remote_input_path, new_name=name)
                    else: log.warning("[DRY] Not uploading '" + path + "' to '" + remote_input_path)

            # The specified remote directory (for re-usage of already uploaded input) does not exist
            else:
                if has_remote_input: raise ValueError("Cannot specify 'has_remote_input' and 'remote_input_path' simultaneously")
                if not self.is_directory(remote_input_path): raise RuntimeError("The remote input directory does not exist: '" + remote_input_path + "'")

        # We have a dictionary
        elif types.is_dictionary(input_path):

            # If a remote input path is not specified
            if remote_input_path is None:

                # Determine the full path to the remote input directory on the remote system
                remote_input_path = fs.join(remote_simulation_path, "in")

                # Create the remote directory
                if not dry: self.create_directory(remote_input_path)
                else: log.warning("[DRY] Not creating directory '" + remote_input_path + "' ...")

                # Upload the local input files to the new remote directory
                for name, path in input_path.items():

                    # Check
                    if has_remote_input:

                        # Is local
                        if fs.is_file(path):

                            # Debugging
                            log.debug("Uploading the '" + path + "' file to '" + remote_input_path + "' under the name '" + name + "'")

                            # Upload
                            if not dry: self.upload(path, remote_input_path, new_name=name)
                            else: log.warning("[DRY] Not uploading '" + path + "' to '" + remote_input_path + "' ...")

                        # Is remote
                        elif self.is_file(path):

                            # Debugging
                            log.debug("Copying the '" + path + "' file from '" + fs.directory_of(path) + "' to '" + remote_input_path + "' under the name '" + name + "'")

                            # Copy
                            if not dry: self.copy_file(path, remote_input_path, new_name=name)
                            else: log.warning("[DRY] Not copying file '" + path + "' to '" + remote_input_path + "' ...")

                        # Not found
                        else: raise ValueError("The input file '" + name + "' is not found remotely or locally")

                    # No checking
                    else:

                        # Debugging
                        log.debug("Uploading the '" + path + "' file to '" + remote_input_path + "' under the name '" + name + "'")

                        # Upload the file, giving it the desired name
                        if not dry: self.upload(path, remote_input_path, new_name=name)
                        else: log.warning("[DRY] Not uploading '" + path + "' to '" + remote_input_path + "' ...")

            # The specified remote directory (for re-usage of already uploaded input) does not exist
            else:
                if has_remote_input: raise ValueError("Cannot specify 'has_remote_input' and 'remote_input_path' simultaneously")
                if not self.is_directory(remote_input_path): raise RuntimeError("The remote input directory does not exist: '" + remote_input_path + "'")

        # Invalid format for arguments.input_path
        else: raise ValueError("Invalid value for 'input_path': must be None, local directory path, list of file paths or SimulationInput object")

        # Return the remote input path
        return remote_input_path

    # -----------------------------------------------------------------

    def create_simulation_object(self, arguments, name, simulation_id, remote_simulation_path, local_ski_path, local_input_path, local_output_path):

        """
        This function ...
        :param arguments:
        :param name:
        :param simulation_id:
        :param remote_simulation_path:
        :param local_ski_path:
        :param local_input_path:
        :param local_output_path:
        :return:
        """

        # Create a new remote simulation object
        simulation = RemoteSimulation(local_ski_path, local_input_path, local_output_path)

        # Determine and set the simulation file path
        simulation_file_path = fs.join(self.local_skirt_host_run_dir, str(simulation_id) + ".sim")

        # Set the host ID and cluster name (if applicable)
        simulation.host_id = self.host_id
        simulation.cluster_name = self.cluster_name

        # Set other attributes
        simulation.path = simulation_file_path
        simulation.id = simulation_id
        simulation.name = name
        simulation.remote_ski_path = arguments.ski_pattern
        simulation.remote_simulation_path = remote_simulation_path
        simulation.remote_input_path = arguments.input_path
        simulation.remote_output_path = arguments.output_path
        simulation.submitted_at = time.timestamp()

        # Return the simulation object
        return simulation

    # -----------------------------------------------------------------

    def adapt_simulation_object(self, simulation, arguments, name, simulation_id, remote_simulation_path, local_ski_path, local_input_path, local_output_path):

        """
        This function ...
        :param simulation:
        :param arguments:
        :param name:
        :param simulation_id:
        :param remote_simulation_path:
        :param local_ski_path:
        :param local_input_path:
        :param local_output_path:
        :return:
        """

        # Set basic properties
        simulation.ski_path = local_ski_path
        simulation.input_path = local_input_path
        simulation.output_path = local_output_path

        # Set the host ID and cluster name (if applicable)
        simulation.host_id = self.host_id
        simulation.cluster_name = self.cluster_name

        # Set other attributes
        simulation.id = simulation_id
        simulation.name = name
        simulation.remote_ski_path = arguments.ski_pattern
        simulation.remote_simulation_path = remote_simulation_path
        simulation.remote_input_path = arguments.input_path
        simulation.remote_output_path = arguments.output_path
        simulation.submitted_at = time.timestamp()

    # -----------------------------------------------------------------

    def schedule(self, arguments, name, scheduling_options, local_ski_path, jobscript_dir_path=None, dry=False, save_jobscript=True):

        """
        This function ...
        :param arguments:
        :param name:
        :param scheduling_options:
        :param local_ski_path:
        :param jobscript_dir_path:
        :param dry:
        :param save_jobscript:
        :return:
        """

        # Inform the suer
        log.info("Scheduling simulation '" + name + "' on the remote host ...")

        # Jobscript directory path is given
        if jobscript_dir_path is None:

            # Determine the jobscript path
            local_simulation_path = fs.directory_of(local_ski_path)
            jobscript_dir_path = local_simulation_path

        # Verify the scheduling options
        scheduling_options = self._verify_scheduling_options(scheduling_options, arguments, jobscript_dir_path)

        # Now get the options
        nodes = scheduling_options.nodes
        ppn = scheduling_options.ppn
        mail = scheduling_options.mail
        full_node = scheduling_options.full_node
        walltime = scheduling_options.walltime
        local_jobscript_path = scheduling_options.local_jobscript_path

        modules = []
        # module spider mympirun
        #modules.append("vsc-mympirun/3.4.3-intel-2016b-Python-2.7.12")

        local_skirt_version = introspection.skirt_main_version()
        if local_skirt_version == 7:
            modules.append("iimpi/2016b") # this loads GCC, icc, impi, python, iimpi, vsc-base, vsc-mympirun etc.
            modules.append("vsc-mympirun") # this loads mympirun
        elif local_skirt_version == 8:
            modules.append("intel")
            modules.append("CMake/3.10.1-GCCcore-6.4.0")
            modules.append("vsc-mympirun")
        else: raise ValueError("Invalid version")

        # Create a job script next to the (local) simulation's ski file
        jobscript_name = fs.name(local_jobscript_path)
        #jobscript = JobScript(local_jobscript_path, arguments, self.host.cluster, self.skirt_path,
        #                      self.host.mpi_command, modules, walltime, nodes, ppn, name=name, mail=mail,
        #                      full_node=full_node, bind_to_cores=self.host.force_process_binding)

        # Determine remote jobscript location
        remote_simulation_path = fs.directory_of(arguments.ski_pattern)  # NEW, to avoid having to pass this as an argument
        remote_jobscript_path = fs.join(remote_simulation_path, jobscript_name)

        # Set header lines
        header_lines = []
        header_lines.append("To submit manualy, upload this file to the remote filesystem in the directory:")
        header_lines.append("'" + remote_simulation_path + "'")
        header_lines.append("with the name '" + jobscript_name + "' and enter the following commmands:")
        header_lines.append("cd '" + remote_simulation_path + "' # navigate to the simulation directory")
        header_lines.append("module swap cluster/[cluster_name] # swap to cluster [cluster_name], if desired")
        header_lines.append("qsub " + jobscript_name)
        header_lines.append("when the job has been submitted, copy the job ID that is returned [XXXX] and run the following command (locally):")
        header_lines.append("pts set_postponed_job_id " + self.host_id + " " + name + " XXXX")
        header_lines.append("")

        # Create the job script
        jobscript = SKIRTJobScript(name, arguments, self.host.id, self.host.cluster, self.skirt_path, self.host.mpi_command, walltime,
                                   modules, mail=mail, bind_to_cores=self.host.force_process_binding,
                                   extra_header_lines=header_lines, remote=self) #, hyperthreading=self.use_hyperthreading_skirt)

        # Save the job script locally
        if save_jobscript: jobscript.saveto(local_jobscript_path)
        else: jobscript.path = local_jobscript_path

        # Upload and submit
        job_id = None
        if dry: log.warning("[DRY] Dry mode is enabled, run job '" + jobscript_name + "' by locating the job script file at '" + local_jobscript_path + "' and following the instructions therein")
        else:
            job_id = self.upload_and_submit_job(local_jobscript_path, remote_jobscript_path, remote_simulation_path)
            log.debug("The job ID for the simulation is '" + str(job_id) + "'")

        # Return the job ID
        return job_id

    # -----------------------------------------------------------------

    def upload_and_submit_job(self, local_jobscript_path, remote_jobscript_path, remote_simulation_path):

        """
        This function ...
        :return: 
        """

        # Debugging
        log.debug("Uploading and submitting job from script '" + remote_jobscript_path + "'")

        # Copy the job script to the remote simulation directory
        self.upload(local_jobscript_path, remote_simulation_path)

        ## Swap clusters
        # Then, swap to the desired cluster and launch the job script
        # output = subprocess.check_output("module swap cluster/" + self._clustername + "; qsub " + self._path, shell=True, stderr=subprocess.STDOUT)

        # Submit the job script to the remote scheduling system
        # output = self.execute("qsub " + remote_jobscript_path, contains_extra_eof=True)
        output = self.execute("qsub " + remote_jobscript_path)

        # The queue number of the submitted job is used to identify this simulation
        job_id = int(output[0].split(".")[0])
        return job_id

    # -----------------------------------------------------------------

    def upload_and_submit_worker_job(self, local_jobscript_path, local_data_path, remote_path):

        """
        This function ...
        :param local_jobscript_path:
        :param local_data_path:
        :return:
        """

        # Debugging
        log.debug("Uploading and submitting ")

        # Upload

        # Execute the command
        command = "wsub -batch " + remote_jobscript_path + " -data " + remote_data_path
        output = self.execute(command)

        # Get the number of work items
        nitems = int(output[0].split("work items: ")[1])

        # Get the queue number
        job_id = int(output[1].split(".")[0])
        return job_id

    # -----------------------------------------------------------------

    def schedule_multisim(self, arguments, scheduling_options, jobscripts_path):

        """
        This function ...
        :param arguments:
        :param scheduling_options:
        :return:
        """

        # Inform the suer
        log.info("Scheduling a job of " + str(len(arguments)) +  " simulations on the remote host ...")

        arguments_of_first_simulation = arguments[0] # for verifying the scheduling options (ALL SIMULATION HAVE THE SAME PARALLELIZATION MODE IN THIS MULTISIM FUNCTION)

        # Verify the scheduling options
        scheduling_options = self._verify_scheduling_options(scheduling_options, arguments_of_first_simulation, jobscripts_path)

        # Now get the options
        nodes = scheduling_options.nodes
        ppn = scheduling_options.ppn
        mail = scheduling_options.mail
        full_node = scheduling_options.full_node
        walltime = scheduling_options.walltime
        local_jobscript_path = scheduling_options.local_jobscript_path

        # Jobscript name
        jobscript_name = fs.name(local_jobscript_path)
        jobscript = MultiJobScript(local_jobscript_path, )

        # Copy the job script to the remote simulation directory
        remote_simulation_path = fs.directory_of(arguments.ski_pattern)  # NEW, to avoid having to pass this as an argument
        remote_jobscript_path = fs.join(remote_simulation_path, jobscript_name)
        self.upload(local_jobscript_path, remote_simulation_path)

        # SUBMIT THE JOB
        output = self.execute("qsub " + remote_jobscript_path)

        # The queue number of the submitted job is used to identify this simulation
        job_id = int(output[0].split(".")[0])

        # Return the job ID
        return job_id

    # -----------------------------------------------------------------

    def start(self, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Inform the user
        log.info("Starting simulation on the remote host ...")

        # Send the command to the remote machine using a screen session so that we can safely detach from the remote shell
        command = arguments.to_command(self.scheduler, skirt_path=self.skirt_path, mpirun_path=self.host.mpi_command,
                                       bind_to_cores=self.host.force_process_binding, to_string=True, remote=self)
        self.execute("screen -d -m " + command, output=False)

        # Generate a new simulation ID based on the ID's currently in use
        simulation_id = self._new_simulation_id()

        # Return the simulation ID
        return simulation_id

    # -----------------------------------------------------------------

    def _simulation_ids_in_use(self):

        """
        This function ...
        :return:
        """

        # Check the contents of the local run directory to see which simulation id's are currently in use
        current_ids = []
        for name in fs.files_in_path(self.local_skirt_host_run_dir, extension="sim", returns="name"):

            # Get the simulation ID and add it to the list
            current_ids.append(int(name))

        # Return the list of currently used ID's
        return current_ids

    # -----------------------------------------------------------------

    def _new_simulation_id(self):

        """
        This function ...
        :return:
        """

        # Get a list of the ID's currently in use
        current_ids = self._simulation_ids_in_use()

        # Return the lowest missing integer
        return numbers.lowest_missing_integer(current_ids)

    # -----------------------------------------------------------------

    def _new_simulation_ids(self, count):

        """
        This function ...
        :param count:
        :return:
        """

        # Get a list of the ID's currently in use
        current_ids = self._simulation_ids_in_use()

        # Initialize a list to contain the new ID's
        new_ids = []

        # Sort the current simulation ID's and find the lowest 'missing' integer number
        if len(current_ids) > 0:
            current_ids = sorted(current_ids)
            for index in range(max(current_ids)):
                if current_ids[index] != index:
                    new_ids.append(index)
                    if len(new_ids) == count: return new_ids

            # Complement with new ID's
            max_id = max(new_ids)
            missing = count - len(new_ids)

            for index in range(max_id+1, max_id+1+missing):

                new_ids.append(index)

            return new_ids

        # If no simulation ID's are currently in use, return a list of the integers from 0 to count-1
        else: return range(count)

    # -----------------------------------------------------------------

    def retrieve(self, retrieve_crashed=None, check_crashed=False, check_data=False):

        """
        This function ...
        :param retrieve_crashed:
        :param check_crashed:
        :param check_data:
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Initialize a list to contain the simulations that have been retrieved
        simulations = []

        # Loop over the different entries of the status list
        for path, simulation_status in self.get_status():

            # Skip already fully analysed simulations
            if simulation_status == analysed_name: continue

            # If a simulation is only partly analysed, add it to the list of retrieved simulations (again) so that
            # the analysis can be completed
            elif simulation_status.startswith(analysed_name):

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                # Add the retrieved simulation to the list
                simulations.append(simulation)

            # If a simulation has been retrieved earlier, but is not yet analysed, also add it to the list of retrieved
            # simulations (again) so that its results can be analysed
            elif simulation_status == retrieved_name:

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                # Add the retrieved simulation to the list
                simulations.append(simulation)

            # Finished simulations
            elif simulation_status == finished_name:

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                # Retrieve the simulation
                self.retrieve_simulation(simulation)

                # Remove the simulation from the remote
                simulation.remove_from_remote(self)

                # Add the simulation to the list of retrieved simulations
                simulations.append(simulation)

            # Crashed simulations
            elif crashed_name in simulation_status:

                # Do we need to bother looking?
                if retrieve_crashed is None and not check_crashed: continue

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                # Check the simulation ID
                if retrieve_crashed and simulation.id not in retrieve_crashed:

                    # If check crashed for their output
                    if check_crashed and "writing" in simulation_status: # already in the writing phase

                        # Give message to the user
                        log.success("Simulation crashed but already in the writing phase. Trying to retrieve the present output and continuing from there")

                    # Don't check crashed simulations
                    else: continue
                
                # Retrieve the simulation
                self.retrieve_simulation(simulation)

                # Don't remove from the remote
                log.warning("Because the simulation '" + simulation.name + "' crashed, the remote simulation directories will be kept for manual inspection")

                # Check the data of the simulation (whether it is valid)
                if check_data:

                    # Load the simulation data
                    data = SimulationData.from_directory(simulation.output_path)

                    # Check whether valid
                    if not data.valid:

                        # Give error message
                        log.error("The simulation data has invalid files. Not counting this simulation as succesfully retrieved: analysis will not be performed")
                        log.error("Summary of the retrieved data at '" + simulation.output_path + "':")

                        # Show a summary of the simulation data
                        data.show()

                        # Don't add the simulation to the list
                        continue

                # Add the simulation to the list of retrieved simulations
                simulations.append(simulation)

        # Return the list of retrieved simulations
        return simulations

    # -----------------------------------------------------------------

    def retrieve_simulation(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Debug info
        log.debug("Retrieving simulation " + str(simulation.name) + " with id " + str(simulation.id) + " ...")

        # If retrieve file types are not defined, download the complete output directory
        if simulation.retrieve_types is None or simulation.retrieve_types == "None": self.retrieve_simulation_all(simulation)

        # If retrieve file types are defined, download these files seperately to the local filesystem
        else: self.retrieve_simulation_types(simulation)

        # Check if there are any output files in the local output directory
        if fs.is_empty(simulation.output_path): raise ValueError("Retrieval was unsuccesful: no output files in '" + simulation.output_path + "'")

        # If retrieval was succesful, add this information to the simulation file
        simulation.retrieved = True
        simulation.save()

        # Debug info
        log.debug("Successfully retrieved the necessary simulation output")

    # -----------------------------------------------------------------

    def retrieve_simulation_all(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Debug info
        log.debug("Retrieve file types are not defined, retrieving complete remote output directory ...")
        log.debug("Local output directory: " + simulation.output_path)

        from ...magic.core.fits import is_valid

        # Check whether the output directory exists; if not, create it
        if not fs.is_directory(simulation.output_path): fs.create_directory(simulation.output_path)

        # Check whether the local output directory is not empty
        if not fs.is_empty(simulation.output_path):

            # Same number of files?
            #if fs.nfiles_in_path(simulation.output_path) == self.nfiles_in_path(simulation.remote_output_path):

            # Check the exact filenames
            copy_paths = []
            for remote_filepath, filename in self.files_in_path(simulation.remote_output_path, returns=["path", "name"], extensions=True):

                # Determine local filepath
                local_filepath = fs.join(simulation.output_path, filename)

                # Is the file present locally?
                if fs.is_file(local_filepath):

                    # Does the filesize match?
                    if fs.file_nbytes(local_filepath) != self.file_nbytes(remote_filepath):

                        log.warning("Output file '" + filename + "' was incompletely retrieved: downloading again ...")
                        fs.remove_file(local_filepath)
                        copy_paths.append(remote_filepath)

                    # For a FITS file, ALSO check whether it's valid
                    elif local_filepath.endswith(".fits") and not is_valid(local_filepath):

                        log.warning("Output file '" + filename + "' was incompletely retrieved: downloading again ...")
                        fs.remove_file(local_filepath)
                        copy_paths.append(remote_filepath)

                else: copy_paths.append(remote_filepath)

            # If the files have all been download already
            if len(copy_paths) == 0: log.success("All output files of simulation '" + simulation.name + "' have already been retrieved succesfully")
            else:

                # Debugging
                log.debug("Retrieving files: " + ", ".join(copy_paths))

                # Download
                self.download(copy_paths, simulation.output_path)

        # Download the simulation output
        else: self.download(simulation.remote_output_path, simulation.output_path)

    # -----------------------------------------------------------------

    def retrieve_simulation_types(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Debug info
        log.debug("Retrieve file types are defined, retrieving specific output files ...")
        log.debug("Local output directory: " + simulation.output_path)

        from ...magic.core.fits import is_valid

        # Create a list for the paths of the files that have to be copied to the local filesystem
        copy_paths = []

        # Loop over the files that are present in the remoute output directory
        for filename in self.files_in_path(simulation.remote_output_path, extensions=True):

            # Determine the full path to the output file
            filepath = fs.join(simulation.remote_output_path, filename)

            # Check whether the file has to be retrieved
            if needs_retrieval(filename, simulation.retrieve_types):

                local_filepath = fs.join(simulation.output_path, filename)

                if not fs.is_file(local_filepath): copy_paths.append(filepath)

                else:

                    if fs.file_nbytes(local_filepath) != self.file_nbytes(filepath):

                        log.warning("Output file '" + filename + "' was incompletely retrieved: downloading again ...")
                        fs.remove_file(local_filepath)
                        copy_paths.append(filepath)

                    elif local_filepath.endswith(".fits") and not is_valid(local_filepath):

                        log.warning("Output file '" + filename + "' was incompletely retrieved: downloading again ...")
                        fs.remove_file(local_filepath)
                        copy_paths.append(filepath)

                    else: log.warning("The output file '" + filename + "' has already been (completely) retrieved: skipping ...")

        # All already retrieved
        if len(copy_paths) == 0: log.success("All necessary output files have already been retrieved succesfully")
        else:

            # Debugging
            log.debug("Retrieving files: " + ", ".join(copy_paths))
            log.debug("Local output directory: " + simulation.output_path)

            # Check whether the output directory exists; if not, create it
            if not fs.is_directory(simulation.output_path): fs.create_directory(simulation.output_path)

            # Download the list of files to the local output directory
            self.download(copy_paths, simulation.output_path)

    # -----------------------------------------------------------------

    def _get_simulation_status_not_scheduler(self, simulation, screen_states=None):

        """
        This function ...
        :param simulation:
        :param screen_states:
        :return:
        """

        # The name of the ski file (the simulation prefix)
        ski_name = simulation.prefix()

        # Determine the path to the remote log file
        remote_log_file_path = simulation.remote_log_file_path

        # Check whether the simulation has already been analysed
        if simulation.analysed: simulation_status = analysed_name

        # Partly analysed
        elif simulation.analysed_any:

            analysed = []
            if simulation.analysed_all_extraction: analysed.append("extraction")
            if simulation.analysed_all_plotting: analysed.append("plotting")
            if simulation.analysed_all_misc: analysed.append("misc")
            if simulation.analysed_batch: analysed.append("batch")
            if simulation.analysed_scaling: analysed.append("scaling")
            if simulation.analysed_all_extra: analysed.append("extra")

            if len(analysed) > 0: simulation_status = "analysed: " + ", ".join(analysed)
            else: simulation_status = "analysed: started"

        # Check whether the simulation has already been retrieved
        elif simulation.retrieved: simulation_status = retrieved_name

        # Check whether the simulation has finished
        elif simulation.finished: simulation_status = "finished"

        # Get the simulation status from the remote log file if not yet retrieved
        else: simulation_status = self.status_from_log_file(remote_log_file_path, simulation.handle, ski_name, screen_states=screen_states, remote_ski_path=simulation.remote_ski_path)

        # Return the simulation status
        return simulation_status

    # -----------------------------------------------------------------

    def _get_simulation_status_scheduler(self, simulation, jobs_status, python_session):

        """
        This function ...
        :param simulation:
        :param python_session:
        :return:
        """

        # The name of the ski file (the simulation prefix)
        ski_name = simulation.prefix()

        # The path to the simulation log file
        remote_log_file_path = simulation.remote_log_file_path

        # Check if the simulation has already been analysed
        if simulation.analysed: simulation_status = analysed_name

        # Partly analysed
        elif simulation.analysed_any:

            analysed = []
            if simulation.analysed_all_extraction: analysed.append("extraction")
            if simulation.analysed_all_plotting: analysed.append("plotting")
            if simulation.analysed_all_misc: analysed.append("misc")
            if simulation.analysed_batch: analysed.append("batch")
            if simulation.analysed_scaling: analysed.append("scaling")
            if simulation.analysed_all_extra: analysed.append("extra")

            if len(analysed) > 0: simulation_status = "analysed: " + ", ".join(analysed)
            else: simulation_status = "analysed: started"

        # Check if the simulation has already been retrieved
        elif simulation.retrieved: simulation_status = retrieved_name

        # Simulation is finished
        elif simulation.finished: simulation_status = "finished"

        # No handle?
        elif simulation.handle is None:
            log.warning("Simulation '" + simulation.name + "' has no execution handle")
            simulation_status = "unknown"

        # Simulation is a job
        elif simulation.handle.type == "job":

            # Get the job ID
            job_id = simulation.handle.value

            # Check if the job ID is in the list of queued or running jobs
            if job_id in jobs_status:

                # Check the status of this simulation
                job_status = jobs_status[job_id]

                # This simulation is still queued
                if job_status == 'Q': simulation_status = queued_name

                # This simulation is currently running
                elif job_status == 'R': simulation_status = self.running_status_from_log_file(remote_log_file_path)

                # If the job has been cancelled, check whether some part of the log file was already present
                # (the simulation was running but was aborted) or the log file is not present (the simulation is cancelled)
                elif job_status == "C":

                    if self.is_file(remote_log_file_path): simulation_status = aborted_name
                    else: simulation_status = cancelled_name

                # This simulation has an unknown status, check the log file
                else: simulation_status = self.status_from_log_file_job(remote_log_file_path, ski_name)

            # Job not present in the queue anymore: finished, crashed or aborted
            else: simulation_status = self.status_from_log_file_job(remote_log_file_path, ski_name)

        # Simulation is part of a group of simulations in a job
        elif simulation.handle.type == "group-job":

            # Get the job ID
            job_id = simulation.handle.value

            # Check if the job ID is in the list of queued or running jobs
            if job_id in jobs_status:

                # Check the status of the job
                job_status = jobs_status[job_id]

                # If the job is still queued, the simulation is also still queued
                if job_status == 'Q': simulation_status = queued_name

                # If the job is running
                elif job_status == 'R':

                    # Check if the log file exists
                    if self.is_file(remote_log_file_path):

                        # Get the last two lines of the remote log file
                        output = self.read_last_lines(remote_log_file_path, 2)

                        # Get the last line of the actual simulation
                        if len(output) == 0: return "invalid: cannot read log file" #simulation_status = "invalid: cannot read log file"
                        elif len(output) == 1: last = output[0]
                        elif " Available memory: " in output[1]: last = output[0]
                        else: last = output[1]

                        # Interpret the content of the last line
                        if " Finished simulation " + ski_name in last: simulation_status = finished_name
                        elif " *** Error: " in last: simulation_status = self.crashed_status_from_log_file(remote_log_file_path)
                        else: simulation_status = self.running_status_from_log_file(remote_log_file_path)

                    # The job is running but this simulation does not have a log file yet
                    else: simulation_status = queued_name

                # If the job has been cancelled
                elif job_status == 'C':

                    # Check if the log file exists
                    if self.is_file(remote_log_file_path): simulation_status = self.status_from_log_file_job(remote_log_file_path, ski_name)
                    else: simulation_status = cancelled_name

                # This simulation has an unknown status, check the log file
                else: simulation_status = self.status_from_log_file_job(remote_log_file_path, ski_name)

            # Job not present in the queue anymore
            else: simulation_status = self.status_from_log_file_job(remote_log_file_path, ski_name)

        # Simulation is managed with an SQL database
        elif simulation.handle.type == "sql":

            # Get the simulation queue name
            queue_name = simulation.handle.value

            # Get the queue variable name in the remote python session
            queue_variable_name = "queue__" + queue_name

            # Get the status
            simulation_status = python_session.get_simple_property(queue_variable_name, "select('name=?', (simulation_name,))[0]['status']")

        # Invalid simulation handle
        else: raise RuntimeError("Unrecognized simulation handle")

        # Return the simulation status
        return simulation_status

    # -----------------------------------------------------------------

    def get_jobs_status(self, cluster_name=None):

        """
        This function ...
        :param cluster_name:
        :return:
        """

        # Only for specific cluster
        if cluster_name is not None: return self._get_jobs_status_cluster(cluster_name)

        # For all clusters
        else:

            status = dict()

            # Get status for currently active cluster
            active_cluster_name = self.active_cluster_name
            status_active_cluster = self._get_jobs_status_cluster(active_cluster_name, check_swap=False, swap=False, reswap=False)
            status.update(status_active_cluster)

            # Loop over the clusters
            for clustername in self.cluster_names:
                if clustername == active_cluster_name: continue # already done
                status_cluster = self._get_jobs_status_cluster(clustername, check_swap=False, reswap=False)
                status.update(status_cluster)

            # Swap back to original
            self.swap_cluster(active_cluster_name)

            # Return the status
            return status

    # -----------------------------------------------------------------

    def _get_jobs_status_cluster(self, cluster_name, check_swap=True, swap=True, reswap=True, timeout=5):

        """
        This function ...
        :param cluster_name:
        :param check_swap:
        :param swap:
        :param reswap:
        :param timeout:
        :return:
        """

        # Swap cluster
        if check_swap:
            if cluster_name != self.active_cluster_name:
                self.swap_cluster(cluster_name)
                swapped = True
            else: swapped = False
        elif swap:
            self.swap_cluster(cluster_name)
            swapped = True
        else: swapped = False

        # Create a dictionary that contains the status of the different jobs that are scheduled or running on the cluster
        queue_status = dict()

        # Obtain job status information through the 'qstat' command
        try: output = self.execute("qstat", timeout=timeout)
        except TimeOutReached:
            log.warning("The jobs status could not be obtained on cluster '" + cluster_name + "': time out reached")
            return queue_status

        # Check every line in the output
        for line in output:

            # If this line mentions a job
            if "master" in line:

                # Get the job ID
                jobid = int(line.split(".")[0])

                # Split the line
                splitted_line = line.split(" ")

                # Get the status (Q=queued, R=running)
                if "short" in splitted_line: position = splitted_line.index("short")
                elif "long" in splitted_line: position = splitted_line.index("long")
                else: continue
                jobstatus = splitted_line[position - 1]

                # Add the status of this job to the dictionary
                queue_status[jobid] = jobstatus

        # Swap back to original cluster?
        if swapped and reswap: self.swap_cluster(self.cluster_name)

        # Return the queue status
        return queue_status

    # -----------------------------------------------------------------

    def iterate_simulations(self, simulation_names=None):

        """
        This function ...
        :param simulation_names:
        :return:
        """

        from ..basics.configuration import prompt_proceed

        # Loop over the simulation files
        for path in fs.files_in_path(self.local_skirt_host_run_dir, extension="sim", sort=int):

            # Open the simulation file
            try: simulation = RemoteSimulation.from_file(path)
            except EOFError:

                message = "The simulation file '" + path + "' is not readable: perhaps serialization has failed"

                log.error(message)
                remove = prompt_proceed("remove the simulation file?")
                if remove: fs.remove_file(path)
                log.error("Skipping this simulation ...")
                continue

            # Name in list?
            if simulation_names is not None and simulation.name not in simulation_names: continue

            # Yield the simulation
            yield simulation

    # -----------------------------------------------------------------

    def get_status(self, as_dict=False, returns=("path", "status"), simulation_names=None, simulations=None):

        """
        This function ..
        :param as_dict:
        :param returns:
        :param simulation_names:
        :param simulations:
        :return:
        """

        # Initialize a list to contain the statuses
        entries = []

        # Initialize extra arguments
        jobs_status = None
        session = None
        states = None

        # If the remote has a scheduling system for launching jobs
        if self.scheduler:

            # Get the status of the jobs
            jobs_status = self.get_jobs_status()

            # Start a remote python session
            #if self.has_pts: session = self.load_queues_in_session()
            #else: session = None
            session = None

        # If the remote host does not use a scheduling system
        else: states = self.screen_states() # Get state of screen sessions

        # Get the simulations
        if simulations is None: simulations = self.iterate_simulations(simulation_names=simulation_names)
        else:
            if simulation_names is not None: raise ValueError("Cannot specify simulation names and simulation objects")
            if types.is_dictionary(simulations): simulations = simulations.values()
            elif types.is_sequence(simulations): pass
            else: raise ValueError("Container type of the simulations not recognized")

        # Loop over the simulations for this remote host
        for simulation in simulations:

            # Get the status of the simulation
            simulation_status = self.get_simulation_status(simulation, screen_states=states, jobs_status=jobs_status, session=session)

            # Construct item
            item = []
            for rtype in returns:
                if rtype == "name": item.append(simulation.name)
                elif rtype == "path": item.append(simulation.path)
                elif rtype == "id": item.append(id)
                elif rtype == "status": item.append(simulation_status)
                elif rtype == "simulation": item.append(simulation)
                else: raise ValueError("Unknown return variable '" + rtype + "'")

            # Add the simulation properties to the list
            entries.append(tuple(item))

        # Return the list of simulation properties
        if as_dict:
            if len(returns) != 2: raise ValueError("Cannot return as dict if number of return variables is not 2")
            return OrderedDict(entries)
        else: return entries

    # -----------------------------------------------------------------

    def load_queues_in_session(self):

        """
        This function ...
        :return:
        """

        # Start python session
        session = self.start_python_session(assume_pts=True)

        # Import statements
        session.import_package("SimulationQueue", from_name="pts.core.remote.queue")

        # Load all simulation queues
        for path, name in session.files_in_path(self.pts_run_path, extension="queue", returns=["path", "name"]):

            # Load the queue
            session.define_simple_variable("queue__" + name, "SimulationQueue.from_file('" + path + "'")

        # Return the session
        return session

    # -----------------------------------------------------------------

    def get_simulation_status(self, simulation, screen_states=None, jobs_status=None, session=None):

        """
        This function ...
        :param simulation:
        :param screen_states:
        :param jobs_status:
        :param session:
        :return:
        """

        # Uses scheduler?
        if self.scheduler:

            # Check jobs status
            if jobs_status is None: jobs_status = self.get_jobs_status(cluster_name=simulation.cluster_name)

            # Check session
            #if session is None:
                #if self.has_pts: session = self.load_queues_in_session()
                #else: session = None

            # Get the status
            simulation_status = self._get_simulation_status_scheduler(simulation, jobs_status, session)

        # Doesn't use scheduler
        else:

            # Check whether the handle is defined
            if simulation.handle is None:

                # Warning to get attention
                log.warning("Simulation '" + simulation.name + "' [" + str(simulation.id) + "] on remote host '" + self.host_id + "' doesn't appear to have an execution handle. Assuming it is still running in attached mode through another terminal.")
                simulation_status = running_name

            # Handle is defined
            else: simulation_status = self._get_simulation_status_not_scheduler(simulation, screen_states=screen_states)

        # If the simulation is finished, set the flag if necessary
        if is_finished_status(simulation_status) and not simulation.finished:
            simulation.finished = True
            simulation.save()

        # Return the status
        return simulation_status

    # -----------------------------------------------------------------

    def status_from_log_file(self, file_path, handle, simulation_prefix, screen_states=None, remote_ski_path=None):

        """
        This function ...
        :param file_path:
        :param handle:
        :param simulation_prefix:
        :param screen_states:
        :param remote_ski_path:
        :return:
        """

        # If the log file exists
        if self.is_file(file_path):

            # Get the last two lines of the remote log file
            output = self.read_last_lines(file_path, 2)

            # Get the last line of the actual simulation
            if len(output) == 0: return "invalid: cannot read log file"
            elif len(output) == 1: last = output[0]
            elif " Available memory: " in output[1]: last = output[0]
            else: last = output[1]

            # Interpret the content of the last line
            if " Finished simulation " + simulation_prefix in last: simulation_status = finished_name
            elif " *** Error: " in last: simulation_status = self.crashed_status_from_log_file(file_path)
            else:

                # Screen session
                if handle.type == "screen":

                    screen_name = handle.value

                    # Check whether screen is active
                    if screen_states is not None:
                        if screen_name not in screen_states: active_screen = False
                        else: active_screen = screen_states[screen_name] == "detached" or screen_states[screen_name] == "attached"
                    else: active_screen = self.is_active_screen(screen_name)

                    # Get status of simulation
                    if active_screen: simulation_status = self.running_status_from_log_file(file_path)
                    else: simulation_status = aborted_name

                # Attached terminal session
                elif handle.type == "tty":

                    session_rank = handle.value
                    if session_rank in self.ttys: simulation_status = self.running_status_from_log_file(file_path)
                    else: simulation_status = aborted_name

                # Invalid execution handle
                else: raise ValueError("Invalid execution handle")

        # If the log file does not exist, the simulation has not started yet or has been cancelled
        else:

            # Screen session
            if handle.type == "screen":

                # The simulation has not started or it's screen session has been cancelled
                screen_name = handle.value

                # Check whether screen is active
                if screen_states is not None:
                    if screen_name not in screen_states: active_screen = False
                    else: active_screen = screen_states[screen_name] == "detached" or screen_states[screen_name] == "attached"
                else: active_screen = self.is_active_screen(screen_name)

                # Set status of simulation
                if active_screen:

                    # # CAN ALSO BE CRASHED (OR OUTPUT HAS BEEN DELETED ACCIDENTLY)
                    # screen_script_path = handle.remote_screen_script_path
                    # screen_output_path = handle.remote_screen_output_path
                    #
                    # # IDEALLY WE WOULD LOCATE THE REMOTE PATH OF THE SCREEN SCRIPT
                    # # THIS WILL WORK IN THE FUTURE, JUST ADDED REMOTE SCREEN SCRIPT PATH TO EXECUTION HANDLES FOR FUTURE SIMULATIONS
                    # if screen_script_path is not None and self.is_file(screen_script_path) and screen_output_path is not None and remote_ski_path is not None:
                    #
                    #     #raise NotImplementedError("Not yet implemented")
                    #     screenlog_path = fs.join(screen_output_path, "screenlog.0")
                    #     last_line_starting_simulation = self.get_last_line_containing_memoize(screenlog_path, "Constructing a simulation from ski file")
                    #     if last_line_starting_simulation is None: simulation_status = "queued" # nothing has started yet
                    #     else:
                    #         ski_path = last_line_starting_simulation.split("from ski file '")[1].split("'")[0]
                    #         # Open the screen script
                    #         screen = ScreenScript.from_remote_file(screen_script_path, self)
                    #         #past_simulations = screen.get_simulations_before_ski_path(ski_path)
                    #         past_ski_paths = screen.get_ski_paths_before_ski_path(ski_path)
                    #         if remote_ski_path in past_ski_paths: simulation_status = "invalid: output missing"
                    #         else: simulation_status = "queued"
                    #
                    # # LOOK AT SCREEN OUTPUT TO SEE WHETHER SIMULATION HAS STARTED
                    # elif screen_output_path is not None and remote_ski_path is not None:
                    #
                    #     # VERY SLOW METHOD TO DO FOR EACH QUEUED SIMULATION
                    #
                    #     screenlog_path = fs.join(screen_output_path, "screenlog.0")
                    #     if self.is_file(screenlog_path) and self.contains_line(screenlog_path, "Constructing a simulation from ski file '" + remote_ski_path + "'"): simulation_status = "invalid: output missing"
                    #     else: simulation_status = "queued"
                    #
                    # # NO INFORMATION FROM THE SCREEN: ASSUME STILL QUEUED
                    # else: simulation_status = "queued"

                    simulation_status = queued_name

                # Screen is not active anymore
                else: simulation_status = cancelled_name

            # Attached terminal session
            elif handle.type == "tty":

                session_rank = handle.value
                if session_rank in self.ttys: simulation_status = queued_name
                else: simulation_status = cancelled_name

            # Invalid execution handle
            else: raise ValueError("Invalid execution handle")

        # Return the string that indicates the simulation status
        return simulation_status

    # -----------------------------------------------------------------

    def status_from_log_file_job(self, file_path, simulation_prefix):

        """
        This function ...
        :param file_path:
        :param simulation_prefix:
        :return:
        """

        # Check whether the log file exists
        if self.is_file(file_path):

            # Get the last two lines of the remote log file
            output = self.read_last_lines(file_path, 2)

            # Get the last line of the actual simulation
            if len(output) == 0: return "invalid: cannot read log file"
            if len(output) == 1: last = output[0]
            elif " Available memory: " in output[1]: last = output[0]
            else: last = output[1]

            # Interpret the content of the last line
            if " Finished simulation " + simulation_prefix in last: simulation_status = finished_name
            elif " *** Error: " in last: simulation_status = self.crashed_status_from_log_file(file_path)

            # The simulation cannot be running because we would have seen it in the qstat output
            # So with a partial log file, it must have been aborted
            else: simulation_status = aborted_name

        # If the log file does not exist, the simulation has been cancelled (if it would just be still scheduled
        # we would have encountered its job ID in the qstat output)
        else: simulation_status = cancelled_name

        # Return the string that indicates the simulation status
        return simulation_status

    # -----------------------------------------------------------------

    def crashed_status_from_log_file(self, file_path):

        """
        This function ...
        :param file_path:
        :return:
        """

        # Return string from simulation status
        status_string = str(LogSimulationStatus(file_path, self))
        if crashed_name not in status_string: raise RuntimeError("Something went wrong: status from log file is '" + status_string + "' but simulation is supposed to be crashed")
        return status_string

    # -----------------------------------------------------------------

    def running_status_from_log_file(self, file_path, nlibrary_entries=None):

        """
        This function ...
        :param nlibrary_entries:
        :return:
        """

        # Try to interpret the full status from the last log line
        status_string = self.running_status_from_last_log_line(file_path, nlibrary_entries=nlibrary_entries)
        if status_string is not None: return status_string

        # Return string from simulation status
        status_string = str(LogSimulationStatus(file_path, self))
        #print(status_string)
        if not is_running_status(status_string): raise RuntimeError("Something went wrong: status from log file is '" + status_string + "' but simulation is supposed to be running")
        return status_string

    # -----------------------------------------------------------------

    def running_status_from_last_log_line(self, file_path, nlibrary_entries=None):

        """
        This function ...
        :param file_path:
        :param nlibrary_entries:
        :return:
        """

        # Get the line
        line = self.get_last_line(file_path)

        # Return the status
        return get_status_from_last_log_line(line, nlibrary_entries=nlibrary_entries)

    # -----------------------------------------------------------------

    def _verify_scheduling_options(self, options, arguments, jobscript_dir_path=None):

        """
        This function ...
        :param options:
        :param arguments:
        :param jobscript_dir_path:
        :return:
        """

        # If scheduling options is not defined, create a new SchedulingOptions object
        if options is None: options = SchedulingOptions()

        # Test the presence of the 'nodes' and 'ppn' options
        if options.nodes is None or options.ppn is None:

            # The number of threads per core that should be used
            threads_per_core = self.threads_per_core if self.use_hyperthreading_skirt else 1

            # Get the requirements in number of nodes and ppn
            processors = arguments.parallel.processes * arguments.parallel.threads / threads_per_core
            assert math.ceil(processors) == math.floor(processors) # make sure is integer
            processors = int(processors)
            nodes, ppn = self.get_requirements(processors)

            # Set the nodes and pppn
            options.nodes = nodes
            options.ppn = ppn

        # Check if 'mail' option is defined
        if options.mail is None: options.mail = False

        # Check if 'full_node' option is defined
        if options.full_node is None: options.full_node = True

        # We want to estimate the walltime here if it is not defined in the options
        if options.walltime is None: raise RuntimeError("The walltime is not defined in the scheduling options")

        # Check if job script path is defined
        if options.local_jobscript_path is None:

            # Check whether path to directory where jobscript should be saved locally is given
            if jobscript_dir_path is None: raise ValueError("Jobscript directory path must be specified in this function if the local jobscript path has not been set in the scheduling options")

            # Determine the job script path
            local_jobscript_path = fs.join(jobscript_dir_path, time.unique_name("job") + ".sh")

            # Set the local jobscript path
            options.local_jobscript_path = local_jobscript_path

        # Return the scheduling options
        return options

# -----------------------------------------------------------------
