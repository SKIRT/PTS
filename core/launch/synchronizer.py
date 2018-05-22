#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.synchronizer Contains the RemoteSynchronizer class, which can be used to

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..simulation.simulation import RemoteSimulation
from ..remote.host import find_host_ids, has_simulations, has_tasks
from .analyser import SimulationAnalyser
from ..basics.configurable import Configurable
from ..simulation.remote import SKIRTRemote
from ..remote.remote import Remote
from ..tools import filesystem as fs
from ..basics.log import log
from ..basics.task import Task
from ..tools import formatting as fmt
from ..tools import introspection
from ..simulation.status import LogSimulationStatus
from ..simulation.remote import get_retrieved_simulations, get_retrieved_tasks, get_status_simulations, get_status_tasks
from ..simulation.remote import queued_name, running_name, finished_name, retrieved_name, analysed_name, aborted_name, cancelled_name, crashed_name, unknown_name, invalid_name

# -----------------------------------------------------------------

has_skirt = introspection.skirt_is_present()

# -----------------------------------------------------------------

class RemoteSynchronizer(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(RemoteSynchronizer, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Initialize a list to contain different SKIRTRemote instances for the different remote hosts
        self.remotes = []

        # List of host IDs (when remotes are not necessary)
        self._host_ids = None

        # The simulation results analyser
        self.analyser = None

        # Initialize a list to contain the retrieved simulations
        self.simulations = []

        # Initialize a list to contain the retrieved tasks
        self.tasks = []

        # Initialize a list for the runnning simulations
        self.running_simulations = []

        # Initialize a list for the running tasks
        self.running_tasks = []

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 2. Gather the simulations and tasks
        self.gather()

        # 3. Analyse
        if self.config.analyse: self.analyse()

        # 4. Announce the status of the simulations and tasks
        self.announce()

        # 5. Show the progress of the simulation
        if self.config.show_progress and not self.config.offline: self.show_progress()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(RemoteSynchronizer, self).setup(**kwargs)

        # Create the simulation analyser
        self.analyser = SimulationAnalyser(self.config.analysis)

        # Check flags
        if not self.config.retrieve and self.config.analyse: raise ValueError("Cannot analyse when retrieve is disabled")

        # Offline?
        if self.config.offline:

            # Set host IDS
            if "remotes" in kwargs: self.host_ids = [remote.host_id for remote in kwargs.pop("remotes")]
            elif "host_ids" in kwargs: self.host_ids = kwargs.pop("host_ids")
            elif self.config.host_ids is not None: self.host_ids = self.config.host_ids
            else: self.host_ids = find_host_ids()

        # Not offline
        else:

            # Load the remote instances
            if "remotes" in kwargs: self.remotes = kwargs.pop("remotes")
            elif "host_ids" in kwargs: self.remotes = [Remote(host_id=host_id) for host_id in kwargs.pop("host_id")]
            else:

                # Determine the host IDs
                if self.config.host_ids is not None: host_ids = self.config.host_ids
                else: host_ids = find_host_ids()

                # Loop over the host IDs
                for host_id in host_ids:

                    # If there are currently no simulations corresponding to this host, skip it
                    if (not has_simulations(host_id)) and (not has_tasks(host_id)): continue

                    # Create and setup a remote execution context
                    remote = Remote()
                    if not remote.setup(host_id):
                        log.warning("Remote host '" + host_id + "' is not available: skipping ...")
                        continue

                    # Setup SKIRT remote environment
                    if has_skirt and remote.has_skirt: remote = SKIRTRemote.from_remote(remote)

                    # Add the remote to the list of remote objects
                    self.remotes.append(remote)

        # Set the number of allowed open file handles
        fs.set_nallowed_open_files(self.config.nopen_files)

        # Set analyser options
        self.analyser.config.ignore_missing_data = self.config.ignore_missing_data

    # -----------------------------------------------------------------

    @property
    def nremotes(self):

        """
        This function ...
        :return:
        """

        return len(self.remotes)

    # -----------------------------------------------------------------

    @property
    def has_remotes(self):

        """
        This function ...
        :return:
        """

        return self.nremotes > 0

    # -----------------------------------------------------------------

    @property
    def host_ids(self):

        """
        This function ...
        :return:
        """

        if self._host_ids is not None: return self._host_ids
        else: return [remote.host_id for remote in self.remotes]

    # -----------------------------------------------------------------

    @host_ids.setter
    def host_ids(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        if self.has_remotes: raise ValueError("Cannot set host IDS when there are remote instances")
        self._host_ids = value

    # -----------------------------------------------------------------

    def get_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        #if not self.has_remotes: raise ValueError("No remotes")
        for remote in self.remotes:
            if remote.host_id == host_id: return remote
        #raise ValueError("Remote for host '" + host_id + "' not found")

        # Create new remote (shouldn't happen)
        return Remote(host_id=host_id)

    # -----------------------------------------------------------------

    @property
    def nsimulations(self):

        """
        Thisnfunction ...
        :return:
        """

        return len(self.simulations)

    # -----------------------------------------------------------------

    @property
    def nrunning_simulations(self):

        """
        This function ...
        :return:
        """

        return len(self.running_simulations)

    # -----------------------------------------------------------------

    @property
    def has_single_simulation(self):

        """
        This function ...
        :return:
        """

        return self.nsimulations == 1

    # -----------------------------------------------------------------

    @property
    def has_single_running_simulation(self):

        """
        This function ...
        :return:
        """

        return self.nrunning_simulations == 1

    # -----------------------------------------------------------------

    @property
    def single_simulation(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_simulation: raise ValueError("Not a single simulation")
        return self.simulations[0]

    # -----------------------------------------------------------------

    @property
    def single_running_simulation(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_running_simulation: raise ValueError("Not a single running simulation")
        return self.running_simulations[0]

    # -----------------------------------------------------------------

    @property
    def ntasks(self):

        """
        This function ...
        :return:
        """

        return len(self.tasks)

    # -----------------------------------------------------------------

    @property
    def nrunning_tasks(self):

        """
        This function ...
        :return:
        """

        return len(self.running_tasks)

    # -----------------------------------------------------------------

    @property
    def has_single_task(self):

        """
        This function ...
        :return:
        """

        return self.ntasks == 1

    # -----------------------------------------------------------------

    @property
    def has_single_running_task(self):

        """
        This function ...
        :return:
        """

        return self.nrunning_tasks == 1

    # -----------------------------------------------------------------

    @property
    def single_task(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_task: raise ValueError("Not a single task")
        return self.tasks[0]

    # -----------------------------------------------------------------

    @property
    def single_running_task(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_running_task: raise ValueError("Not a single running task")
        return self.running_tasks[0]

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the synchronizer...")

        # Log out from the remotes
        for remote in self.remotes: remote.logout()

        # Clear the list of remotes
        self.remotes = []

        # Clear the list of simulations
        self.simulations = []

        # Clear the list of tasks
        self.tasks = []

    # -----------------------------------------------------------------

    def get_ids_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        if self.config.ids is None: return None
        elif host_id not in self.config.ids: return None
        else: return self.config.ids[host_id]

    # -----------------------------------------------------------------

    def gather(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Gathering the SKIRT simulations and PTS tasks ...")

        # Retrieve
        if self.config.retrieve and not self.config.offline: self.retrieve()

        # Just load with current status
        else: self.load()

    # -----------------------------------------------------------------

    def load(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading finished SKIRT simulations and PTS tasks ...")

        # Load SKIRT simulations
        if self.config.simulations and has_skirt: self.load_simulations()

        # Load PTS tasks
        if self.config.tasks: self.load_tasks()

    # -----------------------------------------------------------------

    def load_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading SKIRT simulations ...")

        # Loop over the host IDS
        for host_id in self.host_ids:

            # Inform the user
            log.debug("Loading the retrieved simulations of remote host '" + host_id + "' ...")

            # Get retrieved simulations
            self.simulations += get_retrieved_simulations(host_id)

    # -----------------------------------------------------------------

    def load_tasks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading PTS tasks ...")

        # Loop over the host IDS
        for host_id in self.host_ids:

            # Inform the user
            log.debug("Loading the retrieved tasks of remote host '" + host_id + "' ...")

            # Get retrieved tasks
            self.tasks += get_retrieved_tasks(host_id)

    # -----------------------------------------------------------------

    def get_retrieve_crashed_ids_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        if self.config.retrieve_crashed is None: return None
        elif host_id not in self.config.retrieve_crashed: return None
        else: return self.config.retrieve_crashed[host_id]

    # -----------------------------------------------------------------

    def retrieve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Retrieving finished SKIRT simulations and PTS tasks ...")

        # Retrieve SKIRT simulations
        if self.config.simulations and has_skirt: self.retrieve_simulations()

        # Retrieve PTS tasks
        if self.config.tasks: self.retrieve_tasks()

    # -----------------------------------------------------------------

    def retrieve_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Retrieving the output of finished simulations ...")

        # Loop over the different remotes
        for remote in self.remotes:

            # Check whether SKIRT is present on the remote
            if not remote.has_skirt: continue

            # Inform the user
            log.debug("Retrieving the simulations of remote '" + remote.system_name + "' ...")

            # Retrieve simulations
            self.simulations += remote.retrieve(retrieve_crashed=self.get_retrieve_crashed_ids_for_host(remote.host_id),
                                                check_crashed=self.config.check_crashed, check_data=self.config.check_data)

    # -----------------------------------------------------------------

    def retrieve_tasks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Retrieving tasks ...")

        # Loop over the different remotes
        for remote in self.remotes:

            # Inform the user
            log.debug("Retrieving the tasks of remote '" + remote.system_name + "' ...")

            # Retrieve tasks
            self.tasks += remote.retrieve_tasks()

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the output of retrieved SKIRT simulations and PTS tasks ...")

        # Analyse the output of the retrieved simulations
        if self.config.simulations and has_skirt: self.analyse_simulations()

        # Analyse the output of the retrieved tasks
        if self.config.tasks: self.analyse_tasks()

    # -----------------------------------------------------------------

    def analyse_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing simulations ...")

        # Loop over the list of simulations and analyse them
        for simulation in self.simulations:

            # Run the analyser on the simulation
            self.analyser.run(simulation=simulation)

            # Clear the simulation analyser
            self.analyser.clear()

    # -----------------------------------------------------------------

    def analyse_tasks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing tasks ...")

        # Loop over the list of retrieved tasks
        for task in self.tasks:

            # Analyse the task
            task.analyse()

    # -----------------------------------------------------------------

    def announce(self):

        """
        This function ...
        :return:
        """

        # Announce the status of the SKIRT simulations
        if self.config.simulations and has_skirt: self.announce_simulations()

        # Announce the status of the PTS tasks
        if self.config.tasks: self.announce_tasks()

    # -----------------------------------------------------------------

    def announce_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        #log.info("SKIRT simulations:")

        #first = True

        # Loop over the different hosts
        for host_id in self.host_ids:

            # Offline?
            if self.config.offline:

                status = get_status_simulations(host_id)
                remote = None

            # Online
            else:

                # Get remote
                remote = self.get_remote(host_id)

                # Check whether SKIRT is present
                if not remote.has_skirt: continue

                # Get the status of the different simulations
                status = remote.get_status()

            # Show the name of the current remote
            if len(status) > 0: log.info("Simulations on remote '" + host_id + "':")
            print()

            # Get the status of the different simulations
            for path, simulation_status in status:

                #if first:
                #    # Inform the user
                #    log.info("SKIRT simulations:")
                #    first = False

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                prefix = " - "
                tag = "[" + str(simulation.id) + "]"

                # Finished, retrieved and analysed simulation (remote output has already been removed, if requested)
                if simulation_status.startswith(analysed_name):

                    if (self.config.ids is not None and (
                            host_id in self.config.ids and simulation.id in self.config.ids[host_id])) \
                            or (self.config.statuses is not None and analysed_name in self.config.statuses):

                        tag = "[ X ]"

                        # Remove the simulation file
                        fs.remove_file(path)

                    formatter = fmt.green

                # Finished and retrieved simulation (remote output has already been removed, if requested)
                elif simulation_status == retrieved_name:

                    if (self.config.ids is not None and (host_id in self.config.ids and simulation.id in self.config.ids[host_id]))\
                            or (self.config.statuses is not None and retrieved_name in self.config.statuses):

                        tag = "[ X ]"

                        # Remove the simulation file
                        fs.remove_file(path)

                    formatter = fmt.green

                # Finished, but not yet retrieved simulation
                elif simulation_status == finished_name:

                    if (self.config.ids is not None and (
                            host_id in self.config.ids and simulation.id in self.config.ids[host_id])) \
                            or (self.config.statuses is not None and finished_name in self.config.statuses):
                        log.warning(
                            "The simulation with ID " + str(simulation.id) + " has finished, but has not been"
                                                                             " retrieved yet. Deleting it now would mean all simulation output is lost. Run "
                                                                             " 'pts status' again to retrieve the simulation output.")

                    formatter = fmt.blue

                    simulation_status += " (do 'pts status' again to retrieve)"

                # Running simulation
                elif running_name in simulation_status:

                    # ADD TO RUNNNING SIMULATIONS
                    self.running_simulations.append(simulation)

                    if (self.config.ids is not None and (host_id in self.config.ids and simulation.id in self.config.ids[host_id]))\
                            or (self.config.statuses is not None and running_name in self.config.statuses):

                        if remote.host.scheduler:

                            tag = "[ X ]"

                            remote.kill_job(simulation.id)

                            # Remove the simulation file
                            fs.remove_file(path)

                            # Remove the remote input, output and simulation directory
                            if simulation.has_input: remote.remove_directory(simulation.remote_input_path)
                            remote.remove_directory(simulation.remote_output_path)
                            remote.remove_directory(simulation.remote_simulation_path)

                            simulation_status += " -> aborted"

                        else: log.warning("Aborting simulations not running on a host with a scheduling system is not"
                                          " implemented yet. ")

                    formatter = fmt.reset

                # Simulations with invalid state
                elif invalid_name in simulation_status:

                    formatter = fmt.red + fmt.bold

                # Simulations with unknown status (because offline)
                elif simulation_status == unknown_name:

                    formatter = fmt.lightyellow

                # Crashed simulation
                elif crashed_name in simulation_status:

                    if (self.config.ids is not None and (host_id in self.config.ids and simulation.id in self.config.ids[host_id]))\
                            or (self.config.statuses is not None and crashed_name in self.config.statuses):

                        tag = "[ X ]"

                        # Remove the simulation file
                        fs.remove_file(path)

                        # Remove the remote input, output and simulation directory
                        if simulation.has_input: remote.remove_directory(simulation.remote_input_path)
                        remote.remove_directory(simulation.remote_output_path)
                        remote.remove_directory(simulation.remote_simulation_path)

                    formatter = fmt.lightred

                # Cancelled simulation
                elif simulation_status == cancelled_name:

                    if (self.config.ids is not None and (host_id in self.config.ids and simulation.id in self.config.ids[host_id]))\
                            or (self.config.statuses is not None and cancelled_name in self.config.statuses):

                        tag = "[ X ]"

                        # Remove the simulation file
                        fs.remove_file(path)

                        # Remove the remote input, output and simulation directory
                        if simulation.has_input: remote.remove_directory(simulation.remote_input_path)
                        remote.remove_directory(simulation.remote_output_path)
                        remote.remove_directory(simulation.remote_simulation_path)

                    formatter = fmt.lightyellow

                # Aborted simulation
                elif simulation_status == aborted_name:

                    if (self.config.ids is not None and (host_id in self.config.ids and simulation.id in self.config.ids[host_id]))\
                            or (self.config.statuses is not None and aborted_name in self.config.statuses):

                        tag = "[ X ]"

                        # Remove the simulation file
                        fs.remove_file(path)

                        # Remove the remote input, output and simulation directory
                        if simulation.has_input: remote.remove_directory(simulation.remote_input_path)
                        remote.remove_directory(simulation.remote_output_path)
                        remote.remove_directory(simulation.remote_simulation_path)

                    formatter = fmt.lightyellow

                # Queued simulation
                elif simulation_status == queued_name:

                    if (self.config.ids is not None and (host_id in self.config.ids and simulation.id in self.config.ids[host_id]))\
                            or (self.config.statuses is not None and queued_name in self.config.statuses):

                        if remote.host.scheduler:

                            tag = "[ X ]"

                            # Stop the simulation
                            remote.stop_job(simulation.id)

                            # Remove the simulation file
                            fs.remove_file(path)

                            # Remove the remote input, output and simulation directory
                            if simulation.has_input: remote.remove_directory(simulation.remote_input_path)
                            remote.remove_directory(simulation.remote_output_path)
                            remote.remove_directory(simulation.remote_simulation_path)

                            simulation_status += " -> cancelled"

                        else: log.warning("Cancelling simulations not running on a host with a scheduling system is not"
                                          " implemented yet. ")

                    formatter = fmt.reset

                # Other
                else: formatter = fmt.reset

                # Show the status of the current simulation
                print(formatter + prefix + tag + " " + simulation.name + ": " + simulation_status + fmt.reset)

            print()

    # -----------------------------------------------------------------

    def announce_tasks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        #log.info("PTS tasks:")

        first = True

        # Loop over the different hosts
        for host_id in self.host_ids:

            # Offline?
            if self.config.offline:

                status = get_status_tasks(host_id)
                remote = None

            # Online?
            else:

                # Get remote
                remote = self.get_remote(host_id)

                # Get the status of the different tasks
                status = remote.get_task_status()

            # Show the name of the current remote
            if len(status) > 0: log.info("Tasks on remote '" + host_id + "':")
            print()

            # Get the status of the different tasks
            for path, task_status in status:

                if first:
                    # Inform the user
                    log.info("PTS tasks:")
                    first = False

                # Open the task file
                task = Task.from_file(path)

                prefix = " - "
                tag = "[" + str(task.id) + "]"

                # Tasks with invalid state
                if "invalid" in task_status:

                    formatter = fmt.red + fmt.bold

                # Simulations with unknown status (because offline)
                elif task_status == "unknown":

                    formatter = fmt.lightyellow

                elif "crashed" in task_status:

                    formatter = fmt.lightred

                # Finished, retrieved and analysed simulation (remote output has already been removed, if requested)
                elif task_status == "analysed":

                    # ...

                    formatter = fmt.green

                # Retrieved tasks (remote output has already been removed, if requested)
                elif task_status == "retrieved":

                #if (self.config.ids is not None and (
                #        remote.host.id in self.config.ids and simulation.id in self.config.ids[remote.host.id])) \
                #        or (self.config.statuses is not None and "retrieved" in self.config.statuses):
                #    tag = "[ X ]"

                    # Remove the simulation file
                    #fs.remove_file(path)

                    formatter = fmt.green

                # Finished, but not yet retrieved task
                elif task_status == "finished":

                    formatter = fmt.blue

                    task_status += " (do 'pts status' again to retrieve)"

                # Running task
                elif "running" in task_status:

                    formatter = fmt.reset

                    # ADD TO RUNNING TASKS
                    self.running_tasks.append(task)

                # Cancelled task
                elif task_status == "cancelled":

                    formatter = fmt.lightyellow

                # Aborted task
                elif task_status == "aborted":

                    formatter = fmt.lightyellow

                # Queued task
                elif task_status == "queued":

                    formatter = fmt.reset

                else: formatter = fmt.reset

                # Show the status of the current task
                print(formatter + prefix + tag + " " + task.name + ": " + task_status + fmt.reset)

            print()

    # -----------------------------------------------------------------

    def show_progress(self):

        """
        This function ...
        :return:
        """

        # Check
        if not self.has_single_running_simulation:
            log.warning("Cannot show the progress when there are multiple simulations still running")
            return

        # Get the simulation
        simulation = self.single_running_simulation

        # Inform the user
        log.info("Show the progress of the simulation '" + simulation.name + "' ...")

        # Create the status object
        status = LogSimulationStatus(simulation.remote_log_file_path, remote=self.get_remote(simulation.host_id), debug_output=self.config.debug_output)

        # Show the simulation progress
        with log.no_debugging(): success = status.show_progress(simulation.handle)

        # Check whether not crashed
        if not success: raise RuntimeError("The simulation crashed")

# -----------------------------------------------------------------
