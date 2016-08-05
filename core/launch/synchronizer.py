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
from ..basics.host import find_host_ids, has_simulations, has_tasks
from .analyser import SimulationAnalyser
from ..basics.configurable import OldConfigurable
from ..simulation.remote import SkirtRemote
from ..tools import filesystem as fs
from ..tools.logging import log
from ..basics.task import Task
from ..tools import formatting as fmt

# -----------------------------------------------------------------

class RemoteSynchronizer(OldConfigurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(RemoteSynchronizer, self).__init__(config, "core")

        # -- Attributes --

        # Initialize a list to contain different SkirtRemote instances for the different remote hosts
        self.remotes = []

        # The simulation results analyser
        self.analyser = SimulationAnalyser()

        # Initialize a list to contain the retrieved simulations
        self.simulations = []

        # Initialize a list to contain the retrieved tasks
        self.tasks = []

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new RemoteSynchronizer instance
        synchronizer = cls()

        ## Adjust the configuration settings according to the command-line arguments

        # Set the remote name and the delete dictionary
        if hasattr(arguments, "remote"): synchronizer.config.remote = arguments.remote
        if hasattr(arguments, "ids"): synchronizer.config.ids = arguments.ids
        if hasattr(arguments, "status"): synchronizer.config.statuses = arguments.status
        if hasattr(arguments, "relaunch"): synchronizer.config.relaunch = arguments.relaunch

        # Return the new synchronizer
        return synchronizer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Retrieve the simulations and tasks
        self.retrieve()

        # 3. Analyse
        self.analyse()

        # 4. Announce the status of the simulations
        self.announce()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(RemoteSynchronizer, self).setup()

        # Find host ids for which a configuration file has been created by the user
        for host_id in find_host_ids():

            # If a list of remotes is defined and this remote is not in it, skip it
            if self.config.remote is not None and host_id not in self.config.remote: continue

            # If there are currently no simulations corresponding to this host, skip it
            if (not has_simulations(host_id)) and (not has_tasks(host_id)): continue

            # Create a remote SKIRT execution context
            remote = SkirtRemote()

            # Setup the remote execution context
            remote.setup(host_id)

            # Add the remote to the list of remote objects
            self.remotes.append(remote)

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

    def retrieve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Retrieving finished SKIRT simulations and PTS tasks ...")

        # Retrieve SKIRT simulations
        self.retrieve_simulations()

        # Retrieve PTS tasks
        self.retrieve_tasks()

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

            # Inform the user
            log.debug("Retrieving the simulations of remote '" + remote.system_name + "' ...")

            # Retrieve simulations
            self.simulations += remote.retrieve()

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
        self.analyse_simulations()

        # Analyse the output of the retrieved tasks
        self.analyse_tasks()

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
            self.analyser.run(simulation)

            # Clear the analyser
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

            # Loop over the analyser classes that are defined for this task
            for analyser_class in task.analyser_classes:

                # Create an instance of the analyser class
                analyser = analyser_class.for_task(task)

                # Run the analyser
                analyser.run()

            # Set analysed flag to True
            task.analysed = True
            task.save()

    # -----------------------------------------------------------------

    def announce(self):

        """
        This function ...
        :return:
        """

        # Announce the status of the SKIRT simulations
        self.announce_simulations()

        # Announce the status of the PTS tasks
        self.announce_tasks()

    # -----------------------------------------------------------------

    def announce_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("SKIRT simulations:")

        # Loop over the different remotes
        for remote in self.remotes:

            # Get the status of the different simulations
            status = remote.get_status()

            # Show the name of the current remote
            if len(status) > 0: log.info("Simulations on remote '" + remote.host_id + "':")
            print()

            # Get the status of the different simulations
            for path, simulation_status in status:

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                prefix = " - "
                tag = "[" + str(simulation.id) + "]"

                # Finished, retrieved and analysed simulation (remote output has already been removed, if requested)
                if simulation_status == "analysed":

                    if (self.config.ids is not None and (
                            remote.host.id in self.config.ids and simulation.id in self.config.ids[remote.host.id])) \
                            or (self.config.statuses is not None and "retrieved" in self.config.statuses):

                        tag = "[ X ]"

                        # Remove the simulation file
                        fs.remove_file(path)

                    formatter = fmt.green

                # Finished and retrieved simulation (remote output has already been removed, if requested)
                elif simulation_status == "retrieved":

                    if (self.config.ids is not None and (remote.host.id in self.config.ids and simulation.id in self.config.ids[remote.host.id]))\
                            or (self.config.statuses is not None and "retrieved" in self.config.statuses):

                        tag = "[ X ]"

                        # Remove the simulation file
                        fs.remove_file(path)

                    formatter = fmt.green

                # Finished, but not yet retrieved simulation
                elif simulation_status == "finished":

                    if (self.config.ids is not None and (
                            remote.host.id in self.config.ids and simulation.id in self.config.ids[remote.host.id])) \
                            or (self.config.statuses is not None and "finished" in self.config.statuses):
                        log.warning(
                            "The simulation with ID " + str(simulation.id) + " has finished, but has not been"
                                                                             " retrieved yet. Deleting it now would mean all simulation output is lost. Run "
                                                                             " 'pts status' again to retrieve the simulation output.")

                    formatter = fmt.blue

                    simulation_status += " (do 'pts status' again to retrieve)"

                # Running simulation
                elif "running" in simulation_status:

                    if (self.config.ids is not None and (remote.host.id in self.config.ids and simulation.id in self.config.ids[remote.host.id]))\
                            or (self.config.statuses is not None and "running" in self.config.statuses):

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

                # Crashed simulation
                elif simulation_status == "crashed":

                    if (self.config.ids is not None and (remote.host.id in self.config.ids and simulation.id in self.config.ids[remote.host.id]))\
                            or (self.config.statuses is not None and "crashed" in self.config.statuses):

                        tag = "[ X ]"

                        # Remove the simulation file
                        fs.remove_file(path)

                        # Remove the remote input, output and simulation directory
                        if simulation.has_input: remote.remove_directory(simulation.remote_input_path)
                        remote.remove_directory(simulation.remote_output_path)
                        remote.remove_directory(simulation.remote_simulation_path)

                    formatter = fmt.lightred

                # Cancelled simulation
                elif simulation_status == "cancelled":

                    if (self.config.ids is not None and (remote.host.id in self.config.ids and simulation.id in self.config.ids[remote.host.id]))\
                            or (self.config.statuses is not None and "cancelled" in self.config.statuses):

                        tag = "[ X ]"

                        # Remove the simulation file
                        fs.remove_file(path)

                        # Remove the remote input, output and simulation directory
                        if simulation.has_input: remote.remove_directory(simulation.remote_input_path)
                        remote.remove_directory(simulation.remote_output_path)
                        remote.remove_directory(simulation.remote_simulation_path)

                    formatter = fmt.lightyellow

                # Aborted simulation
                elif simulation_status == "aborted":

                    if (self.config.ids is not None and (remote.host.id in self.config.ids and simulation.id in self.config.ids[remote.host.id]))\
                            or (self.config.statuses is not None and "aborted" in self.config.statuses):

                        tag = "[ X ]"

                        # Remove the simulation file
                        fs.remove_file(path)

                        # Remove the remote input, output and simulation directory
                        if simulation.has_input: remote.remove_directory(simulation.remote_input_path)
                        remote.remove_directory(simulation.remote_output_path)
                        remote.remove_directory(simulation.remote_simulation_path)

                    formatter = fmt.lightyellow

                # Queued simulation
                elif simulation_status == "queued":

                    if (self.config.ids is not None and (remote.host.id in self.config.ids and simulation.id in self.config.ids[remote.host.id]))\
                            or (self.config.statuses is not None and "queued" in self.config.statuses):

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
        log.info("PTS tasks:")

        # Loop over the different remotes
        for remote in self.remotes:

            # Get the status of the different tasks
            status = remote.get_task_status()

            # Show the name of the current remote
            if len(status) > 0: log.info("Tasks on remote '" + remote.host_id + "':")
            print()

            # Get the status of the different tasks
            for path, task_status in status:

                # Open the task file
                task = Task.from_file(path)

                prefix = " - "
                tag = "[" + str(task.id) + "]"

                # Tasks with invalid state
                if "invalid" in task_status:

                    formatter = fmt.red + fmt.bold

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
