#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.status Contains the SimulationStatus class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..tools import time
from ..tools.progress import Bar, BAR_FILLED_CHAR, BAR_EMPTY_CHAR
from ..tools.logging import log
from .logfile import get_last_phase, get_nprocesses, get_simulation_phase
from ..basics.handle import ExecutionHandle

# -----------------------------------------------------------------

phase_descriptions = dict()
phase_descriptions["setup"] = "setup of the simulation"
phase_descriptions["stellar"] = "emission of stellar photons"
phase_descriptions["spectra"] = "calculation of dust emission spectra"
phase_descriptions["dust"] = "emission of dust photons"
phase_descriptions["comm"] = "communication"
phase_descriptions["write"] = "writing results"
phase_descriptions["wait"] = "waiting for synchronization"

# -----------------------------------------------------------------

class SimulationStatus(object):

    """
    This class ...
    """

    def __init__(self, log_path, remote=None):

        """
        The constructor ...
        :param log_path:
        :param remote:
        """

        # Set attributes
        self.log_path = log_path
        self.remote = remote

        # The status
        self.status = None

        # The phase
        self.phase = None

        # The simulation phase
        self.simulation_phase = None

        # The stage (if applicable)
        self.stage = None

        # The cycle (if applicable)
        self.cycle = None

        # The progress (if applicable)
        self.progress = None

        # More info
        self.extra = None

        # Get the status
        self.refresh()

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        string = self.status

        # Add phase
        if self.phase is not None: string += ": " + phase_descriptions[self.phase]

        # Add simulation phase
        if self.simulation_phase is not None: string += " in the " + self.simulation_phase.lower() + " phase"

        # Add stage, cycle
        if self.stage is not None: string += " stage " + str(self.stage + 1)
        if self.cycle is not None: string += " [cycle " + str(self.cycle) + "]"

        # Add progress
        if self.progress is not None: string += " " + str(self.progress) + "%"

        # Return the string
        return string

    # -----------------------------------------------------------------

    @property
    def running(self):

        """
        This function ...
        :return:
        """

        return self.status == "running"

    # -----------------------------------------------------------------

    @property
    def finished(self):

        """
        This function ...
        :return:
        """

        return self.status == "finished"

    # -----------------------------------------------------------------

    @property
    def crashed(self):

        """
        This function ...
        :return:
        """

        return self.status == "crashed"

    # -----------------------------------------------------------------

    @property
    def aborted(self):

        """
        This function ...
        :return:
        """

        return self.status == "aborted"

    # -----------------------------------------------------------------

    @property
    def not_started(self):

        """
        This function ...
        :return:
        """

        return self.status == "not started"

    # -----------------------------------------------------------------

    def show_progress(self, process_or_handle, refresh_time=3):

        """
        This function ...
        :param process_or_handle:
        :param refresh_time: time in seconds
        :return:
        """

        last_phase = None
        last_stage = None
        last_cycle = None
        last_extra = None

        # Wait for a bit ?

        # Refresh loop
        while True:

            # Not yet started
            if self.not_started:

                # Execution handle
                if isinstance(process_or_handle, ExecutionHandle):

                    if process_or_handle.type != "screen": raise NotImplementedError("Execution handle must be 'screen'")
                    screen_name = process_or_handle.value

                    # Check the screen status
                    if not self.remote.is_active_screen(screen_name):
                        log.error("The simulation has stopped")
                        name = "screenlog.0"
                        screenlog_path = None
                        if process_or_handle.remote_screen_output_path is not None:
                            screenlog_path = fs.join(process_or_handle.remote_screen_output_path, name)
                            if not self.remote.is_file(screenlog_path): screenlog_path = None
                        if self.remote.is_file(self.log_path):
                            lines = list(self.remote.read_lines(self.log_path))
                            for line in lines: log.error(line)
                        elif screenlog_path is not None:
                            lines = list(self.remote.read_lines(screenlog_path))
                            for line in lines: log.error(line)
                        else: log.error("Log file hasn't been created yet and no screen output found")
                        return False

                # Process
                else:

                    returncode = process_or_handle.poll()
                    if returncode is not None:
                        log.error("The simulation has stopped")
                        if fs.is_file(self.log_path):
                            lines = list(fs.read_lines(self.log_path))
                            for line in lines: log.error(line)
                        else: log.error("Log file hasn't been created yet")
                        return False

                log.info("Waiting for simulation to start ...")
                time.wait(refresh_time)
                self.refresh(process_or_handle)
                continue

            # Finished: break loop
            if self.finished:
                log.success("Simulation finished")
                return True

            # Crashed: break loop
            elif self.crashed:
                log.error("Simulation crashed")
                return False

            # Aborted: break loop
            elif self.aborted:
                log.error("Simulation has been aborted")
                return False

            # New phase
            elif last_phase is None or self.phase != last_phase:

                last_phase = self.phase
                if last_phase is None:
                    time.wait(refresh_time)
                    self.refresh(process_or_handle)
                    continue
                else:
                    if self.simulation_phase is not None: log.info("Starting " + phase_descriptions[last_phase] + " in " + self.simulation_phase.lower() + " phase ...")
                    else: log.info("Starting " + phase_descriptions[last_phase] + " ...")

            # Self absorption phase
            if self.phase == "dust" and self.simulation_phase == "DUST SELF-ABSORPTION":

                total_length = 100

                if self.stage is None:
                    self.refresh_after(1)
                    continue

                if last_stage is None or self.stage != last_stage:
                    log.info("Starting stage " + str(self.stage + 1) + " ...")
                    last_stage = self.stage

                if self.cycle is None:
                    self.refresh_after(1)
                    continue

                if last_cycle is None or self.cycle != last_cycle:

                    log.info("Starting cycle " + str(self.cycle) + " ...")
                    last_cycle = self.cycle

                    # Progress bar
                    with Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                             filled_char=BAR_FILLED_CHAR, expected_size=total_length, every=1, add_datetime=True) as bar:

                        # Loop
                        while True:

                            if self.stage != last_stage:
                                bar.show(100) # make sure it always ends on 100%
                                break
                            if self.cycle != last_cycle:
                                bar.show(100) # make sure it always ends on 100%
                                break
                            if self.progress is None: bar.show(100)
                            else: bar.show(int(self.progress))
                            self.refresh_after(1)

                    self.refresh_after(refresh_time)

                else: self.refresh_after(refresh_time)

            # Stellar emission: show progress bar
            elif self.phase == "stellar" or self.phase == "spectra" or self.phase == "dust":

                total_length = 100

                # Progress bar
                with Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                         filled_char=BAR_FILLED_CHAR, expected_size=total_length, every=1, add_datetime=True) as bar:
                    # Loop
                    while True:
                        if self.phase != last_phase:
                            bar.show(100) # make sure it always ends on 100%
                            break
                        if self.progress is None:
                            bar.show(100)
                        else: bar.show(int(self.progress))
                        self.refresh_after(1)

            # Still the same phase
            else:
                if self.extra is not None:
                    if last_extra is None or last_extra != self.extra:
                        log.info(self.extra)
                    last_extra = self.extra
                self.refresh_after(refresh_time)

    # -----------------------------------------------------------------

    def refresh_after(self, seconds, process_or_handle=None):

        """
        This function ...
        :param seconds:
        :param process_or_handle:
        :return:
        """

        # Wait and refresh
        time.wait(seconds)
        self.refresh(process_or_handle)

    # -----------------------------------------------------------------

    def refresh(self, process_or_handle=None):

        """
        This function ...
        :param process_or_handle:
        :return:
        """

        # Clear everything
        self.status = None
        self.phase = None
        self.simulation_phase = None
        self.stage = None
        self.cycle = None
        self.progress = None
        self.extra = None

        # Check whether the log file exists
        if self.remote is not None: present = self.remote.is_file(self.log_path)
        else: present = fs.is_file(self.log_path)

        # If not exists, not started
        if not present:
            self.status = "not started"
            return

        # Get the log file lines
        if self.remote is not None: lines = list(self.remote.read_lines(self.log_path))
        else: lines = list(fs.read_lines(self.log_path))

        if len(lines) == 0:
            self.status = "invalid: cannot read log file"
            return

        # Check whether finished
        last_two_lines = lines[-2:]
        if len(lines) == 1: last = lines[0]
        elif " Available memory: " in last_two_lines[1]: last = last_two_lines[0]
        else: last = last_two_lines[1]

        # Interpret the content of the last line
        if " Finished simulation " in last:
            self.status = "finished"
            return

        elif " *** Error: " in last:
            self.status = "crashed"
            return

        # Running
        else:

            # Process or handle is passed: check whether the process or screen is still running and not aborted
            if process_or_handle is not None:

                # Execution handle
                if isinstance(process_or_handle, ExecutionHandle):

                    if process_or_handle.type != "screen": raise NotImplementedError("Execution handle must be 'screen'")
                    screen_name = process_or_handle.value

                    if not self.remote.is_active_screen(screen_name):
                        self.status = "aborted"
                        return

                # Process
                else:

                    returncode = process_or_handle.poll()
                    if returncode is not None:
                        #log.error("The simulation has been aborted")
                        #return False
                        self.status = "aborted"
                        return

            # Status is 'running'
            self.status = "running"

            # Get the last phase in the log file
            last_phase, start_index = get_last_phase(lines)

            # Set the phase
            self.phase = last_phase

            # Get the 'real' simulation phase
            #start_line = lines[start_index]
            #simulation_phase = get_simulation_phase(start_line, self.simulation_phase)

            # Set the simulation phase
            #self.simulation_phase = simulation_phase

            # The setup phase
            if self.phase == "setup":

                self.simulation_phase = "SETUP"

                # Loop over the lines in reversed order
                for line in reversed(lines):

                    if "Writing dust cell properties" in line:
                        self.extra = "Writing dust cell properties ..."
                        break
                    elif "Setting the value of the density in the cells" in line:
                        self.extra = "Calculating dust densities ..."
                        break
                    elif "Writing data to plot the dust grid":
                        self.extra = "Writing the dust grid ..."
                        break
                    elif "Subdividing node" in line:
                        self.extra = "Creating the dust grid tree ..."
                        break
                    elif "Starting subdivision" in line:
                        self.extra = "Creating the dust grid tree ..."
                        break
                    elif "Adding dust population" in line:
                        self.extra = "Setting dust populations ..."
                        break
                    elif line.startswith("Reading"):
                        self.extra = "Reading input ..."
                        break

            # Get the progress of the stellar emission
            elif self.phase == "stellar":

                # Set simulation phase
                self.simulation_phase = "STELLAR EMISSION"

                # Loop over the lines in reversed order
                for line in reversed(lines):

                    if "Launched stellar emission photon packages" in line:
                        self.progress = float(line.split("packages: ")[1].split("%")[0])
                        break

                # Set initial value for progress
                else: self.progress = 0

            # Get the progress of the spectra calculation phase
            elif self.phase == "spectra":

                nprocesses = get_nprocesses(lines)
                total_entries = None
                entries_per_process = None

                # Loop over the lines
                for line in lines[start_index-4:]:

                    # Check whether dust self-absorption or dust emission phase
                    if "Starting the dust self-absorption phase" in line: self.simulation_phase = "DUST SELF-ABSORPTION"
                    elif "Starting the dust emission phase" in line: self.simulation_phase = "DUST EMISSION"

                    # If this is the log message that marks the very start of the spectra calculation, record the associated time
                    # If this log message states the total number of library entries that are used, record this number
                    if "Library entries in use" in line:

                        #spectra_start = log_file.contents["Time"][i]

                        # Get the total number of library entries in use and the number of entries per process
                        total_entries = int(line.split("use: ")[1].split(" out of")[0])
                        entries_per_process = total_entries / nprocesses

                        # Initial value of progress
                        self.progress = 0.

                    # Get the progress
                    elif "Calculating emission for" in line:

                        entry = float(line.split()[-1][:-3])

                        # Determine the progress
                        # if self.staggered: fraction = entry / total_entries
                        # else: fraction = (entry - process * entries_per_process) / entries_per_process

                        fraction = entry / total_entries
                        progress = float(fraction) * 100.

                        # Set the progress
                        self.progress = progress

            # Get the progress of a dust emission cycle
            elif self.phase == "dust":

                # Loop over the lines in reversed order
                for line in reversed(lines):

                    if "Launched dust emission photon packages" in line:
                        self.progress = float(line.split("packages: ")[1].split("%")[0])
                        self.simulation_phase = "DUST EMISSION"
                        break

                    elif "Launched last-stage dust self-absorption cycle" in line:
                        self.cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
                        self.stage = 2
                        self.progress = float(line.split("packages: ")[1].split("%")[0])
                        self.simulation_phase = "DUST SELF-ABSORPTION"
                        break

                    elif "Starting the last-stage dust self-absorption cycle" in line:
                        self.cycle = int(line.split("cycle ")[1].split("...")[0])
                        self.stage = 2
                        self.progress = 0.
                        self.simulation_phase = "DUST SELF-ABSORPTION"
                        break

                    elif "Launched second-stage dust self-absorption cycle" in line:
                        self.cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
                        self.stage = 1
                        self.progress = float(line.split("packages: ")[1].split("%")[0])
                        self.simulation_phase = "DUST SELF-ABSORPTION"
                        break

                    elif "Starting the second-stage dust self-absorption cycle" in line:
                        self.cycle = int(line.split("cycle ")[1].split("...")[0])
                        self.stage = 1
                        self.progress = 0.
                        self.simulation_phase = "DUST SELF-ABSORPTION"
                        break

                    elif "Launched first-stage dust self-absorption cycle" in line:
                        self.cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
                        self.stage = 0
                        self.progress = float(line.split("packages: ")[1].split("%")[0])
                        self.simulation_phase = "DUST SELF-ABSORPTION"
                        break

                    elif "Starting the first-stage dust self-absorption cycle" in line:
                        self.cycle = int(line.split("cycle ")[1].split("...")[0])
                        self.stage = 0
                        self.progress = 0.
                        self.simulation_phase = "DUST SELF-ABSORPTION"
                        break

                # No break encountered, set initial value for progress
                else: self.progress = 0.

# -----------------------------------------------------------------
