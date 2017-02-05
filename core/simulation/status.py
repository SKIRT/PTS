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

        # The cycle (if applicable)
        self.cycle = None

        # The progress (if applicable)
        self.progress = None

        # Get the status
        self.refresh()

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        if self.phase is None: return self.status
        elif "self-absorption" in self.phase: return "running: " + str(self.phase) + ", cycle " + str(self.cycle) + "] " + str(self.progress) + "%"
        elif "stellar emission" in self.phase or "dust emission" in self.phase: return "running: " + str(self.phase) + " " + str(self.progress) + "%"
        else: return "running: " + str(self.phase)

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
    def not_started(self):

        """
        This function ...
        :return:
        """

        return self.status == "not started"

    # -----------------------------------------------------------------

    def show_progress(self):

        """
        This function ...
        :return:
        """

        last_phase = None
        last_cycle = None
        while True:

            if self.not_started:

                log.info("Waiting for simulation to start ...")
                time.wait(5)
                self.refresh()
                continue

            if self.finished:
                log.success("Simulation finished")
                break

            elif self.crashed:
                log.success("Simulation crashed")
                break

            elif last_phase is None or self.phase != last_phase:
                last_phase = self.phase
                if last_phase is None:
                    time.wait(5)
                    self.refresh()
                    continue
                else: log.info("Starting " + last_phase + " ...")

            # Stellar emission: show progress bar
            if self.phase == "stellar emission":

                total_length = 100

                # Progress bar
                with Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                         filled_char=BAR_FILLED_CHAR, expected_size=total_length, every=1) as bar:
                    # Loop
                    while True:
                        if self.phase != last_phase: break
                        if self.progress is None:
                            bar.show(100)
                        else:
                            bar.show(int(self.progress))
                        time.wait(1)
                        self.refresh()

            elif "self-absorption" in self.phase:

                total_length = 100

                if self.cycle is None:
                    time.wait(1)
                    self.refresh()
                    continue

                if last_cycle is None or self.cycle != last_cycle:

                    log.info("Starting cycle " + str(self.cycle) + " ...")
                    last_cycle = self.cycle

                    # Progress bar
                    with Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                             filled_char=BAR_FILLED_CHAR, expected_size=total_length, every=1) as bar:

                        # Loop
                        while True:

                            if self.cycle != last_cycle: break
                            if self.progress is None:
                                bar.show(100)
                            else:
                                bar.show(int(self.progress))
                            time.wait(1)
                            self.refresh()

            else:

                # Wait for 5 seconds
                time.wait(5)

                # Refresh
                self.refresh()

    # -----------------------------------------------------------------

    def refresh(self):

        """
        This function ...
        :return:
        """

        # Clear everything
        self.status = None
        self.phase = None
        self.cycle = None
        self.progress = None

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

            self.status = "running"

            # Loop over the log lines
            for line in lines:

                if "Starting setup" in line:

                    self.phase = "setup"
                    self.cycle = None
                    self.progress = 0.

                elif "Starting the stellar emission phase" in line:

                    self.phase = "stellar emission"
                    self.cycle = None
                    self.progress = 0.

                elif "Launched stellar emission photon packages" in line:

                    self.cycle = None
                    self.progress = float(line.split("packages: ")[1].split("%")[0])

                elif "Starting the first-stage dust self-absorption cycle" in line:

                    self.phase = "self-absorption stage 1"
                    self.progress = 0.

                elif "Launched first-stage dust self-absorption cycle" in line:

                    self.cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
                    self.progress = float(line.split("packages: ")[1].split("%")[0])

                elif "Starting the second-stage dust self-absorption cycle" in line:

                    self.phase = "self-absorption stage 2"
                    self.progress = 0.

                elif "Launched second-stage dust self-absorption cycle" in line:

                    self.cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
                    self.progress = float(line.split("packages: ")[1].split("%")[0])

                elif "Starting the last-stage dust self-absorption cycle" in line:

                    self.phase = "self-absorption stage 3"
                    self.progress = 0.

                elif "Launched last-stage dust self-absorption cycle" in line:

                    self.cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
                    self.progress = float(line.split("packages: ")[1].split("%")[0])

                elif "Starting the dust emission phase" in line:

                    self.phase = "dust emission"
                    self.cycle = None
                    self.progress= 0.

                elif "Launched dust emission photon packages" in line:

                    self.progress = float(line.split("packages: ")[1].split("%")[0])

                elif "Starting writing results" in line:

                    self.phase = "writing"
                    self.progress = 0.

# -----------------------------------------------------------------
