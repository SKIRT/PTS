#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.status Contains the SimulationStatus class and derived classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod
from subprocess import Popen

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..tools import time
from ..tools.progress import Bar, BAR_FILLED_CHAR, BAR_EMPTY_CHAR
from ..basics.log import log
from .logfile import get_last_phase, get_nprocesses
from ..basics.handle import ExecutionHandle
from ..tools import terminal
from ..tools import strings

# -----------------------------------------------------------------

phase_descriptions = dict()
phase_descriptions["start"] = "the simulation"
phase_descriptions["setup"] = "setup of the simulation"
phase_descriptions["stellar"] = "emission of stellar photons"
phase_descriptions["spectra"] = "calculation of dust emission spectra"
phase_descriptions["dust"] = "emission of dust photons"
phase_descriptions["comm"] = "communication"
phase_descriptions["write"] = "writing results"
phase_descriptions["wait"] = "waiting for synchronization"

# -----------------------------------------------------------------

skirt_debug_output_prefix = "SKIRT :: [[ "
skirt_debug_output_suffix = " ]]"
ndebug_output_whitespaces = 55

# -----------------------------------------------------------------

similarity_threshold = 0.90

# -----------------------------------------------------------------

class SimulationStatus(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, debug_output=False, ignore_output=None):

        """
        This function ...
        :return:
        """

        # Set attributes
        self.debug_output = debug_output
        self.ignore_output = ignore_output

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

        ##

        # The log lines
        self.log_lines = []

    # -----------------------------------------------------------------

    @property
    def nlines(self):

        """
        This function ...
        :return:
        """

        return len(self.log_lines)

    # -----------------------------------------------------------------

    @property
    def last_line(self):

        """
        This function ...
        :return:
        """

        return self.log_lines[-1]

    # -----------------------------------------------------------------

    @property
    def last_message(self):

        """
        This function ...
        :return:
        """

        if self.nlines == 0: return None
        else: return self.last_line[26:]

    # -----------------------------------------------------------------

    @property
    def first_line(self):

        """
        This function ...
        :return:
        """

        return self.log_lines[0]

    # -----------------------------------------------------------------

    @property
    def first_message(self):

        """
        This function ...
        :return:
        """

        return self.first_line[26:]

    # -----------------------------------------------------------------

    @abstractmethod
    def show_progress(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        pass

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

class LogSimulationStatus(SimulationStatus):

    """
    This class ...
    """

    def __init__(self, log_path, remote=None, debug_output=False, ignore_output=None):

        """
        The constructor ...
        :param log_path:
        :param remote:
        :param debug_output:
        :param ignore_output:
        """

        # Call the constructor of the base class
        super(LogSimulationStatus, self).__init__(debug_output=debug_output, ignore_output=ignore_output)

        # Set attributes
        self.log_path = log_path
        self.remote = remote

        # Number of consecutive similar log lines that are encountered
        self.nsimilar = 0

        # Refresh the status
        self.refresh()

    # -----------------------------------------------------------------

    def show_progress(self, process_or_handle, refresh_time=3, finish_at=None, finish_after=None):

        """
        This function ...
        :param process_or_handle:
        :param refresh_time: time in seconds
        :param finish_at:
        :param finish_after:
        :return:
        """

        last_phase = None
        last_stage = None
        last_cycle = None
        last_extra = None

        # Wait for a bit ?

        #print(process_or_handle, type(process_or_handle))
        #print(process_or_handle.poll, process_or_handle.poll())

        # Refresh loop
        while True:

            # Not yet started
            if self.not_started:

                #print("NOT STARTED")

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

                # Subprocess process
                elif isinstance(process_or_handle, Popen):

                    returncode = process_or_handle.poll()
                    if returncode is not None:
                        log.error("The simulation has stopped")
                        if fs.is_file(self.log_path):
                            lines = list(fs.read_lines(self.log_path))
                            for line in lines: log.error(line)
                        else: log.error("Log file hasn't been created yet")
                        return False

                # # Pexpect spawn object?
                # elif introspection.lazy_isinstance(process_or_handle, "spawn", "pexpect", return_false_if_fail=True):
                #
                #     from ..tools.terminal import fetch_lines
                #     if not process_or_handle.isalive():
                #         log.error("The simulation has stopped")
                #         # Read the lines
                #         for line in fetch_lines(process_or_handle): log.error(line)
                #         return False

                # Not recognized
                else:
                    raise ValueError("Cannot interpret the process or execution handle that is passed")
                    #log.error("Cannot interpret the process or execution handle that is passed")

                log.info("Waiting for simulation to start ...")
                time.wait(refresh_time)
                self.refresh(process_or_handle, finish_at=finish_at, finish_after=finish_after)
                continue

            # Finished: break loop
            if self.finished:
                log.success("Simulation finished")
                return True

            # Crashed: break loop
            elif self.crashed:

                # Check whether FINISH POINT HAS BEEN REACHED
                for i, line in enumerate(reversed(self.log_lines)):

                    if finish_at is not None and finish_at in line:

                        log.success("Reached finish_at point")
                        if process_or_handle is not None and not isinstance(process_or_handle, ExecutionHandle):
                            terminal.kill(process_or_handle.pid)
                        return True

                    if finish_after is not None and finish_after in line and i != 0:

                        log.success("Reached finish_after point")
                        if process_or_handle is not None and not isinstance(process_or_handle, ExecutionHandle):
                            terminal.kill(process_or_handle.pid)
                        return True

                else: # loop finished without break or return
                    log.error("Simulation crashed")
                    return False

            # Aborted: break loop
            elif self.aborted:
                log.error("Simulation has been aborted")
                return False

            # New phase
            elif last_phase is None or self.phase != last_phase:

                #print("NEW PHASE")

                last_phase = self.phase
                if last_phase is None:
                    time.wait(refresh_time)
                    self.refresh(process_or_handle, finish_at=finish_at, finish_after=finish_after)
                    continue
                else:
                    if self.simulation_phase is not None: log.info("Starting " + phase_descriptions[last_phase] + " in " + self.simulation_phase.lower() + " phase ...")
                    else: log.info("Starting " + phase_descriptions[last_phase] + " ...")

            # Self absorption phase
            if self.phase == "dust" and self.simulation_phase == "DUST SELF-ABSORPTION":

                #print("DUST SELF-ABSORPTION")

                total_length = 100

                if self.stage is None:
                    self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)
                    continue

                if last_stage is None or self.stage != last_stage:
                    log.info("Starting stage " + str(self.stage + 1) + " ...")
                    last_stage = self.stage

                if self.cycle is None:
                    self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)
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
                            self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)

                    self.refresh_after(refresh_time, finish_at=finish_at, finish_after=finish_after)

                else: self.refresh_after(refresh_time, finish_at=finish_at, finish_after=finish_after)

            # Stellar emission: show progress bar
            elif self.phase == "stellar" or self.phase == "spectra" or self.phase == "dust":

                #print("STELLAR AND DUST EMISSION, SPECTRA CALCULATION")

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
                        self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)

            # Still the same phase
            else:

                #print("SAME PHASE")

                if self.extra is not None:
                    if last_extra is None or last_extra != self.extra:
                        log.info(self.extra)
                    last_extra = self.extra
                self.refresh_after(refresh_time, finish_at=finish_at, finish_after=finish_after)

    # -----------------------------------------------------------------

    def refresh_after(self, seconds, process_or_handle=None, finish_after=None, finish_at=None):

        """
        This function ...
        :param seconds:
        :param process_or_handle:
        :param finish_after:
        :param finish_at:
        :return:
        """

        # Wait and refresh
        time.wait(seconds)
        self.refresh(process_or_handle, finish_after=finish_after, finish_at=finish_at)

    # -----------------------------------------------------------------

    def refresh(self, process_or_handle=None, finish_at=None, finish_after=None, similar_log_frequency=10):

        """
        This function ...
        :param process_or_handle:
        :param finish_at:
        :param finish_after:
        :param similar_log_frequency:
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

        #if present: print("LOG FILE IS PRESENT")
        #else: print("LOG FILE IS NOT (YET) PRESENT")

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

        # There are new lines
        #print(len(lines), self.nlines)
        if self.debug_output and len(lines) > self.nlines:

            # Get current number of columns of the shell
            total_ncolumns = terminal.ncolumns()
            usable_ncolumns = total_ncolumns - 26 - len(skirt_debug_output_prefix) - len(skirt_debug_output_suffix) - ndebug_output_whitespaces
            if usable_ncolumns < 20: usable_ncolumns = 20
            nnew = len(lines) - self.nlines
            previous_message = ""
            for line in lines[-nnew:]:
                message = line[26:]
                if strings.similarity(message, previous_message) > similarity_threshold:
                    self.nsimilar += 1
                    if self.nsimilar > 10 * similar_log_frequency: # there have been 10 messages allowed through already
                        similar_log_frequency *= 10
                    if self.nsimilar % similar_log_frequency != 0: continue
                self.nsimilar = 0

                # Show only if not want to be ignored
                if self.ignore_output is not None and contains_any(message, self.ignore_output): continue
                message = strings.add_whitespace_or_ellipsis(message, usable_ncolumns, ellipsis_position="center")
                log.debug(skirt_debug_output_prefix + message + skirt_debug_output_suffix)
                previous_message = message

        # SET THE LOG LINES
        self.log_lines = lines

        # Check whether finished
        last_two_lines = lines[-2:]
        if len(lines) == 1: last = lines[0]
        elif " Available memory: " in last_two_lines[1]: last = last_two_lines[0]
        else: last = last_two_lines[1]

        #print("LAST LINE", last)

        #print("LINES", lines)

        #print("FINISH AFTER", finish_after)

        # NEW: CHECK WHETHER USER WANTS TO FINISH AFTER A SPECIFIC LINE
        #if "finish_after" in kwargs:
        if finish_after is not None:
            #finish_after_line = kwargs["finish_after"]
            #is_last = finish_after_line in lines[-1]
            for i, line in enumerate(reversed(lines)):
                if i != 0 and finish_after in line: # if it is not the last line (i = 0), meaning there is at least one line after it
                    self.status = "finished"
                    # KILL THE PROCESS
                    if process_or_handle is not None and not isinstance(process_or_handle, ExecutionHandle):
                        terminal.kill(process_or_handle.pid) # KILL
                    return

        #print("FINISH AT", finish_at)

        #if "finish_at" in kwargs:
        if finish_at is not None:
            #finish_at_line = kwargs["finish_at"]
            for line in reversed(lines):
                if finish_at in line:
                    self.status = "finished"
                    # KILL THE PROCESS
                    if process_or_handle is not None and not isinstance(process_or_handle, ExecutionHandle):
                        terminal.kill(process_or_handle.pid) # KILL
                    return

        # Interpret the content of the last line
        if " Finished simulation " in last:
            #print("FOUND FINISHED SIMULATION IN LAST LOG LINE [" + last + "]")
            self.status = "finished"
            return

        elif " *** Error: " in last:
            #print("FOUND ERROR IN LAST LOG LINE [" + last + "]")
            #print("LINES:")
            #for line in lines: print("  " + line)
            self.status = "crashed"
            return

        # Running
        else:

            # Process or handle is passed: check whether the process or screen is still running and not aborted
            if process_or_handle is not None:

                # Execution handle
                if isinstance(process_or_handle, ExecutionHandle):

                    if process_or_handle.type != "screen": raise NotImplementedError("Execution handle must be 'screen' (here it's '" + str(process_or_handle.type) + "')")
                    screen_name = process_or_handle.value

                    if not self.remote.is_active_screen(screen_name):
                        self.status = "aborted"
                        return

                # Process
                elif isinstance(process_or_handle, Popen):

                    returncode = process_or_handle.poll()
                    if returncode is not None:
                        #log.error("The simulation has been aborted")
                        #return False
                        self.status = "aborted"
                        return
                    else: # Still running, but is it hanging because of a problem?
                        pass
                        #out, err = process_or_handle.communicate()
                        #print(out)
                        #print(err)
                        # e.g. a possible error on which it hangs can be:
                        # A system call failed during shared memory initialization that should
                        # not have.  It is likely that your MPI job will now either abort or
                        # experience performance degradation.
                        # Local host:  druif.ugent.be
                        # System call: ftruncate(2)
                        # Error:       No space left on device (errno 28)

                # # Pexpect spawn object?
                # elif introspection.lazy_isinstance(process_or_handle, "spawn", "pexpect", return_false_if_fail=True):
                #     if not process_or_handle.isalive():
                #         self.status = "aborted"
                #         return

                # Not recognized
                else:
                    #raise ValueError("Cannot interpret the process or execution handle that is passed")
                    log.error("Cannot interpret the process or execution handle that is passed")

            # Status is 'running'
            self.status = "running"

            # Get the phase info
            self.phase, self.simulation_phase, stage, cycle, progress, extra = get_phase_info(lines)

# -----------------------------------------------------------------

class SpawnSimulationStatus(SimulationStatus):

    """
    This function ...
    """

    def __init__(self, child, debug_output=False, ignore_output=None):

        """
        This function ...
        :param child:
        :param debug_output:
        :param ignore_output:
        """

        # Call the constructor of the base class
        super(SpawnSimulationStatus, self).__init__(debug_output=debug_output, ignore_output=ignore_output)

        # Set the child
        self.child = child

        # Progress bar where needed
        self._bar = None

    # -----------------------------------------------------------------

    def show_progress(self, finish_at=None, finish_after=None, refresh_frequency=5, similar_log_frequency=10):

        """
        This function ...
        :param finish_at:
        :param finish_after:
        :param refresh_frequency: number lines after which to refresh
        :param similar_log_frequency:
        :return:
        """

        # Track the progress and show
        return self.track_progress(show=True, finish_at=finish_at, finish_after=finish_after,
                                   refresh_frequency=refresh_frequency, similar_log_frequency=similar_log_frequency)

    # -----------------------------------------------------------------

    def track_progress(self, show=False, finish_at=None, finish_after=None, refresh_frequency=5, similar_log_frequency=10):

        """
        This function ...
        :param show:
        :param finish_at:
        :param finish_after:
        :param refresh_frequency:
        :param similar_log_frequency:
        :return:
        """

        break_after_next_line = False

        # Clear everything
        self.status = None
        self.phase = None
        self.simulation_phase = None
        self.stage = None
        self.cycle = None
        self.progress = None
        self.extra = None

        # Get current number of columns of the shell
        total_ncolumns = terminal.ncolumns()
        usable_ncolumns = total_ncolumns - 26 - len(skirt_debug_output_prefix) - len(
            skirt_debug_output_suffix) - ndebug_output_whitespaces
        if usable_ncolumns < 20: usable_ncolumns = 20

        last_phase = None
        last_stage = None
        last_cycle = None
        last_extra = None

        # Loop over the lines
        nsimilar = 0
        for line in terminal.fetch_lines(self.child):

            # Show the SKIRT line if debug_output is enabled
            if self.debug_output:
                message = line[26:]
                if self.last_message is not None and strings.similarity(message, self.last_message) > similarity_threshold:
                    nsimilar += 1
                    if nsimilar > 10 * similar_log_frequency: # there have been 10 messages allowed through already
                        similar_log_frequency *= 10
                    if nsimilar % similar_log_frequency != 0: continue
                nsimilar = 0
                # Show only if not want to be ignored
                if self.ignore_output is not None and contains_any(message, self.ignore_output): continue
                message = strings.add_whitespace_or_ellipsis(message, usable_ncolumns, ellipsis_position="center")
                log.debug(skirt_debug_output_prefix + message + skirt_debug_output_suffix)

            # Add the line
            self.log_lines.append(line)

            # CHECK WHETHER WE HAVE TO BREAK
            if break_after_next_line:
                self.status = "finished"
                #break
                #return
                log.success("Simulation finished")
                return True

            # NEW: CHECK WHETHER USER WANTS TO FINISH AFTER A SPECIFIC LINE
            #if "finish_after" in kwargs:
            if finish_after is not None:
                #finish_after_line = kwargs["finish_after"]
                #is_last = finish_after_line in lines[-1]
                #for i, line in enumerate(reversed(lines)):

                #if i != 0 and finish_after in line: # if it is not the last line (i = 0), meaning there is at least one line after it
                if finish_after in line:
                    break_after_next_line = True
                    continue

                    #self.status = "finished"
                    # KILL THE PROCESS
                    #if process_or_handle is not None and not isinstance(process_or_handle, ExecutionHandle):
                    #    terminal.kill(process_or_handle.pid) # KILL
                    #return

            #if "finish_at" in kwargs:
            if finish_at is not None:
                #finish_at_line = kwargs["finish_at"]
                if finish_at in line:
                    self.status = "finished"
                    # KILL THE PROCESS
                    #if process_or_handle is not None and not isinstance(process_or_handle, ExecutionHandle):
                    #    terminal.kill(process_or_handle.pid) # KILL
                    #return
                    #break
                    log.success("Simulation finished")
                    return True

            # CONTINUE 'JUST EXPECTING' AND FILLING LOG_LINES UNTIL WE ACTUALLY HAVE TO DO THE WORK OF CHECKING THE PHASE AND STUFF
            if self.nlines % refresh_frequency != 0: continue

            # Debugging
            #if self.debug_output: log.debug("Refreshing the status ...")

            # Interpret the content of the last line
            if " Finished simulation " in line:
                # print("FOUND FINISHED SIMULATION IN LAST LOG LINE [" + last + "]")
                self.status = "finished"
                #return
                #break
                log.success("Simulation finished")
                return True

            elif " *** Error: " in line:
                # print("FOUND ERROR IN LAST LOG LINE [" + last + "]")
                # print("LINES:")
                # for line in lines: print("  " + line)
                self.status = "crashed"
                #return
                #break
                log.error("Simulation crashed")
                return False

            # Running
            else:

                # Check if the child is alive
                if not self.child.isalive():
                    self.status = "aborted"
                    #break
                    log.error("Simulation has been aborted")
                    return False

                # Running
                self.status = "running"

                # Get the info
                self.phase, self.simulation_phase, self.stage, self.cycle, self.progress, self.extra = get_phase_info(self.log_lines)

            #print(last_phase, self.phase)
            #print(last_stage, self.stage)
            #print(last_cycle, self.cycle)
            #print(last_extra, self.extra)

            # # Finished: break loop
            # if self.finished:
            #     log.success("Simulation finished")
            #     return True

            # # Crashed: break loop
            # elif self.crashed:
            #     log.error("Simulation crashed")
            #     return False
            #
            # # Aborted: break loop
            # elif self.aborted:
            #     log.error("Simulation has been aborted")
            #     return False

            # New phase
            #elif last_phase is None or self.phase != last_phase:
            if last_phase is None or self.phase != last_phase:

                last_phase = self.phase
                if last_phase is None: continue
                else:
                    if self.simulation_phase is not None: log.info("Starting " + phase_descriptions[last_phase] + " in " + self.simulation_phase.lower() + " phase ...")
                    else: log.info("Starting " + phase_descriptions[last_phase] + " ...")

            # Self absorption phase
            if self.phase == "dust" and self.simulation_phase == "DUST SELF-ABSORPTION":

                # print("DUST SELF-ABSORPTION")

                total_length = 100

                if self.stage is None:
                    #self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)
                    continue

                if last_stage is None or self.stage != last_stage:
                    log.info("Starting stage " + str(self.stage + 1) + " ...")
                    last_stage = self.stage

                if self.cycle is None:
                    #self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)
                    continue

                if last_cycle is None or self.cycle != last_cycle:

                    log.info("Starting cycle " + str(self.cycle) + " ...")
                    last_cycle = self.cycle

                    # Progress bar
                    if self._bar is None: self._bar = Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                             filled_char=BAR_FILLED_CHAR, expected_size=total_length, every=1, add_datetime=True)

                    if self.stage != last_stage:
                        self._bar.show(100) # make sure it always ends on 100%
                        self._bar.__exit__()

                    if self.cycle != last_cycle:
                        self._bar.show(100) # make sure it always ends on 100%
                        self._bar.__exit__()

                    if self.progress is None:
                        self._bar.show(100)
                        #self._bar.__exit__() # ?

                    else: self._bar.show(int(self.progress))
                    #self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)

                    #self.refresh_after(refresh_time, finish_at=finish_at, finish_after=finish_after)
                    continue

                else: continue
                #self.refresh_after(refresh_time, finish_at=finish_at, finish_after=finish_after)

            # Stellar emission: show progress bar
            elif self.phase == "stellar" or self.phase == "spectra" or self.phase == "dust":

                total_length = 100

                # Progress bar
                if self._bar is None: self._bar = Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                         filled_char=BAR_FILLED_CHAR, expected_size=total_length, every=1,
                         add_datetime=True)
                # Loop
                #while True:

                if self.phase != last_phase:
                    self._bar.show(100)  # make sure it always ends on 100%
                    self._bar.__exit__()

                if self.progress is None:
                    self._bar.show(100)
                    self._bar.__exit__()

                else: self._bar.show(int(self.progress))
                #self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)
                continue

            # Same phase
            else:

                if self.extra is not None:
                    if last_extra is None or last_extra != self.extra:
                        log.info(self.extra)
                    last_extra = self.extra

        # We shouldn't get here
        raise RuntimeError("We shouldn't get here!")

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        import pexpect
        self.child.expect(pexpect.EOF)

# -----------------------------------------------------------------

def get_setup_info(lines):

    """
    This function ...
    :param lines:
    :return:
    """

    # Initialize everything
    simulation_phase = None
    stage = None
    cycle = None
    progress = None
    extra = None

    simulation_phase = "SETUP"

    # Loop over the lines in reversed order
    for line in reversed(lines):

        if "Writing dust cell properties" in line:
            extra = "Writing dust cell properties ..."
            break
        elif "Setting the value of the density in the cells" in line:
            extra = "Calculating dust densities ..."
            break
        elif "Writing data to plot the dust grid" in line:
            extra = "Writing the dust grid ..."
            break
        elif "Subdividing node" in line:
            extra = "Creating the dust grid tree ..."
            break
        elif "Starting subdivision" in line:
            extra = "Creating the dust grid tree ..."
            break
        elif "Adding dust population" in line:
            extra = "Setting dust populations ..."
            break
        elif line.startswith("Reading"):
            extra = "Reading input ..."
            break

    # Return the info
    return simulation_phase, stage, cycle, progress, extra

# -----------------------------------------------------------------

def get_stellar_info(lines):

    """
    This function ...
    :param lines:
    :return:
    """

    # Initialize everything
    simulation_phase = None
    stage = None
    cycle = None
    progress = None
    extra = None

    # Set simulation phase
    simulation_phase = "STELLAR EMISSION"

    # Loop over the lines in reversed order
    for line in reversed(lines):

        if "Launched stellar emission photon packages" in line:
            progress = float(line.split("packages: ")[1].split("%")[0])
            break

    # Set initial value for progress
    else: progress = 0

    # Return the info
    return simulation_phase, stage, cycle, progress, extra

# -----------------------------------------------------------------

def get_spectra_info(lines, start_index):

    """
    This function ...
    :param lines:
    :param start_index:
    :return:
    """

    # Initialize everything
    simulation_phase = None
    stage = None
    cycle = None
    progress = None
    extra = None

    nprocesses = get_nprocesses(lines)
    total_entries = None
    entries_per_process = None

    # Loop over the lines
    for line in lines[start_index - 4:]:

        # Check whether dust self-absorption or dust emission phase
        if "Starting the dust self-absorption phase" in line: simulation_phase = "DUST SELF-ABSORPTION"
        elif "Starting the dust emission phase" in line: simulation_phase = "DUST EMISSION"

        # If this is the log message that marks the very start of the spectra calculation, record the associated time
        # If this log message states the total number of library entries that are used, record this number
        if "Library entries in use" in line:

            # spectra_start = log_file.contents["Time"][i]

            # Get the total number of library entries in use and the number of entries per process
            total_entries = int(line.split("use: ")[1].split(" out of")[0])
            entries_per_process = total_entries / nprocesses

            # Initial value of progress
            progress = 0.

        # Get the progress
        elif "Calculating emission for" in line:

            entry = float(line.split()[-1][:-3])

            # Determine the progress
            # if self.staggered: fraction = entry / total_entries
            # else: fraction = (entry - process * entries_per_process) / entries_per_process

            fraction = entry / total_entries
            progress = float(fraction) * 100.

            # Set the progress
            progress = progress

    # Return the info
    return simulation_phase, stage, cycle, progress, extra

# -----------------------------------------------------------------

def get_dust_info(lines):

    """
    This function ...
    :param lines:
    :return:
    """

    # Initialize everything
    simulation_phase = None
    stage = None
    cycle = None
    progress = None
    extra = None

    # Loop over the lines in reversed order
    for line in reversed(lines):

        if "Launched dust emission photon packages" in line:
            progress = float(line.split("packages: ")[1].split("%")[0])
            simulation_phase = "DUST EMISSION"
            break

        elif "Launched last-stage dust self-absorption cycle" in line:
            cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
            stage = 2
            progress = float(line.split("packages: ")[1].split("%")[0])
            simulation_phase = "DUST SELF-ABSORPTION"
            break

        elif "Starting the last-stage dust self-absorption cycle" in line:
            cycle = int(line.split("cycle ")[1].split("...")[0])
            stage = 2
            progress = 0.
            simulation_phase = "DUST SELF-ABSORPTION"
            break

        elif "Launched second-stage dust self-absorption cycle" in line:
            cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
            stage = 1
            progress = float(line.split("packages: ")[1].split("%")[0])
            simulation_phase = "DUST SELF-ABSORPTION"
            break

        elif "Starting the second-stage dust self-absorption cycle" in line:
            cycle = int(line.split("cycle ")[1].split("...")[0])
            stage = 1
            progress = 0.
            simulation_phase = "DUST SELF-ABSORPTION"
            break

        elif "Launched first-stage dust self-absorption cycle" in line:
            cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
            stage = 0
            progress = float(line.split("packages: ")[1].split("%")[0])
            simulation_phase = "DUST SELF-ABSORPTION"
            break

        elif "Starting the first-stage dust self-absorption cycle" in line:
            cycle = int(line.split("cycle ")[1].split("...")[0])
            stage = 0
            progress = 0.
            simulation_phase = "DUST SELF-ABSORPTION"
            break

    # No break encountered, set initial value for progress
    else: progress = 0.

    # Return the info
    return simulation_phase, stage, cycle, progress, extra

# -----------------------------------------------------------------

def get_phase_info(lines):

    """
    This function ...
    :return: 
    """

    # Initialize everything
    phase = None
    simulation_phase = None
    stage = None
    cycle = None
    progress = None
    extra = None

    # Get the last phase in the log file
    last_phase, start_index = get_last_phase(lines)

    # Set the phase
    phase = last_phase

    # The setup phase
    if phase == "setup": simulation_phase, stage, cycle, progress, extra = get_setup_info(lines)

    # Get the progress of the stellar emission
    elif phase == "stellar": simulation_phase, stage, cycle, progress, extra = get_stellar_info(lines)

    # Get the progress of the spectra calculation phase
    elif phase == "spectra": simulation_phase, stage, cycle, progress, extra = get_spectra_info(lines)

    # Get the progress of a dust emission cycle
    elif phase == "dust": simulation_phase, stage, cycle, progress, extra = get_dust_info(lines)

    # Return the info
    return phase, simulation_phase, stage, cycle, progress, extra

# -----------------------------------------------------------------

def contains_any(string, patterns):

    """
    Thisf unction ...
    :param string:
    :param patterns:
    :return:
    """

    for pattern in patterns:
        if pattern in string: return True
    return False

# -----------------------------------------------------------------
