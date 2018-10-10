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
import sys
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
from ..tools import formatting as fmt
from ..tools.sequences import contains_any
from ..tools.stringify import tostr

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

default_similarity_threshold = 0.90
#new_default_similarity_threshold = 0.73 # decreased till similarity between 'Computing density for cell' messages was reached
new_default_similarity_threshold = 0.6 # even more decreased for computing density for cell limit?

# -----------------------------------------------------------------

increase_similarity_after = 5
increase_similarity_factor = 10

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

        # The simulation phase
        self.simulation_phase = None

        # The phase and previous phase
        self.phase = None
        self.previous_phase = None

        # The stage (if applicable)
        self.stage = None
        self.previous_stage = None

        # The cycle (if applicable)
        self.cycle = None
        self.previous_cycle = None

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
        return len(self.log_lines)

    # -----------------------------------------------------------------

    @property
    def last_line(self):
        return self.log_lines[-1]

    # -----------------------------------------------------------------

    @property
    def last_message(self):
        if self.nlines == 0: return None
        else: return get_message(self.last_line)

    # -----------------------------------------------------------------

    @property
    def first_line(self):
        return self.log_lines[0]

    # -----------------------------------------------------------------

    @property
    def first_message(self):
        return get_message(self.first_line)

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

        # Start with status (e.g. running)
        string = self.status

        # Add phase
        if self.phase is not None: string += ": " + phase_descriptions[self.phase]

        # Add simulation phase
        if self.simulation_phase is not None: string += " in the " + self.simulation_phase.lower() + " phase"

        # Add stage, cycle
        if self.stage is not None: string += " stage " + str(self.stage + 1)
        if self.cycle is not None: string += " [cycle " + str(self.cycle) + "]"

        # Add previous phase, stage cycle
        if self.previous_phase is not None and self.phase == "spectra":
            string += " after the " + phase_descriptions[self.previous_phase] + " phase"
            if self.previous_stage is not None and self.previous_cycle is not None: string += " (stage " + str(self.previous_stage) + ", cycle " + str(self.previous_cycle) + ")"
            elif self.previous_stage is not None: string += " (stage " + str(self.previous_stage) + ")"

        # Add progress
        if self.progress is not None: string += " " + tostr(self.progress, round=True, decimal_places=1) + "%"

        # Return the string
        return string

    # -----------------------------------------------------------------

    @property
    def running(self):
        return self.status == "running"

    # -----------------------------------------------------------------

    @property
    def finished(self):
        return self.status == "finished"

    # -----------------------------------------------------------------

    @property
    def crashed(self):
        return self.status == "crashed"

    # -----------------------------------------------------------------

    @property
    def aborted(self):
        return self.status == "aborted"

    # -----------------------------------------------------------------

    @property
    def not_started(self):
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
        self.current_nsimilar = 0

        # Flag
        self.ignored_previous = False

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
                                #print("END OF STAGE?")
                                bar.show(100) # make sure it always ends on 100%
                                break
                            if self.cycle != last_cycle:
                                #print("END OF CYCLE?")
                                bar.show(100) # make sure it always ends on 100%
                                break
                            if self.progress is None:
                                #print("PROGRESS IS NONE!")
                                bar.show(100)
                            else:
                                #print("PROGRESS: " + str(int(self.progress)))
                                bar.show(int(self.progress))
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
                            #print("END OF PHASE?")
                            bar.show(100) # make sure it always ends on 100%
                            break
                        if self.progress is None:
                            #print("PROGRESS IS NONE!")
                            bar.show(100)
                        else:
                            #print("PROGRESS: " + str(int(self.progress)))
                            bar.show(int(self.progress))
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

    def refresh(self, process_or_handle=None, finish_at=None, finish_after=None, similar_log_frequency=10,
                similarity_threshold=new_default_similarity_threshold):

        """
        This function ...
        :param process_or_handle:
        :param finish_at:
        :param finish_after:
        :param similar_log_frequency:
        :param similarity_threshold:
        :return:
        """

        # Set previous
        self.previous_phase = self.phase
        self.previous_stage = self.stage
        self.previous_cycle = self.cycle

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

        # There are new lines and we are not in the middle of a progress bar
        if self.debug_output and len(lines) > self.nlines and self.progress is None:

            # Get current number of columns of the shell
            total_ncolumns = terminal.ncolumns()
            usable_ncolumns = total_ncolumns - 26 - len(skirt_debug_output_prefix) - len(skirt_debug_output_suffix) - ndebug_output_whitespaces
            if usable_ncolumns < 20: usable_ncolumns = 20
            nnew = len(lines) - self.nlines
            previous_message = ""

            for line in lines[-nnew:]:

                message = get_message(line)

                # Current message is similar to the previous message
                similarity = strings.similarity(message, previous_message)
                #print(similarity, similarity_threshold)
                if similarity > similarity_threshold:

                    self.nsimilar += 1
                    self.current_nsimilar += 1
                    #print(self.nsimilar, increase_similarity_after * similar_log_frequency)

                    if self.nsimilar > increase_similarity_after * similar_log_frequency: # there have been 10 messages allowed through already

                        #print(str(increase_similarity_after) + " messages have been allowed to pass through")
                        similar_log_frequency *= increase_similarity_factor

                    if self.current_nsimilar % similar_log_frequency != 0:

                        #sys.stdout.write("\r" + fmt.blue + skirt_debug_output_prefix + "." * self.nsimilar + skirt_debug_output_suffix + fmt.reset)
                        if self.current_nsimilar > usable_ncolumns:
                            ndots = usable_ncolumns
                            nspaces = 0
                        else:
                            ndots = self.current_nsimilar
                            nspaces = usable_ncolumns - ndots
                        sys.stdout.write(fmt.blue + time.timestamp() + " D " + skirt_debug_output_prefix + "." * ndots + " " * nspaces + skirt_debug_output_suffix + fmt.reset + "\r")
                        sys.stdout.flush()
                        continue

                    else: print("")

                else: self.nsimilar = 0
                self.current_nsimilar = 0

                # Show only if not want to be ignored
                if self.ignore_output is not None and contains_any(message, self.ignore_output):
                    self.ignored_previous = True
                    continue
                if self.ignored_previous and message.startswith("  "): continue # ignore sub-messages
                if not message.startswith("  "): self.ignored_previous = False
                message = strings.add_whitespace_or_ellipsis(message, usable_ncolumns, ellipsis_position="center")
                #log.debug(skirt_debug_output_prefix + message + skirt_debug_output_suffix)
                print(fmt.blue + time.timestamp() + " D " + skirt_debug_output_prefix + message + skirt_debug_output_suffix + fmt.reset)
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

            # Get the info of the phase at which the crash happened
            self.phase, self.simulation_phase, self.stage, self.cycle, self.progress, self.extra = get_phase_info(self.log_lines)

            # Return
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
            self.phase, self.simulation_phase, self.stage, self.cycle, self.progress, self.extra = get_phase_info(lines)

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

    def track_progress(self, show=False, finish_at=None, finish_after=None, refresh_frequency=5,
                       similar_log_frequency=10, similarity_threshold=new_default_similarity_threshold):

        """
        This function ...
        :param show:
        :param finish_at:
        :param finish_after:
        :param refresh_frequency:
        :param similar_log_frequency:
        :param similarity_threshold:
        :return:
        """

        break_after_next_line = False

        # Set previous
        self.previous_phase = self.phase
        self.previous_stage = self.stage
        self.previous_cycle = self.cycle

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
        usable_ncolumns = total_ncolumns - 26 - len(skirt_debug_output_prefix) - len(skirt_debug_output_suffix) - ndebug_output_whitespaces
        if usable_ncolumns < 20: usable_ncolumns = 20

        last_phase = None
        last_stage = None
        last_cycle = None
        last_extra = None

        # Loop over the lines
        nsimilar = 0
        current_nsimilar = 0
        ignored_previous = False
        for line in terminal.fetch_lines(self.child):

            # Show the SKIRT line if debug_output is enabled and we are not in the middle of a progress bar
            if self.debug_output and self.progress is None:

                # Get the log message
                message = get_message(line)

                # The current message is similar to the previous message
                if self.last_message is not None and strings.similarity(message, self.last_message) > similarity_threshold:

                    nsimilar += 1
                    current_nsimilar += 1

                    #print(nsimilar, increase_similarity_after * similar_log_frequency)
                    if nsimilar > increase_similarity_after * similar_log_frequency: # there have been 10 messages allowed through already

                        #print(str(increase_similarity_after) + " messages have been allowed to pass through")
                        similar_log_frequency *= increase_similarity_factor

                    # We currently don't want to show the log messages
                    if current_nsimilar % similar_log_frequency != 0:

                        #sys.stdout.write("\r" + fmt.blue + skirt_debug_output_prefix + "." * nsimilar + skirt_debug_output_suffix + fmt.reset)
                        if current_nsimilar > usable_ncolumns:
                            ndots = usable_ncolumns
                            nspaces = 0
                        else:
                            ndots = current_nsimilar
                            nspaces = usable_ncolumns - ndots
                        #total_length = ndots + nspaces
                        #print(nsimilar, current_nsimilar, total_length)
                        if show:
                            sys.stdout.write(fmt.blue + time.timestamp() + " D " + skirt_debug_output_prefix + "." * ndots + " " * nspaces + skirt_debug_output_suffix + fmt.reset + "\r")
                            sys.stdout.flush()
                        #print(ndots)
                        continue

                    else: print("")

                # The current message is not similar to the previous message
                else: nsimilar = 0

                # Debug output is allowed through, reset current nsimilar
                current_nsimilar = 0

                # Show only if not want to be ignored
                if self.ignore_output is not None and contains_any(message, self.ignore_output):
                    ignored_previous = True
                    continue
                if ignored_previous and message.startswith("  "): continue # ignore sub-messages
                if not message.startswith("  "): ignored_previous = False
                #print(list(line))
                message = strings.add_whitespace_or_ellipsis(message, usable_ncolumns, ellipsis_position="center")
                #log.debug(skirt_debug_output_prefix + message + skirt_debug_output_suffix)
                if show: print(fmt.blue + time.timestamp() + " D " + skirt_debug_output_prefix + message + skirt_debug_output_suffix + fmt.reset)

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

                # Get the info
                self.phase, self.simulation_phase, self.stage, self.cycle, self.progress, self.extra = get_phase_info(self.log_lines)

                # Return fail
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

            # New phase
            #elif last_phase is None or self.phase != last_phase:
            if last_phase is None or self.phase != last_phase:

                last_phase = self.phase
                if last_phase is None: continue
                else:
                    if self.simulation_phase is not None:
                        if show: log.info("Starting " + phase_descriptions[last_phase] + " in " + self.simulation_phase.lower() + " phase ...")
                    elif show: log.info("Starting " + phase_descriptions[last_phase] + " ...")

            # Self absorption phase
            if self.phase == "dust" and self.simulation_phase == "DUST SELF-ABSORPTION":

                # print("DUST SELF-ABSORPTION")

                total_length = 100

                if self.stage is None:
                    #self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)
                    continue

                if last_stage is None or self.stage != last_stage:
                    if show: log.info("Starting stage " + str(self.stage + 1) + " ...")
                    last_stage = self.stage

                if self.cycle is None:
                    #self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)
                    continue

                if last_cycle is None or self.cycle != last_cycle:

                    if show: log.info("Starting cycle " + str(self.cycle) + " ...")
                    last_cycle = self.cycle

                    if show:
                        # Progress bar
                        if self._bar is None: self._bar = Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                                 filled_char=BAR_FILLED_CHAR, expected_size=total_length, every=1, add_datetime=True)

                        if self.stage != last_stage:
                            self._bar.show(100) # make sure it always ends on 100%
                            self._bar.__exit__(None, None, None)

                        if self.cycle != last_cycle:
                            self._bar.show(100) # make sure it always ends on 100%
                            self._bar.__exit__(None, None, None)

                        if self.progress is None:
                            self._bar.show(100)
                            #self._bar.__exit__(None, None, None) # ?

                        else: self._bar.show(int(self.progress))
                    #self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)

                    #self.refresh_after(refresh_time, finish_at=finish_at, finish_after=finish_after)
                    continue

                else: continue
                #self.refresh_after(refresh_time, finish_at=finish_at, finish_after=finish_after)

            # Stellar emission, dust emission or spectra calculation: show progress bar
            elif self.phase == "stellar" or self.phase == "spectra" or self.phase == "dust":

                total_length = 100

                if show:
                    # Progress bar
                    if self._bar is None: self._bar = Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                             filled_char=BAR_FILLED_CHAR, expected_size=total_length, every=1,
                             add_datetime=True)
                    # Loop
                    #while True:

                    if self.phase != last_phase:
                        self._bar.show(100)  # make sure it always ends on 100%
                        self._bar.__exit__(None, None, None)

                    if self.progress is None:
                        self._bar.show(100)
                        self._bar.__exit__(None, None, None)

                    else: self._bar.show(int(self.progress))
                #self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)
                continue

            # Same phase
            else:

                if self.extra is not None:
                    if last_extra is None or last_extra != self.extra:
                        if show: log.info(self.extra)
                    last_extra = self.extra

        # We shouldn't get here
        #raise RuntimeError("We shouldn't get here!")

        # Check whether finished
        for line in reversed(self.log_lines):

            message = get_message(line)
            #print("MESS", message)
            if "Finished simulation" in message:
                if show: log.success("Simulation finished")
                self.status = "finished"
                return True

            elif "Available memory" in message:
                if show: log.success("Simulation finished")
                self.status = "finished"
                return True

            elif " *** Error:" in message:
                if show: log.error("Simulation crashed")
                self.status = "crashed"
                return False

        # No break
        else:
            if show: log.error("Simulation has been aborted")
            self.status = "aborted"
            return False

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
            try:
                progress = float(line.split("packages: ")[1].split("%")[0])
                break
            except: pass # SOMETHING WEIRD WITH THE LINE

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

    #nprocesses = get_nprocesses(lines)
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
            #entries_per_process = total_entries / nprocesses

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
            try:
                progress = float(line.split("packages: ")[1].split("%")[0])
                simulation_phase = "DUST EMISSION"
                break
            except: pass  # SOMETHING WEIRD WITH THE LINE

        elif "Launched last-stage dust self-absorption cycle" in line:
            cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
            stage = 2

            try:
                progress = float(line.split("packages: ")[1].split("%")[0])
                simulation_phase = "DUST SELF-ABSORPTION"
                break
            except: pass # SOMETHING WEIRD WITH THE LINE

        elif "Starting the last-stage dust self-absorption cycle" in line:
            cycle = int(line.split("cycle ")[1].split("...")[0])
            stage = 2
            progress = 0.
            simulation_phase = "DUST SELF-ABSORPTION"
            break

        elif "Launched second-stage dust self-absorption cycle" in line:
            cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
            stage = 1
            try:
                progress = float(line.split("packages: ")[1].split("%")[0])
                simulation_phase = "DUST SELF-ABSORPTION"
                break
            except: pass  # SOMETHING WEIRD WITH THE LINE

        elif "Starting the second-stage dust self-absorption cycle" in line:
            cycle = int(line.split("cycle ")[1].split("...")[0])
            stage = 1
            progress = 0.
            simulation_phase = "DUST SELF-ABSORPTION"
            break

        elif "Launched first-stage dust self-absorption cycle" in line:
            cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
            stage = 0
            try:
                progress = float(line.split("packages: ")[1].split("%")[0])
                simulation_phase = "DUST SELF-ABSORPTION"
                break
            except: pass  # SOMETHING WEIRD WITH THE LINE

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
    #phase = None
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
    elif phase == "spectra": simulation_phase, stage, cycle, progress, extra = get_spectra_info(lines, start_index)

    # Get the progress of a dust emission cycle
    elif phase == "dust": simulation_phase, stage, cycle, progress, extra = get_dust_info(lines)

    # Return the info
    return phase, simulation_phase, stage, cycle, progress, extra

# -----------------------------------------------------------------

def get_message_type_symbol(line):

    """
    Thisf unction ...
    :param line:
    :return:
    """

    if is_regular_line(line): return line[24]
    else: return None

# -----------------------------------------------------------------

def get_message_type(line):

    """
    This function ...
    :param line:
    :return:
    """

    symbol = get_message_type_symbol(line)
    return get_type_for_symbol(symbol)

# -----------------------------------------------------------------

def get_type_for_symbol(symbol):

    """
    This function ...
    :param symbol:
    :return:
    """

    if symbol is None: return None
    elif symbol == " ": return "info"
    elif symbol == "?": return "question"
    elif symbol == "*": return "error"
    elif symbol == "-": return "success"
    elif symbol == "D": return "debug"
    elif symbol == "!": return "warning"
    else:
        log.warning("Unknown message symbol: '" + symbol + "'")
        return None

# -----------------------------------------------------------------

def get_message(line):

    """
    Thisnf unction ...
    :param line: 
    :return: 
    """

    if is_regular_line(line): return line[26:]
    else: return line

# -----------------------------------------------------------------

def get_type_and_message(line):

    """
    This function ...
    :param line:
    :return:
    """

    if is_regular_line(line):
        symbol = line[24]
        message = line[26:]
        return get_type_for_symbol(symbol), message
    else: return None, line

# -----------------------------------------------------------------

def is_regular_line(line):

    """
    This function ...
    :param line:
    :return:
    """

    if line[2] == "/" and line[5] == "/" and line[13] == ":" and line[16]: return True
    else: return False

# -----------------------------------------------------------------

default_ignore_output = ["Adding dust population", "Grain sizes range from", "Grain composition grid",
                         "Reading heat capacity data", "Reading grain composition", "closed.", "Reading SED data",
                         "Reading FITS file", "Frame dimensions:", "Writing grain size information", "created.",
                         "Writing optical dust population", "Writing dust population masses",
                         "Writing combined dust mix properties", "Reading wavelength grid data",
                         "Number of leaf cells of each level:", "    Level", ("wavelengths from", "micron to"),
                         ("grain sizes from", "micron to"), ("temperatures from", "K to")]

# -----------------------------------------------------------------

def show_log_summary(lines, debug_output=False, ignore_output=default_ignore_output, refresh_frequency=5, similar_log_frequency=10,
                     similarity_threshold=new_default_similarity_threshold):

    """
    This function ...
    :param lines:
    :param debug_output:
    :param ignore_output:
    :param refresh_frequency:
    :param similar_log_frequency:
    :param similarity_threshold:
    :return:
    """

    #break_after_next_line = False

    # The status
    status = None

    # The phase
    phase = None

    # The simulation phase
    simulation_phase = None

    # The stage (if applicable)
    stage = None

    # The cycle (if applicable)
    cycle = None

    # The progress (if applicable)
    progress = None

    # More info
    extra = None

    # Progress bar where needed
    _bar = None

    # Get current number of columns of the shell
    total_ncolumns = terminal.ncolumns()
    usable_ncolumns = total_ncolumns - 26 - len(skirt_debug_output_prefix) - len(
        skirt_debug_output_suffix) - ndebug_output_whitespaces
    if usable_ncolumns < 20: usable_ncolumns = 20

    # Initialize
    last_phase = None
    last_stage = None
    last_cycle = None
    last_extra = None

    last_message = None

    # Loop over the lines
    log_lines = []
    nsimilar = 0
    current_nsimilar = 0
    ignored_previous = False
    for line in lines:

        # Get the log message
        message = get_message(line)

        # Show the SKIRT line if debug_output is enabled and we are not in the middle of a progress bar
        if debug_output and progress is None:

            # The current message is similar to the previous message
            if last_message is not None: similarity = strings.similarity(message, last_message)
            else: similarity = None
            #print(similarity, similarity_threshold)
            if last_message is not None and similarity > similarity_threshold:

                nsimilar += 1
                current_nsimilar += 1

                # print(nsimilar, increase_similarity_after * similar_log_frequency)
                if nsimilar > increase_similarity_after * similar_log_frequency:  # there have been 10 messages allowed through already

                    # print(str(increase_similarity_after) + " messages have been allowed to pass through")
                    similar_log_frequency *= increase_similarity_factor

                # We currently don't want to show the log messages
                if current_nsimilar % similar_log_frequency != 0:

                    # sys.stdout.write("\r" + fmt.blue + skirt_debug_output_prefix + "." * nsimilar + skirt_debug_output_suffix + fmt.reset)
                    if current_nsimilar > usable_ncolumns:
                        ndots = usable_ncolumns
                        nspaces = 0
                    else:
                        ndots = current_nsimilar
                        nspaces = usable_ncolumns - ndots
                    total_length = ndots + nspaces
                    # print(nsimilar, current_nsimilar, total_length)
                    sys.stdout.write(fmt.blue + time.timestamp() + " D " + skirt_debug_output_prefix + "." * ndots + " " * nspaces + skirt_debug_output_suffix + fmt.reset + "\r")
                    sys.stdout.flush()
                    # print(ndots)
                    continue

                else: print("")

            # The current message is not similar to the previous message
            else: nsimilar = 0

            # Debug output is allowed through, reset current nsimilar
            current_nsimilar = 0

            # Show only if not want to be ignored
            if ignore_output is not None and contains_any(message, ignore_output):
                ignored_previous = True
                continue
            if ignored_previous and message.startswith("  "): continue  # ignore sub-messages
            if not message.startswith("  "): ignored_previous = False
            # print(list(line))
            message = strings.add_whitespace_or_ellipsis(message, usable_ncolumns, ellipsis_position="center")
            # log.debug(skirt_debug_output_prefix + message + skirt_debug_output_suffix)
            print(fmt.blue + time.timestamp() + " D " + skirt_debug_output_prefix + message + skirt_debug_output_suffix + fmt.reset)

        # Add the line
        log_lines.append(line)
        last_message = message
        #nlines += 1
        nlines = len(log_lines)

        # CHECK WHETHER WE HAVE TO BREAK
        # if break_after_next_line:
        #
        #     status = "finished"
        #     # break
        #     # return
        #     log.success("Simulation finished")
        #     return True

        # CONTINUE 'JUST EXPECTING' AND FILLING LOG_LINES UNTIL WE ACTUALLY HAVE TO DO THE WORK OF CHECKING THE PHASE AND STUFF
        if nlines % refresh_frequency != 0: continue

        # Interpret the content of the last line
        if " Finished simulation " in line:

            # print("FOUND FINISHED SIMULATION IN LAST LOG LINE [" + last + "]")
            status = "finished"
            # return
            # break
            log.success("Simulation finished")
            return True

        elif " *** Error: " in line:

            # print("FOUND ERROR IN LAST LOG LINE [" + last + "]")
            # print("LINES:")
            # for line in lines: print("  " + line)
            status = "crashed"
            # return
            # break
            log.error("Simulation crashed")

            # Get the info
            phase, simulation_phase, stage, cycle, progress, extra = get_phase_info(log_lines)

            # Return fail
            return False

        # Running
        else:

            # Check if the child is alive
            #if not self.child.isalive():
            #    self.status = "aborted"
            #    # break
            #    log.error("Simulation has been aborted")
            #    return False

            # Running
            status = "running"

            # Get the info
            phase, simulation_phase, stage, cycle, progress, extra = get_phase_info(log_lines)

        # New phase
        # elif last_phase is None or self.phase != last_phase:
        if last_phase is None or phase != last_phase:

            last_phase = phase
            if last_phase is None: continue
            else:
                if simulation_phase is not None: log.info("Starting " + phase_descriptions[last_phase] + " in " + simulation_phase.lower() + " phase ...")
                else: log.info("Starting " + phase_descriptions[last_phase] + " ...")

        # Self absorption phase
        if phase == "dust" and simulation_phase == "DUST SELF-ABSORPTION":

            # print("DUST SELF-ABSORPTION")

            total_length = 100

            if stage is None:
                # self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)
                continue

            if last_stage is None or stage != last_stage:
                log.info("Starting stage " + str(stage + 1) + " ...")
                last_stage = stage

            if cycle is None:
                # self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)
                continue

            if last_cycle is None or cycle != last_cycle:

                log.info("Starting cycle " + str(cycle) + " ...")
                last_cycle = cycle

                # Progress bar
                if _bar is None: _bar = Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                                                      filled_char=BAR_FILLED_CHAR, expected_size=total_length, every=1,
                                                      add_datetime=True)

                if stage != last_stage:
                    _bar.show(100)  # make sure it always ends on 100%
                    _bar.__exit__(None, None, None)

                if cycle != last_cycle:
                    _bar.show(100)  # make sure it always ends on 100%
                    _bar.__exit__(None, None, None)

                if progress is None:
                    _bar.show(100)
                    # self._bar.__exit__(None, None, None) # ?

                else: _bar.show(int(progress))
                # self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)

                # self.refresh_after(refresh_time, finish_at=finish_at, finish_after=finish_after)
                continue

            else: continue
            # self.refresh_after(refresh_time, finish_at=finish_at, finish_after=finish_after)

        # Stellar emission, dust emission or spectra calculation: show progress bar
        elif phase == "stellar" or phase == "spectra" or phase == "dust":

            total_length = 100

            # Progress bar
            if _bar is None: _bar = Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                                          filled_char=BAR_FILLED_CHAR, expected_size=total_length, every=1,
                                          add_datetime=True)
            # Loop
            # while True:

            if phase != last_phase:
                _bar.show(100)  # make sure it always ends on 100%
                _bar.__exit__(None, None, None)

            if progress is None:
                _bar.show(100)
                _bar.__exit__(None, None, None)

            else: _bar.show(int(progress))
            # self.refresh_after(1, finish_at=finish_at, finish_after=finish_after)
            continue

        # Same phase
        else:

            if extra is not None:
                if last_extra is None or last_extra != extra: log.info(extra)
                last_extra = extra

# -----------------------------------------------------------------

def get_status_from_last_log_line(line, nlibrary_entries=None):

    """
    This function ...
    :param line:
    :param nlibrary_entries:
    :return:
    """

    # Stellar emission
    if "Launched stellar emission photon packages" in line:

        progress = float(line.split("packages:")[1].split("%")[0])
        return "running: emission of stellar photons in the stellar emission phase " + tostr(progress, round=True, decimal_places=1) + "%"

    # Emission spectra
    elif "Calculating emission for library entry" in line:

        if nlibrary_entries is not None:

            entry = int(line.split("for library entry")[1].split("...")[0])
            progress = float(entry)/nlibrary_entries * 100.
            return "running: calculation of dust emission spectra " + tostr(progress, round=True, decimal_places=1) + "%"

        else: return "running: calculation of dust emission spectra"

    # Dust self-absorption
    elif "Launched" in line and "dust self-absorption cycle" in line and "photon packages" in line:

        # Get stage
        stage_description = line.split("Launched")[1].split("dust self-absorption")[0].strip()
        if stage_description == "first-stage": stage = 1
        elif stage_description == "second-stage": stage = 2
        elif stage_description == "last-stage": stage = 3
        else: raise ValueError("Invalid dust self-absorption stage: " + stage_description)

        # Get cycle and progress
        cycle = int(line.split("cycle")[1].split("photon packages")[0])
        progress = float(line.split("packages:")[1].split("%")[0])

        # Return status
        return "running: emission of dust photons in the dust self-absorption phase stage " + str(stage) + " [cycle " + str(cycle) + "] " + tostr(progress, round=True, decimal_places=1) + "%"

    # Cannot interpret status
    else: return None

# -----------------------------------------------------------------
