#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.logfile Contains the LogFile class, providing a convenient way of extracting
#  information from a simulation's log file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import warnings

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ..tools import time
from ..tools import filesystem as fs
from ..tools.utils import lazyproperty

# -----------------------------------------------------------------

possible_phases = ["setup", "wait", "comm", "stellar", "spectra", "dust", "write"]

# -----------------------------------------------------------------

class LogFile(object):

    """
    This class ...
    """

    def __init__(self, path=None, contents=None, name=None, verbose_logging=None, memory_logging=None):

        """
        The constructor ...
        :param contents:
        :param name:
        :param verbose_logging:
        :param memory_logging:
        :return:
        """

        # Local path is passed
        if path is not None:

            # Set the log file path
            self.path = path

            # Determine the name of the log file
            name = fs.name(self.path)

            # Determine the simulation prefix
            self.prefix = name.split("_")[0]

            # Determine the process rank associated with this log file
            try: self.process = int(name.split("_logP")[1].split(".txt")[0])
            except IndexError: self.process = 0

            # Parse the log file
            self.contents, self.verbose_logging, self.memory_logging = parse(path, return_options=True)

        # Table is passed
        elif contents is not None:

            # Name must be given
            if name is None: raise ValueError("Name must be specified")

            # No path
            self.path = None

            # Determine the simulation prefix
            self.prefix = name.split("_")[0]

            # Determine the process rank associated with this log file
            try: self.process = int(name.split("_logP")[1].split(".txt")[0])
            except IndexError: self.process = 0

            # Set contents
            self.contents = contents
            self.verbose_logging = verbose_logging
            self.memory_logging = memory_logging

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        return cls(path=path)

    # -----------------------------------------------------------------

    @classmethod
    def from_remote_file(cls, path, remote):

        """
        Thisn function ...
        :param path:
        :param remote:
        :return:
        """

        # Determine the name of the log file
        name = fs.name(path)

        # Parse the contents
        contents, verbose_logging, memory_logging = parse_from_lines(remote.read_lines(path), return_options=True)

        # Create and return
        return cls(contents=contents, name=name, verbose_logging=verbose_logging, memory_logging=memory_logging)

    # -----------------------------------------------------------------

    @lazyproperty
    def messages(self):

        """
        This function ...
        :return:
        """

        return list(self.contents["Message"])

    # -----------------------------------------------------------------

    @lazyproperty
    def host(self):

        """
        This function ...
        :return:
        """

        # Loop over all the log file entries
        for i in range(len(self.contents)):

            # Get the current log message
            message = self.contents["Message"][i]

            if "Running on" in message:

                host = message.split("on ")[1].split(" for")[0]

                if len(host.split(".")) == 3: return host.split(".")[1]
                elif len(host.split(".")) == 2: return host.split(".")[0]
                else: return host

        return None

    # -----------------------------------------------------------------

    @lazyproperty
    def logging_options(self):

        """
        This function ...
        :return:
        """

        from ..launch.options import LoggingOptions
        return LoggingOptions(verbose=self.verbose_logging, memory=self.memory_logging)

    # -----------------------------------------------------------------

    @lazyproperty
    def t_0(self):

        """
        This function ...
        :return:
        """

        # Return the time of the first log message
        return self.contents["Time"][0]

    # -----------------------------------------------------------------

    @lazyproperty
    def t_last(self):

        """
        This function ...
        :return:
        """

        # Return the time of the last log message
        return self.contents["Time"][-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def finished_simulation_message(self):

        """
        This function ...
        :return:
        """

        # Loop over all the log file entries
        for i in range(len(self.contents)):

            # Get the current log message
            message = self.contents["Message"][i]

            # Look for the message that indicates the start of the simulation
            if "Finished simulation" in message: return message

        # Return None if the message could not be found
        return None

    # -----------------------------------------------------------------

    @lazyproperty
    def starting_simulation_message(self):

        """
        This function ...
        :return:
        """

        # Loop over all the log file entries
        for i in range(len(self.contents)):

            # Get the current log message
            message = self.contents["Message"][i]

            # Look for the message that indicates the start of the simulation
            if "Starting simulation" in message: return message

        # Return None if the message could not be found
        return None

    # -----------------------------------------------------------------

    @lazyproperty
    def total_runtime(self):

        """
        This function ...
        :return:
        """

        # Cannot determine total runtime
        if self.finished_simulation_message is None: return None

        seconds = float(self.finished_simulation_message.split(" in ")[-1].split(" s")[0])
        return seconds

    # -----------------------------------------------------------------

    @lazyproperty
    def elapsed_time(self):

        """
        This function ...
        :return:
        """

        seconds = (self.t_last - self.t_0).total_seconds()
        return seconds

    # -----------------------------------------------------------------

    @property
    def has_memory(self):

        """
        This function ...
        :return:
        """

        return "Memory" in self.contents.colnames

    # -----------------------------------------------------------------

    @lazyproperty
    def peak_memory(self):

        """
        This function ...
        :return:
        """

        last_message = self.last_message

        if "Peak memory usage" in last_message: return float(last_message.split("Peak memory usage: ")[1].split(" GB")[0])
        else: raise RuntimeError("The log file was aborted before it could state the peak memory usage")

    # -----------------------------------------------------------------

    @lazyproperty
    def setup_peak_memory(self):

        """
        This function ...
        :return:
        """

        peak_memory = None

        # Loop over all log file entries
        for i in range(len(self.contents)):

            # Skip entries corresponding to other phases
            if self.contents["Phase"][i] != "setup": continue

            # Replace peak memory if appropriate
            if peak_memory is None or self.contents["Memory"][i] > peak_memory: peak_memory = self.contents["Memory"][i]

        # Return the peak memory
        return peak_memory

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_peak_memory(self):

        """
        This function ...
        :return:
        """

        peak_memory = None

        # Loop over all log file entries
        for i in range(len(self.contents)):

            # Skip entries corresponding to other phases
            if self.contents["Phase"][i] != "stellar": continue

            # Replace peak memory if appropriate
            if peak_memory is None or self.contents["Memory"][i] > peak_memory: peak_memory = self.contents["Memory"][i]

        # Return the peak memory
        return peak_memory

    # -----------------------------------------------------------------

    @lazyproperty
    def spectra_peak_memory(self):

        """
        This function ...
        :return:
        """

        peak_memory = None

        # Loop over all log file entries
        for i in range(len(self.contents)):

            # Skip entries corresponding to other phases
            if self.contents["Phase"][i] != "spectra": continue

            # Replace peak memory if appropriate
            if peak_memory is None or self.contents["Memory"][i] > peak_memory: peak_memory = self.contents["Memory"][i]

        # Return the peak memory
        return peak_memory

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_peak_memory(self):

        """
        This function ...
        :return:
        """

        peak_memory = None

        # Loop over all log file entries
        for i in range(len(self.contents)):

            # Skip entries corresponding to other phases
            if self.contents["Phase"][i] != "dust": continue

            # Replace peak memory if appropriate
            if peak_memory is None or self.contents["Memory"][i] > peak_memory: peak_memory = self.contents["Memory"][i]

        # Return the peak memory
        return peak_memory

    # -----------------------------------------------------------------

    @lazyproperty
    def writing_peak_memory(self):

        """
        This function ...
        :return:
        """

        peak_memory = None

        # Loop over all log file entries
        for i in range(len(self.contents)):

            # Skip entries corresponding to other phases
            if self.contents["Phase"][i] != "write": continue

            # Replace peak memory if appropriate
            if peak_memory is None or self.contents["Memory"][i] > peak_memory: peak_memory = self.contents["Memory"][i]

        # Return the peak memory
        return peak_memory

    # -----------------------------------------------------------------

    @lazyproperty
    def last_message(self):

        """
        This function ...
        :return:
        """

        return self.contents["Message"][-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def nchunks(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries
        for i in range(len(self.contents)):

            line = self.contents["Message"][i]

            if line.startswith("Using") and line.endswith("chunks."):
                chunks = int(line.split("Using ")[1].split(" chunks")[0])
                return chunks

        # Number of chunks could not be determined
        return None

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_packages(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Skip entries corresponding to other phases
            if not self.contents["Phase"][i] == "stellar": continue

            # Search for the line stating the number of photon packages
            if "photon packages for each of" in self.contents["Message"][i]:

                #print(self.contents["Message"][i])
                # Return the number of stellar photon packages
                #stellar_packages = int(self.contents["Message"][i].split("(")[1].split(" photon")[0]) # doesn't work anymore?
                stellar_packages = int(self.contents["Message"][i].split(" photon packages")[0])
                return stellar_packages

        # If the number of stellar photon packages could not be determined, return None
        return None

    # -----------------------------------------------------------------

    @property
    def npackages(self):
        return self.stellar_packages

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_packages_per_process(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Skip entries corresponding to other phases
            if not self.contents["Phase"][i] == "stellar": continue

            # Search for the line
            if "photon packages per wavelength per process" in self.contents["Message"][i]:

                # Get the number of stellar photon packages per process
                packages = int(self.contents["Message"][i].split("(")[1].split(" photon")[0])
                return packages

        # Not found
        return None

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_packages(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries in reversed order
        for i in reversed(range(len(self.contents))):

            # Skip entries not corresponding to the dust emission phase
            if not self.contents["Phase"][i] == "dust": continue

            # Search for the line stating the number of photon packages
            if "photon packages for each of" in self.contents["Message"][i]:

                # Return the number of dust emission photon packages
                dust_packages = int(self.contents["Message"][i].split("(")[1].split(" photon")[0])
                return dust_packages

        # If the number of dust photon packages could not be determined, return None
        return None

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_packages_per_process(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries in reversed order
        for i in reversed(range(len(self.contents))):

            # Skip entries not corresponding to the dust emission phase
            if not self.contents["Phase"][i] == "dust": continue

            # Search for the line stating the number of photon packages per process
            if "photon packages per wavelength per process" in self.contents["Message"][i]:

                # Return the number of dust emission photon packages per process
                dust_packages = int(self.contents["Message"][i].split("(")[1].split(" photon")[0])
                return dust_packages

        # If the number of dust photon packages per process could not be determined, return None
        return None

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_cells(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Skip entries not corresponding to the setup
            #if self.contents["Phase"][i] != "setup": continue

            line = self.contents["Message"][i]

            if "Absorbed Stellar Luminosity Table" in line:

                # Example: "Absorbed Stellar Luminosity Table is not distributed. Size is (256000,160)"
                ncells = int(line.split(" (")[1].split(",")[0])
                return ncells

        # Not found, try tree dust grid
        return self.dust_cells_tree

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_cells_tree(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Search for the line stating the total number of leafs in the tree
            if "Total number of leaves" in self.contents["Message"][i]:

                # Return the number of leaves
                leaves = int(self.contents["Message"][i].split(": ")[1])
                return leaves

        # If the number of nodes could not be determined, return None
        return None

    # -----------------------------------------------------------------

    @property
    def ndust_cells(self):
        return self.dust_cells

    # -----------------------------------------------------------------

    @lazyproperty
    def uses_transient_heating(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Only look during the setup
            # Skip entries corresponding to other phases
            if self.contents["Phase"][i] != "setup": continue

            # Search for a particular line
            if "Precalculating cached values for transient dust emissivity computations" in self.contents["Message"][i]: return True

        # Not found
        return False

    # -----------------------------------------------------------------

    @lazyproperty
    def selfabsorption(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Search for particular line
            if "Starting the dust self-absorption phase" in self.contents["Message"][i]: return True

        # Not found
        return False

    # -----------------------------------------------------------------

    @lazyproperty
    def write_grid(self):

        """
        This function ...
        :return:
        """

        pass

        #"Will be outputting 3D grid data up to level" # Only for tree dust grids ...

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_file(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Only look during the setup
            # Skip entries corresponding to other phases
            if self.contents["Phase"][i] != "setup": continue

            line = self.contents["Message"][i]

            # Search for a particular line
            if "Reading wavelength grid data from file" in line:
                filepath = line.split("from file ")[1].split("...")[0]
                return filepath

        return None

    # -----------------------------------------------------------------

    @lazyproperty
    def uses_tree(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Only look during the setup
            # Skip entries corresponding to other phases
            if self.contents["Phase"][i] != "setup": continue

            # Search for a particular line
            if "Starting subdivision of level" in self.contents["Message"][i]: return True

        # If the line was not found
        return False

    # -----------------------------------------------------------------

    @lazyproperty
    def tree_nodes(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Only look during the setup
            # Skip entries corresponding to other phases
            if self.contents["Phase"][i] != "setup": continue

            # Search for the line stating the total number of nodes in the tree
            if "Total number of nodes" in self.contents["Message"][i]:

                # Return the number of nodes
                nodes = int(self.contents["Message"][i].split(": ")[1])
                return nodes

        # If the number of nodes could not be determined, return None
        return None

    # -----------------------------------------------------------------

    @property
    def tree_leaf_distribution(self):

        """
        This function ...
        :return:
        """

        from .tree import DustGridTreeDistribution

        # Keep track of the levels and corresponding number of cells
        levels = []
        counts = []

        # Current search level
        level = None

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Triggered
            if "Number of leaf cells of each level" in self.contents["Message"][i]: level = 0

            # Get the number of cells at this level
            elif level is not None:

                level_string = "Level " + str(level)

                # Check whether this level exists
                if level_string in self.contents["Message"][i]:

                    # Get the number of cells for this level
                    cells = int(self.contents["Message"][i].split(level_string + ": ")[1].split(" cells")[0])

                    # Add entries to the appropriate lists
                    levels.append(level)
                    counts.append(cells)

                    level += 1

                # This level does not exist anymore in the tree, return the result as a distribution
                else:
                    #return Distribution.from_counts("Number of cells", counts, levels)
                    return DustGridTreeDistribution.from_counts("Number of cells", counts, levels)

        # If the tree leaf distribution could not be determined, return None
        return None

    # -----------------------------------------------------------------

    @lazyproperty
    def tree_levels(self):

        """
        This function ...
        :return:
        """

        level = None

        # Loop over the log entries
        for i in range(len(self.contents)):

            if "Starting subdivision of level" in self.contents["Message"][i]:
                level = int(self.contents["Message"][i].split("of level ")[1].split("...")[0])

            elif "Construction of the tree finished" in self.contents["Message"][i]: return level

        # If the number of tree levels could not be determined, return None
        return None

    # -----------------------------------------------------------------

    @lazyproperty
    def npopulations(self):

        """
        This function ...
        :return:
        """

        npopulations = 0

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Only look during the setup
            # Skip entries corresponding to other phases
            if self.contents["Phase"][i] != "setup": continue

            # Search for specific messages
            if "Adding dust population" in self.contents["Message"][i]:
                #current_npopulations = int(self.contents["Message"][i].split("#")[1].split(" based on")[0])+1
                #npopulations = current_npopulations
                npopulations += 1

        # Return the number of populations
        return npopulations

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grain_types(self):

        """
        This function ...
        :return:
        """

        types = []

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Only look during the setup
            # Skip entries corresponding to other phases
            if self.contents["Phase"][i] != "setup": continue

            if "Grain composition grid" in self.contents["Message"][i]:
                types.append(self.contents["Message"][i].split("grid (")[1].split(")")[0])

        # Return the dust grain types
        return types

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths(self):

        """
        This function ...
        :return:
        """

        # Loop over the log entries
        for i in range(len(self.contents)):

            # Search for the line stating the number of wavelengths
            if "photon packages for each of" in self.contents["Message"][i]:

                # Return the number of wavelengths
                #wavelengths = int(self.contents["Message"][i].split("for each of ")[1].split(" wavelengths")[0])
                wavelengths = int(self.contents["Message"][i].split("for each of ")[1].split(" ")[0]) # There was a corrupted log line once, and this helped
                return wavelengths

        # If the number of wavelengths is not found, return None
        return None

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):
        return self.wavelengths

    # -----------------------------------------------------------------

    @lazyproperty
    def processes(self):

        """
        This function ...
        :return:
        """

        if self.starting_simulation_message is None: raise ValueError("There is something wrong with the log file")

        # Read number of processes from "Starting simulation ..." message
        if "with" in self.starting_simulation_message:
            processes = int(self.starting_simulation_message.split(' with ')[1].split()[0])
            return processes
        else: return 1

    # -----------------------------------------------------------------

    @property
    def nprocesses(self):
        return self.processes

    # -----------------------------------------------------------------

    @lazyproperty
    def threads(self):

        """
        This function ...
        :return:
        """

        # Indicate
        triggered = False
        max_thread_index = 0

        # Loop over all the log file entries
        for i in range(len(self.contents)):

            # Get the current log message
            message = self.contents["Message"][i]

            # Look for the message that indicates the start of the simulation
            if "Initializing random number generator" in message:

                triggered = True
                max_thread_index = int(message.split("thread number ")[1].split(" with seed")[0])

            # If the interesting log messages have been encountered, but the current log message doesn't state a thread
            # index, return the number of threads
            elif triggered:

                threads = max_thread_index + 1
                return threads

            # If none of the above, we are still in the part before the setup of the Random class

        # We should not get here
        raise ValueError("The number of threads could not be determined from the log file")

    # -----------------------------------------------------------------

    @property
    def nthreads(self):
        return self.threads

    # -----------------------------------------------------------------

    @lazyproperty
    def data_parallel(self):

        """
        This function ...
        :return:
        """

        # Loop over all log file entries
        for i in range(len(self.contents)):

            # Get the current log message
            message = self.contents["Message"][i]

            # Look for a message that indicates whether the Absorbed Stellar Luminosity Table is distributed or not
            if "Absorbed Stellar Luminosity Table" in message:

                if "is not distributed" in message: return False
                elif "is distributed" in message: return True
                else: raise ValueError("Log message truncated")

            # Alternatively, look at the Dust Emission Spectra Table

        # Return false if no messages regarding the distributed-ness of the tables was encountered
        return False

        # Simple implementation (only works if the log file is complete)
        #return "in data parallelization mode" in self.finished_simulation_message

        # Simple implementation (should be enough, but I was working with log files where this was not yet present)
        #return "in data parallelization mode" in self.starting_simulation_message

    # -----------------------------------------------------------------

    @lazyproperty
    def parallelization(self):

        """
        Thisf unction ...
        :return:
        """

        # WARNING: THE NUMBER OF THREADS PER CORE IS NOT KNOWN (1 IS ASSUMED!)
        from .parallelization import Parallelization
        return Parallelization.from_processes_and_threads(self.nprocesses, self.nthreads, data_parallel=self.data_parallel)

# -----------------------------------------------------------------

def parse_from_lines(lines, return_options=False):

    """
    This function ...
    :param lines: lines (can be generator or file handle)
    :param return_options:
    :return:
    """

    # Initialize lists for the columns
    times = []
    phases = []
    messages = []
    types = []
    memories = []

    verbose_logging = None
    memory_logging = None

    # PARSING

    # The current phase
    current_phase = None

    # The phase before that
    previous_phase = None

    # The phase even before that
    previousprevious_phase = None

    # Loop over all lines in the log file
    for line in lines:

        # If the line contains an error, skip it (e.g. when convergence has not been reached after a certain
        # number of dust-selfabsorption cycles)
        if "*** Error:" in line: continue

        # Remove the line ending (if there)
        line = line.rstrip("\n")

        #print(line)

        # Check whether the timestamp is valid
        # 17/09/2017 19:51:29.080
        #print(line)
        if time.has_valid_timestamp(line): t = time.parse_line(line) # Get the date and time information of the current line
        else:
            index = 1
            while True:
                if index > len(line) - 26:
                    index = None
                    break
                if time.has_valid_timestamp(line[index:]): break
                index += 1
            if index is None:
                warnings.warn("Not a valid line: '" + line + "': skipping ...")
                continue
            else:
                line = line[index:]
                t = time.parse_line(line)

        # Add the time
        times.append(t)

        # Check whether the log file was created in verbose logging mode
        if verbose_logging is None: verbose_logging = "[P" in line

        # Check whether the log file was created in memory logging mode
        if memory_logging is None: memory_logging = "GB)" in line

        # Get the memory usage at the current line, if memory logging was enabled for the simulation
        if memory_logging:

            try: memory = float(line.split(" (")[1].split(" GB)")[0])
            except ValueError:
                warnings.warn("Invalid line: '" + line + "': cannot interpret memory")
                memory = None
            memories.append(memory)

        if memory_logging: message = line.split("GB) ")[1]
        elif verbose_logging: message = line.split("] ")[1]
        else: message = line[26:]
        messages.append(message)

        typechar = line[24]
        if typechar == " ": types.append("info")
        elif typechar == "-": types.append("success")
        elif typechar == "!": types.append("warning")
        elif typechar == "*": types.append("error")
        else:
            #warnings.warn("Could not determine the type of log message from '" + line + "'")
            types.append("info")

        # Get the simulation phase
        current_phase, previous_phase, previousprevious_phase = get_phase(line, current_phase, previous_phase, previousprevious_phase)
        phases.append(current_phase)

    # Create the table data structures
    data = [times, phases, messages, types]
    names = ["Time", "Phase", "Message", "Type"]

    # If memory logging was enabled, add the 2 additional columns
    if memory_logging:
        data.append(memories)
        names.append("Memory")

    # Create the table and return it
    table = Table(data=data, names=names, meta={"name": "the contents of the simulation's log file"})
    if return_options: return table, verbose_logging, memory_logging
    else: return table

# -----------------------------------------------------------------

def parse(path, return_options=False):

    """
    This function ...
    :param path:
    :param return_options:
    :return:
    """

    # Open the log file
    with open(path, 'r') as fh: return parse_from_lines(fh, return_options=return_options)

# -----------------------------------------------------------------

def get_phase(line, current, previous, previousprevious):

    """
    This function ...
    :param line:
    :param current:
    :param previous:
    :param previousprevious:
    :return:
    """

    if current == "comm":
        if search_end(line):
            if previous == "wait" and (previousprevious == "stellar" or previousprevious == "write"): return previousprevious, current, None
            elif previous == "stellar" or previous == "write": return previous, current, previous

    # Search for the different simulation phases
    if search_start(line): return "start", current, previous
    elif search_end(line): return None, current, previous
    elif search_setup(line): return "setup", current, previous
    elif search_wait(line): return "wait", current, previous
    elif search_comm(line): return "comm", current, previous
    elif search_stellar(line): return "stellar", current, previous
    elif search_spectra(line): return "spectra", current, previous
    elif search_dust(line): return "dust", current, previous
    elif search_write(line): return "write", current, previous
    else: return current, previous, previousprevious # phase still going on

# -----------------------------------------------------------------

def get_last_phase(lines):

    """
    This function ...
    :return:
    """

    # The current phase
    current_phase = None

    # The phase before that
    previous_phase = None

    # The phase even before that
    previousprevious_phase = None

    # The start index of the current phase
    start_index = 0

    # Loop over the log lines
    for index, line in enumerate(lines):

        # Remember current phase before checking next line
        current_phase_before = current_phase

        # Get the simulation phase
        current_phase, previous_phase, previousprevious_phase = get_phase(line, current_phase, previous_phase, previousprevious_phase)

        # If new phase, set start index
        if current_phase != current_phase_before: start_index = index

    # Return the current phase
    return current_phase, start_index

# -----------------------------------------------------------------

def search_start(line):

    """
    This function ...
    :param line:
    :return:
    """

    return "Starting simulation" in line

# -----------------------------------------------------------------

def search_end(line):

    """
    This function ...
    :param line:
    :return:
    """

    if "Finished setup" in line: return True
    elif "Finished the stellar emission phase" in line: return True
    elif "Finished communication of the absorbed luminosities" in line: return True
    elif "Finished communication of the dust emission spectra" in line: return True
    elif "Finished the" in line and "-stage dust self-absorption cycle" in line: return True
    elif "Finished the dust emission phase" in line: return True
    elif "Finished writing results" in line: return True
    elif "Finished communication of" in line: return True # patch Dries, test why this is necessary
    else: return False

# -----------------------------------------------------------------

def search_setup(line):

    """
    This function ...
    :param line:
    :return:
    """

    if "Starting setup" in line: return True
    elif "Finished communication of the dust densities" in line: return True
    else: return False

# -----------------------------------------------------------------

def search_wait(line):

    """
    This function ...
    :param line:
    :return:
    """

    if "Waiting for other processes to finish the calculation of the dust cell densities" in line: return True
    elif "Waiting for other processes to finish the setup" in line: return True
    elif "Waiting for other processes to finish the stellar emission phase" in line: return True
    elif "Waiting for other processes to finish the emission spectra calculation" in line: return True
    elif "Waiting for other processes to finish this self-absorption cycle" in line: return True
    elif "Waiting for other processes to finish the dust emission phase" in line: return True
    else: return False

# -----------------------------------------------------------------

def search_comm(line):

    """
    This function ...
    :param line:
    :return:
    """

    if "Starting communication of the dust densities" in line: return True
    elif "Starting communication of the absorbed luminosities" in line: return True
    elif "Starting communication of the dust emission spectra" in line: return True
    elif "Starting communication of" in line: return True # patch Dries, test why this is necessary
    else: return False

# -----------------------------------------------------------------

def search_stellar(line):

    """
    This function ...
    :param line:
    :return:
    """

    return "Starting the stellar emission phase" in line

# -----------------------------------------------------------------

def search_spectra(line):

    """
    This function ...
    :param line:
    :return:
    """

    return "Library entries in use" in line

# -----------------------------------------------------------------

def search_dust(line):

    """
    This function ...
    :param line:
    :return:
    """

    return "Dust emission spectra calculated" in line

# -----------------------------------------------------------------

def search_write(line):

    """
    This function ...
    :param line:
    :return:
    """

    return "Starting writing results" in line

# -----------------------------------------------------------------

def get_start_line(lines):

    """
    This function ...
    :param lines:
    :return:
    """

    # Loop over all the lines
    for line in lines:

        # Look for the message that indicates the start of the simulation
        if "Starting simulation" in line: return line

    # Return None if the message could not be found
    return None

# -----------------------------------------------------------------

def get_nprocesses(lines):

    """
    This function ...
    :param lines:
    :return:
    """

    # Get line of the log file which states that the simulation is starting
    start = get_start_line(lines)
    if start is None: return None

    # Read number of processes from "Starting simulation ..." message
    if "with" in start:
        processes = int(start.split(' with ')[1].split()[0])
        return processes
    else: return 1

# -----------------------------------------------------------------

def get_simulation_phase(line, simulation_phase):

    """
    This function ...
    :param line:
    :param simulation_phase: current simulation phase
    :return:
    """

    # Check for marker for the setup phase
    if "Starting setup..." in line: simulation_phase = "SETUP"
    elif "Starting the stellar emission phase..." in line: simulation_phase = "STELLAR EMISSION"
    elif "Starting the dust self-absorption phase..." in line: simulation_phase = "DUST SELF-ABSORPTION"
    elif "Starting the dust emission phase..." in line: simulation_phase = "DUST EMISSION"
    elif "Starting writing results..." in line: simulation_phase = "WRITING"
    else: pass

    return simulation_phase

# -----------------------------------------------------------------

# def has_valid_timestamp(line):
#
#     """
#     This function ...
#     :param line:
#     :return:
#     """
#
#     # 17/09/2017 19:51:29.080
#
#     if line[2] != "/": return False
#     if line[5] != "/": return False
#     if line[10] != " ": return False
#     if line[13] != ":": return False
#     if line[16] != ":": return False
#     if line[19] != ".": return False
#
#     # All checks passed
#     return True

# -----------------------------------------------------------------
