#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package core.logfile

# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS modules
from ..tools import time

# *****************************************************************

class LogFile(object):

    """
    This class ...
    """

    def __init__(self, path):

        """
        The constructor ...
        :return:
        """

        # Set the log file path
        self.path = path

        # Determine the name of the log file
        name = os.path.basename(self.path)

        # Determine the simulation prefix
        self.prefix = name.split("_")[0]

        # Determine the process rank associated with this log file
        try: self.process = int(name.split("_logP")[1].split(".txt")[0])
        except IndexError: self.process = 1

        # Parse the log file
        self.contents = parse(path)

        # Determine the number of processes and threads
        self.processes = get_processes(self.contents)
        self.threads = get_threads(self.contents)

# *****************************************************************

def get_processes(table):

    """
    This function ...
    :param table:
    :return:
    """

    for message in table["Message"]:
        if "Starting simulation" in message:
            if "with" in message: return int(message.split(' with ')[1].split()[0])
            else: return 1
    raise ValueError("Cannot determine the number of processes from the log file")

# *****************************************************************

def get_threads(table):

    """
    This function ...
    :param table:
    :return:
    """

    triggered = False
    max_thread_number = 0
    for message in table["Message"]:
        if "Initializing random number generator" in message:
            triggered = True
            max_thread_number = int(message.split("thread number ")[1].split(" with seed")[0])
        elif triggered: return max_thread_number+1
    raise ValueError("Cannot determine the number of threads from the log file")

# *****************************************************************

def parse(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Initialize lists for the columns
    times = []
    phases = []
    messages = []
    types = []
    memories = []

    verbose_logging = None
    memory_logging = True

    # Open the log file
    with open(path, 'r') as f:

        # The current phase
        current_phase = None

        # Loop over all lines in the log file
        for line in f:

            # Remove the line ending
            line = line[:-1]

            # Get the date and time information of the current line
            t = time.parse(line)
            times.append(t)

            # Check whether the log file was created in verbose logging mode
            if verbose_logging is None: verbose_logging = "[P" in line

            # Get the memory usage at the current line, if memory logging was enabled for the simulation
            if memory_logging:
                try:
                    memory = float(line.split(" (")[1].split(" GB)")[0])
                    memories.append(memory)
                except IndexError:
                    memory_logging = False

            if memory_logging:
                message = line.split("GB) ")[1]
            elif verbose_logging:
                message = line.split("] ")[1]
            else:
                message = line[26:]
            messages.append(message)

            typechar = line[24]
            if typechar == " ": types.append("info")
            elif typechar == "-": types.append("success")
            elif typechar == "!": types.append("warning")
            elif typechar == "*": types.append("error")
            else: raise ValueError("Could not determine the type of log message")

            # Get the simulation phase
            current_phase = get_phase(line, current_phase)
            phases.append(current_phase)

    # Set the contents variable
    if memory_logging:
        return Table([times, phases, messages, types, memories], names=('Time', 'Phase', 'Message', 'Type', 'Memory'), meta={"name": "the contents of the simulation's log file"})
    else:
        return Table([times, phases, messages, types], names=('Time', 'Phase', 'Message', 'Type'), meta={"name": "the contents of the simulation's log file"})

# *****************************************************************

def get_phase(line, current):

    """
    This function ...
    :param line:
    :return:
    """

    # Search for the different simulation phases
    if search_start(line): return "start"
    elif search_end(line): return None
    elif search_setup(line): return "setup"
    elif search_wait(line): return "wait"
    elif search_comm(line): return "comm"
    elif search_stellar(line): return "stellar"
    elif search_spectra(line): return "spectra"
    elif search_dust(line): return "dust"
    elif search_write(line): return "write"
    else: return current

# *****************************************************************

def search_start(line):

    """
    This function ...
    :param line:
    :return:
    """

    if "Starting simulation" in line: return True
    else: return False

# *****************************************************************

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
    else: return False

# *****************************************************************

def search_setup(line):

    """
    This function ...
    :param line:
    :return:
    """

    if "Starting setup" in line: return True
    elif "Finished communication of the dust densities" in line: return True
    else: return False

# *****************************************************************

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

# *****************************************************************

def search_comm(line):

    """
    This function ...
    :param line:
    :return:
    """

    if "Starting communication of the dust densities" in line: return True
    elif "Starting communication of the absorbed luminosities" in line: return True
    elif "Starting communication of the dust emission spectra" in line: return True
    else: return False

# *****************************************************************

def search_stellar(line):

    """
    This function ...
    :param line:
    :return:
    """

    if "Starting the stellar emission phase" in line: return True
    else: return False

# *****************************************************************

def search_spectra(line):

    """
    This function ...
    :param line:
    :return:
    """

    if "Library entries in use" in line: return True
    else: return False

# *****************************************************************

def search_dust(line):

    """
    This function ...
    :param line:
    :return:
    """

    if "Dust emission spectra calculated" in line: return True
    else: return False

# *****************************************************************

def search_write(line):

    """
    This function ...
    :param line:
    :return:
    """

    if "Starting writing results" in line: return True
    else: return False

# *****************************************************************
