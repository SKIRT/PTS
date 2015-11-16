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

class MemoryFile(object):

    """
    This class ...
    """

    def __init__(self, path):

        """
        The constructor ...
        :return:
        """

        # Set the memory file path
        self.path = path

        # Determine the name of the memory file
        name = os.path.basename(self.path)

        # Determine the simulation prefix
        self.prefix = name.split("_")[0]

        # Determine the process rank associated with this memory file
        try: self.process = int(name.split("_memoryP")[1].split(".txt")[0])
        except IndexError: self.process = 1

        # Parse the memory file
        self.contents = parse(path)

# *****************************************************************

def parse(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Initialize lists for the columns
    times = []
    deltas = []
    addresses = []

    multiprocessing = None

    # Open the memory file
    with open(path, 'r') as f:

        # Loop over all lines in the memory file
        for line in f:

            # Remove the line ending
            line = line[:-1]

            # Check whether the simulation was run in multiprocessing mode
            if multiprocessing is None: multiprocessing = "[P" in line

            # Get the memory increment or decrement
            index = 33 if multiprocessing else 26
            #print(line)
            if line[index] == "+": memory = float(line.split("+")[1].split(" GB")[0])
            elif line[index] == "-": memory = -float(line.split("-", 1)[1].split(" GB")[0])
            else: continue
            deltas.append(memory)

            # Get the date and time information of the current line
            t = time.parse(line)
            times.append(t)

            # Get the address of the associated Array
            address = None
            addresses.append(address)

    # Create and return a table with the contents of the memory file
    contents = Table([times, deltas, addresses], names=('Time', 'Delta', 'Address'), meta={"name": "the contents of the simulation's memory file"})
    contents["Delta"].unit = "GB"
    return contents

# *****************************************************************
