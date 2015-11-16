#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package extract.memory Extract memory information from SKIRT or SkirtMemory output

# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from datetime import datetime

# Import astronomical modules
from astropy.table import Table

# *****************************************************************

class MemoryExtractor(object):

    """
    This function ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        pass

    # *****************************************************************

    def run(self, simulation, output_path):

        """
        This function ...
        :return:
        """

        # Obtain the log files and memory files created by the simulation
        log_files = simulation.logfiles()
        memory_files = simulation.memoryfiles()

        # Check whether memory mode has been used for the simulation
        memory_mode = len(memory_files) > 0

        # Initialize lists for the columns
        process_list = []
        phase_list = []
        seconds_list = []
        memory_list = []
        delta_list = []

        # Loop over all log files to determine the earliest recorded time
        t0 = datetime.now()
        for log_file in log_files:
            if log_file.t0 < t0: t0 = log_file.t0

        # Loop over the log files again and fill the column lists
        for log_file in log_files:

            # Get the process rank associated with this log file
            process = log_file.process

            # If memory mode was used
            if memory_mode:

                # Get the associated memory file
                memory_file = memory_files[process]
                assert process == memory_file.process

            # Create two counter variables
            j = k = 0

            # Loop over both the entries of the log file and the memory file
            end_log = False
            end_memory = False if memory_mode else True
            while True:

                # Check whether to end looping over either the log or memory file's contents
                if j == len(log_file.contents): end_log = True
                if k == len(memory_file.contents): end_memory = True

                # Calculate the number of seconds that have passed since the earliest recorded log time
                seconds_log = (log_file.contents["Time"][j] - t0).total_seconds()
                seconds_memory = (memory_files.contents["Time"][k] - t0).total_seconds()

                if (not end_log) and seconds_log < seconds_memory:

                    # Fill in the column lists
                    process_list.append(process)
                    phase_list.append(log_file.contents["Phase"][j])
                    seconds_list.append(seconds_log)
                    memory_list.append(log_file.contents["Memory"][j])
                    delta_list.append(None)

                    j += 1

                elif (not end_memory) and seconds_memory < seconds_log:

                    # Fill in the column lists
                    process_list.append(process)
                    phase_list.append(phase_list[len(phase_list)-1])
                    seconds_list.append(seconds_memory)
                    memory_list.append(None)
                    delta_list.append(memory_file.contents["Delta"][k])

                    k += 1

                # Exit the double loop
                else: break

        # Create a table
        table = Table([process_list, phase_list, seconds_list, memory_list, delta_list], names=('Process rank', 'Simulation phase', 'Simulation time', 'Memory usage', 'Memory (de)allocation'))
        table["Simulation time"].unit = "s"
        table["Memory usage"].unit = "GB"
        table["Memory (de)allocation"].unit = "GB"

        # Write the table to file
        table.write(output_path, format="ascii")

# *****************************************************************

# Calculate the number of wavelengths and dust cells
#Nlambda = skifile.nwavelengths()
#Ncells = skifile.ncells()
#Ndoubles = Nlambda * Ncells

#plt.plot(seconds_list, memory_list)

#totals = np.cumsum(deltas)

#for i in range(len(times)):
#    print times[i], totals[i]

#plt.step(times, totals)
#plt.fill_between(times, totals, color='green')
#plt.bar(times, totals, color='r')
#plt.show()
