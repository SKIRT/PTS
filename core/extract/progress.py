#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.extract.progress Contains the ProgressTable class and the the ProgressExtractor class.
#  The latter class is used for extracting simulation progress from a simulation's log files into a ProgressTable object.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ..basics.log import log

# -----------------------------------------------------------------

class NoProgressData(Exception):

    """
    This class ...
    """

    def __init__(self, message, simulation_name=None):

        """
        Thisf unction ...
        :param message:
        :param simulation_name:
        """

        # Call the base class constructor with the parameters it needs
        super(NoProgressData, self).__init__(message)

        # The simulation name
        self.simulation_name = simulation_name

# -----------------------------------------------------------------

class ProgressTable(Table):

    """
    This function ...
    """

    @classmethod
    def from_columns(cls, process_list, phase_list, seconds_list, progress_list):

        """
        This function ...
        :param process_list:
        :param phase_list:
        :param seconds_list:
        :param progress_list:
        :return:
        """

        names = ["Process rank", "Phase", "Time", "Progress"]
        data = [process_list, phase_list, seconds_list, progress_list]

        # Call the constructor of the base class
        table = cls(data, names=names, masked=True)

        # Set the column units
        table["Time"].unit = "s"
        table["Progress"].unit = "%"

        table.path = None

        return table

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Open the table
        #table = cls.read(path, format="ascii.ecsv")
        table = super(ProgressTable, cls).read(path, format="ascii.ecsv")

        # Set the path
        table.path = path

        # Return the table
        return table

    # -----------------------------------------------------------------

    @classmethod
    def from_remote_file(cls, path, remote):

        """
        This function ...
        :param path:
        :param remote:
        :return:
        """

        # Open the contents
        contents = remote.get_text(path)

        # Open the table
        table = cls.read(contents, format="ascii.ecsv")

        # Return the table
        return table

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Save to the current path
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Write the table in ECSV format
        self.write(path, format="ascii.ecsv")

        # Set the path
        self.path = path

# -----------------------------------------------------------------

def extract_progress_cwd():

    """
    Thisf unction ...
    :return:
    """

    from pts.core.simulation.simulation import createsimulations

    # Create a SkirtSimulation object based on a log file present in the current working directory
    simulation = createsimulations(single=True)

    # Create a new ProgressExtractor instance
    extractor = ProgressExtractor()

    # Run the extractor and get the table
    extractor.run(simulation)
    table = extractor.table

    # Return the progress table
    return table

# -----------------------------------------------------------------

class ProgressExtractor(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # -- Attributes --

        self.log_files = None
        #self.staggered = None
        self.table = None

        # The output path
        self.output_path = None

    # -----------------------------------------------------------------

    def run(self, simulation, output_path=None):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation, output_path=output_path)

        # 2. Perform the extraction
        self.extract()

        # 3. Write the results
        if self.output_path is not None: self.write()

    # -----------------------------------------------------------------

    def setup(self, simulation, output_path=None):

        """
        This function ...
        :param simulation:
        :param output_path:
        :return:
        """

        # Obtain the log files created by the simulation
        self.log_files = simulation.logfiles()

        # Determine whether the emission spectra calculation was performed using a staggered assignment scheme
        # self.staggered = simulation.parameters().staggered()

        # Set the output path
        self.output_path = output_path

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting ...")

        number_of_processes = None

        # Initialize lists for the columns
        process_list = []
        phase_list = []
        seconds_list = []
        progress_list = []

        # Loop over the log files again and fill the column lists
        for log_file in self.log_files:

            # Get the total number of processes
            if number_of_processes is None: number_of_processes = log_file.processes
            else: assert number_of_processes == log_file.processes

            # Get the process rank associated with this log file
            process = log_file.process

            stellar_start = None
            spectra_start = None
            dust_start = None

            first_spectra_phase = True
            last_dust_phase = False

            total_entries = None
            entries_per_process = None

            # Loop over the entries in the log file
            for i in range(len(log_file.contents)):

                # Get the description of the current simulation phase
                phase = log_file.contents["Phase"][i]

                # The current log message
                message = log_file.contents["Message"][i]

                # The log file entries corresponding to the stellar emission phase
                if phase == "stellar":

                    # If this is the log message that marks the very start of the stellar emission phase, record the associated time
                    if "photon packages for" in message:

                        stellar_start = log_file.contents["Time"][i]

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append(phase)

                        # Add the seconds entry
                        seconds_list.append(0.0)

                        # Get the progress and add it to the list
                        progress_list.append(0.0)

                    # If this is one of the log messages that log stellar emission progress
                    elif "Launched stellar emission photon packages" in message:

                        # Add the seconds entry
                        seconds = (log_file.contents["Time"][i] - stellar_start).total_seconds()

                        # Get the progress and add it to the list
                        try: progress = float(message.split("packages: ")[1].split("%")[0])
                        except: continue  # INVALID LINE

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append(phase)

                        # Add the seconds and progress
                        seconds_list.append(seconds)
                        progress_list.append(progress)

                # The log file entries corresponding to the stellar emission phase
                elif phase == "spectra" and first_spectra_phase:

                    # If this is the log message that marks the very start of the spectra calculation, record the associated time
                    # If this log message states the total number of library entries that are used, record this number
                    if "Library entries in use" in message:

                        spectra_start = log_file.contents["Time"][i]

                        # Get the total number of library entries in use and the number of entries per process
                        total_entries = int(message.split("use: ")[1].split(" out of")[0])
                        entries_per_process = total_entries / number_of_processes

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append(phase)

                        # Add the seconds entry
                        seconds_list.append(0.0)

                        # Get the progress and add it to the list
                        progress_list.append(0.0)

                    elif "Calculating emission for" in message:

                        entry = float(message.split()[-1][:-3])

                        # Determine the progress
                        #if self.staggered: fraction = entry / total_entries
                        #else: fraction = (entry - process * entries_per_process) / entries_per_process

                        fraction = entry / total_entries

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append(phase)

                        # Add the seconds entry
                        seconds = (log_file.contents["Time"][i] - spectra_start).total_seconds()
                        seconds_list.append(seconds)

                        # Get the progress and add it to the list
                        progress = float(fraction*100.0)
                        progress_list.append(progress)

                # The log file entries corresponding to the dust emission phase
                # We only want to record the progress of the 'last' dust emission phase
                elif phase == "dust" and last_dust_phase:

                    # If this is the log message that marks the very start of the dust emission phase, record the associated time
                    if "photon packages for" in message:

                        dust_start = log_file.contents["Time"][i]

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append(phase)

                        # Add the seconds entry
                        seconds_list.append(0.0)

                        # Get the progress and add it to the list
                        progress_list.append(0.0)

                    # If this is one of the log messages that log dust emission progress
                    elif "Launched dust emission photon packages" in message:

                        # Add the seconds entry
                        seconds = (log_file.contents["Time"][i] - dust_start).total_seconds()

                        # Get the progress and add it to the list
                        try: progress = float(message.split("packages: ")[1].split("%")[0])
                        except: continue # INVALID LINE

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append(phase)

                        # Add the seconds and progress
                        seconds_list.append(seconds)
                        progress_list.append(progress)

                # Record the end of the spectra calculation (the first log message of the emission phase of the self-absorption cycle)
                elif phase == "dust" and first_spectra_phase:

                    # If this line indicates the end of the dust emission spectra calculation
                    if "Dust emission spectra calculated" in message:

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append("spectra")

                        # Add the seconds entry
                        seconds = (log_file.contents["Time"][i] - spectra_start).total_seconds()
                        seconds_list.append(seconds)

                        # Add 100% progress to the list
                        progress_list.append(100.0)

                        # Indicate that the first spectra phase has already been processed (subsequent spectra phases can be ignored)
                        first_spectra_phase = False

                # Log messages that fall in between phases
                elif phase is None:

                    # The current log message
                    message = log_file.contents["Message"][i]

                    # Look for messages indicating whether this dust photon shooting phase corresponds to
                    # one of the dust self-absorption cycles or the actual dust emission phase
                    if "dust self-absorption cycle" in message: last_dust_phase = False
                    elif "Starting the dust emission phase" in message: last_dust_phase = True

                    elif "Finished the stellar emission phase" in message:

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append("stellar")

                        # Add the seconds entry
                        seconds = (log_file.contents["Time"][i] - stellar_start).total_seconds()
                        seconds_list.append(seconds)

                        # Add 100% progress to the list
                        progress_list.append(100.0)

                    elif "Finished the dust emission phase" in message:

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append("dust")

                        # Add the seconds entry
                        seconds = (log_file.contents["Time"][i] - dust_start).total_seconds()
                        seconds_list.append(seconds)

                        # Add 100% progress to the list
                        progress_list.append(100.0)

        # Create the progress table
        self.table = ProgressTable.from_columns(process_list, phase_list, seconds_list, progress_list)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the table to file
        self.table.saveto(self.output_path)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the table to None
        self.table = None

# -----------------------------------------------------------------
