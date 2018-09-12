#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.extract.timeline Contains the TimeLineTable and TimeLineExtractor classes. The latter class is used
#  for extracting timeline information from a simulation's log files to a TimeLineTable object.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from datetime import datetime

# Import astronomical modules
from ..tools import tables

# Import the relevant PTS classes and modules
from pts.core.tools.utils import lazyproperty
from ..basics.log import log
from ..basics.table import SmartTable

# -----------------------------------------------------------------

def extract_timeline(simulation=None, log_files=None):

    """
    This function ...
    :param simulation:
    :param log_files:
    :return:
    """

    # Create a TimeLineExtractor instance
    extractor = TimeLineExtractor()

    # Run the timeline extractor
    extractor.run(simulation=simulation, log_files=log_files)

    # Return the timeline
    return extractor.table

# -----------------------------------------------------------------

class NoTimingData(Exception):

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
        super(NoTimingData, self).__init__(message)

        # The simulation name
        self.simulation_name = simulation_name

# -----------------------------------------------------------------

class TimeLineTable(SmartTable):

    """
    This function ...
    """

    @classmethod
    def from_columns(cls, process_list, phase_list, start_list, end_list, simulation_phase_list, annotation_list):

        """
        This function ...
        :param process_list:
        :param phase_list:
        :param start_list:
        :param end_list:
        :param simulation_phase_list:
        :param annotation_list:
        :return:
        """

        # Define column names and units
        names = ["Process rank", "Phase", "Start time", "End time", "Simulation phase", "Annotation"]
        units = [None, None, "s", "s", None, None]

        # NEW
        return super(TimeLineTable, cls).from_columns(process_list, phase_list, start_list, end_list, simulation_phase_list, annotation_list, names=names, units=units)

    # -----------------------------------------------------------------

    @lazyproperty
    def nprocesses(self):

        """
        This function ...
        :return:
        """

        return max(self["Process rank"]) + 1

    # -----------------------------------------------------------------

    def duration(self, phase, single=False, only_root=False, simulation_phase=None, annotation_contains=None, annotation=None):

        """
        This function ...
        :param phase:
        :param single:
        :param only_root:
        :param simulation_phase:
        :param annotation_contains:
        :param annotation:
        :return:
        """

        # Keep track of the total amount of time spent in the specified phase
        total = 0.0

        assert self["Process rank"][0] == 0

        processes = []
        skip_to_next_process = False

        # Loop over the table rows
        for i in range(len(self)):

            # Get tbe process rank
            rank = self["Process rank"][i]

            # Only add the contributions from the root process
            if only_root and rank > 0: break

            # Skip to next process
            if skip_to_next_process:
                if rank == processes[-1]: continue
                else: skip_to_next_process = False

            # Add the process rank so that we can check whether we have encounterd each of them at the end
            if rank not in processes: processes.append(rank)

            # Check simulation phase
            if simulation_phase is not None and self["Simulation phase"][i] != simulation_phase: continue

            # Check annotation
            if annotation_contains is not None:
                if self["Annotation"][i] is None: continue
                elif annotation_contains not in self["Annotation"][i]: continue
            if annotation is not None and annotation != self["Annotation"][i]: continue

            # Check whether the current entry corresponds to the desired phase
            if self["Phase"][i] == phase:

                # Get the start and end time for the phase
                start = self["Start time"][i]
                end = self["End time"][i]

                # Calculate the time duration for this phase, returning it if single=True, otherwise add it to the total
                total += end - start
                if single: skip_to_next_process = True

        # Assert that we had every process
        if not only_root: assert sorted(list(processes)) == range(self.nprocesses)

        # Return the (average) total amount of time spent in the specified phase
        if only_root: return total
        else: return total / self.nprocesses

    # -----------------------------------------------------------------

    def duration_without(self, phases, only_root=False):

        """
        This function ...
        :param phases:
        :param only_root:
        :return:
        """

        # If no phases are given, set an empty list
        if phases is None: phases = []

        # Create a list of phases is only one is given
        if isinstance(phases, basestring): phases = [phases]

        # Keep track of the total amount of time spent in phases other than the specified phase
        total = 0.0

        assert self["Process rank"][0] == 0

        processes = []
        skip_to_next_process = False

        # Loop over the table rows
        for i in range(len(self)):

            # Get the process rank
            rank = self["Process rank"][i]

            # Only add the contributions from the root process
            if only_root and rank > 0: break

            # Skip to next process
            if skip_to_next_process:
                if rank == processes[-1]: continue
                else: skip_to_next_process = False

            # Add the process rank so that we can check whether we have encounterd each of them at the end
            if rank not in processes: processes.append(rank)

            # Check whether the current entry corresponds to a phase different from the specified phase
            if self["Phase"][i] not in phases:

                # Get the start and end time for this phase
                start = self["Start time"][i]
                end = self["End time"][i]

                # Add the duration to the total
                total += end - start

        # Assert that we had every process
        if not only_root: assert sorted(processes) == range(self.nprocesses)

        # Return the total amount of time spent in phases other than the specified phase
        if only_root: return total
        else: return total / self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def total(self):

        """
        This function ...
        :return:
        """

        return self.duration_without(None)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_root(self):

        """
        This function ...
        :return:
        """

        return self.duration_without(None, only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration_without(None) * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def setup(self):

        """
        This function ...
        :return:
        """

        return self.duration("setup")

    # -----------------------------------------------------------------

    @lazyproperty
    def setup_root(self):

        """
        This function ...
        :return:
        """

        return self.duration("setup", only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def setup_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration("setup") * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar(self):

        """
        This function ...
        :return:
        """

        return self.duration("stellar", single=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_root(self):

        """
        This function ...
        :return:
        """

        return self.duration("stellar", single=True, only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration("stellar", single=True) * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def spectra(self):

        """
        This function ...
        :return:
        """

        return self.duration("spectra")

    # -----------------------------------------------------------------

    @lazyproperty
    def spectra_root(self):

        """
        This function ...
        :return:
        """

        return self.duration("spectra", only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def spectra_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration("spectra") * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def dust(self):

        """
        This function ...
        :return:
        """

        return self.duration("dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_root(self):

        """
        This function ...
        :return:
        """

        return self.duration("dust", only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration("dust") * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def dustemission_spectra(self):

        """
        This function ...
        :return:
        """

        # Return the duration of the dust spectra calculation phase belonging to the (final) dust emission phase
        return self.duration("spectra", single=True, simulation_phase="DUST EMISSION")

    # -----------------------------------------------------------------

    @lazyproperty
    def dustemission_spectra_root(self):

        """
        This function ...
        :return:
        """

        return self.duration("spectra", single=True, simulation_phase="DUST EMISSION", only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def dustemission_spectra_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration("spectra", single=True, simulation_phase="DUST EMISSION") * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def dustemission_photons(self):

        """
        This function ...
        :return:
        """

        # Return the duration of the dust photon shooting phase belonging to the (final) dust emission phase
        return self.duration("dust", single=True, simulation_phase="DUST EMISSION")

    # -----------------------------------------------------------------

    @lazyproperty
    def dustemission_photons_root(self):

        """
        This function ...
        :return:
        """

        return self.duration("spectra", single=True, simulation_phase="DUST EMISSION", only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def dustemission_photons_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration("spectra", single=True, simulation_phase="DUST EMISSION") * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def writing(self):

        """
        This function ....
        :return:
        """

        return self.duration("write")

    # -----------------------------------------------------------------

    @lazyproperty
    def writing_root(self):

        """
        This function ...
        :return:
        """

        return self.duration("write", only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def writing_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration("write") * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def communication(self):

        """
        This function ...
        :return:
        """

        return self.duration("comm")

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_root(self):

        """
        This function ...
        :return:
        """

        return self.duration("comm", only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration("comm") * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_densities(self):

        """
        This function ...
        :return:
        """

        # Get the communication phase that is the communication of the dust densities during the setup
        return self.duration("comm", single=True, annotation_contains="dust densities")

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_densities_root(self):

        """
        This function ...
        :return:
        """

        return self.duration("comm", single=True, annotation_contains="dust densities", only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_densities_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration("comm", single=True, annotation_contains="dust densities") * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_stellar_absorption(self):

        """
        This function ...
        :return:
        """

        # Get the communication phase that is the (first) communication of the absorbed stellar luminosities
        return self.duration("comm", single=True, annotation_contains="absorbed stellar luminosity table")

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_dust_absorption(self):

        """
        This function ...
        :return:
        """

        # Get the communication phase that is the (first) communication of the absorbed stellar luminosities
        return self.duration("comm", single=True, annotation_contains="absorbed dust luminosity table")

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_emission(self):

        """
        This function ...
        :return:
        """

        # Get the communication phase ...
        return self.duration("comm", single=True, annotation_contains="dust emission spectra table")

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_instruments(self):

        """
        This function ...
        :return:
        """

        # Get the communication phase ...
        return self.duration("comm", single=True, annotation_contains="observed fluxes")

    # -----------------------------------------------------------------

    @lazyproperty
    def waiting(self):

        """
        This function ...
        :return:
        """

        return self.duration("wait")

    # -----------------------------------------------------------------

    @lazyproperty
    def waiting_root(self):

        """
        This function ...
        :return:
        """

        return self.duration("wait", only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def waiting_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration("wait") * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def other(self):

        """
        This function ...
        :return:
        """

        return self.duration(None)

    # -----------------------------------------------------------------

    @lazyproperty
    def other_root(self):

        """
        This function ...
        :return:
        """

        return self.duration(None, only_root=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def other_all_processes(self):

        """
        This function ...
        :return:
        """

        return self.duration(None) * self.nprocesses

    # -----------------------------------------------------------------

    @lazyproperty
    def serial(self):

        """
        This function ...
        :return:
        """

        return self.setup + self.writing + self.other

    # -----------------------------------------------------------------

    @lazyproperty
    def serial_root(self):

        """
        This function ...
        :return:
        """

        return self.setup_root + self.writing_root + self.other_root

    # -----------------------------------------------------------------

    @lazyproperty
    def parallel(self):

        """
        This function ...
        :return:
        """

        return self.stellar + self.spectra + self.dust

    # -----------------------------------------------------------------

    @lazyproperty
    def parallel_root(self):

        """
        This function ...
        :return:
        """

        return self.stellar_root + self.spectra_root + self.dust_root

    # -----------------------------------------------------------------

    @lazyproperty
    def overhead(self):

        """
        This function ...
        :return:
        """

        return self.communication + self.waiting

    # -----------------------------------------------------------------

    @lazyproperty
    def overhead_root(self):

        """
        This function ...
        :return:
        """

        return self.communication_root + self.waiting_root

# -----------------------------------------------------------------

class TimeLineExtractor(object):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        :return:
        """

        # -- Attributes --

        # The list of log files created by the simulation
        self.log_files = None

        # The table containing the timeline information
        self.table = None

        # The output path
        self.output_path = None

    # -----------------------------------------------------------------

    def run(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(*args, **kwargs)

        # 2. Perform the extraction
        self.extract()

        # 3. Write the results
        if self.output_path is not None: self.write()

    # -----------------------------------------------------------------

    def setup(self, *args, **kwargs):

        """
        This funtion ...
        :param args:
        :param kwargs:
        :return:
        """

        # Simulation is passed
        if kwargs.get("simulation", None) is not None:

            # Obtain the log files created by the simulation
            simulation = kwargs.pop("simulation")
            log_files = simulation.logfiles()

        # Log files are passed
        elif kwargs.get("log_files", None) is not None: log_files = kwargs.pop("log_files")

        # One log file is passed
        elif kwargs.get("log_file", None) is not None: log_files = [kwargs.pop("log_file")]

        # Not enough input
        else: raise ValueError("Too little input")

        # Set the log files
        self.log_files = log_files

        # Check whether there are log files
        if not self.has_logfiles: raise ValueError("No log files found in simulation output")

        # Set the output path
        #if len(args) > 0: output_path = args[1]
        #else: output_path = kwargs.pop("output_path", None)
        self.output_path = kwargs.pop("output_path", None)

    # -----------------------------------------------------------------

    @property
    def nlogfiles(self):

        """
        This function ...
        :return:
        """

        return len(self.log_files)

    # -----------------------------------------------------------------

    @property
    def has_logfiles(self):

        """
        This function ...
        :return:
        """

        return self.nlogfiles > 0

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting ...")

        # Initialize lists for the columns
        process_list = []
        phase_list = []
        start_list = []
        end_list = []
        simulation_phase_list = []
        annotation_list = []

        # Loop over all log files to determine the earliest recorded time
        t_0 = datetime.now()
        for log_file in self.log_files:
            if log_file.t_0 < t_0: t_0 = log_file.t_0

        # Keep track of the actual simulation phase (setup, stellar, self-absorption, dust emission, writing)
        simulation_phase = None

        # Keep track of the dust self-absorption stage and cycle
        dust_selfabsorption_stage = None
        dust_selfabsorption_cycle = None

        # Loop over the log files again and fill the column lists
        unique_processes = []
        for log_file in self.log_files:

            # Get the process rank associated with this log file
            process = log_file.process
            unique_processes.append(process)

            # Get the first message
            message = log_file.contents["Message"][0]

            # Keep track of the current phase while looping over the log file entries
            current_phase = log_file.contents["Phase"][0]
            process_list.append(process)
            phase_list.append(current_phase)
            start_list.append((log_file.t_0 - t_0).total_seconds())
            simulation_phase_list.append(get_simulation_phase(message, None))
            annotation_list.append(None)

            # Loop over all log file entries (lines)
            for j in range(len(log_file.contents)):

                # Get the description of the current simulation phase
                phase = log_file.contents["Phase"][j]

                # Get the current message
                message = log_file.contents["Message"][j]

                # Get the simulation phase
                simulation_phase = get_simulation_phase(message, simulation_phase)

                # Check for markers that indicate a dust self-absorption cycle
                if "Starting the" in message and "dust self-absorption cycle" in message:
                    stage_description = message.split("Starting the ")[1].split(" dust self-absorption")[0]
                    if stage_description == "first-stage": dust_selfabsorption_stage = 0
                    elif stage_description == "second-stage": dust_selfabsorption_stage = 1
                    elif stage_description == "last-stage": dust_selfabsorption_stage = 2
                    else: raise ValueError("Invalid dust self-absorption stage: " + stage_description)
                    dust_selfabsorption_cycle = int(message.split("dust self-absorption cycle ")[1].split("...")[0])

                if simulation_phase == "DUST EMISSION":
                    dust_selfabsorption_stage = None
                    dust_selfabsorption_cycle = None

                # If a new phase is not entered
                if phase == current_phase: continue

                # Determine the current time
                seconds = (log_file.contents["Time"][j] - t_0).total_seconds()

                # Check whether the current entry corresponds to the desired phase
                annotation = None
                if simulation_phase == "DUST SELF-ABSORPTION":
                    if phase == "dust": annotation = str(dust_selfabsorption_stage) + "," + str(dust_selfabsorption_cycle)
                    elif phase == "spectra": annotation = str(dust_selfabsorption_stage) + "," + str(dust_selfabsorption_cycle)
                if phase == "comm":
                    communication_of = log_file.contents["Message"][j].split("communication of")[1].split("...")[0].lower()
                    if communication_of.startswith("the "): communication_of = communication_of.split("the ")[1]
                    annotation = communication_of

                # Mark the end of the previous phase
                end_list.append(seconds)

                # Mark the start of the current phase
                process_list.append(process)
                phase_list.append(phase)
                start_list.append(seconds)

                # Set the simulation phase
                simulation_phase_list.append(simulation_phase)

                # Set the annotation
                annotation_list.append(annotation)

                # Update the current phase
                current_phase = phase

            # Add the last recorded time in the log file as the end of the last simulation phase
            end_list.append((log_file.t_last - t_0).total_seconds())

        # Fix for when the number of phases does not correspond between the different processes
        if len(unique_processes) > 1: verify_phases(process_list, phase_list, start_list, end_list, simulation_phase_list, annotation_list)

        # Create the table
        self.table = TimeLineTable.from_columns(process_list, phase_list, start_list, end_list, simulation_phase_list, annotation_list)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the table to file
        tables.write(self.table, self.output_path)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the table to None
        self.table = None

# -----------------------------------------------------------------

def verify_phases(process_list, phase_list, start_list, end_list, simulation_phase_list, annotation_list):

    """
    This function ...
    :param process_list:
    :param phase_list:
    :param start_list:
    :param end_list:
    :param simulation_phase_list:
    :param annotation_list:
    :return:
    """

    # TODO: right now, it does not fix everything ! only the case where the root process misses a phase at the end compared to the other processes

    # Get the number of processes
    nprocs = max(process_list) + 1

    # For each process rank, get the index in the columns where the entries of that process start
    process_indices = [None] * nprocs
    previous_process = -1
    for i in range(len(process_list)):
        if process_list[i] == previous_process: continue
        else:
            process_indices[process_list[i]] = i
            previous_process = process_list[i]
            if previous_process == nprocs - 1: break

    max_number_of_entries = max([process_indices[i+1]-process_indices[i] for i in range(len(process_indices)-1)])

    for entry_index in range(max_number_of_entries):

        # Check whether the phases correspond and if the process ranks are actually the ranks we loop over below
        phase_process_0 = phase_list[entry_index]
        if process_list[entry_index] != 0: # but already an entry of process rank 1
            if process_list[entry_index] != 1: raise RuntimeError("I don't know how to solve this inconsistent set of columns")

            missing_phase = phase_list[process_indices[1]+entry_index]
            missing_end = end_list[process_indices[1]+entry_index]

            missing_simulation_phase = simulation_phase_list[process_indices[1]+entry_index]
            missing_annotation = simulation_phase_list[process_indices[1]+entry_index]

            process_list.insert(entry_index, 0)
            phase_list.insert(entry_index, missing_phase)
            start_list.insert(entry_index, end_list[entry_index-1])
            end_list.insert(entry_index, missing_end)
            simulation_phase_list.insert(entry_index, missing_simulation_phase)
            annotation_list.insert(entry_index, missing_annotation)

        #for rank in range(1,nprocs):

            #print(phase_process_0, phase_list[process_indices[rank]+entry_index])
            #print(rank, process_list[process_indices[rank]+entry_index])

# -----------------------------------------------------------------

def get_simulation_phase(message, simulation_phase):

    """
    This function ...
    :param message:
    :param simulation_phase:
    :return:
    """

    # Check for marker for the setup phase
    if message == "Starting setup...": simulation_phase = "SETUP"
    elif message == "Starting the stellar emission phase...": simulation_phase = "STELLAR EMISSION"
    elif message == "Starting the dust self-absorption phase...": simulation_phase = "DUST SELF-ABSORPTION"
    elif message == "Starting the dust emission phase...": simulation_phase = "DUST EMISSION"
    elif message == "Starting writing results...": simulation_phase = "WRITING"
    else: pass

    return simulation_phase

# -----------------------------------------------------------------
