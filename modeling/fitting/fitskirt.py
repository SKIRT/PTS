#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.fitskirt Contains the FitSKIRTLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import math
import copy
import subprocess
from lxml import etree
from datetime import datetime

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import archive as arch
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.tools import introspection, monitoring
from ...core.tools import parallelization as par
from ...core.simulation.parallelization import Parallelization
from ...core.tools import xml
from ...core.tools import parsing
from ...core.prep.smile import skirt_quantities_to_pts_quantities
from ...core.basics.range import RealRange
from ...core.tools import stringify
from ...core.basics.map import Map
from ...core.tools import time

# -----------------------------------------------------------------

class FitSKIRTDefinition(object):

    """
    This function ...
    """

    def __init__(self, ski_path, fski_path, output_path, input_path=None, name=None):

        """
        The constructor ...
        :param ski_path:
        :param fski_path:
        :param output_path:
        :param input_path:
        :param name:
        """

        # Options for the ski file pattern
        self.ski_path = ski_path
        self.fski_path = fski_path

        # The input and output paths
        self.input_path = input_path
        self.output_path = output_path

        # A name for this simulation
        self.name = name

# -----------------------------------------------------------------

class FitSKIRTArguments(object):

    """
    This class ...
    """

    def __init__(self, definition, parallelization=None):

        """
        The constructor ...
        :param definition:
        :param parallelization:
        """

        print(definition.ski_path)
        print(definition.input_path)

        # Check whether the ski file and the other input are in the same directory
        if fs.directory_of(definition.ski_path) != definition.input_path:

            # Create new temporary input directory
            temp_path = fs.create_directory_in(introspection.pts_temp_dir, time.unique_name("FitSKIRT"))

            # Copy the input and the ski file
            fs.copy_files(fs.files_in_path(definition.input_path), temp_path)
            fs.copy_file(definition.ski_path, temp_path)

            # Set the

        # Set the name of the ski file in the fski file
        fski = FskiFile(definition.fski_path)
        fski.set_ski_name(fs.name(definition.ski_path))
        fski.save()

        # Options for the ski file pattern
        self.fski_path = definition.fski_path
        self.relative = None

        # The input and output paths
        self.input_path = definition.input_path
        self.output_path = definition.output_path

        # Options for parallelization
        self.parallel = Map()
        self.parallel.simulations = None
        self.parallel.threads = parallelization.threads if parallelization is not None else None
        self.parallel.processes = parallelization.processes if parallelization is not None else None

    # -----------------------------------------------------------------

    @property
    def prefix(self):

        """
        This function ...
        :return:
        """

        if not fs.is_file(self.fski_path): raise RuntimeError("Cannot determine the prefix for the fski pattern '" + self.fski_path + "'. Does it define multiple files?")
        return fs.strip_extension(fs.name(self.fski_path))

    # -----------------------------------------------------------------

    def to_command(self, fitskirt_path, mpi_command, scheduler, bind_to_cores=False, threads_per_core=1, to_string=False, remote=None):

        """
        This function ...
        :param fitskirt_path:
        :param mpi_command:
        :param scheduler:
        :param bind_to_cores:
        :param threads_per_core:
        :param to_string:
        :param remote:
        :return:
        """

        # Create the argument list
        arguments = fitskirt_command(fitskirt_path, mpi_command, bind_to_cores, self.parallel.processes, self.parallel.threads, threads_per_core, scheduler, remote)

        # # Parallelization options
        if self.parallel.threads > 0: arguments += ["-t", str(self.parallel.threads)]
        if self.parallel.simulations > 1 and self.parallel.processes <= 1: arguments += ["-s", str(self.parallel.simulations)]
        if self.parallel.dataparallel and self.parallel.processes > 1: arguments += ["-d"]

        # Options for input and output
        if self.input_path is not None: arguments += ["-i", self.input_path]
        if self.output_path is not None: arguments += ["-o", self.output_path]

        # Ski file pattern
        if self.relative: arguments += ["-k"]
        if isinstance(self.fski_path, basestring): arguments += [self.fski_path]
        elif isinstance(self.fski_path, list): arguments += self.fski_path
        else: raise ValueError("The ski pattern must consist of either a string or a list of strings")

        # If requested, convert the argument list into a string
        if to_string:

            # Create the final command string for this simulation
            command = " ".join(arguments)

            # Return the command string
            return command

        # Otherwise, return the list of argument values
        else: return arguments

# -----------------------------------------------------------------

def fitskirt_command(fitskirt_path, mpi_command, bind_to_cores, processes, threads, threads_per_core, scheduler, remote=None):

    """
    This function ...
    :param fitskirt_path:
    :param mpi_command:
    :param bind_to_cores:
    :param processes:
    :param threads:
    :param threads_per_core:
    :param scheduler:
    :param remote:
    :return:
    """

    # Multiprocessing mode
    if processes > 1:

        # Determine the command based on whether or not a scheduling system is used
        if scheduler: command = mpi_command.split()
        else: command = mpi_command.split() + ["-np", str(processes)]

        # If 'process to core' binding must be enabled, add the 'cpus-per-proc' option
        # (see https://www.open-mpi.org/faq/?category=tuning)
        if bind_to_cores:

            # Hyperthreading: threads_per_core will be > 1
            # No hyperthreading: threads_per_core will be 1
            # cores / process = (cores / thread) * (threads / process)
            cores_per_process = int(threads / threads_per_core)

            # Check if cpus-per-proc option is possible
            if remote is None or remote.mpi_knows_cpus_per_proc_option:
                command += ["--cpus-per-proc", str(cores_per_process)] # "CPU'S per process" means "core per process" in our definitions
            else: log.warning("The MPI version on the remote host does not know the 'cpus-per-proc' command. Processes cannot be bound to cores")

        # Add the SKIRT path and return the final command list
        command += [fitskirt_path]
        return command

    # Singleprocessing mode, no MPI command or options
    else: return [fitskirt_path]

# -----------------------------------------------------------------

class FskiFile(object):

    """
    This class ...
    """

    def __init__(self, path):

        """
        This function ...
        :param path:
        """

        # Check the file path
        if not path.lower().endswith(".fski"): raise ValueError("Invalid filename extension for ski file")

        # Set the path to the fski file
        self.path = fs.absolute_path(path)

        # load the XML tree (remove blank text to avoid confusing the pretty printer when saving)
        self.tree = etree.parse(arch.opentext(self.path), parser=etree.XMLParser(remove_blank_text=True))

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function converts the tree into a string
        :return:
        """

        return etree.tostring(self.tree, encoding="UTF-8", xml_declaration=True, pretty_print=False, with_tail=False)

    # -----------------------------------------------------------------

    def saveto(self, filepath, update_path=True):

        """
        This function saves the (possibly updated) contents of the SkiFile instance into the specified file.
        The filename \em must end with ".ski". Saving to and thus replacing the ski file from which this
        SkiFile instance was originally constructed is allowed, but often not the intention.
        :param filepath:
        :param update_path:
        :return:
        """

        if not filepath.lower().endswith(".fski"): raise ValueError("Invalid filename extension for ski file")

        # Update the producer and time attributes on the root element
        root = self.tree.getroot()
        root.set("producer", "Python Toolkit for SKIRT (SkiFile class)")
        root.set("time", datetime.now().strftime("%Y-%m-%dT%H:%M:%S"))

        # Serialize the XML tree
        outfile = open(os.path.expanduser(filepath), "wb")
        outfile.write(etree.tostring(self.tree, encoding="UTF-8", xml_declaration=True, pretty_print=True))
        outfile.close()

        # Update the fski file path
        if update_path: self.path = filepath

    # -----------------------------------------------------------------

    def save(self):

        """
        This function saves the ski file to the original path
        :return:
        """

        self.saveto(self.path)

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function returns a copy (a deep copy) of this ski file
        :return:
        """

        ski = copy.deepcopy(self)
        ski.path = None # set the path to None so this copy won't be involuntarily saved over the original file
        return ski

    # -----------------------------------------------------------------

    def to_string(self):

        """
        This function returns the ski contents as a string
        :return:
        """

        return etree.tostring(self.tree)

    # -----------------------------------------------------------------

    def to_lines(self):

        """
        This function returns the ski contents as a list of strings (lines)
        :return:
        """

        return etree.tostringlist(self.tree)

    # -----------------------------------------------------------------

    def get_fitscheme(self):

        """
        This function returns the FitScheme element
        :return:
        """

        return self.tree.getroot().getchildren()[0]

    # -----------------------------------------------------------------

    def get_simulation(self):

        """
        This function ...
        :return:
        """

        return self.get_unique_base_element("simulation")

    # -----------------------------------------------------------------

    @property
    def ski_name(self):

        """
        This function ...
        :return:
        """

        # Get the adjustable simulation
        simulation = self.get_simulation()

        # Return the value
        return self.get_value(simulation, "skiName")

    # -----------------------------------------------------------------

    def set_ski_name(self, filename):

        """
        This function ...
        :param filename:
        :return:
        """

        # Get the adjustable simulation
        simulation = self.get_simulation()

        # Set the ski name
        self.set_value(simulation, "skiName", filename)

    # -----------------------------------------------------------------

    def get_parameter_ranges(self):

        """
        This function ...
        :return:
        """

        return self.get_unique_base_element("parameterRanges")

    # -----------------------------------------------------------------

    def get_ranges(self):

        """
        This function ...
        :return:
        """

        parameter_ranges = self.get_parameter_ranges()
        range_elements = self.get_children(self.get_child_with_name(parameter_ranges, "ranges"))

        # Create a dictionary for the parameter ranges
        ranges = dict()

        # Loop over the elements
        for element in range_elements:

            # Get properties
            name = self.get_value(element, "label")
            ptype = self.get_value(element, "quantityType")
            pmin = self.get_value(element, "minimumValue")
            pmax = self.get_value(element, "maximumValue")

            # Determine the actual parsing type for the range
            if ptype == "dimless": range_type = "real_range"
            else:
                quantity_type = skirt_quantities_to_pts_quantities[ptype]
                range_type = quantity_type + "_range"

            # Create a string to parse the range
            range_string = pmin + ">" + pmax

            # Parse the range
            parsing_function = getattr(parsing, range_type)
            the_range = parsing_function(range_string)

            # Add the range to the dictionary
            ranges[name] = the_range

        # Return the ranges dictionary
        return ranges

    # -----------------------------------------------------------------

    def get_reference_images(self):

        """
        This function ...
        :return:
        """

        return self.get_unique_base_element("referenceImages")

    # -----------------------------------------------------------------

    def remove_reference_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the image
        component = self.get_image(name)

        # Get the parent
        parent = component.getparent()

        # Remove the reference image
        parent.remove(component)

    # -----------------------------------------------------------------

    def remove_all_reference_images(self):

        """
        This function ...
        :return:
        """

        # Loop over the reference images
        for name in self.get_image_names(): self.remove_reference_image(name)

    # -----------------------------------------------------------------

    def add_reference_image(self, filename, luminosity_ranges, kernel_fwhm, kernel_type="Gaussian", kernel_dimension=6):

        """
        This function ...
        :param filename:
        :param luminosity_ranges:
        :param kernel_fwhm:
        :param kernel_type:
        :param kernel_dimension:
        :return:
        """

        parent = self.get_child_with_name(self.get_reference_images(), "images")

        # Set min and max
        min_luminosities = [range.min for range in luminosity_ranges]
        max_luminosities = [range.max for range in luminosity_ranges]

        # Attributes
        attrs = {"filename": filename, "minLuminosities": stringify.stringify(min_luminosities)[1], "maxLuminosities": stringify.stringify(max_luminosities)[1]}
        image = parent.makeelement("ReferenceImage", attrs)

        kernel_parent = image.makeelement("kernel", {"type": "ConvolutionKernel"})

        # Add children
        kernel_attrs = {"fwhm": repr(kernel_fwhm), "dimension": str(kernel_dimension)}
        kernel = kernel_parent.makeelement("GaussianKernel", kernel_attrs)

        # Add to parent
        kernel_parent.append(kernel)

        image.append(kernel_parent)

        # Add the image
        parent.append(image)

    # -----------------------------------------------------------------

    def get_images(self):

        """
        This function ...
        :return:
        """

        return self.get_children(self.get_child_with_name(self.get_reference_images(), "images"))

    # -----------------------------------------------------------------

    def get_image_names(self):

        """
        This function ...
        :return:
        """

        names = []

        for image in self.get_images():

            name = self.get_value(image, "filename")
            names.append(name)

        return names

    # -----------------------------------------------------------------

    def get_image(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        for image in self.get_images():

            name = self.get_value(image, "filename")
            if name == image_name: return image

    # -----------------------------------------------------------------

    def get_image_min_luminosities(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        return parsing.real_list(self.get_value(self.get_image(image_name), "minLuminosities"))

    # -----------------------------------------------------------------

    def get_image_max_luminosities(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        return parsing.real_list(self.get_value(self.get_image(image_name), "maxLuminositites"))

    # -----------------------------------------------------------------

    def get_image_luminosity_ranges(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        ranges = []

        for min_value, max_value in zip(self.get_image_min_luminosities(image_name), self.get_image_max_luminosities(image_name)):
            ranges.append(RealRange(min_value, max_value))

        return ranges

    # -----------------------------------------------------------------

    def get_image_kernel(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        return xml.get_unique_element(self.get_image(image_name), "kernel")

    # -----------------------------------------------------------------

    def get_image_kernel_type(self, image_name):

        """
        This function ...
        :return:
        """

        return self.get_image_kernel(image_name).tag

    # -----------------------------------------------------------------

    def get_image_fwhm(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        return float(self.get_value(self.get_image_kernel(image_name), "fwhm"))

    # -----------------------------------------------------------------

    def set_image_fwhm(self, image_name, value):

        """
        This function ...
        :param image_name:
        :param value:
        :return:
        """

        self.set_value(self.get_image_kernel(image_name), "fwhm", repr(value))

    # -----------------------------------------------------------------

    def get_image_dimension(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        return int(self.get_value(self.get_image_kernel(image_name), "dimension"))

    # -----------------------------------------------------------------

    def set_image_dimension(self, image_name, value):

        """
        This function ...
        :param image_name:
        :param value:
        :return:
        """

        self.set_value(self.get_image_kernel(image_name), "dimension", str(value))

    # -----------------------------------------------------------------

    def get_optimization(self):

        """
        This function ...
        :return:
        """

        return self.get_unique_base_element("optim")

    # -----------------------------------------------------------------

    def get_population_size(self):

        """
        This function ...
        :return:
        """

        return int(self.get_value(self.get_optimization(), "popsize"))

    # -----------------------------------------------------------------

    def set_population_size(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.set_value(self.get_optimization(), "popsize", str(value))

    # -----------------------------------------------------------------

    def get_ngenerations(self):

        """
        This function ...
        :return:
        """

        return int(self.get_value(self.get_optimization(), "generations"))

    # -----------------------------------------------------------------

    def set_ngenerations(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.set_value(self.get_optimization(), "generations", str(value))

    # -----------------------------------------------------------------

    def get_mutation_rate(self):

        """
        This function ...
        :return:
        """

        return float(self.get_value(self.get_optimization(), "pmut"))

    # -----------------------------------------------------------------

    def set_mutation_rate(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.set_value(self.get_optimization(), "pmut", repr(value))

    # -----------------------------------------------------------------

    def get_crossover_rate(self):

        """
        This function ...
        :return:
        """

        return float(self.get_value(self.get_optimization(), "pcross"))

    # -----------------------------------------------------------------

    def set_crossover_rate(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.set_value(self.get_optimization(), "pcross", repr(value))

    # -----------------------------------------------------------------

    def get_child_with_name(self, element, child_name):

        """
        This function returns the child of the passed element with a given name
        :param element:
        :param child_name:
        :return:
        """

        for child in element.getchildren():
            if child.tag == child_name: return child
        return None

    # -----------------------------------------------------------------

    def get_children(self, element):

        """
        This function ...
        :param element:
        :return:
        """

        return element.getchildren()

    # -----------------------------------------------------------------

    def get_unique_base_element(self, name):

        """
        This function returns the xml tree element with the specified name that is at the base level of the simulation hierarchy
        :param name:
        :return:
        """

        return xml.get_unique_element(self.get_fitscheme(), "//" + name)

    # -----------------------------------------------------------------

    def get_unique_base_element_direct(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return xml.get_unique_element_direct(self.get_fitscheme(), "//" + name)

    # -----------------------------------------------------------------

    def get_value(self, element, name):

        """
        This function returns the string value of a property of an element
        :param element:
        :param name:
        :return:
        """

        if name not in element.attrib: raise ValueError("A property '" + name + "' does not exist for this element")
        return element.get(name)

    # -----------------------------------------------------------------

    def set_value(self, element, name, string):

        """
        This function ...
        :param element:
        :param name:
        :param string:
        :return:
        """

        element.set(name, string)

# -----------------------------------------------------------------

class FitSKIRT(object):

    """
    This function ...
    """

    def __init__(self, path=None, mpi_style="generic"):

        """
        This function ...
        :param path:
        :param mpi_style:
        """

        # Set the FitSKIRT path
        self.path = path if path is not None else ""

        if not self.path.endswith("skirt"): self.path = fs.join(self.path, "fitskirt")
        if self.path != "fitskirt": self.path = os.path.realpath(os.path.expanduser(self.path))

        if self.path == "fitskirt":
            if introspection.skirt_is_present(): self.path = introspection.fitskirt_path
            else: raise EnvironmentError("FitSKIRT is not installed or not in the PATH environment variable")

        # Indicate we are not running yet
        self._process = None

        # Set the MPI style
        self.mpi_style = mpi_style.lower()

    # -----------------------------------------------------------------

    def run(self, definition, parallelization, wait=True, silent=False):

        """
        This function ...
        :param definition:
        :param parallelization:
        :param wait:
        :param silent:
        :return:
        """

        # Create the arguments
        arguments = FitSKIRTArguments(definition, parallelization)

        # Check whether MPI is present on this system if multiple processe are requested
        if arguments.parallel.processes > 1 and not introspection.has_mpi():
            log.warning("No mpirun executable: not running")
            return []

        # Determine the MPI command
        if self.mpi_style == "lsf":
            scheduler = True
            mpi_command = "mpirun -lsf"
        elif self.mpi_style == "generic":
            scheduler = False
            mpi_command = "mpirun"
        else: raise ValueError("Invalid MPI style")

        # Get the command string
        command = arguments.to_command(self.path, mpi_command, scheduler)

        # Debugging
        command_string = " ".join(command)
        log.debug("The command to launch FitSKIRT is: '" + command_string + "'")

        # Launch the FitSKIRT command
        if wait:
            self._process = None
            if silent: subprocess.call(command, stdout=open(os.devnull,'w'), stderr=open(os.devnull,'w'))
            else: subprocess.call(command)
        else: self._process = subprocess.Popen(command, stdout=open(os.path.devnull, 'w'), stderr=subprocess.STDOUT)

        # Return the list of simulations so that their results can be followed up
        #simulations = arguments.simulations(simulation_names=simulation_names)

        # Check whether FitSKIRT has started
        returncode = self._process.poll() if self._process is not None else None
        if wait or returncode is not None: # when wait=True, or returncode is not None, FitSKIRT executable should have finished

            # Check presence of log file
            #if not fs.is_file(simulations.logfilepath()): raise RuntimeError("SKIRT executable has stopped but log file is not present")

            # Check presence of output files
            if fs.is_empty(definition.output_path): raise RuntimeError("FitSKIRT executable has stopped but no output present")

        # Return the list of simulations
        #return simulations

# -----------------------------------------------------------------

class FitSKIRTLauncher(Configurable):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(FitSKIRTLauncher, self).__init__(config, interactive)

        # The FitSKIRT execution context
        self.fitskirt = FitSKIRT()

        # The definition
        self.definition = None

        # The parallelization scheme
        self.parallelization = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the simulation definition (if necessary)
        if self.definition is None: self.create_definition()

        # 2. Set the parallelization scheme
        if not self.has_parallelization: self.set_parallelization()
        else: self.check_parallelization()

        # 3. Launch the simulation
        self.launch()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    @property
    def has_parallelization(self):

        """
        This function ...
        :return:
        """

        return False

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(FitSKIRTLauncher, self).setup(**kwargs)

        # Setup the remote execution context
        #if self.config.remote is not None:
        #    self.remote = SkirtRemote()
        #    self.remote.setup(self.config.remote, self.config.cluster)

        # Create output directory
        if self.config.create_output and not fs.is_directory(self.config.output): fs.create_directory(self.config.output)

    # -----------------------------------------------------------------

    def create_definition(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the FitSKIRT definition ...")

        # Create the simulation definition
        self.definition = FitSKIRTDefinition(self.config.ski, self.config.fski, self.config.output, self.config.input)

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme ...")

        # Get the number of threads per core
        threads_per_core = par.nthreads_per_core()

        # Check whether MPI is available on this system
        if introspection.has_mpi():

            # Pure MPI
            processes = par.ncores() * threads_per_core
            threads = 1
            threads_per_core = 1

        # No MPI available
        else:

            processes = 1
            threads = min(int(math.ceil(monitoring.free_cpus())), 1)

        # Debugging
        log.debug("The number of thread per core is " + str(threads_per_core))
        log.debug("The number of processes is " + str(processes))

        # Set the parallelization scheme
        self.parallelization = Parallelization.from_processes_and_threads(processes, threads, threads_per_core=threads_per_core)

    # -----------------------------------------------------------------

    def check_parallelization(self):

        """
        This function checks whether the parallelization scheme that is asked by the user is possible given the
        number of cores and hyperthreads per core on the remote host.
        Returns:
        """

        # Inform the user
        log.info("Checking the parallelization scheme ...")

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching FitSKIRT ...")

        # Launch remotely or locally
        if self.config.remote is not None: self.launch_remote()
        else: self.launch_local()

    # -----------------------------------------------------------------

    def launch_remote(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def launch_local(self):

        """
        This function ...
        :return:
        """

        # INform the user
        log.info("Launching FitSKIRT locally ...")

        # Run FitSKIRT
        self.fitskirt.run(self.definition, silent=False, wait=True, parallelization=self.parallelization)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
