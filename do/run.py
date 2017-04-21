#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.run Run a PTS command.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import importlib

# Import the relevant PTS classes and modules
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from .commandline import start_and_clear

# -----------------------------------------------------------------

def run_locally(command_name, module_path, class_name, config, input_files, output_files, output, log):

    """
    This function ...
    :param command_name: 
    :param module_path:
    :param class_name:
    :param config:
    :param input_files:
    :param output_files:
    :param output:
    :param log:
    :return: 
    """

    # Start message
    log.start("Starting " + command_name + " ...")

    # Get the class
    cls = introspection.get_class(module_path, class_name)

    # Create the class instance, configure it with the configuration settings
    inst = cls(config)

    # Set input files
    input_dict = {}
    if input_files is not None:

        # Loop over the names of the input variables
        for name in input_files:

            # Get class path
            classpath, filepath = input_files[name]
            modulepath, classname = classpath.rsplit(".", 1)

            # Get input class
            input_module = importlib.import_module(modulepath)
            input_class = getattr(input_module, classname)

            # Open the input file
            input_object = input_class(filepath)

            # Set to input dict
            input_dict[name] = input_object

    # Start
    start_and_clear(command_name, inst.run, **input_dict)

    # Write output files
    if output_files is not None:

        types = dict()

        # Loop over the names of the attributes for output
        for name in output_files:

            # Get filepath
            filepath = output_files[name]

            # Get the output object
            output_object = getattr(inst, name)

            # Set the type
            types[name] = type(output_object).__module__ + "." + type(output_object).__name__

            # Save the output object
            real_filepath = filepath + "." + type(output_object).default_extension
            output_object.saveto(real_filepath)

        # Save the types
        from pts.core.tools import serialization
        types_path = fs.join(output, "types.dat")
        serialization.write_dict(types, types_path)

# -----------------------------------------------------------------

def run_remotely(command_name, config, keep, host_id, log):

    """
    This function ...
    :param command_name:
    :param config:
    :param keep:
    :param host_id:
    :param log:
    :return: 
    """

    # Additional imports
    from pts.core.remote.remote import Remote

    # Start message
    log.start("Starting " + command_name + " on remote host " + host_id + " ...")

    # Debugging
    log.debug("Initializing the remote ...")

    # Initialize the remote execution environment
    remote = Remote(host_id=host_id)

    # Run PTS remotely
    task = remote.run_pts(command_name, config, keep_remote_output=keep)

    # Succesfully submitted
    log.success("Succesfully submitted the PTS job to the remote host")

# -----------------------------------------------------------------
