#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.pts Contains the PTSRemoteLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import importlib
import exceptions
import tempfile

# Import the relevant PTS classes and modules
from ..remote.remote import Remote
from ..tools import introspection
from ..tools.logging import log
from ..basics.configuration import ConfigurationDefinition, DictConfigurationSetter
from ..tools import filesystem as fs
from ..basics.task import Task
from ...do.commandline import start_target

# -----------------------------------------------------------------

def launch_local(pts_command, config_dict, input_dict=None, analysers=None, analysis_info=None):

    """
    This function ...
    :param pts_command:
    :param config_dict:
    :param input_dict:
    :param analysers:
    :param analysis_info:
    :return:
    """

    # Resolve the PTS command
    subproject, command_name, description, class_name, class_module_path, configuration_module_path, configuration_method = find_match(pts_command)

    # Create the configuration
    config = create_configuration(command_name, class_name, configuration_module_path, config_dict, input_dict, description)

    # Get the class
    cls = introspection.get_class(class_module_path, class_name)

    # Create the class instance, configure it with the configuration settings
    inst = cls(config)

    # Start
    start_target(command_name, inst.run)

    # Exact command name
    exact_command_name = subproject + "/" + command_name

    # Create task
    task = Task(exact_command_name, config.to_string())

    # Set analysis info
    task.analysis_info = analysis_info

    # Set analysers
    for analyser_class in analysers: task.add_analyser(analyser_class)

    # Do the analysis
    task.analyse()

# -----------------------------------------------------------------

class PTSRemoteLauncher(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # -- Attributes --

        # Create the remote execution context
        self.remote = Remote()

    # -----------------------------------------------------------------

    def setup(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Setup the remote
        self.remote.setup(host_id)

    # -----------------------------------------------------------------

    @property
    def host(self):

        """
        This function ...
        :return:
        """

        return self.remote.host

    # -----------------------------------------------------------------

    @property
    def host_id(self):

        """
        This function ...
        :return:
        """

        return self.remote.host_id

    # -----------------------------------------------------------------

    def run_detached(self, pts_command, config_dict, input_dict=None, analysers=None, analysis_info=None,
                     keep_remote_output=False, remove_local_output=False):

        """
        This function ...
        :param pts_command:
        :param config_dict:
        :param input_dict:
        :param analysers:
        :param analysis_info:
        :param keep_remote_output:
        :param remove_local_output:
        :return:
        """

        # Initialize
        subproject, exact_command_name, class_name, class_module_path, config = self._initialize(pts_command, config_dict, input_dict)

        # Debugging
        log.debug("Running in detached mode ...")

        # Run PTS remotely
        task = self.remote.run_pts(exact_command_name, config, input_dict=input_dict, keep_remote_output=keep_remote_output, remove_local_output=remove_local_output)

        # Set the analysers
        if analysers is not None:
            for analyser in analysers: task.add_analyser(analyser)

        # Set the analysis info
        if analysis_info is not None: task.analysis_info = analysis_info

        # Save the task
        task.save()

        # Succesfully submitted
        log.success("Succesfully submitted the PTS job to the remote host")

    # -----------------------------------------------------------------

    def run_and_analyse(self, pts_command, config_dict, local_output_path, analysers, analysis_info, use_session=True):

        """
        This function ...
        :param pts_command:
        :param config_dict:
        :param local_output_path:
        :param analysers:
        :param analysis_info:
        :param use_session:
        :return:
        """

        # Run attached
        config = self.run_attached(pts_command, config_dict, return_config=True, use_session=use_session)

        # Create a new Task object
        task = Task(pts_command, config.to_string())

        # Set the host ID and cluster name (if applicable)
        task.host_id = self.host_id
        task.cluster_name = None

        # Generate a new task ID
        task_id = self.remote._new_task_id()

        # Set properties such as the task ID and name and the screen name
        task.id = task_id

        # Set local and remote output path
        task.local_output_path = local_output_path

        # Set analysis info
        task.analysis_info = analysis_info

        # Set analysers
        for analyser_class in analysers: task.add_analyser(analyser_class)

        # Analyse
        task.analyse()

    # -----------------------------------------------------------------

    def run_attached(self, pts_command, config_dict, input_dict=None, return_output_names=None, unpack=False,
                     return_config=False, use_session=True):

        """
        This function ...
        :param pts_command:
        :param config_dict:
        :param input_dict:
        :param return_output_names:
        :param unpack:
        :param return_config:
        :param use_session:
        :return:
        """

        if use_session: return self._run_attached_session(pts_command, config_dict, input_dict, return_output_names, unpack, return_config)
        else: return self._run_attached_no_session(pts_command, config_dict, input_dict, return_output_names, unpack, return_config)

    # -----------------------------------------------------------------

    def _run_attached_session(self, pts_command, config_dict, input_dict=None, return_output_names=None, unpack=False,
                              return_config=False):

        """
        This function ...
        :param pts_command:
        :param config_dict:
        :param input_dict:
        :param return_output_names:
        :param unpack:
        :param return_config:
        :return:
        """

        # Don't write when running in attached mode
        config_dict["write"] = False

        # Initialize
        subproject, exact_command_name, class_name, class_module_path, config = self._initialize(pts_command, config_dict, input_dict)

        # Create a remote temporary directory (for the config and input)
        remote_temp_path = self.remote.new_temp_directory()

        # Debugging
        log.info("Running in attached mode ...")

        # START REMOTE PYTHON SESSION
        #python = self.remote.start_python_session(output_path=remote_temp_path)
        python = self.remote.start_python_session(output_path=remote_temp_path, attached=True)

        # Import the class from which to make an instance
        python.import_package("importlib", show_output=log.is_debug())
        python.send_line("module = importlib.import_module('" + class_module_path + "')", show_output=log.is_debug())  # get the module of the class
        python.send_line("cls = getattr(module, '" + class_name + "')", show_output=log.is_debug())  # get the class

        # Inform the user
        log.start("Starting " + exact_command_name + " ...")

        # Always create a log file while executing remotely
        config.report = True
        config.log_path = remote_temp_path
        config.path = remote_temp_path

        #### UPLOADING THE CONFIG TO THE REMOTE ####

        # Debugging
        log.debug("Saving the configuration file locally ...")

        # Determine path to the temporarily saved local configuration file
        temp_path = tempfile.gettempdir()
        temp_conf_path = fs.join(temp_path, "config.cfg")

        # Save the configuration file to the temporary directory
        config.saveto(temp_conf_path)

        # Debugging
        log.debug("Uploading the configuration file to '" + remote_temp_path + "' ...")

        # Upload the config file
        remote_conf_path = fs.join(remote_temp_path, fs.name(temp_conf_path))
        self.remote.upload_retry(temp_conf_path, remote_temp_path, show_output=log.is_debug())

        # Remove the original config file
        fs.remove_file(temp_conf_path)

        ##### UPLOADING THE INPUT #####

        # If input is given
        if input_dict is not None:

            # Debugging
            log.debug("Uploading the input ...")

            ## Save the input locally

            local_input_filepaths = []

            for name in input_dict:

                # Determine filepath
                path = fs.join(temp_path, name + "." + input_dict[name].default_extension)

                # Save
                input_dict[name].saveto(path)

                # Add the filepath
                local_input_filepaths.append(path)

            #### UPLOAD THE INPUT :

            # Upload the input files
            self.remote.upload_retry(local_input_filepaths, remote_temp_path, show_output=log.is_debug())

            ### LOAD THE INPUT DICT REMOTELY
            python.send_line("input_dict = dict()", show_output=log.is_debug())
            for name in input_dict:

                # Determine the remote filepath
                remote_filepath = fs.join(remote_temp_path, name + "." + input_dict[name].default_extension)

                # Import the class of the filetype remotely
                classpath = str(type(input_dict[name])).split("'")[1].split("'")[0]

                modulepath, classname = classpath.rsplit(".", 1)

                python.send_line("input_module = importlib.import_module('" + modulepath + "')", show_output=log.is_debug())  # get the module of the class
                python.send_line("input_cls = getattr(input_module, '" + classname + "')", show_output=log.is_debug())  # get the class

                # Open the input file
                python.send_line("input_dict['" + name + "'] = input_cls.from_file('" + remote_filepath + "')", show_output=True)
        ###

        # Import the Configuration class remotely
        python.import_package("Configuration", from_name="pts.core.basics.configuration", show_output=log.is_debug())

        # Load the config into the remote python session
        python.send_line("config = Configuration.from_file('" + remote_conf_path + "')", show_output=log.is_debug())

        # Create the class instance, configure it with the configuration settings
        python.send_line("inst = cls(config)", show_output=log.is_debug())

        # Run the instance
        if input_dict is not None: output = python.send_line("inst.run(**input_dict)", show_output=True, timeout=None) # no timeout, this can take a while
        else: output = python.send_line("inst.run()", show_output=True, timeout=None) # no timeout, this can take a while

        # Check the log output
        last_line = output[-1]
        if "Error:" in last_line:
            error_message = last_line.split(": ")[1]
            error_type = last_line.split(":")[0]
            raise getattr(exceptions, error_type)(error_message)

        # Set the output
        output_list = None
        if return_output_names is not None:

            # Initialize output list
            output_list = []

            # Fill in the values in the dict
            for name in return_output_names: output_list.append(python.get_simple_property("inst", name, show_output=log.is_debug()))

        ######

        # End the python session
        del python

        # Return the output(can be None if return_output_names was None)
        if output_list is None:
            if return_config: return config
            else: return None
        if unpack:
            if len(output_list) == 1:
                if return_config: return output_list[0], config
                else: return output_list[0]
            else:
                if return_config: return output_list, config
                else: return output_list
        else:
            if return_config: return dict(zip(return_output_names, output_list)), config
            else: return dict(zip(return_output_names, output_list))

    # -----------------------------------------------------------------

    def _run_attached_no_session(self, pts_command, config_dict, input_dict=None, return_output_names=None, unpack=False,
                                return_config=False):

        """
        This function ...
        :param pts_command:
        :param config_dict:
        :param input_dict:
        :param return_output_names:
        :param unpack:
        :param return_config:
        :return:
        """

        # Don't write when running in attached mode
        config_dict["write"] = False

        # Initialize
        subproject, exact_command_name, class_name, class_module_path, config = self._initialize(pts_command, config_dict, input_dict)

        # Create a remote temporary directory (for the config and input)
        remote_temp_path = self.remote.new_temp_directory()
        remote_output_path = self.remote.create_directory_in(remote_temp_path, "out")

        # Debugging
        log.info("Running in attached mode ...")

        # Always create a log file while executing remotely
        config.report = True
        config.log_path = remote_temp_path
        config.path = remote_temp_path

        #### UPLOADING THE CONFIG TO THE REMOTE ####

        # Debugging
        log.debug("Saving the configuration file locally ...")

        # Determine path to the temporarily saved local configuration file
        temp_path = tempfile.gettempdir()
        temp_conf_path = fs.join(temp_path, "config.cfg")

        # Save the configuration file to the temporary directory
        config.saveto(temp_conf_path)

        # Debugging
        log.debug("Uploading the configuration file to '" + remote_temp_path + "' ...")

        # Upload the config file
        remote_conf_path = fs.join(remote_temp_path, fs.name(temp_conf_path))
        self.remote.upload_retry(temp_conf_path, remote_temp_path, show_output=log.is_debug())

        # Remove the original config file
        fs.remove_file(temp_conf_path)

        ##### UPLOADING THE INPUT #####

        # If input is given
        if input_dict is not None:

            remote_input_path = self.remote.create_directory_in(remote_temp_path, "in")

            # Debugging
            log.debug("Uploading the input ...")

            ## Save the input locally

            local_input_filepaths = []

            for name in input_dict:

                # Determine filepath
                path = fs.join(temp_path, name + "." + input_dict[name].default_extension)

                # Save
                input_dict[name].saveto(path)

                # Add the filepath
                local_input_filepaths.append(path)

            #### UPLOAD THE INPUT :

            # Upload the input files
            self.remote.upload_retry(local_input_filepaths, remote_input_path, show_output=log.is_debug())

            # Set remote input paths
            remote_input_paths = dict()
            remote_input_classes = dict()

            for name in input_dict:

                # Determine the remote filepath
                remote_filepath = fs.join(remote_input_path, name + "." + input_dict[name].default_extension)

                # Import the class of the filetype remotely
                classpath = str(type(input_dict[name])).split("'")[1].split("'")[0]
                #modulepath, classname = classpath.rsplit(".", 1)

                # Set the file path
                remote_input_paths[name] = remote_filepath

                # Set the file class path
                remote_input_classes[name] = classpath

        # No input required
        else: remote_input_path = remote_input_paths = remote_input_classes = None

        # If output is requested
        if return_output_names is not None:

            remote_output_path = self.remote.create_directory_in(remote_temp_path, "output")

            # Set remote output paths
            remote_output_paths = dict()

            for name in return_output_names:

                # Determine remote filepath (but without extension)
                remote_filepath = fs.join(remote_output_path, name)

                # Set the file path
                remote_output_paths[name] = remote_filepath

        # No output required
        else: remote_output_paths = None

        # Construct the PTS command string
        command_string = "pts --configfile " + remote_conf_path + " --input " + remote_input_path
        command_string += " --output " + remote_output_path + " " + exact_command_name

        # Add input file paths
        if remote_input_paths is not None:

            input_file_paths_parts = []

            for name in remote_input_paths:

                filepath = remote_input_paths[name]
                class_path = remote_input_classes[name]
                part = "'" + name + "':('" + class_path + "','" + filepath + "')"
                input_file_paths_parts.append(part)

            # Construct dictionary string
            input_file_paths_string = "{" + ",".join(input_file_paths_parts) + "}"

            # Add option
            command_string += " --input_files " + '"' + input_file_paths_string + '"'

        # Add output file paths
        if remote_output_paths is not None:

            output_file_paths_parts = []

            for name in remote_output_paths:

                filepath = remote_output_paths[name]
                part = "'" + name + "':'" + filepath +"'"
                output_file_paths_parts.append(part)

            # Construct dictionary string
            output_file_paths_string = "{" + ",".join(output_file_paths_parts) + "}"

            # Add option
            command_string += " --output_files " + '"' + output_file_paths_string + '"'

        # Run the command
        output = self.remote.execute(command_string, show_output=True)

        # Check the log output
        last_line = output[-1]
        if "Error:" in last_line:
            error_message = last_line.split(": ")[1]
            error_type = last_line.split(":")[0]
            raise getattr(exceptions, error_type)(error_message)

        # Set the output
        output_list = None
        if return_output_names is not None:

            # Initialize output list
            output_list = []

            # Open the types.dat file
            types_file_path = fs.join(remote_output_path, "types.dat")
            from ..tools import serialization
            types = serialization.load_dict(types_file_path)

            # Fill in the values in the dict
            for name in return_output_names:

                # Search for the file
                paths = self.remote.files_in_path(remote_output_path, exact_name=name)
                if len(paths) > 0: raise RuntimeError("Encountered multiple files with the name '" + name + "': " + str(paths))
                filepath = paths[0]

                # Get the type
                classpath = types[name]
                modulepath, classname = classpath.rsplit(".", 1)

                # Get output class
                output_module = importlib.import_module(modulepath)
                output_class = getattr(output_module, classname)

                # Open the output object
                output_object = output_class(filepath)

                # Add to list
                output_list.append(output_object)

        # Return the output(can be None if return_output_names was None)
        if output_list is None:
            if return_config: return config
            else: return None
        if unpack:
            if len(output_list) == 1:
                if return_config: return output_list[0], config
                else: return output_list[0]
            else:
                if return_config: return output_list, config
                else: return output_list
        else:
            if return_config: return dict(zip(return_output_names, output_list)), config
            else: return dict(zip(return_output_names, output_list))

    # -----------------------------------------------------------------

    def _initialize(self, pts_command, config_dict, input_dict):

        """
        This function ...
        :param pts_command:
        :return:
        """

        # Resolve the PTS command
        subproject, command_name, description, class_name, class_module_path, configuration_module_path, configuration_method = find_match(pts_command)

        # Create the configuration
        config = create_configuration(command_name, class_name, configuration_module_path, config_dict, input_dict, description)

        # Exact command name
        exact_command_name = subproject + "/" + command_name

        # Inform the user about starting the PTS command remotely
        log.start("Starting " + exact_command_name + " on remote host " + self.remote.host_id + " ...")

        # Return
        return subproject, exact_command_name, class_name, class_module_path, config

# -----------------------------------------------------------------

def find_match(pts_command):

    """
    This function ...
    :param pts_command:
    :return:
    """

    # Find possible PTS commands
    scripts = introspection.get_scripts()
    tables = introspection.get_arguments_tables()

    # Find matches
    matches = introspection.find_matches_scripts(pts_command, scripts)
    table_matches = introspection.find_matches_tables(pts_command, tables)

    # No match
    if len(matches) + len(table_matches) == 0:

        #log.error("There is no match")
        #show_all_available(scripts, tables)
        #exit()
        raise ValueError("There is no match with a PTS command for '" + pts_command + "'")

    # If there is a unique match in an existing script, raise error
    elif len(matches) == 1 and len(table_matches) == 0: raise ValueError("This do command cannot be executed remotely")

    # If there is an unique match in a table
    elif len(table_matches) == 1 and len(matches) == 0:

        subproject, index = table_matches[0]
        command_name = tables[subproject]["Command"][index]
        description = tables[subproject]["Description"][index]
        class_path_relative = tables[subproject]["Path"][index]
        class_path = "pts." + subproject + "." + class_path_relative
        module_path, class_name = class_path.rsplit('.', 1)

        configuration_method_table = tables[subproject]["Configuration method"][index]

    # Show possible matches if there are more than just one
    else: raise ValueError("The PTS command '" + pts_command + "' is ambigious")

    configuration_name = tables[subproject]["Configuration"][index]
    if configuration_name == "--": configuration_name = command_name
    configuration_module_path = "pts." + subproject + ".config." + configuration_name

    # Return
    return subproject, command_name, description, class_name, module_path, configuration_module_path, configuration_method_table

# -----------------------------------------------------------------

def create_configuration(command_name, class_name, configuration_module_path, config_dict, input_dict=None, description=None):

    """
    This function ...
    :param command_name:
    :param class_name:
    :param configuration_module_path:
    :param config_dict:
    :param input_dict:
    :param description:
    :return:
    """

    # Find definition
    try:
        configuration_module = importlib.import_module(configuration_module_path)
        definition = getattr(configuration_module, "definition")
    except ImportError:
        log.warning("No configuration definition found for the " + class_name + " class")
        definition = ConfigurationDefinition()  # Create new configuration definition

    ## SET THE INPUT FILENAMES

    if input_dict is not None:

        config_dict["input"] = dict()
        for name in input_dict:
            config_dict["input"] = name + "." + input_dict[name].default_extension  # generate a default filename

    ## CREATE THE CONFIGURATION

    # Create the configuration setter
    if config_dict is None: config_dict = dict()  # no problem if all options are optional
    setter = DictConfigurationSetter(config_dict, command_name, description)

    # Create the configuration from the definition and from the provided configuration dictionary
    config = setter.run(definition)

    # Return the configuration
    return config

# -----------------------------------------------------------------
