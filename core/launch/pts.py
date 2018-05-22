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
from ..basics.log import log
from ..tools import filesystem as fs
from ..basics.task import Task
from ...do.commandline import start_target
from ..basics.configuration import create_configuration_passive
from ..tools import strings
from ..tools import types
from ..basics.configuration import get_definition
from ..basics.map import Map
from ..tools import terminal
from ..tools import conda
from ..tools.stringify import tostr

# -----------------------------------------------------------------

class RemoteInstance(object):

    """
    This class ...
    """

# -----------------------------------------------------------------

def get_definition_for_command(command, cwd=None):

    """
    This function ...
    :param command:
    :param cwd:
    :return:
    """

    # Find possible PTS commands
    scripts = introspection.get_scripts()
    tables = introspection.get_arguments_tables()

    # Find matches
    matches = introspection.find_matches_scripts(command, scripts)
    table_matches = introspection.find_matches_tables(command, tables)

    # No match
    if len(matches) + len(table_matches) == 0: raise ValueError("There is no match with a PTS command for '" + command + "'")

    # If there is a unique match in an existing script, raise error
    elif len(matches) == 1 and len(table_matches) == 0:

        # Get subproject and script name
        subproject, filename = matches[0]
        do_subpath = fs.join(introspection.pts_subproject_dir("do"), subproject)
        filepath = fs.join(do_subpath, filename)

        # Default
        write_config = True

        # Read the lines of the file
        triggered = False
        definition_lines = []
        all_relevant_lines = []
        for line in fs.read_lines(filepath):

            if "ConfigurationDefinition(" in line:
                if not line.startswith("#"): all_relevant_lines.append(line)
                triggered = True
                if "write_config" in line:
                    if line.split("write_config=")[1].startswith("True"): write_config = True
                    elif line.split("write_config=")[1].startswith("False"): write_config = False
                    else: raise IOError("Don't know how to interpret '" + line + "'")
            elif "parse_arguments(" in line: break
            elif triggered:

                if line.strip() == "": continue
                definition_lines.append(line)
                if not line.startswith("#"): all_relevant_lines.append(line)

            elif not line.startswith("#") and line.strip("") != "": all_relevant_lines.append(line)

        # Execute
        namespace = Map()
        introspection.execute_lines(all_relevant_lines, namespace)
        definition = namespace.definition
        definition.write_config = write_config
        definition.add_flag("debug", "enable debug output", letter="d")
        definition.add_flag("brief", "brief output", letter="b")
        definition.add_flag("report", "write a report file")
        # Add config path
        if definition.write_config:
            definition.add_optional("config_path", "directory_path", "directory for the configuration file to be written to (relative to the working directory or absolute) (if None, the output directory is used)")
        return definition

    # If there is an unique match in a table
    elif len(table_matches) == 1 and len(matches) == 0:

        subproject, index = table_matches[0]
        command_name = tables[subproject]["Command"][index]
        description = tables[subproject]["Description"][index]
        class_path_relative = tables[subproject]["Path"][index]
        class_path = "pts." + subproject + "." + class_path_relative
        module_path, class_name = class_path.rsplit('.', 1)
        configuration_method_table = tables[subproject]["Configuration method"][index]

        # Determine configuration module path
        configuration_name = tables[subproject]["Configuration"][index]
        if configuration_name == "--": configuration_name = command_name
        configuration_module_path = "pts." + subproject + ".config." + configuration_name

        # Get the definition
        definition = get_definition(class_name, configuration_module_path, cwd=cwd)

        # Exact command name
        exact_command_name = subproject + "/" + command_name

        # Return the definition
        return definition

    # Show possible matches if there are more than just one
    else: raise ValueError("The PTS command '" + command + "' is ambigious")

# -----------------------------------------------------------------

def execute_pts_local(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    if len(args) == 0: raise ValueError("No arguments were provided")

    # Determine the python path
    python_path = conda.conda_python_path_for_pts()
    pts_main_path = fs.join(introspection.pts_do_dir, "__main__.py")

    # Get options
    show_output = kwargs.pop("show_output", False)
    cwd = kwargs.pop("cwd", None)

    # Create the command
    command = construct_pts_command_string(python_path, pts_main_path, args, kwargs)

    # Execute
    output = terminal.execute(command, show_output=show_output, cwd=cwd)

    # Return the output
    return output

# -----------------------------------------------------------------

def execute_pts_remote(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    if len(args) == 0: raise ValueError("No arguments were provided")

    # Get the remote
    remote = args[0]

    # Get options
    show_output = kwargs.pop("show_output", False)
    cwd = kwargs.pop("cwd", None)

    # Determine python path
    python_path = remote.conda_pts_environment_python_path

    # Create the command
    command = construct_pts_command_string(python_path, remote.pts_main_path, args[1:], kwargs) # strip remote from the args

    # Execute
    output = remote.execute(command, show_output=show_output, cwd=cwd)

    # If errors popped up
    for line in output:
        if "error: unrecognized argument:" in line: raise ValueError("Invalid argument (possibly because of outdated PTS version?): [" + line.split("unrecognized argument: ")[1] + "]")
        if "error: unrecognized arguments:" in line: raise ValueError("Invalid arguments (possibly because of outdated PTS version?): [" + line.split("unrecognized arguments: ")[1] + "]")
        if "error: argument" in line and "invalid" in line: raise ValueError("Invalid argument value: " + line.split("value: ")[1])
    if "Error:" in output[-1]: raise RuntimeError(output[-1])

    # Return the output lines
    return output

# -----------------------------------------------------------------

def construct_pts_command_string(python_path, pts_main_path, args, kwargs):

    """
    Thisf unction ...
    :return:
    """

    # Get PTS command
    if len(args) == 0: raise ValueError("No command specified")
    command = args[0]
    if len(args) > 1: positional = args[1:]
    else: positional = []

    # Create string of positional arguments
    positional_string = " ".join([strings.add_quotes_if_spaces(tostr(item)) for item in positional])

    _definition = None

    # Get optional arguments
    optional_parts = []
    for name in kwargs:

        original_name = name

        # A flag option
        if types.is_boolean_type(kwargs[name]):

            # Get the definition if not yet done
            if _definition is None: _definition = get_definition_for_command(command)

            # Check
            if name not in _definition.flags: raise ValueError("There is no '" + name + "' flag in the configuration definition")

            if kwargs[name]: # value of True

                # If default is True, this option should not be specified on the command line
                if _definition.flags[name].default: continue

                # If the default is False, the option should be added to the command line
                else: pass # (name = name)

            else:  # value of False

                # If default is True, the option should be added to the command line, but with not_
                if _definition.flags[name].default: name = "not_" + name

                # If the default is False, this option should not be specified on the command line
                else: continue

        # Set command-line name
        if strings.is_character(name): option_name = "-" + name
        else: option_name = "--" + name

        if types.is_boolean_type(kwargs[original_name]): optional_parts.append(option_name)
        else:
            #value = strings.add_quotes_if_spaces(str(kwargs[name]))
            value = strings.add_quotes_if_spaces(tostr(kwargs[original_name]))
            optional_parts.append(option_name + " " + value)

    # Create string of optional arguments
    optional_string = " ".join(optional_parts)

    # Create command string
    command = python_path + " " + pts_main_path + " " + command + " " + positional_string + " " + optional_string

    # Return the command
    return command

# -----------------------------------------------------------------

def load_task(host_id, task_id):

    """
    This function ...
    :param host_id: 
    :param task_id: 
    :return: 
    """

    # Determine path
    host_id_run_path = fs.join(introspection.pts_run_dir, host_id)
    task_path = fs.join(host_id_run_path, str(task_id) + ".task")

    # Load the task
    task = Task.from_file(task_path)

    # Return the task
    return task

# -----------------------------------------------------------------

def launch_local(pts_command, config_dict, input_dict=None, analysers=None, analysis_info=None, description=None,
                 cwd=None, debug=False, brief=False):

    """
    This function ...
    :param pts_command:
    :param config_dict:
    :param input_dict:
    :param analysers:
    :param analysis_info:
    :param description:
    :param cwd:
    :param debug:
    :param brief:
    :return:
    """

    # Resolve the PTS command
    subproject, command_name, command_description, class_name, class_module_path, configuration_module_path, configuration_method = find_match(pts_command)

    # Set description
    if description is None: description = command_description

    # Create the configuration
    config = create_configuration_passive(command_name, class_name, configuration_module_path, config_dict, description, cwd=cwd)

    # Set logging options
    if debug: config["debug"] = True
    if brief: config["brief"] = True

    # Set the log level temporarily
    previous_level = log.level
    level = "INFO"
    if config.debug: level = "DEBUG"
    if config.brief: level = "SUCCESS"
    log.setLevel(level)

    # Get the class
    cls = introspection.get_class(class_module_path, class_name)

    # Change working directory
    if cwd is not None: previous_cwd = fs.change_cwd(cwd)
    else: previous_cwd = None

    # Create the class instance, configure it with the configuration settings
    inst = cls(config)

    # Start
    if input_dict is None: input_dict = {}
    start_target(command_name, inst.run, **input_dict)

    # Exact command name
    exact_command_name = subproject + "/" + command_name

    # If analysis has to be performed, create task and analyse
    if analysers is not None:

        # Create task
        task = Task(exact_command_name, config.to_string())

        # Set analysis info
        task.analysis_info = analysis_info

        # Set analysers
        for analyser_class in analysers: task.add_analyser(analyser_class)

        # Do the analysis
        task.analyse()

    # Reset the log level
    log.setLevel(previous_level)

    # Move back to the previous working directory
    if previous_cwd is not None: fs.change_cwd(previous_cwd)

    # Return the instance
    return inst

# -----------------------------------------------------------------

def launch_remote(remote, pts_command, config_dict, input_dict=None, analysers=None, analysis_info=None,
                  return_output_names=None, unpack=False, use_session=False, local_output_path=None):

    """
    This function ...
    :param remote:
    :param pts_command:
    :param config_dict:
    :param input_dict:
    :param analysers:
    :param analysis_info:
    :param return_output_names:
    :param unpack:
    :param use_session:
    :param local_output_path:
    :return:
    """

    # Create the PTS remote launcher
    launcher = PTSRemoteLauncher(remote=remote)

    # Run
    return launcher.run_and_analyse(pts_command, config_dict, local_output_path, analysers, analysis_info,
                                    use_session=use_session, input_dict=input_dict, return_output_names=return_output_names,
                                    unpack=unpack)

# -----------------------------------------------------------------

class PTSRemoteLauncher(object):

    """
    This class ...
    """

    def __init__(self, remote=None):

        """
        The constructor ...
        :param remote:
        :return:
        """

        # -- Attributes --

        # Create the remote execution context
        if remote is not None: self.remote = remote
        else: self.remote = Remote()

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

    def run_and_analyse(self, pts_command, config_dict, local_output_path, analysers, analysis_info, use_session=True,
                        input_dict=None, return_output_names=None, unpack=False):

        """
        This function ...
        :param pts_command:
        :param config_dict:
        :param local_output_path:
        :param analysers:
        :param analysis_info:
        :param use_session:
        :param input_dict:
        :param return_output_names:
        :param unpack:
        :return:
        """

        # Run attached
        if return_output_names is not None: output, config = self.run_attached(pts_command, config_dict, return_config=True, use_session=use_session,
                                                            input_dict=input_dict, return_output_names=return_output_names,
                                                                               unpack=unpack)
        else:
            config = self.run_attached(pts_command, config_dict, return_config=True, use_session=use_session, input_dict=input_dict)
            output = None

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

        # Return output
        if return_output_names is not None: return output

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
        python = self.remote.start_python_session(output_path=remote_temp_path, attached=True)

        # Import the class from which to make an instance
        python.import_package("importlib", show_output=log.is_debug)
        python.send_line("module = importlib.import_module('" + class_module_path + "')", show_output=log.is_debug)  # get the module of the class
        python.send_line("cls = getattr(module, '" + class_name + "')", show_output=log.is_debug)  # get the class

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
        self.remote.upload_retry(temp_conf_path, remote_temp_path, show_output=log.is_debug)

        # Remove the original config file
        fs.remove_file(temp_conf_path)

        ##### UPLOADING THE INPUT #####

        # If input is given
        if input_dict is not None:

            # Debugging
            log.debug("Uploading the input ...")

            # Load the dictionary in the remote session
            python.load_dictionary("input_dict", input_dict, local_temp_path=temp_path, remote_temp_path=remote_temp_path)

        ###

        # Import the Configuration class remotely
        python.import_package("Configuration", from_name="pts.core.basics.configuration", show_output=log.is_debug)

        # Load the config into the remote python session
        python.send_line("config = Configuration.from_file('" + remote_conf_path + "')", show_output=log.is_debug)

        # Create the class instance, configure it with the configuration settings
        python.send_line("inst = cls(config)", show_output=log.is_debug)

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
            for name in return_output_names: output_list.append(python.get_simple_property("inst", name, show_output=log.is_debug))

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
        self.remote.upload_retry(temp_conf_path, remote_temp_path, show_output=log.is_debug)

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
            self.remote.upload_retry(local_input_filepaths, remote_input_path, show_output=log.is_debug)

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
        command_string = "pts --configfile " + remote_conf_path
        if remote_input_path is not None: command_string += " --input " + remote_input_path
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
        :param config_dict:
        :param input_dict:
        :return:
        """

        # Resolve the PTS command
        subproject, command_name, description, class_name, class_module_path, configuration_module_path, configuration_method = find_match(pts_command)

        # Create the configuration
        config = create_configuration_passive(command_name, class_name, configuration_module_path, config_dict, description)

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
    if len(matches) + len(table_matches) == 0: raise ValueError("There is no match with a PTS command for '" + pts_command + "'")

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
