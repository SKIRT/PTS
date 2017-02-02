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
import tempfile

# Import the relevant PTS classes and modules
from ..remote.remote import Remote
from ..tools import introspection
from ..tools.logging import log
from ..basics.configuration import ConfigurationDefinition, DictConfigurationSetter
from ..tools import filesystem as fs

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

    def run_attached(self, pts_command, config_dict, input_dict=None, return_output_names=None, unpack=False, return_config=False):

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

        # START REMOTE PYTHON SESSION
        python = self.remote.start_python_session()

        # Import the class from which to make an instance
        python.import_package("importlib")
        python.send_line("module = importlib.import_module('" + class_module_path + "')")  # get the module of the class
        python.send_line("cls = getattr(module, '" + class_name + "')")  # get the class

        # Inform the user
        log.start("Starting " + exact_command_name + " ...")

        # Create a remote temporary directory (for the config and input)
        remote_temp_path = self.remote.new_temp_directory()

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
        self.remote.upload(temp_conf_path, remote_temp_path)

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
            self.remote.upload(local_input_filepaths, remote_temp_path)

            ### LOAD THE INPUT DICT REMOTELY

            python.send_python_line("input_dict = dict()")

            for name in input_dict:

                # Determine the remote filepath
                remote_filepath = fs.join(remote_temp_path, name + "." + input_dict[name].default_extension)

                # Import the class of the filetype remotely
                classpath = str(type(input_dict[name])).split("'")[1].split("'")[0]

                modulepath, classname = classpath.rsplit(".", 1)

                python.send_line("input_module = importlib.import_module('" + modulepath + "')")  # get the module of the class
                python.send_line("input_cls = getattr(input_module, '" + classname + "')")  # get the class

                # Open the input file
                python.send_line("input_dict['" + name + "'] = input_cls.from_file('" + remote_filepath + "')", show_output=True)
        ###

        # Import the Configuration class remotely
        python.import_package("Configuration", from_name="pts.core.basics.configuration")

        # Load the config into the remote python session
        python.send_python_line("config = Configuration.from_file('" + remote_conf_path + "')")

        # Create the class instance, configure it with the configuration settings
        python.send_python_line("inst = cls(config)")

        # Run the instance
        if input_dict is not None: python.send_line("inst.run(**input_dict)", show_output=True, timeout=None) # no timeout, this can take a while
        else: python.send_line("inst.run()", show_output=True, timeout=None) # no timeout, this can take a while

        # Set the output
        output_list = None
        if return_output_names is not None:

            # Initialize output list
            output_list = []

            # Fill in the values in the dict
            for name in return_output_names: output_list.append(python.get_simple_property("inst", name))

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

    def _initialize(self, pts_command, config_dict, input_dict):

        """
        This function ...
        :param pts_command:
        :return:
        """

        # Resolve the PTS command
        subproject, command_name, description, class_name, class_module_path, configuration_module_path, configuration_method = find_match(pts_command)

        ## GET THE CONFIGURATION DEFINITION

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

        ###

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
