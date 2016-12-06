#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.remote Contains the Remote class.

# -----------------------------------------------------------------

# Import standard modules
import re
import sys
import pexpect
from pexpect import pxssh
import tempfile
import StringIO

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from .host import Host
from .vpn import VPN
from ..tools.logging import log
from ..tools import parsing
from ..tools import filesystem as fs
from ..tools import time
from ..basics.task import Task
from ..tools import introspection

# -----------------------------------------------------------------

connected_remotes = dict()

# -----------------------------------------------------------------

class Remote(object):

    """
    This function ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(Remote, self).__init__()

        # -- Attributes --

        # The SSH interface, an instance of the pxssh class
        self.ssh = pxssh.pxssh()

        # The host instance
        self.host = None

        # The VPN service
        self.vpn = None

        # A flag indicating whether the connection with the remote has been established
        self.connected = False

        # A regular expression object that strips away special unicode characters, used on the remote console output
        self.ansi_escape = re.compile(r'\x1b[^m]*m')

        # Flag that says whether we are in a remote python session
        self.in_python_session = False

    # -----------------------------------------------------------------

    def setup(self, host_id, cluster_name=None):

        """
        This function ...
        :param host_id:
        :param cluster_name:
        :return:
        """

        # Create the host object
        self.host = Host(host_id, cluster_name)

        # If a VPN connection is required for the remote host
        if self.host.requires_vpn: self.connect_to_vpn()

        # Make the connection
        self.login()

        # TODO: swap to cluster here?

        # Load modules
        self.load_modules()

    # -----------------------------------------------------------------

    def __del__(self):

        """
        The destructor ...
        :return:
        """

        # Disconnect from the remote host
        self.logout()

    # -----------------------------------------------------------------

    def load_modules(self):

        """
        This function ...
        :return:
        """

        # Load the necessary modules
        if self.host.modules is not None:

            log.info("Loading necessary modules...")
            self.execute("module load " + " ".join(self.host.modules), output=False)

    # -----------------------------------------------------------------

    def load_installation_modules(self):

        """
        This function ...
        :return:
        """

        # Load the modules
        if self.host.installation_modules is not None:

            log.info("Loading installation modules ...")
            self.execute("module load " + " ".join(self.host.installation_modules), output=False)

    # -----------------------------------------------------------------

    def unload_all_modules(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Unloading all modules ...")

        # Module purge
        self.execute("module purge", output=False)

    # -----------------------------------------------------------------

    def connect_to_vpn(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Connecting to vpn service '" + self.host.vpn.service + "' ...")

        # Connect to the VPN service
        self.vpn = VPN(self.host.vpn.service)
        self.vpn.connect(self.host.vpn.user, self.host.vpn.password, self.host.vpn.secret, self.host.vpn.prompt_time_delay)

    # -----------------------------------------------------------------

    def login(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Logging in to the remote environment on host '" + self.host.id + "' ...")

        # Connect to the remote host
        self.connected = self.ssh.login(self.host.name, self.host.user, self.host.password, port=self.host.port, login_timeout=30)

        # Check whether connection was succesful
        if not self.connected: raise RuntimeError("Connection failed")

        # If the connection was succesful, add the remote to the dictionary of currently connected remotes
        if self.connected: connected_remotes[self.host.id] = self

    # -----------------------------------------------------------------

    def logout(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        if log is not None: log.info("Logging out from the remote environment ...")
        # the conditional statement is because of this error message during destruction at the end of a script:
        # Exception AttributeError: "'NoneType' object has no attribute 'info'" in <bound method SkirtRemote.__del__ of <pts.core.simulation.remote.SkirtRemote object at 0x118628d50>> ignored

        # Disconnect
        if self.connected:

            # End python session if running
            if self.in_python_session: self.end_python_session()

            self.ssh.logout()
            if connected_remotes is not None: del connected_remotes[self.host.id] # the conditional statement is because 'connected_remotes' is set to None during destruction at the end of a script
            self.connected = False

        # Disconnect from the VPN service if necessary
        if self.vpn is not None: self.vpn.disconnect()

    # -----------------------------------------------------------------

    def start_screen(self, name, local_script_path, script_destination, screen_output_path=None,
                     keep_remote_script=False, attached=False):

        """
        This function ...
        :param name:
        :param local_script_path:
        :param script_destination:
        :param screen_output_path:
        :param keep_remote_script:
        :param attached:
        :return:
        """

        # Copy the script to the remote host
        self.upload(local_script_path, script_destination)

        # Rename the remote script
        local_script_name = fs.name(local_script_path)
        remote_script_name = name + ".sh"
        remote_script_path = fs.join(script_destination, remote_script_name)
        self.rename_file(script_destination, local_script_name, remote_script_name)

        # Make the shell script executable
        self.execute("chmod +x " + remote_script_path, output=False)

        # Record the screen output: 'script' command
        #if screen_output_path is not None: self.execute("script " + screen_output_path)

        # If the output directory is given for the screen output
        cwd = None
        if screen_output_path is not None:

            # Get the path to the current working directory
            cwd = self.working_directory

            # Now change the working directory to the screen_output_path
            self.change_cwd(screen_output_path)

        # Create the screen session and execute the batch script
        if attached:
            #self.execute("screen -S " + name + " -L -m " + remote_script_path, output=False, show_output=True)
            self.execute("sh " + remote_script_path, output=False, show_output=True)
        else: self.execute("screen -S " + name + " -L -d -m " + remote_script_path, output=False, timeout=None)

        # Remove the remote shell script
        if not keep_remote_script: self.execute("rm " + remote_script_path, output=False)

        # Navigate back to the previous working directory after changing for the screen output
        if screen_output_path is not None: self.change_cwd(cwd)

    # -----------------------------------------------------------------

    def kill_screen(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Quit the specified screen session
        self.execute("screen -S " + name + " -X quit", output=False)

    # -----------------------------------------------------------------

    def kill_job(self, id):

        """
        This function ...
        :param id:
        :return:
        """

        # Stop the job with the specified ID
        self.execute("qdel " + str(id), output=False)

    # -----------------------------------------------------------------

    def screen_state(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Execute the 'screen -ls' command
        output = self.execute("screen -ls")

        # Loop over the different active screen sessions
        for entry in output:

            # Check if the specified screen name corresponds to the current entry
            if name in entry:

                # Check the state of this screen session
                if "Detached" in entry: return "detached"
                elif "Attached" in entry: return "attached"
                else: raise ValueError("Screen " + name + " in unkown state")

        # If the screen name was not encountered, return None (the screen session has finished or has been aborted)
        return None

    # -----------------------------------------------------------------

    def is_active_screen(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        state = self.screen_state(name)
        return state == "detached" or state == "attached"

    # -----------------------------------------------------------------

    def is_attached_screen(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.screen_state(name) == "attached"

    # -----------------------------------------------------------------

    def is_detached_screen(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.screen_state(name) == "detached"

    # -----------------------------------------------------------------

    def get_requirements(self, processors):

        """
        This function calculates the required amount of nodes and processors per node, given a certain number of
        processors.
        :param processors:
        :return:
        """

        # Calculate the necessary amount of nodes
        nodes = processors // self.cores_per_node + (processors % self.cores_per_node > 0)

        # Determine the number of processors per node
        ppn = processors if nodes == 1 else self.cores_per_node

        # Return the number of nodes and processors per node
        return nodes, ppn

    # -----------------------------------------------------------------

    def launch_pts_command(self, command, arguments, show_output=True):

        """
        This function ...
        :param command:
        :param arguments:
        :param show_output:
        :return:
        """

        # Convert the list of commands back to a single string
        portions = []
        for argument in arguments:
            if " " in argument: portions.append("'" + argument + "'")
            else: portions.append(argument)
        argument_string = " ".join(portions)

        # Execute the command
        self.execute("pts " + command + " " + argument_string, show_output=show_output)

    # -----------------------------------------------------------------

    def run_pts(self, command, config, input_dict=None, keep_remote_output=False, remove_local_output=False):

        """
        This function ...
        :param command:
        :param config:
        :param input_dict:
        :param keep_remote_output:
        :param remove_local_output:
        :return:
        """

        # Create a name for this PTS session
        unique_session_name = time.unique_name(command.replace("/", "-"))

        # Create a remote temporary directory
        remote_temp_path = self.new_temp_directory()

        ##

        ## CHANGE THE LOG PATH TO A REMOTE PATH AND CHANGE THE CWD

        # Get the local working directory path
        local_cwd = config.path

        # Always create a log file while executing remotely
        config.report = True
        config.log_path = remote_temp_path
        config.path = remote_temp_path

        # Don't save the config remotely again (we will upload it explicitly)
        config.config_path = None

        config.write_config = False

        ##

        # DETERMINE REMOTE OUTPUT PATH

        # If the 'output' parameter is specified in the configuration, get the full absolute path to the corresponding output directory
        if "output" in config:

            # Determine the local output path that is desired by the user
            local_output_path = fs.absolute_or_in(config.output, local_cwd)

            # Create the local output path if necessary
            if not fs.is_directory(local_output_path): fs.create_directory(local_output_path)

            # Set the remote output path
            remote_output_path = fs.join(config.path, "out")

            # Change the out path for running on the remote
            config.output = remote_output_path

        # else, the output directory is the remote working directory
        # THIS IS FOR 'OLDER' CONFIGURABLE DERIVED CLASSES THAT DON'T USE THE 'OUTPUT_PATH' and 'OUTPUT_PATH_FILE' MECHANICS ...
        else:

            # Set local output path and remote output path
            local_output_path = local_cwd
            remote_output_path = config.path

        ###

        # Debugging
        log.debug("Saving the configuration file locally ...")

        # Determine path to the temporarily saved local configuration file
        temp_path = tempfile.gettempdir()
        temp_conf_path = fs.join(temp_path, unique_session_name + ".cfg")

        # Save the configuration file to the temporary directory
        config.save(temp_conf_path)

        # Debugging
        log.debug("Uploading the configuration file to '" + remote_temp_path + "' ...")

        # Upload the config file
        remote_conf_path = fs.join(remote_temp_path, fs.name(temp_conf_path))
        self.upload(temp_conf_path, remote_temp_path)

        # Remove the original config file
        fs.remove_file(temp_conf_path)

        #### UPLOAD THE INPUT

        if input_dict is not None:

            # Create remote input path
            input_path = fs.join(remote_temp_path, "input")
            self.create_directory(input_dict)

            for name in input_dict:

                filename = name + "." + input_dict[name].default_extension

                local_filename = fs.join(temp_path, filename) # local temporary path

                # TODO: save the input item here

                filepath = fs.join(input_path, filename) # remote temporary path

                # TODO: UPLOAD THE INPUT item here

                # TODO: delete the local input item

        ####

        # Debugging
        log.debug("Creating a script for remote execution ...")

        # Determine the path to the remote equivalent of this file
        remote_main_path = fs.join(self.pts_package_path, "do", "__main__.py")

        # Create a bash script
        temp_script_path = fs.join(temp_path, unique_session_name + ".sh")

        # Write the lines to the script file
        with open(temp_script_path, 'w') as script_file:

            #script_file.write("#!/usr/bin/env python\n")
            #script_file.write("# -*- coding: utf8 -*-\n")
            #script_file.write("\n")
            script_file.write("python " + remote_main_path + " --configfile " + remote_conf_path + " " + command + "\n")

        # Execute the script
        # name, local_script_path, script_destination, screen_output_path=None, keep_remote_script=False
        self.start_screen(unique_session_name, temp_script_path, remote_temp_path, screen_output_path=remote_temp_path, keep_remote_script=keep_remote_output)

        # Remove the local script file
        fs.remove_file(temp_script_path)

        # Remove the remote temporary directory
        if keep_remote_output: log.info("Remote output will be placed in '" + remote_output_path + "'")

        # Create a new Task object
        task = Task(command, config.to_string())

        # Generate a new task ID
        task_id = self._new_task_id()

        # Determine the path to the task file
        task_file_path = fs.join(self.local_pts_host_run_dir, str(task_id) + ".task")
        task.path = task_file_path

        # Set properties such as the task ID and name and the screen name
        task.id = task_id
        task.remote_temp_pts_path = remote_temp_path
        task.name = unique_session_name
        task.screen_name = unique_session_name
        task.remote_screen_output_path = remote_temp_path

        # Set local and remote output path
        task.local_output_path = local_output_path
        task.remote_output_path = remote_output_path

        # Other
        task.remove_remote_output = not keep_remote_output
        task.remove_local_output = remove_local_output

        # Save the task
        task.save()

        # Return the task
        return task

    # -----------------------------------------------------------------

    def execute(self, command, output=True, expect_eof=True, contains_extra_eof=False, show_output=False, timeout=None):

        """
        This function ...
        :param command:
        :param output:
        :param expect_eof:
        :param contains_extra_eof:
        :param show_output:
        :param timeout:
        :return:
        """

        # Send the command
        self.ssh.sendline(command)

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        # Otherwise, assure that the logfile is set to 'None'
        if show_output: self.ssh.logfile = sys.stdout
        else: self.ssh.logfile = None

        # Retrieve the output if requested
        eof = self.ssh.prompt(timeout=timeout)

        # If an extra EOF is used before the actual output line (don't ask me why but I encounter this on the HPC UGent infrastructure), do prompt() again
        if contains_extra_eof: eof = self.ssh.prompt()

        # If the command could not be sent, raise an error
        if not eof and expect_eof and not contains_extra_eof: raise RuntimeError("The command could not be sent")

        # Set the log file back to 'None'
        self.ssh.logfile = None

        # Ignore the first and the last line (the first is the command itself, the last is always empty)
        if output:
            # Trial and error to get it right for HPC UGent login nodes; don't know what is happening
            if contains_extra_eof: return self.ssh.before.replace('\x1b[K', '').split("\r\n")[1:-1]
            else: return self.ansi_escape.sub('', self.ssh.before).replace('\x1b[K', '').split("\r\n")[1:-1]

    # -----------------------------------------------------------------

    def execute_lines(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        output = kwargs.pop("output", True)
        show_output = kwargs.pop("show_output", False)
        timeout = kwargs.pop("timeout", None)

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        # Otherwise, assure that the logfile is set to 'None'
        if show_output: self.ssh.logfile = sys.stdout
        else: self.ssh.logfile = None

        # Loop over the lines
        for line in args:

            if isinstance(line, basestring):

                # Send the command
                self.ssh.sendline(line)

            elif isinstance(line, tuple):

                # Expect
                self.ssh.expect(line[0])
                self.ssh.sendline(line[1])

            else: raise ValueError("Lines must be strings or tuples")

        # Synchronize
        eof = self.ssh.prompt(timeout=timeout)

        # If the command could not be sent, raise an error
        if not eof: raise RuntimeError("The commands could not be sent")

        # Set the log file back to 'None'
        self.ssh.logfile = None

        # Return the output
        if output: return self.ansi_escape.sub('', self.ssh.before).replace('\x1b[K', '').split("\r\n")[1:-1]

    # -----------------------------------------------------------------

    def in_python_virtual_environment(self):

        """
        This function ...
        :return:
        """

        # python -c 'import sys; print (sys.real_prefix)'

        lines = []
        lines.append("import sys")
        lines.append("answer = False")
        lines.append("answer = True if hasattr(sys, 'real_prefix') else answer")
        lines.append("answer = True if 'conda' in sys.prefix.lower() else answer")
        lines.append("answer = True if 'canopy' in sys.prefix.lower() else answer")
        lines.append("answer = True if 'epd' in sys.prefix.lower() else answer")
        lines.append("answer = True if 'enthought' in sys.prefix.lower() else answer")
        lines.append("print answer")

        # Determine the command
        command = 'python -c "' + "; ".join(lines) + '"'

        # Run the command
        output = self.execute(command)

        if output[0] == "True": return True
        elif output[0] == "False": return False
        else: raise RuntimeError("Something went wrong: answer from remote = '" + str(output))

    # -----------------------------------------------------------------

    def start_python_session(self, assume_pts=True):

        """
        This function ...
        :param assume_pts: assume PTS is present, import some basic PTS tools
        :return:
        """

        if self.in_python_session:
            log.warning("Already in a python session")
            return

        # Initiate a python session
        self.ssh.sendline("python")
        self.ssh.expect("\r\n>>>")

        # Set flag
        self.in_python_session = True

        if assume_pts:

            # Import standard PTS tools
            self.import_python_package("filesystem", as_name="fs", from_name="pts.core.tools")

            # Set logging level to match that of local PTS
            if log.is_debug():

                self.import_python_package("setup_log", from_name="pts.core.tools.logging")
                self.send_python_line("log = setup_log('DEBUG')")

            else: self.import_python_package("log", from_name="pts.core.tools.logging")

    # -----------------------------------------------------------------

    def send_python_line(self, line, output=False, in_loop=False, show_output=False, timeout=30):

        """
        This function ...
        :param line:
        :param output:
        :param in_loop:
        :param show_output:
        :param timeout: default is 30 seconds, use None to have no timeout
        :return:
        """

        if not self.in_python_session:
            log.warning("Not in a remote python session")
            return

        # Show output
        if show_output: self.ssh.logfile = sys.stdout
        else: self.ssh.logfile = None

        # Send line and expect
        self.ssh.sendline(line)
        if in_loop: self.ssh.expect("\r\n...", timeout=timeout)
        else: self.ssh.expect("\r\n>>>", timeout=timeout)

        # Set the log file back to 'None'
        self.ssh.logfile = None

        # Check and return output
        output_lines = self.ssh.before.split("\r\n")[1:]

        # Check for errors
        if len(output_lines) > 0 and "Traceback (most recent call last)" in output_lines[0]:

            # Print error message and traceback
            log.error("Something went wrong during remote execution on host " + self.host_id)
            log.error("Command: " + line)
            for line in output_lines: log.error(line)
            exit()

        if output: return output_lines

    # -----------------------------------------------------------------

    def do_python_loop(self, top_statement, in_loop_lines, show_output=False):

        """
        This function ...
        :param top_statement:
        :param in_loop_lines:
        :param show_output:
        :return:
        """

        # Top statement
        self.send_python_line(top_statement, in_loop=True)

        # In-loop lines
        for line in in_loop_lines: self.send_python_line("    " + line, in_loop=True, show_output=show_output)

        # Finish the loop
        self.send_python_line("")

    # -----------------------------------------------------------------

    def remove_python_variable(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        self.send_python_line("del " + name)

    # -----------------------------------------------------------------

    def import_python_package(self, name, as_name=None, from_name=None):

        """
        This function ...
        :param name:
        :param as_name:
        :param from_name:
        :return:
        """

        command = ""

        if from_name is not None:
            command += "from " + from_name + " "

        command += "import " + name

        if as_name is not None:
            command += " as " + as_name

        # Execute the import command
        #output = self.send_python_line(command, output=True, show_output=log.is_debug())
        output = self.send_python_line(command, output=True)

        # If output is given, this is normally not so good
        if len(output) > 0:

            # Check output
            last_line = output[-1]
            if "cannot import" in last_line: log.warning(last_line)

    # -----------------------------------------------------------------

    def define_simple_python_variable(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        command = name + " = " + str(value)
        self.send_python_line(command)

    # -----------------------------------------------------------------

    def get_python_attributes(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_simple_python_variable("vars(" + name + ")")

    # -----------------------------------------------------------------

    def get_simple_python_variable(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        output = self.send_python_line(name, output=True)

        #if len(output) < 1: raise NameError("No such variable: '" + name + "'")
        if len(output) == 0: return None
        elif len(output) > 1: raise RuntimeError("Unexpected output: " + str(output))

        return eval(output[0])

    # -----------------------------------------------------------------

    def get_simple_python_property(self, variable, name):

        """
        This function ...
        :param variable:
        :param name:
        :return:
        """

        try: return self.get_simple_python_variable(variable + "." + name)
        except NameError: raise AttributeError("Variable '" + variable + "' has no attribute '" + name + "'")

    # -----------------------------------------------------------------

    def get_python_string(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        output = self.send_python_line(name, output=True)
        return output[0][1:-1]

    # -----------------------------------------------------------------

    def python_variables(self):

        """
        This function ...
        :return:
        """

        # Get the complete dir() list of variables
        variables = self.get_simple_python_variable("dir()")

        # '__builtins__', '__doc__', '__name__', '__package__'
        variables.remove("__builtins__")
        variables.remove("__doc__")
        variables.remove("__name__")
        variables.remove("__package__")

        # Return the list of user-defined variables
        return variables

    # -----------------------------------------------------------------

    def end_python_session(self):

        """
        This function ...
        :return:
        """

        if not self.in_python_session:
            log.warning("Not in a remote python session")
            return

        self.ssh.sendline("exit()")
        self.ssh.prompt()

        # Set flag
        self.in_python_session = False

    # -----------------------------------------------------------------

    def execute_python_interactive(self, lines, show_output=False):

        """
        This function ...
        :param lines:
        :param show_output:
        :return:
        """

        # If already in python session
        if self.in_python_session:

            # Send the lines consecutively
            for line in lines: self.send_python_line(line, show_output=show_output)

        # Not yet in python session
        else:

            # Don't import PTS stuff
            self.start_python_session(assume_pts=False)

            # Send the lines consecutively
            for line in lines: self.send_python_line(line, show_output=show_output)

            # End the python session
            self.end_python_session()

    # -----------------------------------------------------------------

    def execute_python_script(self, script_path, show_output=False):

        """
        This function ...
        :param script_path:
        :param show_output:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def rename_file(self, directory, old_name, new_name):

        """
        This function ...
        :param directory:
        :param old_name:
        :param new_name:
        :return:
        """

        if self.in_python_session: self.send_python_line("fs.rename_file('" + directory + "', '" + old_name + "', '" + new_name + "')")
        else:

            # Determine the old and new file path
            old_path = fs.join(directory, old_name)
            new_path = fs.join(directory, new_name)

            # Use the 'mv' command to rename the file
            self.execute("mv " + old_path + " " + new_path)

    # -----------------------------------------------------------------

    def remove_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        if self.in_python_session: self.send_python_line("fs.remove_directory('" + path + "')")
        else: self.execute("rm -rf " + path, output=False)

    # -----------------------------------------------------------------

    def remove_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        if self.in_python_session: self.send_python_line("fs.remove_file('" + path + "')")
        else: self.execute("rm " + path, output=False)

    # -----------------------------------------------------------------

    def change_cwd(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        if self.in_python_session: raise RuntimeError("Cannot change the working directory while in a remote python interactive session.")
        else:

            # Try to change the directory, give an error if this fails
            output = self.execute("cd " + path)
            if len(output) > 0 and "No such file or directory" in output[0]: raise RuntimeError("The directory does not exist")

    # -----------------------------------------------------------------

    def directories_in_path(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        if self.in_python_session: return self.get_simple_python_variable("fs.directories_in_path('" + path + "')")
        else:

            # Get the path to the current working directory
            working_directory = self.working_directory

            # Change the working directory to the provided path
            self.change_cwd(path)

            # List the directories in the provided path
            output = self.execute("for i in $(ls -d */); do echo ${i%%/}; done")

            # Change the working directory back to the original working directory
            self.change_cwd(working_directory)

            # Return the list of directories
            return output

    # -----------------------------------------------------------------

    def files_in_path(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        if self.in_python_session: return self.get_simple_python_variable("fs.files_in_path('" + path + "')")
        else:

            # Get the path to the current working directory
            working_directory = self.working_directory

            # Change the working directory to the provided path
            self.change_cwd(path)

            # List the files in the provided path
            output = self.execute("for f in *; do [[ -d $f ]] || echo $f; done")

            # Change the working directory back to the original working directory
            self.change_cwd(working_directory)

            # Return the list of directories
            return output

    # -----------------------------------------------------------------

    def to_home_directory(self):

        """
        This function ...
        """

        if self.in_python_session: raise RuntimeError("Cannot navigate to other directories whilst in a remote python interactive session")

        # Navigate to the home directory
        self.execute("cd ~", output=False)

    # -----------------------------------------------------------------

    def create_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        if self.in_python_session: self.send_python_line("fs.create_directory('" + path + "')")
        else:

            # Create the remote directory
            self.execute("mkdir " + path, output=False)

    # -----------------------------------------------------------------

    def create_directories(self, *paths):

        """
        This function ...
        :return:
        """

        if self.in_python_session: self.send_python_line("fs.create_directories(*" + str(paths) + ")")
        else:

            # Create the remote directories
            self.execute("mkdir " + " ".join(paths), output=False)

    # -----------------------------------------------------------------

    def append_line(self, filepath, line):

        """
        This function ...
        :param filepath:
        :param line:
        :return:
        """

        command = "echo '" + line + "' >>" + filepath
        self.execute(command, output=False)

    # -----------------------------------------------------------------

    def append_lines(self, filepath, lines):

        """
        This function ...
        :param filepath:
        :param lines:
        :return:
        """

        for line in lines: self.append_line(filepath, line)

    # -----------------------------------------------------------------

    def download(self, origin, destination, timeout=60, new_name=None, compress=False, show_output=False):

        """
        This function ...
        :param origin:
        :param destination:
        :param timeout:
        :param new_name:
        :param compress:
        :param show_output:
        :return:
        """

        # If debugging is enabled, always show the scp output
        if log.is_debug(): show_output = True

        # Construct the command string
        copy_command = "scp "
        if compress: copy_command += "-C "

        # Add the host address
        copy_command += self.host.user + "@" + self.host.name + ":"

        # If the origin is a string, we assume it represents a single file path or directory path
        if isinstance(origin, basestring):

            # Check if the origin represents a file
            #if self.is_file(origin): copy_command += origin.replace(" ", "\\\ ") + " "
            if self.is_file(origin): copy_command += "'" + origin + "' "

            # Check if it represents a directory
            #elif self.is_directory(origin): copy_command += origin.replace(" ", "\\ ") + "/* " + "-r "
            #elif self.is_directory(origin): copy_command += origin.replace(" ", "\\\ ") + "/* "
            elif self.is_directory(origin): copy_command += "'" + origin + "/*' "

            # The origin does not exist
            else: raise ValueError("The specified path " + origin + " does not represent an existing directory or file on the remote host")

        # If the origin is a list, we assume it contains multiple file paths
        elif isinstance(origin, list):

            # Check whether the files exist remotely
            for file_path in origin:
                if not self.is_file(file_path): raise ValueError("The file " + file_path + " does not exist on the remote host")

            # Escape possible space characters
            #origin = [path.replace(" ", "\\\ ") for path in origin]
            origin = ["'" + path + "'" for path in origin]

            # Add a quotation mark character because the seperate file paths are going to be separated by spaces
            # (the command is going to be of the form scp username@ip.of.server.copyfrom:"file1.log file2.log" "~/yourpathtocopy")
            copy_command += '"'

            # Add the file paths to the command string
            copy_command += " ".join(origin)

            # Add another quotation mark to identify the end of the filepath list
            copy_command += '" '

        # Add the destination path to the command
        #copy_command += destination.replace(" ", "\\\ ") + "/"
        copy_command += "'" + destination + "/'"
        if new_name is not None: copy_command += new_name

        # Debugging
        log.debug("Copy command: " + copy_command)

        # Create the pexpect child instance
        child = pexpect.spawn(copy_command, timeout=timeout)
        if self.host.password is not None:
            child.expect(['password: '])
            child.sendline(self.host.password)

        # If the output does not have to be shown on the console, create a temporary file where the output is written to
        if not show_output:

            # Temporary file for output of the scp command
            #temp_file_path = tempfile.mktemp()
            #temp_file = open(temp_file_path, 'w')

            # New way: use a string stream
            temp_file = StringIO.StringIO()
            child.logfile = temp_file

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        else: child.logfile = sys.stdout

        # Execute the command and get the output
        child.expect(pexpect.EOF, timeout=None)
        child.close()

        if not show_output:

            # Retrieve file contents -- this will be
            # 'First line.\nSecond line.\n'
            stdout = child.logfile.getvalue()

            # Raise an error if something went wrong
            if child.exitstatus != 0: raise RuntimeError(stdout)

            # Get the output lines
            lines = stdout.split("\n")

            # Check for messages that signal an error that has occured
            for line in lines:
                if "not a regular file" in line: raise ValueError(line)

            # Debugging: show the output of the scp command
            log.debug("Copy stdout: " + str(" ".join(lines)))

    # -----------------------------------------------------------------

    def upload(self, origin, destination, timeout=60, new_name=None, compress=False, show_output=False):

        """
        This function ...
        :param origin:
        :param destination:
        :param timeout:
        :param new_name:
        :param compress:
        :param show_output:
        :return:
        """

        # If debugging is enabled, always show the scp output
        if log.is_debug(): show_output = True

        # Construct the command string
        copy_command = "scp "
        if compress: copy_command += "-C "

        # If the origin is a string, we assume it represents a single file path or directory path
        if isinstance(origin, basestring):

            # Check if the origin represents a file
            #if fs.is_file(origin): copy_command += origin.replace(" ", "\\\ ") + " "
            if fs.is_file(origin): copy_command += "'" + origin + "' "

            # Check if it represents a directory
            #elif fs.is_directory(origin): copy_command += "-r " + origin.replace(" ", "\\\ ") + "/ "
            elif fs.is_directory(origin): copy_command += "-r '" + origin + "/' "

            # The origin does not exist
            else: raise ValueError("The specified path " + origin + " does not represent an existing directory or file")

        # If the origin is a list, we assume it contains multiple file paths
        elif isinstance(origin, list):

            # Check whether the files exist locally
            for file_path in origin:
                if not fs.is_file(file_path): raise ValueError("The file " + file_path + " does not exist")

            origin = ["'" + path + "'" for path in origin]

            # Add the file paths to the command string
            copy_command += " ".join(origin) + " "

        # Invalid argument
        else: raise ValueError("The origin must be a string or a list of strings")

        # Add the host address and the destination directory
        #copy_command += self.host.user + "@" + self.host.name + ":" + destination.replace(" ", "\\\ ") + "/"
        copy_command += self.host.user + "@" + self.host.name + ":'" + destination + "/'"
        if new_name is not None: copy_command += new_name

        # Debugging
        log.debug("Copy command: " + copy_command)

        # Create the pexpect child instance
        child = pexpect.spawn(copy_command, timeout=timeout)
        if self.host.password is not None:
            child.expect(['password: '])
            child.sendline(self.host.password)

        # If the output does not have to be shown on the console, create a temporary file where the output is written to
        if not show_output:

            # Temporary file for output of the scp command
            #temp_file_path = tempfile.mktemp()
            #temp_file = open(temp_file_path, 'w')

            # New way: use a string stream
            temp_file = StringIO.StringIO()
            child.logfile = temp_file

        # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
        else: child.logfile = sys.stdout

        # Execute the command and get the output
        try:
            child.expect(pexpect.EOF, timeout=None)
        except pexpect.EOF:
            pass
        child.close()

        if not show_output:

            # Retrieve file contents -- this will be
            # 'First line.\nSecond line.\n'
            stdout = child.logfile.getvalue()

            # Raise an error if something went wrong
            if child.exitstatus != 0: raise RuntimeError(stdout)

            # Get the output lines
            lines = stdout.split("\n")

            # Check for messages that signal an error that has occured
            for line in lines:
                if "not a regular file" in line: raise ValueError(line)

            # Debugging: show the output of the scp command
            log.debug("Copy stdout: " + str(" ".join(lines)))

    # -----------------------------------------------------------------

    def read_text_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Expand the path to absolute form
        path = self.expand_user_path(path)

        # Load the text file into a variable
        self.execute("value='cat " + path + "'")

        # Print the variable to the console, and obtain the output
        return self.execute('echo "$($value)"')

    # -----------------------------------------------------------------

    def install(self, private=False, key_password=None):

        """
        This function ...
        :param private:
        :param key_password:
        :return:
        """

        # Navigate to the home directory
        self.execute("cd ~", output=False)

        # Create the SKIRT directory
        self.create_directory("SKIRT")

        # In the SKIRT directory, create the necessary subdirectories
        self.execute("cd SKIRT", output=False)
        self.create_directories("git", "run", "release")

        # Clone the SKIRT repository
        if private:
            output = self.execute("git clone git@github.ugent.be:SKIRT/SKIRT.git git", expect_eof=False)
            self.ssh.expect(['id_rsa: '])
            self.ssh.sendline(key_password)

        else: self.execute("git clone https://github.com/SKIRT/SKIRT.git git", output=False)

        # Compile the SKIRT code
        self.execute("./makeSKIRT.sh", output=False)

        # Put SKIRT in the PATH environment variable

    # -----------------------------------------------------------------

    @property
    def system_name(self):

        """
        This function ...
        :return:
        """

        return self.host.system_name

    # -----------------------------------------------------------------

    @property
    def home_directory(self):

        """
        This function ...
        :return:
        """

        # If we are in a python session
        if self.in_python_session:

            self.import_python_package("expanduser", from_name="os.path")
            return self.get_python_string("expanduser('~')")

        # If we are not in a python session
        else:

            # Find out the path to the user's home directory and return it
            output = self.execute("echo $HOME")
            return output[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def session_temp_directory(self):

        """
        This function ...
        :return:
        """

        # Determine the path to a new temporary directory
        path = fs.join(self.pts_temp_path, time.unique_name("session"))

        # Create the directory from within the remote python session
        if self.in_python_session:

            # Create from python
            self.send_python_line("fs.create_directory('" + path + "')")

        # Create the directory using a bash command
        else: self.create_directory(path)

        # Return the path to the new temporary directory
        return path

    # -----------------------------------------------------------------

    def new_temp_directory(self):

        """
        This function ...
        :return:
        """

        # Generate the path to a new unique temporary directory
        path = fs.join(self.pts_temp_path, time.unique_name("new"))

        # Create the directory from within the remote python session
        if self.in_python_session:

            # Create from python
            self.send_python_line("fs.create_directory('" + path + "')")

        # Create the directory using a bash command
        else: self.create_directory(path)

        # Return the path to the new temporary directory
        return path

    # -----------------------------------------------------------------

    def expand_user_path(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        if not path.startswith("~"): return path
        else: return fs.join(self.home_directory, path.split("~/")[1])

    # -----------------------------------------------------------------

    def find_executable(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the output of the 'which' command
        output = self.execute("which " + name)

        if len(output) == 0: return None

        # Only one line is expected
        return output[0]

    # -----------------------------------------------------------------

    def is_executable(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.find_executable(name) is not None

    # -----------------------------------------------------------------

    @property
    def working_directory(self):

        """
        This function ...
        :return:
        """

        # Find out the path to the current working directory and return it
        output = self.execute("echo $PWD")
        return output[0]

    # -----------------------------------------------------------------

    @property
    def free_cores(self):

        """
        This function ...
        :return:
        """

        # The calculation doesn't work for multiple nodes
        if self.nodes > 1: raise ValueError("This function cannot be called for systems with multiple nodes")

        # Determine the number of free cores on the only node
        return self.cores_per_node * (1.0 - self.cpu_load)

    # -----------------------------------------------------------------

    @property
    def free_sockets(self):

        """
        This function ...
        :return:
        """

        # The calculation doesn't work for multiple nodes
        if self.nodes > 1: raise ValueError("This function cannot be called for systems with multiple nodes")

        # Determine the number of free cores on the only node
        return self.sockets_per_node * (1.0 - self.cpu_load)

    # -----------------------------------------------------------------

    @property
    def free_memory(self):

        """
        This function ...
        :return:
        """

        # Use the 'free' command to get information about the virtual memory usage
        #output = self.execute("free -t | grep 'Total'")
        output = self.execute("free -t | grep buffers/cache")
        splitted = output[0].split(":")[1].split()

        # Calculate the free amount of memory in gigabytes
        free = float(splitted[1]) / 1e6

        # Return the free amount of virtual memory in gigabytes
        return free

    # -----------------------------------------------------------------

    @property
    def free_space(self):

        """
        This function ...
        :return:
        """

        # Use the 'df' command to obtain information about the free disk space
        output = self.execute("df -lh")

        total = 0.0
        used = 0.0
        free = 0.0

        for entry in output[1:]:

            if not entry.startswith("/dev/"): continue

            splitted = entry.split()

            total += float(splitted[1].split("G")[0]) if "G" in splitted[1] else float(splitted[1].split("T")[0]) * 1e3
            used += float(splitted[2].split("G")[0]) if "G" in splitted[2] else float(splitted[2].split("T")[0]) * 1e3
            free += float(splitted[3].split("G")[0]) if "G" in splitted[3] else float(splitted[3].split("T")[0]) * 1e3

        # Return the amount of free memory in gigabytes
        return free

    # -----------------------------------------------------------------

    @property
    def scheduler(self):

        """
        This property ...
        :return:
        """

        return self.host.scheduler

    # -----------------------------------------------------------------

    @property
    def use_hyperthreading(self):

        """
        This function ...
        :return:
        """

        return self.host.use_hyperthreading

    # -----------------------------------------------------------------

    @property
    def multi_node_communication(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, check whether multi node communication is possible based on
        # the configuration of the current cluster
        if self.scheduler: return self.host.cluster.multi_node_communication

        # If no scheduler is used, raise an error (this function should not get called)
        else: raise RuntimeError("This function should only be called when using a remote with a scheduling system")

    # -----------------------------------------------------------------

    @property
    def tty(self):

        """
        This function ...
        :return:
        """

        output = self.execute("tty")
        return int(output.split("pts/")[1])

    # -----------------------------------------------------------------

    @property
    def ttys(self):

        """
        This function ...
        :return:
        """

        output = self.execute("who | awk '{print $2,$NF}' |grep -v '(:[0-9]'")

        session_numbers = []
        for line in output: session_numbers.append(int(line.split("pts/")[1].split(" ")[0]))
        return session_numbers

    # -----------------------------------------------------------------

    @property
    def virtual_memory_per_node(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the amount of virtual memory memory per node is defined
        # in the configuration (this is in GB)
        if self.scheduler: return self.host.cluster.memory

        # If no scheduler is used, assume the number of nodes is 1 and get the total virtual memory (total swap)
        else:

            output = self.execute("free -t | grep Swap")
            splitted = output[0].split(":")[1].split()

            # Calculate the free amount of memory in gigabytes
            total_swap = float(splitted[0]) / 1e6

            # Return the free amount of virtual memory in gigabytes
            return total_swap

    # -----------------------------------------------------------------

    @property
    def nodes(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the number of nodes is defined in the host configuration
        if self.scheduler: return self.host.cluster.nodes

        # If no scheduling system is used, assume the system is only concised of one node
        else: return 1

    # -----------------------------------------------------------------

    @property
    def cores_per_node(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the number of cores on the computing nodes is defined in the configuration
        if self.scheduler: return self.host.cluster.cores_per_socket * self.host.cluster.sockets_per_node

        # If no scheduler is used, the computing node is the actual node we are logged in to
        else:

            # Use the 'lscpu' command to obtain the total number of CPU's (=hardware threads!)
            output = self.execute("lscpu | grep '^CPU(s)'")
            cpus = int(float(output[0].split(":")[1]))

            # Return the number of physical cores
            return cpus / self.threads_per_core

    # -----------------------------------------------------------------

    @property
    def threads_per_core(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the number of threads per core is defined in the configuration
        if self.scheduler: return self.host.cluster.threads_per_core

        # If no scheduler is used, the computing node is the actual node we are logged in to
        else:

            # Use the 'lscpu' command to get the number of hardware threads per core
            output = self.execute("lscpu | grep '^Thread(s) per core'")
            threads_per_core = int(float(output[0].split(":")[1]))

            # Return the amount of hyperthreads or 'hardware' threads per physical core
            return threads_per_core

    # -----------------------------------------------------------------

    @property
    def sockets_per_node(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the number of sockets per node is defined in the configuration
        if self.scheduler: return self.host.cluster.sockets_per_node

        # If no scheduler is used, the computing node is the actual node we are logged in to
        else:

            # Use the 'lscpu' command to get the number of NUMA domains
            output = self.execute("lscpu | grep '^Socket(s):'")
            nsockets = int(output[0])

            # Return the number of sockets
            return nsockets

    # -----------------------------------------------------------------

    @property
    def cores_per_socket(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the number of cores per socket is defined in the configuration
        if self.scheduler: return self.host.cluster.cores_per_socket

        # If no scheduler is used, the computing node is the actual node we are logged in to
        else:

            # Use the 'lscpu' command
            output = self.execute("lscpu | grep '^Core(s) per socket:'")
            ncores = int(output[0])

            # Return the number of cores per socket
            return ncores

    # -----------------------------------------------------------------

    @property
    def numa_domains(self):

        """
        This function returns the number of NUMA domains (per node)
        :return:
        """

        # If the remote host uses a scheduling system, the number of numa domains is defined in the configuration
        if self.scheduler: return self.host.cluster.numa_domains_per_node

        # If no scheduler is used, the computing node is the actual node we are logged in to
        else:

            # Use the 'lscpu' command to get the number of NUMA domains
            output = self.execute("lscpu | grep '^NUMA node(s):'")
            numa_domains = int(output[0])

            # Return the number of NUMA domains
            return numa_domains

    # -----------------------------------------------------------------

    @property
    def numa_cpus(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the CPU ranks
        cpus = [[] for _ in range(self.numa_domains)]

        # Only works for up to 10 domains!
        output = self.execute("lscpu")

        for line in output:

            # Skip irrelevant lines
            if not line.startswith("NUMA node") or line.startswith("NUMA node(s)"): continue

            # Get the NUMA domain
            domain = int(line.split("NUMA node")[1].split(" CPU")[0])

            # Get the list of CPU's
            string = line.split(": ")[1].strip()
            cpu_list = parsing.integer_list(string)

            # Set the CPU list for the current NUMA domain
            cpus[domain] = cpu_list

        # Return the list of CPU lists for each NUMA domain
        return cpus

    # -----------------------------------------------------------------

    def cpus_for_numa_domain(self, domain_index):

        """
        This function ...
        :param domain_index:
        :return:
        """

        return self.numa_cpus[domain_index]

    # -----------------------------------------------------------------

    @property
    def hpc_ugent_node_status(self):

        """
        This function ...
        :return:
        """

        lines = ["from vsc.jobs.pbs.nodes import collect_nodeinfo"]
        lines.append("node_list, state_list, types = collect_nodeinfo()")
        lines.append("print node_list")
        lines.append("print state_list")
        lines.append("print types")

        # For interpreting 'types':
        #template = "%sppn=%s, physmem=%sGB, swap=%sGB, vmem=%sGB, local disk=%sGB"
        #for typ, nodes in sorted(types.items(), key=lambda x: len(x[1]), reverse=True):
            # most frequent first
            #cores, phys, swap, disk = typ
            #txt.append(template % (offset, cores, phys, swap, phys + swap, disk))

        output = self.execute_python_interactive(lines)

        return output

    # -----------------------------------------------------------------

    @property
    def cpu_load(self):

        """
        This function ...
        :return:
        """

        # Use the 'top' command to get the current CPU load
        output = self.execute("top -b -n1 | grep 'Cpu(s)' | awk '{print $2 + $4}'")

        # Convert the output into the fraction of CPU that is used
        load = float(output[0]) / 100.0

        # Return the current CPU load
        return load

    # -----------------------------------------------------------------

    @property
    def memory_load(self):

        """
        This function ...
        :return:
        """

        # Use the 'free' command to get information about the virtual memory usage
        output = self.execute("free -t | grep 'Total'")
        splitted = output[0].split(":")[1].split()

        # Calculate the total and used amount of memory in gigabytes
        total = float(splitted[0]) / 1e6
        used = float(splitted[1]) / 1e6
        free = float(splitted[2]) / 1e6

        # Return the fraction of virtual memory that is currently in use
        return used / total

    # -----------------------------------------------------------------

    @property
    def platform(self):

        """
        This fucntion ...
        :return:
        """

        lines = []
        lines.append("import platform")
        lines.append("answer = None")
        lines.append("answer = 'MacOS' if platform.system() == 'Darwin' else answer")
        lines.append("answer = 'Windows' if platform.system() == 'Windows' else answer")
        lines.append("answer = 'Linux' if platform.system() == 'Linux' else answer")
        lines.append("print answer")

        # Determine the command
        command = 'python -c "' + "; ".join(lines) + '"'

        # Run the command
        output = self.execute(command)

        return output[0]

    # -----------------------------------------------------------------

    @property
    def architecture(self):

        """
        This function ...
        :return:
        """

        # architecture
        lines = []
        lines.append("import platform")
        lines.append("answer = None")
        lines.append("answer = '64bit' if platform.architecture()[0] == '64bit' else answer")
        lines.append("print answer")

        # Determine the command
        command = 'python -c "' + "; ".join(lines) + '"'

        # Run the command
        output = self.execute(command)

        return output[0]

    # -----------------------------------------------------------------

    @property
    def python_distributions(self):

        """
        This function ...
        :return:
        """

        # Dictionary
        packages = dict()

        # pip list
        output = self.execute("pip list")

        # Loop over the lines in the output
        for entry in output:

            # Known messages that corrupt the output
            if "You are using pip version" in entry: continue
            if "You should consider upgrading via the" in entry: continue

            # Get name and version
            name, version = entry.split(" (")
            version = version[:-1]

            # Add the package name with its version
            packages[name] = version

        # Return the dictionary of python packages with their version numbers
        return packages

    # -----------------------------------------------------------------

    @property
    def installed_python_packages(self):

        """
        This function ...
        :return:
        """

        # Check if in python session
        was_in_python_session = self.in_python_session

        # Start python session
        if not was_in_python_session: self.start_python_session()

        # Import PTS introspection tools
        self.import_python_package("introspection", from_name="pts.core.tools")

        # Get packages
        packages = self.get_simple_python_variable("introspection.installed_python_packages()")

        # End python session again
        if not was_in_python_session: self.end_python_session()

        # Return the dictionary of python packages with their version numbers
        return packages

    # -----------------------------------------------------------------

    def is_present_package(self, package):

        """"
        This function ...
        :param package:
        """

        # Check if in python session
        was_in_python_session = self.in_python_session

        # Start python session
        if not was_in_python_session: self.start_python_session()

        # Import PTS introspection tools
        self.import_python_package("introspection", from_name="pts.core.tools")

        # Check if present
        present = self.get_simple_python_variable("introspection.is_present_package('" + package + "')")

        # End python session again
        if not was_in_python_session: self.end_python_session()

        # Return
        return present

    # -----------------------------------------------------------------

    def file_or_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        if self.in_python_session: return self.get_simple_python_variable("fs.file_or_directory('" + path + "')")

        else:

            # Launch a bash command to check whether the path exists as either a file or a directory on the remote file system
            output = self.execute("if [ -f " + path + " ]; then echo file; elif [ -d " + path + " ]; then echo directory; else echo None; fi")

            # Return the result
            if output[0] == "None": return None
            else: return output[0]

    # -----------------------------------------------------------------

    def evaluate_boolean_expression(self, expression):

        """
        This function ...
        :param expression:
        :return:
        """

        # The commmand
        command = "if [ " + expression + " ]; then echo True; else echo False; fi"

        # Launch a bash command to check whether the path exists as a directory on the remote file system
        output = self.execute(command)

        # Return the result
        return output[0] == "True"

    # -----------------------------------------------------------------

    def is_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Call the corresponding PTS function
        if self.in_python_session: return self.get_simple_python_variable("fs.is_directory('" + path + "')")

        # Launch a bash command to check whether the path exists as a directory on the remote file system
        else: return self.evaluate_boolean_expression("-d " + path)

    # -----------------------------------------------------------------

    def is_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Call the corresponding PTS function
        if self.in_python_session: return self.get_simple_python_variable("fs.is_file('" + path + "')")

        # Launch a bash command to check whether the path exists as a regular file on the remote file system
        else: return self.evaluate_boolean_expression("-f " + path)

    # -----------------------------------------------------------------

    def is_subdirectory(self, path, parent_path):

        """
        This function ...
        :param path:
        :param parent_path:
        :return:
        """

        if self.in_python_session: return self.get_simple_python_variable("fs.is_subdirectory('" + path + "', '" + parent_path + "')")
        else:

            if not self.is_directory(path): raise ValueError("Not a directory: " + path)

            path = self.expand_user_path(path)
            parent_path = self.expand_user_path(parent_path)
            return path.startswith(parent_path)

    # -----------------------------------------------------------------

    @property
    def host_id(self):

        """
        This function ...
        :return:
        """

        return self.host.id

    # -----------------------------------------------------------------

    @property
    def cluster_name(self):

        """
        This function ...
        :return:
        """

        return self.host.cluster_name

    # -----------------------------------------------------------------

    @property
    def skirt_root_path(self):

        """
        This function ...
        :return:
        """

        path = self.expand_user_path("~/SKIRT")
        return path

    # -----------------------------------------------------------------

    @property
    def skirt_repo_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.skirt_root_path, "git")

    # -----------------------------------------------------------------

    @property
    def pts_root_path(self):

        """
        This function ...
        :return:
        """

        path = self.expand_user_path("~/PTS")
        return path

    # -----------------------------------------------------------------

    @property
    def pts_package_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.pts_root_path, "pts")

    # -----------------------------------------------------------------

    @property
    def pts_temp_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.pts_root_path, "temp")
        if not self.is_directory(path): self.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @lazyproperty
    def local_pts_host_run_dir(self):

        """
        This function ...
        :return:
        """

        path = fs.join(introspection.pts_run_dir, self.host.id)
        if not fs.is_directory(path): fs.create_directory(path, recursive=True)
        return path

    # -----------------------------------------------------------------

    def retrieve_tasks(self):

        """
        This function ...
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Initialize a list to contain the tasks that have been retrieved
        tasks = []

        # Loop over the different entries of the status list
        for path, task_status in self.get_task_status():

            # If a task has been retrieved earlier, but is not yet analysed, also add it to the list of retrieved
            # tasks (again) so that its results can be analysed
            if task_status == "retrieved":

                # Open the task file
                task = Task.from_file(path)

                # Add the task to the list
                tasks.append(task)

            # Finished task
            elif task_status == "finished":

                # Open the task file
                task = Task.from_file(path)

                # CHECK WHETHER OUTPUT HAS TO BE DOWNLOADED
                if task.local_output_path is not None:

                    # Debug info
                    log.debug("Retrieving the complete remote output directory ...")
                    log.debug("Local output directory: " + task.local_output_path)
                    log.debug("Remote output directory: " + task.remote_output_path)

                    # Check whether the output directory exists; if not, create it
                    if not fs.is_directory(task.local_output_path): fs.create_directory(task.local_output_path)

                    # Download the PTS task output
                    self.download(task.remote_output_path, task.local_output_path)

                # Local output path not defined
                else: log.warning("Local output path not defined: remote PTS task output will not be retrieved (look on the remote filesystem for the results in '" + task.remote_output_path + "')")

                # Add the retrieved task to the list
                tasks.append(task)

                # If retrieval was succesful, add this information to the task file
                task.retrieved = True
                task.save()

                ## REMOVE REMOTE OUTPUT IF REQUESTED
                if task.remove_remote_output:

                    # Remove the temporary PTS directory if it contains the output directory
                    if self.is_subdirectory(task.remote_output_path, task.remote_temp_pts_path): self.remove_directory(task.remote_temp_pts_path)
                    else:
                        # Remove the output directory and the temporary directory seperately
                        self.remove_directory(task.remote_output_path)
                        self.remove_directory(task.remote_temp_pts_path)

        # Return the list of retrieved tasks
        return tasks

    # -----------------------------------------------------------------

    def get_task_status(self):

        """
        This function ..
        :return:
        """

        # Initialize a list to contain the statuses
        entries = []

        # If the remote host does not use a scheduling system
        if not self.scheduler:

            # Search for task files in the local PTS run/host_id directory
            for path in fs.files_in_path(self.local_pts_host_run_dir, extension="task", sort=int):

                # Open the task file
                task = Task.from_file(path)

                # Check whether the task has already been analysed
                if task.analysed: task_status = "analysed"

                # Check whether the task has already been retrieved
                elif task.retrieved: task_status = "retrieved"

                # Not yet retrieved, could be any other state ...
                else:

                    # Get screen name and output path
                    screen_name = task.screen_name
                    output_path = task.remote_output_path
                    log_output_path = task.remote_log_path

                    # Check whether the remote output path exists
                    if not self.is_directory(output_path): task_status = "invalid: remote output directory has been deleted"

                    # The remote output path exists
                    else:

                        # Check whether the log file is present
                        log_path = None
                        for filename in self.files_in_path(log_output_path):
                            if "log" in filename:
                                log_path = fs.join(log_output_path, filename)
                                break

                        # Check whether the report file exists
                        if log_path is not None:

                            if self.is_active_screen(screen_name): task_status = "running"
                            else:

                                # Get the last two lines of the remote log file
                                output = self.execute("tail -2 " + log_path)
                                last_line = output[1]

                                # Check whether the last line states that the task has finished
                                if "Finished" in last_line: task_status = "finished"
                                elif "Error:" in last_line: task_status = "crashed: " + last_line
                                else: task_status = "aborted"

                        # If the log file does not exist, the task has not started yet or has been cancelled
                        else:

                            # The task has not started or it's screen session has been cancelled
                            if self.is_active_screen(screen_name): task_status = "queued"
                            else: task_status = "cancelled"

                # Add the task properties to the list
                entries.append((path, task_status))

        # If the remote has a scheduling system for launching jobs
        #else: raise NotImplementedError("Running PTS tasks on remote schedulers is not supported yet")

        # Return the list of task properties
        return entries

    # -----------------------------------------------------------------

    def _task_ids_in_use(self):

        """
        This function ...
        :return:
        """

        # Check the contents of the local run directory to see which task id's are currently in use
        current_ids = []
        for name in fs.files_in_path(self.local_pts_host_run_dir, extension="task", returns="name"):

            # Get the task ID and add it to the list
            current_ids.append(int(name))

        # Return the list of currently used ID's
        return current_ids

    # -----------------------------------------------------------------

    def _new_task_id(self):

        """
        This function ...
        :return:
        """

        # Get a list of the ID's currently in use
        current_ids = self._task_ids_in_use()

        # Sort the current task ID's and find the lowest 'missing' integer number
        if len(current_ids) > 0:
            current_ids = sorted(current_ids)
            task_id = max(current_ids) + 1
            for index in range(max(current_ids)):
                if current_ids[index] != index:
                    task_id = index
                    break

            # Return the task ID
            return task_id

        # If no task ID's are currently in use, return 0
        else: return 0

# -----------------------------------------------------------------
