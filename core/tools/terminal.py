#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.terminal Provides functions for interacting with the terminal (launching commands and getting output).

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import sys
import pexpect
#import subprocess

# -----------------------------------------------------------------

def execute(command, output=True, show_output=False, timeout=None, expect=None, cwd=None):

    """
    This function ...
    :return:
    """

    # WITH SUBPROCESS:
    #output = subprocess.check_output(command, shell=True, stdout=sys.stdout, stderr=sys.stderr, cwd=cwd)
    #if output: return output.split("\n")[:-1]

    # WITH PEXPECT:
    child = pexpect.spawn(command, timeout=timeout, cwd=cwd)

    # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
    # Otherwise, assure that the logfile is set to 'None'
    if show_output: child.logfile = sys.stdout
    else: child.logfile = None

    # Expect
    if expect is not None: child.expect(expect)
    else: child.expect(pexpect.EOF)

    # Set the log file back to 'None'
    child.logfile = None

    # Ignore the first and the last line (the first is the command itself, the last is always empty)
    if output: return child.before.replace('\x1b[K', '').split("\r\n")[1:-1]

# -----------------------------------------------------------------

def execute_lines(*args, **kwargs):

    """
    This function ...
    :return:
    """

    # Get arguments
    output = kwargs.pop("output", True)
    show_output = kwargs.pop("show_output", False)
    timeout = kwargs.pop("timeout", None)
    cwd = kwargs.pop("cwd", None)

    # Execute first line
    assert isinstance(args[0], basestring)
    child = pexpect.spawn(args[0], timeout=timeout, cwd=cwd)

    # If the output has to be shown on the console, set the 'logfile' to the standard system output stream
    # Otherwise, assure that the logfile is set to 'None'
    if show_output: child.logfile = sys.stdout
    else: child.logfile = None

    # Loop over the lines
    for line in args[1:]:

        # Just a command where completion is expected
        if isinstance(line, basestring):

            # Send the command
            child = child.sendline(line)
            #child.expect()
            #child.expect("$", timeout=timeout)

        # Tuple: something is expected and must be filled in
        elif isinstance(line, tuple):

            # Expect
            if len(line) == 3 and line[2]:

                # index = self.ssh.expect([self.ssh.PROMPT, line[0]]) # this is not working, why?
                index = child.expect(["$", line[0]], timeout=timeout)
                if index == 0: pass
                elif index == 1: child.sendline(line[1])
                # eof = self.ssh.prompt()
            else:

                #self.ssh.expect(line[0])
                #self.ssh.sendline(line[1])
                child.expect(line[0], timeout=timeout)
                child.sendline(line[1])

        # Invalid
        else: raise ValueError("Lines must be strings or tuples")

    # Expect
    child.expect(pexpect.EOF)

    # Set the log file back to 'None'
    child.logfile = None

    # Return the output
    if output: return child.before.split("\r\n")[1:-1]

# -----------------------------------------------------------------
