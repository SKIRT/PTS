#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
# Main module for the do package
# -----------------------------------------------------------------

## \package pts.do.__main__ Execute one of the scripts in this package from any current directory.
#
# Before proceeding, ensure that your login script includes the extra lines as described in the Installation Guide
# and that you have logged in again after that change.
#
# To execute a python script named \c example.py residing in this directory, enter "pts example" in a Terminal window.
# This will work regardless of your current directory.
# You can include command line arguments as well, as in "pts example arg1 arg2".
# You can shorten the script name to the first few letters as long as there are no two scripts with matching names.
# For example "pts exa arg1 arg2" would still execute \c example.py, assuming that only one script has a name
# starting with the string "exa".
#
# Use "ipts" rather than "pts" to enter interactive python mode (with the >>> prompt) after executing the script.
# This is useful for testing or experimenting with pts functionality: the script imports the relevant
# pts module(s) and initializes some key objects that then can be used from the interactive python prompt.
#

# -----------------------------------------------------------------

# Import standard modules
import sys
import argparse

# Import the relevant PTS modules
from pts.core.tools import introspection, parsing
from pts.core.tools import filesystem as fs
from pts.do.run import run_configurable, no_match, no_input, show_version, multiple_matches

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser(prog="pts")
parser.add_argument("do_command", type=str, help="the name of the PTS do command (preceeded by the subproject name and a slash if ambigious; i.e. 'subproject/do_command')", default=None, nargs='?')
parser.add_argument("--version", "-v", action="store_true", help="show the PTS version")
parser.add_argument("--interactive", action="store_true", help="use interactive mode for the configuration")
parser.add_argument("--arguments", action="store_true", help="use argument mode for the configuration")
parser.add_argument("--configfile", type=str, help="use a configuration file")
parser.add_argument("--rerun", action="store_true", help="use the last used configuration")
parser.add_argument("--remote", type=str, help="launch the PTS command remotely")
parser.add_argument("--keep", action="store_true", help="keep the remote output")
parser.add_argument("--input", type=str, help="the name/path of the input directory")
parser.add_argument("--output", type=str, help="the name/path of the output directory")
parser.add_argument("--input_files", type=parsing.string_tuple_dictionary, help="dictionary of (class_path, input_file_path) where the key is the input variable name")
parser.add_argument("--output_files", type=parsing.string_string_dictionary, help="dictionary of output file paths where the key is the variable (attribute) name")
parser.add_argument("options", nargs=argparse.REMAINDER, help="options for the specific do command")

# -----------------------------------------------------------------

# Find possible PTS commands
scripts = introspection.get_scripts()
tables = introspection.get_arguments_tables()

# -----------------------------------------------------------------

# No input?
if len(sys.argv) == 1: no_input(parser, scripts, tables)

# -----------------------------------------------------------------

# Parse the command-line arguments
args = parser.parse_args()

# -----------------------------------------------------------------

if args.version: show_version()

# -----------------------------------------------------------------

# Check input and output options, should be directories
if args.input is not None and not fs.is_directory(args.input): raise ValueError("Input path should be an existing directory")
if args.output is not None and not fs.is_directory(args.output): raise ValueError("Output path should be an existing directory")

# -----------------------------------------------------------------

# Get the name of the do script
script_name = args.do_command

# Construct clean arguments list
sys.argv = ["pts", args.do_command] + args.options

# -----------------------------------------------------------------

# Find matches
matches = introspection.find_matches_scripts(script_name, scripts)
table_matches = introspection.find_matches_tables(script_name, tables)

# -----------------------------------------------------------------

# No match
if len(matches) + len(table_matches) == 0: no_match(script_name, scripts, tables)

# If there is a unique match in an existing script, return it
elif len(matches) == 1 and len(table_matches) == 0: #run_script(matches, args) # TEMPORARILY DOESN'T WORK BECAUSE NOT ALL DO SCRIPTS ARE FUTURE PROOF

    if args.remote is not None: raise ValueError("This do command cannot be executed remotely")

    match = matches[0]

    # Execute the matching script, after adjusting the command line arguments so that it appears that the script was executed directly
    target = fs.join(introspection.pts_do_dir, match[0], match[1])
    sys.argv[0] = target
    del sys.argv[1]
    print("Executing: " + match[0] + "/" + match[1] + " " + " ".join(sys.argv[1:]))

    # Execute the script
    exec open(target)

# If there is an unique match in a table
elif len(table_matches) == 1 and len(matches) == 0: run_configurable(table_matches, args, tables)

# Show possible matches if there are more than just one
else: multiple_matches(matches, table_matches, tables)

# -----------------------------------------------------------------
