#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import os
import imp
import inspect
from collections import defaultdict

# -----------------------------------------------------------------

pts_package_dir = os.path.dirname(inspect.getfile(inspect.currentframe()))

# Create an empty dictionary to contain the required modules together with the places of use
modules = defaultdict(list)

# Recursively loop over all files inside this working directory
for directory, subdirs, files in os.walk(pts_package_dir):
    
    # Loop over all files in the (sub)directory
    for filename in files:

        # If the file is not a python script, skip it
        if not filename.endswith(".py"): continue
        
        # Determine the full path to this file
        filepath = os.path.join(directory, filename)
        
        # Read the lines of the script file
        for line in open(filepath, 'r'):
            
            # Look for an 'import yyy' or 'from yyy import zzz' statement
            if line.startswith("import ") or (line.startswith("from ") and "import" in line): 
                
                # Get the name of the module
                module = line.split()[1].split(".")[0]

                # Get the path of the script, relative to the 'PTS/git' directory
                rel_filepath = filepath.split("PTS/pts/")[1]
                    
                # Add the module name to the list
                if module: modules[module].append(rel_filepath)

# -----------------------------------------------------------------

# List all required modules
for module, files in modules.items(): 
    
    # Check whether the module is present on this system
    try: 
        imp.find_module(module)
        message = "present"
    except ImportError:
        found = False
        message = "not found!"
    
    # List the name of the module
    print(module, ":", message)
    
    # List the files where this module is used
    for file in files: print("  - ", file)

# -----------------------------------------------------------------
