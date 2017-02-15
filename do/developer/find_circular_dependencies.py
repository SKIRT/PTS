#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.find_circular_dependencies Find circular dependencies.

## THIS SCRIPT DOES NOT WORK YET!

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import introspection

# -----------------------------------------------------------------

def find_links_to(module_path, modules, dependencies_dict, level=0, already_processed=None):

    """
    This function ...
    :param level:
    """

    #links = []

    #print(module_path, modules)

    if already_processed is None: already_processed = []

    #print(already_processed)

    for path in modules:

        if path in already_processed: continue

        #if path == module_path: links.append((level, path))
        #else: links += find_links_to(module_path, dependencies_dict[path], dependencies_dict, level=level+1, already_processed=already_processed)

        links = find_links_to(module_path, dependencies_dict[path], dependencies_dict, level=level + 1,
                               already_processed=already_processed)
        if links: return module_path, path

        already_processed.append(path)

    #return links
    return False

# -----------------------------------------------------------------

modules = introspection.get_internal_dependencies()

# Loop over the modules
for path in modules:

    # Loop over the dependencies
    #for dep_path in modules[path]:

    links = find_links_to(path, modules[path], modules)

    if links: print(links)

    #for link in links:

    #    print(link)

# -----------------------------------------------------------------