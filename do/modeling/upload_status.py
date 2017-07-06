#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.upload_status Upload the status page.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import webbrowser

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.modeling.html.generator import HTMLGenerator
from pts.core.remote.host import Host
from pts.core.tools import introspection
from pts.core.remote.mounter import RemoteMounter

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition()
definition.add_flag("generate", "first (re)generate the HTML", False)
definition.add_flag("check", "check by opening in the browser", False)
config = parse_arguments("upload_status", definition)

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

# Determine the filename and the url
filename = environment.galaxy_name + ".html"
url = "http://users.ugent.be/~sjversto/" + filename

# -----------------------------------------------------------------

# Generate the HTML
if config.generate:

    # Generate the HTML
    generator = HTMLGenerator()
    generator.config.path = modeling_path
    generator.run()

# -----------------------------------------------------------------

# Get account
username, password = introspection.get_account("ugent.be")

# Create host
host = Host("www", name="files.ugent.be", user=username, password=password, mount_point=username + "/www/users", protocol="smb")

# Mount
mounter = RemoteMounter()
mount_path = mounter.mount(host)

# -----------------------------------------------------------------

# Get the status page
fs.copy_file(environment.html_status_path, mount_path, new_name=filename)

# -----------------------------------------------------------------

# Unmount
mounter.unmount(host)

# -----------------------------------------------------------------

# CHECK
if config.check: webbrowser.open(url, new=2)

# -----------------------------------------------------------------
