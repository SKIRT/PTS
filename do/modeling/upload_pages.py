#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.upload_pages Upload the HTML pages.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.modeling.html.all import AllPagesGenerator
from pts.core.remote.host import load_host, Host
from pts.core.remote.mounter import RemoteMounter
from pts.modeling.core.environment import verify_modeling_cwd
from pts.core.tools import browser

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition()
definition.add_flag("generate", "first (re)generate the HTML", False)
definition.add_flag("check", "check by opening in the browser", False)
definition.add_optional("credentials", "username_password", "username and password")
config = parse_arguments("upload_status", definition)

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

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
    generator = AllPagesGenerator()
    generator.config.path = modeling_path
    generator.run()

# -----------------------------------------------------------------

# Credentials are given: create temporary host runtime object
if config.credentials is not None:

    # Create
    host = Host("www")

    # Set properties
    host.mount_point = "sjversto/www/users"
    host.protocol = "smb"
    host.name = "files.ugent.be"

    # Set username and password
    host.user = config.credentials.username
    host.password = config.credentials.password

# Load the host
else: host = load_host("www")

# -----------------------------------------------------------------

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
if config.check: browser.open_url(url)

# -----------------------------------------------------------------
