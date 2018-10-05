#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.backup_website Backup your personal website locally.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.remote.mounter import RemoteMounter
from pts.core.tools import archive

# -----------------------------------------------------------------

# Mount
mounter = RemoteMounter()
mount_path = mounter.mount("www")

# -----------------------------------------------------------------

new_name = "www_backup"
backup_path = fs.join(fs.home, new_name)

# -----------------------------------------------------------------

# Copy contents
fs.copy_directory(mount_path, fs.home, new_name=new_name)

# -----------------------------------------------------------------

# Unmount
mounter.unmount("www")

# -----------------------------------------------------------------

# Create zip
archive.compress_directory_in_place(backup_path, remove=True)

# -----------------------------------------------------------------
