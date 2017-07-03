#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.config Configuration settings and paths for the eagle package.
#
# This module defines configuration settings and paths for the eagle package on the host computing system.
# For example, after executing "import eagle.config" the path to the SKIRT executable for
# the system on which the script is running can be found as the value of "eagle.config.skirt_path".
# The following table lists the configuration variables defined in this module.
#
#<TABLE>
#<TR><TD><B>Variable name</B></TD>  <TD><B>Description of value</B></TD></TR>
#<TR><TD>skirt_path</TD>            <TD>The absolute path to the SKIRT executable</TD></TR>
#<TR><TD>eagledata_path</TD>        <TD>A dictionary containing key-value pairs providing the absolute path to the
#                                       eagle data directory (value) for each relevant eagle simulation (key)</TD></TR>
#<TR><TD>database_path</TD>         <TD>The absolute path to the directory containing the SKIRT-runs database</TD></TR>
#<TR><TD>backup_path</TD>           <TD>The absolute path to the directory containing backups of the database</TD></TR>
#<TR><TD>templates_path</TD>        <TD>The absolute path to the directory containing ski file templates</TD></TR>
#<TR><TD>jobs_path</TD>             <TD>The absolute path to the directory containing batch job scripts and logs</TD></TR>
#<TR><TD>results_path</TD>          <TD>The absolute path to the directory containing the SKIRT results</TD></TR>
#<TR><TD>collections_path</TD>      <TD>The absolute path to the directory containing the collections with
#                                       statistics on sets of SKIRT simulation results</TD></TR>
#<TR><TD>plots_path</TD>            <TD>The absolute path to the directory containing the plots showing
#                                       statistics on sets of SKIRT simulation results</TD></TR>
#</TABLE>
#
# The module also offers some utility functions for simple tasks such as obtaining a time stamp.
#

# -----------------------------------------------------------------

import datetime
import os.path

# -----------------------------------------------------------------

## This function returns a string representing the current time and date in the format "YYYY-MM-DD--hh-mm-ss".
# This format ensures proper collation (strings sort in date/time order), is easy to read for a human user,
# and can be used as part of a filename on any platform (since there are no nasty characters).
def timestamp():
    return datetime.datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

## This function returns an absolute version of the specified path. The original path may be absolute, relative to
# a user's home folder, or relative to the current working directory.
def absolutepath(path):
    return os.path.realpath(os.path.expanduser(path))

# -----------------------------------------------------------------

# provide a dictionary of configuration values for the Cosma6 cluster in Durham
configuration = {
    'skirt_path': absolutepath("~/SKIRT8/release/SKIRT/main/skirt"),
    'eagledata_path': {
        # Lowres simulations
        'RefL0100N1504':    "/cosma5/data/Eagle/ScienceRuns/Planck1/L0100N1504/PE/REFERENCE/data/",
        'RefL0050N0752':    "/cosma5/data/Eagle/ScienceRuns/Planck1/L0050N0752/PE/REFERENCE/data",
        'AGNdT9L0050N0752': "/cosma5/data/Eagle/ScienceRuns/Planck1/L0050N0752/PE/S15_AGNdT9/data",
        'RefL0025N0376':    "/cosma5/data/Eagle/ScienceRuns/Planck1/L0025N0376/PE/REFERENCE/data",
        # Hires simulations:
        'RefL0025N0752':    "/cosma5/data/Eagle/ScienceRuns/Planck1/L0025N0752/PE/REFERENCE/data",
        'RecalL0025N0752':  "/cosma5/data/Eagle/ScienceRuns/Planck1/L0025N0752/PE/RECALIBRATED/data",
    },
    'database_path':    "/cosma6/data/dp004/pcamps/Eagle/Database",
    'backup_path':      "/cosma6/data/dp004/pcamps/Eagle/Backup",
    'templates_path':   "/cosma6/data/dp004/pcamps/Eagle/Templates",
    'jobs_path':        "/cosma6/data/dp004/pcamps/Eagle/Jobs",
    'results_path':     "/cosma6/data/dp004/pcamps/Eagle/Results",
    'collections_path': "/cosma6/data/dp004/pcamps/Eagle/Collections",
    'plots_path':       "/cosma6/data/dp004/pcamps/Eagle/Plots",
}

# -----------------------------------------------------------------

# copy the configuration settings to the module's namespace
globals().update(configuration)

# -----------------------------------------------------------------

# get the user's account information for the public EAGLE database, if available
public_eagle_database_username = None
public_eagle_database_password = None
try:
    public_eagle_database_username, public_eagle_database_password = \
        open(os.path.join(database_path,"public_eagle_database_account_info.txt")).readline().split()
except Exception:
    pass

# -----------------------------------------------------------------
