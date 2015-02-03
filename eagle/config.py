#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package eagle.config Configuration settings and paths for the eagle package.
#
# This module defines configuration settings and paths for the eagle package as the values of
# variables in its outermost namespace. For example, after executing "import eagle.config" the name
# of the host on which the script is running can be found as the value of "eagle.config.hostname".
#
# The following table lists the variables defined by this module.
#
#<TABLE>
#<TR><TD><B>Variable name</B></TD>  <TD><B>Description of value</B></TD></TR>
#<TR><TD>hostname</TD>              <TD>The name of the host on which the Python script is running</TD></TR>
#<TR><TD>skirt_path</TD>            <TD>The absolute path to the SKIRT executable</TD></TR>
#<TR><TD>eagledata_path</TD>        <TD>A dictionary containing key-value pairs providing the absolute path to the
#                                       eagle data directory (value) for each relevant eagle simulation (key)</TD></TR>
#<TR><TD>catalogs_path</TD>         <TD>The absolute path to the directory containing the catalog files corresponding
#                                       to each of the the eagle snapshots used with this package</TD></TR>
#<TR><TD>database_path</TD>         <TD>The absolute path to the directory containing the SKIRT-runs database</TD></TR>
#<TR><TD>backup_path</TD>           <TD>The absolute path to the directory containing backups of the database</TD></TR>
#<TR><TD>templates_path</TD>        <TD>The absolute path to the directory containing ski file templates</TD></TR>
#<TR><TD>results_path</TD>          <TD>The absolute path to the directory containing the SKIRT results</TD></TR>
#<TR><TD>queue</TD>                 <TD>The name of the queue to which jobs should be submitted, or None</TD></TR>
#<TR><TD>default_eaglesim</TD>      <TD>The identifier of the eagle simulation currently in use (must be one of the
#                                       keys in the \em eagledata_path dictionary)</TD></TR>
#<TR><TD>default_redshift</TD>      <TD>The redshift of the snapshot currently in use</TD></TR>
#</TABLE>
#
# In addition this module offers some utility functions for simple tasks such as obtaining a time stamp.
#

# -----------------------------------------------------------------

import datetime
import os.path
import socket

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

# get the name of the user logged in on the terminal controlling this process
username = os.getlogin()

# get the name of the current host so we can let the configuration settings depend on it
if socket.gethostname().find('.')>=0:
    hostname=socket.gethostname()
else:
    hostname=socket.gethostbyaddr(socket.gethostname())[0]

# -----------------------------------------------------------------

# configuration for the COSMA cluster in Durham
if "cosma" in hostname:
    skirt_path = absolutepath("~/SKIRT/release/SKIRTmain/skirt")
    eagledata_path = { 'Ref100':  "/cosma5/data/Eagle/ScienceRuns/Planck1/L0100N1504/PE/REFERENCE/data/",
                       'Ref25':   "/cosma5/data/Eagle/ScienceRuns/Planck1/L0025N0752/PE/REFERENCE/data",
                       'Recal25': "/cosma5/data/Eagle/ScienceRuns/Planck1/L0025N0752/PE/RECALIBRATED/data",
                       'Ref12':   "/cosma5/data/Eagle/ScienceRuns/Planck1/L0012N0188/PE/REFERENCE/data" }
    catalogs_path =  "/cosma5/data/Eagle/SkirtAnalysis/Catalogs"
    database_path =  "/cosma5/data/Eagle/SkirtAnalysis/Database"
    backup_path =    "/cosma5/data/Eagle/SkirtAnalysis/Backup"
    templates_path = "/cosma5/data/Eagle/SkirtAnalysis/Templates"
    results_path =   "/cosma5/data/Eagle/SkirtAnalysis/Results"
    queue = "cosma5"
    default_eaglesim = 'Ref12'
    default_redshift = 0

# -----------------------------------------------------------------

# configuration for Peter's desktop at work
elif "obiwan" in hostname:
    skirt_path = absolutepath("~/SKIRT/release/SKIRTmain/skirt")
    eagledata_path = { 'Ref100': absolutepath("~/EAGLE/Snapshots/L0100N1504REF"),
                       'Ref25': absolutepath("~/EAGLE/Snapshots/L0025N0752REF"),
                       'Recal25': absolutepath("~/EAGLE/Snapshots/L0025N0752RECAL"),
                       'Ref12': absolutepath("~/EAGLE/Snapshots/L0012N0188REF") }
    catalogs_path =  absolutepath("~/EAGLE/Catalogs")
    database_path =  absolutepath("~/EAGLE/Database")
    backup_path =    absolutepath("~/EAGLE/Backup")
    templates_path = absolutepath("~/EAGLE/Templates")
    results_path =   absolutepath("~/EAGLE/Results")
    queue = None
    default_eaglesim = 'Ref12'
    default_redshift = 0

# -----------------------------------------------------------------

else:
    raise ValueError("Unknown host: " + hostname)
