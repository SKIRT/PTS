#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.zip Transparently opening files from a ZIP archive.
#
# The functions in this module allow opening a file contained in a ZIP archive without first extracting the file.

# -----------------------------------------------------------------

import os
import os.path
import zipfile
import StringIO

# -----------------------------------------------------------------

## This function opens a text file in read-only mode. If a file exists at the specified path, it is simply opened.
# Otherwise, the function looks for a ZIP archive with the same name as the directory in which the file would have
# resided, but with the ".zip" extension added, and it attempts to open a file with the same name from the archive.
# In both cases, the function returns a read-only file-like object that offers sequential access, i.e. it provides
# only the following methods: read(), readline(), readlines(), \_\_iter\_\_(), next().
#
def opentext(filepath):
    # if the specified file exists, simply open it
    if os.path.isfile(filepath):
        return open(filepath, 'r')

    # otherwise, try a zip archive with the same name as the directory in which the file would have resided
    directory,filename = os.path.split(filepath)
    zippath = directory + ".zip"
    return zipfile.ZipFile(zippath,'r').open(filename,'r')

## This function opens a binary file in read-only mode. If a file exists at the specified path, it is simply opened.
# Otherwise, the function looks for a ZIP archive with the same name as the directory in which the file would have
# resided, but with the ".zip" extension added, and it attempts to open a file with the same name from the archive.
# In both cases, the function returns a read-only file-like object that offers full random access functionality.
#
# In case the file is opened from a ZIP archive, the complete file contents is loaded into a memory buffer. This is
# necessary to enable random access to the decompressed data stream.
#
def openbinary(filepath):
    # if the specified file exists, simply open it
    if os.path.isfile(filepath):
        return open(filepath, 'rb')

    # otherwise, try a zip archive with the same name as the directory in which the file would have resided
    directory,filename = os.path.split(filepath)
    zippath = directory + ".zip"
    return StringIO.StringIO(zipfile.ZipFile(zippath,'r').read(filename))

## This function returns a list of the names of the files in the specified directory. If the directory does not exist,
# the function looks for a ZIP archive with the same name as the directory, but with the ".zip" extension added,
# and it lists the names of the files in that archive instead.
#
def listdir(dirpath):
    # if the specified directory exists, list its contents
    if os.path.isdir(dirpath):
        return os.listdir(dirpath)

    # otherwise, try a zip archive with the same name as the directory
    zippath = dirpath + ".zip"
    return zipfile.ZipFile(zippath,'r').namelist()

# -----------------------------------------------------------------

