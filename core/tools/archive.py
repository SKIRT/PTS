#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.archive Compressing, decompressing and transparently opening files from ZIP/GZ/BZ2 archives.
#
# The functions in this module allow opening a file contained in a ZIP archive without first extracting the file.
# It also allows compressing and decompressing .gz and .bz2 files.

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import zipfile
import StringIO
import shutil
import gzip
import bz2
import subprocess

# Import the relevant PTS classes and modules
from . import filesystem as fs

# -----------------------------------------------------------------

extensions = ["bz2", "zip", "gz"]

# -----------------------------------------------------------------

# Signatures for different compressed file formats
signatures = dict()
signatures["\x1f\x8b\x08"] = "gz"
signatures["\x42\x5a\x68"] = "bz2"
signatures["\x50\x4b\x03\x04"] = "zip"
signatures["\x52\x61\x72\x21\x1a\x07\x00"] = "rar4"    # 52 61 72 21 1A 07 00
signatures["\x52\x61\x72\x21\x1a\x07\x01\x00"] = "rar5" # 52 61 72 21 1A 07 01 00

# -----------------------------------------------------------------

max_signature_length = max(len(x) for x in signatures)

# -----------------------------------------------------------------

def is_compressed(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    return compression_type(filepath) is not None

# -----------------------------------------------------------------

def compression_type(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    #
    file_start = fs.read_start(filepath, max_signature_length)

    for magic, filetype in signatures.items():
        if file_start.startswith(magic):
            return filetype

    return None

# -----------------------------------------------------------------

def fix_extension(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    ct = compression_type(filepath)
    if ct is None: raise ValueError("Not a compressed file")

    filename = fs.name(filepath)

    # Rename if necessary
    if not filename.endswith(ct):
        new_filepath = filepath + "." + ct
        fs.rename_file_path(filepath, new_filepath)
    else: new_filepath = filepath

    # Return the filepath
    return filepath

# -----------------------------------------------------------------

def is_archive(filepath):

    """
    This function ...
    :param filepath: 
    :return: 
    """

    for extension in extensions:
        if filepath.endswith("." + extension): return True
    return False

# -----------------------------------------------------------------

def bare_name(filepath):

    """
    This function ...
    :param filepath: 
    :return: 
    """

    for extension in extensions:
        if filepath.endswith("." + extension): return filepath.split("." + extension)[0]
    return filepath

# -----------------------------------------------------------------

def decompress_directory_in_place(filepath, remove=False, into_root=False):

    """
    This function ...
    :param filepath:
    :param remove:
    :param into_root:
    :return:
    """

    from ..basics.log import log

    # Inform the user
    log.info("Decompressing '" + filepath + "' ...")

    # Tar.gz
    if filepath.endswith(".tar.gz"):

        # Determine the path of the directory
        new_path = filepath.split(".tar.gz")[0]
        dir_path = fs.directory_of(new_path)

        # Debugging
        log.debug("New path: '" + new_path + "'")
        log.debug("Decompressing in directory '" + dir_path + "' ...")

        # Decompress
        command = "tar -zxvf " + filepath + " --directory " + dir_path
        log.debug("Decompress command: '" + command + "'")
        subprocess.call(command, shell=True)

    else: raise NotImplementedError("Not implemented yet")

    if into_root:

        for path in fs.files_in_path(new_path):
            fs.move_file(path, dir_path)
        fs.remove_directory(new_path)
        new_path = dir_path

    # Remove the file
    if remove: fs.remove_file(filepath)

    # Return the new path
    return new_path

# -----------------------------------------------------------------

def compress_directory_in_place(directory_path, remove=False):

    """
    This function ...
    :param directory_path:
    :param remove:
    :return:
    """

    from ..basics.log import log

    # Determine filename
    name = fs.name(directory_path)
    filename = name + ".tar.gz"
    filepath = fs.join(fs.directory_of(directory_path), filename)

    # Inform the user
    log.info("Compressing '" + directory_path + "' ...")
    command = "tar -zcvf " + filepath + " " + directory_path
    log.debug("Decompress command: '" + command + "'")
    subprocess.call(command, shell=True)

    # Remove the file
    if remove: fs.remove_directory(directory_path)

    # Return the new path
    return filepath

# -----------------------------------------------------------------

def decompress_file_in_place(path, remove=False):

    """
    This function ...
    :param path:
    :param remove:
    :return:
    """

    from ..basics.log import log

    # Inform the user
    log.info("Decompressing '" + path + "' ...")

    # Check extension
    if path.endswith(".bz2"):
        new_path = path.rstrip(".bz2")
        decompress_bz2(path, new_path)
    elif path.endswith(".gz"):
        new_path = path.rstrip(".gz")
        if new_path.endswith(".tar"): new_path = new_path.split(".tar")[0]
        decompress_gz(path, new_path)
    elif path.endswith(".zip"):
        new_path = path.rstrip(".zip")
        decompress_zip(path, new_path)
    else: raise ValueError("Unrecognized archive type (must be bz2, gz [or tar.gz] or zip)")

    # Remove the original file if requested
    if remove: fs.remove_file(path)

    # Return the new path
    return new_path

# -----------------------------------------------------------------

def decompress_file(path, new_path):

    """
    This funtion ...
    :param path:
    :param new_path:
    :return:
    """

    if path.endswith(".bz2"): decompress_bz2(path, new_path)
    elif path.endswith(".gz"): decompress_gz(path, new_path)
    elif path.endswith(".xz"): decompress_xz(path, new_path)
    elif path.endswith(".zip"): decompress_zip(path, new_path)
    else: raise ValueError("Unrecognized archive type (must be bz2, gz or zip)")

# -----------------------------------------------------------------

def decompress_files(filepaths, remove=False):

    """
    This function ...
    :param filepaths:
    :param remove:
    :return:
    """

    # Initialize a list for the decompressed file paths
    new_paths = []

    # Loop over the files
    for filepath in filepaths:

        # Get the name of the file
        filename = fs.name(filepath)

        # Get directory of the file
        path = fs.directory_of(filepath)

        # Strip the bz2 extension
        newfilename = fs.strip_extension(filename)

        # Determine path to new file
        newfilepath = fs.join(path, newfilename)

        # Decompress this file
        decompress_file(filepath, newfilepath)

        # If succesful, add the new path to the list
        new_paths.append(newfilepath)

    # Remove original files if requested
    if remove: fs.remove_files(filepaths)

    # Return the list of new file paths
    return new_paths

# -----------------------------------------------------------------

def decompress_zip(zip_path, new_path):

    """
    This function decompresses a .zip file
    :param zip_path:
    :param new_path:
    :return:
    """

    # If directory is specified
    if fs.is_directory(new_path):
        name = fs.name(zip_path).rstrip(".zip")
        new_path = fs.join(new_path, name)

    with zipfile.ZipFile(zip_path, 'w') as myzip:
        myzip.write(new_path)

# -----------------------------------------------------------------

def decompress_gz(gz_path, new_path):

    """
    This function decompresses a .gz file
    :param gz_path:
    :param new_path:
    :return:
    """

    # If directory is specified
    if fs.is_directory(new_path):
        name = fs.name(gz_path).rstrip(".gz")
        if name.endswith(".tar"): name = name.split(".tar")[0]
        new_path = fs.join(new_path, name)

    # Decompress
    with gzip.open(gz_path, 'rb') as f_in:
        with open(new_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

# -----------------------------------------------------------------

def decompress_xz(xz_path, new_path):

    """
    This function ...
    :param xz_path:
    :param new_path:
    :return:
    """

    # If directory is specified
    if fs.is_directory(new_path):
        name = fs.name(xz_path).rstrip(".xz")
        if name.endswith(".tar"): name = name.split(".tar")[0]
        new_path = fs.join(new_path, name)

    # Decompress
    with gzip.open(xz_path, 'rb') as f_in:
        with open(new_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

# -----------------------------------------------------------------

def decompress_bz2(bz2_path, new_path):

    """
    This function decompresses a .bz2 file
    :param bz2_path:
    :param new_path:
    :return:
    """

    # If directory is specified
    if fs.is_directory(new_path):
        name = fs.name(bz2_path).rstrip(".bz2")
        new_path = fs.join(new_path, name)

    # Decompress, create decompressed new file
    with open(new_path, 'wb') as new_file, bz2.BZ2File(bz2_path, 'rb') as file:
        for data in iter(lambda: file.read(100 * 1024), b''):
            new_file.write(data)

# -----------------------------------------------------------------

## Get lines
def get_lines(filepath):

    # File exists
    if fs.is_file(filepath): return fs.get_lines(filepath)

    # Look for compressed file
    directory, filename = os.path.split(filepath)
    zippath = directory + ".zip"
    if not fs.is_file(zippath): raise IOError("File is not found")

    # Return the lines from the compressed file
    with zipfile.ZipFile(zippath, 'r').open(filename, 'r') as fh: return fs.get_lines_filehandle(fh)

# -----------------------------------------------------------------

## This function opens a text file in read-only mode. If a file exists at the specified path, it is simply opened.
# Otherwise, the function looks for a ZIP archive with the same name as the directory in which the file would have
# resided, but with the ".zip" extension added, and it attempts to open a file with the same name from the archive.
# In both cases, the function returns a read-only file-like object that offers sequential access, i.e. it provides
# only the following methods: read(), readline(), readlines(), \_\_iter\_\_(), next().
#
def opentext(filepath):

    # If the specified file exists, simply open it
    if fs.is_file(filepath): return open(filepath, 'r')

    # otherwise, try a zip archive with the same name as the directory in which the file would have resided
    directory,filename = os.path.split(filepath)
    zippath = directory + ".zip"
    if not fs.is_file(zippath): raise IOError("File is not found")

    # Return the filehandle
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

## This function returns True if the specified file exists at the specified path and/or inside a ZIP archive with
# the same name as the directory in which the file would have resided, but with the ".zip" extension added.
# Otherwise the function returns False.
#
def isfile(filepath):
    # if the specified file exists, we are done
    if os.path.isfile(filepath):
        return True

    # otherwise, try a zip archive with the same name as the directory in which the file would have resided
    directory,filename = os.path.split(filepath)
    zippath = directory + ".zip"
    if os.path.isfile(zippath):
        try:
            zipfile.ZipFile(zippath,'r').getinfo(filename)
            return True
        except KeyError:
            pass
    return False

## This function returns a sorted list of the names of the files in the specified directory and/or in the ZIP archive
# with the same name as the directory, but with the ".zip" extension added. If both the directory and the archive
# exist, the two lists are merged while removing duplicates. The returned list is optionally limited to
# filenames that end in the string (or strings) specified as the second argument to this function.
#
def listdir(dirpath, endswith=None):
    filenames = [ ]
    # if the specified directory exists, list its contents
    if os.path.isdir(dirpath):
        filenames += os.listdir(dirpath)
    # if the corresponding zip archive exists, list its contents as well
    zippath = dirpath + ".zip"
    if os.path.isfile(zippath):
        filenames += zipfile.ZipFile(zippath,'r').namelist()
    # if requested, filter the names
    if endswith != None:
        filenames = filter(lambda fn: fn.endswith(endswith), filenames)
    # remove duplicates and sort the list
    return sorted(list(set(filenames)))

# -----------------------------------------------------------------

