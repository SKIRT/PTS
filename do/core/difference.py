#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.difference Report differences, if any, between the 'out' and 'ref' directories of a test case.
#
# This script reports the differences, if any, between the contents of the 'out' and 'ref' directories of a test case.
# The report is more extensive than the single line in the test report produced by the skirttestsuite class.
#
# If there are no arguments, the script compares the 'out' and 'ref' subdirectories in the current directory.
# Otherwise the first argument is interpreted as the relative or absolute path of a directory, and the script
# compares the 'out' and 'ref' subdirectories of that directory.

# -----------------------------------------------------------------

# Import standard modules
import sys
import os.path
import filecmp

# Import the relevant PTS classes and modules
from pts.core.test.skirt import equalfiles

# -----------------------------------------------------------------

def report_difference(casedirpath):
    # get the directories to be compared
    refpath = os.path.join(casedirpath, "ref")
    outpath = os.path.join(casedirpath, "out")
    if not os.path.isdir(refpath):
        print "Test case has no reference directory"
        return
    if not os.path.isdir(refpath):
        print "Test case has no output directory"
        return

    # check for recursive subdirectories
    if len(filter(lambda fn: os.path.isdir(fn), os.listdir(refpath))) > 0:
        print "Reference directory contains a subdirectory"
        return
    if len(filter(lambda fn: os.path.isdir(fn), os.listdir(outpath))) > 0:
        print "Output directory contains a sub directory"
        return

    # verify list of filenames
    dircomp = filecmp.dircmp(outpath, refpath, ignore=['.DS_Store'])
    if (len(dircomp.left_only) > 0):
        print "Output contains " + str(len(dircomp.left_only)) + " extra file(s)"
    if (len(dircomp.right_only) > 0):
        print "Output misses " + str(len(dircomp.right_only)) + " file(s)"

    # compare common files
    matches, mismatches, errors = filecmp.cmpfiles(outpath, refpath, dircomp.common, shallow=False)
    for filename in matches:
        print "Output file matches: " + filename
    for filename in mismatches + errors:
        if equalfiles(os.path.join(outpath, filename), os.path.join(refpath, filename)):
            print "Output file matches: " + filename
        else:
            print "Output file differs: " + filename + "     <-------"

# -----------------------------------------------------------------

# get the command-line argument specifying the target directory, if any
argument = sys.argv[1] if len(sys.argv) > 1 else ""
casedirpath = os.path.realpath(os.path.expanduser(argument))

# perform the comparison
print "Reporting on test case results in " + casedirpath + "..."
report_difference(casedirpath)
print "End of report"

# -----------------------------------------------------------------
