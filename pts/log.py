#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.log Support for logging
#
# The class in this module supports logging of time-stamped messages to the console and to a report file.

# -----------------------------------------------------------------

import datetime
import os
import os.path

# -----------------------------------------------------------------
#  Utility functions
# -----------------------------------------------------------------

## This private utility function returns the specified message decorated with a timestamp.
def timestamp(message, delimiter="   "):
    stampedmessage = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S.%f")[:-3] + delimiter + message
    return stampedmessage

# -----------------------------------------------------------------
#  SkirtTestSuite class
# -----------------------------------------------------------------

## An instance of the Log class provides logging services to the console, with an optional copy to a text file.
class Log:

    ## The constructor takes two optional arguments: \em reportpath specifies the path of the directory to
    # the file which will receive a copy of all log messages; \em reportname specifies the name of the
    # log file (without filename extension; a timestamp is automatically added to the filename). If the
    # \em reportpath argument is empty or missing, no report file is created.
    def __init__(self, reportpath="", reportname=""):
        self._doreport = True if reportpath else False
        if self._doreport:
            reportpath = os.path.realpath(os.path.expanduser(reportpath))
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d--%H-%M-%S")
            reportfilepath = os.path.join(reportpath, "report_" + reportname + "_" + timestamp + ".txt")
            self._report = open(reportfilepath, 'w')

    ## This function logs an informational message.
    def info(self, message):
        stampedmessage = timestamp(message)
        if self._doreport:
            self._report.write(stampedmessage + "\n")
            self._report.flush()
            os.fsync(self._report.fileno())
        print stampedmessage

    ## This function logs a warning message.
    def warning(self, message):
        stampedmessage = timestamp(message, delimiter=" ! ")
        if self._doreport:
            self._report.write(stampedmessage + "\n")
            self._report.flush()
            os.fsync(self._report.fileno())
        print "\033[35m" + stampedmessage + "\033[0m"

    ## This function logs a success message.
    def success(self, message):
        stampedmessage = timestamp(message, delimiter=" - ")
        if self._doreport:
            self._report.write(stampedmessage + "\n")
            self._report.flush()
            os.fsync(self._report.fileno())
        print "\033[32m" + stampedmessage + "\033[0m"

    ## This function logs an error message.
    def error(self, message):
        stampedmessage = timestamp(message, delimiter=" * ")
        if self._doreport:
            self._report.write(stampedmessage + "\n")
            self._report.flush()
            os.fsync(self._report.fileno())
        print "\033[31m" + stampedmessage + "\033[0m"

    ## This function closes the log file, if it was created.
    def finish(self):
        if self._doreport:
            self._report.close()

# -----------------------------------------------------------------
