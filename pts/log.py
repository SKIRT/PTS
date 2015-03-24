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

# Import standard modules
import datetime
import os
import os.path

# -----------------------------------------------------------------

## This private utility function returns the specified message decorated with a timestamp.
def timestamp(message, delimiter="   "):
    stampedmessage = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S.%f")[:-3] + delimiter + message
    return stampedmessage

# -----------------------------------------------------------------

# Define the different logging levels
levels = {"info": 0, "warning": 1, "success": 2, "error": 3}

# -----------------------------------------------------------------

## An instance of the Log class provides logging services to the console, with an optional copy to a text file.
class Log(object):

    ## The constructor takes two optional arguments:
    #
    #  - reportpath: specifies the path of the directory to the file which will receive a copy of all log messages.
    #                If the reportpath argument is empty or missing, no report file is created.
    #  - reportname: specifies the name of the log file (without filename extension; a timestamp is automatically
    #                added to the filename).
    #  - maxlevel: the maximum logging level. The default is "info".
    #
    def __init__(self, reportpath="", reportname="", maxlevel="info"):

        # Check whether a report file should be created
        self._doreport = True if reportpath else False

        # Set the maximum logging level
        self._maxlevel = levels[maxlevel]

        # Open the report file if necessary
        if self._doreport:

            reportpath = os.path.realpath(os.path.expanduser(reportpath))
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d--%H-%M-%S")
            reportfilepath = os.path.join(reportpath, "report_" + reportname + "_" + timestamp + ".txt")
            self._report = open(reportfilepath, 'w')

    ## This function logs an informational message.
    def info(self, message):

        # Do not do anything if the maximum logging level is above the "info" level
        if levels["info"] < self._maxlevel: return

        # Add a time stamp to the message
        stampedmessage = timestamp(message)

        # If required, write the stamped message to the report file
        if self._doreport:

            self._report.write(stampedmessage + "\n")
            self._report.flush()
            os.fsync(self._report.fileno())

        # Print the stamped message
        print stampedmessage

    ## This function logs a warning message.
    def warning(self, message):

        # Do not do anything if the maximum logging level is above the "warning" level
        if levels["warning"] < self._maxlevel: return

        # Add a time stamp to the message
        stampedmessage = timestamp(message, delimiter=" ! ")

        # If required, write the stamped message to the report file
        if self._doreport:

            self._report.write(stampedmessage + "\n")
            self._report.flush()
            os.fsync(self._report.fileno())

        # Print the stamped message
        print "\033[35m" + stampedmessage + "\033[0m"

    ## This function logs a success message.
    def success(self, message):

        # Do not do anything if the maximum logging level is above the "success" level
        if levels["success"] < self._maxlevel: return

        # Add a time stamp to the message
        stampedmessage = timestamp(message, delimiter=" - ")

        # If required, write the stamped message to the report file
        if self._doreport:

            self._report.write(stampedmessage + "\n")
            self._report.flush()
            os.fsync(self._report.fileno())

        # Print the stamped message
        print "\033[32m" + stampedmessage + "\033[0m"

    ## This function logs an error message.
    def error(self, message):

        # Do not do anything if the maximum logging level is above the "error" level
        if levels["warning"] < self._maxlevel: return

        # Add a time stamp to the message
        stampedmessage = timestamp(message, delimiter=" * *** Error: ")

        # If required, write the stamped message to the report file
        if self._doreport:

            self._report.write(stampedmessage + "\n")
            self._report.flush()
            os.fsync(self._report.fileno())

        # Print the stamped message
        print "\033[31m" + stampedmessage + "\033[0m"

    ## This function closes the log file, if necessary
    def finish(self):

        # If a report was written, close the corresponding file
        if self._doreport: self._report.close()

# -----------------------------------------------------------------
