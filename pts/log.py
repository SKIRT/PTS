#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

import datetime
import os
import os.path

# -----------------------------------------------------------------
#  Utility functions
# -----------------------------------------------------------------

def timestamp(message):
    stampedmessage = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S ") + message
    return stampedmessage

# -----------------------------------------------------------------
#  SkirtTestSuite class
# -----------------------------------------------------------------

class Log:
    
    def __init__(self, reportpath="", reportname=""):

        self._doreport = True if reportpath else False
        if self._doreport:
            reportpath = os.path.realpath(os.path.expanduser(reportpath))
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d--%H-%M-%S")
            reportfilepath = os.path.join(reportpath, "report_" + reportname + "_" + timestamp + ".txt")
            self._report = open(reportfilepath, 'w')
    
    def info(self, message):
        
        stampedmessage = timestamp(message)
        
        if self._doreport:
        
            self._report.write(stampedmessage + "\n")
            self._report.flush()    
            os.fsync(self._report.fileno())
        
        print stampedmessage

    def warning(self, message):
        
        stampedmessage = timestamp(message)
        
        if self._doreport:
        
            self._report.write(stampedmessage + "\n")
            self._report.flush()
            os.fsync(self._report.fileno())
        
        print "\033[35m" + stampedmessage + "\033[0m"

    def success(self, message):
        
        stampedmessage = timestamp(message)
        
        if self._doreport:
        
            self._report.write(stampedmessage + "\n")
            self._report.flush()
            os.fsync(self._report.fileno())
        
        print "\033[32m" + stampedmessage + "\033[0m"
    
    def error(self, message):
        
        stampedmessage = timestamp(message)
        
        if self._doreport:
            
            self._report.write(stampedmessage + "\n")
            self._report.flush()
            os.fsync(self._report.fileno())
        
        print "\033[31m" + stampedmessage + "\033[0m"

    def finish(self):
        
        if self._doreport:
            self._report.close()

# -----------------------------------------------------------------
