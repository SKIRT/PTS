#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.looper Functions to perform some given work for all eligible SKIRT-run records
#
# This module offers functions to perform some given work (by invoking a callback function)
# for all SKIRT-run records that are at a given workflow stage and have the 'scheduled' status.

# -----------------------------------------------------------------

# Import standard modules
import os.path

# Import the relevant PTS classes and modules
from . import config as config
from .filelock import FileLock
from .database import Database

# -----------------------------------------------------------------

## This function repeatedly invokes the provided callback function for a SKIRT-run record that is
# at the given workflow stage and has the 'scheduled' status, until there are no more such records.
# Each new record is selected arbitrarily from the eligible records after the previous invocation
# completes. Access to the database is locked so that multiple processes can "pull" work simultaneously.
# Before invoking the callback function, the status of the SKIRT-run record is set to 'running'.
# If the callback function returns normally, the status is updated to 'succeeded'; if an exception
# is raised; the status is updated to 'failed'.
#
# The callback function is passed a single argument containing the SKIRT-run database record
# to be handled. This record offers dictionary-style access to the database fields.
#
def loop(stage, callback):

    while(True):
        # get a record to be processed; return if none are available
        db = Database()
        with FileLock(os.path.join(config.database_path, "database_lock_file.txt"), timeout=10, delay=1):
            with db.transaction():
                records = db.select("stage = ? and status='scheduled'", (stage,))
                if len(records) < 1: return
                record = records[0]
                runid = record["runid"]
                db.updatestatus((runid,), 'running')
        db.close()

        try:
            # invoke the callback function
            print "Processing {} for SKIRT-run {}...".format(stage, runid)
            callback(record)

            # set the runstatus of the database record to 'succeeded'
            db = Database()
            with db.transaction():
                db.updatestatus((runid,), 'succeeded')
            db.close()

        except:
            # set the runstatus of the database record to 'failed'
            db = Database()
            with db.transaction():
                db.updatestatus((runid,), 'failed')
            db.close()
            raise

# -----------------------------------------------------------------
