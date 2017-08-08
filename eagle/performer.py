#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.performer Functions to perform some given work for all eligible SKIRT-run records
#
# This module offers functions to perform some given work (by invoking a callback function)
# for all SKIRT-run records that are at a given workflow stage and have the 'scheduled' status,
# or for all explicitly specified SKIRT-run records.

# -----------------------------------------------------------------

# Import standard modules
import os.path
import time

# Import the relevant PTS classes and modules
from ..core.basics.log import log
from . import config as config
from .filelock import FileLock
from .database import Database
from .skirtrun import runids_in_range

# -----------------------------------------------------------------

## This function repeatedly invokes the provided callback function for a SKIRT-run record that is
# at the given workflow stage and has the 'scheduled' status, until there are no more such records,
# or until the specified run time (given in seconds) has been surpassed.
# Each new record is selected arbitrarily from the eligible records after the previous invocation
# completes. Access to the database is locked so that multiple processes can "pull" work simultaneously.
# Before invoking the callback function, the status of the SKIRT-run record is set to 'running'.
# If the callback function returns normally, the status is updated to 'succeeded'; if an exception
# is raised; the status is updated to 'failed'.
#
# The callback function is passed a single argument containing the SKIRT-run database record
# to be handled. This record offers dictionary-style access to the database fields.
#
def loop(callback, stage, runtime, chunksize=1):
    # loop until runtime has been surpassed
    starttime = time.time()
    while (time.time()-starttime)<runtime:
        # get a chunk of records to be processed; return if none are available
        db = Database()
        with FileLock(os.path.join(config.database_path, "database_lock_file.txt"), timeout=500, delay=5):
            with db.transaction():
                allrecords = db.select("stage = ? and status='scheduled'", (stage,))
                if len(allrecords) < 1: return
                chunkrecords = allrecords[:chunksize]
                db.updatestatus(chunkrecords, 'running')
        db.close()

        try:
            # invoke the callback function
            for record in chunkrecords:
                log.info("Processing {} for SKIRT-run {}...".format(stage, record['runid']))
                callback(record)

            # set the runstatus of the database records to 'succeeded'
            db = Database()
            with db.transaction():
                db.updatestatus(chunkrecords, 'succeeded')
            db.close()

        except:
            # set the runstatus of the database records to 'failed'
            db = Database()
            with db.transaction():
                db.updatestatus(chunkrecords, 'failed')
            db.close()
            raise

# -----------------------------------------------------------------

## This function invokes the provided callback function for each of the specified SKIRT-run records
# (a comma-seperated list of run-ids and/or run-id ranges expressed as two run-ids with a dash in between).
# No locking occurs and the process is performed regardless of the stage and status of the SKIRT-run records.
# Before invoking the callback function, the stage of the SKIRT-run record is set to the specified stage
# with a status of 'running'. If the callback function returns normally, the status is updated to 'succeeded';
# if an exception is raised; the status is updated to 'failed'.
#
# The callback function is passed a single argument containing the SKIRT-run database record
# to be handled. This record offers dictionary-style access to the database fields.
#
def force(callback, stage, runidspec):

    for runid in runids_in_range(runidspec):
        # get a record to be processed; return if none are available
        db = Database()
        with db.transaction():
            records = db.select("runid = ?", (runid,))
            if len(records) != 1: raise ValueError("Runid {} is not in database".format(runid))
            record = records[0]
            db.updatestage((runid,), stage, 'running')
        db.close()

        try:
            # invoke the callback function
            log.info("Processing {} for SKIRT-run {}...".format(stage, runid))
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
