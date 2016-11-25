#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.database Functions to help create, update and query the SKIRT-runs database.
#
# This module describes the SKIRT-runs database, documents how to access it from production code, and
# defines a class to help manage, update and query the database. Some of the functions are intended
# for interactive use during database maintenance operations (e.g. adding an extra field) and thus
# serve as an example rather than as a polished tool for use in production code.
#
# <B>Objective</B>
#
# The SKIRT-runs database keeps track of SKIRT simulations performed on EAGLE simulation data.
# Each record in the database represents a single SKIRT simulation and holds information on:
#  - the current processing stage (e.g. extract, simulate) and status (e.g. running, succeeded)
#  - the EAGLE data fed to this SKIRT simulation (e.g. snapshot, galaxy ID)
#  - the SKIRT parameters used for the simulation (e.g. ski template)
#
# <B>Fields</B>
#
# The following table lists the fields (columns) for each record (row) in the main (and so far only)
# table in the database.
#
#<TABLE>
#<TR><TD><B>Field</B></TD>  <TD><B>Type</B></TD>  <TD><B>Description of value</B></TD></TR>
#<TR><TD>runid</TD>         <TD>Integer</TD>
#                           <TD>An automatically assigned unique identifier for this SKIRT run. This value is used
#                               as the name of a directory containing the SKIRT in/out data for this run. These
#                               directories are located in the directory indicated by the variable \em results_path
#                               provided by the eagle.config package.</TD></TR>
#<TR><TD>label</TD>         <TD>Text</TD>
#                           <TD>A short user-specified identifier used to indicate records that belong together,
#                               facilitating querying or updating those records as a set.</TD></TR>
#<TR><TD>stage</TD>         <TD>Text</TD>
#                           <TD>An identifier for the current stage of this SKIRT run. The value should be one of
#                               "insert", "extract", "simulate", "observe", "store", "completed", or "removed".
#                               </TD></TR>
#<TR><TD>status</TD>        <TD>Text</TD>
#                           <TD>An identifier for the execution status of the current stage of this SKIRT run.
#                               The value should be one of "scheduled", "running", "failed", or "succeeded".</TD></TR>
#<TR><TD>statusdate</TD>    <TD>Date</TD>
#                           <TD>The time and date on which this record's stage and/or status were most recently
#                               changed.</TD></TR>
#<TR><TD>eaglesim</TD>      <TD>Text</TD>
#                           <TD>The name of the EAGLE simulation from which data is obtained for this SKIRT run,
#                               for example 'RefL0100N1504'. This value is used as a key into the \em eagledata_path
#                               dictionary provided by the eagle.config package.</TD></TR>
#<TR><TD>snaptag</TD>       <TD>Integer</TD>
#                           <TD>The tag (0-28) of the EAGLE snapshot (indicating redshift 20 through 0) from which
#                               data is obtained for this SKIRT run.
#                               This value is used to determine the appropriate filenames in the EAGLE data.</TD></TR>
#<TR><TD>galaxyid</TD>      <TD>Integer</TD>
#                           <TD>The identifier (GalaxyID) in the public EAGLE database for the galaxy
#                               used for this SKIRT run.</TD></TR>
#<TR><TD>groupnr</TD>       <TD>Integer</TD>
#                           <TD>The group number identifying the particles of this galaxy in the snapshot.</TD></TR>
#<TR><TD>subgroupnr</TD>    <TD>Integer</TD>
#                           <TD>The subgroup number identifying the particles of this galaxy in the snapshot.</TD></TR>
#<TR><TD>starmass</TD>      <TD>Real</TD>
#                           <TD>The intrinsic stellar mass of this galaxy, in solar masses.</TD></TR>
#<TR><TD>copx</TD>          <TD>Real</TD>
#                           <TD>The x-coordinate of the center of potential of this galaxy, in cMpc.</TD></TR>
#<TR><TD>copy</TD>          <TD>Real</TD>
#                           <TD>The y-coordinate of the center of potential of this galaxy, in cMpc.</TD></TR>
#<TR><TD>copz</TD>          <TD>Real</TD>
#                           <TD>The z-coordinate of the center of potential of this galaxy, in cMpc.</TD></TR>
#<TR><TD>skitemplate</TD>   <TD>Text</TD>
#                           <TD>The name (without extension) of the ski file used as a template for this SKIRT run.
#                               The actual ski file is automatically derived from the template. The templates are
#                               located in the directory indicated by the variable \em templates_path
#                               provided by the eagle.config package.</TD></TR>
#<TR><TD>numpp</TD>         <TD>Real</TD>
#                           <TD>The number of photon packages to be launched for this SKIRT run.</TD></TR>
#<TR><TD>deltamax</TD>      <TD>Real</TD>
#                           <TD>The maximum mass fraction of the dust grid used for this SKIRT run.</TD></TR>
#</TABLE>
#
# <B>Implementation</B>
#
# The SKIRT-runs database is implemented using the sqlite3 package (part of standard Python), which is based on
# the SQLite library version 3.7.3. See www.sqlite.org for information on the supported SQL dialect and its syntax.
#

# -----------------------------------------------------------------

import os.path
import shutil
import sqlite3

from . import config as config

# -----------------------------------------------------------------

## This variable holds an enumeration for the supported database field types
fieldtype_enum = ('text', 'numeric')

## These variables hold enumerations for the values of the stage and status fields
stage_enum = ('insert', 'extract', 'simulate', 'observe', 'store', 'completed', 'removed')
status_enum = ('scheduled', 'running', 'failed', 'succeeded')

# -----------------------------------------------------------------

## This function copies the complete database to the backup location (under a time-stamped name).
# You should call this function before opening the database with the purpose of updating it.
def backup():
    original = os.path.join(config.database_path, "SKIRT-runs.db")
    timestamp = config.timestamp()
    backup = os.path.join(config.backup_path, "SKIRT-runs_backup_"+timestamp+".db")
    shutil.copyfile(original, backup)

# -----------------------------------------------------------------
#  Database class
# -----------------------------------------------------------------

## This class represents the SKIRT-runs database and provides functions to help manage, update and query the database.
# Some of these functions are intended for interactive use during database maintenance operations (e.g. adding an
# extra field) and thus serve as an example rather than as a polished tool for use in production code.
#
class Database:

    # -------- opening and closing --------

    ## The constructor opens the SKIRT-runs database and internally retains a connection to it. Once the connection is
    # no longer needed, you should close it by calling the close() function on the database object. For example:
    #
    #\verbatim
    #import eagle.database
    #db = eagle.database.Database()
    # ...
    #db.close()
    #\endverbatim
    #
    # The optional filepath argument allows opening a SKIRT-runs database with a nonstandard location and name.
    # In almost all cases this argument should be omitted so that the constructor opens the standard database.
    def __init__(self, filepath=None):
        if filepath==None:
            filepath = os.path.join(config.database_path, "SKIRT-runs.db")
        self._con = sqlite3.connect(filepath, timeout=60)
        self._con.row_factory = sqlite3.Row
        self._con.text_factory = sqlite3.OptimizedUnicode

    ## This function closes the connection to the database. Any uncommited changes are lost.
    # You can no longer use the database after calling this function.
    def close(self):
        self._con.close()
        self._con = None;

    # -------- querying --------

    ## This function prints information on the table declaration(s) in the SKIRT-runs database.
    def showinfo(self):
        print "Info for database SKIRT-runs"
        for table in self._con.execute('''select name from sqlite_master
                                          where type='table' and not name glob 'sqlite*' order by name'''):
            tablename = table['name']
            print "  Table: " + tablename
            for index in self._con.execute('''select name from sqlite_master
                                              where type='index' and tbl_name=? order by name''', (tablename,)):
                print "    Index: " + index['name']
            for col in self._con.execute("pragma table_info('" + tablename + "')"):
                print "    Column: " + col['name'] + " - " + col['type']

    ## This function returns the largest run-id currently in the database (i.e. the run-id of the most
    # recently inserted row).
    def maxrunid(self):
        cursor = self._con.execute("select max(runid) from skirtruns")
        id = cursor.fetchone()[0]
        return id if id is not None else 0

    ## This function returns a sequence of row objects representing the set of database records selected by
    # the specified SQL \em where expression. The optional \em params argument provides a sequence of values
    # that will replace the corresponding question marks in the \em where expression.
    # Each row object contains all fields present in the database. To extract a value, use the field name as
    # an index to the row object; never rely on the order of the fields because this may change over time.
    # For example:
    #
    #\verbatim
    #for row in db.select("label=?", (label, )):
    #    print row["galaxyid"], row["stage"], row["status"]
    #\endverbatim
    #
    def select(self, where, params=None):
        cursor = self._con.execute("select * from skirtruns where " + where, params if params!=None else ())
        return cursor.fetchall()

    ## This function prints the contents of a sequence of row objects for visual inspection. If the \em refetch
    # argument is False or missing, the function simply prints the contents of the specified row objects.
    # For example:
    #
    #\verbatim
    #db.show(db.select("stage='insert'"))
    #\endverbatim
    #
    # If the \em refetch argument is True, the function retrieves the data from the data base based on the run-id's
    # specified in the \em rows argument. In this case, the \em rows argument provides a sequence of runid's;
    # each item in the sequence must be either a number or an object that returns a number when indexed with
    # the 'runid' string. For example:
    #
    #\verbatim
    #db.show((4,78,23), refetch=True)
    #\endverbatim
    #
    def show(self, rows, refetch=False):
        if refetch:
            runids = rows
            rows = [ ]
            for runid in runids:
                rows += self.select("runid=?", (cleanrunid(runid),))

        if len(rows)==0:
            print "There are no records to be shown"
        else:
            print rows[0].keys()
            for row in rows: print row
            print "{} records were shown".format(len(rows))

    # -------- updating --------

    ## This function returns an object that can be used as the context manager of a with block, automatically
    # committing any changes at the end of the block, or performing a rollback if an error is raised within the block.
    # For example:
    #
    #\verbatim
    #with db.transaction():
    #    for item in itemlist: db.insert(...item...)
    #\endverbatim
    #
    def transaction(self):
        return self._con;

    ## This function commits any changes to the database. In production code you should prefer the use of a with block
    # (see the context() function). So this function is primarily for interactive use.
    def commit(self):
        return self._con.commit();

    ## This function inserts a new record into the database with the specified field values. The function
    # automatically determines values for runid (unique serial nr), stage('insert'), status ('succeeded'),
    # and statusdate (now). The change is \em not committed.
    def insert(self, label, eaglesim, snaptag, galaxyid, groupnr, subgroupnr, starmass, copx, copy, copz,
                     skitemplate, numpp, deltamax):
        stage = 'insert'
        status = 'succeeded'
        statusdate = config.timestamp()
        self._con.execute('''insert into skirtruns (label, stage, status, statusdate, eaglesim, snaptag,
                                                    galaxyid, groupnr, subgroupnr, starmass, copx, copy, copz,
                                                    skitemplate, numpp, deltamax)
                             values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',
                           (label, stage, status, statusdate, eaglesim, snaptag,
                            galaxyid, groupnr, subgroupnr, starmass, copx, copy, copz,
                            skitemplate, numpp, deltamax)
                          )

    ## This function updates the value of a field in the records specified through a list of run-id's.
    # The change is \em not committed.
    # The \em runids argument provides a sequence of run-id's; each item in the sequence must be either a number
    # or an object that returns a number when indexed with the 'runid' string. It is not possible to update the runid
    # field. To update the stage, status and statusdate fields, use the updatestatus() or updatestatus() functions.
    def updatefield(self, runids, fieldname, value):
        if fieldname in ('stage','status'):
            raise ValueError("Use the updatestage() or updatestatus() functions to update these fields")
        for runid in runids:
            self._con.execute("update skirtruns set " + fieldname + " = ? where runid = ?", (value,cleanrunid(runid)))

    ## This function updates the value of the stage and status fields, and stores a new timestamp in the statusdate
    # field, for the specified records. The default new value for the status field is 'scheduled', but another value
    # can be specified. The change is \em not committed.
    # The \em runids argument provides a sequence of runid's; each item in the sequence must be either a number
    # or an object that returns a number when indexed with the 'runid' string.
    def updatestage(self, runids, stage, status='scheduled'):
        if not stage in stage_enum: raise ValueError("Unsupported stage value: " + stage)
        if not status in status_enum: raise ValueError("Unsupported status value: " + status)
        statusdate = config.timestamp()
        for runid in runids:
            self._con.execute("update skirtruns set stage = ?, status = ?, statusdate = ? where runid = ?",
                               (stage,status,statusdate,cleanrunid(runid)))

    ## This function updates the value of the status field (without touching the stage field), and stores a new
    # timestamp in the statusdate field, for the specified records. The change is \em not committed.
    # The \em runids argument provides a sequence of runid's; each item in the sequence must be either a number
    # or an object that returns a number when indexed with the 'runid' string.
    def updatestatus(self, runids, status):
        if not status in status_enum: raise ValueError("Unsupported status value: " + status)
        statusdate = config.timestamp()
        for runid in runids:
            self._con.execute("update skirtruns set status = ?, statusdate = ? where runid = ?",
                               (status,statusdate,cleanrunid(runid)))

    ## This function updates all fields of the database record for the specified run-id to the values contained
    # in the specified row object. The row object should be obtained through the select() function from a database
    # with the same fields as the target database. The value of the 'runid' field in the row is ignored.
    # The function returns True if any of the field values were modified, and False if no changes were made
    # (because all fields already had the proper values). However, in both cases the change is \em not committed.
    def updaterow(self, runid, newrow):
        # get the row for this run-id
        cursor = self._con.execute("select * from skirtruns where runid = ?", (runid,))
        oldrows = cursor.fetchall()
        if len(oldrows) != 1: raise ValueError("The specified run-id does not match a database record: " + str(runid))
        oldrow = oldrows[0]

        # update each field if needed
        modified = False
        for fieldname in oldrow.keys():
            if fieldname != 'runid' and oldrow[fieldname] != newrow[fieldname]:
                self._con.execute("update skirtruns set " + fieldname + " = ? where runid = ?",
                                  (newrow[fieldname],runid))
                modified = True
        return modified

    # -------- maintaining --------

    ## This function creates an empty skirtruns table in a fresh SKIRT-runs database. It is \em not intended for use in
    # production code. The change is automatically committed.
    def createtable(self):
        # use "integer primary key autoincrement" for the automatically assigned unique record id
        # use "text" affinity for strings and dates;
        # use "numeric" affinity for integers and reals (except the unique record id)
        self._con.execute('''create table skirtruns (
                             runid integer primary key autoincrement,
                             label text,
                             stage text,
                             status text,
                             statusdate text,
                             eaglesim text,
                             snaptag numeric,
                             galaxyid numeric,
                             groupnr numeric,
                             subgroupnr numeric,
                             starmass numeric,
                             copx numeric,
                             copy numeric,
                             copz numeric,
                             skitemplate text,
                             numpp numeric,
                             deltamax numeric
                          )''')

    ## This function adds a field to the main table in the SKIRT-runs database. It is \em not intended for use in
    # production code. The change is automatically committed. The field type should be 'text' for strings and dates,
    # or 'numeric' for integers and reals. The initial value should have an appropriate type, or be omitted.
    def addfield(self, fieldname, fieldtype='text', initialvalue=None):
        if not fieldtype in fieldtype_enum: raise ValueError("unsupported field type: "+fieldtype)
        self._con.execute("alter table skirtruns add column " + fieldname + " " + fieldtype)
        if initialvalue!=None:
            self._con.execute("update skirtruns set " + fieldname + " = ?", (initialvalue,))
            self._con.commit()

    ## This function executes an arbitrary SQL statement on the SKIRT-runs database and returns an sqlite3.Cursor
    # object representing the result. It is \em not intended for use in production code.
    def execute(self, sql, params=None):
        return self._con.execute(sql, params if params!=None else ())

# -----------------------------------------------------------------

## This helper function extracts a run-id from its argument, which must be either a number
# or an object that returns a number when indexed with the 'runid' string.
def cleanrunid(runid):
    try: clean = int(runid)
    except: clean = 0
    if clean<=0: clean = int(runid['runid'])
    return clean

# -----------------------------------------------------------------
