#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package eagle.database Functions to help create, update and query the SKIRT-runs database.
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
#  - the current status of the SKIRT simulation (e.g. to be scheduled, running or finished)
#  - the EAGLE data fed to this SKIRT simulation (e.g. snapshot, galaxy ID)
#  - the SKIRT parameters used for the simulation
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
#<TR><TD>username</TD>      <TD>Text</TD>
#                           <TD>The name of the user responsible for inserting this SKIRT run into the database,
#                               for example "pcamps".</TD></TR>
#<TR><TD>label</TD>         <TD>Text</TD>
#                           <TD>A short user-specified identifier used to indicate records that belong together,
#                               facilitating querying or updating those records as a group.</TD></TR>
#<TR><TD>runstatus</TD>     <TD>Text</TD>
#                           <TD>An identifier for the current status of this SKIRT run. In the current implementation,
#                               the value should be one of "inserted", "scheduled", "running", "failed",
#                               "completed", "archived", or "removed".</TD></TR>
#<TR><TD>statusdate</TD>    <TD>Date</TD>
#                           <TD>The time and date on which this record's run status was most recently changed.</TD></TR>
#<TR><TD>queue</TD>         <TD>Text</TD>
#                           <TD>A short identifier for the queue or computing system on which the SKIRT run is being
#                               performed or has been performed. The value of this field is relevant only when the
#                               runstatus is 'scheduled' or later.</TD></TR>
#<TR><TD>eaglesim</TD>      <TD>Text</TD>
#                           <TD>A short name for the EAGLE simulation from which data is obtained for this SKIRT run,
#                               for example 'Ref100Mpc'. This value is used as a key into the \em eagledata_path
#                               dictionary provided by the eagle.config package.</TD></TR>
#<TR><TD>redshift</TD>      <TD>Real</TD>
#                           <TD>The redshift of the EAGLE snapshot from which data is obtained for this SKIRT run.
#                               This value is used to determine the appropriate filenames in the EAGLE data.</TD></TR>
#<TR><TD>galaxyid</TD>      <TD>Integer</TD>
#                           <TD>The identifier (GalaxyID) in the public EAGLE database for the galaxy
#                               used for this SKIRT run.</TD></TR>
#<TR><TD>skitemplate</TD>   <TD>Text</TD>
#                           <TD>The name (without extension) of the ski file used as a template for this SKIRT run.
#                               The actual ski file is automatically derived from the template. The templates are
#                               located in the directory indicated by the variable \em templates_path
#                               provided by the eagle.config package.</TD></TR>
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
import eagle.config as config

# -----------------------------------------------------------------

## This variable holds an enumeration for the supported database field types
fieldtype_enum = ('text', 'numeric')

## This variable holds an enumeration for the values of the runstatus field
runstatus_enum = ('inserted', 'scheduled', 'running', 'failed', 'completed', 'archived', 'removed')

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
        self._con = sqlite3.connect(filepath)
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

    ## This function returns a sequence of row objects representing the set of database records selected by
    # the specified SQL \em where expression. The optional \em params argument provides a sequence of values
    # that will replace the corresponding question marks in the \em where expression.
    # Each row object contains all fields present in the database. To extract a value, use the field name as
    # an index to the row object; never rely on the order of the fields because this may change over time.
    # For example:
    #
    #\verbatim
    #for row in db.select("label=?", (label, )):
    #    print row["galaxyid"], row["runstatus"]
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
    #db.show(db.select("runstatus='inserted'"))
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
    # automatically determines values for runid (unique serial nr), username (from config), runstatus ('inserted')
    # and statusdate (now). The change is \em not committed.
    def insert(self, label, eaglesim, redshift, galaxyid, skitemplate):
        username = config.username
        runstatus = runstatus_enum[0]
        statusdate = config.timestamp()
        self._con.execute('''insert into skirtruns (username, runstatus, statusdate,
                                    label, eaglesim, redshift, galaxyid, skitemplate)
                             values (?,?,?,?,?,?,?,?)''',
                           (username, runstatus, statusdate, label, eaglesim, redshift, galaxyid, skitemplate)
                          )

    ## This function updates the value of a field in the records specified through a list of run-id's.
    # The change is \em not committed.
    # The \em runids argument provides a sequence of run-id's; each item in the sequence must be either a number
    # or an object that returns a number when indexed with the 'runid' string. It is not possible to update the
    # runid field. To update the runstatus and statusdate fields, use the updatestatus() function.
    def updatefield(self, runids, fieldname, value):
        if 'status' in fieldname: raise ValueError("Use the updatestatus() function to update status fields")
        for runid in runids:
            self._con.execute("update skirtruns set " + fieldname + " = ? where runid = ?", (value,cleanrunid(runid)))

    ## This function updates the value of the runstatus field, and stores a new timestamp in
    # in the statusdate field, in the records specified through a list of run-id's.
    # The change is \em not committed.
    # The \em runids argument provides a sequence of runid's; each item in the sequence must be either a number
    # or an object that returns a number when indexed with the 'runid' string.
    def updatestatus(self, runids, runstatus):
        if not runstatus in runstatus_enum:
            raise ValueError("Unsupported runstatus value: " + runstatus)
        statusdate = config.timestamp()
        for runid in runids:
            self._con.execute("update skirtruns set runstatus = ?, statusdate = ? where runid = ?",
                               (runstatus,statusdate,cleanrunid(runid)))

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
                             username text,
                             label text,
                             runstatus text,
                             statusdate text,
                             queue text,
                             eaglesim text,
                             redshift numeric,
                             galaxyid numeric,
                             skitemplate text
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
