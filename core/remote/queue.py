#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.queue Contains the SimulationQueue class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import sqlite3

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs

# -----------------------------------------------------------------

## This variable holds an enumeration for the supported database field types
fieldtype_enum = ('text', 'numeric')

## This variable holds an enumeration for the values of the runstatus field
runstatus_enum = ('inserted', 'scheduled', 'running', 'failed', 'completed', 'archived', 'removed')

# -----------------------------------------------------------------

# Initializing the queue
#import eagle.database
#db = eagle.database.Database()
#db.createtable()
#db.close()
#

# Adding to queue
# backup the data base
#database.backup()

# insert a new record into the database for each selected galaxy
#db = database.Database()
#with db.transaction():
#    for record in records:
#        db.insert(label, eaglesim, 0, record["galaxyid"], skitemplate)
#db.close()

## OTHER STUFF
# # construct the index on galayid, and add a galaxy-id property
#self._ids = { }       # key = galaxy id, value = index in collection
#self._info["galaxy_id"] = self._info["skirt_run_id"].copy()
#db = Database()
#index = 0
#for runid in self._info["skirt_run_id"]:
#    galaxy_id = db.select("runid=?", (runid,))[0]["galaxyid"]
#    self._ids[galaxy_id] = index
#    self._info["galaxy_id"][index] = galaxy_id
#    index+=1
#db.close()

# # set the runstatus of the database record to 'running'
#db = Database()
#with db.transaction():
#    db.updatestatus((runid,), 'running')
#db.close()

# set the runstatus of the database record to 'completed'
#db = Database()
#with db.transaction():
#    db.updatestatus((runid,), 'completed')
#db.close()

# set the runstatus of the database record to 'failed'
#db = Database()
#with db.transaction():
#    db.updatestatus((runid,), 'failed')
#db.close()
#raise


# open the database
#db = Database()

# get the records corresponding to the run-id sequence from the database
#records = db.select("runid in (" + ",".join(map(str,runids)) + ")")

# update the records to indicate that a run has been scheduled and submit the job;
# do this within a single transaction context to ensure that the scheduled job sees the updated records
#print "Submitting job to queue", config.queue, "for run-ids", runids
#with db.transaction():
#    db.updatestatus(runids, 'scheduled')
#    db.updatefield(runids, 'queue', config.queue)
#    subprocess.call(("bsub",), stdin=open(jobscriptname))

#print "Please manually perform scheduled run-ids", runids
#with db.transaction():
#            db.updatestatus(runids, 'scheduled')
#            db.updatefield(runids, 'queue', config.hostname)

# close the database
#db.close()

# -----------------------------------------------------------------

class SimulationQueue(object):

    """
    This class ...
    """

    def __init__(self, name, path):

        """
        This function ...
        :param name:
        :param path:
        :return:
        """

        #self.remote = remote

        # The queue name
        self.name = name

        #self.pid = None
        #self.pip_file_name = None

        # Make connection
        self._con = sqlite3.connect(path)
        self._con.row_factory = sqlite3.Row
        self._con.text_factory = sqlite3.OptimizedUnicode

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        name = fs.strip_extension(fs.name(path))
        return cls(name, path)

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        self.close()

    # -----------------------------------------------------------------

    def close(self):

        """
        This function closes the connection to the database. Any uncommited changes are lost.
        You can no longer use the database after calling this function.
        :return:
        """

        self._con.close()
        self._con = None

    # -----------------------------------------------------------------

    ## This function creates an empty skirtruns table in a fresh SKIRT-runs database. It is \em not intended for use in
    # production code. The change is automatically committed.
    def createtable(self):
        # use "integer primary key autoincrement" for the automatically assigned unique record id
        # use "text" affinity for strings and dates;
        # use "numeric" affinity for integers and reals (except the unique record id)
        #self._con.execute('''create table skirtruns (
        #                     runid integer primary key autoincrement,
        #                     username text,
        #                     label text,
        #                     runstatus text,
        #                     statusdate text,
        #                     queue text,
        #                     eaglesim text,
        #                     redshift numeric,
        #                     galaxyid numeric,
        #                     skitemplate text
        #                  )''')

        self._con.execute('''create table simulations (
                          runid integer primary key autoincrement,
                          name text,
                          ski_path text,
                          input_path text,
                          output_path text,
                          status text,
                          nprocesses integer,
                          nthreads integer,
                          data_parallel integer,

                          )''')

    # -----------------------------------------------------------------

    ## This function adds a field to the main table in the SKIRT-runs database. It is \em not intended for use in
    # production code. The change is automatically committed. The field type should be 'text' for strings and dates,
    # or 'numeric' for integers and reals. The initial value should have an appropriate type, or be omitted.
    def addfield(self, fieldname, fieldtype='text', initialvalue=None):
        if not fieldtype in fieldtype_enum: raise ValueError("unsupported field type: " + fieldtype)
        self._con.execute("alter table skirtruns add column " + fieldname + " " + fieldtype)
        if initialvalue != None:
            self._con.execute("update skirtruns set " + fieldname + " = ?", (initialvalue,))
            self._con.commit()

    ## This function executes an arbitrary SQL statement on the SKIRT-runs database and returns an sqlite3.Cursor
    # object representing the result. It is \em not intended for use in production code.
    def execute(self, sql, params=None):
        return self._con.execute(sql, params if params != None else ())

    # -----------------------------------------------------------------

    def maxrunid(self):

        """
        This function returns the largest run-id currently in the database (i.e. the run-id of the most
        recently inserted row).
        """

        cursor = self._con.execute("select max(runid) from skirtruns")
        return cursor.fetchone()[0]

    # -----------------------------------------------------------------

    # -------- updating --------

    ## This function returns an object that can be used as the context manager of a with block, automatically
    # committing any changes at the end of the block, or performing a rollback if an error is raised within the block.
    # For example:
    #
    # \verbatim
    # with db.transaction():
    #    for item in itemlist: db.insert(...item...)
    # \endverbatim
    #
    def transaction(self):

        return self._con

    # -----------------------------------------------------------------

    def commit(self):

        """
        This function commits any changes to the database. In production code you should prefer the use of a with block
        (see the context() function). So this function is primarily for interactive use.
        :return:
        """

        return self._con.commit()

    # -----------------------------------------------------------------

    def insert(self, label, eaglesim, redshift, galaxyid, skitemplate,
               align=1, seed=1, numpp=5e5, deltamax=3e-6, fdust=0.3, fpdr=0.15, keep=1):

        """
        This function inserts a new record into the database with the specified field values. The function
        automatically determines values for runid (unique serial nr), username (from config), runstatus ('inserted')
        and statusdate (now). The change is \em not committed.
        :param label:
        :param eaglesim:
        :param redshift:
        :param galaxyid:
        :param skitemplate:
        :param align:
        :param seed:
        :param numpp:
        :param deltamax:
        :param fdust:
        :param fpdr:
        :param keep:
        :return:
        """

        username = config.username
        runstatus = runstatus_enum[0]
        statusdate = config.timestamp()
        self._con.execute('''insert into skirtruns (username, runstatus, statusdate,
                                                    label, eaglesim, redshift, galaxyid, skitemplate,
                                                    align, seed, numpp, deltamax, fdust, fpdr, keep)
                             values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',
                          (username, runstatus, statusdate,
                           label, eaglesim, redshift, galaxyid, skitemplate,
                           align, seed, numpp, deltamax, fdust, fpdr, keep))

# -----------------------------------------------------------------
