#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.dbadapters Pyevolve have a feature in which you can save the statistics of every
#  generation in a database, file or call an URL with the statistics as param.
#  You can use the database to plot evolution statistics graphs later. In this
#  module, you'll find the adapters above cited.
#
# warning: the use the of a DB Adapter can reduce the performance of the Genetic Algorithm.
# seealso:
#   Method :meth:`GSimpleGA.GSimpleGA.setDBAdapter`
#      DB Adapters are set in the GSimpleGA Class.

# -----------------------------------------------------------------

# Import standard modules
from abc import ABCMeta, abstractmethod
import types
import datetime

# Import other evolve modules
from . import statistics
from . import constants
from . import utils

# Import the relevant PTS classes and modules
from ...core.basics.log import log

# -----------------------------------------------------------------

class DataBaseAdapter(object):

    """
    DataBaseAdapter Class - The base class for all DB Adapters
    If you want to create your own DB Adapter, you must subclass this
    class.
    :param frequency: the the generational dump frequency
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, frequency, identify, name):

        """
        The class constructor
        """

        self.statsGenFreq = frequency

        if identify is None: self.identify = datetime.datetime.now().strftime("%d/%m/%y-%H:%M")
        else: self.identify = identify

        # Set its own name
        self.name = name

    # -----------------------------------------------------------------

    def setIdentify(self, identify):

        """
        Sets the identify of the statistics
        :param identify: the id string
        """

        if identify is None: self.identify = datetime.datetime.now().strftime("%d/%m/%y-%H:%M")
        else: self.identify = identify

    # -----------------------------------------------------------------

    def getIdentify(self):

        """
        Return the statistics identify
        :rtype: identify string
        """

        return self.identify

    # -----------------------------------------------------------------

    def getStatsGenFreq(self):

        """
        Returns the frequency of statistical dump
        :rtype: the generation interval of statistical dump
        """

        return self.statsGenFreq

    # -----------------------------------------------------------------

    def setStatsGenFreq(self, statsGenFreq):

        """
        Set the frequency of statistical dump
        :param statsGenFreq: the generation interval of statistical dump
        """

        self.statsGenFreq = statsGenFreq

    # -----------------------------------------------------------------

    @abstractmethod
    def open(self, ga_engine):

        """
        This method is called one time to do the initialization of
        the DB Adapter
        :param ga_engine: the GA Engine
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def commit_and_close(self):

        """
        This method is called at the end of the evolution, to closes the
        DB Adapter and commit the changes
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def insert(self, ga_engine):

        """
        Insert the stats
        :param ga_engine: the GA Engine
        """

        pass

# -----------------------------------------------------------------

class PopulationsFile(DataBaseAdapter):

    """
    This class ...
    """

    def __init__(self, filepath=constants.CDefPopulationsFileName, identify=None,
                 frequency=constants.CDefPopulationsStatsGenFreq, reset=True, name=constants.CDefPopulationsName):

        """
        This function ...
        :param filepath: 
        :param identify: 
        :param frequency:
        :param reset:
        """

        # Call the constructor of the base class
        super(PopulationsFile, self).__init__(frequency, identify, name)

        # Set properties
        self.filepath = filepath
        self.reset = reset

        # THe file handle
        self.handle = None

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return: 
        """

        ret = "PopulationsFile Adapter [File='%s', identify='%s']" % (self.filepath, self.getIdentify())
        return ret

    # -----------------------------------------------------------------

    def open(self, ga_engine):

        """
        Open the populations file
        """

        # Debugging
        log.debug("Opening the " + self.name + " ...")

        # Debugging
        #log.debug("Opening the populations file to dump genome properties '%s'", self.filepath)
        open_mode = 'w' if self.reset else 'a'

        self.handle = open(self.filepath, open_mode)

    # -----------------------------------------------------------------

    def close(self):

        """ 
        Closes the populations file handle
        """

        log.debug("Closing the populations file [%s]", self.filepath)

        if self.handle: self.handle.close()

    # -----------------------------------------------------------------

    def commit_and_close(self):

        """
        Commits and closes
        """

        self.close()

    # -----------------------------------------------------------------

    def insert(self, engine):

        """
        Inserts the genomes into the populations file
        :param engine: the GeneticEngine
        """

        # Get the generation index (0 for intitial -> n for Generation n-1)
        generation = engine.getCurrentGeneration()

        # Get the internal population (the survivors (newborns + elitism))
        population = engine.get_population()

        # Loop over the population's individuals
        for key in population.keys:

            line = [self.getIdentify(), str(generation)]

            # Get the individual
            individual = population[key]

            # Add entry to the list
            key = str(key)
            genes_string = str(individual.genomeList)  # .genomeList works for all G1D genomes (list, binary string, ...)

            # Add to the line
            line.append(key)
            line.append(genes_string)

            # Loop over the entries
            self.handle.write(" ".join(line) + "\n")

            # FLUSH: IMPORTANT SINCE WE WANT TO GET GENOME INFORMATION FROM INTIIAL GENERATION WHEN WE HAVE JUST SCORED IT AND GENERATED GENERATION0
            self.handle.flush()

# -----------------------------------------------------------------

class DBFileCSV(DataBaseAdapter):

    """
    DBFileCSV Class - Adapter to dump statistics in CSV format
    Inheritance diagram for :class:`DBAdapters.DBFileCSV`:
    .. inheritance-diagram:: DBAdapters.DBFileCSV

    Example:
      >>> adapter = DBFileCSV(filename="file.csv", identify="run_01", frequency = 1, reset=True)
      :param filename: the CSV filename
      :param identify: the identify of the run
      :param frequency: the generational dump frequency
      :param reset: if True, the file old data will be overwrite with the new
    """

    def __init__(self, filename=constants.CDefCSVFileName, identify=None,
                frequency=constants.CDefCSVFileStatsGenFreq, reset=True, name=constants.CDefCSVName):

        """
        The creator of DBFileCSV Class
        """

        # Call the constructor of the base class
        super(DBFileCSV, self).__init__(frequency, identify, name)

        # The CSV module
        self.csvmod = None

        self.filename = filename
        self.csvWriter = None
        self.fHandle = None
        self.reset = reset

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        The string representation of adapter
        """

        ret = "DBFileCSV DB Adapter [File='%s', identify='%s']" % (self.filename, self.getIdentify())
        return ret

    # -----------------------------------------------------------------

    def open(self, ga_engine):

        """
        Open the CSV file or creates a new file
        :param ga_engine: the GA Engine
        .. versionchanged:: 0.6
         The method now receives the *ga_engine* parameter.
        """

        # Debugging
        log.debug("Opening the " + self.name + " ...")

        if self.csvmod is None:
         log.debug("Loading the csv module...")
         self.csvmod = utils.importSpecial("csv")

        # Debugging
        log.debug("Opening the CSV file to dump statistics '%s'", self.filename)
        open_mode = 'w' if self.reset else 'a'

        self.fHandle = open(self.filename, open_mode)
        self.csvWriter = self.csvmod.writer(self.fHandle, delimiter=',')

    # -----------------------------------------------------------------

    def close(self):

        """
        Closes the CSV file handle
        """

        log.debug("Closing the " + self.name + " ...")
        #log.debug("Closing the CSV file [%s]", self.filename)

        if self.fHandle: self.fHandle.close()

    # -----------------------------------------------------------------

    def commit_and_close(self):

        """
        Commits and closes
        """

        self.close()

    # -----------------------------------------------------------------

    def insert(self, ga_engine):

        """
        Inserts the stats into the CSV file
        :param ga_engine: the GA Engine
        .. versionchanged:: 0.6
         The method now receives the *ga_engine* parameter.
        """

        stats = ga_engine.getStatistics()
        generation = ga_engine.getCurrentGeneration()
        line = [self.getIdentify(), generation]
        line.extend(stats.asTuple())
        self.csvWriter.writerow(line)

# -----------------------------------------------------------------

class DBURLPost(DataBaseAdapter):

    """
    DBURLPost Class - Adapter to call an URL with statistics
    Inheritance diagram for :class:`DBAdapters.DBURLPost`:
    .. inheritance-diagram:: DBAdapters.DBURLPost
    Example:
      >>> dbadapter = DBURLPost(url="http://localhost/post.py", identify="test")
    
    The parameters that will be sent is all the statistics described in the :class:`Statistics.Statistics`
    class, and the parameters:
    
    **generation**
      The generation of the statistics
    
    **identify**
      The id specified by user
    
    .. note:: see the :class:`Statistics.Statistics` documentation.
    
    :param url: the URL to be used
    :param identify: the identify of the run
    :param frequency: the generational dump frequency
    :param post: if True, the POST method will be used, otherwise GET will be used.
    
    """

    def __init__(self, url, identify=None,
                frequency=constants.CDefURLPostStatsGenFreq, post=True, name=constants.CDefURLPostName):

        """
        The creator of the DBURLPost Class.
        """

        super(DBURLPost, self).__init__(frequency, identify, name)
        self.urllibmod = None

        self.url = url
        self.post = post

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        The string representation of adapter
        """

        ret = "DBURLPost DB Adapter [URL='%s', identify='%s']" % (self.url, self.getIdentify())
        return ret

    # -----------------------------------------------------------------

    def open(self, ga_engine):

        """
        Load the modules needed
        :param ga_engine: the GA Engine
        .. versionchanged:: 0.6
         The method now receives the *ga_engine* parameter.
        """

        # Debugging
        log.debug("Opening the " + self.name + " ...")

        if self.urllibmod is None:
         log.debug("Loading urllib module...")
         self.urllibmod = utils.importSpecial("urllib")

    # -----------------------------------------------------------------

    def insert(self, ga_engine):

        """ Sends the data to the URL using POST or GET
        :param ga_engine: the GA Engine
        .. versionchanged:: 0.6
         The method now receives the *ga_engine* parameter.
        """

        log.debug("Sending http request to %s.", self.url)
        stats = ga_engine.getStatistics()
        response = None
        params = stats.internalDict.copy()
        params["generation"] = ga_engine.getCurrentGeneration()
        params["identify"] = self.getIdentify()
        if self.post:  # POST
         response = self.urllibmod.urlopen(self.url, self.urllibmod.urlencode(params))
        else:  # GET
         response = self.urllibmod.urlopen(self.url + "?%s" % (self.urllibmod.urlencode(params)))
        if response:
         response.close()

# -----------------------------------------------------------------

class DBSQLite(DataBaseAdapter):

   """
   DBSQLite Class - Adapter to dump data in SQLite3 database format
   Inheritance diagram for :class:`DBAdapters.DBSQLite`:
   .. inheritance-diagram:: DBAdapters.DBSQLite
   Example:
      >>> dbadapter = DBSQLite(identify="test")

   When you run some GA for the first time, you need to create the database, for this, you
   must use the *resetDB* parameter:

      >>> dbadapter = DBSQLite(identify="test", resetDB=True)

   This parameter will erase all the database tables and will create the new ones.
   The *resetDB* parameter is different from the *resetIdentify* parameter, the *resetIdentify*
   only erases the rows with the same "identify" name.
   """

   def __init__(self, dbname=constants.CDefSQLiteDBName, identify=None, resetDB=False,
                resetIdentify=True, frequency=constants.CDefSQLiteStatsGenFreq,
                commit_freq=constants.CDefSQLiteStatsCommitFreq, name=constants.CDefSQLiteName):

        """
        The creator of the DBSQLite Class
        :param dbname: the database filename
        :param identify: the identify if the run
        :param resetDB: if True, the database structure will be recreated
        :param resetIdentify: if True, the identify with the same name will be overwrite with new data
        :param frequency: the generational dump frequency
        :param commit_freq: the commit frequency
        """

        # Call the constructor of the base class
        super(DBSQLite, self).__init__(frequency, identify, name)

        # Set attributes
        self.sqlite3mod = None
        self.connection = None
        self.resetDB = resetDB
        self.resetIdentify = resetIdentify
        self.dbName = dbname
        self.typeDict = {types.FloatType: "real"}
        self.cursorPool = None
        self.commitFreq = commit_freq

        self.named_individuals = False

   # -----------------------------------------------------------------

   def __repr__(self):

        """
        The string representation of adapter
        """

        ret = "DBSQLite DB Adapter [File='%s', identify='%s']" % (self.dbName, self.getIdentify())
        return ret

   # -----------------------------------------------------------------

   def open(self, ga_engine):

        """
        Open the database connection
        :param ga_engine: the GA Engine
        .. versionchanged:: 0.6
         The method now receives the *ga_engine* parameter.
        """

        # Debugging
        log.debug("Opening the " + self.name + " ...")

        if self.sqlite3mod is None:
         log.debug("Loading sqlite3 module ...")
         self.sqlite3mod = utils.importSpecial("sqlite3")

        log.debug("Opening the " + self.name + " from '%s'", self.dbName)
        self.connection = self.sqlite3mod.connect(self.dbName)

        temp_stats = statistics.Statistics()

        # Set named individuals flag
        self.named_individuals = ga_engine.named_individuals

        # NEW
        self.createStructure(temp_stats)

        # Reset
        if self.resetDB: self.resetStructure(statistics.Statistics())

        # Create structure
        #self.createStructure(temp_stats)

        # Reset
        if self.resetIdentify: self.resetTableIdentify()

   # -----------------------------------------------------------------

   def commit_and_close(self):

        """
        Commit changes on database and closes connection
        """

        self.commit()
        self.close()

   # -----------------------------------------------------------------

   def close(self):

        """
        Close the database connection
        """

        # Debugging
        log.debug("Closing the " + self.name + " ...")

        if self.cursorPool:

         self.cursorPool.close()
         self.cursorPool = None

        self.connection.close()

   # -----------------------------------------------------------------

   def commit(self):

        """
        Commit changes to database
        """

        # Debugging
        log.debug("Commiting changes to the " + self.name + " ...")

        self.connection.commit()

   # -----------------------------------------------------------------

   def getCursor(self):

        """
        Return a cursor from the pool
        :rtype: the cursor
        """

        if not self.cursorPool:

         log.debug("Creating new cursor for the " + self.name + " ...")
         self.cursorPool = self.connection.cursor()
         return self.cursorPool

        else: return self.cursorPool

   # -----------------------------------------------------------------

   def createStructure(self, stats):

        """
        Create table using the Statistics class structure
        :param stats: the statistics object:
        """

        # Inform the user
        log.debug("Creating structure for the " + self.name + " ...")

        c = self.getCursor()
        pstmt = "create table if not exists %s(identify text, generation integer, " % (constants.CDefSQLiteDBTable)
        for k, v in stats.items():
         pstmt += "%s %s, " % (k, self.typeDict[type(v)])
        pstmt = pstmt[:-2] + ")"
        log.debug("Creating table %s: %s.", constants.CDefSQLiteDBTable, pstmt)
        c.execute(pstmt)

        # Create command
        if self.named_individuals: pstmt = """create table if not exists %s (identify text, generation integer, individual text, fitness real, raw real)""" % (constants.CDefSQLiteDBTablePop)
        else: pstmt = """create table if not exists %s (identify text, generation integer, individual integer, fitness real, raw real)""" % (constants.CDefSQLiteDBTablePop)

        log.debug("Creating table %s: %s.", constants.CDefSQLiteDBTablePop, pstmt)

        c.execute(pstmt)
        self.commit()

   # -----------------------------------------------------------------

   def resetTableIdentify(self):

        """
        Delete all records on the table with the same Identify
        """

        # Debugging
        log.debug("Resetting the " + self.name + " for the identifier '" + self.getIdentify() + "' ...")

        c = self.getCursor()
        stmt = "delete from %s where identify = ?" % (constants.CDefSQLiteDBTable)
        stmt2 = "delete from %s where identify = ?" % (constants.CDefSQLiteDBTablePop)

        log.debug("Erasing data from the tables with the run identifier = %s", self.getIdentify())

        try:
            c.execute(stmt, (self.getIdentify(),))
            c.execute(stmt2, (self.getIdentify(),))

        except self.sqlite3mod.OperationalError, expt:
         if str(expt).find("no such table") >= 0:
            print "\n ## The DB Adapter can't find the tables ! Consider enable the parameter resetDB ! ##\n"

        self.commit()

   # -----------------------------------------------------------------

   def resetStructure(self, stats):

        """
        Deletes de current structure and calls createStructure
        :param stats: the statistics object
        """

        # Debugging
        log.debug("Resetting structure of the " + self.name + ", dropping table and creating new empty table ...")

        # Drop
        c = self.getCursor()
        c.execute("drop table if exists %s" % (constants.CDefSQLiteDBTable,))
        c.execute("drop table if exists %s" % (constants.CDefSQLiteDBTablePop,))

        # Commit, create structure
        self.commit()
        self.createStructure(stats)

   # -----------------------------------------------------------------

   def insert(self, ga_engine):

        """
        Inserts the statistics data to database
        :param ga_engine: the GA Engine
        .. versionchanged:: 0.6
         The method now receives the *ga_engine* parameter.
        """

        # Inform the user
        log.debug("Inserting data into the " + self.name + " ...")

        # Get stats, and population and generation information
        stats = ga_engine.getStatistics()
        population = ga_engine.get_population()
        generation = ga_engine.getCurrentGeneration()

        c = self.getCursor()
        pstmt = "insert into %s values (?, ?, " % (constants.CDefSQLiteDBTable)

        for i in xrange(len(stats)): pstmt += "?, "
        pstmt = pstmt[:-2] + ")"
        c.execute(pstmt, (self.getIdentify(), generation) + stats.asTuple())

        pstmt = "insert into %s values(?, ?, ?, ?, ?)" % (constants.CDefSQLiteDBTablePop,)
        tups = []

        # Named
        if self.named_individuals:

         # Loop over the individuals
         for name in population.names:

            ind = population[name]
            tups.append((self.getIdentify(), generation, name, ind.fitness, ind.score))

        # Not named
        else:

         for i in xrange(len(population)):
            ind = population[i]
            tups.append((self.getIdentify(), generation, i, ind.fitness, ind.score))

        # Execute SQL commands
        c.executemany(pstmt, tups)

        # Commit?
        self.commit()
        #if (generation % self.commitFreq == 0): self.commit()

# -----------------------------------------------------------------

class DBXMLRPC(DataBaseAdapter):

    """
    DBXMLRPC Class - Adapter to dump statistics to a XML Remote Procedure Call
    Inheritance diagram for :class:`DBAdapters.DBXMLRPC`:
    .. inheritance-diagram:: DBAdapters.DBXMLRPC
    
    Example:
      >>> adapter = DBXMLRPC(url="http://localhost:8000/", identify="run_01", frequency = 1)
    
      :param url: the URL of the XML RPC
      :param identify: the identify of the run
      :param frequency: the generational dump frequency
    
    
    .. note:: The XML RPC Server must implement the *insert* method, wich receives
             a python dictionary as argument.
    
    Example of an server in Python: ::
    
      import xmlrpclib
      from SimpleXMLRPCServer import SimpleXMLRPCServer
    
      def insert(l):
          print "Received statistics: %s" % l
    
      server = SimpleXMLRPCServer(("localhost", 8000), allow_none=True)
      print "Listening on port 8000..."
      server.register_function(insert, "insert")
      server.serve_forever()
    
    .. versionadded:: 0.6
      The :class:`DBXMLRPC` class.
    """

    def __init__(self, url, identify=None, frequency=constants.CDefXMLRPCStatsGenFreq, name="XML remote procedure"):
    
        """
        The creator of DBXMLRPC Class
        """

        # Call the constructor of the base class
        super(DBXMLRPC, self).__init__(frequency, identify, name)
        self.xmlrpclibmod = None

        self.url = url
        self.proxy = None

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        The string representation of adapter
        """

        ret = "DBXMLRPC DB Adapter [URL='%s', identify='%s']" % (self.url, self.getIdentify())
        return ret

    # -----------------------------------------------------------------

    def open(self, ga_engine):

        """
        Open the XML RPC Server proxy
        :param ga_engine: the GA Engine
        .. versionchanged:: 0.6
         The method now receives the *ga_engine* parameter.
        """

        # Debugging
        log.debug("Opening the " + self.name + " ...")

        if self.xmlrpclibmod is None:
             log.debug("Loding the xmlrpclib module...")
             self.xmlrpclibmod = utils.importSpecial("xmlrpclib")

        log.debug("Opening the XML RPC Server Proxy on %s", self.url)
        self.proxy = self.xmlrpclibmod.ServerProxy(self.url, allow_none=True)

    # -----------------------------------------------------------------

    def insert(self, ga_engine):

        """
        Calls the XML RPC procedure
        :param ga_engine: the GA Engine
        .. versionchanged:: 0.6
         The method now receives the *ga_engine* parameter.
        """

        # Inform the user
        log.debug("Inserting data into the " + self.name + " ...")

        stats = ga_engine.getStatistics()
        generation = ga_engine.getCurrentGeneration()
        di = stats.internalDict.copy()
        di.update({"identify": self.getIdentify(), "generation": generation})
        self.proxy.insert(di)

# -----------------------------------------------------------------

class DBVPythonGraph(DataBaseAdapter):

   """
   The DBVPythonGraph Class - A DB Adapter for real-time visualization using VPython
   Inheritance diagram for :class:`DBAdapters.DBVPythonGraph`:
   .. inheritance-diagram:: DBAdapters.DBVPythonGraph
   .. note:: to use this DB Adapter, you **must** install VPython first.
   Example:
      adapter = DBVPythonGraph(identify="run_01", frequency = 1)
      ga_engine.setDBAdapter(adapter)

   :param identify: the identify of the run
   :param genmax: use the generations as max value for x-axis, default False
   :param frequency: the generational dump frequency

   .. versionadded:: 0.6
      The *DBVPythonGraph* class.
   """

   def __init__(self, identify=None, frequency=20, genmax=False, name="graph"):

        """
        This function ...
        :param identify:
        :param frequency:
        :param genmax:
        """

        # Call the constructor of the base class
        super(DBVPythonGraph, self).__init__(frequency, identify, name)

        # Set properties
        self.genmax = genmax
        self.vtkGraph = None
        self.curveMin = None
        self.curveMax = None
        self.curveDev = None
        self.curveAvg = None

   # -----------------------------------------------------------------

   def makeDisplay(self, title_sec, x, y, ga_engine):

      """
      Used internally to create a new display for VPython.
      :param title_sec: the title of the window
      :param x: the x position of the window
      :param y: the y position of the window
      :param ga_engine: the GA Engine

      :rtype: the window (the return of gdisplay call)
      """
      #title = "Pyevolve v.%s - %s - id [%s]" % (__version__, title_sec, self.identify)
      title = "PTS/evolve"
      if self.genmax:
         disp = self.vtkGraph.gdisplay(title=title, xtitle='Generation', ytitle=title_sec,
                                       xmax=ga_engine.getGenerations(), xmin=0., width=500,
                                       height=250, x=x, y=y)
      else:
         disp = self.vtkGraph.gdisplay(title=title, xtitle='Generation', ytitle=title_sec,
                                       xmin=0., width=500, height=250, x=x, y=y)
         return disp

   # -----------------------------------------------------------------

   def open(self, ga_engine):

      """
      Imports the VPython module and creates the four graph windows
      :param ga_engine: the GA Engine
      """

      log.debug("Opening the " + self.name + " ...")
      #log.debug("Loading visual.graph (VPython) module...")

      if self.vtkGraph is None:
         self.vtkGraph = utils.importSpecial("visual.graph").graph

      display_rawmin = self.makeDisplay("Raw Score (min)", 0, 0, ga_engine)
      display_rawmax = self.makeDisplay("Raw Score (max)", 0, 250, ga_engine)
      display_rawdev = self.makeDisplay("Raw Score (std. dev.)", 500, 0, ga_engine)
      display_rawavg = self.makeDisplay("Raw Score (avg)", 500, 250, ga_engine)

      self.curveMin = self.vtkGraph.gcurve(color=self.vtkGraph.color.red, gdisplay=display_rawmin)
      self.curveMax = self.vtkGraph.gcurve(color=self.vtkGraph.color.green, gdisplay=display_rawmax)
      self.curveDev = self.vtkGraph.gcurve(color=self.vtkGraph.color.blue, gdisplay=display_rawdev)
      self.curveAvg = self.vtkGraph.gcurve(color=self.vtkGraph.color.orange, gdisplay=display_rawavg)

   # -----------------------------------------------------------------

   def insert(self, ga_engine):

        """
        Plot the current statistics to the graphs
        :param ga_engine: the GA Engine
        """

        # Inform the user
        log.debug("Inserting data into the " + self.name + " ...")

        stats = ga_engine.getStatistics()
        generation = ga_engine.getCurrentGeneration()

        self.curveMin.plot(pos=(generation, stats["rawMin"]))
        self.curveMax.plot(pos=(generation, stats["rawMax"]))
        self.curveDev.plot(pos=(generation, stats["rawDev"]))
        self.curveAvg.plot(pos=(generation, stats["rawAve"]))

# -----------------------------------------------------------------

class DBMySQLAdapter(DataBaseAdapter):

   """ DBMySQLAdapter Class - Adapter to dump data in MySql database server
   Inheritance diagram for :class:`DBAdapters.DBMySQLAdapter`:
   .. inheritance-diagram:: DBAdapters.DBMySQLAdapter
   Example:
      >>> dbadapter = DBMySQLAdapter("pyevolve_username", "password", identify="run1")

   or

      >>> dbadapter = DBMySQLAdapter(user="username", passwd="password",
      ...                            host="mysqlserver.com.br", port=3306, db="pyevolve_db")

   When you run some GA for the first time, you need to create the database, for this, you
   must use the *resetDB* parameter as True.

   This parameter will erase all the database tables and will create the new ones.
   The *resetDB* parameter is different from the *resetIdentify* parameter, the *resetIdentify*
   only erases the rows with the same "identify" name, and *resetDB* will drop and recreate
   the tables.

   :param user: mysql username (must have permission to create, drop, insert, etc.. on tables
   :param passwd: the user password on MySQL server
   :param host: the hostname, default is "localhost"
   :param port: the port, default is 3306
   :param db: the database name, default is "pyevolve"
   :param identify: the identify if the run
   :param resetDB: if True, the database structure will be recreated
   :param resetIdentify: if True, the identify with the same name will be overwrite with new data
   :param frequency: the generational dump frequency
   :param commit_freq: the commit frequency
   """

   def __init__(self, user, passwd, host=constants.CDefMySQLDBHost, port=constants.CDefMySQLDBPort,
                db=constants.CDefMySQLDBName, identify=None, resetDB=False, resetIdentify=True,
                frequency=constants.CDefMySQLStatsGenFreq, commit_freq=constants.CDefMySQLStatsCommitFreq,
                name=constants.CDefMySQLName):

      """
      The creator of the DBMySQLAdapter Class
      """

      # Callt he cosntructor of the base class
      super(DBMySQLAdapter, self).__init__(frequency, identify, name)

      # Attributes
      self.mysqldbmod = None
      self.connection = None
      self.resetDB = resetDB
      self.resetIdentify = resetIdentify
      self.db = db
      self.host = host
      self.port = port
      self.user = user
      self.passwd = passwd
      self.typeDict = {types.FloatType: "DOUBLE(14,6)"}
      self.cursorPool = None
      self.commitFreq = commit_freq

      self.named_individuals = False

   # -----------------------------------------------------------------

   def __repr__(self):

      """
      The string representation of adapter
      """

      ret = "DBMySQLAdapter DB Adapter [identify='%s', host='%s', username='%s', db='%s']" % (self.getIdentify(),
            self.host, self.user, self.db)

      return ret

   # -----------------------------------------------------------------

   def open(self, ga_engine):

      """
      Open the database connection
      :param ga_engine: the GA Engine
      .. versionchanged:: 0.6
         The method now receives the *ga_engine* parameter.
      """

      # Debugging
      log.debug("Opening the " + self.name + " ...")

      if self.mysqldbmod is None:
         log.debug("Loading MySQLdb module...")
         self.mysqldbmod = utils.importSpecial("MySQLdb")

      # Set the named_individuals flag
      self.named_individuals = ga_engine.named_individuals

      log.debug("Opening database, host=%s", self.host)
      self.connection = self.mysqldbmod.connect(host=self.host, user=self.user,
                                                passwd=self.passwd, db=self.db,
                                                port=self.port)
      temp_stats = statistics.Statistics()
      self.createStructure(temp_stats)

      if self.resetDB: self.resetStructure(statistics.Statistics())

      if self.resetIdentify: self.resetTableIdentify()

   # -----------------------------------------------------------------

   def commit_and_close(self):

      """
      Commit changes on database and closes connection
      """

      self.commit()
      self.close()

   # -----------------------------------------------------------------

   def close(self):

      """
      Close the database connection
      """

      log.debug("Closing the " + self.name + " ...")

      if self.cursorPool:
         self.cursorPool.close()
         self.cursorPool = None
      self.connection.close()

   # -----------------------------------------------------------------

   def commit(self):

      """
      Commit changes to database
      """

      # Debugging
      log.debug("Commiting changes to the " + self.name + " ...")

      self.connection.commit()

   # -----------------------------------------------------------------

   def getCursor(self):

      """
      Return a cursor from the pool
      :rtype: the cursor
      """

      if not self.cursorPool:

         log.debug("Creating new cursor for database ...")
         self.cursorPool = self.connection.cursor()
         return self.cursorPool

      else: return self.cursorPool

   # -----------------------------------------------------------------

   def createStructure(self, stats):

      """
      Create table using the Statistics class structure
      :param stats: the statistics object
      """

      # INform the user
      log.debug("Creating structure of the " + self.name + " ...")

      c = self.getCursor()
      pstmt = "create table if not exists %s(identify VARCHAR(80), generation INTEGER, " % (constants.CDefMySQLDBTable)
      for k, v in stats.items():
         pstmt += "%s %s, " % (k, self.typeDict[type(v)])
      pstmt = pstmt[:-2] + ")"
      log.debug("Creating table %s: %s", constants.CDefSQLiteDBTable, pstmt)
      c.execute(pstmt)

      if self.named_individuals:
         pstmt = """create table if not exists %s(identify VARCHAR(80), generation INTEGER,
                           individual VARCHAR(80), fitness DOUBLE(14,6), raw DOUBLE(14,6))""" % (constants.CDefMySQLDBTablePop)
      else:
         pstmt = """create table if not exists %s(identify VARCHAR(80), generation INTEGER,
                 individual INTEGER, fitness DOUBLE(14,6), raw DOUBLE(14,6))""" % (constants.CDefMySQLDBTablePop)

      log.debug("Creating table %s: %s", constants.CDefMySQLDBTablePop, pstmt)

      # Execute
      c.execute(pstmt)

      self.commit()

   # -----------------------------------------------------------------

   def resetTableIdentify(self):

      """
      Delete all records on the table with the same Identify
      """

      # Inform the user
      log.debug("Resetting the " + self.name + " for identifier '" + self.getIdentify() + "' ...")

      c = self.getCursor()

      stmt = "delete from %s where identify = '%s'" % (constants.CDefMySQLDBTable, self.getIdentify())
      stmt2 = "delete from %s where identify = '%s'" % (constants.CDefMySQLDBTablePop, self.getIdentify())

      log.debug("Erasing data from the tables with the run identifier = %s", self.getIdentify())
      c.execute(stmt)
      c.execute(stmt2)

      self.commit()

   # -----------------------------------------------------------------

   def resetStructure(self, stats):

      """
      Deletes de current structure and calls createStructure
      :param stats: the statistics object
      """

      # Debugging
      log.debug("Reseting structure of the " + self.name + ", dropping table and creating new empty table ...")

      # Drop
      c = self.getCursor()
      c.execute("drop table if exists %s" % (constants.CDefMySQLDBTable,))
      c.execute("drop table if exists %s" % (constants.CDefMySQLDBTablePop,))

      # Commit, create structure
      self.commit()
      self.createStructure(stats)

   # -----------------------------------------------------------------

   def insert(self, ga_engine):

      """ Inserts the statistics data to database
      :param ga_engine: the GA Engine
      .. versionchanged:: 0.6
         The method now receives the *ga_engine* parameter.
      """

      # Inform the user
      log.debug("Inserting data into the " + self.name + " ...")

      # Get data
      stats = ga_engine.getStatistics()
      population = ga_engine.get_population()
      generation = ga_engine.getCurrentGeneration()

      c = self.getCursor()
      pstmt = "insert into " + constants.CDefMySQLDBTable + " values (%s, %s, "
      for i in xrange(len(stats)):
         pstmt += "%s, "
      pstmt = pstmt[:-2] + ")"
      c.execute(pstmt, (self.getIdentify(), generation) + stats.asTuple())

      pstmt = "insert into " + constants.CDefMySQLDBTablePop + " values(%s, %s, %s, %s, %s)"

      tups = []

      if self.named_individuals:
         for name in population.names:
            ind = population[name]
            tups.append((self.getIdentify(), generation, name, ind.fitness, ind.score))
      else:
         for i in xrange(len(population)):
            ind = population[i]
            tups.append((self.getIdentify(), generation, i, ind.fitness, ind.score))

      c.executemany(pstmt, tups)

      # Commit?
      #self.commit()
      if generation % self.commitFreq == 0: self.commit()

# -----------------------------------------------------------------
