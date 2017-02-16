#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.plot.heatmap Contains the HeatMapPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import matplotlib.pyplot as plt
import sqlite3

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.plot import Plot
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class HeatMapPlotter(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(HeatMapPlotter, self).__init__(config)

        # The database connection
        self.database = None

        # The plot
        self.plot = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        if len(identify_list) == 1 and not popGraph:

            if options.compare_raw or options.compare_fitness:
                parser.error("You can't use this graph type with only one identify !")

            conn = sqlite3.connect(options.dbfile)
            conn.row_factory = sqlite3.Row
            c = conn.cursor()

            if options.genrange:
                genrange = options.genrange.split(":")
                ret = c.execute("select * from statistics where identify = ? and generation between ? and ?",
                                (options.identify, genrange[0], genrange[1]))
            else:
                ret = c.execute("select * from statistics where identify = ?", (options.identify,))

            pop = ret.fetchall()

            ret.close()
            conn.close()

            if len(pop) <= 0:
                print("No statistic data found for the identify '%s' !" % (options.identify,))
                exit()

            print("%d generations found !" % (len(pop),))


        elif len(identify_list) > 1 and not popGraph:

            pop = []
            if (not options.compare_raw) and (not options.compare_fitness):
                parser.error("You can't use many ids with this graph type !")

            conn = sqlite3.connect(options.dbfile)
            conn.row_factory = sqlite3.Row
            c = conn.cursor()
            for item in identify_list:
                if options.genrange:
                    genrange = options.genrange.split(":")
                    ret = c.execute("select * from statistics where identify = ? and generation between ? and ?",
                                    (item, genrange[0], genrange[1]))
                else:
                    ret = c.execute("select * from statistics where identify = ?", (item,))
                fetchall = ret.fetchall()
                if len(fetchall) > 0:
                    pop.append(fetchall)

            ret.close()
            conn.close()

            if len(pop) <= 0:
                print("No statistic data found for the identify list '%s' !" % (options.identify,))
                exit()

            print("%d identify found !" % (len(pop),))


        if options.pop_heatmap_raw:

            if options.outfile:
                graph_pop_heatmap_raw(pop, options.minimize, options.colormap,
                                      options.outfile + "." + options.extension)
            else:
                graph_pop_heatmap_raw(pop, options.minimize, options.colormap)

        if options.pop_heatmap_fitness:

            if options.outfile:
                graph_pop_heatmap_fitness(pop, options.minimize, options.colormap,
                                          options.outfile + "." + options.extension)
            else:
                graph_pop_heatmap_fitness(pop, options.minimize, options.colormap)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(HeatMapPlotter, self).setup(**kwargs)

        # Get the database
        if "database" in kwargs:

            self.database = kwargs.pop("database")

        elif self.config.database is not None:

            self.database = sqlite3.connect(self.config.database)

        elif fs.is_file(fs.join(self.config.path, "database.db")):

            # Load database from file
            self.database = sqlite3.connect(fs.join(self.config.path, "database.db"))

        else: raise ValueError("Database not found and not specified")

    # -----------------------------------------------------------------

    def initialize(self):

        """
        This function ...
        :return:
        """

        print("Loading database and creating graph...")

        identify_list = options.identify.split(",")
        identify_list = map(str.strip, identify_list)

        pop = None

        conn = sqlite3.connect(options.dbfile)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        if options.genrange:
            genrange = options.genrange.split(":")
            ret = c.execute(
                "select distinct generation from population where identify = ? and generation between ? and ?",
                (options.identify, genrange[0], genrange[1]))
        else:
            ret = c.execute("select distinct generation from population where identify = ?", (options.identify,))

        generations = ret.fetchall()
        if len(generations) <= 0:
            print("No generation data found for the identify '%s' !" % (options.identify,))
            exit()

        pop = []
        for gen in generations:
            pop_tmp = []

            if options.lindrange:
                individual_range = options.lindrange.split(":")
                ret = c.execute("""
                                 select *  from population
                                 where identify = ?
                                 and generation = ?
                                 and individual between ? and ?
                                 """, (options.identify, gen[0], individual_range[0], individual_range[1]))
            else:
                ret = c.execute("""
                                 select *  from population
                                 where identify = ?
                                 and generation = ?
                                 """, (options.identify, gen[0]))

            ret_fetch = ret.fetchall()
            for it in ret_fetch:
                if options.pop_heatmap_raw:
                    pop_tmp.append(it["raw"])
                else:
                    pop_tmp.append(it["fitness"])
            pop.append(pop_tmp)

        ret.close()
        conn.close()

        if len(pop) <= 0:
            print("No statistic data found for the identify '%s' !" % (options.identify,))
            exit()

        print("%d generations found !" % (len(pop),))

        popGraph = True

    # -----------------------------------------------------------------

    def graph_pop_heatmap_raw(self, pop, minimize, colormap="jet", filesave=None):

        """
        This function ...
        :param pop:
        :param minimize:
        :param colormap:
        :param filesave:
        :return:
        """

        pylab.imshow(pop, aspect="auto", interpolation="gaussian", cmap=matplotlib.cm.__dict__[colormap])
        pylab.title("Plot of pop. raw scores along the generations")
        pylab.xlabel('Population')
        pylab.ylabel('Generations')
        pylab.grid(True)
        pylab.colorbar()

        if filesave:
            pylab.savefig(filesave)
            print("Graph saved to %s file !" % (filesave,))
        else:
            pylab.show()

    # -----------------------------------------------------------------

    def graph_pop_heatmap_fitness(self, pop, minimize, colormap="jet", filesave=None):

        """
        This function ...
        :param pop:
        :param minimize:
        :param colormap:
        :param filesave:
        :return:
        """

        pylab.imshow(pop, aspect="equal", interpolation="gaussian", cmap=matplotlib.cm.__dict__[colormap])
        pylab.title("Plot of pop. fitness scores along the generations")
        pylab.xlabel('Population')
        pylab.ylabel('Generations')
        pylab.grid(True)
        pylab.colorbar()

        if filesave:
            pylab.savefig(filesave)
            print("Graph saved to %s file !" % (filesave,))
        else:
            pylab.show()

# -----------------------------------------------------------------
