#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.check_database Check the database created by the optimizer.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import sqlite3

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter
from pts.core.tools.logging import setup_log
from pts.core.basics.range import IntegerRange

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition(write_config=False)
definition.add_required("filename", "file_path", "database file name/path")
definition.add_required("run_id", "string", "optimization run ID")
definition.add_required("check", "nonnegative_integer", "")
setter = InteractiveConfigurationSetter("check_database", add_cwd=False, add_logging=False)
config = setter.run(definition)

# Create logger
log = setup_log(level="DEBUG")

# -----------------------------------------------------------------

# Connect to the database
conn = sqlite3.connect(config.filename)

# Creating rows
conn.row_factory = sqlite3.Row

# Create cursor
c = conn.cursor()

# -----------------------------------------------------------------

identifier = config.run_id

# -----------------------------------------------------------------

# Range of generations
generation_range = IntegerRange(0, 1)

# Range of individuals
individual_range = IntegerRange(1, 2)

# -----------------------------------------------------------------

# First check
if config.index == 0:

    # Select multiple generations
    ret = c.execute("select distinct generation from population where identify = ? and generation between ? and ?", (identifier, generation_range.min, generation_range.max))
    generations = ret.fetchall()

    # Select just one generation
    #ret = c.execute("select distinct generation from population where identify = ?", (identifier,))
    #generations = ret.fetchall()

    if len(generations) == 0: raise RuntimeError("No generations found in the range")

    # The data
    pop = []

    # Loop over the generations
    for gen in generations:

        print("generation", gen)

        pop_tmp = []

        #if options.lindrange:

        ret = c.execute("""
                   select *  from population
                   where identify = ?
                   and generation = ?
                   and individual between ? and ?
                   """, (identifier, gen[0], individual_range.min, individual_range.max))

        #else:
        #    ret = c.execute("""
        #               select *  from population
        #               where identify = ?
        #               and generation = ?
        #               """, (options.identify, gen[0]))

        ret_fetch = ret.fetchall()

        if len(ret_fetch) == 0: raise RuntimeError("No individuals found in the range")

        # Loop over the individuals
        for it in ret_fetch:

            print("individual", it)

            # Get scores
            raw = it["raw"]
            fitness = it["fitness"]

            print("raw", raw)
            print("fitness", fitness)

            # Add the score
            pop_tmp.append(fitness)

        pop.append(pop_tmp)

    # Close
    ret.close()

    # Show the data
    print(pop)

# -----------------------------------------------------------------

elif config.index == 1:

    #if options.genrange:
        #genrange = options.genrange.split(":")
    ret = c.execute("select * from statistics where identify = ? and generation between ? and ?", (identifier, generation_range.min, generation_range.max))
    #else: ret = c.execute("select * from statistics where identify = ?", (options.identify,))

    pop = ret.fetchall()

    ret.close()

    print("Generations:", pop)

# -----------------------------------------------------------------

elif config.index == 2: pass

# -----------------------------------------------------------------

conn.close()

# -----------------------------------------------------------------
