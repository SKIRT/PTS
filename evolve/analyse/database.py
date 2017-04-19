#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.analyse.database Contains functions for analysing the database

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import sqlite3

# -----------------------------------------------------------------

def load_database(path):

    """
    This function ...
    :param path: 
    :return: 
    """

    # Connect to the database
    conn = sqlite3.connect(path)

    # Creating rows
    conn.row_factory = sqlite3.Row

    # Create cursor
    c = conn.cursor()

    # Return the cursor
    return c

# -----------------------------------------------------------------

def get_generations(path, run_id):

    """
    This function ...
    :param path: 
    :param run_id: 
    :return: 
    """

    cursor = load_database(path)

    # Select multiple generations
    ret = cursor.execute("select distinct generation from population where identify = ?", run_id)
    generations = ret.fetchall()

    # Return the generation numbers
    return [gen[0] for gen in generations]

# -----------------------------------------------------------------

def _get_individuals(path, run_id, generation, individual_range=None):

    """
    This function ...
    :param path: 
    :param run_id: 
    :param generation: 
    :return: 
    """

    # Get the cursor
    cursor = load_database(path)

    if individual_range is not None:
        ret = cursor.execute("""
                             select *  from population
                             where identify = ?
                             and generation = ?
                             and individual between ? and ?
                             """, (run_id, generation, individual_range.min, individual_range.max))
    else:
        ret = cursor.execute("""
                           select *  from population
                           where identify = ?
                           and generation = ?
                           """, (run_id, generation))

    ret_fetch = ret.fetchall()

    if len(ret_fetch) == 0: raise RuntimeError("No individuals found in the range")

    return ret_fetch

# -----------------------------------------------------------------

def get_scores(path, run_id, generation):

    """
    This function ...
    :param path: 
    :param run_id: 
    :param generation: 
    :return: 
    """

    scores = []

    # Loop over the individuals
    for it in _get_individuals(path, run_id, generation):

        # Get score
        raw = it["raw"]

        # Add the score
        scores.append(raw)

    # Return the scores
    return scores

# -----------------------------------------------------------------

def get_fitnesses(path, run_id, generation):

    """
    This function ...
    :param path: 
    :param run_id: 
    :param generation: 
    :return: 
    """

    fitnesses = []

    # Loop over the individuals
    for it in _get_individuals(path, run_id, generation):

        # Get score
        fitness = it["fitness"]

        # Add the score
        fitnesses.append(fitness)

    # Return the scores
    return fitnesses

# -----------------------------------------------------------------
