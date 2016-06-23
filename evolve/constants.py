#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.Consts Pyevolve have defaults in all genetic operators, settings and etc,
#  this is an issue to help the user in the API use and minimize the source code needed to make simple things.
#  In the module :mod:`Consts`, you will find those defaults settings. You are encouraged to see the constants,
#  but not to change directly on the module, there are methods for this.
#

# -----------------------------------------------------------------

# Import other evolve modules
import scaling
import selectors
import initializators
import mutators
import crossovers
from tree import GTreeGP

# -----------------------------------------------------------------

# Required python version 2.5+
CDefPythonRequire = (2, 5)

# Types of sort
# - raw: uses the "score" attribute
# - scaled: uses the "fitness" attribute
sortType = {
    "raw": 0,
    "scaled": 1
}

# Optimization type
# - Minimize or Maximize the Evaluator Function
minimaxType = {"minimize": 0,
               "maximize": 1
               }

CDefESCKey = 27

CDefImportList = {"visual.graph": "you must install VPython !",
                  "csv": "csv module not found !",
                  "urllib": "urllib module not found !",
                  "sqlite3": "sqlite3 module not found, are you using Jython or IronPython ?",
                  "xmlrpclib": "xmlrpclib module not found !",
                  "MySQLdb": "MySQLdb module not found, you must install mysql-python !",
                  "pydot": "Pydot module not found, you must install Pydot to plot graphs !"}

####################
# Defaults section #
####################

# - Tournament selector
CDefTournamentPoolSize = 2

# - Scale methods defaults
CDefScaleLinearMultiplier = 1.2
CDefScaleSigmaTruncMultiplier = 2.0
CDefScalePowerLawFactor = 1.0005
CDefScaleBoltzMinTemp = 1.0
CDefScaleBoltzFactor = 0.05
# 40 temp. = 500 generations
CDefScaleBoltzStart = 40.0

# - Population Defaults
CDefPopSortType = sortType["scaled"]
CDefPopMinimax = minimaxType["maximize"]
CDefPopScale = scaling.LinearScaling

# - GA Engine defaults
CDefGAGenerations = 100
CDefGAMutationRate = 0.02
CDefGACrossoverRate = 0.9
CDefGAPopulationSize = 80
CDefGASelector = selectors.GRankSelector
CDefGAElitismReplacement = 1

# - This is general used by integer/real ranges defaults
CDefRangeMin = 0
CDefRangeMax = 100

# - G1DBinaryString defaults
CDefG1DBinaryStringMutator = mutators.G1DBinaryStringMutatorFlip
CDefG1DBinaryStringCrossover = crossovers.G1DBinaryStringXSinglePoint
CDefG1DBinaryStringInit = initializators.G1DBinaryStringInitializator
CDefG1DBinaryStringUniformProb = 0.5

# - G2DBinaryString defaults
CDefG2DBinaryStringMutator = mutators.G2DBinaryStringMutatorFlip
CDefG2DBinaryStringCrossover = crossovers.G2DBinaryStringXUniform
CDefG2DBinaryStringInit = initializators.G2DBinaryStringInitializator
CDefG2DBinaryStringUniformProb = 0.5

# - GTree defaults
CDefGTreeInit = initializators.GTreeInitializatorInteger
CDefGGTreeMutator = mutators.GTreeMutatorIntegerRange
CDefGTreeCrossover = crossovers.GTreeCrossoverSinglePointStrict

# - GTreeGP defaults
CDefGTreeGPInit = initializators.GTreeGPInitializator
CDefGGTreeGPMutator = mutators.GTreeGPMutatorSubtree
CDefGTreeGPCrossover = crossovers.GTreeGPCrossoverSinglePoint

# - G1DList defaults
CDefG1DListMutIntMU = 2
CDefG1DListMutIntSIGMA = 10

CDefG1DListMutRealMU = 0
CDefG1DListMutRealSIGMA = 1

CDefG1DListMutator = mutators.G1DListMutatorSwap
CDefG1DListCrossover = crossovers.G1DListCrossoverSinglePoint
CDefG1DListInit = initializators.G1DListInitializatorInteger
CDefG1DListCrossUniformProb = 0.5

# SBX Crossover defaults
# Crossover distribution index for SBX
# Larger Etac = similar to parents
# Smaller Etac = far away from parents
CDefG1DListSBXEtac = 10
CDefG1DListSBXEPS = 1.0e-14

# - G2DList defaults
CDefG2DListMutIntMU = 2
CDefG2DListMutIntSIGMA = 10

CDefG2DListMutRealMU = 0
CDefG2DListMutRealSIGMA = 1

CDefG2DListMutator = mutators.G2DListMutatorSwap
CDefG2DListCrossover = crossovers.G2DListCrossoverUniform
CDefG2DListInit = initializators.G2DListInitializatorInteger
CDefG2DListCrossUniformProb = 0.5

# Gaussian Gradient
CDefGaussianGradientMU = 1.0
CDefGaussianGradientSIGMA = (1.0 / 3.0)  # approx. +/- 3-sigma is +/- 10%

# - DB Adapters SQLite defaults
CDefSQLiteDBName = "pyevolve.db"
CDefSQLiteDBTable = "statistics"
CDefSQLiteDBTablePop = "population"
CDefSQLiteStatsGenFreq = 1
CDefSQLiteStatsCommitFreq = 300

# - DB Adapters MySQL defaults
CDefMySQLDBName = "pyevolve"
CDefMySQLDBTable = "statistics"
CDefMySQLDBTablePop = "population"
CDefMySQLDBHost = "localhost"
CDefMySQLDBPort = 3306
CDefMySQLStatsGenFreq = 1
CDefMySQLStatsCommitFreq = 300

# - DB Adapters URL Post defaults
CDefURLPostStatsGenFreq = 100

# - DB Adapters CSV File defaults
CDefCSVFileName = "pyevolve.csv"
CDefCSVFileStatsGenFreq = 1

# - DB Adapter XML RPC
CDefXMLRPCStatsGenFreq = 20

# Util Consts
CDefBroadcastAddress = "255.255.255.255"
nodeType = {"TERMINAL": 0, "NONTERMINAL": 1}

CDefGPGenomes = [GTreeGP]

# Migration Consts
CDefGenMigrationRate = 20
CDefMigrationNIndividuals = 3
CDefGenMigrationReplacement = 3

CDefNetworkIndividual = 1
CDefNetworkInfo = 2

# -----------------------------------------------------------------
