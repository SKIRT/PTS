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
import Scaling
import Selectors
import Initializators
import Mutators
import Crossovers
from GTree import GTreeGP

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
CDefPopScale = Scaling.LinearScaling

# - GA Engine defaults
CDefGAGenerations = 100
CDefGAMutationRate = 0.02
CDefGACrossoverRate = 0.9
CDefGAPopulationSize = 80
CDefGASelector = Selectors.GRankSelector
CDefGAElitismReplacement = 1

# - This is general used by integer/real ranges defaults
CDefRangeMin = 0
CDefRangeMax = 100

# - G1DBinaryString defaults
CDefG1DBinaryStringMutator = Mutators.G1DBinaryStringMutatorFlip
CDefG1DBinaryStringCrossover = Crossovers.G1DBinaryStringXSinglePoint
CDefG1DBinaryStringInit = Initializators.G1DBinaryStringInitializator
CDefG1DBinaryStringUniformProb = 0.5

# - G2DBinaryString defaults
CDefG2DBinaryStringMutator = Mutators.G2DBinaryStringMutatorFlip
CDefG2DBinaryStringCrossover = Crossovers.G2DBinaryStringXUniform
CDefG2DBinaryStringInit = Initializators.G2DBinaryStringInitializator
CDefG2DBinaryStringUniformProb = 0.5

# - GTree defaults
CDefGTreeInit = Initializators.GTreeInitializatorInteger
CDefGGTreeMutator = Mutators.GTreeMutatorIntegerRange
CDefGTreeCrossover = Crossovers.GTreeCrossoverSinglePointStrict

# - GTreeGP defaults
CDefGTreeGPInit = Initializators.GTreeGPInitializator
CDefGGTreeGPMutator = Mutators.GTreeGPMutatorSubtree
CDefGTreeGPCrossover = Crossovers.GTreeGPCrossoverSinglePoint

# - G1DList defaults
CDefG1DListMutIntMU = 2
CDefG1DListMutIntSIGMA = 10

CDefG1DListMutRealMU = 0
CDefG1DListMutRealSIGMA = 1

CDefG1DListMutator = Mutators.G1DListMutatorSwap
CDefG1DListCrossover = Crossovers.G1DListCrossoverSinglePoint
CDefG1DListInit = Initializators.G1DListInitializatorInteger
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

CDefG2DListMutator = Mutators.G2DListMutatorSwap
CDefG2DListCrossover = Crossovers.G2DListCrossoverUniform
CDefG2DListInit = Initializators.G2DListInitializatorInteger
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
