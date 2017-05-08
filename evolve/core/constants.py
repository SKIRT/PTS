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
#

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
#minimaxType = {"minimize": 0,
#               "maximize": 1
#               }

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
CDefPopMinimax = "maximize" #minimaxType["maximize"]
from .scaling import LinearScaling
CDefPopScale = LinearScaling

# - GA Engine defaults
CDefGAGenerations = 100
CDefGAMutationRate = 0.02
CDefGACrossoverRate = 0.9
CDefGAPopulationSize = 80
from .selectors import GRankSelector
CDefGASelector = GRankSelector
CDefGAElitismReplacement = 1

# - This is general used by integer/real ranges defaults
CDefRangeMin = 0
CDefRangeMax = 100

# - G1DBinaryString defaults
from .mutators import G1DBinaryStringMutatorFlip
CDefG1DBinaryStringMutator = G1DBinaryStringMutatorFlip
from .crossovers import G1DBinaryStringXSinglePoint
CDefG1DBinaryStringCrossover = G1DBinaryStringXSinglePoint
from .initializators import G1DBinaryStringInitializator
CDefG1DBinaryStringInit = G1DBinaryStringInitializator
CDefG1DBinaryStringUniformProb = 0.5

# - G2DBinaryString defaults
from .mutators import G2DBinaryStringMutatorFlip
CDefG2DBinaryStringMutator = G2DBinaryStringMutatorFlip
from .crossovers import G2DBinaryStringXUniform
CDefG2DBinaryStringCrossover = G2DBinaryStringXUniform
from .initializators import G2DBinaryStringInitializator
CDefG2DBinaryStringInit = G2DBinaryStringInitializator
CDefG2DBinaryStringUniformProb = 0.5

# - GTree defaults
from .initializators import GTreeInitializatorInteger
CDefGTreeInit = GTreeInitializatorInteger
from .mutators import GTreeMutatorIntegerRange
CDefGGTreeMutator = GTreeMutatorIntegerRange
from .crossovers import GTreeCrossoverSinglePointStrict
CDefGTreeCrossover = GTreeCrossoverSinglePointStrict

# - GTreeGP defaults
from .initializators import GTreeGPInitializator
CDefGTreeGPInit = GTreeGPInitializator
from .mutators import GTreeGPMutatorSubtree
CDefGGTreeGPMutator = GTreeGPMutatorSubtree
from .crossovers import GTreeGPCrossoverSinglePoint
CDefGTreeGPCrossover = GTreeGPCrossoverSinglePoint

# - G1DList defaults
CDefG1DListMutIntMU = 2
CDefG1DListMutIntSIGMA = 10

CDefG1DListMutRealMU = 0
CDefG1DListMutRealSIGMA = 1

from .mutators import G1DListMutatorSwap
CDefG1DListMutator = G1DListMutatorSwap
from .crossovers import G1DListCrossoverSinglePoint
CDefG1DListCrossover = G1DListCrossoverSinglePoint
from .initializators import G1DListInitializatorInteger
CDefG1DListInit = G1DListInitializatorInteger
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

from .mutators import G2DListMutatorSwap
CDefG2DListMutator = G2DListMutatorSwap
from .crossovers import G2DListCrossoverUniform
CDefG2DListCrossover = G2DListCrossoverUniform
from .initializators import G2DListInitializatorInteger
CDefG2DListInit = G2DListInitializatorInteger
CDefG2DListCrossUniformProb = 0.5

# Gaussian Gradient
CDefGaussianGradientMU = 1.0
CDefGaussianGradientSIGMA = (1.0 / 3.0)  # approx. +/- 3-sigma is +/- 10%

# - DB Adapters SQLite defaults
CDefSQLiteName = "SQLite database"
CDefSQLiteDBName = "pyevolve.db"
CDefSQLiteDBTable = "statistics"
CDefSQLiteDBTablePop = "population"
CDefSQLiteStatsGenFreq = 1
CDefSQLiteStatsCommitFreq = 300

# - DB Adapters MySQL defaults
CDefMySQLName = "MySQL database"
CDefMySQLDBName = "pyevolve"
CDefMySQLDBTable = "statistics"
CDefMySQLDBTablePop = "population"
CDefMySQLDBHost = "localhost"
CDefMySQLDBPort = 3306
CDefMySQLStatsGenFreq = 1
#CDefMySQLStatsCommitFreq = 300
CDefMySQLStatsCommitFreq = 1

# - DB Adapters URL Post defaults
CDefURLPostName = "URL"
CDefURLPostStatsGenFreq = 100

# - NEW: DB Adapters for populations file
CDefPopulationsName = "populations file"
CDefPopulationsFileName = "populations.dat"
CDefPopulationsStatsGenFreq = 1

# - DB Adapters CSV File defaults
CDefCSVName = "statistics file"
CDefCSVFileName = "pyevolve.csv"
CDefCSVFileStatsGenFreq = 1

# - DB Adapter XML RPC
CDefXMLRPCStatsGenFreq = 20

# Util Consts
CDefBroadcastAddress = "255.255.255.255"
nodeType = {"TERMINAL": 0, "NONTERMINAL": 1}

from .tree import GTreeGP
CDefGPGenomes = [GTreeGP]

# Migration Consts
CDefGenMigrationRate = 20
CDefMigrationNIndividuals = 3
CDefGenMigrationReplacement = 3

CDefNetworkIndividual = 1
CDefNetworkInfo = 2

# -----------------------------------------------------------------
