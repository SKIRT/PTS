#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.scalingtest A class representing a SKIRT scaling test
#

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path

# Import the relevant PTS class
from pts.skifile import SkiFile
from pts.skirtsimulation import SkirtSimulation
from pts.skirtexec import SkirtExec
from pts.log import Log
from pts.plotscaling import plottimes, plotspeedups, ploteffs

class ScalingTest:

    # -----------------------------------------------------------------

    ## The constructor accepts the following arguments:
    #  - Blabla
    #
    def __init__(self, path, system, mode):
        
        # Set the path
        self._path = path
        
        # Set the output and result paths
        self._outpath = os.path.join(self._path, "out")
        self._respath = os.path.join(self._path, "res")
        
        # Create the output directories if they do not already exist
        try: os.mkdir(self._outpath)
        except OSError: pass
        try: os.mkdir(self._respath)
        except OSError: pass
        
        # Set the system name
        self._system = system
        
        # Set the mode
        self._mode = mode
        
        # Create the logging mechanism
        self._log = Log()
        
        # Create skirt execution context
        self._skirt = SkirtExec(log=self._log)
    
    def run(self, processes, threads):

        # The ski file to be used for the scaling test
        skifile = os.path.join(self._path, "scaling.ski")
        
        # Log the system name, the test mode and the version of SKIRT used for this test
        self._log.info("Starting parallel scaling benchmark for " + self._system + " in " + self._mode + " mode.")        
        self._log.info("Using " + self._skirt.version())

        # Create a file containing the results of the scaling test
        resultfilename = os.path.join(self._respath, self._system + "_" + self._mode + "_" + str(processes) + "_" + str(threads) + ".dat")
        resultfile = open(resultfilename, "w")
        
        # Write a header containing useful information about this test to the results file
        resultfile.write("# Parallel scaling benchmark results for " + self._system + " in " + self._mode + " mode\n")
        resultfile.write("# Using " + self._skirt.version() + "\n")
        resultfile.write("# Column 1: Number of processes p\n")
        resultfile.write("# Column 2: Number of threads per process t\n")
        resultfile.write("# Column 3: Total number of threads (t*p)\n")
        resultfile.write("# Column 4: Execution time for the setup (s)\n")
        resultfile.write("# Column 5: Execution time for the stellar emission phase (s)\n")
        resultfile.write("# Column 6: Execution time for the writing phase (s)\n")
        resultfile.write("# Column 7: Execution time for the simulation (s)\n")

        # Make a list of the different configurations (number of processes, number of threads) to be used for this test
        conflist = self._makelist(processes, threads)
                
        # Perform the simulations
        for configuration in conflist:
            
            self._log.info("Running simulation with " + str(configuration[0]) + " process(es) consisting of " + str(configuration[1]) + " thread(s)")
            self._runandlog(skifile, configuration[0], configuration[1], resultfile)

        # Close the results file
        resultfile.close()
        
        # End with some log messages
        self._log.info("Finished parallel scaling benchmark.")
        self._log.info("Results were written to " + resultfilename)
        
    def plot(self):
        
        # Get a list of result files (i.e. all *.dat files in the current directory)
        filenames = sorted(filter(lambda fn: fn.endswith(".dat"), os.listdir(self._respath)), key=self._totalthreads)

        # Generate the plots
        plottimes(filenames, "scaling_times.pdf", self._respath, xlim=(0,40))
        plotspeedups(filenames, "scaling_speedups.pdf", self._respath, xlim=(0,40))
        ploteffs(filenames, "scaling_effs.pdf", self._respath, xlim=(0,40))

    ## This function runs the simulation once with the specified number of threads,
    # and writes the timing results to the specified file object
    def _runandlog(self, skifile, processes, threads, resultfile):
    
        simulation = self._skirt.execute(skipattern=skifile, outpath=self._outpath, threads=threads, processes=processes)[0]
        if simulation.status() != "Finished": raise ValueError("Simulation " + simulation.status())

        for line in open(simulation.logfilepath()):
            
            if "Finished setup in" in line:
                setuptime = float(line.split(" in ")[1].split()[0])
                
            elif "Finished the stellar emission phase in" in line:
                stellartime = float(line.split(" in ")[1].split()[0])
                    
            elif "Finished writing results in" in line:
                writingtime = float(line.split(" in ")[1].split()[0])

            elif "Finished simulation" in line:
                simulationtime = float(line.split(" in ")[1].split()[0])

        resultfile.write(str(processes) + " " + str(threads) + " " + str(processes*threads) + " " + str(setuptime) 
                                        + " " + str(stellartime) + " " + str(writingtime) + " " + str(simulationtime) + "\n")
        
    ## This function returns the total number of cores for a benchmark result file name,
    # which is assumed be formatted as <system-name>_<#hard-cores>_<#extra-cores>.dat
    def _totalthreads(self, filename):
        
        segments = filename.split("_")
        threads = segments[2].split(".")[0]
        return int(threads)
        
    def _makelist(self, processes, threads):

        list = []
        
        for nprocs in (1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 36, 40, 48, 56):
            
            # Stop when nprocs exceeds the maximum number of processes
            if nprocs > processes: break
            
            # For one processes, determine the thread counts if we are in multithreading or hybrid mode
            if nprocs == 1 and self._mode != "mpi":
                
                for nthreads in (1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 36, 40, 48, 56):
                    
                    # Stop when nthreads exceeds the maximum number of threads
                    if nthreads > threads: break
                    
                    # Add this configuration 3 times
                    list.append((nprocs,nthreads))
                    list.append((nprocs,nthreads))
                    list.append((nprocs,nthreads))
            
            else:
                
                # Add this configuration (nprocs,threads) 3 times
                list.append((nprocs,threads))
                list.append((nprocs,threads))
                list.append((nprocs,threads))
        
        return list