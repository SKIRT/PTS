#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.sphanalyze Analyzing SPH output data.
#
# SKIRT can import particle lists output from SPH simulations to build
# star or dust density distributions. To help optimize this process we
# need to understand how the particles are distributed over the domain,
# taking into account their smoothing length (or rather the "support").
#
# For example, say we create some cuboidal grid over the domain and
# make a list of all particles affecting each cell in the grid. If the
# number of particles affecting each cell is substantially lower than
# the total number of particles, we can accordingly speed up the
# density calculations.
#
# The functions in this module prepares some statistic on SPH data files
# to assist in building or fine-tuning the appropriate algorithms.

# -----------------------------------------------------------------

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------

## Create a one-dimensional grid such that each cell holds approximately
#  the same number of particles.
#
# parameters:
# - coords: unsorted array of particle coordinates to be fitted in the grid
# - cmin,cmax: lower and upper limits of the grid
# - gridsize: number of cells in the grid
#
# returns: array of length gridsize+1 containing grid separation points
#
def makegrid(coords, cmin, cmax, gridsize):
    # determine the particle distribution at a decent resolution
    nbins = gridsize*100
    binwidth = (cmax-cmin)/nbins
    bins = np.zeros(nbins, 'int')
    for coord in coords:
        bins[int( (coord-cmin)/binwidth )] += 1

    # determine grid separation points based on the cumulative distribution
    grid = np.zeros(gridsize+1)
    grid[0] = cmin
    percell = len(coords) / gridsize  # target number of particles per cell
    cumul = 0       # cumulative number of particles in processed bins
    gridindex = 1   # index of the next grid separation point to be filled
    for binindex in range(nbins):
        cumul += bins[binindex]
        if cumul > percell*gridindex:
            grid[gridindex] = cmin + (binindex+1)*binwidth
            gridindex += 1
            if gridindex >= gridsize:
                break
    grid[gridsize] = cmax
    return grid

# -----------------------------------------------------------------

## Determine the range of indices in a one-dimensional grid
#  affected by a certain particle.
#
# parameters:
# - grid: the separation points of the grid; length gridsize+1
# - center: position of the center of the particle
# - radius: radius of the particle
#
# returns: a tuple (i1,i2) with 0 <= i1 < i2 <= gridsize;
#     going from lower to higher indices in the grid:
#        - i1 is the index of the first cell overlapping the particle
#        - i2 is the index of the first subsequent cell no longer overlapping
#             the particle, thus following Python slice semantics
#
def getindexrange(grid, center, radius):
    # process arguments
    gridsize = len(grid)-1
    c1 = center - radius
    c2 = center + radius

    # set defaults so that outliers are treated reasonably
    i1 = i2 = gridsize

    # look for the first affected cell
    for index in range(1, gridsize+1):
        if c1 < grid[index]:
            i1 = index-1
            break

    # look for the first subsequent unaffected cell
    for index in range(i1+1, gridsize+1):
        if c2 < grid[index]:
            i2 = index
            break

    return (i1,i2)

# -----------------------------------------------------------------

## Determine whether an axis-aligned bounding box intersects with a sphere.
#  Algorithm due to Jim Arvo in "Graphics Gems" (1990).
#
# parameters:
# - xmin, xmax, ymin, ymax, zmin, zmax: ll and ur corner of the bounding box
# - xc, yc, zc, r: center and radius of the sphere
#
# returns true if the bounding box and the sphere intersect, false otherwise
#
def intersects(xmin, xmax, ymin, ymax, zmin, zmax,  xc, yc, zc, r):
    squaredist = r*r

    if xc < xmin:
        squaredist -= (xc - xmin)**2
    elif xc > xmax:
        squaredist -= (xc - xmax)**2
    if yc < ymin:
        squaredist -= (yc - ymin)**2
    elif yc > ymax:
        squaredist -= (yc - ymax)**2
    if zc < zmin:
        squaredist -= (zc - zmin)**2
    elif zc > zmax:
        squaredist -= (zc - zmax)**2

    return squaredist > 0.;

# -----------------------------------------------------------------

## Analyze the specified data file
#
# parameters:
# - filename: name of the SPH data file
# - columns: tuple of 0-based column indices for x,y,z,h (smoothing length)
# - scale: multiplication factor that scales all data points
# - supfac: multiplication factor to obtain the support radius from h
#           (i.e. maximum effective range, assuming spherical distribution)
# - gridsize: number of cells in each direction (i.e. gridsize**3 in total)
# - plotfile: name of the plot output file; if missing plot is shown in window
# - basicplot: if true, plot histograms for the particle radii and positions
# - circleplot: if true, plot a circle for each particle position and radius
# - cellplot: if true, create a smart grid and plot a histogram of the
#             number of particles affecting each grid cell
# - accurate: if false, a particle is considered to affect all cells that
#             intersect the particle's enclosing rectangle;
#             if true, only cells that effectively intersect the particle's
#             sphere are counted - but this is A LOT SLOWER
#
def analyze(filename, columns=(0,1,2,3), scale=1, supfac=1, gridsize=9, plotfile=None,
            basicplot=False, circleplot=False, cellplot=False, accurate=False):
    # read the data
    print "---------------------------------------"
    print "Reading data file", filename
    data = np.loadtxt(filename, usecols=columns) * scale
    x,y,z,r = data[:,0], data[:,1], data[:,2], data[:,3]*supfac

    # print number of particles
    print "Number of particles: {0}".format(len(r))

    # print domain size
    xmin, xmax = min(x-r), max(x+r)
    ymin, ymax = min(y-r), max(y+r)
    zmin, zmax = min(z-r), max(z+r)
    print "Domain range: [{0:.0f}..{1:.0f}] [{2:.0f}..{3:.0f}] [{4:.0f}..{5:.0f}]".format(xmin,xmax,ymin,ymax,zmin,zmax)

    # print maximum particle radius as a fraction of domain size
    rmax = max(r)
    xperc = 100*rmax/(xmax-xmin)
    yperc = 100*rmax/(ymax-ymin)
    zperc = 100*rmax/(zmax-zmin)
    print "Maximum particle radius: {0:.0f} ({1:.0f}% of X size, {2:.0f}% of Y size, {3:.0f}% of Z size)" \
                        .format(rmax, xperc,yperc,zperc)

    # produce basic plot if requested
    if basicplot:
        plt.figure(1, figsize=(12,12))
        plt.suptitle("Data file: {0}\nTotal of {1} particles".format(filename, len(r)), fontsize=16)

        # plot a histogram for the particle radii
        plt.subplot(221)
        plt.hist(r, bins=100)
        plt.title("Particle radius")

        # plot a histogram for the particle positions in each dimension
        plt.subplot(222)
        plt.hist(x, bins=100, range=(xmin,xmax))
        plt.title("Particle X coordinate")
        plt.subplot(223)
        plt.hist(y, bins=100, range=(ymin,ymax))
        plt.title("Particle Y coordinate")
        plt.subplot(224)
        plt.hist(z, bins=100, range=(zmin,zmax))
        plt.title("Particle Z coordinate")

        # show or save the figure
        if plotfile==None:
            plt.show()
        else:
            plt.savefig(plotfile)
            plt.close()

    # produce circle plot if requested
    if circleplot:
        plt.figure(1, figsize=(15,15))
        plt.suptitle("Data file: {0}\nTotal of {1} particles".format(filename, len(r)), fontsize=16)

        # plot a histogram for the particle radii
        plt.subplot(221)
        plt.hist(r, bins=100, log=True)
        plt.title("Particle radius")

        # plot projections of the particle positions/radii on each coordinate plane
        alpha = 0.01
        if len(r) < 20000: alpha=0.02
        if len(r) < 2000: alpha=0.05
        s = r * r;   # size is in square points so depends on scale of the plot --> not accurate
        plt.subplot(222, aspect='equal')
        plt.scatter(x, y, s, edgecolors='none', alpha=alpha)
        plt.title("Projection on XY plane")
        plt.subplot(223, aspect='equal')
        plt.scatter(x, z, s, edgecolors='none', alpha=alpha)
        plt.title("Projection on XZ plane")
        plt.subplot(224, aspect='equal')
        plt.scatter(y, z, s, edgecolors='none', alpha=alpha)
        plt.title("Projection on YZ plane")

        # show or save the figure
        if plotfile==None:
            plt.show()
        else:
            plt.savefig(plotfile)
            plt.close()

    # create smart grid and produce cell plot if requested
    if cellplot:
        # create smart grid
        print "Creating smart {0}x{0}x{0} grid ({1} cells)".format(gridsize, gridsize**3)
        xgrid = makegrid(x, xmin, xmax, gridsize)
        ygrid = makegrid(y, ymin, ymax, gridsize)
        zgrid = makegrid(z, zmin, zmax, gridsize)

        # create array to count number of particles affecting each cell
        partsincell = np.zeros((gridsize,gridsize,gridsize),'int')
        # create array to count number of cells affected by each particle
        cellsbypart = np.zeros(len(r),'int')

        # count cells (potentially) affected by each particle, and vice versa
        for index in range(len(r)):
            # get the overall index range
            xi1,xi2 = getindexrange(xgrid, x[index], r[index])
            yi1,yi2 = getindexrange(ygrid, y[index], r[index])
            zi1,zi2 = getindexrange(zgrid, z[index], r[index])

            # perform the counting
            if (accurate):
                # perform an accurate intersection test for each cell
                for xi in range(xi1,xi2):
                    for yi in range(yi1,yi2):
                        for zi in range(zi1,zi2):
                            if (intersects(xgrid[xi], xgrid[xi+1], ygrid[yi], ygrid[yi+1], zgrid[zi], zgrid[zi+1],  \
                                           x[index], y[index], z[index], r[index]) ):
                                partsincell[xi,yi,zi] += 1
                                cellsbypart[index] += 1
            else:
                # count all cells intersecting the enclosing rectangle
                partsincell[xi1:xi2,yi1:yi2,zi1:zi2] += 1
                cellsbypart[index] += (xi2-xi1)*(yi2-yi1)*(zi2-zi1)

            # report progress
            if index%50000==0:
                print " .. processed {0} particles".format(index)

        # print basic statistic
        print "Maximum {0} particles affecting a single cell".format( max(partsincell.flatten()) )

        # start the figure
        plt.figure(1, figsize=(12,12))
        plt.suptitle("Data file: {0}\nTotal of {1} particles -- Smart {2}x{2}x{2} grid ({3} cells)" \
                        .format(filename, len(r), gridsize, gridsize**3), fontsize=16)

        # plot a histogram for the number of particles affecting each cell
        plt.subplot(211)
        plt.hist(partsincell.flatten(), bins=100)
        plt.title("Number of particles affecting a single cell")

        # plot a histogram for the number of cells affected by each particle
        plt.subplot(212)
        plt.hist(cellsbypart, bins=100)
        plt.title("Number of cells affected by a single particle")

        # show or save the figure
        if plotfile==None:
            plt.show()
        else:
            plt.savefig(plotfile)
            plt.close()

    print "-----"

# -----------------------------------------------------------------
