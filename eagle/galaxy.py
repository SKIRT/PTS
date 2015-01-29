#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package eagle.galaxy Extracting data from EAGLE simulation results.
#
# EAGLE is a cosmological hydrodynamical simulation conducted by the VIRGO consortium. The simulation
# includes detailed recipes for star formation, gas cooling, stellar evolution, and feedback from supernovae and AGN.
# Using more than 10 billion particles, it numerically resolves thousands of galaxies in a representative cosmological
# volume that also contains groups and clusters.
#
# The EAGLE simulation output is stored in a (large) set of data files in the HDF5 format, documented at the
# <a href="http://www.hdfgroup.org/HDF5/">HFD5 home page</a>. The output is organized in \em snapshots, where
# each snapshot represents the state of the universe at a particular time (or equivalently, redshift).
#
# The classes and functions in this module allow extracting information relevant for SKIRT from the EAGLE output.
# Specifically, there are provisions for selecting gravitationally bound systems based on certain criteria
# (such as total mass or number of particles), and for writing the particles contained in such systems to a text
# file ready for SKIRT import.
#
# This module has to deal with three units systems: EAGLE snapshot units (documented through hdf5 attributes
# in the snapshot files), SKIRT import units (documented in the SKIRT SPH classes), and presentation units (used
# for displaying galaxy properties and specifying selection criteria in this module's Galaxy class).
# The following table lists the most important units in each system.
#
#<TABLE>
#<TR><TD><B>Physical Quantity</B></TD>  <TD><B>EAGLE snapshot</B></TD>
#                                       <TD><B>SKIRT import</B></TD>
#                                       <TD><B>Presentation</B></TD></TR>
#<TR><TD>position, size</TD>            <TD>\f$\textrm{Mpc}\,a\,h^{-1}\f$</TD>
#                                       <TD>\f$\textrm{pc}\f$</TD>
#                                       <TD>\f$\textrm{kpc}\f$</TD></TR>
#<TR><TD>mass</TD>                      <TD>\f$10^{10}\,\textrm{M}_\odot\,h^{-1}\f$</TD>
#                                       <TD>\f$\textrm{M}_\odot\f$</TD>
#                                       <TD>\f$10^9\,\textrm{M}_\odot\f$</TD></TR>
#<TR><TD>velocity</TD>                  <TD>\f$\textrm{km/s}\,a^{1/2}\f$</TD>
#                                       <TD>--</TD>
#                                       <TD>--</TD></TR>
#<TR><TD>time, age</TD>                 <TD>--</TD>
#                                       <TD>year</TD>
#                                       <TD>--</TD></TR>
#<TR><TD>temperature</TD>               <TD>K</TD>
#                                       <TD>--</TD>
#                                       <TD>--</TD></TR>
#</TABLE>
#
# Note the corrections for cosmological scale factor \f$a\f$ and hubble parameter \f$h\f$ in the EAGLE snapshot units.

# -----------------------------------------------------------------

import os
import os.path
import numpy as np
import h5py

import eagle.config as config
from pts.geometry import Transform
import eagle.starformation as sf

# -----------------------------------------------------------------

## This class represents the set of files containing a single EAGLE snapshot (i.e. the state of the universe
# at a particular redshift). The constructor locates all files, verifies that the snapshot is complete,
# and cashes some relevant information. The other functions return or export data extracted from the snapshot.
#
class Snapshot:

    ## The constructor locates all files containing particle data for the snapshot, verifies that the data is
    # consistent and complete, and cashes some relevant information. The arguments are:
    #  - eaglesim: the EAGLE simulation for this snaphot; this value is used as a key into the
    #               \em config.eagledata_path dictionary to obtain the path to the eagle data directory;
    #  - redshift: the redshift for this snapshot, with an accuracy of at least 3 digits after the decimal point;
    #    the default value is redshift zero.
    #
    # The specified redshift is converted into a filename pattern (formatted as zNNNpFFF where NNN is the portion
    # before the decimal point and FFF is the fractional portion), which is used to filter the files in the directory.
    # Any subdirectories with a name matching the redshift pattern are searched as well (recursively).
    def __init__(self, eaglesim, redshift=0):

        # format the redshift pattern
        intportion = int(redshift)
        fracportion = int((redshift-intportion)*1000+0.5)
        pattern = "z{0:0>3}p{1:0>3}".format(intportion,fracportion)

        # make a list of all hdf5 file paths that match the pattern
        self.directory = config.eagledata_path[eaglesim]
        paths = findhdf5files(self.directory, pattern)
        if len(paths) == 0: raise ValueError("No matching HDF5 files found")

        # open one of the files, arbitrarily, and retrieve the relevant header information
        hdf = h5py.File(paths[0], 'r')
        self.runlabel = hdf["Header"].attrs["RunLabel"]
        self.hubbleparam = hdf["Header"].attrs["HubbleParam"]
        self.boxsize = toparsec(hdf["Header"].attrs["BoxSize"], self.hubbleparam, 1)
        self.redshift = hdf["Header"].attrs["Redshift"]
        self.expansionfactor = hdf["Header"].attrs["ExpansionFactor"]

         # set up attribute dictionaries
        self.constants          = sf.readAttrs(hdf, field='Constants')
        self.header          = sf.readAttrs(hdf, field='Header')
        self.runpars          = sf.readAttrs(hdf, field='RuntimePars')

        # extract star formation law parameters
        self.schmidtpars      = sf.schmidtParameters(self.constants, self.runpars)

        # resampling parameters
        self.mcl          = 1.5e4     # in M_solar
        self.mpdr          = 1.5e5 # in M_solar
        self.agethresh           = 1e8     # in yr

        hdf.close()

        # verify that all files belong to the same snapshot
        for path in paths:
            hdf = h5py.File(path, 'r')
            if self.runlabel != hdf["Header"].attrs["RunLabel"]: raise ValueError("Run label does not match")
            if self.redshift != hdf["Header"].attrs["Redshift"]: raise ValueError("Redshift does not match")
            hdf.close()

        # look for a file that contains particle data, and retrieve the total number of such files
        self.numfiles = 0
        for path in paths:
            hdf = h5py.File(path, 'r')
            if "PartType0" in hdf.keys() and "PartType4" in hdf.keys():
                self.numfiles = hdf["Header"].attrs["NumFilesPerSnapshot"]
                hdf.close()
                break
            hdf.close()

        # build an ordered list of paths for all files containing particle data
        # the appropriate index is derived from the filename (portion between the last two periods)
        self.paths = [ None for x in range(self.numfiles) ]
        for path in paths:
            hdf = h5py.File(path, 'r')
            if "PartType0" in hdf.keys() and "PartType4" in hdf.keys():
                index = int(path.rsplit('.',2)[1])
                self.paths[index] = path
            hdf.close()

        # verify that the list of files is complete
        missing = [ x for x in self.paths if x==None ]
        if len(missing)>0: raise ValueError("There are " + str(len(missing)) + " missing particle files")

    ## This function returns the total number of star particles in the snapshot, as a long integer
    def numstarparticles(self):
        result = long(0)
        for path in self.paths:
            hdf = h5py.File(path, 'r')
            result += hdf["PartType4/ParticleIDs"].shape[0]
            hdf.close()
        return result

    ## This function returns the total number of gas particles in the snapshot, as a long integer
    def numgasparticles(self):
        result = long(0)
        for path in self.paths:
            hdf = h5py.File(path, 'r')
            result += hdf["PartType0/ParticleIDs"].shape[0]
            hdf.close()
        return result

    ## This function prints some information on the snapshot to the console.
    def printinfo(self):
        print "Directory: " + self.directory
        print "Run label: " + self.runlabel
        print "Box size:  {0:.1f} Mpc".format(self.boxsize/1e6)
        print "Redshift:  {0:.3f}".format(self.redshift)
        print "Expansion: {0:.1f}%".format(self.expansionfactor*100)
        print "Files:     {0}".format(self.numfiles)
        print "Star particles: {0:11,}".format(self.numstarparticles())
        print "Gas particles: {0:12,}".format(self.numgasparticles())

    ## This function saves a cross reference table for the snapshot in a file named "crossreference.dat"
    # in the \em config.database_path directory. For each particle subgroup with at least the specified
    # minimum number of star particles AND at least the specified minimum number of gas particles,
    # the cross reference identifies the hdf5 data files holding the corresponding particle information.
    # Each text line in the file represents a table entry and has the following columns:
    # group number, subgroup number, number of star particles, number of gas particles, list of file indices.
    def exportcrossreference(self, minparticles = 1000):
        # loop over all star and gas particles and count group membership, creating a dictionary
        # with key = ( group number, subgroup number )
        # and value = [ number of star particles, number of gas particles, set of file indices ]
        print "Building cross reference table..."
        subgroups = { }
        for index,path in zip(range(len(self.paths)),self.paths):
            if index%10==0: print " opening file", index+1, "of", len(self.paths)
            hdf = h5py.File(path, 'r')
            for groupnumber,subnumber in zip(hdf["PartType4/GroupNumber"][:], hdf["PartType4/SubGroupNumber"][:]):
                if subnumber < 1073741824:   # (1073741824 means "not in subgroup")
                    key = (abs(groupnumber), subnumber)
                    if not key in subgroups: subgroups[key] = [0, 0, set() ]
                    subgroups[key][0] += 1
                    subgroups[key][2].add(index)
            for groupnumber,subnumber in zip(hdf["PartType0/GroupNumber"][:],hdf["PartType0/SubGroupNumber"][:]):
                if subnumber < 1073741824:   # (1073741824 means "not in subgroup")
                    key = (abs(groupnumber), subnumber)
                    if not key in subgroups: subgroups[key] = [0, 0, set() ]
                    subgroups[key][1] += 1
                    subgroups[key][2].add(index)
            hdf.close()
        print " done."

        # save cross reference entries for all sufficiently large subgroups
        print "Saving cross reference table..."
        out = open(os.path.join(config.database_path,"crossreference.dat"),'w')
        for key in sorted(subgroups.keys()):
            value = subgroups[key]
            if value[0] >= minparticles and value[1] >= minparticles:
                line = str(key[0]) + ' ' + str(key[1]) + ' ' + str(value[0]) + ' ' + str(value[1])
                for index in sorted([i for i in value[2]]): line += ' ' + str(index)
                out.write(line + '\n')
        out.close()
        print " done."

    ## This function saves catalog information for the snapshot in a file named "catalog.dat" in the
    # \em config.database_path directory.
    # The function expects a file named "crossreference.dat" to be present in the \em config.database_path directory,
    # in the format as exported by the exportcrossreference() function. The created catalog contains summary data
    # for each subgroup listed in the cross reference. This data can be used to easily calculate properties such as
    # total mass, barycenter, or spin factor. Each line in the catalog file represents a catalog entry (i.e. subgroup)
    # and contains a fixed number of columns, as listed in the following table. All data is in snapshot units.
    #
    #<TABLE>
    #<TR><TD><B>column index</B></TD> <TD><B>field description</B></TD></TR>
    #<TR><TD>0</TD>     <TD>group number</TD></TR>
    #<TR><TD>1</TD>     <TD>subgroup number</TD></TR>
    #<TR><TD>2</TD>     <TD>number of star particles</TD></TR>
    #<TR><TD>3</TD>     <TD>number of gas particles</TD></TR>
    #<TR><TD>4</TD>     <TD>\f$\sum_* m_i \f$ (stellar mass)</TD></TR>
    #<TR><TD>5</TD>     <TD>\f$\sum_\textrm{gas} m_i \f$ (gas mass)</TD></TR>
    #<TR><TD>6</TD>     <TD>\f$\min (r_{i,x})\f$</TD></TR>
    #<TR><TD>7</TD>     <TD>\f$\min (r_{i,y})\f$</TD></TR>
    #<TR><TD>8</TD>     <TD>\f$\min (r_{i,z})\f$</TD></TR>
    #<TR><TD>9</TD>     <TD>\f$\max (r_{i,x})\f$</TD></TR>
    #<TR><TD>10</TD>    <TD>\f$\max (r_{i,y})\f$</TD></TR>
    #<TR><TD>11</TD>    <TD>\f$\max (r_{i,z})\f$</TD></TR>
    #<TR><TD>12</TD>    <TD>\f$\sum m_i r_{i,x}\f$</TD></TR>
    #<TR><TD>13</TD>    <TD>\f$\sum m_i r_{i,y}\f$</TD></TR>
    #<TR><TD>14</TD>    <TD>\f$\sum m_i r_{i,z}\f$</TD></TR>
    #<TR><TD>15</TD>    <TD>\f$\sum m_i v_{i,x}\f$</TD></TR>
    #<TR><TD>16</TD>    <TD>\f$\sum m_i v_{i,y}\f$</TD></TR>
    #<TR><TD>17</TD>    <TD>\f$\sum m_i v_{i,z}\f$</TD></TR>
    #<TR><TD>18</TD>    <TD>\f$\sum m_i(\vec{r}_i\times\vec{v}_i)_x\f$</TD></TR>
    #<TR><TD>19</TD>    <TD>\f$\sum m_i(\vec{r}_i\times\vec{v}_i)_y\f$</TD></TR>
    #<TR><TD>20</TD>    <TD>\f$\sum m_i(\vec{r}_i\times\vec{v}_i)_y\f$</TD></TR>
    #</TABLE>
    #
    def exportcatalog(self):
        # load the cross reference table into a dictionary with as key (group number, subgroup number) and as
        # value a Galaxy instance with zeroed initial fields.
        print "Loading cross reference table..."
        catalog = { }
        for line in open(os.path.join(config.database_path,"crossreference.dat"),'r'):
            fields = line.split()
            key = (int(fields[0]), int(fields[1]))
            catalog[key] = Galaxy(self, fields[0:4])
        print " done:", len(catalog), "entries."

        # loop over all star and gas particles and add their information to the corresponding galaxy, if present
        print "Building catalog..."
        for index,path in zip(range(len(self.paths)),self.paths):
            if index%10==0: print " opening file", index+1, "of", len(self.paths)
            hdf = h5py.File(path, 'r')
            stars = hdf["PartType4"]
            for groupnumber,subnumber,mass,r,v in zip(stars["GroupNumber"][:], stars["SubGroupNumber"][:],
                                                  stars["Mass"][:], stars["Coordinates"][:], stars["Velocity"][:]):
                key = (abs(groupnumber), subnumber)
                if key in catalog:
                    catalog[key].addparticle(mass, 0., r, v)
            gas = hdf["PartType0"]
            for groupnumber,subnumber,mass,r,v in zip(gas["GroupNumber"][:], gas["SubGroupNumber"][:],
                                                  gas["Mass"][:], gas["Coordinates"][:], gas["Velocity"][:]):
                key = (abs(groupnumber), subnumber)
                if key in catalog:
                    catalog[key].addparticle(0., mass, r, v)
            hdf.close()
        print " done."
        print ""

        # save the catalog entries
        print "Saving catalog..."
        out = open(os.path.join(config.database_path,"catalog.dat"),'w')
        for key in sorted(catalog.keys()):
            out.write(catalog[key].catalogline() + '\n')
        out.close()
        print " done."

    ## This function returns a GalaxyList object holding all galaxies in the snapshot listed in the catalog.
    # The function expects a file named "catalog.dat" to be present in the \em config.database_path directory,
    # in the format as exported by the exportcatalog() function.
    def galaxies(self):
        # build a galaxy list from the catalog
        galaxylist = GalaxyList()
        for line in open(os.path.join(config.database_path,"catalog.dat"),'r'):
            galaxylist.addgalaxy(Galaxy(self, line.split()))
        return galaxylist

# -----------------------------------------------------------------

## This class represents a single particle subgroup (or "galaxy") in a snapshot, identified by group number
# and subgroup number. It serves to accumulate the information needed for a catalog entry, and to extract
# information about the galaxy from the catalog entry and/or from underlying snapshot data files.
class Galaxy:

    ## This constructor initializes all property fields for the galaxy. The first argument provides a reference
    # to the snapshot containing this galaxy. The second argument is a list of field values used to initialize
    # the galaxy properties. The first four values must be: group number, subgroup number, number of star particles,
    # and number of gas particles. If there are no additional values, the other properties are cleared.
    # If there are additional values, they must correspond to the fields of a catalog entry.
    def __init__(self, snapshot, fields):
        # remember reference to snapshot
        self.snapshot = snapshot
        # initialize the properties that must always be provided as arguments to the constructor
        self.groupnumber, self.subgroupnumber, self.numstarparticles, self.numgasparticles = map(int, fields[0:4])
        # initialize the other properties, either from the constructor arguments or to zero
        # note: these fields are in snapshot units!
        if len(fields) > 4:
            self.sum_m_star, self.sum_m_gas,  \
            self.min_r_x, self.min_r_y, self.min_r_z,  \
            self.max_r_x, self.max_r_y, self.max_r_z,  \
            self.sum_mr_x, self.sum_mr_y, self.sum_mr_z,  \
            self.sum_mv_x, self.sum_mv_y, self.sum_mv_z,  \
            self.sum_mrv_x, self.sum_mrv_y, self.sum_mrv_z = map(float, fields[4:])
        else:
            self.sum_m_star = self.sum_m_gas = 0.
            self.min_r_x = self.min_r_y = self.min_r_z = float('inf')
            self.max_r_x = self.max_r_y = self.max_r_z = float('-inf')
            self.sum_mr_x = self.sum_mr_y = self.sum_mr_z = 0.
            self.sum_mv_x = self.sum_mv_y = self.sum_mv_z = 0.
            self.sum_mrv_x = self.sum_mrv_y = self.sum_mrv_z = 0.

    ## This function adds the statistics for a single star or gas particle to the catalog properties for this entry.
    # The r and v arguments must be (x,y,z)-tuples.
    def addparticle(self, m_star, m_gas, r, v):
        self.sum_m_star += m_star
        self.sum_m_gas += m_gas
        self.min_r_x = min(self.min_r_x, r[0])
        self.min_r_y = min(self.min_r_y, r[1])
        self.min_r_z = min(self.min_r_z, r[2])
        self.max_r_x = max(self.max_r_x, r[0])
        self.max_r_y = max(self.max_r_y, r[1])
        self.max_r_z = max(self.max_r_z, r[2])
        m = m_star + m_gas
        self.sum_mr_x += m * r[0]
        self.sum_mr_y += m * r[1]
        self.sum_mr_z += m * r[2]
        self.sum_mv_x += m * v[0]
        self.sum_mv_y += m * v[1]
        self.sum_mv_z += m * v[2]
        self.sum_mrv_x += m * (r[1]*v[2]-r[2]*v[1])
        self.sum_mrv_y += m * (r[2]*v[0]-r[0]*v[2])
        self.sum_mrv_z += m * (r[0]*v[1]-r[1]*v[0])

    ## This function returns a string containing the properties for this catalog entry in a format appropriate
    # for saving into the catalog text file.
    def catalogline(self):
        return (21*"{:1.9g} ").format( \
            self.groupnumber, self.subgroupnumber, self.numstarparticles, self.numgasparticles,
            self.sum_m_star, self.sum_m_gas,
            self.min_r_x, self.min_r_y, self.min_r_z,
            self.max_r_x, self.max_r_y, self.max_r_z,
            self.sum_mr_x, self.sum_mr_y, self.sum_mr_z,
            self.sum_mv_x, self.sum_mv_y, self.sum_mv_z,
            self.sum_mrv_x, self.sum_mrv_y, self.sum_mrv_z )

    ## This function returns the 2-tuple (group number, subgroup number) identifying the galaxy in the snapshot
    def key(self):
        return (self.groupnumber, self.subgroupnumber)

    ## This function returns the total number of particles (star + gas).
    def numparticles(self):
        return self.numstarparticles + self.numgasparticles

    ## This function returns the stellar mass of the galaxy \f$\sum_* m_i \f$, in solar mass units.
    def starmass(self):
        return tosolar(self.sum_m_star, self.snapshot.hubbleparam)

    ## This function returns the gas mass of the galaxy \f$\sum_\textrm{gas} m_i \f$, in solar mass units.
    def gasmass(self):
        return tosolar(self.sum_m_gas, self.snapshot.hubbleparam)

    ## This function returns the total mass of the galaxy \f$\sum m_i \f$, taking into account stars and gas,
    # in solar mass units.
    def totalmass(self):
        return tosolar(self.sum_m_star + self.sum_m_gas, self.snapshot.hubbleparam)

    ## This function returns the galaxy's center of mass (COM) position \f$\frac{\sum m_i \vec{r}_i}{\sum m_i}\f$,
    # taking into account stars and gas, as a 3-tuple with the (x,y,z) coordinates in parsec
    def barycenter(self):
        mass = self.sum_m_star + self.sum_m_gas
        return ( toparsec(self.sum_mr_x/mass, self.snapshot.hubbleparam, self.snapshot.expansionfactor),
                 toparsec(self.sum_mr_y/mass, self.snapshot.hubbleparam, self.snapshot.expansionfactor),
                 toparsec(self.sum_mr_z/mass, self.snapshot.hubbleparam, self.snapshot.expansionfactor) )

    ## This function returns the galaxy's center of mass (COM) velocity \f$\frac{\sum m_i \vec{v}_i}{\sum m_i}\f$,
    # taking into account stars and gas, as a 3-tuple with the (x,y,z) components in km/s.
    def centralvelocity(self):
        mass = self.sum_m_star + self.sum_m_gas
        return ( tokms(self.sum_mv_x/mass, self.snapshot.expansionfactor),
                 tokms(self.sum_mv_y/mass, self.snapshot.expansionfactor),
                 tokms(self.sum_mv_z/mass, self.snapshot.expansionfactor) )

    ## This function returns an estimate for the radius of the galaxy, in parsec, determined as the largest
    # star or gas particle distance from the center of mass along one the coordinate axes.
    def radius(self):
        mass = self.sum_m_star + self.sum_m_gas
        bx,by,bz = self.sum_mr_x/mass, self.sum_mr_y/mass, self.sum_mr_z/mass
        maxr = np.max(np.abs(np.array(( bx-self.min_r_x, bx-self.max_r_x,
                                        by-self.min_r_y, by-self.max_r_y,
                                        bz-self.min_r_z, bz-self.max_r_z ))))
        return toparsec(maxr, self.snapshot.hubbleparam, self.snapshot.expansionfactor)

    ## This function returns the galaxy's angular momentum, taking into account stars and gas,
    # as a 3-tuple with (x,y,z) components in units of \f$\textrm{M}_\odot\,\textrm{pc}\,\textrm{km/s}\f$.
    # The angular momentum is calculated as
    # \f[\vec{L}=\sum m_i (\vec{r}_i-\vec{R}_\textrm{com})\times(\vec{v}_i-\vec{V}_\textrm{com})
    #           =\sum m_i (\vec{r}_i\times\vec{v}_i) - (\sum m_i)(\vec{R}_\textrm{com}\times\vec{V}_\textrm{com})\f]
    # where \f$\vec{R}_\textrm{com}=\frac{\sum m_i \vec{r}_i}{\sum m_i}\f$
    # and \f$\vec{V}_\textrm{com}=\frac{\sum m_i \vec{v}_i}{\sum m_i}\f$.
    def angularmomentum(self):
        M = self.sum_m_star + self.sum_m_gas
        MRV1 = np.array((self.sum_mrv_x, self.sum_mrv_y, self.sum_mrv_z))
        MRV2 = cross((self.sum_mr_x, self.sum_mr_y, self.sum_mr_z), (self.sum_mv_x, self.sum_mv_y, self.sum_mv_z)) / M
        L = MRV1 - MRV2
        h = self.snapshot.hubbleparam
        a = self.snapshot.expansionfactor
        return tosolar(toparsec(tokms(L, a), h, a), h)

    ## This function returns the direction of the galaxy's angular momentum, taking into account stars and gas,
    # as a 3-tuple with normalized direction cosines. If the angular momentum is too small to allow normalization,
    # the function returns (0,0,0).
    def rotationaxis(self):
        L = self.angularmomentum()
        norm = np.sqrt(np.sum(L*L))
        if norm < 1e-12: return (0,0,0)
        return L/norm

    ## This function returns a dimensionless measure for the galaxy's spin, calculated as
    # \f[\lambda=\frac{|\vec{L}|}{MR\,V}=\frac{|\vec{L}|}{\sqrt{GM^3R}}\f]
    # where M stands for the galaxy's total stellar and gas mass, R stands for the galaxy's estimated radius,
    # V stands for the circle velocity at radius R, and G is the gravitational constant.
    # The circle velocity is estimated using \f$V=\sqrt{GM/R}\f$.
    def spin(self):
        M = self.totalmass()         # mass: Msun
        L = self.angularmomentum()
        L = np.sqrt(np.sum(L*L))     # angular momentum: Msun pc km/s
        R = self.radius()            # radius: pc
        G = 4.302e-3                 # gravitational constant: pc Msun^-1 (km/s)^2
        return L / np.sqrt(G*M*M*M*R)  # spin parameter: dimensionless

    ## This function prints a header line for the information lines printed by the printinfo() function
    def printinfotitle(self):
        print "group subgr  #star part.  #gas part.  star mass  gas mass    radius      spin"
        print "                                       (10^9 solar mass)      (kpc)"

    ## This function prints some information on the galaxy on a single line
    def printinfo(self):
        print " {0:4} {1:3}  {2:11,} {3:11,}   {4:9.1f} {5:9.1f}  {6:9.1f}  {7:9.3f}"  \
            .format(self.groupnumber, self.subgroupnumber,
                    self.numstarparticles, self.numgasparticles,
                    self.starmass()/1e9, self.gasmass()/1e9, self.radius()/1e3, self.spin())

    ## This function exports particle data for the galaxy to two text files that can be imported in SKIRT.
    # The function expects a file named "crossreference.dat" to be present in the \em config.database_path directory,
    # in the format as exported by the exportcrossreference() function. The exported files are placed in the
    # specified directory and are named "galaxy_G_S_stars.dat" and "galaxy_G_S_gas.dat", where G and S are replaced
    # by the galaxy's group number and subgroup number. The file format is as described for SKIRT SPH import.
    def export(self, directory=""):
        # get the list of relevant snapshot files
        fileindices = [ ]
        for line in open(os.path.join(config.database_path,"crossreference.dat"),'r'):
            fields = line.split()
            if self.groupnumber==int(fields[0]) and self.subgroupnumber==int(fields[1]):
                fileindices = map(int, fields[4:])
        if len(fileindices)==0: raise ValueError("No cross reference info found")
        print "Exporting galaxy ({0},{1}) from {2} files...".format(  \
                    self.groupnumber, self.subgroupnumber, len(fileindices))

        # ---- get the particle data

        # initialise star and gas dictionaries
        sdat        = {}
        gdat        = {}
        yngstars    = {}
        hiiregions  = {}
        sdat['r']        = np.column_stack([[],[],[]])
        sdat['v']        = np.column_stack([[],[],[]])
        sdat['h']        = np.array([])
        sdat['im']       = np.array([])
        sdat['m']        = np.array([])
        sdat['Z']        = np.array([])
        sdat['born']     = np.array([])
        sdat['rho_born'] = np.array([])
        gdat['r']        = np.column_stack([[],[],[]])
        gdat['v']        = np.column_stack([[],[],[]])
        gdat['h']        = np.array([])
        gdat['m']        = np.array([])
        gdat['Z']        = np.array([])
        gdat['T']        = np.array([])
        gdat['rho']      = np.array([])
        gdat['sfr']      = np.array([])

        # read particle informaton
        for index in fileindices:

            hdf   = h5py.File(self.snapshot.paths[index], 'r')
            stars = hdf["PartType4"]
            gas   = hdf["PartType0"]

            # index for subhalo stars
            insubhalo = (stars["GroupNumber"][:] == self.groupnumber) & (stars["SubGroupNumber"][:] == self.subgroupnumber)
            sdat['r']        = np.concatenate((sdat['r'],        stars["Coordinates"][:][insubhalo]), axis=0)
            sdat['h']        = np.concatenate((sdat['h'],        stars["SmoothingLength"][:][insubhalo]), axis=0)
            sdat['im']       = np.concatenate((sdat['im'],       stars["InitialMass"][:][insubhalo]), axis=0)
            sdat['m']        = np.concatenate((sdat['m'],        stars["Mass"][:][insubhalo]), axis=0)
            sdat['v']        = np.concatenate((sdat['v'],        stars["Velocity"][:][insubhalo]), axis=0)
            sdat['Z']        = np.concatenate((sdat['Z'],        stars["SmoothedMetallicity"][:][insubhalo]), axis=0)
            sdat['born']     = np.concatenate((sdat['born'],     stars["StellarFormationTime"][:][insubhalo]), axis=0)
            sdat['rho_born'] = np.concatenate((sdat['rho_born'], stars["BirthDensity"][:][insubhalo]), axis=0)

            # index for subhalo gas
            insubhalo = (gas["GroupNumber"][:] == self.groupnumber) & (gas["SubGroupNumber"][:] == self.subgroupnumber)
            gdat['r']        = np.concatenate((gdat['r'],   gas["Coordinates"][:][insubhalo]), axis=0)
            gdat['h']        = np.concatenate((gdat['h'],   gas["SmoothingLength"][:][insubhalo]), axis=0)
            gdat['m']        = np.concatenate((gdat['m'],   gas["Mass"][:][insubhalo]), axis=0)
            gdat['v']        = np.concatenate((gdat['v'],   gas["Velocity"][:][insubhalo]), axis=0)
            gdat['Z']        = np.concatenate((gdat['Z'],   gas["SmoothedMetallicity"][:][insubhalo]), axis=0)
            gdat['T']        = np.concatenate((gdat['T'],   gas["Temperature"][:][insubhalo]), axis=0)
            gdat['rho']      = np.concatenate((gdat['rho'], gas["Density"][:][insubhalo]), axis=0)
            gdat['sfr']      = np.concatenate((gdat['sfr'], gas["StarFormationRate"][:][insubhalo]), axis=0)

            hdf.close()

        # convert units
        sdat['r']        = toparsec(sdat['r'], self.snapshot.hubbleparam, self.snapshot.expansionfactor)
        sdat['r']        = periodicCorrec(sdat['r'], self.snapshot.boxsize)
        sdat['h']        = toparsec(sdat['h'], self.snapshot.hubbleparam, self.snapshot.expansionfactor)
        sdat['im']       = tosolar(sdat['im'], self.snapshot.hubbleparam)
        sdat['Z']        = sdat['Z']
        sdat['t']        = age(sdat['born']) - age(self.snapshot.expansionfactor)
        sdat['rho_born'] *= 6.7699e-31
        gdat['r']        = toparsec(gdat['r'], self.snapshot.hubbleparam, self.snapshot.expansionfactor)
        gdat['r']        = periodicCorrec(gdat['r'], self.snapshot.boxsize)
        gdat['h']        = toparsec(gdat['h'], self.snapshot.hubbleparam, self.snapshot.expansionfactor)
        gdat['m']        = tosolar(gdat['m'], self.snapshot.hubbleparam)
        gdat['Z']        = gdat['Z']
        gdat['T']        = gdat['T']
        gdat['sfr']      = gdat['sfr']
        gdat['rho']      = togcm3(gdat['rho'], self.snapshot.hubbleparam, self.snapshot.expansionfactor)

        # density conversion from g cm^-3 to M_sun Mpc^-3
        densconv = ((self.snapshot.constants['CM_PER_MPC']/1.e6)**3) / self.snapshot.constants['SOLAR_MASS']

        # calculate new fields
        sdat['P']        = sf.getPtot(sdat['rho_born'], self.snapshot.schmidtpars)
        #sdat['sfr_born'] = sf.getSFR(sdat['rho_born'], sdat['m'], self.snapshot.schmidtpars)
        sdat['h_mapp']   = (self.snapshot.mpdr / (4 * np.pi * sdat['rho_born'] * densconv))**(1/3.)

        gdat['P']        = sf.getPtot(gdat['rho'].copy(), self.snapshot.schmidtpars)
        gdat['h_mapp']   = (self.snapshot.mpdr / (4 * np.pi * gdat['rho'] * densconv))**(1/3.)

        # ---- convert to Local Galactic Coordinates (LGC)

        # calculate stellar centre of mass and translational velocity using shrinking aperture technique
        com, v_bar = shrinkingCentroid(sdat['r'], sdat['m'], sdat['v'])

        # find unit rotation axis vector, choosing to use only stellar information and an aperture of 30 kpc
        n_rot = rotAxis(sdat['r'], sdat['v'], sdat['m'], com, v_bar, apt = 3.e4, aptfrac = 0.08)

        # set up transform object
        transf = Transform()
        bx, by, bz = com[0], com[1], com[2]
        transf.translate(-bx, -by, -bz)
        a, b, c = n_rot[0], n_rot[1], n_rot[2]
        v = np.sqrt(b*b+c*c)
        if v > 0.3:
            transf.rotateX(c/v, -b/v)
            transf.rotateY(v, -a)
        else:
            v = np.sqrt(a*a+c*c)
            transf.rotateY(c/v, -a/v)
            transf.rotateX(v, -b)

        # transform coordinates
        sdat['r'],w = transf.transform_vec(sdat['r'][:,0],sdat['r'][:,1],sdat['r'][:,2], np.ones(sdat['r'].shape[0]))
        gdat['r'],w = transf.transform_vec(gdat['r'][:,0],gdat['r'][:,1],gdat['r'][:,2], np.ones(gdat['r'].shape[0]))

        # ---- resample star forming regions

        # set up GALAXEV array
        bcstars = np.column_stack([[],[],[],[],[],[],[]])

        # set up MAPPINGS-III array
        mapstars = np.column_stack([[],[],[],[],[],[],[],[],[]])

        # set up dust array
        dust = np.column_stack([[],[],[],[],[],[],[]])

        # index for particles to resample
        issf = gdat['sfr'] > 0.
        isyoung = sdat['t'] < self.snapshot.agethresh

        # append older stars to GALAXEV array
        bcstars = np.concatenate((bcstars, np.column_stack([sdat['r'], sdat['h'], sdat['im'], sdat['Z'], sdat['t']])[~isyoung]), axis=0)

        # append non-SF gas data to dust array
        dust = np.concatenate((dust, np.column_stack([gdat['r'], gdat['h'], gdat['m'], gdat['Z'], gdat['T']])[~issf].copy()), axis=0)

        # resample stars
        if isyoung.any():
            for k in sdat.keys():
                sdat[k] = sdat[k][isyoung].copy()

            # calculate SFR at birth of young star particles in M_sun / yr
            sdat['sfr']       = sf.getSFR(sdat['rho_born'], sdat['im'], self.snapshot.schmidtpars)

            resmpstr = sf.stochResamp(sdat['sfr'], sdat['im'], delm=self.snapshot.mcl, thresh_age=self.snapshot.agethresh)
            ms, ts, sfrs, idxs, mdiffs = resmpstr
            isinfant = ts < 1.e7

            yngstars['r']     = sdat['r'][idxs][~isinfant]
            yngstars['h']     = sdat['h'][idxs][~isinfant]
            yngstars['im']    = ms[~isinfant]
            yngstars['Z']     = sdat['Z'][idxs][~isinfant]
            yngstars['t']     = ts[~isinfant]

            bcstars = np.concatenate((bcstars, np.column_stack([yngstars['r'], yngstars['h'], yngstars['im'], yngstars['Z'], yngstars['t']])), axis=0)

            hiiregions['r']     = sdat['r'][idxs][isinfant]
            hiiregions['h']     = sdat['h'][idxs][isinfant]
            hiiregions['SFR']   = ms[isinfant] * 1e-7               # Assume constant SFR over 10 Myr, so sfr =  mcl / 10 Myr
            hiiregions['Z']     = sdat['Z'][idxs][isinfant]
            hiiregions['P']     = sdat['P'][idxs][isinfant] * 0.1   # Convert to Pa for output
            hiiregions['logC']  = 0.6*np.log10(ms[isinfant]) + 0.4*np.log10(hiiregions['P']) - 0.4*np.log10(self.snapshot.constants['BOLTZMANN'])
            hiiregions['h_mapp'] = sdat['h_mapp'][idxs][isinfant]

            # Set 0.2 as fiducial value (Jonsson (2010))
            hiiregions['fPDR']  = np.array([0.2]*hiiregions['P'].size)

            # append to MAPPINGSIII array
            mapstars = np.concatenate((mapstars, np.column_stack([hiiregions['r'], hiiregions['h_mapp'], hiiregions['SFR'],
                                                                  hiiregions['Z'], hiiregions['logC'], hiiregions['P'],
                                                                  hiiregions['fPDR']])), axis=0)
            # T is set to 0 K as T unavailable for resampled star particles
            sdat['T'] = np.zeros(sdat['im'].shape[0])

            # add unspent young star particle material to dust array.
            dust = np.concatenate((dust, np.column_stack([sdat['r'], sdat['h'], sdat['im'] - mdiffs, sdat['Z'], sdat['T']]).copy()), axis=0)

        # resample gas
        if issf.any():
            for k in gdat.keys():
                gdat[k] = gdat[k][issf].copy()

            resmpstr = sf.stochResamp(gdat['sfr'], gdat['m'], delm=self.snapshot.mcl, thresh_age=self.snapshot.agethresh)
            ms, ts, sfrs, idxs, mdiffs = resmpstr
            isinfant = ts < 1.e7

            yngstars['r']       = gdat['r'][idxs][~isinfant]
            yngstars['h']       = gdat['h'][idxs][~isinfant]
            yngstars['im']      = ms[~isinfant]
            yngstars['Z']       = gdat['Z'][idxs][~isinfant]
            yngstars['t']       = ts[~isinfant]

            bcstars = np.concatenate((bcstars, np.column_stack([yngstars['r'], yngstars['h'], yngstars['im'], yngstars['Z'], yngstars['t']])), axis=0)

            hiiregions['r']     = gdat['r'][idxs][isinfant]
            hiiregions['h']     = gdat['h'][idxs][isinfant]
            hiiregions['SFR']   = ms[isinfant] * 1e-7               # Assume constant SFR over 10 Myr, so sfr =  mcl / 10 Myr
            hiiregions['Z']     = gdat['Z'][idxs][isinfant]
            hiiregions['P']     = gdat['P'][idxs][isinfant] * 0.1     # convert to Pa
            hiiregions['logC']  = 0.6*np.log10(ms[isinfant]) + 0.4*np.log10(hiiregions['P']) - 0.4*np.log10(self.snapshot.constants['BOLTZMANN'])
            hiiregions['h_mapp'] = gdat['h_mapp'][idxs][isinfant]

             # Set 0.2 as fiducial value a la Jonsson (2010)
            hiiregions['fPDR']   = np.array([0.2]*hiiregions['P'].size)

            # append to MAPPINGSIII array
            mapstars = np.concatenate((mapstars, np.column_stack([hiiregions['r'], hiiregions['h_mapp'], hiiregions['SFR'],
                                                                  hiiregions['Z'], hiiregions['logC'], hiiregions['P'],
                                                                  hiiregions['fPDR']])), axis=0)

            # add unspent SF gas material to dust array
            dust = np.concatenate((dust, np.column_stack([gdat['r'], gdat['h'], gdat['m'] - mdiffs, gdat['Z'], gdat['T']]).copy()), axis=0)

        # ---- write output files

        # open output files
        starsfilename = "galaxy_{0}_{1}_stars.dat".format(self.groupnumber, self.subgroupnumber)
        gasfilename = "galaxy_{0}_{1}_gas.dat".format(self.groupnumber, self.subgroupnumber)
        hiifilename = "galaxy_{0}_{1}_hii.dat".format(self.groupnumber, self.subgroupnumber)
        starsfile = open(os.path.join(directory,starsfilename), 'w')
        starsfile.write('# SPH Star Particles\n')
        starsfile.write('# Extracted from EAGLE HDF5 snapshot to SKIRT6 format\n')
        starsfile.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) t(yr)\n')
        gasfile = open(os.path.join(directory,gasfilename), 'w')
        gasfile.write('# SPH Gas Particles\n')
        gasfile.write('# Extracted from EAGLE HDF5 snapshot to SKIRT6 format\n')
        gasfile.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) T(K) SFR(Msun/yr)\n')
        hiifile = open(os.path.join(directory,hiifilename), 'w')
        hiifile.write('# SPH Hii Particles\n')
        hiifile.write('# Extracted from EAGLE HDF5 snapshot to SKIRT6 format\n')
        hiifile.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) SFR(Msun/yr) Z(0-1) logC P(Pa) f_PDR\n')

        # save particle data
        np.savetxt(starsfile, bcstars, fmt=['%f']*7)
        np.savetxt(gasfile, dust, fmt=['%f']*7)
        np.savetxt(hiifile, mapstars, fmt=['%f']*7+['%e','%f'])

        # close output files
        starsfile.close()
        gasfile.close()
        hiifile.close()

# -----------------------------------------------------------------

## This class holds a list of galaxy objects. It provides facilities to remove galaxies from the list
# based on certain criteria, and to pretty-print information on the remaining items.
class GalaxyList:

    ## This constructor creates an empty galaxy list.
    def __init__(self):
        self.galaxies = [ ]

    ## This function adds a galaxy object to the list.
    def addgalaxy(self, galaxy):
        self.galaxies += [ galaxy ]

    ## This function returns the specified galaxy object, if it is in the list; otherwise returns None.
    def galaxy(self, groupnumber, subgroupnumber):
        for galaxy in self.galaxies:
            if galaxy.groupnumber == groupnumber and galaxy.subgroupnumber == subgroupnumber:
                return galaxy
        return None;

    ## This function removes all galaxies from the list that have a total mass above the specified value,
    # in \f$10^9\,\textrm{M}_\odot\f$.
    def remove_starmass_above(self, mass):
        self.galaxies = [ galaxy for galaxy in self.galaxies if galaxy.starmass() < mass*1e9 ]

    ## This function removes all galaxies from the list that have a total mass above the specified value,
    # in \f$10^9\,\textrm{M}_\odot\f$.
    def remove_gasmass_above(self, mass):
        self.galaxies = [ galaxy for galaxy in self.galaxies if galaxy.gasmass() < mass*1e9 ]

    ## This function removes all galaxies from the list that have a total mass above the specified value,
    # in \f$10^9\,\textrm{M}_\odot\f$.
    def remove_totalmass_above(self, mass):
        self.galaxies = [ galaxy for galaxy in self.galaxies if galaxy.totalmass() < mass*1e9 ]

    ## This function removes all galaxies from the list that have an estimated radius above the specified value,
    # in kpc.
    def remove_radius_above(self, radius):
        self.galaxies = [ galaxy for galaxy in self.galaxies if galaxy.radius() < radius*1e3 ]

    ## This function removes all galaxies from the list that have an estimated spin measure below the specified value.
    def remove_spin_below(self, spin):
        self.galaxies = [ galaxy for galaxy in self.galaxies if galaxy.spin() > spin ]

    ## This function prints some basic properties for all galaxies currently in the list. The output has
    # one line per galaxy and is sorted on descending total number of particles.
    def printinfo(self):
        if len(self.galaxies) > 0:
            print "There are {0} galaxies in the list".format(len(self.galaxies))
            self.galaxies.sort(key=Galaxy.numparticles, reverse=True)
            self.galaxies[0].printinfotitle()
            for galaxy in self.galaxies:
                galaxy.printinfo()
        else:
            print "There are no galaxies in the list"

    ## This function exports particle data for all galaxies currently in the list, as described for the
    # Galaxy.export() function.
    def export(self):
        print "Exporting {0} galaxies...".format(len(self.galaxies))
        self.galaxies.sort(key=Galaxy.numparticles, reverse=True)
        for galaxy in self.galaxies:
            galaxy.export()
        print " done."

# -----------------------------------------------------------------

## This private helper function recursively locates HDF5 files in a nested directory structure,
# as long as the directory or file name contains the specified pattern string.
def findhdf5files(dirpath, pattern):
    result = []
    for name in os.listdir(dirpath):
        if name.find(pattern) >= 0:
            path = os.path.join(dirpath,name);
            if os.path.isdir(path):
                result += findhdf5files(path, pattern)
            if os.path.isfile(path) and name.lower().endswith(".hdf5") and name.find("particles") >= 0:
                result += [ path ]
    return result

# -----------------------------------------------------------------

## This private helper function computes the cross product of two vectors. The incoming vectors must be given as
# 3-tuples (x,y,z); the result is returned as a numpy array with three elements so that it can be used with
# component-wise vector calculations.
def cross(a, b):
    return np.array((a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]))

# -----------------------------------------------------------------

## This private helper function converts a length/distance/size value from EAGLE snapshot units to parsec
def toparsec(eaglevalue, hubbleparam, expansionfactor):
    return eaglevalue * (1e6 * expansionfactor / hubbleparam)

## This private helper function converts a mass value from EAGLE snapshot units to solar masses
def tosolar(eaglevalue, hubbleparam):
    return eaglevalue * (1e10 / hubbleparam)

## This private helper function converts a velocity value from EAGLE snapshot units to km/s
def tokms(eaglevalue, expansionfactor):
    return eaglevalue * np.sqrt(expansionfactor)

## This private helper function converts a current density value from EAGLE snapshot units to g m^-3
def togcm3(eaglevalue, hubbleparam, expansionfactor):
    return eaglevalue * (6.7699e-31 * (expansionfactor**-3) * (hubbleparam**2))

## This private helper function returns the age of a star (in yr) given the universe expansion factor
# when the star was born (in range 0-1)
def age(R):
    H0 = 2.3e-18
    OmegaM0 = 0.27
    yr = 365.25 * 24 * 3600
    T0 = 13.7e9
    return T0 - (2./3./H0/np.sqrt(1-OmegaM0)) * np.arcsinh(np.sqrt( (1/OmegaM0-1)*R**3 )) / yr

## This private helper function returns the periodicity corrected coordinates input as a (N,3)
# numpy array, and takes the box size (in units of crds) and a test length in units of box size
def periodicCorrec(crds, boxsize, testfact = 0.5):
    for i in range(3):
        crd = crds[:,i]
        booldx = np.abs(crd - crd.min()) > boxsize * testfact
        if booldx.any():
            crd[booldx] = crd[booldx]  - boxsize
    return crds

# This private helper function returns the centre of mass or the centre of mass and mean velocity
# from input particle data, in the units of crds and vels respectively
def centroid(crds, masses, vels = []):
    moments = (crds * np.column_stack([masses]*3)).sum(axis = 0)
    M = masses.sum()
    if vels.any():
        momenta = (vels * np.column_stack([masses]*3)).sum(axis = 0)
        return moments/M, momenta/M
    else:
        return moments/M

# This private helper function returns the periodicity corrected coordinates input as a (N,3)
# numpy array, and takes the box size (in units of crds) and a test length in units of box size.
# Finds the centroid of a given set of particles, each time reducing the
# maximum distance between the previously found centroid to particles used
# to calculate the next centroid.
def shrinkingCentroid(crds, masses, vels, thresh=200, shrinkfactor=1.2):
    N = np.inf       # N set high initially to consider all particles
    while N >= thresh:
        C, C_vel = centroid(crds, masses, vels)

        # define new aperture size as eps
        eps = ((((crds-C)**2).sum(axis = -1))**0.5).max()/float(shrinkfactor)

        # index for particles within new aperture size
        shiftcrds = crds - C
        boolidx = ((shiftcrds**2).sum(axis = -1))**0.5 < eps

        N = boolidx.sum()
        crds = crds[boolidx].copy()
        masses = masses[boolidx].copy()
        vels = vels[boolidx].copy()

    return C, C_vel

# This private helper function returns the unit vector pointing in the direction of the rotation
# axis for input particle data, input CoM and input mean velocity. apt specifies an aperture to
# consider particles within in the units of pos, and aptfrac defines and inner radius within which to
# exclude particles in units of apt.
def rotAxis(crds, vels, mass, com, v_bar, apt = 3e4, aptfrac = 0.08):
    # put in centre of mass and rest frame
    pos = crds - com
    v_rel = vels - v_bar

    # calculate apertures
    disp = (pos**2).sum(axis=-1)
    outapt = disp  < apt ** 2
    inapt = disp  > (aptfrac * apt) ** 2
    totapt = inapt * outapt

    # calculate J vectors in arbitrary units
    Js = np.cross(pos[totapt], v_rel[totapt]) * np.column_stack([mass[totapt]]*3)

    # calculate net J vector and normalise to unit vector
    J = Js.sum(axis = 0)
    n_vect = J * (np.dot(J, J).sum()) ** -0.5
    return n_vect

# -----------------------------------------------------------------
