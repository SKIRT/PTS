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
import eagle.starformation as sf
from eagle.connection import Connection
from pts.geometry import Transform

# -----------------------------------------------------------------

## This class represents the set of files containing a single EAGLE snapshot (i.e. the state of the universe
# at a particular redshift). The constructor locates all files, verifies that the snapshot is complete,
# and cashes some relevant information. The other functions return or export data extracted from the snapshot.
#
class Snapshot:

    ## The constructor locates all files containing particle data for the snapshot, verifies that the data is
    # consistent and complete, and cashes some relevant information. The arguments are:
    #  - eaglesim: the EAGLE simulation for this snaphot; this value is used as a key into the
    #    \em config.eagledata_path dictionary to obtain the path to the eagle data directory,
    #    the default value is the value of the config.default_eaglesim field;
    #  - redshift: the redshift for this snapshot, with an accuracy of at least 3 digits after the decimal point;
    #    the default value is the value of the config.default_redshift field.
    #
    # The specified redshift is converted into a filename pattern (formatted as zNNNpFFF where NNN is the portion
    # before the decimal point and FFF is the fractional portion), which is used to filter the files in the directory.
    # Any subdirectories with a name matching the redshift pattern are searched as well (recursively).
    def __init__(self, eaglesim=config.default_eaglesim, redshift=config.default_redshift):
        # remember the constructor arguments
        self.eaglesim = eaglesim
        self.orig_redshift = redshift

        # format the redshift pattern
        intportion = int(redshift)
        fracportion = int((redshift-intportion)*1000+0.5)
        self.pattern = "z{0:0>3}p{1:0>3}".format(intportion,fracportion)

        # make a list of all hdf5 file paths that match the pattern
        self.directory = config.eagledata_path[eaglesim]
        paths = findhdf5files(self.directory, self.pattern)
        if len(paths) == 0: raise ValueError("No matching HDF5 files found")

        # open one of the files, arbitrarily, and retrieve the relevant header information
        hdf = h5py.File(paths[0], 'r')
        self.runlabel = hdf["Header"].attrs["RunLabel"]
        self.hubbleparam = hdf["Header"].attrs["HubbleParam"]
        self.boxsize = toparsec(hdf["Header"].attrs["BoxSize"], self.hubbleparam, 1)
        self.redshift = hdf["Header"].attrs["Redshift"]
        self.expansionfactor = hdf["Header"].attrs["ExpansionFactor"]

        # set up attribute dictionaries
        self.constants = sf.readAttrs(hdf, field='Constants')
        self.header = sf.readAttrs(hdf, field='Header')
        self.runpars = sf.readAttrs(hdf, field='RuntimePars')

        # extract star formation law parameters
        self.schmidtpars = sf.schmidtParameters(self.constants, self.runpars)

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
        print "Simulation: " + self.eaglesim
        print "Box size:  {0:.1f} Mpc".format(self.boxsize/1e6)
        print "Redshift:  {0:.3f}".format(self.redshift)
        print "Expansion: {0:.1f}%".format(self.expansionfactor*100)
        print "Files:     {0}".format(self.numfiles)
        print "Star particles: {0:11,}".format(self.numstarparticles())
        print "Gas particles: {0:12,}".format(self.numgasparticles())

    ## This function returns the absolute file path for the catalog corresponding to this snapshot.
    # The catalog resides in the \em config.catalogs_path directory, and its file name contains the
    # simulation code and redshift pattern for this snapshot.
    def catalogfilepath(self):
        return os.path.join(config.catalogs_path, self.eaglesim + "_" + self.pattern + "_catalog.dat")

    ## This function builds and saves a catalog for the snapshot in a file with the path returned by the
    # catalogfilepath() function. The function first queries the public EAGLE database to obtain a list
    # of galaxies that have at least the specified stellar mass within an aperture of 30 kpc. It then
    # builds a catalog including a text line for each of these galaxies with the following columns:
    # galaxy ID, group number, subgroup number, number of star particles, number of gas particles, total stellar mass,
    # total gass mass, list of indices of the hdf5 data files holding the corresponding particle information.
    # The masses are not limited to an aperture and expressed in \f$\textrm{M}_\odot\f$.
    def exportcatalog(self, minstarmass):

        # query the public EAGLE database
        print "Querying public EAGLE database..."
        # (currently we always query the snapshot for redshift zero)
        if self.orig_redshift!=0: raise ValueError("Nonzero redshift not supported at this time")
        con = Connection("camps", password="4nAPEqcs")
        records = con.execute_query('''
            SELECT
                gal.GalaxyID as galaxyid,
                gal.GroupNumber as groupnumber,
                gal.SubGroupNumber as subgroupnumber
            FROM
                {0}_SubHalo as gal,
                {0}_Aperture as ape
            WHERE
                ape.Mass_Star > {1} and
                gal.SnapNum = 28 and
                ape.ApertureSize = 30 and
                gal.GalaxyID = ape.GalaxyID
            '''.format(config.eagledatabase_name[self.eaglesim], minstarmass))
        print "There will be {} galaxies in this catalog".format(len(records))

        # build a dictionary with a key for every reported galaxy and placeholders for the values
        # key = ( group number, subgroup number )
        # value = [ galaxy ID, # star particles, # gas particles, stellar mass, gass mass, set of file indices ]
        print "Building catalog..."
        subgroups = { }
        for record in records:
            key = (int(record['groupnumber']), int(record['subgroupnumber']))
            subgroups[key] = [ int(record['galaxyid']), 0, 0, 0., 0., set() ]

        # loop over all star and gas particles and track group membership, filling in the values of the dictionary
        for index,path in zip(range(len(self.paths)),self.paths):
            if index%10==0: print " opening file", index+1, "of", len(self.paths)
            hdf = h5py.File(path, 'r')
            stars = hdf["PartType4"]
            for groupnumber,subnumber,mass in zip(stars["GroupNumber"][:],stars["SubGroupNumber"][:],stars["Mass"][:]):
                if subnumber < 1073741824:   # (1073741824 means "not in subgroup")
                    key = (abs(groupnumber), subnumber)
                    if key in subgroups:
                        subgroup = subgroups[key]
                        subgroup[1] += 1
                        subgroup[3] += mass
                        subgroup[5].add(index)
            gas = hdf["PartType0"]
            for groupnumber,subnumber,mass in zip(gas["GroupNumber"][:],gas["SubGroupNumber"][:],gas["Mass"][:]):
                if subnumber < 1073741824:   # (1073741824 means "not in subgroup")
                    key = (abs(groupnumber), subnumber)
                    if key in subgroups:
                        subgroup = subgroups[key]
                        subgroup[2] += 1
                        subgroup[4] += mass
                        subgroup[5].add(index)
            hdf.close()
        print " done."

        # save catalog entries
        print "Saving catalog..."
        out = open(self.catalogfilepath(),'w')
        for key in sorted(subgroups.keys()):
            value = subgroups[key]
            line = str(value[0])
            for index in range(0,2): line += ' ' + str(key[index])
            for index in range(1,3): line += ' ' + str(value[index])
            for index in range(3,5): line += ' ' + str(tosolar(value[index],self.hubbleparam))
            for index in sorted([i for i in value[5]]): line += ' ' + str(index)
            out.write(line + '\n')
        out.close()
        print " done."

    ## This function returns a GalaxyList object holding all galaxies in the snapshot listed in the catalog.
    # The function expects that the catalog has already been cosntructed by the exportcatalog() function.
    def galaxies(self):
        # build a galaxy list from the catalog
        galaxylist = GalaxyList()
        for line in open(self.catalogfilepath(),'r'):
            galaxylist.addgalaxy(Galaxy(self, line.split()))
        return galaxylist

# -----------------------------------------------------------------

## This class holds a list of Galaxy objects. It provides facilities to remove galaxies from the list
# based on certain criteria, and to pretty-print information on the remaining items.
class GalaxyList:

    ## This constructor creates an empty galaxy list.
    def __init__(self):
        self.galaxies = { }     # key: galaxy id; value: Galaxy object

    ## This function adds a galaxy object to the list.
    def addgalaxy(self, galaxy):
        self.galaxies[ galaxy.galaxyid ] = galaxy

    ## This function returns the Galaxy object with the specified galaxy id, or None if it is not in the list.
    def galaxy(self, galaxyid):
        return self.galaxies[galaxyid] if galaxyid in self.galaxies else None

    ## This function prints some basic properties for all galaxies currently in the list. The output has
    # one line per galaxy and is sorted on descending total number of particles.
    def printinfo(self):
        if len(self.galaxies) > 0:
            print "There are {0} galaxies in the list".format(len(self.galaxies))
            self.galaxies.itervalues().next().printinfotitle()
            for galaxy in sorted(self.galaxies.values(), key=Galaxy.numparticles, reverse=True):
                galaxy.printinfo()
        else:
            print "There are no galaxies in the list"

# -----------------------------------------------------------------

## This class represents a single galaxy (a particle subgroup) in a snapshot. It serves to extract
# information about the galaxy from the catalog entry and/or from underlying snapshot data files.
class Galaxy:

    ## This constructor initializes all property fields for the galaxy. The first argument provides a reference
    # to the snapshot containing this galaxy. The second argument is a list of field values used to initialize
    # the galaxy properties, corrresponding to the columns of the catalog entry for the galaxy:
    # galaxy id, group number, subgroup number, number of star particles, number of gas particles, total stellar mass,
    # total gass mass, list of indices of the hdf5 data files holding the corresponding particle information.
    # The masses are expressed in \f$\textrm{M}_\odot\f$.
    def __init__(self, snapshot, fields):
        self.snapshot = snapshot
        self.galaxyid, self.groupnumber, self.subgroupnumber, self.numstarparticles, self.numgasparticles \
                = map(int, fields[0:5])
        self.starmass, self.gasmass = map(float, fields[5:7])
        self.fileindices = map(int, fields[7:])
        if len(self.fileindices)==0: raise ValueError("No cross reference info provided")

    ## This function returns the total number of particles in the galaxy
    def numparticles(self):
        return self.numstarparticles + self.numgasparticles

    ## This function returns a formatted string identifying the galaxy which can be used as a filename prefix.
    def prefix(self):
        return "{}_{}".format(self.snapshot.eaglesim, self.galaxyid)

    ## This function prints a header line for the information lines printed by the printinfo() function
    def printinfotitle(self):
        print "galaxyid group subgr  #star part.  #gas part.  star mass  gas mass"
        print "------------------------------------------------------------------"

    ## This function prints some information on the galaxy on a single line
    def printinfo(self):
        print "{:8} {:5} {:4}  {:11,} {:11,}   {:9.1e} {:9.1e}"  \
            .format(self.galaxyid, self.groupnumber, self.subgroupnumber,
                    self.numstarparticles, self.numgasparticles, self.starmass, self.gasmass)

    ## This function exports particle data for the galaxy to three text files that can be imported in SKIRT.
    # The exported files are placed in the specified directory and are named "SIM_GID_stars.dat",
    # "SIM_GID_hii.dat", and "SIM_GID_gas.dat" where SIM and GID are replaced respectively
    # by the name of the snapshot in which the galaxy resides, and the galaxy ID identifying the galaxy
    # in the public EAGLE database. The file format is as described for SKIRT SPH import.
    # In addition, the function creates a text file named "SIM_GID_info.txt", which contains relevant statistics
    # including particle numbers and various total masses. The contents is documented in the file.
    def export(self, directory=""):
        print "Exporting galaxy ({0},{1}) from {2} files...".format(  \
                    self.groupnumber, self.subgroupnumber, len(self.fileindices))

        # ---- get the particle data

        # initialise star and gas dictionaries
        sdat        = {}
        gdat        = {}
        yngstars    = {}
        hiiregions  = {}
        sdat['r']        = np.column_stack([[],[],[]])
        sdat['h']        = np.array([])
        sdat['im']       = np.array([])
        sdat['m']        = np.array([])
        sdat['v']        = np.column_stack([[],[],[]])
        sdat['Z']        = np.array([])
        sdat['born']     = np.array([])
        sdat['rho_born'] = np.array([])
        gdat['r']        = np.column_stack([[],[],[]])
        gdat['h']        = np.array([])
        gdat['m']        = np.array([])
        gdat['v']        = np.column_stack([[],[],[]])
        gdat['Z']        = np.array([])
        gdat['T']        = np.array([])
        gdat['rho']      = np.array([])
        gdat['sfr']      = np.array([])

        # read particle informaton
        for index in self.fileindices:

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
        sdat['m']        = tosolar(sdat['m'], self.snapshot.hubbleparam)
        sdat['Z']        = sdat['Z']
        sdat['t']        = age(sdat['born']) - age(self.snapshot.expansionfactor)
        sdat['rho_born'] *= 6.7699e-31
        gdat['r']        = toparsec(gdat['r'], self.snapshot.hubbleparam, self.snapshot.expansionfactor)
        gdat['r']        = periodicCorrec(gdat['r'], self.snapshot.boxsize)
        gdat['h']        = toparsec(gdat['h'], self.snapshot.hubbleparam, self.snapshot.expansionfactor)
        gdat['m']        = tosolar(gdat['m'], self.snapshot.hubbleparam)
        gdat['Z']        = gdat['Z']
        gdat['T']        = gdat['T']
        gdat['rho']      = togcm3(gdat['rho'], self.snapshot.hubbleparam, self.snapshot.expansionfactor)
        gdat['sfr']      = gdat['sfr']

        # density conversion from g cm^-3 to M_sun Mpc^-3
        densconv = ((self.snapshot.constants['CM_PER_MPC']/1.e6)**3) / self.snapshot.constants['SOLAR_MASS']

        # calculate the ISM pressure
        sdat['P']        = sf.getPtot(sdat['rho_born'], self.snapshot.schmidtpars)
        gdat['P']        = sf.getPtot(gdat['rho'], self.snapshot.schmidtpars)

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

        # apply 30kpc aperture (i.e. remove all particles outside the aperture)
        applyAperture(sdat, 30e3)
        applyAperture(gdat, 30e3)

        # ---- gather statistics about data as it is read from the hdf5 snapshot

        info = { }
        info["original_particles_stars"] = len(sdat['m'])
        info["original_initial_mass_stars"] = sdat['im'].sum()
        info["original_mass_stars"] = sdat['m'].sum()
        info["original_particles_gas"] = len(gdat['m'])
        info["original_mass_gas"] = gdat['m'].sum()
        info["original_mass_baryons"] = info["original_mass_stars"] + info["original_mass_gas"]

        # ---- initialize statistics about the exported data

        info["exported_particles_old_stars"] = 0
        info["exported_initial_mass_old_stars"] = 0
        info["exported_mass_old_stars"] = 0

        info["exported_particles_non_star_forming_gas"] = 0
        info["exported_mass_non_star_forming_gas"] = 0

        info["exported_particles_young_stars_from_stars"] = 0
        info["exported_initial_mass_young_stars_from_stars"] = 0
        info["exported_mass_young_stars_from_stars"] = 0

        info["exported_particles_hii_regions_from_stars"] = 0
        info["exported_initial_mass_hii_regions_from_stars"] = 0
        info["exported_mass_hii_regions_from_stars"] = 0

        info["exported_particles_unspent_gas_from_stars"] = 0
        info["exported_mass_unspent_gas_from_stars"] = 0

        info["exported_particles_young_stars_from_gas"] = 0
        info["exported_initial_mass_young_stars_from_gas"] = 0
        info["exported_mass_young_stars_from_gas"] = 0

        info["exported_particles_hii_regions_from_gas"] = 0
        info["exported_initial_mass_hii_regions_from_gas"] = 0
        info["exported_mass_hii_regions_from_gas"] = 0

        info["exported_particles_unspent_gas_from_gas"] = 0
        info["exported_mass_unspent_gas_from_gas"] = 0

        # ---- resample star forming regions

        # set up GALAXEV array
        bcstars = np.column_stack([[],[],[],[],[],[],[]])

        # set up MAPPINGS-III array
        mapstars = np.column_stack([[],[],[],[],[],[],[],[],[]])

        # set up dust array
        dust = np.column_stack([[],[],[],[],[],[],[]])

        # index for particles to resample
        issf = gdat['sfr'] > 0.
        isyoung = sdat['t'] < 1e8   # 100 Myr

        # append older stars to GALAXEV array
        if (~isyoung).any():
            bcstars = np.concatenate((bcstars, np.column_stack([sdat['r'], sdat['h'], sdat['im'], sdat['Z'], sdat['t']])[~isyoung]), axis=0)
            info["exported_particles_old_stars"] = np.count_nonzero(~isyoung)
            info["exported_initial_mass_old_stars"] = sdat['im'][~isyoung].sum()
            info["exported_mass_old_stars"] = sdat['m'][~isyoung].sum()

        # append non-SF gas data to dust array
        if (~issf).any():
            dust = np.concatenate((dust, np.column_stack([gdat['r'], gdat['h'], gdat['m'], gdat['Z'], gdat['T']])[~issf].copy()), axis=0)
            info["exported_particles_non_star_forming_gas"] = np.count_nonzero(~issf)
            info["exported_mass_non_star_forming_gas"] = gdat['m'][~issf].sum()

        # resample stars
        if isyoung.any():
            for k in sdat.keys():
                sdat[k] = sdat[k][isyoung].copy()

            # calculate SFR at birth of young star particles in M_sun / yr
            sdat['sfr']       = sf.getSFR(sdat['rho_born'], sdat['im'], self.snapshot.schmidtpars)

            ms, ts, idxs, mdiffs = sf.stochResamp(sdat['sfr'], sdat['im'])
            isinfant = ts < 1e7

            if (~isinfant).any():
                yngstars['r']  = sdat['r'][idxs][~isinfant]
                yngstars['h']  = sdat['h'][idxs][~isinfant]
                yngstars['im'] = ms[~isinfant]
                yngstars['Z']  = sdat['Z'][idxs][~isinfant]
                yngstars['t']  = ts[~isinfant]
                bcstars = np.concatenate((bcstars, np.column_stack([yngstars['r'], yngstars['h'], yngstars['im'], yngstars['Z'], yngstars['t']])), axis=0)
                info["exported_particles_young_stars_from_stars"] = np.count_nonzero(~isinfant)
                info["exported_initial_mass_young_stars_from_stars"] = ms[~isinfant].sum()
                info["exported_mass_young_stars_from_stars"] = info["exported_initial_mass_young_stars_from_stars"]

            if (isinfant).any():
                hiiregions['r']     = sdat['r'][idxs][isinfant]
                hiiregions['h']     = sdat['h'][idxs][isinfant]
                hiiregions['SFR']   = ms[isinfant] * 1e-7               # Assume constant SFR over 10 Myr, so sfr =  mcl / 10 Myr
                hiiregions['Z']     = sdat['Z'][idxs][isinfant]
                hiiregions['P']     = sdat['P'][idxs][isinfant] * 0.1   # Convert to Pa for output
                hiiregions['logC']  = 0.6*np.log10(ms[isinfant]) + 0.4*np.log10(hiiregions['P']) - 0.4*np.log10(self.snapshot.constants['BOLTZMANN'])
                hiiregions['h_mapp']  = (ms[isinfant] / (sdat['rho_born'][idxs][isinfant] * densconv))**(1/3.)
                hiiregions['fPDR']  = 1 - (ts[isinfant]/1e7)            # Covering fraction goes from 1 to 0 over HII region lifetime

                # randomly shift the positions of the HII regions
                sf.stochShiftPos(hiiregions['r'], hiiregions['h'], hiiregions['h_mapp'])

                # append to MAPPINGSIII array
                mapstars = np.concatenate((mapstars, np.column_stack([hiiregions['r'], hiiregions['h_mapp'], hiiregions['SFR'],
                                                                      hiiregions['Z'], hiiregions['logC'], hiiregions['P'],
                                                                      hiiregions['fPDR']])), axis=0)
                info["exported_particles_hii_regions_from_stars"] = np.count_nonzero(isinfant)
                info["exported_initial_mass_hii_regions_from_stars"] = ms[isinfant].sum()
                info["exported_mass_hii_regions_from_stars"] = info["exported_initial_mass_hii_regions_from_stars"]

            # T is set to 0 K as T unavailable for resampled star particles
            sdat['T'] = np.zeros(sdat['im'].shape[0])

            # add unspent young star particle material to dust array
            mass = sdat['im'] - mdiffs
            dust = np.concatenate((dust, np.column_stack([sdat['r'], sdat['h'], mass, sdat['Z'], sdat['T']]).copy()), axis=0)
            info["exported_particles_unspent_gas_from_stars"] = len(mass)
            info["exported_mass_unspent_gas_from_stars"] = mass.sum()

        # resample gas
        if issf.any():
            for k in gdat.keys():
                gdat[k] = gdat[k][issf].copy()

            ms, ts, idxs, mdiffs = sf.stochResamp(gdat['sfr'], gdat['m'])
            isinfant = ts < 1.e7

            if (~isinfant).any():
                yngstars['r']  = gdat['r'][idxs][~isinfant]
                yngstars['h']  = gdat['h'][idxs][~isinfant]
                yngstars['im'] = ms[~isinfant]
                yngstars['Z']  = gdat['Z'][idxs][~isinfant]
                yngstars['t']  = ts[~isinfant]
                bcstars = np.concatenate((bcstars, np.column_stack([yngstars['r'], yngstars['h'], yngstars['im'], yngstars['Z'], yngstars['t']])), axis=0)
                info["exported_particles_young_stars_from_gas"] = np.count_nonzero(~isinfant)
                info["exported_initial_mass_young_stars_from_gas"] = ms[~isinfant].sum()
                info["exported_mass_young_stars_from_gas"] = info["exported_initial_mass_young_stars_from_gas"]

            if (isinfant).any():
                hiiregions['r']     = gdat['r'][idxs][isinfant]
                hiiregions['h']     = gdat['h'][idxs][isinfant]
                hiiregions['SFR']   = ms[isinfant] * 1e-7               # Assume constant SFR over 10 Myr, so sfr =  mcl / 10 Myr
                hiiregions['Z']     = gdat['Z'][idxs][isinfant]
                hiiregions['P']     = gdat['P'][idxs][isinfant] * 0.1     # convert to Pa
                hiiregions['logC']  = 0.6*np.log10(ms[isinfant]) + 0.4*np.log10(hiiregions['P']) - 0.4*np.log10(self.snapshot.constants['BOLTZMANN'])
                hiiregions['h_mapp'] = (ms[isinfant] / (gdat['rho'][idxs][isinfant] * densconv))**(1/3.)
                hiiregions['fPDR']  = 1 - (ts[isinfant]/1e7)            # Covering fraction goes from 1 to 0 over HII region lifetime

                # randomly shift the positions of the HII regions
                sf.stochShiftPos(hiiregions['r'], hiiregions['h'], hiiregions['h_mapp'])

                # append to MAPPINGSIII array
                mapstars = np.concatenate((mapstars, np.column_stack([hiiregions['r'], hiiregions['h_mapp'], hiiregions['SFR'],
                                                                      hiiregions['Z'], hiiregions['logC'], hiiregions['P'],
                                                                      hiiregions['fPDR']])), axis=0)
                info["exported_particles_hii_regions_from_gas"] = np.count_nonzero(isinfant)
                info["exported_initial_mass_hii_regions_from_gas"] = ms[isinfant].sum()
                info["exported_mass_hii_regions_from_gas"] = info["exported_initial_mass_hii_regions_from_gas"]

            # add unspent SF gas material to dust array; use negative temperature to indicate that it is not a physical value
            mass = gdat['m'] - mdiffs
            dust = np.concatenate((dust, np.column_stack([gdat['r'], gdat['h'], mass, gdat['Z'], -gdat['T']]).copy()), axis=0)
            info["exported_particles_unspent_gas_from_gas"] = len(mass)
            info["exported_mass_unspent_gas_from_gas"] = mass.sum()

        # ---- make some sums and write statistics in info file

        info["exported_particles_young_stars"] = info["exported_particles_young_stars_from_stars"] + info["exported_particles_young_stars_from_gas"]
        info["exported_initial_mass_young_stars"] = info["exported_initial_mass_young_stars_from_stars"] + info["exported_initial_mass_young_stars_from_gas"]
        info["exported_mass_young_stars"] = info["exported_mass_young_stars_from_stars"] + info["exported_mass_young_stars_from_gas"]

        info["exported_particles_stars"] = info["exported_particles_old_stars"] + info["exported_particles_young_stars"]
        info["exported_initial_mass_stars"] = info["exported_initial_mass_old_stars"] + info["exported_initial_mass_young_stars"]
        info["exported_mass_stars"] = info["exported_mass_old_stars"] + info["exported_mass_young_stars"]

        info["exported_particles_hii_regions"] = info["exported_particles_hii_regions_from_stars"] + info["exported_particles_hii_regions_from_gas"]
        info["exported_initial_mass_hii_regions"] = info["exported_initial_mass_hii_regions_from_stars"] + info["exported_initial_mass_hii_regions_from_gas"]
        info["exported_mass_hii_regions"] = info["exported_mass_hii_regions_from_stars"] + info["exported_mass_hii_regions_from_gas"]

        info["exported_particles_unspent_gas"] = info["exported_particles_unspent_gas_from_stars"] + info["exported_particles_unspent_gas_from_gas"]
        info["exported_mass_unspent_gas"] = info["exported_mass_unspent_gas_from_stars"] + info["exported_mass_unspent_gas_from_gas"]

        info["exported_particles_gas"] = info["exported_particles_non_star_forming_gas"] + info["exported_particles_unspent_gas"]
        info["exported_mass_gas"] = info["exported_mass_non_star_forming_gas"] + info["exported_mass_unspent_gas"]
        info["exported_mass_baryons"] = info["exported_mass_stars"] + info["exported_mass_hii_regions"] + info["exported_mass_gas"]

        infofilename = self.prefix() + "_info.txt"
        infofile = open(os.path.join(directory,infofilename), 'w')
        infofile.write('# Statistics for SPH particles extracted from EAGLE HDF5 snapshot to SKIRT6 format\n')
        infofile.write('# Masses are expressed in solar mass units\n')
        maxkeylen = max(map(len,info.keys()))
        for key in sorted(info.keys()):
            valueformat = "d" if "_particles_" in key else ".9e"
            infofile.write( ("{0:"+str(maxkeylen)+"} = {1:15"+valueformat+"}\n").format(key, info[key]) )
        infofile.close()

        # ---- write output files

        # open output files
        starsfilename = self.prefix() + "_stars.dat"
        gasfilename = self.prefix() + "_gas.dat"
        hiifilename = self.prefix() + "_hii.dat"
        starsfile = open(os.path.join(directory,starsfilename), 'w')
        starsfile.write('# SPH Star Particles\n')
        starsfile.write('# Extracted from EAGLE HDF5 snapshot to SKIRT6 format\n')
        starsfile.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) t(yr)\n')
        gasfile = open(os.path.join(directory,gasfilename), 'w')
        gasfile.write('# SPH Gas Particles\n')
        gasfile.write('# Extracted from EAGLE HDF5 snapshot to SKIRT6 format\n')
        gasfile.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) T(K)\n')
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

## This private helper function recursively locates HDF5 files in a nested directory structure,
# as long as the directory or file name contains the specified pattern string.
def findhdf5files(dirpath, pattern):
    result = []
    for name in os.listdir(dirpath):
        if pattern in name and "particle" in name and not "snip" in name:
            path = os.path.join(dirpath,name);
            if os.path.isdir(path):
                result += findhdf5files(path, pattern)
            if os.path.isfile(path) and name.lower().endswith(".hdf5"):
                result += [ path ]
    return result

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

# This private helper function applies a spherical aperture to a dictionary of particle data, i.e. it
# adjusts the dictionary so that the particles outside the aperture are removed from each array.
def applyAperture(data, radius):
    x,y,z = data['r'].T
    inside = (x*x+y*y+z*z) <= (radius*radius)
    for key in data:
        data[key] = data[key][inside]

# -----------------------------------------------------------------
