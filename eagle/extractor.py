#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.extractor Extracting data for a given EAGLE galaxy from the simulation snapshot.
#
# The EAGLE simulation output is stored in a (large) set of data files in the HDF5 format, documented at the
# <a href="http://www.hdfgroup.org/HDF5/">HFD5 home page</a>. The output is organized in \em snapshots, where
# each snapshot represents the state of the universe at a particular time (or equivalently, redshift).
#
# The function in this module allows extracting information relevant for SKIRT from the EAGLE output.
# The function converts physical quantities from EAGLE snapshot units (documented through hdf5 attributes
# in the snapshot files) to SKIRT import units (documented in the SKIRT SPH classes).
# The following table lists some of the units in each system.
#
#<TABLE>
#<TR><TD><B>Physical Quantity</B></TD>  <TD><B>EAGLE snapshot</B></TD>
#                                       <TD><B>SKIRT import</B></TD></TR>
#<TR><TD>position, size</TD>            <TD>\f$\textrm{Mpc}\,a\,h^{-1}\f$</TD>
#                                       <TD>\f$\textrm{pc}\f$</TD></TR>
#<TR><TD>mass</TD>                      <TD>\f$10^{10}\,\textrm{M}_\odot\,h^{-1}\f$</TD>
#                                       <TD>\f$\textrm{M}_\odot\f$</TD></TR>
#<TR><TD>velocity</TD>                  <TD>\f$\textrm{km/s}\,a^{1/2}\f$</TD>
#                                       <TD>--</TD></TR>
#<TR><TD>time, age</TD>                 <TD>--</TD>
#                                       <TD>year</TD></TR>
#<TR><TD>temperature</TD>               <TD>K</TD>
#                                       <TD>--</TD></TR>
#</TABLE>
#
# Note the corrections for cosmological scale factor \f$a\f$ and hubble parameter \f$h\f$ in the EAGLE snapshot units.

# -----------------------------------------------------------------

import os.path
import numpy as np
import h5py
import read_eagle       # EAGLE-specific package by must be seperately installed

from ..core.tools.geometry import Transform
from . import config as config
from .skirtrun import SkirtRun

# -----------------------------------------------------------------

## This function extracts information relevant for SKIRT from the EAGLE output for the galaxy described
# by the specified SKIRT-runs database record. It places the resulting files in the "in" folder of the
# appropriate SkirtRun data structure. The function uses the following fields in the specified record:
# runid, eaglesim, snaptag, galaxyid, groupnr, subgroupnr, copx, copy, copz.
#
# The exported files are named "SIM_GID_stars.dat", "SIM_GID_hii.dat", and "SIM_GID_gas.dat", where
# SIM and GID are replaced respectively by the name of the simulation in which the galaxy resides and by the
# identifier of the galaxy in the public EAGLE database. The file format is as described for SKIRT SPH import.
# In addition, the function creates a text file named "SIM_GID_info.txt", which contains relevant statistics
# including particle numbers and various total masses. The contents is documented in the file.
def extract(record):

    # ---- get the particle data

    # initialise star and gas dictionaries
    sdat        = {}
    gdat        = {}
    yngstars    = {}
    hiiregions  = {}

    # open snapshot and read relevant field attributes
    sfn = snapfilename(record["eaglesim"], record["snaptag"])
    snapshot = read_eagle.EagleSnapshot(sfn)
    params = fieldAttrs(sfn, "Header")
    params.update(fieldAttrs(sfn, "Constants"))
    params.update(fieldAttrs(sfn, "RuntimePars"))
    hubbleparam = params["HubbleParam"]
    expansionfactor = params["ExpansionFactor"]
    schmidtparams = schmidtParameters(params)

    # convert center of potential to snapshot units
    copx = record["copx"] * hubbleparam
    copy = record["copy"] * hubbleparam
    copz = record["copz"] * hubbleparam

    # specify (2*250kpc)^3 physical volume about galaxy centre
    delta = 0.25 * hubbleparam / expansionfactor
    snapshot.select_region(copx-delta, copx+delta, copy-delta, copy+delta, copz-delta, copz+delta)

    # read star particle informaton
    insubhalo = (snapshot.read_dataset(4, "GroupNumber") == record["groupnr"]) & \
                (snapshot.read_dataset(4, "SubGroupNumber") == record["subgroupnr"])
    sdat['r']        = snapshot.read_dataset(4, "Coordinates") [insubhalo]
    sdat['h']        = snapshot.read_dataset(4, "SmoothingLength") [insubhalo]
    sdat['im']       = snapshot.read_dataset(4, "InitialMass") [insubhalo]
    sdat['m']        = snapshot.read_dataset(4, "Mass") [insubhalo]
    sdat['v']        = snapshot.read_dataset(4, "Velocity") [insubhalo]
    sdat['Z']        = snapshot.read_dataset(4, "SmoothedMetallicity") [insubhalo]
    sdat['born']     = snapshot.read_dataset(4, "StellarFormationTime") [insubhalo]
    sdat['rho_born'] = snapshot.read_dataset(4, "BirthDensity") [insubhalo]

    # read gas particle informaton
    insubhalo = (snapshot.read_dataset(0, "GroupNumber") == record["groupnr"]) & \
                (snapshot.read_dataset(0, "SubGroupNumber") == record["subgroupnr"])
    gdat['r']        = snapshot.read_dataset(0, "Coordinates") [insubhalo]
    gdat['h']        = snapshot.read_dataset(0, "SmoothingLength") [insubhalo]
    gdat['m']        = snapshot.read_dataset(0, "Mass") [insubhalo]
    gdat['v']        = snapshot.read_dataset(0, "Velocity") [insubhalo]
    gdat['Z']        = snapshot.read_dataset(0, "SmoothedMetallicity") [insubhalo]
    gdat['T']        = snapshot.read_dataset(0, "Temperature") [insubhalo]
    gdat['rho']      = snapshot.read_dataset(0, "Density") [insubhalo]
    gdat['sfr']      = snapshot.read_dataset(0, "StarFormationRate") [insubhalo]

    # convert units
    sdat['r']        = periodicCorrec(sdat['r'], params["BoxSize"])
    sdat['r']        = toparsec(sdat['r'], hubbleparam, expansionfactor)
    sdat['h']        = toparsec(sdat['h'], hubbleparam, expansionfactor)
    sdat['im']       = tosolar(sdat['im'], hubbleparam)
    sdat['m']        = tosolar(sdat['m'], hubbleparam)
    sdat['t']        = age(sdat['born']) - age(expansionfactor)
    sdat['rho_born'] *= 6.7699e-31
    gdat['r']        = periodicCorrec(gdat['r'], params["BoxSize"])
    gdat['r']        = toparsec(gdat['r'], hubbleparam, expansionfactor)
    gdat['h']        = toparsec(gdat['h'], hubbleparam, expansionfactor)
    gdat['m']        = tosolar(gdat['m'], hubbleparam)
    gdat['rho']      = togcm3(gdat['rho'], hubbleparam, expansionfactor)

    # remember density conversion from g cm^-3 to M_sun Mpc^-3
    densconv = ((params['CM_PER_MPC']/1.e6)**3) / params['SOLAR_MASS']

    # calculate the ISM pressure
    sdat['P']        = getPtot(sdat['rho_born'], schmidtparams)
    gdat['P']        = getPtot(gdat['rho'], schmidtparams)

    # calculate stellar center of mass and translational velocity using shrinking aperture technique
    com, v_bar = shrinkingCentroid(sdat['r'], sdat['m'], sdat['v'])

    # find unit rotation axis vector, using only stellar information and an aperture of 30 kpc
    n_rot = rotAxis(sdat['r'], sdat['v'], sdat['m'], com, v_bar, apt=30e3, aptfrac=0.08)

    # translate to center of mass and line up with angular momentum vector
    transf = Transform()
    transf.translate(-com[0], -com[1], -com[2])
    a, b, c = n_rot
    v = np.sqrt(b*b+c*c)
    if v > 0.3:
        transf.rotateX(c/v, -b/v)
        transf.rotateY(v, -a)
    else:
        v = np.sqrt(a*a+c*c)
        transf.rotateY(c/v, -a/v)
        transf.rotateX(v, -b)
    sdat['r'],w = transf.transform_vec(sdat['r'][:,0],sdat['r'][:,1],sdat['r'][:,2], np.ones(sdat['r'].shape[0]))
    gdat['r'],w = transf.transform_vec(gdat['r'][:,0],gdat['r'][:,1],gdat['r'][:,2], np.ones(gdat['r'].shape[0]))

    # apply 30kpc aperture (i.e. remove all particles outside the aperture)
    applyAperture(sdat, 30e3)
    applyAperture(gdat, 30e3)

    # ---- gather statistics about the data as read from the snapshot

    # information identifying the SKIRT-run record and the galaxy
    info = { }
    info["skirt_run_id"] = record["runid"]
    info["galaxy_id"] = record["galaxyid"]

    # information about the particles
    info["original_particles_stars"] = len(sdat['m'])
    info["original_initial_mass_stars"] = sdat['im'].sum()
    info["original_mass_stars"] = sdat['m'].sum()
    info["original_particles_gas"] = len(gdat['m'])
    info["original_mass_gas"] = gdat['m'].sum()
    info["original_mass_baryons"] = info["original_mass_stars"] + info["original_mass_gas"]

    # information about the direction of the stellar angular momentum axis
    info["original_rotation_axis_x"] = n_rot[0]
    info["original_rotation_axis_y"] = n_rot[1]
    info["original_rotation_axis_z"] = n_rot[2]

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

    info["exported_particles_negative_gas_from_stars"] = 0
    info["exported_particles_negative_gas_from_gas"] = 0
    info["exported_mass_negative_gas_from_stars"] = 0
    info["exported_mass_negative_gas_from_gas"] = 0

    info["exported_particles_unspent_gas_from_gas"] = 0
    info["exported_mass_unspent_gas_from_gas"] = 0

    # ---- resample star forming regions

    # set the "standard" constant covering fraction (see Camps+ 2016)
    f_PDR = 0.1

    # seed the random generator so that a consistent pseudo-random sequence is used for each particular galaxy
    np.random.seed(int(record["galaxyid"]))

    # define HII region age constants (in years)
    young_age = 1e8     # 100 Myr  --> particles below this age are resampled
    infant_age = 1e7    # 10 Myr   --> resampled particles below this age are converted to HII regions
                        #              resampled particles above this age are converted young stars
                        #              <==> lifetime of an HII region

    # set up GALAXEV array
    bcstars = np.column_stack([[],[],[],[],[],[],[]])

    # set up MAPPINGS-III array
    mapstars = np.column_stack([[],[],[],[],[],[],[],[],[]])

    # set up dust array
    dust = np.column_stack([[],[],[],[],[],[],[]])

    # index for particles to resample
    issf = gdat['sfr'] > 0.
    isyoung = sdat['t'] < young_age

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
        sdat['sfr']       = getSFR(sdat['rho_born'], sdat['im'], schmidtparams)

        ms, ts, idxs, mdiffs = stochResamp(sdat['sfr'], sdat['im'])
        isinfant = ts < infant_age

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
            hiiregions['SFR']   = ms[isinfant] / infant_age          # Assume constant SFR over HII region lifetime
            hiiregions['Z']     = sdat['Z'][idxs][isinfant]
            hiiregions['P']     = sdat['P'][idxs][isinfant] * 0.1   # Convert to Pa for output
            hiiregions['logC']  = 0.6*np.log10(ms[isinfant]) + 0.4*np.log10(hiiregions['P']) - 0.4*np.log10(params['BOLTZMANN']) + 0.4
            hiiregions['fPDR']  = np.zeros_like(ts[isinfant]) + f_PDR  # Covering fraction is set to constant value

            # calculate the HII region smoothing length from the mass of the surrounding PDR region,
            # estimated to be 10 times as massive (see Jonsson et al. 2010, MNRAS 403, 17-44),
            # using SKIRT's standard smoothing kernel mass/size normalization: rho = 8/pi * M/h^3;
            # and randomly shift the positions of the HII regions within a similarly enlarged range
            hiiregions['h_mapp'] = (10*ms[isinfant] / (np.pi/8 * sdat['rho_born'][idxs][isinfant] * densconv))**(1/3.)
            stochShiftPos(hiiregions['r'], hiiregions['h'], hiiregions['h_mapp'])

            # append to MAPPINGSIII array
            mapstars = np.concatenate((mapstars, np.column_stack([hiiregions['r'], hiiregions['h_mapp'], hiiregions['SFR'],
                                                                  hiiregions['Z'], hiiregions['logC'], hiiregions['P'],
                                                                  hiiregions['fPDR']])), axis=0)
            info["exported_particles_hii_regions_from_stars"] = np.count_nonzero(isinfant)
            info["exported_initial_mass_hii_regions_from_stars"] = ms[isinfant].sum()
            info["exported_mass_hii_regions_from_stars"] = info["exported_initial_mass_hii_regions_from_stars"]

            # append to dust array with negative mass to compensate for the mass of the surrounding PDR region,
            # considered to be 10 times as massive; use zero temperature as T is unavailable for resampled star particles
            dust = np.concatenate((dust, np.column_stack([hiiregions['r'], hiiregions['h_mapp']*3.,
                                                         -10*ms[isinfant], hiiregions['Z'],
                                                         np.zeros(hiiregions['Z'].shape[0])]).copy()), axis=0)
            info["exported_particles_negative_gas_from_stars"] = np.count_nonzero(isinfant)
            info["exported_mass_negative_gas_from_stars"] = 10*ms[isinfant].sum()

        # add unspent young star particle material to dust array
        # use zero temperature as T is unavailable for resampled star particles
        mass = sdat['im'] - mdiffs
        dust = np.concatenate((dust, np.column_stack([sdat['r'], sdat['h'], mass, sdat['Z'], np.zeros(sdat['Z'].shape[0])]).copy()), axis=0)
        info["exported_particles_unspent_gas_from_stars"] = len(mass)
        info["exported_mass_unspent_gas_from_stars"] = mass.sum()

    # resample gas
    if issf.any():
        for k in gdat.keys():
            gdat[k] = gdat[k][issf].copy()

        ms, ts, idxs, mdiffs = stochResamp(gdat['sfr'], gdat['m'])
        isinfant = ts < infant_age

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
            hiiregions['SFR']   = ms[isinfant] / infant_age          # Assume constant SFR over HII region lifetime
            hiiregions['Z']     = gdat['Z'][idxs][isinfant]
            hiiregions['P']     = gdat['P'][idxs][isinfant] * 0.1     # convert to Pa
            hiiregions['logC']  = 0.6*np.log10(ms[isinfant]) + 0.4*np.log10(hiiregions['P']) - 0.4*np.log10(params['BOLTZMANN']) + 0.4
            hiiregions['fPDR']  = np.zeros_like(ts[isinfant]) + f_PDR  # Covering fraction is set to constant value

            # calculate the HII region smoothing length from the mass of the surrounding PDR region,
            # estimated to be 10 times as massive (see Jonsson et al. 2010, MNRAS 403, 17-44),
            # using SKIRT's standard smoothing kernel mass/size normalization: rho = 8/pi * M/h^3;
            # and randomly shift the positions of the HII regions within a similarly enlarged range
            hiiregions['h_mapp'] = (10*ms[isinfant] / (np.pi/8 * gdat['rho'][idxs][isinfant] * densconv))**(1/3.)
            stochShiftPos(hiiregions['r'], hiiregions['h'], hiiregions['h_mapp'])

            # append to MAPPINGSIII array
            mapstars = np.concatenate((mapstars, np.column_stack([hiiregions['r'], hiiregions['h_mapp'], hiiregions['SFR'],
                                                                  hiiregions['Z'], hiiregions['logC'], hiiregions['P'],
                                                                  hiiregions['fPDR']])), axis=0)
            info["exported_particles_hii_regions_from_gas"] = np.count_nonzero(isinfant)
            info["exported_initial_mass_hii_regions_from_gas"] = ms[isinfant].sum()
            info["exported_mass_hii_regions_from_gas"] = info["exported_initial_mass_hii_regions_from_gas"]

            # append to dust array with negative mass to compensate for the mass of the surrounding PDR region,
            # considered to be 10 times as massive; use negative temperature to indicate that it is not a physical value
            dust = np.concatenate((dust, np.column_stack([hiiregions['r'], hiiregions['h_mapp']*3,
                                                         -10*ms[isinfant], hiiregions['Z'], -gdat['T'][idxs][isinfant]]).copy()), axis=0)
            info["exported_particles_negative_gas_from_gas"] = np.count_nonzero(isinfant)
            info["exported_mass_negative_gas_from_gas"] = 10*ms[isinfant].sum()

        # add unspent SF gas material to dust array; use negative temperature to indicate that it is not a physical value
        mass = gdat['m'] - mdiffs
        dust = np.concatenate((dust, np.column_stack([gdat['r'], gdat['h'], mass, gdat['Z'], -gdat['T']]).copy()), axis=0)
        info["exported_particles_unspent_gas_from_gas"] = len(mass)
        info["exported_mass_unspent_gas_from_gas"] = mass.sum()

    # ---- make some sums and write the statistics and output files

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

    info["exported_particles_negative_gas"] = info["exported_particles_negative_gas_from_stars"] + info["exported_particles_negative_gas_from_gas"]
    info["exported_mass_negative_gas"] = info["exported_mass_negative_gas_from_stars"] + info["exported_mass_negative_gas_from_gas"]

    info["exported_particles_gas"] = info["exported_particles_non_star_forming_gas"] + info["exported_particles_unspent_gas"] + info["exported_particles_negative_gas"]
    info["exported_mass_gas"] = info["exported_mass_non_star_forming_gas"] + info["exported_mass_unspent_gas"] # - info["exported_mass_negative_gas"]
    info["exported_mass_baryons"] = info["exported_mass_stars"] + info["exported_mass_hii_regions"] + info["exported_mass_gas"]

    # create the appropriate SKIRT-run directories
    skirtrun = SkirtRun(record["runid"], create=True)
    filepathprefix = os.path.join(skirtrun.inpath(), "{}_{}_".format(record["eaglesim"], record["galaxyid"]))

    # write the statistics file
    infofile = open(filepathprefix + "info.txt", 'w')
    infofile.write('# Statistics for SPH particles extracted from EAGLE HDF5 snapshot to SKIRT6 format\n')
    infofile.write('# Masses are expressed in solar mass units\n')
    maxkeylen = max(map(len,info.keys()))
    for key in sorted(info.keys()):
        valueformat = "d" if "_particles_" in key or "_id" in key else ".9e"
        infofile.write( ("{0:"+str(maxkeylen)+"} = {1:15"+valueformat+"}\n").format(key, info[key]) )
    infofile.close()

    # ---- write output files

    # open output files
    starsfile = open(filepathprefix + "stars.dat", 'w')
    starsfile.write('# SPH Star Particles\n')
    starsfile.write('# Extracted from EAGLE HDF5 snapshot to SKIRT6 format\n')
    starsfile.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) t(yr)\n')
    gasfile = open(filepathprefix + "gas.dat", 'w')
    gasfile.write('# SPH Gas Particles\n')
    gasfile.write('# Extracted from EAGLE HDF5 snapshot to SKIRT6 format\n')
    gasfile.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) T(K)\n')
    hiifile = open(filepathprefix + "hii.dat", 'w')
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

## This private helper function returns the absolute path to the first EAGLE snapshot file
# corresponding to the given EAGLE simulation name and snapshot tag
def snapfilename(eaglesim, snaptag):

    # snapshot filename segment corresponding to the snapshot tag
    snapname = { 0 : "000_z020p000",
                 1 : "001_z015p132",
                 2 : "002_z009p993",
                 3 : "003_z008p988",
                 4 : "004_z008p075",
                 5 : "005_z007p050",
                 6 : "006_z005p971",
                 7 : "007_z005p487",
                 8 : "008_z005p037",
                 9 : "009_z004p485",
                10 : "010_z003p984",
                11 : "011_z003p528",
                12 : "012_z003p017",
                13 : "013_z002p478",
                14 : "014_z002p237",
                15 : "015_z002p012",
                16 : "016_z001p737",
                17 : "017_z001p487",
                18 : "018_z001p259",
                19 : "019_z001p004",
                20 : "020_z000p865",
                21 : "021_z000p736",
                22 : "022_z000p615",
                23 : "023_z000p503",
                24 : "024_z000p366",
                25 : "025_z000p271",
                26 : "026_z000p183",
                27 : "027_z000p101",
                28 : "028_z000p000" } [snaptag]

    return os.path.join(config.eagledata_path[eaglesim],
                        "particledata_{0}/eagle_subfind_particles_{0}.0.hdf5".format(snapname))

# -----------------------------------------------------------------

## This private helper function reads a hdf5 file field's attributes into a python dictionary.
def fieldAttrs(filename, fieldname):
    fileobj = h5py.File(filename, 'r')
    fieldobj = fileobj[fieldname]
    fieldkeys = list(fieldobj.attrs)
    result = { }
    for key in fieldkeys:
        result[key] = fieldobj.attrs[str(key)]
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

# -----------------------------------------------------------------

## This private helper function returns the periodicity corrected coordinates input as a (N,3)
# numpy array, and takes the box size (in units of crds) and a test length in units of box size
def periodicCorrec(crds, boxsize, testfact = 0.5):
    if len(crds)>0:
        for i in range(3):
            crd = crds[:,i]
            booldx = np.abs(crd - crd.min()) > boxsize * testfact
            if booldx.any():
                crd[booldx] = crd[booldx] - boxsize
    return crds

# -----------------------------------------------------------------

## This private helper function returns the centre of mass or the centre of mass and mean velocity
# from input particle data, in the units of crds and vels respectively
def centroid(crds, masses, vels):
    moments = (crds * np.column_stack([masses]*3)).sum(axis = 0)
    M = masses.sum()
    if vels.any():
        momenta = (vels * np.column_stack([masses]*3)).sum(axis = 0)
        return moments/M, momenta/M
    else:
        return moments/M

# -----------------------------------------------------------------

## This private helper function returns the periodicity corrected coordinates input as a (N,3)
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

# -----------------------------------------------------------------

## This private helper function returns the unit vector pointing in the direction of the rotation
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
    norm2 = np.dot(J, J).sum()
    if norm2 > 0: return J * norm2 ** -0.5
    else: return np.array((0,0,1))

# -----------------------------------------------------------------

## This private helper function applies a spherical aperture to a dictionary of particle data, i.e. it
# adjusts the dictionary so that the particles outside the aperture are removed from each array.
def applyAperture(data, radius):
    x,y,z = data['r'].T
    inside = (x*x+y*y+z*z) <= (radius*radius)
    if inside.any():
        for key in data:
            data[key] = data[key][inside]
    else:
        for key in data:
            shape = list(data[key].shape)
            shape[0] = 0
            data[key] = np.zeros(shape)

# -----------------------------------------------------------------

## This private helper function reads the Schmidt parameters into a python structure.
def schmidtParameters(params):

    # extract relevent unit conversions
    CM_PER_MPC   = params['CM_PER_MPC']
    GAMMA        = params['GAMMA']
    GRAV         = params['GRAVITY']
    K_B          = params['BOLTZMANN']
    M_PROTON     = params['PROTONMASS']
    M_SUN        = params['SOLAR_MASS']
    SEC_PER_YEAR = params['SEC_PER_YEAR']

     # extract relevent runtime parameters used to create EAGLE snapshot
    GammaEff     = params['EOS_Jeans_GammaEffective']
    InitH        = params['InitAbundance_Hydrogen']
    RhoHi        = params['SF_SchmidtLawHighDensThresh_HpCM3']
    RhoNorm      = params['EOS_NormPhysDens_HpCM3']
    SchmidtCoeff = params['SF_SchmidtLawCoeff_MSUNpYRpKPC2']
    SchmidtExp   = params['SF_SchmidtLawExponent']
    SchmidtExpHi = params['SF_SchmidtLawHighDensExponent']
    T_JeansNorm  = params['EOS_Jeans_TempNorm_K']

    # Normalisation in cgs units
    Norm_cgs     = SchmidtCoeff * pow(pow(CM_PER_MPC / 1.e6, 2) / M_SUN , SchmidtExp - 1) / (1.e6 * SEC_PER_YEAR)

    # High density Threshold
    RhoHi_cgs    = RhoHi * M_PROTON / InitH

    # Density normalisation in cgs
    RhoNorm_cgs  = RhoNorm * M_PROTON / InitH

    # Min total Pressure
    P_totc       = RhoNorm * T_JeansNorm * K_B / (InitH * 1.22)

    # Pressure at high density Schmidt law break
    PBreak_cgs   = P_totc * (RhoHi/RhoNorm) ** GammaEff

    # Assume f_g = 1
    NormHi_cgs   = Norm_cgs * (GAMMA * PBreak_cgs / GRAV) ** ((SchmidtExp - SchmidtExpHi) * 0.5)

    # tuple of universal SF parameters
    sfparams     = RhoNorm_cgs, RhoHi_cgs, P_totc, PBreak_cgs, GammaEff

    # tuples of high and low pressure SF parameters
    sf_lo        = Norm_cgs, GAMMA/GRAV, SchmidtExp
    sf_hi        = NormHi_cgs, GAMMA/GRAV, SchmidtExpHi

    return sfparams, sf_lo, sf_hi

# -----------------------------------------------------------------

## This private helper function obtains the SFR of gas from which star particles formed.
#
# Inputs:
#  - rho_form: gas density at formation of star particle
#  - mass: mass of star particle
#  - schmidtpars: parameters for implementing Schmidt law from schmidtParameters()
#
# Outputs:
#  - SFR = Star formation rate for gas particle in input mass units per year
#
def getSFR(rho_form, mass, schmidtpars):

    # unpack universal SF law parameters
    RhoNorm_cgs, RhoHi_cgs, P_totc, PBreak_cgs, GammaEff = schmidtpars[0]

    # Pressure at star formation
    P_form = P_totc * (rho_form / RhoNorm_cgs) ** GammaEff

    # unpack high and low pressure SF law parameters
    sf_lo, sf_hi = schmidtpars[1:]

    # calculate SFR
    if type(rho_form) == np.ndarray:
        hidx = rho_form > RhoHi_cgs
        SFR  = np.zeros(rho_form.size)
        if np.any(hidx):
            SFR[hidx] = mass[hidx] * sf_hi[0] * (sf_hi[1] * P_form[hidx]) ** ((sf_hi[2] - 1) * 0.5)
        if np.any(-hidx):
            SFR[-hidx] = mass[-hidx] * sf_lo[0] * (sf_lo[1] * P_form[-hidx]) ** ((sf_lo[2] - 1) * 0.5)
    else:
        if rho_form > RhoHi_cgs:
            SFR = mass * sf_hi[0] * (sf_hi[1] * P_form) ** ((sf_hi[2] - 1) * 0.5)
        else:
            SFR = mass * sf_lo[0] * (sf_lo[1] * P_form) ** ((sf_lo[2] - 1) * 0.5)

     # return SFR converted to input mass units per year from per second
    return np.array(SFR) * 3.15569e7

# -----------------------------------------------------------------

## This private helper function obtains the ambient pressure of gas from which star particles formed.
#
# Inputs:
#  - rho: gas density of star forming particle
#  - schmidtpars: parameters for implementing Schmidt law from schmidtParameters()
#
# Outputs:
#  - P_tot: Ambient pressure from polytropic effective EoS (Schaye & Dalla Vecchia (2004))
#
def getPtot(rho, schmidtpars):
    RhoNorm_cgs, RhoHi_cgs, P_totc, PBreak_cgs, GammaEff = schmidtpars[0]
    P_form = P_totc * (rho / RhoNorm_cgs) ** GammaEff
    return P_form

# -----------------------------------------------------------------

## This private helper function samples star forming gas particles into a number of sub-particles.
#
# Inputs:
#  - sfr: star formation rate in solar masses per yr
#  - m_gas: particle mass in solar masses
#
# Outputs:
#  - nested arrays with a list of subparticles for each parent input particle:
#     - ms: sub-particle stellar masses in solar masses
#     - ts: lookback times of sub-particle formation
#     - idxs: index of the sub-particle's parent particle in input array
#  - mdiffs: mass of parent particles locked up in new stars; this can be subtracted from the parent gas
#            particles for mass conservation
#
def stochResamp(sfr, m_gas):

    # mass resampling parameters (see Kennicutt & Evans 2012 section 2.5)
    m_min = 700         # minimum mass of sub-particle in M_solar
    m_max = 1e6         # maximum mass of sub-particle in M_solar
    alpha = 1.8         # exponent of power-law mass function
    alpha1 = 1. - alpha

    # age resampling parameters
    thresh_age = 1e8    # period over which to resample in yr (100 Myr)

    # initialise lists for output
    ms   = [[]]
    ts   = [[]]
    idxs = [[]]
    mdiffs = []

    # for each parent particle, determine the star-forming sub-particles
    for i in range(sfr.size):
        sfri = sfr[i]
        mi = m_gas[i]

        # determine the maximum number of sub-particles based on the minimum sub-particle mass
        N = int(max(1,np.ceil(mi/m_min)))

        # generate random sub-particle masses from a power-law distribution between min and max values
        X = np.random.random(N)
        m = (m_min**alpha1 + X*(m_max**alpha1-m_min**alpha1))**(1./alpha1)

        # limit and normalize the list of sub-particles to the total mass of the parent
        mlim = m[np.cumsum(m)<=mi]
        if len(mlim)<1: mlim = m[:1]
        m = mi/mlim.sum() * mlim
        N = len(m)

        # generate random decay lookback time for each sub-particle
        X = np.random.random(N)               # X in range (0,1]
        t = thresh_age + mi/sfri * np.log(1-X)

        # determine mask for sub-particles that form stars by present day
        issf = t > 0.

        # add star-forming sub-particles to the output lists
        ms.append(m[issf])
        ts.append(t[issf])
        idxs.append([i]*np.count_nonzero(issf))
        mdiffs.append(m[issf].sum())

    # convert sub-particle lists into numpy arrays
    ms     = np.hstack(ms)
    ts     = np.hstack(ts)
    idxs   = np.hstack(idxs).astype(int)
    mdiffs = np.array(mdiffs)

    return ms, ts, idxs, mdiffs

# -----------------------------------------------------------------

## This private helper function randomly shifts the positions of HII region sub-particles
# within the smoothing sphere of their parent.
#
# Arguments:
#  - r: parent positions; updated by this function to the shifted positions
#  - h: the smoothing lengths of the parents
#  - h_mapp: the smoothing lengths of the sub-particles
#
def stochShiftPos(r, h, h_mapp):
    # the offset sampling smoothing length is determined so that in the limit of infinite particles,
    # the light distribution is the same as the parent particle kernel;
    # assuming Gaussian kernels this means h_sampling**2 + h_mapp**2 = h**2.
    h_sampling = np.sqrt(np.maximum(0,h*h - h_mapp*h_mapp))

    # sample the offset from a scaled gaussian that resembles a cubic spline kernel
    # (see the documentation of the SPHDustDistribution class in SKIRT)
    r[:,0] += h_sampling * np.random.normal(scale=0.29, size=h_sampling.shape)
    r[:,1] += h_sampling * np.random.normal(scale=0.29, size=h_sampling.shape)
    r[:,2] += h_sampling * np.random.normal(scale=0.29, size=h_sampling.shape)

# -----------------------------------------------------------------
