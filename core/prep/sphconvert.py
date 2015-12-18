#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.sphconvert Converting SPH output data to SKIRT input format.
#
# The functions in this module allow converting SPH data files in text column
# format to the SKIRT input format. Currently supported are:
#  - EAGLE old column text format (compatible with SKIRT5)
#  - AWAT column text format
#  - DOLAG column text format
#
# There is a separate function for star and gas particles, for each format.
# The arguments for each function are:
#  - infile: the name of the input file in foreign format
#  - outfile: the name of the output file in SKIRT6 format (file is overwritten)

# -----------------------------------------------------------------

# Import standard modules
import math as math
import numpy as np

# -----------------------------------------------------------------
#  EAGLE column text format
# -----------------------------------------------------------------

## EAGLE star particles:
# - incoming:  x(kpc) y(kpc) z(kpc) t(yr) h(kpc) Z(0-1) M(Msun)
# - outgoing:  x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) t(yr)
def convert_stars_EAGLE(infile, outfile):
    x,y,z,t,h,Z,M = np.loadtxt(infile, unpack=True)
    fid = open(outfile, 'w')
    fid.write('# SPH Star Particles\n')
    fid.write('# Converted from EAGLE SKIRT5 output format into SKIRT6 format\n')
    fid.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) t(yr)\n')
    np.savetxt(fid, np.transpose((x*1e3,y*1e3,z*1e3,h*1e3,M,Z,t)), fmt="%1.9g")
    fid.close()

## EAGLE gas particles:
# - incoming:  x(kpc) y(kpc) z(kpc) SFR(?) h(kpc) Z(0-1) M(Msun)
# - outgoing:  x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1)
def convert_gas_EAGLE(infile, outfile):
    x,y,z,SFR,h,Z,M = np.loadtxt(infile, unpack=True)
    fid = open(outfile, 'w')
    fid.write('# SPH Gas Particles\n')
    fid.write('# Converted from EAGLE SKIRT5 output format into SKIRT6 format\n')
    fid.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1)\n')
    np.savetxt(fid, np.transpose((x*1e3,y*1e3,z*1e3,h*1e3,M,Z)), fmt="%1.9g")
    fid.close()

# -----------------------------------------------------------------
#  AWAT column text format
# -----------------------------------------------------------------

## AWAT star particles:
# - incoming:  x y z vx vy vz M ms0 mzHe mzC mzN mzO mzNe mzMg mzSi mzFe mzZ Z ts id flagfd rho h ...
# -    units:  x,y,z,h (100kpc); M (1e12 Msun); ts(0.471Gyr) with t = (1Gyr-ts)
# - outgoing:  x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) t(yr)
def convert_stars_AWAT(infile, outfile):
    x,y,z,M,Z,ts,h = np.loadtxt(infile, usecols=(0,1,2,6,17,18,22), unpack=True)
    fid = open(outfile, 'w')
    fid.write('# SPH Star Particles\n')
    fid.write('# Converted from AWAT output format into SKIRT6 format\n')
    fid.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) t(yr)\n')
    np.savetxt(fid, np.transpose((x*1e5,y*1e5,z*1e5,h*1e5,M*1e12,Z,1e9-ts*0.471e9)), fmt="%1.9g")
    fid.close()

## AWAT gas particles:
# - incoming:  x y z vx vy vz M rho u mzHe mzC mzN mzO mzNe mzMg mzSi mzFe mzZ id flagfd h myu nhp Temp ...
# -    units:  x,y,z,h (100kpc); M (1e12 Msun); mzZ (Msun) so that Z=mzZ/(M*1e12)
# - outgoing:  x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1)
def convert_gas_AWAT(infile, outfile):
    x,y,z,M,mzZ,h = np.loadtxt(infile, usecols=(0,1,2,6,17,20), unpack=True)
    fid = open(outfile, 'w')
    fid.write('# SPH Gas Particles\n')
    fid.write('# Converted from AWAT output format into SKIRT6 format\n')
    fid.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1)\n')
    np.savetxt(fid, np.transpose((x*1e5,y*1e5,z*1e5,h*1e5,M*1e12,mzZ/(M*1e12))), fmt="%1.9g")
    fid.close()

# -----------------------------------------------------------------
#  DOLAG column text format
# -----------------------------------------------------------------

# return the age of a star (in yr) given the universe expansion factor when the star was born (in range 0-1)
@np.vectorize
def age(R):
    H0 = 2.3e-18
    OmegaM0 = 0.27
    yr = 365.25 * 24 * 3600
    T0 = 13.7e9
    return T0 - (2./3./H0/np.sqrt(1-OmegaM0)) * np.arcsinh(np.sqrt( (1/OmegaM0-1)*R**3 )) / yr

# return the radius of a particle (in kpc) given its mass (in Msun) and density (in Msun/kpc3)
@np.vectorize
def radius(M,rho):
    return (M/rho*3/4/math.pi*64)**(1./3.)

## DOLAG star particles:
# - incoming:  id x y z vx vy vz M R
# -    units:  x,y,z (kpc); M (Msun); R (0-1); assume Z=0.02 & h=1kpc; calculate t(R)
# - outgoing:  x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) t(yr)
def convert_stars_DOLAG(infile, outfile):
    x,y,z,M,R = np.loadtxt(infile, usecols=(1,2,3,7,8), unpack=True)
    fid = open(outfile, 'w')
    fid.write('# SPH Star Particles\n')
    fid.write('# Converted from DOLAG output format into SKIRT6 format\n')
    fid.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1) t(yr)\n')
    np.savetxt(fid, np.transpose((x*1e3,y*1e3,z*1e3,np.ones_like(x)*1e3,M,np.ones_like(x)*0.02,age(R))), fmt="%1.9g")
    fid.close()

## DOLAG gas particles:
# - incoming:  id x y z vx vy vz M rho T cf u sfr
# -    units:  x,y,z (kpc); M (Msun); assume Z=0.02; calculate h(M,rho)
# - outgoing:  x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1)
def convert_gas_DOLAG(infile, outfile):
    x,y,z,M,rho = np.loadtxt(infile, usecols=(1,2,3,7,8), unpack=True)
    fid = open(outfile, 'w')
    fid.write('# SPH Gas Particles\n')
    fid.write('# Converted from DOLAG output format into SKIRT6 format\n')
    fid.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1)\n')
    np.savetxt(fid, np.transpose((x*1e3,y*1e3,z*1e3,radius(M,rho)*1e3,M,np.ones_like(x)*0.02)), fmt="%1.9g")
    fid.close()

# -----------------------------------------------------------------
#  ULB column text format
# -----------------------------------------------------------------

## ULB gas particles:
# - incoming:  x y z M h rho vx vy vz ...
# -    units:  x,y,z,h (100AU); M (Msun)
# - outgoing:  x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1)
def convert_gas_ULB(infile, outfile):
    PARSEC = 3.08568e16   # 1 parsec (in m)
    AU = 1.496e11         # 1 AU (in m)
    CONV = (100. * AU) / PARSEC
    x,y,z,M,h = np.loadtxt(infile, usecols=(0,1,2,3,4), unpack=True)
    fid = open(outfile, 'w')
    fid.write('# SPH Gas Particles\n')
    fid.write('# Converted from ULB output format into SKIRT6 format\n')
    fid.write('# Columns contain: x(pc) y(pc) z(pc) h(pc) M(Msun) Z(0-1)\n')
    np.savetxt(fid, np.transpose((x*CONV,y*CONV,z*CONV,5*h*CONV,M,np.zeros_like(M)+0.02)), fmt="%1.9g")  # inflated h!
    fid.close()

# -----------------------------------------------------------------
