#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.mappings Contains the Mappings class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ...core.data.sed import SED
from ...core.tools import introspection, tables
from ...core.tools import filesystem as fs
from ...core.units.parsing import parse_unit as u
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

# Determine the path to the Bruzual-Charlot SED directory
bc_path = fs.join(introspection.skirt_repo_dir, "dat", "SED", "BruzualCharlot")
chabrier_path = fs.join(bc_path, "chabrier") # chabrier IMF

# -----------------------------------------------------------------

class BruzualCharlot(object):
    
    """
    This class ...
    """

    def __init__(self, metallicity, age):

        """
        The constructor ...
        :param metallicity:
        :param age:
        """

        # Set properties
        self.metallicity = metallicity
        self.age = age

    # -----------------------------------------------------------------

    @lazyproperty
    def sed(self):

        """
        This function ...
        :return:
        """

        # Create the intrinsic SED and return it
        return create_mappings_sed(self.metallicity, self.pressure, self.compactness, self.covering_factor, self.sfr)

    # -----------------------------------------------------------------

    @classmethod
    def dust_mass_for_luminosity(cls, metallicity, compactness, pressure, covering_factor, luminosity, wavelength):

        """
        This function ...
        :param metallicity:
        :param compactness:
        :param pressure:
        :param covering_factor:
        :param luminosity:
        :param wavelength:
        :return:
        """



    # -----------------------------------------------------------------

    def luminosity_at(self, wavelength, unit="W/micron"):

        """
        This function ...
        :param wavelength:
        :param unit:
        :return:
        """

        return self.sed.photometry_at(wavelength).to(unit)

    # -----------------------------------------------------------------

    def luminosity_for_filter(self, fltr, unit="W/micron"):

        """
        This function ...
        :param fltr:
        :param unit:
        :return:
        """

        #luminosity = filter.integrate(self.sed["Wavelength"], self.sed["Luminosity"])
        luminosity = fltr.convolve(self.sed.wavelengths(unit="micron", asarray=True), self.sed.photometry(unit="W/micron", asarray=True)) # also in W/micron
        luminosity = luminosity * u("W/micron")

        # Return the luminosity in the desired unit
        return luminosity.to(unit)

# -----------------------------------------------------------------

def create_bruzual_charlot_sed(metallicity, age):

    """
    This function ...
    :param metallicity:
    :param age:
    :return:
    """

    // find the appropriate SED from interpolating in the BC library
    int mL, mR;
    double hZ = 0.0;
    if (Z<=_Zv[0])
        mL = mR = 0;
    else if (Z>=_Zv[NZ-1])
        mL = mR = NZ-1;
    else
    {
        mL = NR::locate_clip(_Zv,Z);
        mR = mL+1;
        double ZL = _Zv[mL];
        double ZR = _Zv[mR];
        hZ = (Z-ZL)/(ZR-ZL);
    }
    int pL, pR;
    double ht = 0.0;
    if (t<=_tv[0])
        pL = pR = 0;
    else if (t>=_tv[Nt-1])
        pL = pR = Nt-1;
    else
    {
        pL = NR::locate_clip(_tv,t);
        pR = pL+1;
        double tL = _tv[pL];
        double tR = _tv[pR];
        ht = (t-tL)/(tR-tL);
    }
    const Array& jLLv = _jvv(pL,mL);
    const Array& jLRv = _jvv(pL,mR);
    const Array& jRLv = _jvv(pR,mL);
    const Array& jRRv = _jvv(pR,mR);
    Array jv(Nlambda);
    for (int k=0; k<Nlambda; k++)
        jv[k] = (1.0-ht)*(1.0-hZ)*jLLv[k]
                + (1.0-ht)*hZ*jLRv[k]
                + ht*(1.0-hZ)*jRLv[k]
                + ht*hZ*jRRv[k];

    // resample to the possibly redshifted simulation wavelength grid,
    // convert emissivities to luminosities (i.e. multiply by the wavelength bins),
    // multiply by the mass of the population (in solar masses),
    // and return the result
    return NR::resample<NR::interpolate_loglog>(_lambdagrid->lambdav()*(1-z), _lambdav, jv)
                                    * _lambdagrid->dlambdav() * M;

# -----------------------------------------------------------------

def load_bruzual_charlot_data():

    """
    This function ...
    :return:
    """

    # // number of items in the library read by the constructor
    # const int Nlambda = 1221;
    # const int Nt = 221;
    # const int NZ = 6;

    // local constants for units
    const double Lsun = Units::Lsun();
    const double Angstrom = 1e-10;

    // Prepare the vectors for the Bruzual & Charlot library SEDs
    _lambdav.resize(Nlambda);
    _tv.resize(Nt);
    _Zv.resize(NZ);
    _jvv.resize(Nt,NZ,Nlambda);

    // Fill the metallicity vector
    vector<QString> Zcodev(NZ);
    _Zv[0] = 0.0001;     Zcodev[0] = "m22";
    _Zv[1] = 0.0004;     Zcodev[1] = "m32";
    _Zv[2] = 0.004;      Zcodev[2] = "m42";
    _Zv[3] = 0.008;      Zcodev[3] = "m52";
    _Zv[4] = 0.02;	     Zcodev[4] = "m62";
    _Zv[5] = 0.05;	     Zcodev[5] = "m72";

    // Read the wavelength, age and emissivity vectors from the Bruzual & Charlot library
    for (int m=0; m<NZ; m++)
    {
        QString bcfilename = FilePaths::resource("SED/BruzualCharlot/chabrier/bc2003_lr_"
                                                 + Zcodev[m] + "_chab_ssp.ised_ASCII");
        ifstream bcfile(bcfilename.toLocal8Bit().constData());
        if (! bcfile.is_open()) throw FATALERROR("Could not open the data file " + bcfilename);

        find<Log>()->info("Reading SED data from file " + bcfilename + "...");
        int iNt, iNlambda;
        bcfile >> iNt;
        if (iNt != Nt)
            throw FATALERROR("iNt is not equal to Nt");
        for (int p=0; p<Nt; p++)
        {
            double t;
            bcfile >> t;
            _tv[p] = t;     // age in file in yr, we want in yr
        }
        string dummy;
        for (int l=0; l<6; l++)
            getline(bcfile,dummy); // skip six lines...
        bcfile >> iNlambda;
        if (iNlambda != Nlambda)
            throw FATALERROR("iNlambda is not equal to Nlambda");
        for (int k=0; k<Nlambda; k++)
        {
            double lambda;
            bcfile >> lambda;
            _lambdav[k] = lambda * Angstrom;   // lambda in file in A, we want in m
        }
        for (int p=0; p<Nt; p++)
        {
            Array& jv = _jvv(p,m);
            bcfile >> iNlambda;
            if (iNlambda != Nlambda)
                throw FATALERROR("iNlambda is not equal to Nlambda");
            for (int k=0; k<Nlambda; k++)
            {
                double j;
                bcfile >> j;
                jv[k] = j * Lsun/Angstrom;   // emissivity in file in Lsun/A, we want in W/m.
            }
            int idummy;
            bcfile >> idummy;
            for (int k=0; k<idummy; k++)
            {
                double dummy;
                bcfile >> dummy;
            }
        }
        bcfile.close();
        find<Log>()->info("File " + bcfilename + " closed.");
    }

    // cache the simulation's wavelength grid
    _lambdagrid = find<WavelengthGrid>();

# -----------------------------------------------------------------

def find_nearest(array, value):

    """
    This function ...
    :param array:
    :param value:
    :return:
    """

    idx = (np.abs(array-value)).argmin()
    return array[idx]

# -----------------------------------------------------------------

def locate_clip(array, value):

    """
    This function ...
    :param array: actually a list
    :param value:
    :return:
    """

    #n = array.size
    n = len(array)
    if value < array[0]: return 0
    return locate_basic_impl(array, value, n-1)

# -----------------------------------------------------------------

def locate_basic_impl(xv, x, n):

    """
    This function ...
    :param xv:
    :param x:
    :param n:
    :return:
    """

    jl = -1
    ju = n

    while ju - jl > 1:

        jm = (ju + jl) >> 1
        if x < xv[jm]: ju = jm
        else: jl = jm

    return jl

# -----------------------------------------------------------------
