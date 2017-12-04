#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.data.sun Contains the Sun class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from scipy import interpolate

# Import astronomical modules
#from astropy.units import spectral

# Import the relevant PTS classes and modules
from ...core.data.sed import SED
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ..units.parsing import parse_unit as u

# -----------------------------------------------------------------

# Determine the path to the Sun SED file
sun_sed_path = fs.join(introspection.skirt_repo_dir, "dat", "SED", "Sun", "SunSED.dat")

# -----------------------------------------------------------------

#speed_of_light = constants.c
#solar_luminosity = 3.846e26 * u("W") # 3.828×10^26 W

Lsun = 3.839e26 * u("W")                 # solar luminosity without solar neutrino radiation

# -----------------------------------------------------------------

# Effective wavelengths (in m)
eff_wavelengths = {"FUV": 152e-9 * u("m"),
                   "NUV": 231e-9 * u("m"),
                   "U": 365e-9 * u("m"),
                   "B": 445e-9 * u("m"),
                   "V": 551e-9 * u("m"),
                   "R": 658e-9 * u("m"),
                   "I": 806e-9 * u("m"),
                   "J": 1.22e-6 * u("m"),
                   "H": 1.63e-6 * u("m"),
                   "K": 2.19e-6 * u("m"),
                   "SDSS u": 354e-9 * u("m"),
                   "SDSS g": 477e-9 * u("m"),
                   "SDSS r": 623e-9 * u("m"),
                   "SDSS i": 763e-9 * u("m"),
                   "SDSS z": 913e-9 * u("m"),
                   "IRAC1": 3.56e-6 * u("m"),
                   "IRAC2": 4.51e-6 * u("m"),
                   "WISE1": 3.35e-9 * u("m"),
                   "WISE2": 4.60e-6 * u("m")}

# -----------------------------------------------------------------

#	filter	Vega	AB	Apparent	AB mag of Vega
#1	Buser U (= Johnson)	5.59	6.32	-25.99	0.734
#2	Azusienis & Straizys B (= Johnson)	5.45	5.36	-26.13	-0.098
#3	Azusienis & Straizys V (= Johnson)	4.78	4.79	-26.80	0.011
#4	Bessell R	4.46	4.65	-27.12	0.193
#5	Bessell I	4.11	4.55	-27.47	0.443
#6	F300W	6.09	7.52	-25.50	1.435
#7	F450W	5.32	5.25	-26.26	-0.071
#8	F555W	4.85	4.85	-26.74	-0.001
#9	F606W	4.66	4.75	-26.92	0.096
#10	F702W	4.32	4.58	-27.26	0.263
#11	F814W	4.15	4.57	-27.43	0.417
#12	CFHT U	5.57	6.38	-26.01	0.809
#13	CFHT B	5.49	5.32	-26.09	-0.170
#14	CFHT V	4.81	4.81	-26.77	0.000
#15	CFHT R	4.44	4.64	-27.14	0.199
#16	CFHT I	4.06	4.54	-27.52	0.480
#17	KPNO U	5.59	6.32	-25.99	0.729
#18	KPNO B	5.49	5.43	-26.09	-0.059
#19	KPNO V	4.79	4.79	-26.79	0.008
#20	KPNO R	4.47	4.66	-27.11	0.185
#21	KPNO I	4.11	4.55	-27.47	0.439
#22	Koo & Kron U	5.58	6.29	-26.00	0.714
#23	Koo & Kron J	5.31	5.26	-26.27	-0.047
#24	Koo & Kron F	4.58	4.69	-27.01	0.117
#25	Koo & Kron N	4.11	4.53	-27.47	0.421
#26	SDSS u'	5.46	6.45	-26.12	0.995
#27	SDSS g'	5.22	5.14	-26.36	-0.087
#28	SDSS r'	4.50	4.65	-27.08	0.147
#29	SDSS i'	4.16	4.54	-27.42	0.376
#30	SDSS z'	4.01	4.52	-27.58	0.518
#31	ACS old z	3.99	4.52	-27.59	0.530
#32	FIS555	4.84	4.84	-26.74	0.002
#33	FIS606	4.63	4.72	-26.96	0.096
#34	FIS702	4.32	4.59	-27.26	0.261
#35	FIS814	4.12	4.53	-27.47	0.417
#36	LRIS B	5.46	5.42	-26.12	-0.043
#37	LRIS V	4.82	4.83	-26.76	0.005
#38	LRIS R	4.46	4.63	-27.12	0.174
#39	LRIS Rs	4.33	4.59	-27.25	0.256
#40	LRIS I	4.04	4.53	-27.54	0.484
#41	LRIS Z	4.00	4.52	-27.58	0.520
#42	SPH Un	5.43	6.49	-26.15	1.063
#43	SPH G	5.21	5.11	-26.37	-0.098
#44	SPH Rs	4.39	4.61	-27.19	0.173
#45	12k B	5.45	5.34	-26.13	-0.111
#46	12k R	4.39	4.60	-27.20	0.218
#47	12k I	4.10	4.53	-27.49	0.435
#48	12k V	4.85	4.85	-26.73	-0.001
#49	uh8K i	4.06	4.53	-27.52	0.461
#50	ACS B435	5.49	5.40	-26.09	-0.083
#51	ACS V606	4.67	4.75	-26.91	0.080
#52	ACS SDSS i	4.14	4.54	-27.44	0.393
#53	ACS I814	4.11	4.53	-27.47	0.424
#54	ACS SDSS z	4.00	4.52	-27.58	0.522
#55	Bessell U	5.55	6.36	-26.03	0.807
#56	Bessell B	5.45	5.36	-26.13	-0.096
#57	Bessell V	4.80	4.82	-26.78	0.012
#58	Bessell J	3.67	4.57	-27.91	0.897
#59	Bessell H	3.33	4.71	-28.25	1.379
#60	Bessell K	3.29	5.19	-28.29	1.882
#61	KPNO J	3.66	4.57	-27.92	0.905
#62	KPNO H	3.33	4.71	-28.26	1.387
#63	KPNO K	3.29	5.18	-28.29	1.891
#64	2500	7.96	9.80	-23.62	1.841
#65	2800	6.67	8.23	-24.91	1.557
#66	APM Bj	5.29	5.21	-26.29	-0.073
#67	FOCA UV	10.50	12.39	-21.08	1.893
#68	DEIMOS R	4.44	4.62	-27.14	0.186
#69	Galex FUV	13.97	16.42	-17.61	2.457
#70	Galex NUV	8.45	10.31	-23.13	1.863
#71	SDSS u z=0.1	5.83	6.77	-25.75	0.938
#72	SDSS g z=0.1	5.46	5.36	-26.12	-0.091
#73	SDSS r z=0.1	4.53	4.67	-27.05	0.144
#74	SDSS i z=0.1	4.12	4.48	-27.46	0.362
#75	SDSS z z=0.1	3.90	4.42	-27.68	0.520
#76	NIRI J	3.64	4.57	-27.94	0.927
#77	NIRI H	3.33	4.71	-28.25	1.380
#78	NIRI K	3.29	5.18	-28.29	1.890

# -----------------------------------------------------------------

# SOLAR MAGNITUDES IN DIFFERENT BANDS:
# WHICH KINDS OF MAGNITUDES??
# Filter	B&M	here	difference
#U	5.61	5.56	0.05
#B	5.48	5.45	0.03
#V	4.83	4.80	0.03
#R	4.42	4.46	-0.04
#I	4.08	4.10	-0.02
#J	3.64	3.66	-0.02
#H	3.32	3.32	0.00
#K	3.28	3.28	0.00

# -----------------------------------------------------------------

LsunK = 10**((34.1-5.19)/2.5)  # solar luminosity in K band expressed in W/Hz  (AB magnitude is 5.19)
LsunH = 10**((34.1-4.71)/2.5)  # solar luminosity in H band expressed in W/Hz  (AB magnitude is 4.71)

# -----------------------------------------------------------------

class Sun(object):
    
    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Load the intrinsic SED of the sun
        self.sed = SED.from_text_file(sun_sed_path, wavelength_unit="micron", photometry_unit="W/micron", skiprows=4) # wavelength in micron, luminosity in W/micron

        # The total luminosity
        self.luminosity = Lsun

    # -----------------------------------------------------------------

    def total_luminosity(self, unit="W"):

        """
        This function ...
        :param unit:
        :return:
        """

        return Lsun.to(unit)

    # -----------------------------------------------------------------

    def total_luminosity_as_unit(self):

        """
        This function ...
        :return:
        """

        # Create and return the new unit
        return u("Lsun", self.total_luminosity())

    # -----------------------------------------------------------------

    def luminosity_for_wavelength(self, wavelength, unit="W/micron", density=False, density_strict=False):

        """
        This function ...
        :param wavelength: 
        :param unit:
        :param density:
        :param density_strict:
        :return: 
        """

        wavelengths = self.sed.wavelengths(unit="micron", asarray=True)
        luminosities = self.sed.photometry(unit="W/micron", asarray=True)

        interpolated = interpolate.interp1d(wavelengths, luminosities, kind='linear')
        value = interpolated(wavelength.to("micron").value) * u("W/micron")

        # Return the spectral luminosity
        return value.to(unit, wavelength=wavelength, density=density, density_strict=density_strict) #, equivalencies=spectral()) # using equivalencies invokes an error because astropy sees conversion between W/micron and W/Hz, which it sees as an ENERGY

    # -----------------------------------------------------------------

    def luminosity_for_filter(self, fltr, unit="W/micron", density=False, density_strict=False):

        """
        This function ...
        :param fltr:
        :param unit:
        :param density:
        :param density_strict:
        :return:
        """

        # Convole the Sun SED over the filter transmission curve
        luminosity = fltr.convolve(self.sed.wavelengths(unit="micron", asarray=True), self.sed.photometry(unit="W/micron", asarray=True)) # also in W/micron
        luminosity = luminosity * u("W/micron")

        # Return the spectral luminosity
        return luminosity.to(unit, wavelength=fltr.wavelength, density=density, density_strict=density_strict) #, equivalencies=spectral()) # using equivalencies invokes an error because astropy sees conversion between W/micron and W/Hz, which it sees as an ENERGY

    # -----------------------------------------------------------------

    def luminosity_for_filter_as_unit(self, fltr, density=False, density_strict=False):

        """
        This function ...
        :param fltr:
        :param density:
        :param density_strict:
        :return:
        """

        # Determine a name for this unit
        unit_name = "Lsun_" + fltr.band

        # Create and return the new unit
        from astropy.units import Unit
        return Unit(unit_name, represents=self.luminosity_for_filter(fltr, density=density, density_strict=density_strict))

# -----------------------------------------------------------------
