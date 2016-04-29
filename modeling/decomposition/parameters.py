#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.decomposition.parameters Contains the DecompositionParameters class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...core.basics.map import Map
from ...magic.basics.skygeometry import SkyCoordinate

# -----------------------------------------------------------------

class DecompositionParameters(object):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Basic properties
        self.galaxy_name = None
        self.center = None
        self.major = None
        self.ellipticity = None
        self.position_angle = None
        self.distance = None
        self.distance_error = None
        self.inclination = None
        self.i1_fluxdensity = None
        self.i1_error = None
        self.i2_fluxdensity = None
        self.i2_error = None

        # Bulge properties
        self.bulge = Map()
        self.bulge.f = None
        self.bulge.fluxdensity = None
        self.bulge.q = None
        self.bulge.PA = None
        self.bulge.Re = None
        self.bulge.n = None

        # Disk properties
        self.disk = Map()
        self.disk.f = None
        self.disk.fluxdensity = None
        self.disk.q = None
        self.disk.PA = None
        self.disk.mu0 = None
        self.disk.hr = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new class instance
        parameters = cls()

        ra = None
        dec = None

        # Read the parameter file
        with open(path, 'r') as parameter_file:

            # Loop over all lines in the file
            for line in parameter_file:

                splitted = line.split(": ")

                # Bulge parameters
                if splitted[0] == "Bulge":

                    if splitted[1] == "Relative contribution": parameters.bulge.f = float(splitted[2])
                    elif splitted[1] == "IRAC 3.6um flux density": parameters.bulge.fluxdensity = get_quantity(splitted[2])
                    elif splitted[1] == "Axial ratio": parameters.bulge.q = float(splitted[2])
                    elif splitted[1] == "Position angle": parameters.bulge.PA = get_angle(splitted[2])
                    elif splitted[1] == "Effective radius": parameters.bulge.Re = get_quantity(splitted[2])
                    elif splitted[1] == "Sersic index": parameters.bulge.n = float(splitted[2])

                # Disk parameters
                elif splitted[0] == "Disk":

                    if splitted[1] == "Relative contribution": parameters.disk.f = float(splitted[2])
                    elif splitted[1] == "IRAC 3.6um flux density": parameters.disk.fluxdensity = get_quantity(splitted[2])
                    elif splitted[1] == "Axial ratio": parameters.disk.q = float(splitted[2])
                    elif splitted[1] == "Position angle": parameters.disk.PA = get_angle(splitted[2])
                    elif splitted[1] == "Central surface brightness": parameters.disk.mu0 = get_quantity(splitted[2])
                    elif splitted[1] == "Exponential scale length": parameters.disk.hr = get_quantity(splitted[2])

                # Other parameters
                elif len(splitted) == 2:

                    if splitted[0] == "Name": parameters.galaxy_name = splitted[1]
                    elif splitted[0] == "Center RA": ra = get_quantity(splitted[1])
                    elif splitted[0] == "Center DEC": dec = get_quantity(splitted[1])
                    elif splitted[0] == "Major axis length": parameters.major = get_quantity(splitted[1])
                    elif splitted[0] == "Ellipticity": parameters.ellipticity = float(splitted[1])
                    elif splitted[0] == "Position angle": parameters.position_angle = get_angle(splitted[1])
                    elif splitted[0] == "Distance": parameters.distance = get_quantity(splitted[1])
                    elif splitted[0] == "Distance error": parameters.distance_error = get_quantity(splitted[1])
                    elif splitted[0] == "Inclination": parameters.inclination = get_angle(splitted[1])
                    elif splitted[0] == "IRAC 3.6um flux density": parameters.i1_fluxdensity = get_quantity(splitted[1])
                    elif splitted[0] == "IRAC 3.6um flux density error": parameters.i1_error = get_quantity(splitted[1])
                    elif splitted[0] == "IRAC 4.5um flux density": parameters.i2_fluxdensity = get_quantity(splitted[1])
                    elif splitted[0] == "IRAC 4.5um flux density error": parameters.i2_error = get_quantity(splitted[1])

        # Add the center coordinate
        parameters.center = SkyCoordinate(ra=ra, dec=dec)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create the parameters file
        with open(path, 'w') as parameter_file:

            # Add general info
            print("Name:", self.galaxy_name, file=parameter_file)
            print("Center RA:", str(self.center.ra.to("deg").value) + " deg", file=parameter_file)
            print("Center DEC:", str(self.center.dec.to("deg").value) + " deg", file=parameter_file)
            print("Major axis length:", str(self.major), file=parameter_file)
            print("Ellipticity:", self.ellipticity, file=parameter_file)
            print("Position angle:", str(self.position_angle.to("deg").value) + " deg", file=parameter_file)
            print("Distance:", str(self.distance), file=parameter_file)
            print("Distance error:", str(self.distance_error), file=parameter_file)
            print("Inclination:", str(self.inclination.to("deg").value) + " deg", file=parameter_file)
            print("IRAC 3.6um flux density:", self.i1_fluxdensity, file=parameter_file)
            print("IRAC 3.6um flux density error:", self.i1_error, file=parameter_file)
            print("IRAC 4.5um flux density:", self.i2_fluxdensity, file=parameter_file)
            print("IRAC 4.5um flux density error:", self.i2_error, file=parameter_file)

            # print("Model type:", parameters.model_type, file=parameter_file)
            # print("Number of components:", parameters.number_of_components, file=parameter_file)
            # print("Quality:", parameters.quality, file=parameter_file)

            # Add components parameters
            for component in ["bulge", "disk"]:

                # print(component.title() + ": Interpretation:", parameters[component].interpretation, file=parameter_file)
                print(component.title() + ": Relative contribution:", self[component].f, file=parameter_file)
                print(component.title() + ": IRAC 3.6um flux density:", self[component].fluxdensity, file=parameter_file)
                print(component.title() + ": Axial ratio:", self[component].q, file=parameter_file)
                print(component.title() + ": Position angle:", str(self[component].PA.to("deg").value) + " deg", file=parameter_file)  # (degrees ccw from North)

                if component == "bulge":

                    print(component.title() + ": Effective radius:", str(self[component].Re), file=parameter_file)
                    print(component.title() + ": Sersic index:", self[component].n, file=parameter_file)

                elif component == "disk":

                    print(component.title() + ": Central surface brightness:", self[component].mu0, file=parameter_file)  # (mag/arcsec2)
                    print(component.title() + ": Exponential scale length:", str(self[component].hr), file=parameter_file)

# -----------------------------------------------------------------

def get_quantity(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    splitted = entry.split()
    value = float(splitted[0])
    try: unit = splitted[1]
    except IndexError: unit = default_unit

    # Create a quantity object and return it
    if unit is not None: value = value * Unit(unit)
    return value

# -----------------------------------------------------------------

def get_angle(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    splitted = entry.split()
    value = float(splitted[0])
    try: unit = splitted[1]
    except IndexError: unit = default_unit

    # Create an Angle object and return it
    if unit is not None: value = Angle(value, unit)
    return value

# -----------------------------------------------------------------
