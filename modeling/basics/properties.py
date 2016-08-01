#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.basics.properties Contains the GalaxyProperties class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import parsing
from ...magic.basics.skygeometry import SkyCoordinate

# -----------------------------------------------------------------

class GalaxyProperties(object):
    
    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        self.name = None
        self.ngc_id = None
        self.center = None
        self.major = None
        self.ellipticity = None
        self.position_angle = None
        self.distance = None
        self.distance_error = None
        self.inclination = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :return:
        """

        # Create the properties instance
        properties = cls()

        ra = None
        dec = None

        # Read the parameter file
        with open(path, 'r') as parameter_file:

            # Loop over all lines in the file
            for line in parameter_file:

                splitted = line.split(": ")

                if splitted[0] == "Name":
                    properties.name = splitted[1]
                elif splitted[0] == "Center RA":
                    ra = parsing.quantity(splitted[1])
                elif splitted[0] == "Center DEC":
                    dec = parsing.quantity(splitted[1])
                elif splitted[0] == "Major axis length":
                    properties.major = parsing.quantity(splitted[1])
                elif splitted[0] == "Ellipticity":
                    properties.ellipticity = float(splitted[1])
                elif splitted[0] == "Position angle":
                    properties.position_angle = parsing.angle(splitted[1])
                elif splitted[0] == "Distance":
                    properties.distance = parsing.quantity(splitted[1])
                elif splitted[0] == "Distance error":
                    properties.distance_error = parsing.quantity(splitted[1])
                elif splitted[0] == "Inclination":
                    properties.inclination = parsing.angle(splitted[1])
                elif splitted[0] == "IRAC 3.6um flux density":
                    properties.i1_fluxdensity = parsing.quantity(splitted[1])
                elif splitted[0] == "IRAC 3.6um flux density error":
                    properties.i1_error = parsing.quantity(splitted[1])
                elif splitted[0] == "IRAC 4.5um flux density":
                    properties.i2_fluxdensity = parsing.quantity(splitted[1])
                elif splitted[0] == "IRAC 4.5um flux density error":
                    properties.i2_error = parsing.quantity(splitted[1])

        # Add the center coordinate
        properties.center = SkyCoordinate(ra=ra, dec=dec)

        # Return the properties instance
        return properties

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
            print("Name:", self.name, file=parameter_file)
            print("NGC ID:", self.ngc_id, file=parameter_file)
            print("Center RA:", str(self.center.ra.to("deg").value) + " deg", file=parameter_file)
            print("Center DEC:", str(self.center.dec.to("deg").value) + " deg", file=parameter_file)
            print("Major axis length:", str(self.major), file=parameter_file)
            print("Ellipticity:", self.ellipticity, file=parameter_file)
            print("Position angle:", str(self.position_angle.to("deg").value) + " deg", file=parameter_file)
            print("Distance:", str(self.distance), file=parameter_file)
            print("Distance error:", str(self.distance_error), file=parameter_file)
            print("Inclination:", str(self.inclination.to("deg").value) + " deg", file=parameter_file)

            #print("IRAC 3.6um flux density:", properties.i1_fluxdensity, file=parameter_file)
            #print("IRAC 3.6um flux density error:", properties.i1_error, file=parameter_file)
            #print("IRAC 4.5um flux density:", properties.i2_fluxdensity, file=parameter_file)
            #print("IRAC 4.5um flux density error:", properties.i2_error, file=parameter_file)

# -----------------------------------------------------------------
