#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.emissionlines Contains the EmissionLine and EmissionLines classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..units.parsing import parse_unit as u
from ..tools import strings

# -----------------------------------------------------------------

# define some relevant emission lines as a list of (label, center) tuples
#   Galaxies in the universe, Sparke & Gallagher, table 1.7
#   Inami et al, ApJ 2013
#   http://classic.sdss.org/dr4/algorithms/linestable.html
linedefs = [(0.1215, 0.1198, 0     , r"A$\mathrm{1}$", ("A", "1")),
            (0.2800, 0.2773, 0.2825, r"A$\mathrm{2}$", ("A", "2")),
            (0.3660, 0.3650, 0     , r"", None),
            (0.3730, 0.3715, 0.3745, r"O$\mathrm{II}$", ("O", "II")),
            (0.3875, 0.3850, 0.3910, r"He$\mathrm{I}$", ("He", "I")),
            (0.3970, 0.3950, 0.3987, r"Ca$\mathrm{II}$", ("Ca", "II")),
            (0.4104, 0.4090, 0.4125, r"H$\Delta$", ("H", "delta")),
            (0.4341, 0.4312, 0.4355, r"H$\gamma$", ("H", "gamma")),
            (0.4862, 0.4835, 0.4880, r"H$\beta$", ("H", "beta")),
            (0.4960, 0.4940, 0     , r"O$\mathrm{III}$", ("O", "III")),
            (0.5008, 0.4980, 0.5036, r"O$\mathrm{III}$", ("O", "III")),
            (0.5875, 0.5838, 0.5910, r"Na$\mathrm{I}$", ("Na", "I")),
            (0.6310, 0.6250, 0.6330, r"X$\mathrm{1}$", ("X", "1")),
            (0.6525, 0     , 0     , r"", None),
            (0.6565, 0.6545, 0.6607, r"H$\alpha$", ("H", "alpha")),
            (0.6719, 0.6699, 0.6780, r"S$\mathrm{II}$", ("S", "II")),
            (0.7137, 0.7110, 0.7160, r"X$\mathrm{2}$", ("X", "2")),
            (0.7755, 0.7720, 0.7785, r"X$\mathrm{3}$", ("X", "3")),
            (0.9070, 0.9010, 0.9100, r"X$\mathrm{4}$", ("X", "4")),
            (0.9550, 0.9500, 0.9600, r"X$\mathrm{5}$", ("X", "5")),

            (1.085, 1.077, 1.090, r"Y$\mathrm{1}$", ("Y", "1")),
            (1.282, 1.276, 1.288, r"Y$\mathrm{2}$", ("Y", "2")),
            (1.874, 1.862, 1.883, r"Y$\mathrm{3}$", ("Y", "3")),
            (2.163, 2.123, 2.172, r"Y$\mathrm{4}$", ("Y", "4")),
            (2.630, 2.610, 2.652, r"Y$\mathrm{5}$", ("Y", "5")),

            ( 4.05,  4.01,  4.09, r"Z$\mathrm{1}$", ("Z", "1")),
            ( 9.00,  8.84,  9.11, r"Z$\mathrm{2}$", ("Z", "2")),
            (10.51, 10.28, 10.60, r"S$\mathrm{IV}$", ("S", "IV")),
            (11.25, 0    , 0    , r"Z$\mathrm{3}$", ("Z", "3")),
            (12.80, 0    , 0    , r"Ne$\mathrm{II}$", ("Ne", "II")),
            (15.56, 15.20, 15.96, r"Ne$\mathrm{III}$", ("Ne", "III")),
            (18.20, 0    , 0    , r"Z$\mathrm{4}$", ("Z", "4")),
            (18.70, 18.50, 18.90, r"S$\mathrm{III}$", ("S", "III")),
            (33.67, 33.15, 34.20, r"S$\mathrm{III}$", ("S", "III")),
            (34.80, 34.20, 35.30, r"Si$\mathrm{II}$", ("Si", "II")),

            (51.85, 51.00, 53.40, r"O$\mathrm{III}$", ("O", "III")),
            (88.40, 86.30, 91.00, r"O$\mathrm{III}$", ("O", "III")),
            (157.5, 153.0, 161.0, r"C$\mathrm{II}$", ("C", "II")),
            (160.0, 0    , 0    , r"", None),
            (370.0, 359.5, 382.5, r"C$\mathrm{I}$", ("C", "I")),
            (604.3, 591.7, 617.8, r"C$\mathrm{I}$", ("C", "I")),
           ]

# -----------------------------------------------------------------

def get_identifiers():

    """
    This function ...
    :return:
    """

    ids = []

    # Create the lines
    for center, left, right, label, definition in linedefs:

        if definition is None: continue

        # Get ID
        group = definition[0]
        spec = definition[1]
        id = (group, spec)
        ids.append(id)

    # Return the IDs
    return ids

# -----------------------------------------------------------------

def get_id_strings():

    """
    This function ...
    :return:
    """

    return [group + spec for group, spec in get_identifiers()]

# -----------------------------------------------------------------

class EmissionLine(object):

    """
    This class ...
    """

    def __init__(self, group, specification, center, left, right, label):

        """
        This function ...
        """

        # Set attributes
        self.group = group
        self.specification = specification
        self.center = center
        self.left = left
        self.right = right
        self.label = label

    # -----------------------------------------------------------------

    @property
    def identifier(self):

        """
        This function ...
        :return:
        """

        return (self.group, self.specification)

    # -----------------------------------------------------------------

    @classmethod
    def from_string(cls, string):

        """
        This function ...
        :param string:
        :return:
        """

        # Loop over the linedefs
        for center, left, right, label, definition in linedefs:

            # Skip lines witout definition
            if definition is None: continue

            group = definition[0]
            spec = definition[1]

            # Generate aliases for this line
            aliases = list(strings.generate_from_two_parts(group, spec, connectors=(" ", "-", ".", "_", "")))

            # If the string matches one of the aliases
            if string in aliases: return cls(group, spec, center * u("micron"), left * u("micron"), right * u("micron"), label)

        raise ValueError("No line found that corresponds with the string '" + string + "'")

# -----------------------------------------------------------------

class EmissionLines(list):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Call the constructor of the base class
        super(EmissionLines, self).__init__()

        # Create the lines
        for center, left, right, label, definition in linedefs:

            if definition is not None:
                group = definition[0]
                spec = definition[1]
            else: group = spec = None

            # Create the line
            line = EmissionLine(group, spec, center * u("micron"), left * u("micron"), right * u("micron"), label)

            # Add the line
            self.append(line)

# -----------------------------------------------------------------
