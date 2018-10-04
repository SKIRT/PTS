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
linedefs = [(0.1215, 0.1198, 0     , r"A$\mathrm{1}$", ("A", "1"), "A1"),
            (0.2800, 0.2773, 0.2825, r"A$\mathrm{2}$", ("A", "2"), "A2"),
            (0.3660, 0.3650, 0     , r"", None, "u1"),
            (0.3730, 0.3715, 0.3745, r"O$\mathrm{II}$", ("O", "II"), "OII"),
            (0.3875, 0.3850, 0.3910, r"He$\mathrm{I}$", ("He", "I"), "HeI"),
            (0.3970, 0.3950, 0.3987, r"Ca$\mathrm{II}$", ("Ca", "II"), "CaII"),
            (0.4104, 0.4090, 0.4125, r"H$\Delta$", ("H", "delta"), "Hdelta"),
            (0.4341, 0.4312, 0.4355, r"H$\gamma$", ("H", "gamma"), "Hgamma"),
            (0.4862, 0.4835, 0.4880, r"H$\beta$", ("H", "beta"), "Hbeta"),
            (0.4960, 0.4940, 0     , r"O$\mathrm{III}$", ("O", "III"), "OIIIa"),
            (0.5008, 0.4980, 0.5036, r"O$\mathrm{III}$", ("O", "III"), "OIIIb"),
            (0.5875, 0.5838, 0.5910, r"Na$\mathrm{I}$", ("Na", "I"), "NaI"),
            (0.6310, 0.6250, 0.6330, r"X$\mathrm{1}$", ("X", "1"), "X1"),
            (0.6525, 0     , 0     , r"", None, "u2"),
            (0.6565, 0.6545, 0.6607, r"H$\alpha$", ("H", "alpha"), "Halpha"),
            (0.6719, 0.6699, 0.6780, r"S$\mathrm{II}$", ("S", "II"), "SII"),
            (0.7137, 0.7110, 0.7160, r"X$\mathrm{2}$", ("X", "2"), "X2"),
            (0.7755, 0.7720, 0.7785, r"X$\mathrm{3}$", ("X", "3"), "X3"),
            (0.9070, 0.9010, 0.9100, r"X$\mathrm{4}$", ("X", "4"), "X4"),
            (0.9550, 0.9500, 0.9600, r"X$\mathrm{5}$", ("X", "5"), "X5"),

            (1.085, 1.077, 1.090, r"Y$\mathrm{1}$", ("Y", "1"), "Y1"),
            (1.282, 1.276, 1.288, r"Y$\mathrm{2}$", ("Y", "2"), "Y2"),
            (1.874, 1.862, 1.883, r"Y$\mathrm{3}$", ("Y", "3"), "Y3"),
            (2.163, 2.123, 2.172, r"Y$\mathrm{4}$", ("Y", "4"), "Y4"),
            (2.630, 2.610, 2.652, r"Y$\mathrm{5}$", ("Y", "5"), "Y5"),

            ( 4.05,  4.01,  4.09, r"Z$\mathrm{1}$", ("Z", "1"), "Z1"),
            ( 9.00,  8.84,  9.11, r"Z$\mathrm{2}$", ("Z", "2"), "Z2"),
            (10.51, 10.28, 10.60, r"S$\mathrm{IV}$", ("S", "IV"), "SIV"),
            (11.25, 0    , 0    , r"Z$\mathrm{3}$", ("Z", "3"), "Z3"),
            (12.80, 0    , 0    , r"Ne$\mathrm{II}$", ("Ne", "II"), "NeII"),
            (15.56, 15.20, 15.96, r"Ne$\mathrm{III}$", ("Ne", "III"), "NeIII"),
            (18.20, 0    , 0    , r"Z$\mathrm{4}$", ("Z", "4"), "Z4"),
            (18.70, 18.50, 18.90, r"S$\mathrm{III}$", ("S", "III"), "SIIIa"),
            (33.67, 33.15, 34.20, r"S$\mathrm{III}$", ("S", "III"), "SIIIb"),
            (34.80, 34.20, 35.30, r"Si$\mathrm{II}$", ("Si", "II"), "SiII"),

            (51.85, 51.00, 53.40, r"O$\mathrm{III}$", ("O", "III"), "OIIIc"),
            (88.40, 86.30, 91.00, r"O$\mathrm{III}$", ("O", "III"), "OIIId"),
            (157.5, 153.0, 161.0, r"C$\mathrm{II}$", ("C", "II"), "CII"),
            (160.0, 0    , 0    , r"", None, "u3"),
            (370.0, 359.5, 382.5, r"C$\mathrm{I}$", ("C", "I"), "CIa"),
            (604.3, 591.7, 617.8, r"C$\mathrm{I}$", ("C", "I"), "CIb"),
           ]

# -----------------------------------------------------------------

def get_identifiers():

    """
    This function ...
    :return:
    """

    ids = []

    # Get the identifiers
    for center, left, right, label, definition, identifier in linedefs: ids.append(identifier)

    # Return the IDs
    return ids

# -----------------------------------------------------------------

important_lines = ["Halpha", "Hbeta", "Hgamma", "Hdelta", "A1", "A2"]
strong_lines = ["A1", "OII", "Hbeta", "OIIIb", "Halpha", "SII", "X4", "X5", "Y1", "Y3", "Z1"] # in order of appearance

# -----------------------------------------------------------------

class EmissionLine(object):

    """
    This class ...
    """

    def __init__(self, id, group, specification, center, left, right, label):

        """
        This function ...
        :param id:
        :param group:
        :param specification:
        :param center:
        :param left:
        :param right:
        :param label:
        """

        # Set attributes
        self.id = id
        self.group = group
        self.specification = specification
        self.center = center
        self.left = left
        self.right = right
        self.label = label

    # -----------------------------------------------------------------

    @property
    def identifier(self):
        return self.id

    # -----------------------------------------------------------------

    @classmethod
    def from_string(cls, string):

        """
        This function ...
        :param string:
        :return:
        """

        # Loop over the linedefs
        for center, left, right, label, definition, identifier in linedefs:

            # Get group and spec
            group = definition[0] if definition is not None else None
            spec = definition[1] if definition is not None else None

            # Check for exact match with identifier
            if identifier == string: return cls(identifier, group, spec, center * u("micron"), left * u("micron"), right * u("micron"), label)

            # Skip lines without definition
            if definition is None: continue

            # Generate aliases for this line
            aliases = list(strings.generate_from_two_parts(group, spec, connectors=(" ", "-", ".", "_", "")))
            # WARNING: if there are multiple lines with the same group and spec, only the first line is returned!!

            # If the string matches one of the aliases
            if string in aliases: return cls(identifier, group, spec, center * u("micron"), left * u("micron"), right * u("micron"), label)

        # No match found
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
        for center, left, right, label, definition, identifier in linedefs:

            # Get group and spec
            group = definition[0] if definition is not None else None
            spec = definition[1] if definition is not None else None

            # Create the line
            line = EmissionLine(identifier, group, spec, center * u("micron"), left * u("micron"), right * u("micron"), label)

            # Add the line
            self.append(line)

# -----------------------------------------------------------------
