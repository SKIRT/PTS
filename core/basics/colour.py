#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.colour Contains the Colour class.

# Compatibility between python 2 and 3
from __future__ import print_function

# -----------------------------------------------------------------

def hex_to_rgb(hex):

    """
    This function ...
    #FFFFFF" -> [255,255,255]
    """

    # Pass 16 to the integer function for change of base

    return [int(hex[i:i+2], 16) for i in range(1,6,2)]

# -----------------------------------------------------------------

def rgb_to_hex(rgb):

  """
  [255,255,255] -> "#FFFFFF"
  """

  # Components need to be integers for hex to make sense
  rgb = [int(x) for x in rgb]

  return "#"+"".join(["0{0:x}".format(v) if v < 16 else
            "{0:x}".format(v) for v in rgb])

# -----------------------------------------------------------------

def parse_colour(string):

    """
    This function ...
    :param string:
    :return:
    """

    from ..tools import types

    if types.is_string_type(string):

        if string.startswith("#"): return Colour.from_hex(string)
        elif list(string).count(",") == 2:
            red = float(string.split(",")[0])
            green = float(string.split(",")[1])
            blue = float(string.split(",")[2])
            return Colour.from_rgb(red, green, blue)
        else: return Colour.from_name(string)

    elif types.is_sequence(string): return Colour.from_rgb(string[0], string[1], string[2])
    elif isinstance(string, Colour): return string
    else: raise ValueError("Invalid input")

# -----------------------------------------------------------------

predefined = dict()
predefined["black"] = ("#000000", (0,0,0))
predefined["white"] = ("#FFFFFF", (255,255,255))
predefined["red"] = ("#FF0000",	(255,0,0))
predefined["lime"] = ("#00FF00", (0,255,0))
predefined["blue"] = ("#0000FF", (0,0,255))
predefined["yellow"] = ("#FFFF00", (255,255,0))
predefined["cyan"] = ("#00FFFF", (0,255,255))
predefined["magenta"] = ("#FF00FF", (255,0,255))
predefined["silver"] = ("#C0C0C0", (192,192,192))
predefined["gray"] = ("#808080", (128,128,128))
predefined["maroon"] = ("#800000", (128,0,0))
predefined["olive"] = ("#808000", (128,128,0))
predefined["green"] = ("#008000", (0,128,0))
predefined["purple"] = ("#800080", (128,0,128))
predefined["teal"] = ("#008080", (0,128,128))
predefined["navy"] = ("#000080", (0,0,128))

# -----------------------------------------------------------------

class Colour(object):

    """
    This function ...
    """

    def __init__(self, red, green, blue):

        """
        This function ...
        :param red:
        :param green:
        :param blue:
        """

        self.red = red
        self.green = green
        self.blue = blue

    # -----------------------------------------------------------------

    @classmethod
    def from_rgb(cls, red, green, blue):

        """
        This function ...
        :param red:
        :param green:
        :param blue:
        :return:
        """

        return cls(red, green, blue)

    # -----------------------------------------------------------------

    @classmethod
    def from_hex(cls, hex):

        """
        This function ...
        :param hex:
        :return:
        """

        red, green, blue = hex_to_rgb(hex)
        return cls.from_rgb(red, green, blue)

    # -----------------------------------------------------------------

    @classmethod
    def from_name(cls, name):

        """
        :param name:
        :return:
        """

        if name.lower() in predefined: return cls.from_hex(predefined[name.lower()][0])
        else: raise ValueError("Colour '" + name + "' not recognized")

    # -----------------------------------------------------------------

    @property
    def rgb(self):

        """
        This function ...
        :return:
        """

        return self.red, self.green, self.blue

    # -----------------------------------------------------------------

    @property
    def hex(self):

        """
        This function ...
        :return:
        """

        return rgb_to_hex([self.red, self.green, self.blue])

    # -----------------------------------------------------------------

    @property
    def hex1(self):

        """
        This function ...
        :return:
        """

        return self.hex[1:3]

    # -----------------------------------------------------------------

    @property
    def hex2(self):

        """
        This function ...
        :return:
        """

        return self.hex[3:5]

    # -----------------------------------------------------------------

    @property
    def hex3(self):

        """
        This function ...
        :return:
        """

        return self.hex[5:]

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        return self.hex

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        # TODO: doesn't really work: only shows red or?
        string = '\e]4;1;rgb:' + self.hex1 + '/' + self.hex2 + '/' + self.hex3 + '\e\\\e[31m██ = ' + self.hex + '\e[m' + '\e]104\a'
        return string.replace("\e", "\033")

# -----------------------------------------------------------------
