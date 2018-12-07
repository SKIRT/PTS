#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.colour Contains the Colour class.

# Compatibility between python 2 and 3
from __future__ import print_function

# Import standard modules
from collections import OrderedDict
from matplotlib import colors as mcolors

# Get matplotlib colors
mpl_colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

# Sort colors by hue, saturation, value and name.
#by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name) for name, color in colors.items())
#mpl_color_names = [name for hsv, name in by_hsv]

# -----------------------------------------------------------------

normal = "\033[38;5;%sm"
bold = "\033[1;38;5;%sm"
reset = "\033[0m"

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

# 16 basic HTML color names (HTML 4.01 specification)
predefined = OrderedDict()
predefined["black"] = ("#000000", (0,0,0))
predefined["white"] = ("#FFFFFF", (255,255,255))
predefined["red"] = ("#FF0000",	(255,0,0))
predefined["lime"] = ("#00FF00", (0,255,0))
predefined["blue"] = ("#0000FF", (0,0,255))
predefined["yellow"] = ("#FFFF00", (255,255,0))
predefined["cyan"] = ("#00FFFF", (0,255,255))
predefined["aqua"] = predefined["cyan"]
predefined["magenta"] = ("#FF00FF", (255,0,255))
predefined["fuchsia"] = predefined["magenta"]
predefined["silver"] = ("#C0C0C0", (192,192,192))
predefined["gray"] = ("#808080", (128,128,128))
predefined["maroon"] = ("#800000", (128,0,0))
predefined["olive"] = ("#808000", (128,128,0))
predefined["green"] = ("#008000", (0,128,0))
predefined["purple"] = ("#800080", (128,0,128))
predefined["teal"] = ("#008080", (0,128,128))
predefined["navy"] = ("#000080", (0,0,128))

# X11 colours
predefined["alice_blue"] = ("#F0F8FF", None)
predefined["antique_white"] = ("#FAEBD7", None)
#predefined["aqua"] = ("#00FFFF", None)
predefined["aquamarine"] = ("#7FFFD4", None)
predefined["azure"] = ("#F0FFFF", None)
predefined["beige"] = ("#F5F5DC", None)
predefined["bisque"] = ("#FFE4C4", None)
#predefined["black"] = ("#000000", None)
predefined["blanched_almond"] = ("#FFEBCD", None)
#predefined["blue"] = ("#0000FF", None)
predefined["blue_violet"] = ("#8A2BE2", None)
predefined["brown"] = ("#A52A2A", None)
predefined["burlywood"] = ("#DEB887", None)
predefined["cadet_blue"] = ("#5F9EA0", None)
predefined["chartreuse"] = ("#7FFF00", None)
predefined["chocolate"] = ("#D2691E", None)
predefined["coral"] = ("#FF7F50", None)
predefined["cornflower"] = ("#6495ED", None)
predefined["cornsilk"] = ("#FFF8DC", None)
predefined["crimson"] = ("#DC143C", None)
#predefined["cyan"] = ("#00FFFF", None)
predefined["dark_blue"] = ("#00008B", None)
predefined["dark_cyan"] = ("#008B8B", None)
predefined["dark_goldenrod"] = ("#B8860B", None)
predefined["dark_gray"] = ("#A9A9A9", None)
predefined["dark_green"] = ("#006400", None)
predefined["dark_khaki"] = ("#BDB76B", None)
predefined["dark_magenta"] = ("#8B008B", None)
predefined["dark_olive_green"] = ("#556B2F", None)
predefined["dark_orange"] = ("#FF8C00", None)
predefined["dark_orchid"] = ("#9932CC", None)
predefined["dark_red"] = ("#8B0000", None)
predefined["dark_salmon"] = ("#E9967A", None)
predefined["dark_sea_green"] = ("#8FBC8F", None)
predefined["dark_slate_blue"] = ("#483D8B", None)
predefined["dark_slate_gray"] = ("#2F4F4F", None)
predefined["dark_turquoise"] = ("#00CED1", None)
predefined["dark_violet"] = ("#9400D3", None)
predefined["deep_pink"] = ("#FF1493", None)
predefined["deep_sky_blue"] = ("#00BFFF", None)
predefined["dim_gray"] = ("#696969", None)
predefined["dodger_blue"] = ("#1E90FF", None)
predefined["firebrick"] = ("#B22222", None)
predefined["floral_white"] = ("#FFFAF0", None)
predefined["forest_green"] = ("#228B22", None)
#predefined["fuchsia"] = ("#FF00FF", None)
predefined["gainsboro"] = ("#DCDCDC", None)
predefined["ghost_white"] = ("#F8F8FF", None)
predefined["gold"] = ("#FFD700", None)
predefined["goldenrod"] = ("#DAA520", None)
#predefined["gray"] = ("#BEBEBE", None)
predefined["web_gray"] = ("#808080", None)
#predefined["green"] = ("#00FF00", None)
predefined["web_green"] = ("#008000", None)
predefined["green_yellow"] = ("#ADFF2F", None)
predefined["honeydew"] = ("#F0FFF0", None)
predefined["hot_pink"] = ("#FF69B4", None)
predefined["indian_red"] = ("#CD5C5C", None)
predefined["indigo"] = ("#4B0082", None)
predefined["ivory"] = ("#FFFFF0", None)
predefined["khaki"] = ("#F0E68C", None)
predefined["lavender"] = ("#E6E6FA", None)
predefined["lavender_blush"] = ("#FFF0F5", None)
predefined["lawn_green"] = ("#7CFC00", None)
predefined["lemon_chiffon"] = ("#FFFACD", None)
predefined["light_blue"] = ("#ADD8E6", None)
predefined["light_coral"] = ("#F08080", None)
predefined["light_cyan"] = ("#E0FFFF", None)
predefined["light_goldenrod"] = ("#FAFAD2", None)
predefined["light_gray"] = ("#D3D3D3", None)
predefined["light_green"] = ("#90EE90", None)
predefined["light_pink"] = ("#FFB6C1", None)
predefined["light_salmon"] = ("#FFA07A", None)
predefined["light_sea_green"] = ("#20B2AA", None)
predefined["light_sky_blue"] = ("#87CEFA", None)
predefined["light_slate_gray"] = ("#778899", None)
predefined["light_steel_blue"] = ("#B0C4DE", None)
predefined["light_yellow"] = ("#FFFFE0", None)
#predefined["lime"] = ("#00FF00", None)
predefined["lime_green"] = ("#32CD32", None)
predefined["linen"] = ("#FAF0E6", None)
#predefined["magenta"] = ("#FF00FF", None)
#predefined["maroon"] = ("#B03060", None)
predefined["web_maroon"] = ("#7F0000", None)
predefined["medium_aquamarine"] = ("#66CDAA", None)
predefined["medium_blue"] = ("#0000CD", None)
predefined["medium_orchid"] = ("#BA55D3", None)
predefined["medium_purple"] = ("#9370DB", None)
predefined["medium_sea_green"] = ("#3CB371", None)
predefined["medium_slate_blue"] = ("#7B68EE", None)
predefined["medium_spring_green"] = ("#00FA9A", None)
predefined["medium_turquoise"] = ("#48D1CC", None)
predefined["medium_violet_red"] = ("#C71585", None)
predefined["midnight_blue"] = ("#191970", None)
predefined["mint_cream"] = ("#F5FFFA", None)
predefined["misty_rose"] = ("#FFE4E1", None)
predefined["moccasin"] = ("#FFE4B5", None)
predefined["navajo_white"] = ("#FFDEAD", None)
predefined["navy_blue"] = ("#000080", None)
predefined["old_lace"] = ("#FDF5E6", None)
#predefined["olive"] = ("#808000", None)
predefined["olive_drab"] = ("#6B8E23", None)
predefined["orange"] = ("#FFA500", None)
predefined["orange_red"] = ("#FF4500", None)
predefined["orchid"] = ("#DA70D6", None)
predefined["pale_goldenrod"] = ("#EEE8AA", None)
predefined["pale_green"] = ("#98FB98", None)
predefined["pale_turquoise"] = ("#AFEEEE", None)
predefined["pale_violet_red"] = ("#DB7093", None)
predefined["papaya_whip"] = ("#FFEFD5", None)
predefined["peach_puff"] = ("#FFDAB9", None)
predefined["peru"] = ("#CD853F", None)
predefined["pink"] = ("#FFC0CB", None)
predefined["plum"] = ("#DDA0DD", None)
predefined["powder_blue"] = ("#B0E0E6", None)
#predefined["purple"] = ("#A020F0", None)
predefined["web_purple"] = ("#7F007F", None)
predefined["rebecca_purple"] = ("#663399", None)
#predefined["red"] = ("#FF0000", None)
predefined["rosy_brown"] = ("#BC8F8F", None)
predefined["royal_blue"] = ("#4169E1", None)
predefined["saddle_brown"] = ("#8B4513", None)
predefined["salmon"] = ("#FA8072", None)
predefined["sandy_brown"] = ("#F4A460", None)
predefined["sea_green"] = ("#2E8B57", None)
predefined["seashell"] = ("#FFF5EE", None)
predefined["sienna"] = ("#A0522D", None)
#predefined["silver"] = ("#C0C0C0", None)
predefined["sky_blue"] = ("#87CEEB", None)
predefined["slate_blue"] = ("#6A5ACD", None)
predefined["slate_gray"] = ("#708090", None)
predefined["snow"] = ("#FFFAFA", None)
predefined["spring_green"] = ("#00FF7F", None)
predefined["steel_blue"] = ("#4682B4", None)
predefined["tan"] = ("#D2B48C", None)
#predefined["teal"] = ("#008080", None)
predefined["thistle"] = ("#D8BFD8", None)
predefined["tomato"] = ("#FF6347", None)
predefined["turquoise"] = ("#40E0D0", None)
predefined["violet"] = ("#EE82EE", None)
predefined["wheat"] = ("#F5DEB3", None)
#predefined["white"] = ("#FFFFFF", None)
predefined["white_smoke"] = ("#F5F5F5", None)
#predefined["yellow"] = ("#FFFF00", None)
predefined["yellow_green"] = ("#9ACD3", None)

# -----------------------------------------------------------------

def get_colour_names():

    """
    Thisf unction ...
    :return:
    """

    return predefined.keys()

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
        elif name in mpl_colors: return cls.from_hex(mpl_colors[name])
        else: raise ValueError("Colour '" + name + "' not recognized")

    # -----------------------------------------------------------------

    @property
    def rgb(self):
        return self.red, self.green, self.blue

    # -----------------------------------------------------------------

    @property
    def hex(self):
        return rgb_to_hex([self.red, self.green, self.blue])

    # -----------------------------------------------------------------

    @property
    def hex1(self):
        return self.hex[1:3]

    # -----------------------------------------------------------------

    @property
    def hex2(self):
        return self.hex[3:5]

    # -----------------------------------------------------------------

    @property
    def hex3(self):
        return self.hex[5:]

    # -----------------------------------------------------------------

    @property
    def hex_slashes(self):
        return self.hex1 + "/" + self.hex2 + "/" + self.hex3

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        return self.hex

    # -----------------------------------------------------------------

    # def __repr__(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # ## TODO: doesn't really work: only shows red or?
    #     # #string = '\e]4;1;rgb:' + self.hex1 + '/' + self.hex2 + '/' + self.hex3 + '\e\\\e[31m██ = ' + self.hex + '\e[m' + '\e]104\a'
    #     # #return string.replace("\e", "\033")
    #     #
    #     # #string = normal + self.hex1 + "/" + self.hex2 + "/" + self.hex3 + "██ = " + self.hex + reset
    #     #
    #     # #i = 0
    #     # #string = (normal + "%s" + reset) % (i, self.hex_slashes)
    #     # string = (normal + "%s" + reset) % self.hex_slashes
    #     # return string

# -----------------------------------------------------------------
