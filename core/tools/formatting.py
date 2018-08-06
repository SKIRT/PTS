#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.formatting Formatting text in the terminal.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import traceback

# Import the relevant PTS classes and modules
from . import stringify
from . import filesystem as fs
from ..basics.log import log
from .sequences import equal_sizes, any_empty, all_none
from .strings import printed_length

# -----------------------------------------------------------------

# SOURCE: http://misc.flogisoft.com/bash/tip_colors_and_formatting

# -----------------------------------------------------------------

# Set
bold = "\033[1m"
dim = "\033[2m"
underlined = "\033[4m"
blink = "\033[5m"
inverted = "\033[7m"
hidden = "\033[8m"

# -----------------------------------------------------------------

# Reset
reset = "\033[0m"
#reset_bold = "\033[21m" # Doesn't work on MacOS
reset_bold = "\033[22m" # 22 resets both dim and bold
reset_dim = "\033[22m" # 22 resets both dim and bold
reset_underlined = "\033[24m"
reset_blink = "\033[25m"
reset_inverted = "\033[27m"
reset_hidden = "\033[28m"

# -----------------------------------------------------------------

# Text colours
default_text = "\033[39m"
black = "\033[30m"
red = "\033[31m"
green = "\033[32m"
yellow = "\033[33m"
blue = "\033[34m"
magenta = "\033[35m"
cyan = "\033[36m"
lightgray = "\033[37m"
darkgray = "\033[90m"
lightred = "\033[91m"
lightgreen = "\033[92m"
lightyellow = "\033[93m"
lightblue = "\033[94m"
lightmagenta = "\033[95m"
lightcyan = "\033[96m"
white = "\033[97m"

# -----------------------------------------------------------------

# Background colours
default_background = "\033[49m"
black_background = "\033[40m"
red_background = "\033[41m"
green_background = "\033[42m"
yellow_background = "\033[43m"
blue_background = "\033[44m"
magenta_background = "\033[45m"
cyan_background = "\033[46m"
lightgray_background = "\033[47m"
darkgray_background = "\033[100m"
lightred_background = "\033[101m"
lightgreen_background = "\033[102m"
lightyellow_background = "\033[103m"
lightblue_background = "\033[104m"
lightmagenta_background = "\033[105m"
lightcyan_background = "\033[106m"
white_background = "\033[107m"

# -----------------------------------------------------------------

def center_text_around(text, ch, length=50):

    """
    This function ...
    :param text:
    :param ch:
    :param length:
    :return:
    """

    spaced_text = '%s' % text
    result = spaced_text.center(length, ch)
    return result

# -----------------------------------------------------------------

def print_empty():

    """
    This function ...
    :return:
    """

    print("")

# -----------------------------------------------------------------

def print_filled(ch, length=50, prefix=""):

    """
    This function ...
    :param ch:
    :param length:
    :param prefix:
    :return:
    """

    result = ch * length
    print(prefix + result)

# -----------------------------------------------------------------

def print_border(ch, length=50, prefix=""):

    """
    This function ...
    :param ch:
    :param length:
    :param prefix:
    :return:
    """

    text = " " * (length - 2)
    result = center_text_around(text, ch, length)
    print(prefix + result)

# -----------------------------------------------------------------

def print_centered_around(text, ch, length=50, prefix=""):

    """
    This function ...
    :param text:
    :param ch:
    :param length:
    :param prefix:
    :return:
    """

    result = center_text_around(text, ch, length)
    print(prefix + result)

# -----------------------------------------------------------------

def print_centered_around_border(text, ch, length=50, prefix=""):

    """
    This function ...
    :param text:
    :param ch:
    :param length:
    :param prefix:
    :return:
    """

    without_border = center_text_around(text, " ", length=length-2)
    result = center_text_around(without_border, ch, length=length)
    print(prefix + result)

# -----------------------------------------------------------------

def print_files_in_path(path):

    """
    This function ...
    :param path:
    :return:
    """

    if not log.is_debug: return
    in_path = fs.files_in_path(path, returns="name", extensions=True)
    if len(in_path) == 0: log.debug("No files in '" + path + "'")
    else:
        log.debug(str(len(in_path)) + " files in '" + path + "':")
        print("")
        print(stringify.stringify_list_fancy(in_path, lines_prefix="  ")[1])
        print("")

# -----------------------------------------------------------------

def print_files_in_list(lst, name, only_name=False):

    """
    This function ...
    :param lst:
    :param name:
    :param only_name:
    :return:
    """

    if not log.is_debug: return
    if len(lst) == 0: log.debug("No files in '" + name + "'")
    else:
        log.debug(str(len(lst)) + " files in '" + name + "':")
        print("")
        if only_name: strings = [fs.name(path) for path in lst]
        else: strings = lst
        print(stringify.stringify_list_fancy(strings, lines_prefix="  ")[1])
        print("")

# -----------------------------------------------------------------

def print_directories_in_path(path):

    """
    This function ...
    :param path:
    :return:
    """

    if not log.is_debug: return
    in_path = fs.directories_in_path(path, returns="name")
    if len(in_path) == 0: log.debug("No directories in '" + path + "'")
    else:
        log.debug(str(len(in_path)) + " directories in '" + path + "':")
        print("")
        print(stringify.stringify_list_fancy(in_path, lines_prefix="  ")[1])
        print("")

# -----------------------------------------------------------------

def print_directories_in_list(lst, name):

    """
    This function ...
    :param lst:
    :param name:
    :return:
    """

    if not log.is_debug: return
    if len(lst) == 0: log.debug("No directories in '" + name + "'")
    else:
        log.debug(str(len(lst)) + " directories in '" + name +"':")
        print("")
        print(stringify.stringify_list_fancy(lst, lines_prefix="  ")[1])
        print("")

# -----------------------------------------------------------------

def print_dictionary(dictionary, bullet="-"):

    """
    This function ...
    :param dictionary:
    :param bullet:
    :return: 
    """

    print("")
    for label in dictionary: print(" " + bullet + " " + bold + label + reset_bold + ": " + stringify.tostr(dictionary[label]))
    print("")

# -----------------------------------------------------------------

#def get_formatted_columns()

def print_table(table, scientific=None, ndecimal_places=3, round=False):

    """
    This function ...
    :param table:
    :param scientific:
    :param ndecimal_places:
    :param round:
    :return:
    """

    from .stringify import tostr

    # Get column names and units
    names = table.colnames
    units = [table[colname].unit for colname in names]
    nrows = len(table)

    # Print in columns
    with print_in_columns() as print_row:

        # Set the column names and units
        column_names = []
        column_units = []

        # Get column names and units
        for name in names: column_names.append(bold + name + reset_bold)
        for unit in units: column_units.append("[" + tostr(unit) + "]" if (unit is not None and unit != "") else "")

        # Show the header
        print_row(*column_names)
        if not all_none(column_units): print_row(*column_units)

        # Initialize columns
        #if path is not None: columns = [[] for _ in column_names]
        #else: columns = None

        # Loop over the rows
        for index in range(nrows):

            # Set parts
            parts = []
            #parts.append(" - ")

            # Add info
            for name, unit in zip(names, units):

                #info = self.get_info_for_simulation(simulation_name)

                if hasattr(table[name], "mask") and table[name].mask[index]: string = "--"
                else:

                    value = table[name][index]
                    #if unit is not None: value = value.to(unit).value
                    string = tostr(value, scientific=scientific, decimal_places=ndecimal_places, round=round)

                parts.append(string)

            # Print the row
            #print_row(*parts, color=color)
            print_row(*parts)

# -----------------------------------------------------------------

def print_sequence(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    print("")
    for item in sequence: print(stringify.tostr(item))
    print("")

# -----------------------------------------------------------------

def get_color_code(color):

    """
    THis function ...
    :param color:
    :return:
    """

    #if color is not None: print(color, list(color))

    if color is None: return default_text
    elif color == "black": return black
    elif color == "red": return red
    elif color == "green": return green
    elif color == "yellow": return yellow
    elif color == "blue": return blue
    elif color == "magenta": return magenta
    elif color == "cyan": return cyan
    elif color == "lightgray": return lightgray
    elif color == "darkgray": return darkgray
    elif color == "lightred": return lightred
    elif color == "lightgreen": return lightgreen
    elif color == "lightyellow": return lightyellow
    elif color == "lightblue": return lightblue
    elif color == "lightmagenta": return lightmagenta
    elif color == "lightcyan": return lightcyan
    elif color == "white": return white
    else: raise ValueError("Invalid color: " + color)

# -----------------------------------------------------------------

def get_background_color_code(color):

    """
    This function ...
    :param color:
    :return:
    """

    if color is None: return default_background
    elif color == "black": return black_background
    elif color == "red": return red_background
    elif color == "green": return green_background
    elif color == "yellow": return yellow_background
    elif color == "blue": return blue_background
    elif color == "magenta": return magenta_background
    elif color == "cyan": return cyan_background
    elif color == "lightgray": return lightgray_background
    elif color == "darkgray": return darkgray_background
    elif color == "lightred": return lightred_background
    elif color == "lightgreen": return lightgreen_background
    elif color == "lightyellow": return lightyellow_background
    elif color == "lightblue": return lightblue_background
    elif color == "lightmagenta": return lightmagenta_background
    elif color == "lightcyan": return lightcyan_background
    elif color == "white": return white_background
    else: raise ValueError("Invalid background color: " + color)

# -----------------------------------------------------------------

def colored_sequence(sequence, colors, delimiter=",", background_colors=None):

    """
    This function ...
    :param sequence:
    :param colors:
    :param delimiter:
    :param background_colors:
    :return:
    """

    # Same color?
    if isinstance(colors, basestring): colors = [colors] * len(sequence)
    elif colors is None: colors = [None] * len(sequence)

    # Same background color?
    if isinstance(background_colors, basestring): background_colors = [background_colors] * len(sequence)
    elif background_colors is None: background_colors = [None] * len(sequence)

    parts = []

    for item, color, background_color in zip(sequence, colors, background_colors):

        # Get codes
        code = get_color_code(color)
        background_code = get_background_color_code(background_color)

        part = background_code + code + stringify.tostr(item) + reset
        parts.append(part)

    # Return the color coded sequence string
    return delimiter.join(parts)

# -----------------------------------------------------------------

def print_columns(*columns, **kwargs):

    """
    This function ...
    :param columns:
    :return:
    """

    delimiter = kwargs.pop("delimiter", "  ")
    indent = kwargs.pop("indent", "")
    tostr_kwargs = kwargs.pop("tostr_kwargs", {})
    colors = kwargs.pop("colors", None)
    bolds = kwargs.pop("bold", None)

    # Check sizes
    if not equal_sizes(*columns): raise ValueError("Columns must have equal lengths")

    if any_empty(*columns): raise ValueError("One or more columns are empty")

    # Convert all to strings
    string_columns = []

    max_lengths = []

    for column in columns:

        new_column = [stringify.tostr(entry, **tostr_kwargs) for entry in column]
        lengths = [printed_length(string) for string in new_column]

        #print(new_column, max(lengths))

        max_lengths.append(max(lengths))
        string_columns.append(new_column)

    ncolumns = len(string_columns)
    nrows = len(columns[0])
    if colors is not None and len(colors) != nrows: raise ValueError("Colors must be defined per row")
    if bolds is not None and len(bolds) != nrows: raise ValueError("Bold flag must be defined per row")

    # Loop over the rows and print them
    for i in range(nrows):

        row = ""

        for j in range(ncolumns): # number of columns

            part = string_columns[j][i]
            row += part

            if j != ncolumns - 1:

                nextra = max_lengths[j] - printed_length(part)
                spaces = " " * nextra + delimiter

                row += spaces

        string = indent + row
        if colors is not None:
            color = colors[i]
            color_code = get_color_code(color)
            string = color_code + string + reset
        if bolds is not None:
            bold_row = bolds[i]
            if bold_row: string = bold + string + reset_bold

        # Print row
        print(string)

# -----------------------------------------------------------------

class print_in_columns(object):

    """
    This function ...
    """

    def __init__(self, ncolumns=None, delimiter=" ", indent="", tostr_kwargs=None):

        """
        This function ...
        :param ncolumns:
        :param delimiter:
        :param indent:
        :param tostr_kwargs:
        """

        # Set columns
        self.columns = None
        if ncolumns is not None: self.initialize_columns(ncolumns)

        self.colors = []
        self.bold = []
        self.delimiter = delimiter
        self.indent = indent
        self.tostr_kwargs = tostr_kwargs if tostr_kwargs is not None else {}

    # -----------------------------------------------------------------

    def initialize_columns(self, ncolumns):

        """
        This function ...
        :param ncolumns:
        :return:
        """

        self.columns = [[] for _ in range(ncolumns)]

    # -----------------------------------------------------------------

    def __enter__(self):

        """
        This function ...
        :return:
        """

        return self

    # -----------------------------------------------------------------

    @property
    def ncolumns(self):

        """
        This function ...
        :return:
        """

        return len(self.columns)

    # -----------------------------------------------------------------

    def __call__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Initialize the columns if the number of columns was not yet defined
        if self.columns is None: self.initialize_columns(len(args))

        # Loop over the columns
        for index in range(self.ncolumns):

            arg = args[index] if len(args) > index else ""
            self.columns[index].append(arg)

        # Set row color and bold flag
        color = kwargs.pop("color", None)
        bold = kwargs.pop("bold", False)
        self.colors.append(color)
        self.bold.append(bold)

    # -----------------------------------------------------------------

    def __exit__(self, exc_type, exc_value, tb):

        """
        This function ...
        :param exc_type:
        :param exc_value:
        :param tb:
        :return:
        """

        #print(exc_type, exc_value, traceback)
        if exc_type is not None:
            traceback.print_tb(tb)
            log.error(str(exc_type(exc_value)))
            return

        # Print
        print_columns(*self.columns, delimiter=self.delimiter, indent=self.indent, tostr_kwargs=self.tostr_kwargs, colors=self.colors, bold=self.bold)

# -----------------------------------------------------------------

class itemize(object):

    """
    This function ...
    """

    def __init__(self, bullet="-", spacing=1):

        """
        This function ...
        :param bullet:
        """

        self.bullet = bullet
        self.spacing = spacing

    # -----------------------------------------------------------------

    def __enter__(self):

        """
        This function ...
        :return:
        """

        print("")
        return self

    # -----------------------------------------------------------------

    @property
    def spaces(self):

        """
        This function ...
        :return:
        """

        return " " * self.spacing

    # -----------------------------------------------------------------

    def __call__(self, string):

        """
        This function ...
        :param args:
        :return:
        """

        print(self.spaces + self.bullet + self.spaces + string)

    # -----------------------------------------------------------------

    def __exit__(self, exc_type, exc_value, traceback):

        """
        This function ...
        :param exc_type:
        :param exc_value:
        :param traceback:
        :return:
        """

        print("")

# -----------------------------------------------------------------
