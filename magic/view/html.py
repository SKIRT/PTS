#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.view.html Creating FITS viewers in HTML.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools.stringify import stringify_dict, stringify_list
from ...core.tools import html
from ...core.tools import strings

# -----------------------------------------------------------------

# Find the JS9 directory
js9_path = fs.join(fs.home(), "JS9", "js9")
if not fs.is_directory(js9_path): raise IOError("The JS9 installation directory is not found")

# -----------------------------------------------------------------

support_filename = "js9support.css"
main_filename = "js9.css"
prefs_filename = "js9prefs.js"
support_min_filename = "js9support.min.js"
min_filename = "js9.min.js"
plugins_filename = "js9plugins.js"

# -----------------------------------------------------------------

support_filepath = fs.join(js9_path, support_filename)
main_filepath = fs.join(js9_path, main_filename)
prefs_filepath = fs.join(js9_path, prefs_filename)
support_min_filepath = fs.join(js9_path, support_min_filename)
min_filepath = fs.join(js9_path, min_filename)
plugins_filepath = fs.join(js9_path, plugins_filename)

# -----------------------------------------------------------------

meta = """
<meta http-equiv="X-UA-Compatible" content="IE=Edge;chrome=1" > 
<meta name="viewport" content="width=device-width, initial-scale=1">
"""

# -----------------------------------------------------------------

meta_info = dict()
#meta_info["http-equiv"] = "X-UA-Compatible"
#meta_info["content"] = "IE=Edge;chrome=1"
meta_info["viewport"] = "width=device-width, initial-scale=1"

# -----------------------------------------------------------------

scripts = """
<link type="text/css" rel="stylesheet" href="js9support.css">
<link type="text/css" rel="stylesheet" href="js9.css">
<script type="text/javascript" src="js9prefs.js"></script>
<script type="text/javascript" src="js9support.min.js"></script>
<script type="text/javascript" src="js9.min.js"></script>
<script type="text/javascript" src="js9plugins.js"></script>"""

# -----------------------------------------------------------------

css_scripts = [support_filepath, main_filepath]
javascripts = [prefs_filepath, support_min_filepath, min_filepath, plugins_filepath]

# -----------------------------------------------------------------

js9_elements = """
<div class="JS9Menubar"></div>
    <div class="JS9"></div>"""

# -----------------------------------------------------------------

scales = ["log", "linear", "histeq", "power", "sqrt", "squared", "asinh", "sinh"]
colormaps = ["grey", "heat", "cool", "viridis", "magma", "sls", "a", "b", "red", "green", "blue", "inferno", "plasma", "i8", "aips0", "bb", "he", "hsv", "rainbow", "standard", "staircase", "color"]
zooms = ["toFit"]

# -----------------------------------------------------------------

body_settings = dict()
body_settings["onload"] = "init();"

# -----------------------------------------------------------------

class JS9Image(object):

    """
    This function ...
    """

    def __init__(self, name, path, settings=None, display=None):

        """
        This function ...
        :param name
        :param path:
        :param settings:
        :param display:
        """

        self.name = name.replace(" ", "")
        self.path = path
        self.settings = settings
        self.display = display

    # -----------------------------------------------------------------

    def load(self):

        """
        This function ...
        :return:
        """

        string = 'JS9.Load("' + self.path + '"'
        if self.settings is not None: string += ", {" + stringify_dict(self.settings, identity_symbol=":", quote_key=False, quote_character='"')[1] + "}"
        if self.display is not None: string += ', {display:"' + self.display + '"}'

        # End
        string += ');'
        return string

    # -----------------------------------------------------------------

    def preload(self):

        """
        This function ...
        :return:
        """

        string = 'JS9.Preload("' + self.path + '"'
        if self.settings is not None: string += ", {" + stringify_dict(self.settings, identity_symbol=":", quote_key=False, quote_character='"')[1] + "}"
        if self.display is not None: string += ', {display:"' + self.display + '"}'

        # End
        string += ');'
        return string

# -----------------------------------------------------------------

class JS9Loader(object):

    """
    This function ...
    """

    def __init__(self, text, image, button=False):

        """
        This function ...
        :param text:
        :param image:
        :param button:
        """

        # The text for the link
        self.text = text

        # The image
        self.image = image

        # Flag
        self.button = button

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, text, name, path, settings=None, display=None, button=False):

        """
        This function ...
        :param text:
        :param path:
        :param settings:
        :param display:
        :param button:
        :param name:
        :return:
        """

        # Create the image
        image = JS9Image(name, path, settings, display)

        # Return
        return cls(text, image, button=button)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        if self.button:
            buttonid = self.image.name + "Loader"
            load_html = self.image.load()
            return html.button(buttonid, self.text, load_html, quote_character=strings.other_quote_character(self.text, load_html))
        else: return "<a href='javascript:" + self.image.load() + "'>" + self.text + "</a>"

# -----------------------------------------------------------------

class JS9Preloader(object):

    """
    This function ...
    """

    def __init__(self):

        """
        This function ...
        """

        self.images = []

    # -----------------------------------------------------------------

    def add_image(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        self.images.append(image)

    # -----------------------------------------------------------------

    @property
    def has_images(self):

        """
        This function ...
        :return:
        """

        return len(self.images) > 0

    # -----------------------------------------------------------------

    def add_path(self, name, path, settings=None, display=None):

        """
        This function ...
        :param path:
        :param kwargs:
        :param display:
        :return:
        """

        image = JS9Image(name, path, settings, display)
        self.add_image(image)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        string = '<script type="text/javascript">\n'
        string += '    function init()\n'
        string += '    {\n'

        # Preload all images
        for image in self.images: string += image.preload() + "\n"

        string += "    }\n"
        string += "</script>\n"
        return string

# -----------------------------------------------------------------

class JS9Spawner(object):

    """
    This function ...
    """

    def __init__(self, text, image, button=False, menubar=True):

        """
        This function ...
        :param text:
        :param image:
        :param button:
        :param menubar:
        """

        self.text = text
        self.image = image
        self.button = button
        self.menubar = menubar

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, text, name, path, settings=None, button=False, menubar=True):

        """
        This function ...
        :param text:
        :param name:
        :param path:
        :param settings:
        :param button:
        :param menubar:
        :return:
        """

        # Generate display id
        display_id = "JS9" + name.replace(" ", "")

        # Create image
        image = JS9Image(name, path, settings, display_id)

        # Create
        return cls(text, image, button=button, menubar=menubar)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This fcuntion ...
        :return:
        """

        # // make a reasonably unique id for the JS9 divs
        # id = "divJS9_" + ndiv++;

        #code = """
        #  // make up the html with the unique id
        #  //html = sprintf("<p><div class='JS9Menubar' id='%sMenubar'><\/div><div class='JS9' id='%s'><\/div>", id, id);
        #  html = "<p>"

        function_name = "spawn_" + self.image.display

        if self.button:
            buttonid = self.image.name.replace(" ", "") + "Spawner"
            #load_html = self.image.load()
            #return html.button(buttonid, self.text, load_html, quote_character=strings.other_quote_character(self.text, load_html))
            click_code = html.button(buttonid, self.text, function_name + "()")
        #else: return "<a href='javascript:" + self.image.load() + "'>" + self.text + "</a>"
        else: click_code = "<a href='javascript:" + function_name + "()'>" + self.text + "</a>"

        code = '<script type="text/javascript">'
        code += "\n"

        code += 'function ' + function_name + '()'
        code += "\n"
        code += "{"
        code += "\n"

        spawn_code = "<p>"

        menubar_id = self.image.display + "Menubar"

        if self.menubar: spawn_code += "<div class='JS9Menubar' id='" + menubar_id + "'></div>"
        spawn_code += "<div class='JS9' id='" + self.image.display + "'></div>"
        # code += """
        #   // append to end of page
        #   $(html).appendTo($("body"));
        #       // create the new JS9 display, with associated plugins
        #       JS9.AddDivs(id);
        #       // just a standard load to that display
        #       JS9.Load(file, myopts, {display: id});
        #       break;
        #     default:
        #       alert("unknown load type: "+loadtype);"""

        code += 'html = "' + spawn_code + '";'

        spawn_div_name = self.image.name + "placeholder"

        code += "\n"
        #code += '$(html).appendTo($("body"));'
        code += '$(html).appendTo($("#' + spawn_div_name + '"));'
        code += "\n"
        code += "JS9.AddDivs(" + self.image.display + ");"
        code += "\n"
        code += self.image.load()
        code += "\n"

        #code += "\n"
        code += "}"

        code += "\n"

        code += "</script>"

        code += click_code

        code += "<br>"
        code += "<div id='" + spawn_div_name + "'></div>"

        #code += "break;"

        # Return the code
        return code

# -----------------------------------------------------------------

class JS9Viewer(object):

    """
    This function ...
    """

    def __init__(self, id, width=None, height=None):

        """
        This function ...
        :param id:
        :param width:
        :param height:
        """

        self.id = id
        self.width = width
        self.height = height

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        string = '<div class="JS9" id="' + self.id + '"'

        if self.width is not None: string += ' data-width="' + str(self.width) + 'px"'
        if self.height is not None: string += ' data-height="' + str(self.height) + 'px"'

        string += '></div>'
        return string

# -----------------------------------------------------------------

class JS9Menubar(object):

    """
    This function ...
    """

    def __init__(self, id, width=None, background_color=None, displays=None):

        """
        This function ...
        :param id:
        :param width:
        :param background_color:
        """

        self.id = id
        self.width = width
        self.background_color = background_color
        self.displays = displays

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        string = '<div class="JS9Menubar" id="' + self.id + '"'
        if self.width is not None: string += ' data-width="' + str(self.width) + 'px"'
        if self.background_color is not None: string += ' data-backgroundColor="' + self.background_color + '"'
        if self.displays is not None: string += ' data-displays="' + stringify_list(self.displays)[1] + '"'

        # End
        string += "></div>"
        return string

# -----------------------------------------------------------------

class JS9Console(object):

    """
    Thisn function ...
    """

    def __init__(self, id):

        """
        This function ...
        :param id:
        """

        self.id = id

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        string = '<div class="JS9Console" id="' + self.id + '" ></div>'
        return string

# -----------------------------------------------------------------

class JS9Magnifier(object):

    """
    This class ...
    """

    def __init__(self, id):

        """
        This function ...
        :param id:
        """

        self.id = id

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        string = '<div class="JS9Magnifier" id="' + self.id + '"></div>'
        return string

# -----------------------------------------------------------------
