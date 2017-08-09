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
js9_filename = "js9.js"
min_filename = "js9.min.js"
plugins_filename = "js9plugins.js"

# -----------------------------------------------------------------

support_filepath = fs.join(js9_path, support_filename)
main_filepath = fs.join(js9_path, main_filename)
prefs_filepath = fs.join(js9_path, prefs_filename)
support_min_filepath = fs.join(js9_path, support_min_filename)
min_filepath = fs.join(js9_path, min_filename)
#min_filepath = fs.join(js9_path, js9_filename)
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

    def load(self, regions=None):

        """
        This function ...
        :param regions:
        :return:
        """

        string = 'JS9.Load("' + self.path + '"'
        if self.settings is not None: string += ", {" + stringify_dict(self.settings, identity_symbol=":", quote_key=False, quote_character='"')[1] + "}"
        if self.display is not None: string += ', {display:"' + self.display + '"}'

        # End
        string += ');'

        # Add regions
        if regions is not None:
            string += "\n"
            #string += 'function sleep(ms) { return new Promise(resolve => setTimeout(resolve, ms));}'
            #string += "\n"
            #string += "sleep(5000);"
            string += "var region_id = JS9.AddRegions('" + regions + "'"
            #string += "setTimeout(JS9.AddRegions, 5000, '" + regions + "'"
            #if self.display is not None: string += ', {display:"' + self.display + '"}'
            string += ');'
            #string += "window.alert(region_id);"

        # Return
        return string

    # -----------------------------------------------------------------

    def preload(self, regions=None):

        """
        This function ...
        :param regions:
        :return:
        """

        string = 'JS9.Preload("' + self.path + '"'
        if self.settings is not None: string += ", {" + stringify_dict(self.settings, identity_symbol=":", quote_key=False, quote_character='"')[1] + "}"
        if self.display is not None: string += ', {display:"' + self.display + '"}'

        # End
        string += ');'

        # Add regions
        if regions is not None:
            string += "\n"
           # string += 'function sleep(ms) { return new Promise(resolve => setTimeout(resolve, ms));}'
           # string += "\n"
            #string += "sleep(5000);"
            string += "var region_id = JS9.AddRegions('" + regions + "'"
            #string += "setTimeout(JS9.AddRegions, 5000, '" + regions + "'"
            #if self.display is not None: string += ', {display:"' + self.display + '"}'
            string += ');'
            #string += "window.alert(region_id);"

        # Return
        return string

# -----------------------------------------------------------------

class JS9Loader(object):

    """
    This function ...
    """

    def __init__(self, text, image, button=False, regions=None):

        """
        This function ...
        :param text:
        :param image:
        :param button:
        :param regions:
        """

        # The text for the link
        self.text = text

        # The image
        self.image = image

        # Flag
        self.button = button

        # Regions
        self.regions = regions

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, text, name, path, settings=None, display=None, button=False, regions=None):

        """
        This function ...
        :param text:
        :param path:
        :param settings:
        :param display:
        :param button:
        :param name:
        :param regions:
        :return:
        """

        # Create the image
        image = JS9Image(name, path, settings, display)

        # Return
        return cls(text, image, button=button, regions=regions)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        if self.button:
            buttonid = self.image.name + "Loader"
            load_html = self.image.load(regions=self.regions)
            return html.button(buttonid, self.text, load_html, quote_character=strings.other_quote_character(self.text, load_html))
        else: return "<a href='javascript:" + self.image.load(regions=self.regions) + "'>" + self.text + "</a>"

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
        self.regions = []

    # -----------------------------------------------------------------

    def add_image(self, image, regions=None):

        """
        This function ...
        :param image:
        :param regions:
        :return:
        """

        self.images.append(image)
        self.regions.append(regions)

    # -----------------------------------------------------------------

    @property
    def has_images(self):

        """
        This function ...
        :return:
        """

        return len(self.images) > 0

    # -----------------------------------------------------------------

    def add_path(self, name, path, settings=None, display=None, regions=None):

        """
        This function ...
        :param path:
        :param kwargs:
        :param display:
        :param regions:
        :return:
        """

        image = JS9Image(name, path, settings, display)
        self.add_image(image, regions=regions)
        return image

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
        for image, regions in zip(self.images, self.regions): string += image.preload(regions=regions) + "\n"

        string += "    }\n"
        string += "</script>\n"
        return string

# -----------------------------------------------------------------

class JS9Spawner(object):

    """
    This function ...
    """

    def __init__(self, text, image, button=False, menubar=True, colorbar=False, width=None, height=None, regions=None,
                 add_placeholder=True, background_color="white", replace=False):

        """
        This function ...
        :param text:
        :param image:
        :param button:
        :param menubar:
        :param colorbar:
        :param regions:
        :param background_color:
        """

        self.text = text
        self.image = image
        self.button = button
        self.menubar = menubar
        self.colorbar = colorbar
        self.width = width
        self.height = height
        self.regions = regions
        self.add_placeholder = add_placeholder
        self.background_color = background_color
        self.replace = replace

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, text, name, path, settings=None, button=False, menubar=True, colorbar=False, width=None,
                  height=None, regions=None, add_placeholder=True, background_color="white", replace=False):

        """
        This function ...
        :param text:
        :param name:
        :param path:
        :param settings:
        :param button:
        :param menubar:
        :param colorbar:
        :param width:
        :param height:
        :param regions:
        :param add_placeholder:
        :param background_color:
        :param replace:
        :return:
        """

        # Generate display id
        display_id = "JS9" + name.replace(" ", "")

        # Create image
        image = JS9Image(name, path, settings, display_id)

        # Create
        return cls(text, image, button=button, menubar=menubar, colorbar=colorbar, width=width, height=height,
                   regions=regions, add_placeholder=add_placeholder, background_color=background_color, replace=replace)

    # -----------------------------------------------------------------

    @property
    def display_id(self):

        """
        This function ...
        :return:
        """

        return self.image.display

    # -----------------------------------------------------------------

    @property
    def spawn_div_name(self):

        """
        This function ...
        :return:
        """

        return self.image.name + "placeholder"

    # -----------------------------------------------------------------

    @property
    def placeholder_start(self):

        """
        This function ...
        :return:
        """

        return "<div id='" + self.spawn_div_name + "'>"

    # -----------------------------------------------------------------

    @property
    def placeholder_end(self):

        """
        This function ...
        :return:
        """

        return "</div>"

    # -----------------------------------------------------------------

    @property
    def placeholder(self):

        """
        This function ...
        :return:
        """

        return self.placeholder_start + self.placeholder_end

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
        function_call = function_name + "()"

        #function_call += ";load_region_" + self.image.name.replace(" ", "") + "()"

        if self.button:
            buttonid = self.image.name.replace(" ", "") + "Spawner"
            #load_html = self.image.load()
            #return html.button(buttonid, self.text, load_html, quote_character=strings.other_quote_character(self.text, load_html))
            click_code = html.button(buttonid, self.text, function_call)
        #else: return "<a href='javascript:" + self.image.load() + "'>" + self.text + "</a>"
        else: click_code = "<a href='javascript:" + function_call + "'>" + self.text + "</a>"

        code = '<script type="text/javascript">'
        code += "\n"

        code += 'function ' + function_name + '()'
        code += "\n"
        code += "{"
        code += "\n"

        # Make spawn code
        spawn_code = make_spawn_code(self.image, menubar=self.menubar, colorbar=self.colorbar, width=self.width,
                                     background_color=self.background_color)

        #print(spawn_code)

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

        # Add html code
        if self.replace: code += replace_div(self.spawn_div_name, spawn_code)
        else: code += add_to_div(self.spawn_div_name, spawn_code)

        # Add div ID to JS9
        code += "JS9.AddDivs('" + self.image.display + "');"

        code += "\n"

        # Add load image code
        code += self.image.load(regions=self.regions)
        #code += self.image.load()
        code += "\n"

        #code += "\n"
        code += "}"

        code += "\n"

        # if self.regions is not None:
        #     code += "function load_region_" + self.image.name.replace(" ", "") + "()"
        #     code += "\n"
        #     code += "{"
        #     code += "var region_id = JS9.AddRegions('" + self.regions + "');"
        #     code += "\n"
        #     code += "}"

        code += "</script>"

        if self.replace: code += self.placeholder_start

        code += click_code

        if self.replace: code += self.placeholder_end

        if self.add_placeholder and not self.replace:
            code += "<br>"
            code += self.placeholder

        #code += "break;"

        # Return the code
        return code

# -----------------------------------------------------------------

class JS9Viewer(object):

    """
    This function ...
    """

    def __init__(self, id, width=None, height=None, resize=True, scrolling=True):

        """
        This function ...
        :param id:
        :param width:
        :param height:
        :param resize:
        :param scrolling:
        """

        self.id = id
        self.width = width
        self.height = height
        self.resize = resize
        self.scrolling = scrolling

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        string = '<div class="JS9" id="' + self.id + '"'

        if self.width is not None: string += ' data-width="' + str(self.width) + 'px"'
        if self.height is not None: string += ' data-height="' + str(self.height) + 'px"'

        string += ' resize=' + str(int(self.resize))
        string += ' scrolling=' + str(int(self.scrolling))

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

class JS9Colorbar(object):

    """
    This class ...
    """

    def __init__(self, id, width=None):

        """
        This function ...
        :param id:
        :param width:
        """

        self.id = id
        self.width = width

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        string = '<div class="JS9Colorbar" id="' + self.id + '"'

        if self.width is not None: string += ' data-width="' + str(self.width) + 'px"'

        string += "></div>"
        return string

# -----------------------------------------------------------------

class JS9Panner(object):

    """
    This function ...
    """

    def __init__(self, id, width=None, height=None):

        """
        This function ...
        :param id:
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

        string = '<div class="JS9Panner" id="' + self.id + '"'

        if self.width is not None: string += ' data-width="' + str(self.width) + 'px"'
        if self.height is not None: string += ' data-height="' + str(self.height) + 'px"'

        string += "></div>"
        return string

# -----------------------------------------------------------------

class JS9Window(object):

    """
    This function ...
    """

    def __init__(self, name, width=None, height=None, background_color="white", menubar=True, colorbar=False, menubar_position="top", colorbar_position="bottom", resize=True):

        """
        This function ...
        :param name:
        :param width:
        :param height:
        :param background_color:
        :param menubar:
        :param colorbar:
        :param resize:
        """

        # Create view
        self.view = JS9Viewer(name, width=width, height=height, resize=resize)

        # Set names
        # Set menu name
        # menu_name = display_name + "_menu"
        menu_name = name + "Menubar"
        # displays = [display_name]
        displays = None
        bar_name = name + "Colorbar"

        # Add menu bar
        if menubar: self.menu = JS9Menubar(menu_name, displays=displays, width=width, background_color=background_color)
        else: self.menu = None

        # Add color bar
        if colorbar: self.color = JS9Colorbar(bar_name, width=width)
        else: self.color = None

        # Positions
        self.menubar_position = menubar_position
        self.colorbar_position = colorbar_position

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        string = ""

        # Add menu bar (if requested)
        # if name in self.menus: string += str(self.menus[name])

        # string += newline

        # Add view (if not dynamic)
        # if name in self.views: string += str(self.views[name])

        if self.menu is not None and self.menubar_position == "top":
            string += str(self.menu)
            string += html.newline

        if self.color is not None and self.colorbar_position == "top":
            string += str(self.color)
            string += html.newline

        string += str(self.view)
        #string += html.newline

        if self.color is not None and self.colorbar_position == "bottom":
            string += html.newline
            string += str(self.color)

        if self.menu is not None and self.menubar_position == "bottom":
            string += html.newline
            string += str(self.menu)

        # Return the html string
        return string

# -----------------------------------------------------------------

def make_load_region(region, display=None, changeable=True, movable=True, resizable=True, rotatable=True,
                     removable=True, zoomable=True, lock_x=False, lock_y=False, lock_rotation=False,
                     quote_character='"'):

    """
    This function ...
    :param region:
    :param display:
    :return:
    """

    string = ""

    properties = dict()
    properties["changeable"] = str(changeable).lower()
    properties["movable"] = str(movable).lower()
    properties["resizable"] = str(resizable).lower()
    properties["rotatable"] = str(rotatable).lower()
    properties["removable"] = str(removable).lower()
    properties["zoomable"] = str(zoomable).lower()
    properties["lockMovementX"] = str(lock_x).lower()
    properties["lockMovementY"] = str(lock_y).lower()
    properties["lockRotation"] = str(lock_rotation).lower()

    string += "var region_id = JS9.AddRegions('" + str(region) + "'"

    string += ', {' + stringify_dict(properties, identity_symbol=":", quote_key=False, quote_value=False, quote_character=quote_character)[1] + '}'

    if display is not None:
        if quote_character == '"': string += ', {display:"' + display + '"}'
        elif quote_character == "'": string += ", {display:'" + display + "'}"
        else: raise ValueError("Invalid quote character: " + quote_character)
    string += ');'
    string += "\n"

    #string += "window.alert(region_id);"

    return string

# -----------------------------------------------------------------

def make_load_region_function(name, region, display=None):

    """
    This function ...
    :param name:
    :param region:
    :param display:
    :return:
    """

    if " " in name: raise ValueError("Name cannot contain spaces")
    string = "function " + name + "()"
    string += "\n"
    string += "{"

    for line in make_load_region(region, display=display).split("\n"):
        string += "    " + line + "\n"

    string += "}"
    return string

# -----------------------------------------------------------------

def make_load_regions(regions, display=None):

    """
    This function ...
    :param regions:
    :param display:
    :return:
    """

    string = ""

    for region in regions:

        string += "JS9.AddRegions('" + str(region) + "'"
        if display is not None: string += ', {display:"' + display + '"}'
        string += ');'
        string += "\n"

    return string

# -----------------------------------------------------------------

def make_load_regions_function(name, regions, display=None):

    """
    This function ...
    :param name:
    :param regions:
    :param display:
    :return:
    """

    if " " in name: raise ValueError("Name cannot contain spaces")
    string = "function " + name + "()"
    string += "\n"
    string += "{"

    for line in make_load_regions(regions, display=display).split("\n"):
        string += "    " + line + "\n"

    string += "}"
    return string

# -----------------------------------------------------------------

# Start for run analysis on change
"""
<script type="text/javascript">
    var aname, im;
    var lastim, lastreg;
    var ncall = 0;
    // this is the callback for all region changes
    JS9.Regions.opts.onchange = "runMyAnalysis";
    // called when the function changes to redo the last display
    function redo(){
      if( lastim && lastreg ){
        runMyAnalysis(lastim, lastreg);
      }
    }"""

# -----------------------------------------------------------------

def make_synchronize_regions(indicator_id, display_ids, ellipses, ndecimals=3):

    """
    This function ...
    :param indicator_id:
    :param display_ids:
    :param ellipses: default, original ellipses:
    :param ndecimals:
    :return:
    """

    function_name = "synchronizeRegions"

    code = ""

    code += "var aname, im;\n"
    code += "var lastim, lastreg;\n"
    code += "var ncall = 0;\n"

    x_radii = dict()
    y_radii = dict()
    for display_id in ellipses:
        x_radii[display_id] = ellipses[display_id].radius.x
        y_radii[display_id] = ellipses[display_id].radius.y

    code += "var x_radii = {" + stringify_dict(x_radii, quote_key=False, quote_value=False, identity_symbol=":")[1] + "};\n"
    code += "var y_radii = {" + stringify_dict(y_radii, quote_key=False, quote_value=False, identity_symbol=":")[1] + "};\n"

    code += 'JS9.Regions.opts.onchange = "' + function_name + '";'
    code += "\n"

    #code += """// called when the function changes to redo the last display
    code += "function redo()\n"
    code += "{\n"
    code += "    if( lastim && lastreg )\n"
    code += "    {\n"
    code += "        " + function_name + "(lastim, lastreg);\n"
    code += "    }\n"
    code += "}\n"

    code += "\n\n"

    code += html.round_to_decimals_function

    code += "\n\n"

    existing_function_code = ""
    existing_function_code += "function isExistingDisplayWithImage(displayid)\n"
    existing_function_code += "{\n"
    existing_function_code += "    return JS9.IsDisplay(displayid);\n"
    existing_function_code += "}\n"

    code += existing_function_code

    code += "\n"
    code += 'function ' + function_name + '(im, xreg)\n'
    code += "{\n"

    code += "    lastim = im;\n"
    code += "    lastreg = xreg;\n"

    code += "    var x_radius = lastreg.r1;\n"
    code += "    var y_radius = lastreg.r2;\n"

    #code += "    window.alert(x_radius);\n"
    #code += "    window.alert(y_radius);\n"

    #code += "    window.alert(lastim.display.id);\n"

    #new_factor = 2.0

    code += "    var x_factor = x_radius / x_radii[lastim.display.id];\n"
    code += "    var y_factor = y_radius / y_radii[lastim.display.id];\n"

    #code += "   window.alert('factor = ' + String(x_factor));\n"
    #code += "   window.alert('factor = ' + String(y_factor));\n"

    # Check whether values are close
    code += "    var tolerance = 1e4;\n"
    code += "    var relDifference = x_factor / y_factor - 1.;\n"
    code += "    if(relDifference > tolerance)\n"
    code += "    {\n"
    code += "        throw 'The relative difference is ' + String(relDifference);\n"
    code += "    }\n"
    code += "    var factor = x_factor;\n"

    #code += "   window.alert('factor = ' + String(factor));\n"

    code += "    var rounded = roundToDecimals(factor, " + str(ndecimals) + ");\n"

    #code += "    window.alert(rounded);\n"

    code += "    var new_text = 'Factor: ' + String(rounded);\n"

    #text = "Factor: " + str(new_factor)
    #code += '$("div#' + indicator_id + '").text("' + text + '");'
    code += '    $("div#' + indicator_id + '").text(new_text);\n'

    # LOOP OVER THE OTHER IMAGES
    #code += '    for(var i = 0; i < ' + len(display_ids) + ' ; i++)'

    code += "    var displayIds = " + str(display_ids) + ";\n"

    #code += "    window.alert(displayIds);\n"

    # for display_id in display_ids:
    #
    #     im_temp_name = display_id.lower() + "image"
    #     code += "    try\n    {\n"
    #     code += "        var " + im_temp_name + " = JS9.GetImage({display: '" + display_id + "'});\n"
    #     code += "        window.alert('" + display_id + " image is found');\n"
    #     code += "    }\n"
    #     code += "    catch(err)\n{\n"
    #     code += "        window.alert('" + display_id + " image is NOT found');"
    #     code += "\n    }"

    #code += "    window.alert(displayIds.length);\n"

    # existing_function_code = ""
    # existing_function_code += "function isExistingDisplay(display_id)\n"
    # existing_function_code += "{\n"
    # existing_function_code += "    try\n"
    # existing_function_code += "    {\n"
    # existing_function_code += "        var im_temp_name = JS9.GetImage({display: display_id});\n"
    # existing_function_code += "        return true;\n"
    # existing_function_code += "    }\n"
    # existing_function_code += "    catch(err)\n"
    # existing_function_code += "    {\n"
    # existing_function_code += "        return false;\n"
    # existing_function_code += "    }\n"
    # existing_function_code += "}\n"

    #code += existing_function_code
    #code += existing_function_code

    code += "\n"

    code += "    for(var i = 0; i < displayIds.length ; i++)\n"
    code += "    {\n"
    code += "         var displayid = displayIds[i];\n"
    #code += "        window.alert(i);\n"
    #code += "        var isCurrent = displayIDs[i] == lastim.display.id;\n"
    #code += "        window.alert(isCurrent);\n"
    #code += "        if (isCurrent == true) { window.alert(displayIDs[i] + ' = current'); }\n"
    #code += "        else { window.alert(displayIDs[i] + ' = not current'); }\n"

    #code += "        window.alert(displayIds[i] == lastim.display.id);\n"

    #code += "        if(displayid == lastim.display.id)\n"
    #code += "        if(displayid != lastim.display.id)\n"
    code += "        var exists = isExistingDisplayWithImage(displayid);\n"
    code += "        window.alert('Display ' + displayid + ': ' + String(exists));\n"
    code += "        if(displayid != lastim.display.id && exists)\n"
    #code += "        {\n"
    #code += "            window.alert(displayIds[i] + ' = current');\n"
    #code += "        }\n"
    #code += "        else\n"
    code += "        {\n"
    #code += "            window.alert(displayIds[i] + ' = not current');\n"

    code += "            var x_radius_i = x_radii[displayid];\n"
    code += "            var y_radius_i = y_radii[displayid];\n"

    code += "            x_radius_i = x_radius_i * factor;\n"
    code += "            y_radius_i = y_radius_i * factor;\n"

    code += "            // Get the region\n"

    code += "            var xreg = JS9.GetRegions('all')[0];\n"
    code += "            var regid = xreg.id;\n"

    code += "            JS9.ChangeRegions(regid, {'r1': x_radius_i, 'r2': y_radius_i}, {display: displayid});\n"

    code += "        }\n"
    code += "    }\n"
    code += "\n}"

    # Return the code
    return code

# -----------------------------------------------------------------

def make_spawn_code(image, menubar=True, colorbar=False, width=None, background_color="white"):

    """
    This function ...
    :param image:
    :param menubar:
    :param colorbar:
    :param width:
    :param background_color:
    :return:
    """

    spawn_code = "<p>"

    menubar_id = image.display + "Menubar"
    colorbar_id = image.display + "Colorbar"

    if menubar:

        menubar = JS9Menubar(menubar_id, width=width, background_color=background_color)
        # spawn_code += "<div class='JS9Menubar' id='" + menubar_id + "'></div>"
        spawn_code += str(menubar)

    if colorbar:

        colorbar = JS9Colorbar(colorbar_id, width=width)
        # spawn_code += "<div class='JS9Colorbar' id = '" + colorbar_id + "'></div>"
        spawn_code += str(colorbar)

    spawn_code += '<div class="JS9" id="' + image.display + '"></div>'
    return spawn_code

# -----------------------------------------------------------------

def add_to_div(div_id, text):

    """
    This function ...
    :param div_id:
    :param text:
    :return:
    """

    code = ""
    code += 'html = "' + strings.make_single_quoted(text) + '";'
    code += "\n"
    # code += '$(html).appendTo($("body"));'
    code += '$(html).appendTo($("#' + div_id + '"));'
    code += "\n"

    # Return the code
    return code

# -----------------------------------------------------------------

def replace_div(div_id, text):

    """
    This function ...
    :param div_id:
    :param text:
    :return:
    """

    code = ""
    code += 'html = "' + strings.make_single_quoted(text) + '";'
    code += "\n"
    #code += '$("#' + div_id + '").update(html);'
    code += '$("#' + div_id + '").html(html);'
    code += "\n"

    # Return the code
    return code

# -----------------------------------------------------------------
