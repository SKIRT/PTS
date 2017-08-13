#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.status Contains the ModelingStatus class.

# ORIGINAL PACKAGE:
# simpletable.py - v0.1 2014-07-31 Matheus Vieira Portela

# This module provides simple classes and interfaces to generate simple HTML
# tables based on Python native types, such as lists.

# Author's website: http://matheusvportela.wordpress.com/

#__version__ = '0.3'
#__date__    = '2014-08-20'
#__author__  = 'Matheus Vieira Portela'

### CHANGES ###
# 2014-07-31: v0.1 MVP:
#   - First version
# 2014-08-05: v0.2 MVP:
#   - Method for defining header rows
#   - SimpleTable method to create a SimpleTable from lists
#   - Method to create a table from a simple list of elements and a column size
# 2014-08-20: v0.3 MVP:
#   - Enable SimplePage to accept a list of tables
#   - Enable SimplePage to iterate over its tables

### REFERENCES ###
# Decalage HTML.py module: http://www.decalage.info/python/html

# EXAMPLE:

# test_data = [str(x) for x in range(20)]
# formatted_data = simpletable.fit_data_to_columns(test_data, 5)
# table = simpletable.SimpleTable(formatted_data)
# html_page = simpletable.HTMLPage(table)
# html_page.save("test_page.html")

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .stringify import tostr, stringify_dict
from . import numbers
from . import filesystem as fs
from . import types
from . import time
from . import strings
from ..basics.table import SmartTable

# -----------------------------------------------------------------

line = "<hr>"
newline = "<br>"
item = "<li>"

# -----------------------------------------------------------------

def make_line(css_class=None):

    """
    This function ...
    :return: 
    """

    if css_class is None: return line
    else: return "<hr class='" + css_class + "'>"

# -----------------------------------------------------------------

def make_newline(css_class=None):

    """
    This function ...
    :param css_class:
    :return:
    """

    if css_class is None: return line
    else: return "<br class='" + css_class + "'>"

# -----------------------------------------------------------------

def make_item(css_class=None):

    """
    This function ...
    :param css_class:
    :return:
    """

    if css_class is None: return item
    else: return "<li class='" + css_class + "'>"

# -----------------------------------------------------------------

underline_template = "<span style='text-decoration: underline;'>{text}</span>"
bold_template = "<span style='font-weight:bold'>{text}</span>"
fontsize_template = "<span style='font-size:{size}px'>{text}</span>"
small_template = "<small>{text}</small>"
center_template = "<div style='text-align:center;'>{text}</div>"

# -----------------------------------------------------------------

def make_css(contents):

    """
    This function ...
    :return:
    """

    #code = """<style>
    #body {background-color: powderblue;}
    #h1   {color: blue;}
    #p    {color: red;}
    #</style>"""

    # #lines.append('    <style type="text/css">\n%s\n</style>' % self.css)

    # Add code
    code = "<style type='text/css'>\n"
    code += contents

    # Return
    code += "</style>\n"
    return code

# -----------------------------------------------------------------

def make_body_settings(settings):

    """
    This function ...
    :param settings:
    :return:
    """

    # body
    # {
    #     background - color: powderblue;
    # }

    contents = "body\n"
    contents += "{\n"

    for name in settings:
        value = settings[name]
        contents += "    " + name + ": " + tostr(value) + ";\n"

    contents += "}\n"

    # Make the css code and return
    #return make_css(contents)
    return contents

# -----------------------------------------------------------------

def make_page_width(width):

    """
    This function ...
    :param width:
    :return:
    """

    settings = dict()
    settings["max-width"] = str(width) + "px"
    css_settings = make_body_settings(settings)
    return css_settings

# -----------------------------------------------------------------

def center(html):

    """
    This function ...
    :param html:
    :return:
    """

    return center_template.format(text=html)

# -----------------------------------------------------------------

def updated_footing():

    """
    This function ...
    :return:
    """

    # Generate footing
    footing = center(small_template.format(text="Last updated on " + time.pretty_time()))
    return footing

# -----------------------------------------------------------------

def button(id, text, onclick, quote_character='"'):

    """
    This function ...
    :param id:
    :param text:
    :param onclick:
    :param quote_character:
    :return:
    """

    return '<input id=' + quote_character + id + quote_character + ' onclick=' + quote_character + onclick + quote_character + ' type=' + quote_character + 'button' + quote_character + ' value=' + quote_character + text + quote_character + '>'

# -----------------------------------------------------------------

link_stylesheet_header_template = '<link rel="stylesheet" type="text/css" href="{url}">'

# -----------------------------------------------------------------

page_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{title}</title>
{head}
</head>
<body>
<div class="{style}">
{body}
</div>
</body>
</html>
"""

# -----------------------------------------------------------------

def ordered_list(items, css_class=None):

    """
    This function ...
    :param items:
    :param css_class:
    :return:
    """

    html = ""
    if css_class is not None: html += '<ol class="' + css_class + '">'
    else: html += '<ol>'

    for item in items: html += '<li> ' + item
    html += '</ol>'

    return html

# -----------------------------------------------------------------

def unordered_list(items, css_class=None):

    """
    This function ...
    :param items:
    :param css_class:
    :return:
    """

    html = ""
    if css_class is not None: html += '<ul class="' + css_class + '">'
    else: html += '<ul>'

    for item in items: html += '<li> ' + item
    html += '</ul>'

    return html

# -----------------------------------------------------------------

def color(text, color):

    """
    This function ...
    :param text:
    :param color:
    :return:
    """

    return '<font color="' + color + '">' + text + '</font>'

# -----------------------------------------------------------------

def bold(text):

    """
    This function ...
    :param text:
    :return:
    """
    return '<b>' + text + '</b>'

# -----------------------------------------------------------------

def dictionary(dct, css_class=None, key_color=None, value_color=None, bold_keys=False, bold_values=False, identity_symbol=": "):

    """
    This function ...
    :param dct:
    :param css_class:
    :param key_color:
    :param value_color:
    :param bold_keys:
    :param bold_values:
    :param identity_symbol:
    :return:
    """

    items = []
    for name in dct:

        key = tostr(name)
        if key_color is not None: key = color(key, key_color)
        if bold_keys: key = bold(key)

        value = tostr(dct[name], round=True, html=True, decimal_places=1)
        if value_color is not None: value = color(value, value_color)
        if bold_values: value = bold(value)

        item = key + identity_symbol + value
        items.append(item)

    # Return
    return unordered_list(items, css_class=css_class)

# -----------------------------------------------------------------

def image(url, alttext=None, height=None, width=None, hover=None):

    """
    This function ...
    :param url:
    :param alttext:
    :param height:
    :param width:
    :param hover:
    :return:
    """

    code = '<img src="' + url + '"'
    if alttext is not None: code += ' alt="' + alttext + '"'
    if height is not None: code += " height=" + str(height) + "px"
    if width is not None: code += " width=" + str(width) + "px"
    if hover is not None: code += ' onmouseover="this.src=' + "'" + hover + "';" + '" ' + 'onmouseout="this.src=' + "'" + url + "'" + ';"'
    code += ">"
    return code

# -----------------------------------------------------------------

def get_image_url(image_code):

    """
    This function ...
    :param image_code:
    :return:
    """

    quote_character = image_code.split("src=")[1][0]
    return image_code.split("src=" + quote_character)[1].split(quote_character)[0]

# -----------------------------------------------------------------

def image_preview(url, text, title=None):

    """
    This function ...
    :param url:
    :param text:
    :param title:
    :return:
    """

    code = '<a href="' + url + '" class="preview"'
    if title is not None: code += 'title="' + title + '"'
    code += '>'
    code += text
    code += '</a>'
    return code

# -----------------------------------------------------------------

def hyperlink(url, text=None):

    """
    This function ...
    :param url:
    :param text:
    :return:
    """

    return '<a href="' + url + '">' + text + '</a>'

# -----------------------------------------------------------------

def mailto(address, text=None):

    """
    This function ...
    :param address:
    :param text:
    :return:
    """

    if text is None: text = address
    return hyperlink("mailto:" + address, text)

# -----------------------------------------------------------------

round_to_two_decimals_function = """
function roundToTwo(number)
{    
    return +(Math.round(number + "e+2") + "e-2");
}"""

# -----------------------------------------------------------------

round_to_decimals_function = """
function roundToDecimals(number, precision)
{
    return +(Math.round(number + "e+"+ precision) + "e-" + precision);
}"""

# -----------------------------------------------------------------

sleep_function = "function sleep(ms) { return new Promise(resolve => setTimeout(resolve, ms));}"

# -----------------------------------------------------------------

other_sleep_function = """
function sleep(milliseconds) {
  var start = new Date().getTime();
  for (var i = 0; i < 1e7; i++) {
    if ((new Date().getTime() - start) > milliseconds){
      break;
    }
  }
}"""

# -----------------------------------------------------------------

get_css_rule_function = """
function getCSSRule(ruleName)
{
    ruleName = ruleName.toLowerCase();
    var result = null;
    var find = Array.prototype.find;

    find.call(document.styleSheets, styleSheet => {
        result = find.call(styleSheet.cssRules, cssRule => {
            return cssRule instanceof CSSStyleRule 
                && cssRule.selectorText.toLowerCase() == ruleName;
        });
        return result != null;
    });
    return result;
}
"""

# -----------------------------------------------------------------

file_exists_function = """
function fileExists(image_url)
{
    var http = new XMLHttpRequest();
    window.alert(image_url);
    http.open('HEAD', image_url, false);
    http.send();
    return http.status != 404;
}
"""

# -----------------------------------------------------------------

image_exists_function = """
// The "callback" argument is called with either true or false
// depending on whether the image at "url" exists or not.
function imageExists_callback(url, allImages, i, dark_path, callback) 
{
    var img = new Image();
    img.onload = function() { callback(true, allImages, i, dark_path); };
    img.onerror = function() { callback(false, allImages, i, dark_path); };
    img.src = url;
}"""

# -----------------------------------------------------------------

def make_change_theme_function_template(button_id):

    return """
    function change()
    {
        var elem = document.getElementById('""" + button_id +  """');
        if (elem.value=="Light theme")
        {
            lightTheme();
            elem.value = "Dark theme";
        }   
        else
        { 
            darkTheme();
            elem.value = "Light theme";
        }
    }
    """

# -----------------------------------------------------------------

dark_theme_function = """
function darkTheme()
{
    darkBackground();
    darkText();
    darkImages();
    darkTables();
    darkExtra();
}
"""

# -----------------------------------------------------------------

light_theme_function = """
function lightTheme()
{
    lightBackground();
    lightText();
    lightImages();
    lightTables();
    lightExtra();
}
"""

# -----------------------------------------------------------------

dark_background_function = """
function darkBackground()
{
    document.body.style.backgroundColor = "black";
}
"""

# -----------------------------------------------------------------

dark_text_function = """
function darkText()
{
    document.body.style.color = "white";
}
"""

# -----------------------------------------------------------------

light_background_function = """
function lightBackground()
{
    document.body.style.backgroundColor = "white";
}
"""

# -----------------------------------------------------------------

light_text_function = """
function lightText()
{
    document.body.style.color = "black";
}
"""

# -----------------------------------------------------------------

split_by_last_dot_function = """
var splitByLastDot = function(text)
{
    var index = text.lastIndexOf('.');
    return [text.slice(0, index), text.slice(index + 1)]
}
"""

# -----------------------------------------------------------------

# //allImages = document.getElementsByTagName('img'); // faster as global??

dark_images_function = """
function darkImages()
{
    var allImages = document.getElementsByTagName('img');
    for(var i = 0; i < allImages.length ; i++)
    {     
        var result = splitByLastDot(allImages[i].src);
        var filename = result[0];
        var extension = result[1];
        
        //window.alert("filename: " + filename);
        
        var dark_path = filename + "_dark." + extension;
        
        imageExists_callback(dark_path, allImages, i, dark_path, function(existing, allImages, i, dark_path){
            if(existing == true) { allImages[i].src = dark_path; }
            });
                                
        // didn't work, something I tried myself
        //allImages[i].onerror = function(){
        //    allImages[i].src = filename + "." + extension;
        //    }
    }
}
"""

# -----------------------------------------------------------------

light_images_function = """
function lightImages()
{
    var allImages = document.getElementsByTagName('img');
    for(var i = 0; i < allImages.length ; i++)
    {            
        var result = splitByLastDot(allImages[i].src);
        var filename = result[0];
        var extension = result[1];
        
        var is_dark = filename.endsWith("_dark");
        
        var default_path = filename.split("_dark")[0] + "." + extension;   
        allImages[i].src = default_path;
    }
}
"""

# -----------------------------------------------------------------

dark_tables_function = """
function darkTables()
{
    var hovertable = getCSSRule('table.hovertable tr:hover');
    //window.alert(hovertable);
    hovertable.style["background-color"] = "#3b3d3f";
}
"""

# -----------------------------------------------------------------

light_tables_function = """
function lightTables()
{
    var hovertable = getCSSRule('table.hovertable tr:hover');
    //window.alert(hovertable);
    hovertable.style["background-color"] = "lightgrey";
}
"""

# -----------------------------------------------------------------

def make_theme_button(classes=None):

    """
    This function ...
    :param classes
    :return:
    """

    button_id = "theme"

    #code = """
    #<input id="b1" onclick="change()" type="button" value="Dark theme">
    code = button(button_id, "Dark theme", "change()")
    code += "\n"
    code += "<script>"

    code += "\n"
    code += get_css_rule_function
    code += "\n"
    code += file_exists_function
    code += "\n"
    code += image_exists_function
    code += "\n"
    code += split_by_last_dot_function
    code += "\n"
    code += make_change_theme_function_template(button_id)
    code += "\n"
    code += dark_theme_function
    code += "\n"
    code += light_theme_function
    code += "\n"
    code += dark_background_function
    code += "\n"
    code += dark_text_function
    code += "\n"
    code += light_background_function
    code += "\n"
    code += light_text_function
    code += "\n"
    code += dark_images_function
    code += "\n"
    code += light_images_function
    code += "\n"
    code += dark_tables_function
    code += "\n"
    code += light_tables_function

    # EXTRA
    if classes is not None:

        code += "\n"
        dark_function_code = "function darkExtra()\n{\n"
        light_function_code = "function lightExtra()\n{\n"

        # Loop over the different classes
        for class_name in classes:

            class_code = ""

            class_code += "    var allDivs = document.getElementsByTagName('div');"
            class_code += "\n"
            class_code += "    for(var i = 0; i < allDivs.length ; i++)\n    {"
            class_code += "\n"
            class_code += "        var thing = allDivs[i];"
            class_code += "\n"
            class_code += "        var currentClass = thing.className;"
            class_code += "\n"

            #class_code += "        window.alert(currentClass);\n"

            # Check class
            class_code += "        if (currentClass == '" + class_name + "')\n"
            class_code += "        {\n"

            #class_code += "            window.alert(currentClass);\n"

            dark_code = class_code
            light_code = class_code

            property_name = classes[class_name]

            ## DIDN'T WORK?
            # sequence = property_name.split(".")
            #
            # property_sequence_string = "[" + "][".join("'" + s + "'"  for s in sequence) + "]"
            #
            # dark_code += "            thing" + property_sequence_string + " = 'black';\n"
            # light_code += "            thing" + property_sequence_string + " = 'white';\n"

            ## NEW: only for one level!
            dark_code += '            thing.setAttribute("' + property_name + '", "black");\n'
            light_code += '            thing.setAttribute("' + property_name + '", "white");\n'

            dark_code += "        }\n        else {}"
            light_code += "        }\n        else {}"

            # end of for loop
            dark_code += "\n    }"
            light_code += "\n    }"

            dark_function_code += dark_code
            light_function_code += light_code

        # End function
        dark_function_code += "\n"
        dark_function_code += "}"

        # End function
        light_function_code += "\n"
        light_function_code += "}"

    else:

        dark_function_code = "function darkExtra()\n{}"
        light_function_code = "function lightExtra()\n{}"

    code += dark_function_code
    code += "\n\n"
    code += light_function_code

    code += "\n</script>"

    # Return the code
    return code

# -----------------------------------------------------------------

def make_script_function(name, script):

    """
    This function ...
    :param name:
    :param script:
    :return:
    """

    if " " in name: raise ValueError("Name cannot contain spaces")
    string = "function " + name + "()"
    string += "\n"
    string += "{"
    string += "\n"
    for line in script.split("\n"):
        string += "    " + line + "\n"
    string += "}"
    return string

# -----------------------------------------------------------------

def make_script_button(id, text, script, function_name):

    """
    This function ...
    :param id:
    :param text:
    :param script:
    :param function_name:
    """

    # code = button(id, text, function_name + "()")
    # code += "\n"
    # code += "<script>"
    # code += "\n"
    # code += make_script_function(function_name, script)
    # code += "\n</script>"

    code = "<script>"
    code += "\n"
    code += make_script_function(function_name, script)
    code += "\n</script>"
    code += "\n"
    code += button(id, text, function_name + "()")

    # Return the code
    return code

# -----------------------------------------------------------------

def make_slider(id, min_value, max_value, stepsize, default_value=None, action_function=None, text=None, text_id=None):

    """
    This function ...
    :param id:
    :param min_value:
    :param max_value:
    :param stepsize:
    :param default_value:
    :param action_function:
    :param text:
    :param text_id:
    :return:
    """

    """
    <input type="range" min="0" max="50" value="0" step="5" onchange="showValue(this.value)" />
    <span id="range">0</span>
    <script type="text/javascript">
    function showValue(newValue)
    {
        document.getElementById("range").innerHTML=newValue;
    }
    </script>"""

    code = "<input id='" + id + "' type='range' min='" + str(min_value) + "' max='" + str(max_value) + "'"

    if default_value is not None: code += " value='" + str(default_value) + "'"

    code += " step='" + str(stepsize) + "'"

    if action_function is not None: function_name = "processValue"
    else: function_name = "showValue"

    code += " onchange='" + function_name + "(this.value)'"

    code += "/>\n"

    if text_id is None: text_id = "range"

    if text is None:
        if default_value is None: default_text = "0"
        else: default_text = str(default_value)
    else:
        if default_value is None: default_text = text + "0"
        else: default_text = text + str(default_value)

    if text is None: text = ""

    code += "<span id='" + text_id + "'>\n"
    code += default_text + "\n"
    code += "</span>\n"
    
    #code += "\n"

    code += '<script type="text/javascript">\n'

    if action_function is not None:

        code += "function processValue(newValue)\n"
        code += "{\n"
        code += "    showValue(newValue);\n"
        code += "    " + action_function + "(newValue);\n"
        code += "}\n"
        code += "\n"

    code += "function showValue(newValue)\n"
    code += "{\n"
    code += "    document.getElementById('" + text_id + "').innerHTML = '" + strings.make_double_quoted(text) + "' + String(newValue);\n"
    code += "}\n"
    code += "</script>\n"

    # Return the code
    return code

# -----------------------------------------------------------------

custom_slider_css = """
@mixin rangeThumb {
  width: 18px;
  height: 18px;
  margin: -8px 0  0;
  border-radius: 50%;
  background: #37adbf;
  cursor: pointer;
  border: 0 !important;
}

@mixin rangeTrack {
  width: 100%;
  height: 2px;
  cursor: pointer;
  background: #b2b2b2;
}

.range {
  position: relative;
  width: 550px;
  height: 5px;
}

.range input {
  width: 100%;
  position: absolute;
  top: 2px;
  height: 0;
  -webkit-appearance: none;

  // Thumb
  &::-webkit-slider-thumb {
    -webkit-appearance: none; // needed again for Chrome & Safari
    @include rangeThumb;
  }

  &::-moz-range-thumb {
    @include rangeThumb;
  }

  &::-ms-thumb {
    @include rangeThumb;
  }

  // Track
  &::-webkit-slider-runnable-track {
    @include rangeTrack;
  }

  &::-moz-range-track {
    @include rangeTrack;
  }

  &::-ms-track {
    @include rangeTrack;
  }

  &:focus { // override outline/background on focus
    background: none;
    outline: none;
  }

  &::-ms-track { // A little somethin' somethin' for IE
    width: 100%;
    cursor: pointer;
    background: transparent;
    border-color: transparent;
    color: transparent;
  }
}

// Labels below slider
.range-labels {
  margin: 18px -41px 0;
  padding: 0;
  list-style: none;
  
  li {
    position: relative;
    float: left;
    width: 90.25px;
    text-align: center;
    color: #b2b2b2;
    font-size: 14px;
    cursor: pointer;
    
    &::before {
      position: absolute;
      top: -25px;
      right: 0;
      left: 0;
      content: "";
      margin: 0 auto;
      width: 9px;
      height: 9px;
      background: #b2b2b2;
      border-radius: 50%;
    }
  }
  
  .active {
    color: #37adbf;
  }
  
  .selected::before {
    background: #37adbf;
  }
  
  .active.selected::before {
    display: none;
  }
}"""

# -----------------------------------------------------------------

custom_slider_javascript = """
var sheet = document.createElement('style'),  
  $rangeInput = $('.range input'),
  prefs = ['webkit-slider-runnable-track', 'moz-range-track', 'ms-track'];

document.body.appendChild(sheet);

var getTrackStyle = function (el) {  
  var curVal = el.value,
      val = (curVal - 1) * 16.666666667,
      style = '';
  
  // Set active label
  $('.range-labels li').removeClass('active selected');
  
  var curLabel = $('.range-labels').find('li:nth-child(' + curVal + ')');
  
  curLabel.addClass('active selected');
  curLabel.prevAll().addClass('selected');
  
  // Change background gradient
  for (var i = 0; i < prefs.length; i++) {
    style += '.range {background: linear-gradient(to right, #37adbf 0%, #37adbf ' + val + '%, #fff ' + val + '%, #fff 100%)}';
    style += '.range input::-' + prefs[i] + '{background: linear-gradient(to right, #37adbf 0%, #37adbf ' + val + '%, #b2b2b2 ' + val + '%, #b2b2b2 100%)}';
  }

  return style;
}

$rangeInput.on('input', function () {
  sheet.textContent = getTrackStyle(this);
});

// Change input value on label click
$('.range-labels li').on('click', function () {
  var index = $(this).index();
  
  $rangeInput.val(index + 1).trigger('input');
  
});"""

# -----------------------------------------------------------------

def make_custom_slider(id, options, default=None, action_function=None):

    """
    This function ...
    :param id:
    :param options:
    :param default:
    :param action_function:
    :return:
    """

    """
    <div class="range">
      <input type="range" min="1" max="7" steps="1" value="1">
    </div>
    
    <ul class="range-labels">
      <li class="active selected">Today</li>
      <li>2 days</li>
      <li>3 days</li>
      <li>4 days</li>
      <li>5 days</li>
      <li>6 days</li>
      <li>7 days</li>
    </ul>"""

    noptions = len(options)

    if default is not None: default_index = options.index(default) + 1
    else: default_index = 1

    code = "<div class='range'>"

    code += "  <input id='" + id + "' type='range' min='1' max='" + str(noptions) + "' steps='1' value='" + str(default_index) + "'"
    if action_function is not None: code += " onchange='" + action_function + "(this.value)'"
    code += ">\n"

    code += "</div>\n"

    code += "<ul class='range-labels'>"

    for index, option in enumerate(options):

        if index + 1 == default_index: code += " <li class='active selected'>" + str(option) + "</li>\n"
        else: code += " <li>" + str(option) + "</li>\n"

    code += "</ul>\n"

    # Return the code
    return code

# -----------------------------------------------------------------

class SimpleTableCell(object):

    """A table class to create table cells.
    Example:
    cell = SimpleTableCell('Hello, world!')
    """

    def __init__(self, text, header=False, bgcolor=None, tostr_kwargs=None):

        """
        Table cell constructor.
        Keyword arguments:
        text -- text to be displayed
        header -- flag to indicate this cell is a header cell.
        """

        self.text = text
        self.header = header
        self.bgcolor = bgcolor
        self.tostr_kwargs = tostr_kwargs if tostr_kwargs is not None else {}

    # -----------------------------------------------------------------

    def __str__(self):

        """Return the HTML code for the table cell."""

        #if self.header: return '<th>%s</th>' % (self.text)
        #else: return '<td>%s</td>' % (self.text)

        if self.header:
            if self.bgcolor is not None: return "<th bgcolor='" + self.bgcolor + "'>" + tostr(self.text) + "</th>"
            else: return "<th>" + tostr(self.text, **self.tostr_kwargs) + "</th>"
        else:
            if self.bgcolor is not None: return "<td bgcolor='" + self.bgcolor + "'>" + tostr(self.text) + "</td>"
            else: return "<td>" + tostr(self.text, **self.tostr_kwargs) + "</td>"

# -----------------------------------------------------------------

class SimpleTableRow(object):

    """A table class to create table rows, populated by table cells.

    Example:
    # Row from list
    row = SimpleTableRow(['Hello,', 'world!'])

    # Row from SimpleTableCell
    cell1 = SimpleTableCell('Hello,')
    cell2 = SimpleTableCell('world!')
    row = SimpleTableRow([cell1, cell2])
    """

    def __init__(self, cells, header=False, bgcolors=None, tostr_kwargs=None):

        """Table row constructor.

        Keyword arguments:
        cells -- iterable of SimpleTableCell (default None)
        header -- flag to indicate this row is a header row.
                  if the cells are SimpleTableCell, it is the programmer's
                  responsibility to verify whether it was created with the
                  header flag set to True.
        """

        # Set cells
        if isinstance(cells[0], SimpleTableCell): self.cells = cells
        else:
            if bgcolors is not None: self.cells = [SimpleTableCell(cell, header=header, bgcolor=bgcolor, tostr_kwargs=tostr_kwargs) for cell, bgcolor in zip(cells, bgcolors)]
            else: self.cells = [SimpleTableCell(cell, header=header, tostr_kwargs=tostr_kwargs) for cell in cells]

        # Set header
        self.header = header

    # -----------------------------------------------------------------

    def __str__(self):

        """Return the HTML code for the table row and its cells as a string."""

        row = []

        row.append('<tr>')
        for cell in self.cells: row.append(str(cell))
        row.append('</tr>')
        
        return '\n'.join(row)

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        Iterate through row cells
        """

        for cell in self.cells:
            yield cell

    # -----------------------------------------------------------------

    def add_cell(self, cell):

        """
        Add a SimpleTableCell object to the list of cells.
        """

        self.cells.append(cell)

    # -----------------------------------------------------------------

    def add_cells(self, cells):

        """
        Add a list of SimpleTableCell objects to the list of cells.
        """

        for cell in cells: self.cells.append(cell)

# -----------------------------------------------------------------

class SimpleTable(object):

    """
    A table class to create HTML tables, populated by HTML table rows.

    Example:
    # Table from lists
    table = SimpleTable([['Hello,', 'world!'], ['How', 'are', 'you?']])

    # Table with header row
    table = SimpleTable([['Hello,', 'world!'], ['How', 'are', 'you?']],
                      header_row=['Header1', 'Header2', 'Header3'])

    # Table from SimpleTableRow
    rows = SimpleTableRow(['Hello,', 'world!'])
    table = SimpleTable(rows)
    """

    def __init__(self, rows, header=None, css_class=None, bgcolors=None, tostr_kwargs=None, subheader=None):

        """
        Table constructor.
        Keyword arguments:
        rows -- iterable of SimpleTableRow
        header:
        css_class -- table CSS class
        subheader_row:
        """

        if isinstance(rows[0], SimpleTableRow): self.rows = rows
        else:
            if bgcolors is not None: self.rows = [SimpleTableRow(row, bgcolors=bcolors, tostr_kwargs=tostr_kwargs) for row, bcolors in zip(rows, bgcolors)]
            else: self.rows = [SimpleTableRow(row, tostr_kwargs=tostr_kwargs) for row in rows]

        # Check
        if subheader is not None and header is None: raise ValueError("Cannot specify subheader but not header")

        # Set header row
        if header is None: self.header_row = None
        else:
            if subheader is not None:
                titles = []
                for top, sub in zip(header, subheader):
                    if sub is None: title = top + newline + " "
                    else: title = top + newline + sub
                    titles.append(title)
                self.header_row = SimpleTableRow(titles, header=True, tostr_kwargs=tostr_kwargs)
            else: self.header_row = SimpleTableRow(header, header=True, tostr_kwargs=tostr_kwargs)

        # Set CSS class
        self.css_class = css_class

    # -----------------------------------------------------------------

    @classmethod
    def rasterize(cls, cells, ncolumns=2, header=None, css_class=None, tostr_kwargs=None, subheader=None):

        """
        This function ...
        :param cells:
        :param ncolumns:
        :param header_row:
        :param css_class:
        :param tostr_kwargs:
        :param header:
        :return:
        """

        if len(cells) == 0: raise ValueError("No cells given")

        rows = []

        cellsiter = iter(list(cells) + [""] * ncolumns)

        nrows = numbers.round_up_to_int(len(cells) / float(ncolumns))

        #print("nrows", nrows)

        for i in range(nrows):

            row = []
            for j in range(ncolumns): row.append(cellsiter.next())
            rows.append(row)

        # Create and return
        return cls(rows, header=header, css_class=css_class, tostr_kwargs=tostr_kwargs, subheader=subheader)

    # -----------------------------------------------------------------

    @classmethod
    def from_cells(cls, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        return cls.rasterize(args, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_table(cls, table, css_class=None, bgcolors=None, tostr_kwargs=None, column_names=None, extra_column=None,
                   extra_column_label=None):

        """
        This function ...
        :param table:
        :param css_class:
        :param bgcolors:
        :param tostr_kwargs:
        :param column_names:
        :param extra_column:
        :param extra_column_label:
        :return:
        """

        # Set column names
        if column_names is None: column_names = table.column_names

        # Check
        if extra_column_label is not None and extra_column is None: raise ValueError("Extra column label is specified but extra column is not given")

        # Get rows
        #rows = table.as_tuples(add_units=False)
        rows = table.as_lists(add_units=False)
        if extra_column is not None:
            for i in range(len(rows)): rows[i].append(extra_column[i])

        # Set header
        header = column_names
        if extra_column is not None:
            label = extra_column_label if extra_column_label is not None else ""
            header.append(label)

        # Set subheader
        subheader = table.unit_strings
        if extra_column is not None: subheader.append("")

        # Create the HTML table
        return cls(rows, header=header, subheader=subheader, css_class=css_class, bgcolors=bgcolors, tostr_kwargs=tostr_kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_composites(cls, composites, css_class=None, bgcolors=None, tostr_kwargs=None, labels=None, label="-",
                        extra_column=None, extra_column_label=None):

        """
        This function ...
        :param composites:
        :param css_class:
        :param bgcolors:
        :param tostr_kwargs:
        :param labels:
        :param label:
        :param extra_column:
        :param extra_column_label:
        :return:
        """

        # Create the table
        table = SmartTable.from_composites(*composites, labels=labels, label=label)

        # Create the HTML table
        return cls.from_table(table, css_class=css_class, bgcolors=bgcolors, tostr_kwargs=tostr_kwargs,
                              column_names=table.descriptions, extra_column=extra_column, extra_column_label=extra_column_label)

    # -----------------------------------------------------------------

    @classmethod
    def from_composite(cls, composite, css_class=None, bgcolors=None, tostr_kwargs=None, key_label="Property", value_label="Value"):

        """
        This function ...
        :param composite:
        :param css_class:
        :param bgcolors:
        :param tostr_kwargs:
        :param key_label:
        :param value_label:
        :return:
        """

        return cls(composite.as_tuples(), header=[key_label, value_label], css_class=css_class)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        Return the HTML code for the table as a string.
        """

        table = []

        if self.css_class: table.append('<table class=%s>' % self.css_class)
        else: table.append('<table>')

        # Add header and subheader
        if self.header_row is not None: table.append(str(self.header_row))

        for row in self.rows: table.append(str(row))

        table.append('</table>')
        
        return '\n'.join(table)

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        Iterate through table rows
        """

        for row in self.rows:
            yield row

    # -----------------------------------------------------------------

    def add_row(self, row):

        """
        Add a SimpleTableRow object to the list of rows.
        """

        self.rows.append(row)

    # -----------------------------------------------------------------

    def add_rows(self, rows):

        """
        Add a list of SimpleTableRow objects to the list of rows.
        """

        for row in rows: self.rows.append(row)

# -----------------------------------------------------------------

class HTMLPage(object):

    """
    A class to create HTML pages containing CSS and tables.
    """

    def __init__(self, title, css=None, css_path=None, style=None, encoding="utf-8", language="en", head=None, footing=None, body_settings=None, javascript_path=None):

        """
        HTML page constructor.
        Keyword arguments:
        tables -- List of SimpleTable objects
        css -- Cascading Style Sheet specification that is appended before the
               table string
        encoding -- Characters encoding. Default: UTF-8
        :param language:
        :param head:
        :param footing:
        :param body_settings:
        :param javascript_path:
        """

        # Set properties
        self.title = title
        self.css = css
        self.css_path = css_path
        self.style = style
        self.encoding = encoding
        self.language = language
        self.head = head
        self.footing = footing
        self.body_settings = body_settings
        self.javascript_path = javascript_path

        # The contents
        self.contents = []

        # The path of the HTML file
        self.path = None

    # -----------------------------------------------------------------

    def __iadd__(self, element):

        """
        This function ...
        :param element:
        :return:
        """

        self.contents.append(element)
        return self

    # -----------------------------------------------------------------

    def __str__(self):

        """
        Returns the HTML page as a string.
        """

        lines = []

        lines.append("<!DOCTYPE html>")
        lines.append('<html lang="' + self.language +'">')
        lines.append("<head>")

        # Add css path(s)
        if self.css_path is not None:
            if types.is_string_sequence(self.css_path):
                for url in self.css_path: lines.append('    <link rel="stylesheet" type="text/css" href="{url}">'.format(url=url))
            elif types.is_string_type(self.css_path): lines.append('    <link rel="stylesheet" type="text/css" href="{url}">'.format(url=self.css_path))
            else: raise ValueError("Invalid type for css_path")

        # Add custom css
        if self.css is not None: lines.append(make_css(self.css))

        # Add javascript path(s)
        if self.javascript_path is not None:
            if types.is_string_sequence(self.javascript_path):
                for url in self.javascript_path: lines.append('    <script type="text/javascript" src="{url}"></script>'.format(url=url))
            elif types.is_string_type(self.javascript_path): lines.append('    <script type="text/javascript" src="{url}"></script>'.format(url=self.javascript_path))
            else: raise ValueError("Invalid type for javascript_path")

        # Set encoding
        #lines.append('    <meta charset="UTF-8">')
        lines.append('    <meta http-equiv="Content-Type" content="text/html;charset=%s">' % self.encoding)

        # Set title
        lines.append('    <title>{title}</title>'.format(title=self.title))

        # Set other head information
        if self.head is not None: lines.append(self.head)

        #
        lines.append("</head>")

        if self.body_settings is not None: body_start = "<body " + stringify_dict(self.body_settings, identity_symbol="=", quote_key=False)[1] + ">"
        else: body_start = "<body>"
        lines.append(body_start)

        # Start of styled div
        if self.style is not None: lines.append('<div class="' + self.style + '">')

        # Add content
        for element in self.contents: lines.append(str(element))

        # Footing
        if self.footing is not None: lines.append(str(self.footing))

        # End of styled div
        if self.style is not None: lines.append('</div>')

        # End of body and page
        lines.append("</body>")
        lines.append("</html>")

        # Create the page and return
        return '\n'.join(lines)

    # -----------------------------------------------------------------

    def saveto(self, filepath, update_path=True):

        """
        Save HTML page to a file using the proper encoding:
        :param update_path:
        """

        # Write each line
        #with codecs.open(filepath, 'w', self.encoding) as outfile:
        #    for line in str(self): outfile.write(line)

        fs.write_text(filepath, str(self))

        # Set filepath
        self.path = filepath

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        if self.path is None: raise ValueError("Path not defined")
        self.saveto(self.path)

# -----------------------------------------------------------------

def fit_data_to_columns(data, num_cols):

    """
    Format data into the configured number of columns in a proper format to
    generate a SimpleTable.

    Example:
    test_data = [str(x) for x in range(20)]
    fitted_data = fit_data_to_columns(test_data, 5)
    table = SimpleTable(fitted_data)
    """

    num_iterations = len(data)/num_cols

    if len(data)%num_cols != 0:
        num_iterations += 1

    return [data[num_cols*i:num_cols*i + num_cols] for i in range(num_iterations)]

# -----------------------------------------------------------------
