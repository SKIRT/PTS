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
from . import sequences

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
italic_template = "<span style='font-style:italic'>{text}</span>"
fontsize_template = "<span style='font-size:{size}px'>{text}</span>"
small_template = "<small>{text}</small>"
center_template = "<div style='text-align:center;'>{text}</div>"
strikethrough_template = "<span style='text-decoration:line-through;'>{text}</span>"

# -----------------------------------------------------------------

def font_size(text, size):

    """
    This function ...
    :param text:
    :param size:
    :return:
    """

    return '<font size = "' + str(size) + '">' + text + '</font>'

# -----------------------------------------------------------------

def strikethrough(text):

    """
    This function ...
    :param text:
    :return:
    """

    return strikethrough_template.format(text=text)

# -----------------------------------------------------------------

def italic(text):
    
    """
    This function ...
    :param text: 
    :return: 
    """

    return italic_template.format(text=text)

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

def big(text):

    """
    This function ...
    :param text:
    :return:
    """

    return '<big>' + text + '</big>'

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

def image(url, alttext=None, height=None, width=None, hover=None, class_name=None):

    """
    This function ...
    :param url:
    :param alttext:
    :param height:
    :param width:
    :param hover:
    :param class_name:
    :return:
    """

    code = '<img src="' + url + '"'
    if class_name is not None: code += ' class="' + class_name + '"'
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
    if title is not None: code += ' title="' + title + '"'
    code += '>'
    code += text
    code += '</a>'
    return code

# -----------------------------------------------------------------

def anchor(name, text=None):

    """
    This function ...
    :param name:
    :param text:
    :return:
    """

    code = '<a name="' + name + '"'
    code += ">"

    if text is not None: code += text

    code += "</a>"
    return code

# -----------------------------------------------------------------

def anchor_link(name, text):

    """
    This function ...
    :param name:
    :param text:
    :return:
    """

    return hyperlink("#" + name, text)

# -----------------------------------------------------------------

def hyperlink(url, text):

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
        if (styleSheet.cssRules == null) { return null; }
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
}"""

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

def make_dark_theme_function(background=True, text=True, images=True, tables=True, extra=True):

    """
    This function ...
    :param background:
    :param text:
    :param images:
    :param tables:
    :param extra:
    :return:
    """

    code = "function darkTheme()\n"
    code += "{\n"

    if background: code += "    darkBackground();\n"
    if text: code += "    darkText();\n"
    if images: code += "    darkImages();\n"
    if tables: code += "    darkTables();\n"
    if extra: code += "    darkExtra();\n"

    code += "}"
    return code

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

def make_light_theme_function(background=True, text=True, images=True, tables=True, extra=True):

    """
    This function ...
    :return:
    """

    code = "function lightTheme()\n"
    code += "{\n"

    if background: code += "    lightBackground();\n"
    if text: code += "    lightText();\n"
    if images: code += "    lightImages();\n"
    if tables: code += "    lightTables();\n"
    if extra: code += "    lightExtra();\n"

    code += "}"
    return code

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
        
        //window.alert(String(allImages[i].classList));
        //window.alert(String(allImages[i].classList.contains('invertable')));
        
        imageExists_callback(dark_path, allImages, i, dark_path, function(existing, allImages, i, dark_path)
        {    
            if(existing == true) { allImages[i].src = dark_path;}
        });
        
        if (!allImages[i].src.includes("dark") && allImages[i].classList.contains('invertable'))
        {    
            allImages[i].classList.remove('invertable');
            allImages[i].classList.add('inverted');
             
            //window.alert(String(allImages[i].classList));
        }
                                
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
        
        //window.alert(String(allImages[i].classList))
        //window.alert(String(allImages[i].classList.contains('inverted')));
        
        if (allImages[i].classList.contains('inverted'))
        {
            allImages[i].classList.remove('inverted');
            allImages[i].classList.add('invertable');
        }
    }
}
"""

# -----------------------------------------------------------------

dark_tables_function = """
function darkTables()
{
    var hovertable = getCSSRule('table.hovertable tr:hover, table.sortable tr:hover');
    if (hovertable == null) { return; }
    //window.alert(hovertable);
    hovertable.style["background-color"] = "#3b3d3f";
}
"""

# -----------------------------------------------------------------

light_tables_function = """
function lightTables()
{
    var hovertable = getCSSRule('table.hovertable tr:hover, table.sortable tr:hover');
    if (hovertable == null) { return; }
    //window.alert(hovertable);
    hovertable.style["background-color"] = "lightgrey";
}
"""

# -----------------------------------------------------------------

def make_theme_button(classes=None, images=True):

    """
    This function ...
    :param classes:
    :param images:
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
    #code += dark_theme_function
    code += make_dark_theme_function(images=images)
    code += "\n"
    #code += light_theme_function
    code += make_light_theme_function(images=images)
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

def make_script(script):

    """
    This function ...
    :param script:
    :return:
    """

    code = "<script>"
    code += "\n"
    code += script
    code += "\n</script>"
    return code

# -----------------------------------------------------------------

def make_script_function_script(function_name, script):

    """
    This function ...
    :param function_name:
    :param script:
    :return:
    """

    code = "<script>"
    code += "\n"
    code += make_script_function(function_name, script)
    code += "\n</script>"

    return code

# -----------------------------------------------------------------

def make_script_button(id, text, script, function_name, quote_character='"'):

    """
    This function ...
    :param id:
    :param text:
    :param script:
    :param function_name:
    :param quote_character:
    """

    code = make_script_function_script(function_name, script)
    code += "\n"
    code += button(id, text, function_name + "();", quote_character=quote_character)

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

def make_basic_custom_slider(id, options, default=None, action_function=None, pass_value=True):

    """
    This function ...
    :param id:
    :param options:
    :param default:
    :param action_function:
    :param pass_value:
    :return:
    """

    noptions = len(options)

    if default is not None: default_index = options.index(default) + 1
    else: default_index = 1

    code = ""
    code += "  <input id='" + id + "' type='range' min='1' max='" + str(noptions) + "' steps='1' value='" + str(default_index) + "'"
    # if action_function is not None: code += " onchange='" + action_function + "(this.value)'"
    if action_function is not None:
        if pass_value: code += " onchange='" + action_function + "(this.value)' oninput='" + action_function + "(this.value)'"
        else: code += " onchange='" + action_function + "()' oninput='" + action_function + "()'"
    code += ">\n"

    # Return the code
    return code

# -----------------------------------------------------------------

def make_custom_slider(id, options, default=None, action_function=None, pass_value=True):

    """
    This function ...
    :param id:
    :param options:
    :param default:
    :param action_function:
    :param pass_value:
    :return:
    """

    # EXAMPLE:
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
    #if action_function is not None: code += " onchange='" + action_function + "(this.value)'"
    if action_function is not None:
        if pass_value: code += " onchange='" + action_function + "(this.value)' oninput='" + action_function+ "(this.value)'"
        else: code += " onchange='" + action_function + "()' oninput='" + action_function+ "()'"
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

def make_image_slider(image_id, urls, labels, default, width=None, height=None, basic=True, extra_urls=None, img_class=None, extra_img_class=None):

    """
    This function ...
    :param image_id:
    :param urls:
    :param labels:
    :param default:
    :param width:
    :param height:
    :param basic:
    :param extra_urls:
    :param img_class:
    :param extra_img_class:
    :return:
    """

    # Set slider ID
    slider_id = image_id + "Slider"

    # Set span ID
    span_id = image_id + "Span"

    # Set function name
    function_name = "update" + image_id

    # Get the default URL
    #default_index = labels.index(default)
    #default_url = urls[default_index]

    code = ""

    # Add slider
    if basic: code += center(make_basic_custom_slider(slider_id, labels, default=default, action_function=function_name))
    else: code += center(make_custom_slider(slider_id, labels, default=default, action_function=function_name))

    code += newline

    code += '\n'
    code += "<center>"
    code += '<span align="center" id="' + span_id + '">' + str(default) + "</span>"
    code += "</center>"
    code += "\n"

    code += newline

    # Code for image holder
    code += "<center>"
    code += '<img id="' + image_id + '"'
    if img_class is not None: code += ' class="' + img_class + '"'
    if width is not None: code += ' width="' + str(width)  + 'px"'
    if height is not None: code += ' height="' + str(height) + 'px"'
    code += '>'

    if extra_urls is not None:

        extra_image_id = image_id + "_extra"
        code += '<img id="' + extra_image_id + '"'
        if extra_img_class is not None: code += ' class="' + img_class + '"'
        if width is not None: code += ' width="' + str(width) + 'px"'
        if height is not None: code += ' height="' + str(height) + 'px"'
        code += '>'

    else: extra_image_id = None

    code += "</center>"

    # Make script
    script = ""

    labels_variable_name = "labels_" + image_id
    urls_variable_name = "urls_" + image_id

    string_labels = [str(label) for label in labels]
    script += 'var ' + labels_variable_name + ' = [' + tostr(string_labels, add_quotes=True) + '];\n'
    script += 'var ' + urls_variable_name + ' = [' + tostr(urls, add_quotes=True) + '];\n'

    if extra_urls is not None:
        extra_urls_variable_name = "urls_" + extra_image_id
        script += 'var ' + extra_urls_variable_name + ' = [' + tostr(extra_urls, add_quotes=True) + '];\n'
    else: extra_urls_variable_name = None

    script += "\n"

    # get_url_function_name = 'get' + image_id + 'URL'
    # script += 'function ' + get_url_function_name + '(label)\n'
    # script += '{\n'
    # script += '    var index = labels.indexOf("label");\n'
    # script += '    return urls[index];\n'
    # script += '}\n'

    script += "\n"

    script += 'var val = document.getElementById("' + slider_id + '").value;\n'
    script += 'var index = val - 1;\n'
    script += 'document.getElementById("' + span_id + '").innerHTML = ' + labels_variable_name + '[index];\n'
    #script += 'document.getElementById("' + image_id + '").src = ' + get_url_function_name + '(val);\n'
    script += 'document.getElementById("' + image_id + '").src = ' + urls_variable_name + '[index];\n'

    if extra_urls is not None: script += 'document.getElementById("' + extra_image_id + '").src = ' + extra_urls_variable_name + '[index];\n'

    script += "\n"
    script += 'function ' + function_name + '(newVal)\n'
    script += '{\n'
    script += '    var newIndex = newVal - 1;\n'
    script += '    document.getElementById("' + span_id + '").innerHTML = ' + labels_variable_name + '[newIndex];\n'
    script += '    document.getElementById("' + image_id + '").src = ' + urls_variable_name + '[newIndex];\n'

    if extra_urls is not None: script += '    document.getElementById("' + extra_image_id + '").src = ' + extra_urls_variable_name + '[newIndex];\n'

    script += '\n}'

    # Add javascript
    code += make_script(script)

    # Return the HTML code
    return code

# -----------------------------------------------------------------

def make_usable(name):

    """
    This function ...
    :return:
    """

    return name.replace(" ", "").replace(".", "")  # CANNOT CONTAIN SPACES AND DOTS!!

# -----------------------------------------------------------------

def make_if_statement(labels, values, indent=""):

    """
    This function ...
    :param labels:
    :param values:
    :param indent:
    :return:
    """

    code = ""

    code += indent + "if ("

    comparison_strings = []
    for label, value in zip(labels, values):
        string = label + " == " + str(value)
        comparison_strings.append(string)

    code += " && ".join(comparison_strings)

    code += ")\n"

    return code

# -----------------------------------------------------------------

def make_else_statement(labels, values, indent=""):

    """
    This function ...
    :param labels:
    :param values:
    :param indent:
    :return:
    """

    code = ""

    code += indent + '} else if ('

    comparison_strings = []
    for label, value in zip(labels, values):
        string = label + " == " + str(value)
        comparison_strings.append(string)

    code += " && ".join(comparison_strings)

    code += ")\n"

    return code

# -----------------------------------------------------------------

def make_multi_image_sliders(image_names, urls, regulator_names, labels, default, width=None,
                             height=None, basic=True, img_class=None, label_images=None, table_class=None):

    """
    This function ...
    :param image_names:
    :param urls:
    :param regulator_ids:
    :param labels:
    :param default:
    :param width:
    :param height:
    :param basic:
    :param img_class:
    :param label_images:
    :param table_class:
    :return:
    """

    # Set labels and default values
    if types.is_sequence(labels): labels = {regulator_name: labels for regulator_name in regulator_names}
    if not types.is_dictionary(default): default = {regulator_name: default for regulator_name in regulator_names}

    # Check label images (if specified)
    if label_images is not None:
        for regulator_name in regulator_names:
            if regulator_name not in label_images: raise ValueError("Regulator name '" + regulator_name + "' not in label_images dictionary")
            else:
                for label in labels[regulator_name]:
                    if label not in label_images[regulator_name]: raise ValueError("Label '" + str(label) + "' not in label_images subdictionary for regulator name '" + regulator_name + "'")

    regulator_ids = [make_usable(regulator_name) for regulator_name in regulator_names]

    # Transform labels and default dicts
    labels = {make_usable(regulator_name): label for regulator_name, label in labels.items()}
    default = {make_usable(regulator_name): value for regulator_name, value in default.items()}

    # Transform label images
    if label_images is not None:
        label_images = {make_usable(regulator_name): [image_paths_for_labels[label] for label in labels[make_usable(regulator_name)]] for regulator_name, image_paths_for_labels in label_images.items()}
        # IS NOW A DICT OF LISTS

    code = ""

    #function_names = dict()
    function_name = "update_images"

    span_ids = dict()
    slider_ids = dict()

    #left_cell = ""
    left_cells = ["" for _ in range(2)]

    label_image_ids = dict()

    # Loop over the regulators
    for index, regulator_id in enumerate(regulator_ids):

        function_name_for_regulator = "update_for_" + regulator_id

        # Set slider ID
        slider_id = regulator_id + "Slider"

        # Set span ID
        span_id = regulator_id + "Span"

        # Set function name
        #function_name = "update" + regulator_id

        # Get the default URL
        # default_index = labels.index(default)
        # default_url = urls[default_index]

        which_cell = index % 2

        left_cells[which_cell] += regulator_id.upper()
        left_cells[which_cell] += newline + newline

        # Add slider
        if basic: left_cells[which_cell] += center(make_basic_custom_slider(slider_id, labels[regulator_id], default=default[regulator_id], action_function=function_name_for_regulator, pass_value=False))
        else: left_cells[which_cell] += center(make_custom_slider(slider_id, labels[regulator_id], default=default[regulator_id], action_function=function_name_for_regulator, pass_value=False))

        left_cells[which_cell] += newline

        left_cells[which_cell] += '\n'
        left_cells[which_cell] += "<center>"
        left_cells[which_cell] += '<span align="center" id="' + span_id + '">' + str(default[regulator_id]) + "</span>"
        left_cells[which_cell] += "</center>"
        left_cells[which_cell] += "\n"

        if label_images:

            label_image_id = "label_image_" + regulator_id
            label_image_ids[regulator_id] = label_image_id

            label_image_height = 100
            label_image_width = None

            left_cells[which_cell] += "<center>"
            left_cells[which_cell] += '<img id="' + label_image_id + '"'
            #if img_class is not None: left_cell += ' class="' + img_class + '"'
            if label_image_width is not None: left_cells[which_cell] += ' width="' + str(label_image_width) + 'px"'
            if label_image_height is not None: left_cells[which_cell] += ' height="' + str(label_image_height) + 'px"'
            left_cells[which_cell] += '>'
            left_cells[which_cell] += "</center>"
            left_cells[which_cell] += "\n"

        left_cells[which_cell] += newline

        span_ids[regulator_id] = span_id
        slider_ids[regulator_id] = slider_id

        # remember function names
        #function_names[regulator_id] = function_name

    # Create image IDs
    #image_ids = [make_usable(image_name) for image_name in image_names]
    image_ids = []
    categories_and_names_for_image_id = dict()
    for category in image_names:
        for name in image_names[category]:
            image_id = category + "___" + make_usable(name)
            image_ids.append(image_id)
            categories_and_names_for_image_id[image_id] = (category, name)

    right_cells = ["" for _ in range(4)]

    # Loop over the images
    for index, image_id in enumerate(image_ids):

        which_cell = index % 4
        right_cells[which_cell] += font_size(image_id.upper(), 1)
        right_cells[which_cell] += newline

        # Code for image holder
        #right_cell += "<center>"
        right_cells[which_cell] += '<img id="' + image_id + '"'
        if img_class is not None: right_cells[which_cell] += ' class="' + img_class + '"'
        if width is not None: right_cells[which_cell] += ' width="' + str(width) + 'px"'
        if height is not None: right_cells[which_cell] += ' height="' + str(height) + 'px"'
        right_cells[which_cell] += '>'
        #code += "</center>"

        #if numbers.is_even(index): right_cell += newline
        #if index % 4 == 0 and index != 0: right_cell += newline
        right_cells[which_cell] += newline

    # All cells
    all_cells = left_cells + right_cells

    # Make and add table
    #table = SimpleTable.one_row(left_cell, right_cell, css_class=table_class)
    table = SimpleTable.one_row(*all_cells, css_class=table_class)
    code += str(table)

    # Make script
    script = ""

    # # Set labels
    labels_variable_names = dict()
    for regulator_id in regulator_ids:

         labels_variable_name = "labels_" + regulator_id

         label_images_variable_name = "label_image_paths_" + regulator_id

         labels_variable_names[regulator_id] = labels_variable_name

         string_labels = [str(label) for label in labels[regulator_id]]
         script += 'var ' + labels_variable_name + ' = [' + tostr(string_labels, add_quotes=True, delimiter=', ') + '];\n'

         if label_images is not None:
            script += 'var ' + label_images_variable_name + ' = [' + tostr(label_images[regulator_id], add_quotes=True, delimiter=", ") + '];\n'

         script += "\n"

    # script += 'var val = document.getElementById("' + slider_id + '").value;\n'
    # script += 'var index = val - 1;\n'
    # script += 'document.getElementById("' + span_id + '").innerHTML = ' + labels_variable_name + '[index];\n'
    # # script += 'document.getElementById("' + image_id + '").src = ' + get_url_function_name + '(val);\n'
    # script += 'document.getElementById("' + image_id + '").src = ' + urls_variable_name + '[index];\n'
    # script += "\n"

    # MOVED TO BELOW DEFINITION OF FUNCTIONS???
    # # NEW
    # if label_images is not None:
    #     for regulator_id in regulator_ids:
    #         function_name_regulator = "update_" + regulator_id
    #         script += function_name_regulator + '();\n'
    #         script += "\n"
    #
    # # Initialize images
    # script += function_name + '();\n'
    # script += "\n"

    # Make functions for each image
    for regulator_id in regulator_ids:

        function_name_for_regulator = "update_for_" + regulator_id
        function_name_regulator = "update_" + regulator_id

        slider_id = slider_ids[regulator_id]

        label_images_variable_name = "label_image_paths_" + regulator_id

        if label_images is not None:

            label_image_id = label_image_ids[regulator_id]

            script += 'function ' + function_name_regulator + '()\n'
            script += '{\n'
            script += '    var val_' + regulator_id + ' = document.getElementById("' + slider_id + '").value;\n'
            script += '    var index_' + regulator_id + ' = val_' + regulator_id + ' - 1;\n'
            script += '    var path_' + regulator_id + ' = ' + label_images_variable_name + '[index_' + regulator_id + '];\n'
            script += '    document.getElementById("' + label_image_id + '").src = path_' + regulator_id + ';\n'
            script += '}\n\n'

        script += 'function ' + function_name_for_regulator + '()\n'
        script += '{\n'
        if label_images is not None: script += '    ' + function_name_regulator + "();\n"
        script += '    ' + function_name + "();\n"
        script += '}\n\n'

    # # Make function for each regulator
    # for regulator_id in regulator_ids:
    #
    #     # Get function name
    #     function_name = function_names[regulator_id]
    #
    #     script += 'function ' + function_name + '(newVal)\n'
    #     script += '{\n'
    #     #script += '    var newIndex = newVal - 1;\n'
    #     #script += '    document.getElementById("' + span_id + '").innerHTML = ' + labels_variable_name + '[newIndex];\n'
    #     #script += '    document.getElementById("' + image_id + '").src = ' + urls_variable_name + '[newIndex];\n'
    #     script += '\n}'

    # Make the image change functions
    # Loop over the images, change them
    dependencies = dict()
    for image_id in image_ids:

        image_change_function_name = "update_" + image_id

        category, name = categories_and_names_for_image_id[image_id]

        keys_list = []
        for levels_dict in urls[category][name]:
            keys = levels_dict.keys()
            #keys_list.append(list(sorted(keys)))
            keys_list.append(keys)

        regulator_names_image = sequences.get_all_equal_value(keys_list)
        regulator_ids_image = [make_usable(regulator_name) for regulator_name in regulator_names_image]
        dependencies[image_id] = regulator_ids_image

        regulator_ids_value_names = [strings.replace_first_digit_by_word(regulator_id) + "_value" for regulator_id in regulator_ids_image]

        script += 'function ' + image_change_function_name + '(' + tostr(regulator_ids_value_names) + ')\n'
        script += '{\n'

        # Make if statements
        first = True
        for levels_dict in urls[category][name]:

            levels_values = []
            for rni in regulator_names_image:
                sigma_level = levels_dict[rni]
                levels_values.append(sigma_level)

            if first: script += make_if_statement(regulator_ids_value_names, levels_values, indent="    ")
            else: script += make_else_statement(regulator_ids_value_names, levels_values, indent="    ")

            first = False

            # Path or URL
            path = urls[category][name][levels_dict]

            script += '    {\n'

            #script += '    document.getElementById("' + span_id + '").innerHTML = ' + labels_variable_name + '[newIndex];\n'
            script += '        document.getElementById("' + image_id + '").src = "' + path + '";\n'

        script += '    } else {\n'
        script += '        window.alert("ERROR: ' + ", ".join([name + ': " + String(' + value_name + ') + "' for name, value_name in zip(regulator_ids_image, regulator_ids_value_names)]) + '");\n'
        script += '    }\n'

        script += '}\n'
        script += '\n'

    # Make the function
    script += 'function ' + function_name + '()\n'
    script += '{\n'

    # Get values for regulators
    for regulator_id in regulator_ids:

        labels_variable_name = labels_variable_names[regulator_id]
        span_id = span_ids[regulator_id]
        slider_id = slider_ids[regulator_id]

        script += '    var val_' + regulator_id + ' = document.getElementById("' + slider_id + '").value;\n'
        script += '    var index_' + regulator_id + ' = val_' + regulator_id + ' - 1;\n'
        script += '    var label_' + regulator_id + ' = ' + labels_variable_name + '[index_' + regulator_id + '];\n'
        script += '    document.getElementById("' + span_id + '").innerHTML = label_' + regulator_id + ';\n'
        script += '\n'

    # Update all images
    for image_id in image_ids:

        image_change_function_name = "update_" + image_id

        script += '    arguments_' + image_id + ' = [];\n'
        for regulator_id in dependencies[image_id]:
            script += '    arguments_' + image_id + '.push(label_' + regulator_id + ');\n'

        # UPDATE THE IMAGE
        #script += '    ' + image_change_function_name + '(arguments_' + image_id + ');\n'
        script += '    ' + image_change_function_name + '.apply(this, arguments_' + image_id + ');\n'
        script += "\n"

    script += '\n}'

    # NEW
    if label_images is not None:
        for regulator_id in regulator_ids:
            function_name_regulator = "update_" + regulator_id
            script += function_name_regulator + '();\n'
            script += "\n"

    # Initialize images
    script += function_name + '();\n'
    script += "\n"

    # Add javascript
    code += make_script(script)

    # Return the HTML code
    return code

# -----------------------------------------------------------------

class SimpleTableCell(object):

    """
    A table class to create table cells.
    Example:
    cell = SimpleTableCell('Hello, world!')
    """

    def __init__(self, text, header=False, bgcolor=None, tostr_kwargs=None, strike=False, bold=False, italic=False, wrap=True):

        """
        Table cell constructor.
        Keyword arguments:
        text -- text to be displayed
        header -- flag to indicate this cell is a header cell.
        :param wrap:
        """

        self.text = text
        self.header = header
        self.bgcolor = bgcolor
        self.tostr_kwargs = tostr_kwargs if tostr_kwargs is not None else {}
        self.strike = strike
        self.bold = bold
        self.italic = italic
        self.wrap = wrap

    # -----------------------------------------------------------------

    def __str__(self):

        """
        Return the HTML code for the table cell.
        """

        #if self.header: return '<th>%s</th>' % (self.text)
        #else: return '<td>%s</td>' % (self.text)

        if self.header:

            text = tostr(self.text, **self.tostr_kwargs)
            if self.strike: text = strikethrough(text)
            if self.bold: text = bold(text)
            if self.italic: text = italic(text)
            if self.bgcolor is not None: return "<th bgcolor='" + self.bgcolor + "'>" + text + "</th>"
            else: return "<th>" + text + "</th>"

        else:
            text = tostr(self.text, **self.tostr_kwargs)
            if self.strike: text = strikethrough(text)
            if self.bold: text = bold(text)
            if self.italic: text = italic(text)
            if self.bgcolor is not None: return "<td bgcolor='" + self.bgcolor + "'>" + text + "</td>"
            else: return "<td>" + text + "</td>"

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

    def __init__(self, cells, header=False, bgcolors=None, tostr_kwargs=None, strike=False, bold=False, italic=False, wrap=True):

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
            #if bgcolors is not None: self.cells = [SimpleTableCell(cell, header=header, bgcolor=bgcolor, tostr_kwargs=tostr_kwargs) for cell, bgcolor in zip(cells, bgcolors)]
            #else: self.cells = [SimpleTableCell(cell, header=header, tostr_kwargs=tostr_kwargs) for cell in cells]

            self.cells = []
            for index in range(len(cells)):

                bgcolor = bgcolors[index] if bgcolors is not None else None
                # strikei = strike[index] if strike is not None else None
                # boldi = bold[index] if bold is not None else None
                # italici = italic[index] if italic is not None else None

                celli = SimpleTableCell(cells[index], header=header, tostr_kwargs=tostr_kwargs, bgcolor=bgcolor, strike=strike, bold=bold, italic=italic, wrap=wrap)
                self.cells.append(celli)

        # Set wrap property
        self.wrap = wrap

        # Set header
        self.header = header

    # -----------------------------------------------------------------

    def __str__(self):

        """Return the HTML code for the table row and its cells as a string."""

        row = []

        if self.wrap: start = "<tr>"
        else: start = '<tr style="white-space:nowrap;">'

        row.append(start)
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

    def __init__(self, rows, header=None, css_class=None, bgcolors=None, tostr_kwargs=None, subheader=None,
                 strike_rows=None, bold_rows=None, italic_rows=None, wrap=True):

        """
        Table constructor.
        Keyword arguments:
        rows -- iterable of SimpleTableRow
        header:
        css_class -- table CSS class
        subheader_row:
        :param strike_rows:
        :param bold_rows:
        :param italic_rows:
        :param wrap:
        """

        # Construct the rows
        if isinstance(rows[0], SimpleTableRow): self.rows = rows
        else:
            #if bgcolors is not None: self.rows = [SimpleTableRow(row, bgcolors=bcolors, tostr_kwargs=tostr_kwargs) for row, bcolors in zip(rows, bgcolors)]
            #else: self.rows = [SimpleTableRow(row, tostr_kwargs=tostr_kwargs) for row in rows]

            self.rows = []
            for index in range(len(rows)):

                bcolors = bgcolors[index] if bgcolors is not None else None
                strike = strike_rows[index] if strike_rows is not None else False
                bold = bold_rows[index] if bold_rows is not None else False
                italic = italic_rows[index] if italic_rows is not None else False

                row = SimpleTableRow(rows[index], bgcolors=bcolors, strike=strike, bold=bold, italic=italic, wrap=wrap, tostr_kwargs=tostr_kwargs)
                self.rows.append(row)

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
    def one_row(cls, *cells, **kwargs):

        """
        This function ...
        :param cells:
        :param kwargs:
        :return:
        """

        # Make table
        return cls([cells], **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def one_column(cls, *cells, **kwargs):

        """
        Thisf unction ...
        :param cells:
        :param kwargs:
        :return:
        """

        # Make table
        return cls([[cell] for cell in cells], **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def rasterize(cls, cells, ncolumns=2, header=None, css_class=None, tostr_kwargs=None, subheader=None, wrap=True):

        """
        This function ...
        :param cells:
        :param ncolumns:
        :param header_row:
        :param css_class:
        :param tostr_kwargs:
        :param header:
        :param wrap:
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
        return cls(rows, header=header, css_class=css_class, tostr_kwargs=tostr_kwargs, subheader=subheader, wrap=wrap)

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
    def from_table(cls, table, css_class=None, bgcolors=None, tostr_kwargs=None, column_names=None, extra_columns=None,
                   extra_column_labels=None, strike_rows=None, bold_rows=None, italic_rows=None, wrap=True):

        """
        This function ...
        :param table:
        :param css_class:
        :param bgcolors:
        :param tostr_kwargs:
        :param column_names:
        :param extra_columns:
        :param extra_column_labels:
        :param strike_rows:
        :param bold_rows:
        :param italic_rows:
        :param wrap:
        :return:
        """

        # Set column names
        if column_names is None: column_names = table.column_names

        # Check
        if extra_column_labels is not None and extra_columns is None: raise ValueError("Extra column labels are specified but extra columns are not given")

        # Get rows
        #rows = table.as_tuples(add_units=False)
        rows = table.as_lists(add_units=False)
        if extra_columns is not None:
            for column in extra_columns:
                for i in range(len(rows)): rows[i].append(column[i])

        # Set header
        header = column_names
        if extra_columns is not None:
            for column_label in extra_column_labels:
                label = column_label if column_label is not None else ""
                header.append(label)

        # Set subheader
        subheader = table.unit_strings
        if extra_columns is not None:
            for column in extra_columns: subheader.append("")

        # Create the HTML table
        return cls(rows, header=header, subheader=subheader, css_class=css_class, bgcolors=bgcolors,
                   tostr_kwargs=tostr_kwargs, strike_rows=strike_rows, bold_rows=bold_rows, italic_rows=italic_rows, wrap=wrap)

    # -----------------------------------------------------------------

    @classmethod
    def from_composites(cls, composites, css_class=None, bgcolors=None, tostr_kwargs=None, labels=None, label="-",
                        extra_columns=None, extra_column_labels=None, strike_rows=None, bold_rows=None, italic_rows=None,
                        wrap=True):

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
        :param strike_rows:
        :param bold_rows:
        :param italic_rows:
        :param wrap:
        :return:
        """

        # Create the table
        table = SmartTable.from_composites(*composites, labels=labels, label=label)

        # Create the HTML table
        return cls.from_table(table, css_class=css_class, bgcolors=bgcolors, tostr_kwargs=tostr_kwargs,
                              column_names=table.descriptions, extra_columns=extra_columns,
                              extra_column_labels=extra_column_labels, strike_rows=strike_rows, bold_rows=bold_rows,
                              italic_rows=italic_rows, wrap=wrap)

    # -----------------------------------------------------------------

    @classmethod
    def from_composite(cls, composite, css_class=None, bgcolors=None, tostr_kwargs=None, key_label="Property",
                       value_label="Value"):

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

        return cls(composite.as_tuples(), header=[key_label, value_label], css_class=css_class, bgcolors=bgcolors, tostr_kwargs=tostr_kwargs)

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

    def __init__(self, title, css=None, css_path=None, style=None, encoding="utf-8", language="en", head=None,
                 footing=None, body_settings=None, javascript=None, javascript_path=None, javascript_path_body=None,
                 javascript_header=None):

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
        :param javascript:
        :param javascript_path:
        :param javascript_path_body:
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
        self.javascript = javascript
        self.javascript_header = javascript_header
        self.javascript_path = javascript_path
        self.javascript_path_body = javascript_path_body

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

        # Add javascript in header
        if self.javascript_header is not None:
            lines.append(make_script(self.javascript_header))

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

        # Add custom javascript
        if self.javascript is not None: lines.append(make_script(self.javascript))

        # Add javascript body scripts
        if self.javascript_path_body is not None:

            if types.is_string_sequence(self.javascript_path_body):
                for url in self.javascript_path_body:
                    lines.append('<script src="' + url + '"></script>')
            elif types.is_string_type(self.javascript_path_body): lines.append('<script src="' + self.javascript_path_body + '"></script>')
            else: raise ValueError("Invalid value for 'javascript_path_body'")

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
