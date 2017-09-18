#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.view View an image with JS9.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.view.html import javascripts, css_scripts, JS9Preloader, JS9Window, scales, colormaps, zooms
from pts.magic.view.html import body_settings, make_replace_infs_by_nans, make_load_regions, make_load_region
from pts.core.tools import html
from pts.core.tools import filesystem as fs
from pts.core.tools import browser
from pts.magic.region.list import load_region_list

# -----------------------------------------------------------------

default_scale = "log"
default_colormap = "viridis"
default_zoom = "toFit"

# -----------------------------------------------------------------

# Cretae configuration definition
definition = ConfigurationDefinition()
definition.add_required("image", "file_path", "image path")
definition.add_optional("regions", "file_path", "regions file path")
definition.add_optional("mask", "file_path", "mask file path")
definition.add_optional("page_width", "positive_integer", "page width (in pixels)", 600)
definition.add_optional("with", "positive_integer", "image width", 500)
definition.add_optional("height", "positive_integer", "image height")
definition.add_flag("menubar", "show menu bar", True)
definition.add_flag("colorbar", "show color bar", True)
definition.add_flag("resize", "allow resize", True)
definition.add_flag("scrolling", "allow scrolling", True)
definition.add_flag("onload", "display the masks and regions when the image is loaded", False) # doesn't work yet

# Additional features
definition.add_flag("panner", "add a panner", False)
definition.add_flag("magnifier", "add a magnifier", False)
definition.add_flag("combine_panner_and_magnifier", "combine panner and manifier next to each other", True)

# Regions settings
definition.add_flag("movable", "movable regions", True)
definition.add_flag("rotatable", "rotatable regions", True)
definition.add_flag("removable", "removable regions", True)
definition.add_flag("resizable", "resizable regions", True)

# View settings
definition.add_optional("scale", "string", "scale", default_scale, choices=scales)
definition.add_optional("colormap", "string", "color map", default_colormap, choices=colormaps)
definition.add_optional("zoom", "string", "zoom function", default_zoom, choices=zooms)

# Create config
config = parse_arguments("view", definition)

# -----------------------------------------------------------------

stylesheet_url = "http://users.ugent.be/~sjversto/stylesheet.css"
background_color = "white"
page_style = "ugentstyle"

# -----------------------------------------------------------------

# Create CSS for the page width
css = html.make_page_width(config.page_width)

# Get name
name = fs.strip_extension(fs.name(config.image))
title = name

# Create list of css scripts
css_paths = css_scripts[:]
css_paths.append(stylesheet_url)

# Create the page
page = html.HTMLPage(title, css=css, body_settings=body_settings, style=page_style, css_path=css_paths, javascript_path=javascripts, footing=html.updated_footing())

# Make theme button
classes = dict()
classes["JS9Menubar"] = "data-backgroundColor"
classes["JS9Menubar JS9Plugin"] = "data-backgroundColor"
page += html.center(html.make_theme_button(classes=classes))
page += html.newline

# -----------------------------------------------------------------

image_name = html.make_usable(name)
display_name = "JS9" + image_name

# -----------------------------------------------------------------

# Make settings
settings = dict()
settings["scale"] = config.scale
settings["colormap"] = config.colormap
settings["zoom"] = config.zoom

# Create preloader
preloader = JS9Preloader()
image = preloader.add_path(name, config.image, settings=settings, display=display_name)

# Create window
window = JS9Window(display_name, width=config.width, height=config.height,
                   background_color=background_color, menubar=config.menubar,
                   colorbar=config.colorbar, resize=config.resize, scrolling=config.scrolling,
                   panner=config.panner, magnifier=config.magnifier, combine_panner_and_magnifier=config.combine_panner_and_magnifier)
view = window.view

# -----------------------------------------------------------------

## INFS

# Create nan/infs replacer button
button_id = image_name + "nansinfs"
replace_function_name = "replace_infs_nans_" + image_name
replace_nans_infs = make_replace_infs_by_nans(display_name)

# Create the button
infs_button = html.make_script_button(button_id, "Replace infs", replace_nans_infs, replace_function_name)

# Add the button
page += html.center(infs_button)
page += html.newline

# -----------------------------------------------------------------

## REGIONS

if config.regions is not None:

    # Load the regions
    regions = load_region_list(config.regions)

    # Regions button
    region_button_id = image_name + "regionsbutton"
    load_region_function_name = "load_regions_" + image_name
    load_region = make_load_regions(regions, display=display_name, movable=config.movable, rotatable=config.rotatable, removable=config.removable, resizable=config.resizable, quote_character="'")

    # Create region loader
    regions_button = html.make_script_button(region_button_id, "Load regions", load_region, load_region_function_name)

    # Add the button
    page += html.center(regions_button)
    page += html.newline

# -----------------------------------------------------------------

# Add the window
page += window

# -----------------------------------------------------------------

# Add the preloader
page += preloader

# -----------------------------------------------------------------

# Open the page
browser.open_page(page)

# -----------------------------------------------------------------
