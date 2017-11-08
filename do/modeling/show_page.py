#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.show_page Show any page

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments, prompt_finish
from pts.core.tools import browser
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.html.component import data_page_filename, photometry_page_filename, components_page_filename
from pts.modeling.html.component import preparation_page_filename, maps_page_filename, model_page_filename
from pts.modeling.html.component import fitting_page_filename, seds_page_filename, datacubes_page_filename
from pts.modeling.html.component import fluxes_page_filename, images_page_filename, attenuation_page_filename
from pts.modeling.html.component import colours_page_filename, heating_page_filename
from pts.modeling.truncation.component import html_name, ellipse_page_filename, significance_page_filename

# -----------------------------------------------------------------

# Strings to identify different pages
pages = ["index", "status", "data", "preparation", "components", "photometry", "maps", "model", "fitting", "seds", "datacubes", "fluxes", "images", "attenuation", "colours", "heating"]
maps_pages = ["all", "summary", "old", "young", "ionizing", "dust", "clip", "selection"]
truncation_pages = ["ellipse", "levels"]

# -----------------------------------------------------------------

# Set the page names
page_names = []
page_names.extend(pages)
page_names.extend("maps/" + name for name in maps_pages)
page_names.extend("truncation/" + name for name in truncation_pages)

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("page", "string", "page to show", choices=page_names)
config = parse_arguments("show_page", definition)

# -----------------------------------------------------------------

# Get modeling path
environment = load_modeling_environment_cwd()

# -----------------------------------------------------------------

# Determine the path to the HTML file
if config.page.startswith("maps/"):

    page = config.page.split("maps/")[1]

    # Determine the page path
    if page == "all": page_path = environment.all_maps_html_page_path
    elif page == "summary": page_path = environment.maps_summary_html_page_path
    elif page == "old": page_path = environment.old_maps_html_page_path
    elif page == "young": page_path = environment.young_maps_html_page_path
    elif page == "ionizing": page_path = environment.ionizing_maps_html_page_path
    elif page == "dust": page_path = environment.dust_maps_html_page_path
    elif page == "clip": page_path = environment.clip_maps_html_page_path
    elif page == "selection": page_path = environment.maps_selection_html_page_path
    else: raise ValueError("Invalid maps page name")

    # Set directory
    directory = environment.maps_html_path

# Truncation page
elif config.page.startswith("truncation/"):

    page = config.page.split("truncation/")[1]
    truncation_html_path = fs.create_directory_in(environment.truncation_path, html_name)

    # Determine the page path
    if page == "ellipse": page_path = fs.join(truncation_html_path, ellipse_page_filename)
    elif page == "levels": page_path = fs.join(truncation_html_path, significance_page_filename)
    else: raise ValueError("Invalid truncation page name")

    # Set directory
    directory = environment.truncation_html_path

# Other
else:

    # Determine the page path
    if config.page == "index": page_path = environment.html_index_path
    elif config.page == "status": page_path = environment.html_status_path
    elif config.page == "data": page_path = fs.join(environment.html_path, data_page_filename)
    elif config.page == "preparation": page_path = fs.join(environment.html_path, photometry_page_filename)
    elif config.page == "components": page_path = fs.join(environment.html_path, components_page_filename)
    elif config.page == "photometry": page_path = fs.join(environment.html_path, preparation_page_filename)
    elif config.page == "maps": page_path = fs.join(environment.html_path, maps_page_filename)
    elif config.page == "model": page_path = fs.join(environment.html_path, model_page_filename)
    elif config.page == "fitting": page_path = fs.join(environment.html_path, fitting_page_filename)
    elif config.page == "seds": page_path = fs.join(environment.html_path, seds_page_filename)
    elif config.page == "datacubes": page_path = fs.join(environment.html_path, datacubes_page_filename)
    elif config.page == "fluxes": page_path = fs.join(environment.html_path, fluxes_page_filename)
    elif config.page == "images": page_path = fs.join(environment.html_path, images_page_filename)
    elif config.page == "attenuation": page_path = fs.join(environment.html_path, attenuation_page_filename)
    elif config.page == "colours": page_path = fs.join(environment.html_path, colours_page_filename)
    elif config.page == "heating": page_path = fs.join(environment.html_path, heating_page_filename)
    else: raise ValueError("Invalid page name")

    directory = None

# -----------------------------------------------------------------

# Check whether page exists
if not fs.is_file(page_path): raise ValueError("The page is not present")

# -----------------------------------------------------------------

# # Setup localhost server
# with browser.serve_local_host():
#     # Open
#     # SHOULD BE AN OPEN AND WAIT FOR CLOSING FUNCTION?
#     browser.open_path(page_path)

if directory is not None: fs.change_cwd(directory)

import subprocess
from pts.core.tools import introspection

python_path = introspection.python_executable_path()

command = [python_path, "-m", "SimpleHTTPServer"]
process = subprocess.Popen(command)

#thread = browser.start_localhost()
browser.open_path(page_path)
prompt_finish()

# Kill the server
process.kill()

# -----------------------------------------------------------------
