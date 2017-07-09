#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plotting.model Plot the model.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from IPython.display import display
#from ipywidgets import embed_snippet, embed_minimal_html
import matplotlib.pyplot as plt
#from IPython.terminal.embed import embed_snippet
from ipyvolume.transferfunction import TransferFunctionWidgetJs3
from ipyvolume.style import dark, light, minimal
from ipyvolume.style import create as create_style

import json
import ipywidgets

# Import astronomical modules
from astropy.io.fits import Header
import ipyvolume.pylab as p3
import ipyvolume
from ipyvolume.embed import embed_html, template, get_state, add_referring_widgets

# Import the relevant PTS classes and modules
from pts.core.tools.logging import log
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.modeling.basics.models import load_2d_model, load_3d_model
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.modeling.basics.models import DeprojectionModel3D
from pts.modeling.basics.properties import GalaxyProperties

# -----------------------------------------------------------------

body_template = """<script src="https://unpkg.com/jupyter-js-widgets@~2.0.20/dist/embed.js"></script>
<script type="application/vnd.jupyter.widget-state+json">
{json_data}
</script>
{widget_views}"""

# -----------------------------------------------------------------

widget_view_template = """<script type="application/vnd.jupyter.widget-view+json">
{{
    "model_id": "{model_id}"
}}
</script>"""

# -----------------------------------------------------------------

def xyz(shape=128, limits=[-3, 3], spherical=False, sparse=True, centers=False):

    """
    This function ...
    :param shape:
    :param limits:
    :param spherical:
    :param sparse:
    :param centers:
    :return:
    """

    dim = 3

    try: shape[0]
    except: shape = [shape] * dim

    try: limits[0][0]
    except: limits = [limits] * dim

    if centers: v = [slice(vmin+(vmax-vmin)/float(N)/2, vmax-(vmax-vmin)/float(N)/4, (vmax-vmin)/float(N)) for (vmin, vmax), N in zip(limits, shape)]
    else: v = [slice(vmin, vmax+(vmax-vmin)/float(N)/2, (vmax-vmin)/float(N-1)) for (vmin, vmax), N in zip(limits, shape)]

    if sparse: x, y, z = np.ogrid.__getitem__(v)
    else: x, y, z = np.mgrid.__getitem__(v)

    # RETURN
    if spherical:

        r = np.linalg.norm([x, y, z])
        theta = np.arctan2(y, x)
        phi = np.arccos(z / r)
        return x, y, z, r, theta, phi

    else: return x, y, z

# -----------------------------------------------------------------

def determine_model_limits(components, unit, symmetric=False):

    """
    This function ...
    :param components:
    :param unit:
    :param symmetric:
    :return:
    """

    x_min = 0.0
    x_max = 0.0
    y_min = 0.0
    y_max = 0.0
    z_min = 0.0
    z_max = 0.0

    # print("")
    # Determine limits, loop over components
    for name in components:

        component = components[name]
        x_min_scalar = component.xmin.to(unit).value
        x_max_scalar = component.xmax.to(unit).value
        y_min_scalar = component.ymin.to(unit).value
        y_max_scalar = component.ymax.to(unit).value
        z_min_scalar = component.zmin.to(unit).value
        z_max_scalar = component.zmax.to(unit).value

        if x_min_scalar < x_min: x_min = x_min_scalar
        if x_max_scalar > x_max: x_max = x_max_scalar
        if y_min_scalar < y_min: y_min = y_min_scalar
        if y_max_scalar > y_max: y_max = y_max_scalar
        if z_min_scalar < z_min: z_min = z_min_scalar
        if z_max_scalar > z_max: z_max = z_max_scalar

        # print(name + " limits: ")
        # print("")
        # print(" - x: " + str(component.xrange))
        # print(" - y: " + str(component.yrange))
        # print(" - z: " + str(component.zrange))
        # print("")

    minvalue = min(x_min, y_min, z_min)
    maxvalue = max(x_max, y_max, z_max)

    # Define limits
    if symmetric: limits = [minvalue, maxvalue]
    else: limits = [[x_min, x_max], [y_min, y_max], [z_min, z_max]]

    # Return the limits
    return limits

# -----------------------------------------------------------------

minimal_test = create_style("minimal", {
		'background-color': 'white',
        'box' : {
            'color': 'pink',
            'visible': False,
        },
        'axes': {
            'color': 'black',
            'visible': False,
            'x': {
                'color': '#f00',
                'label': {
                    'color': '#0f0'
                },
                'ticklabel': {
                    'color': '#00f'
                },
            },
            'y': {
                'color': '#0f0',
                'label': {
                    'color': '#00f'
                },
                'ticklabel': {
                    'color': '#f00'
                }
            },
            'z': {
                    'color': '#00f',
                    'label': {
                        'color': '#f00'
                    },
                    'ticklabel': {
                        'color': '#0f0'
                    }
            }
        }
})

# -----------------------------------------------------------------

def plot_galaxy_components(components, draw=True, show=True, shape=128, unit="pc", width=700, height=800, style='light', **kwargs):

    """
    This function ....
    :param components:
    :param draw:
    :param show:
    :param shape:
    :param unit:
    :param width:
    :param height:
    :param style:
    :param kwargs:
    :return:
    """

    # Detemrine the limits
    limits = determine_model_limits(components, unit, symmetric=True)

    # Debugging
    log.debug("Plot limits: " + str(limits))

    # Create coordinate data
    x, y, z, r, theta, phi = xyz(shape=shape, limits=limits, spherical=True)
    data = r * 0

    # Loop over the components
    for name in components:

        # Debugging
        log.debug("Computing the density of the " + name + " component ...")

        component = components[name]
        density = component.density_function(normalize=True)(x, y, z)
        data += density

    # DRAW FIGURE
    if draw:

        # :param lighting: boolean, to use lighting or not, if set to false, lighting parameters will be overriden
        # :param data_min: minimum value to consider for data, if None, computed using np.nanmin
        # :param data_max: maximum value to consider for data, if None, computed using np.nanmax
        # :param tf: transfer function (see ipyvolume.transfer_function, or use the argument below)
        # :param stereo: stereo view for virtual reality (cardboard and similar VR head mount)
        # :param width: width of rendering surface
        # :param height: height of rendering surface
        # :param ambient_coefficient: lighting parameter
        # :param diffuse_coefficient: lighting parameter
        # :param specular_coefficient: lighting parameter
        # :param specular_exponent: lighting parameter
        # :param downscale: downscale the rendering for better performance, for instance when set to 2, a 512x512 canvas will show a 256x256 rendering upscaled, but it will render twice as fast.
        # :param level: level(s) for the where the opacity in the volume peaks, maximum sequence of length 3
        # :param opacity: opacity(ies) for each level, scalar or sequence of max length 3
        # :param level_width: width of the (gaussian) bumps where the opacity peaks, scalar or sequence of max length 3
        # :param kwargs: extra argument passed to Volume and default transfer function

        # DEFAULT:

        # lighting=False, data_min=None, data_max=None, tf=None, stereo=False,
        # width=400, height=500,
        # ambient_coefficient=0.5, diffuse_coefficient=0.8,
        # specular_coefficient=0.5, specular_exponent=5,
        # downscale=1,
        # level=[0.1, 0.5, 0.9], opacity=[0.01, 0.05, 0.1], level_width=0.1,

        level = [0.2]
        #opacity = [0.05, 0.0, 0.0]
        opacity = [0.08, 0.0, 0.0]
        level_width = 0.2
        level_width = [level_width] * 3

        kwargs = dict()
        kwargs["width"] = width
        kwargs["height"] = height
        kwargs["stereo"] = False
        kwargs["level"] = level
        kwargs["opacity"] = opacity
        kwargs["level_width"] = level_width
        kwargs["downscale"] = 1

        # Create transfer function arguments
        tf_kwargs = {}

        # Clip off lists
        min_length = min(len(level), len(level_width), len(opacity))
        level = list(level[:min_length])
        opacity = list(opacity[:min_length])
        level_width = list(level_width[:min_length])
        # append with zeros
        while len(level) < 3:
            level.append(0)
        while len(opacity) < 3:
            opacity.append(0)
        while len(level_width) < 3:
            level_width.append(0)
        for i in range(1,4):
            tf_kwargs["level"+str(i)] = level[i-1]
            tf_kwargs["opacity"+str(i)] = opacity[i-1]
            tf_kwargs["width"+str(i)] = level_width[i-1]
        tf = TransferFunctionWidgetJs3(**tf_kwargs)

        # Set the transfer function
        kwargs["tf"] = tf

        # Set style
        if style == "dark": kwargs["style"] = dark
        elif style == "light": kwargs["style"] = light
        elif style == "minimal": kwargs["style"] = minimal
        else: raise ValueError("Invalid style: " + style)

        # Create the volume plot
        vol = ipyvolume.quickvolshow(data=data.T, **kwargs)
        #vol = p3.volshow(data=data, **kwargs)

        # SHOW?
        if show:
            #p3.volshow()
            vol.show()
        return vol

    # ONLY RETURN THE DATA
    else: return data

# -----------------------------------------------------------------

def generate_html(widgets, drop_defaults=False, all=False, title="ipyvolume embed example", template=template, widget_view_template=widget_view_template, as_dict=False, only_body=False, **kwargs):

    """
    This function ...
    :param widgets:
    :param drop_defaults:
    :param all:
    :param title:
    :param external_json:
    :param template:
    :param widget_view_template:
    :param as_dict:
    :param only_body:
    :param kwargs:
    :return:
    """

    try:
        widgets[0]
    except (IndexError, TypeError):
        widgets = [widgets]

    # collect the state of all relevant widgets
    state = {}
    previous = 0 + ipyvolume.serialize.performance
    try:
        # we cannot serialize binary buffers yet into the json format, so go back to json
        # style, and afterwards set it back
        ipyvolume.serialize.performance = 0
        if all:
            state = ipywidgets.Widget.get_manager_state(drop_defaults=drop_defaults)["state"]
        for widget in widgets:
            if not all:
                get_state(widget, state, drop_defaults=drop_defaults)
        # it may be that other widgets refer to the collected widgets, such as layouts, include those as well
        while add_referring_widgets(state):
            pass
    finally:
        ipyvolume.serialize.performance = previous

    values = dict(extra_script_head="", body_pre="", body_post="")
    values.update(kwargs)
    widget_views = ""
    for widget in widgets:
        widget_views += widget_view_template.format(**dict(model_id=widget.model_id))
    json_data = dict(version_major=1, version_minor=0, state=state)
    values.update(dict(title=title,
              json_data=json.dumps(json_data),
                   widget_views=widget_views))

    # Return
    if as_dict: return values
    elif only_body: return body_template.format(**values)
    else:
        html_code = template.format(**values)
        return html_code

# -----------------------------------------------------------------
