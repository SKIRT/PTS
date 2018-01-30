#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.modeling.maps.selection import ComponentMapsSelection
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.decomposition.decomposition import scalelength_scaleheight_ratios, degeyter_ratio
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
selection = ComponentMapsSelection.from_modeling_path(modeling_path)

# -----------------------------------------------------------------

solar_metallicity = 0.02

# -----------------------------------------------------------------

default_sfr = 1.
default_fuv_attenuation = -2.5 * np.log(1./2.) # half of the light in the FUV band is attenuated
default_ionizing_contribution = 0.3

# -----------------------------------------------------------------

definition = definition.copy()

# The name
definition.add_required("name", "string", "name of the model")

# Stellar maps
definition.add_required("old_stars_map", "string", "choice of old stars map", choices=selection.old_map_names)
definition.add_required("young_stars_map", "string", "choice of young stars map", choices=selection.young_map_names)
definition.add_required("ionizing_stars_map", "string", "choice of ionizing stars map", choices=selection.ionizing_map_names)

# Flags
definition.add_flag("bulge", "add bulge", True)
definition.add_flag("old", "add old stars", True)
definition.add_flag("young", "add young stars", True)
definition.add_flag("ionizing", "add ionizing stars", True)
definition.add_flag("additional", "add additional stellar component(s)", False)

# Output directory
definition.add_optional("output", "directory_path", "output directory")

# Template
default_stellar_template = "BruzualCharlot"
definition.add_optional("default_old_bulge_template", "string", "SED template for the old stellar bulge", default_stellar_template)
definition.add_optional("default_old_disk_template", "string", "SED template for the old stellar disk", default_stellar_template)
definition.add_optional("default_young_template", "string", "SED template for the young stellar population", default_stellar_template)

# Metallicities
definition.add_optional("default_old_bulge_metallicity", "real", "default metallicity for the old stellar bulge", solar_metallicity)
definition.add_optional("default_old_disk_metallicity", "real", "default metallicity for the old stellar disk", solar_metallicity)
definition.add_optional("default_young_metallicity", "real", "default metallicity for the young stellar component", solar_metallicity)
definition.add_optional("default_ionizing_metallicity", "real", "default metallicity for the ionizing stellar component", solar_metallicity)

# Ages
definition.add_optional("default_old_bulge_age", "positive_real", "default age of the old stellar bulge (in Gyr)", 8.)
definition.add_optional("default_old_disk_age", "positive_real", "default age of the old stellar disk (in Gyr)", 8.)
definition.add_optional("default_young_age", "positive_real", "default age of the young stellar component (in Gyr)", 0.1)

# Scale heights
definition.add_optional("scalelength_to_scaleheight", "real", "ratio of scalelength to scaleheight", default=degeyter_ratio, suggestions=scalelength_scaleheight_ratios)
definition.add_optional("young_scaleheight_ratio", "real", "ratio of the young stellar scaleheight to the old stellar scaleheight", 0.5)
definition.add_optional("ionizing_scaleheight_ratio", "real", "ratio of the ionizing scaleheight to the old stellar scaleheight", 0.25)

## OTHER
definition.add_optional("default_sfr", "real", "default value for the SFR", default_sfr)
definition.add_optional("fuv_attenuation", "real", "estimate of the average FUV band attenuation", default_fuv_attenuation)

# Star formation (Mappings)
definition.add_optional("default_ionizing_compactness", "real", "compactness", 6.)
definition.add_optional("default_ionizing_pressure", "quantity", "pressure", "1e12 K/m3", convert_default=True)
definition.add_optional("default_covering_factor", "real", "covering factor", 0.2)

# -----------------------------------------------------------------

# ADVANCED: FLAG TO USE DEFAULTS
definition.add_flag("use_defaults", "use the defaults (do not prompt for component parameters)", False)

# -----------------------------------------------------------------
