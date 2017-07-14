#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
from pts.modeling.maps.component import get_old_stellar_map_names, get_young_stellar_map_names, get_ionizing_stellar_map_names
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

default_sfr = 1.
default_fuv_attenuation = -2.5 * np.log(1./2.) # half of the light in the FUV band is attenuated
default_ionizing_contribution = 0.5

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The name
definition.add_required("name", "string", "name of the model")

# Stellar maps
definition.add_required("old_stars_map", "string", "choice of old stars map", choices=get_old_stellar_map_names(modeling_path))
definition.add_required("young_stars_map", "string", "choice of young stars map", choices=get_young_stellar_map_names(modeling_path))
definition.add_required("ionizing_stars_map", "string", "choice of ionizing stars map", choices=get_ionizing_stellar_map_names(modeling_path))

# Flags
definition.add_flag("bulge", "add bulge", True)
definition.add_flag("old", "add old stars", True)
definition.add_flag("young", "add young stars", True)
definition.add_flag("ionizing", "add ionizing stars", True)
definition.add_flag("additional", "add additional stellar component(s)", True)

# Output directory
definition.add_optional("output", "directory_path", "output directory")

# Template
default_stellar_template = "BruzualCharlot"
definition.add_optional("default_old_bulge_template", "string", "SED template for the old stellar bulge", default_stellar_template)
definition.add_optional("default_old_disk_template", "string", "SED template for the old stellar disk", default_stellar_template)
definition.add_optional("default_young_template", "string", "SED template for the young stellar population", default_stellar_template)

# Metallicities
definition.add_optional("default_old_bulge_metallicity", "real", "default metallicity for the old stellar bulge", 0.03)
definition.add_optional("default_old_disk_metallicity", "real", "default metallicity for the old stellar disk", 0.03)
definition.add_optional("default_young_metallicity", "real", "default metallicity for the young stellar component", 0.03)
definition.add_optional("default_ionizing_metallicity", "real", "default metallicity for the ionizing stellar component", 0.03) # XU KONG et al. 2000

# Ages
definition.add_optional("default_old_bulge_age", "positive_real", "default age of the old stellar bulge (in Gyr)", 8.)
definition.add_optional("default_old_disk_age", "positive_real", "default age of the old stellar disk (in Gyr)", 8.)
definition.add_optional("default_young_age", "positive_real", "default age of the young stellar component (in Gyr)", 0.1)

# Scale heights
degeyter_ratio = 8.26
mosenkov_ratio = None
definition.add_optional("scalelength_to_scaleheight", "real", "ratio of scalelength to scaleheight", default=mosenkov_ratio, suggestions=[mosenkov_ratio, degeyter_ratio])
definition.add_optional("young_scaleheight_ratio", "real", "ratio of the young stellar scaleheight to the old stellar scaleheight", 0.5)
definition.add_optional("ionizing_scaleheight_ratio", "real", "ratio of the ionizing scaleheight to the old stellar scaleheight", 0.25)

## OTHER
definition.add_optional("default_sfr", "real", "default value for the SFR", default_sfr)
definition.add_optional("fuv_attenuation", "real", "estimate of the average FUV band attenuation", default_fuv_attenuation)
definition.add_optional("fuv_ionizing_contribution", "real", "relative contribution to the FUV flux from ionizing stars", default_ionizing_contribution)

# Star formation (Mappings)
definition.add_optional("default_ionizing_compactness", "real", "compactness", 6.)
definition.add_optional("default_ionizing_pressure", "quantity", "pressure", "1e12 K/m3", convert_default=True)
definition.add_optional("default_covering_factor", "real", "covering factor", 0.2)

# -----------------------------------------------------------------
