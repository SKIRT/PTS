#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
#from pts.modeling.component.component import get_cache_host_id
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

definition = definition.copy()

#definition.add_optional("image_names", "string_list", "The image names for which to make the truncated images", ["GALEX FUV", "SDSS g", "IRAC I1", "MIPS 24mu", "Pacs red", "SPIRE PLW"])
definition.add_optional("factor_range", "real_range", "range (min,max) of values for the factor to which the disk ellipse should be scaled to obtain the truncation region", "0.2>1.5", convert_default=True)
definition.add_optional("factor_nvalues", "integer", "the number of values for the factor", 10)
#definition.add_optional("best_factor", "real", "the best estimate for the value of the factor", 0.82)

definition.add_optional("nbins", "positive_integer", "number of bins for signal to noise curve", 20)

# Cache
#cache_host_id = get_cache_host_id(modeling_path)
#if cache_host_id is not None: definition.add_flag("cache", "cache FITS files to the remote host storage", False)
#else: definition.add_fixed("cache", "caching not possible since cache host id not defined", False)

# Value
definition.add_optional("truncated_value", "real", "value to replace truncated pixels with", default="nan", convert_default=True)

# -----------------------------------------------------------------
