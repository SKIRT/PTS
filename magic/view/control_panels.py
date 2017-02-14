# HIDE_DEPENDS
from __future__ import absolute_import
from .source_panel import SourcePanel
from .plot_panel import PlotPanel
from .phot_panel import PhotPanel
from .stats_panel import StatsPanel
from .color_panel import ColorPanel
#from .fits_faker_panel.fits_faker_panel import FitsFakerPanel

control_panels_to_load = [("Source", SourcePanel),
                          ("Color", ColorPanel),
                          ("Plot", PlotPanel),
                          ("Stats", StatsPanel),
                          ("Phot", PhotPanel),
                          #("Faker", FitsFakerPanel)
                          ]
