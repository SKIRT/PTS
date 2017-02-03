#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.imports Contains the ImportsChecker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import importlib

# Import the relevant PTS classes and modules
from ..tools import introspection
from ..tools.logging import log
from ..basics.configurable import Configurable

# -----------------------------------------------------------------

class ImportsChecker(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        :param config:
        """

        # Call the constructor of the base class
        super(ImportsChecker, self).__init__(config)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Check internal imports
        self.check_internal()

        # 3. Check external imports
        self.check_external()

        # 4. Show
        if self.config.show: self.show()

        # 5. Write
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def check_internal(self):

        """
        This function ...
        :return:
        """

        # Get import lines for all modules
        import_lines = introspection.get_internal_imports()

        # Loop over the PTS modules
        for filepath in import_lines:

            lines = import_lines[filepath]

            for line in lines:

                which, unresolved = introspection.get_modules(line, filepath, return_unresolved=True)

                #print(which)

                # Loop over the modules that were found
                for module_path in which:

                    # Only look for internal imports
                    if "PTS/pts" not in module_path: continue

                    # print(module_path)
                    pythonic_path = module_path.split("PTS/")[1].replace("/", ".")[:-3]
                    # print(pythonic_path)

                    subproject = pythonic_path.split("pts.")[1].split(".")[0]

                    if pythonic_path.split(subproject + ".")[1].split(".")[0] == "config": continue

                    # Skip module initialization files, they are empty anyways
                    if pythonic_path.endswith("__init__"): continue
                    if pythonic_path.endswith("__main__"): continue
                    if pythonic_path.endswith("run_queue"): continue
                    if pythonic_path.endswith("enable_qch_mathjax"): continue
                    if pythonic_path.endswith("eagle.config"): continue
                    if pythonic_path.endswith("eagle.collections"): continue
                    if pythonic_path.endswith("eagle.database"): continue
                    if pythonic_path.endswith("eagle.galaxy"): continue
                    if pythonic_path.endswith("eagle.plotresults"): continue
                    if pythonic_path.endswith("eagle.runner"): continue
                    if pythonic_path.endswith("eagle.scheduler"): continue
                    if pythonic_path.endswith("eagle.skirtrun"): continue
                    if pythonic_path.endswith("fit2BB_Md"): continue

                    try:
                        module = importlib.import_module(pythonic_path)
                    except ImportError:
                        log.warning("Importing module '" + pythonic_path + "' failed")
                        continue

                    # Loop over the names imported for this module
                    if which[module_path] is not None:

                        for name in which[module_path]:
                            ok = hasattr(module, name)
                            if not ok: log.warning("Name '" + name + "' could not be imported from module '" + pythonic_path + "'")

                if len(unresolved) > 0:
                    print("Unresolved imports:")
                    print(unresolved)

    # -----------------------------------------------------------------

    def check_external(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------
