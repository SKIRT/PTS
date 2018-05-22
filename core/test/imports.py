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
from ..basics.log import log
from ..basics.configurable import Configurable
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class ImportsChecker(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ImportsChecker, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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

                    # Skip
                    filename = fs.strip_extension(fs.name(module_path))
                    if introspection.skip_module(filename, module_path): continue

                    try: module = importlib.import_module(pythonic_path)
                    except ImportError:
                        log.warning("Importing module '" + pythonic_path + "' failed")
                        continue

                    # Loop over the names imported for this module
                    if which[module_path] is not None:

                        for name in which[module_path]:
                            ok = hasattr(module, name)
                            if not ok: log.warning("Name '" + name + "' could not be imported from module '" + pythonic_path + "'")

                if len(unresolved) == 1:
                    log.warning("Unresolved import: '" + unresolved[0][1] + "' in '" + unresolved[0][0] + "'")
                elif len(unresolved) > 1:
                    log.warning("Unresolved imports:")
                    for module, name in unresolved:
                        log.warning(" - '" + name + "' in '" + module + "'")

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
