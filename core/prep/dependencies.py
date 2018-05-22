#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.dependencies Contains the DependenciesChecker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from pts.do.commandline import show_all_available, show_possible_matches
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

class DependenciesChecker(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DependenciesChecker, self).__init__(*args, **kwargs)

        # Dependencies
        self.dependencies = None
        self.dependencies_for_subproject = dict()
        self.data = dict()

        # Present / not present
        self.present_dependencies = []
        self.not_present_dependencies = []

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the dependencies
        self.get_dependencies()

        # 3. Check with installed packages
        self.check_dependencies()

        # 4. Show
        if self.config.show: self.show()

        # 5. Write
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DependenciesChecker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def get_dependencies(self):

        """
        This function ...
        :return:
        """

        # If no script name is given, list all dependencies of PTS and the
        # PTS modules that use them
        if self.config.script is None:

            if self.config.subprojects:

                self.dependencies = defaultdict(set)

                directories_for_subproject = defaultdict(list)

                # Loop over the subdirectories of the 'do' directory
                for path, name in fs.directories_in_path(introspection.pts_do_dir, returns=["path", "name"]):
                    subproject = name
                    directories_for_subproject[subproject].append(path)

                # Loop over the other directories in 'pts' (other than 'do' and 'doc')
                for path, name in fs.directories_in_path(introspection.pts_package_dir, returns=["path", "name"],
                                                         exact_not_name=["do", "doc"]):
                    subproject = name
                    directories_for_subproject[subproject].append(path)

                encountered_internal_modules = set()

                for subproject in directories_for_subproject:

                    # print(subproject, directories_for_subproject[subproject])

                    # List the dependencies of the matching script
                    dependencies_for_this = defaultdict(set)

                    for dir_path in directories_for_subproject[subproject]:

                        for file_path in fs.files_in_path(dir_path, extension="py", recursive=True):
                            # if subproject == "dustpedia": print(file_path)

                            introspection.add_dependencies(dependencies_for_this, file_path,
                                                           encountered_internal_modules)

                    self.dependencies_for_subproject[subproject] = dependencies_for_this.keys()

                    for dependency in dependencies_for_this:
                        for name in dependencies_for_this[dependency]:
                            self.dependencies[dependency].add(name)

            # No subprojects are specified, list all PTS dependencies
            else: self.dependencies = introspection.get_all_dependencies()

        # If a script name is given
        else:

            scripts = introspection.get_scripts()
            tables = introspection.get_arguments_tables()

            # Find matching 'do' commands (actual scripts or tabulated commands)
            matches = introspection.find_matches_scripts(self.config.script, scripts)
            table_matches = introspection.find_matches_tables(self.config.script, tables)

            # List the dependencies of the matching script
            self.dependencies = defaultdict(set)

            # No match
            if len(matches) + len(table_matches) == 0:
                show_all_available(scripts, tables)
                exit()

            # More matches
            elif len(matches) + len(table_matches) > 1:
                show_possible_matches(matches, table_matches, tables)
                exit()

            # Exactly one match from existing do script
            elif len(matches) == 1 and len(table_matches) == 0:

                # Determine the full path to the matching script
                script_path = fs.join(introspection.pts_do_dir, matches[0][0], matches[0][1])

                introspection.add_dependencies(self.dependencies, script_path, set())

            # Exactly one match from tabulated command
            elif len(table_matches) == 1 and len(matches) == 0:

                # from pts.core.tools import logging
                # configuration

                # Path to class module

                table_match = table_matches[0]
                subproject = table_match[0]
                index = table_match[1]

                relative_class_module_path = tables[subproject]["Path"][index].replace(".", "/").rsplit("/", 1)[0] + ".py"
                class_module_path = fs.join(introspection.pts_subproject_dir(subproject), relative_class_module_path)
                logging_path = fs.join(introspection.pts_package_dir, "core", "basics", "log.py")

                command_name = tables[subproject]["Command"][index]
                configuration_name = tables[subproject]["Configuration"][index]
                if configuration_name == "--": configuration_name = command_name
                configuration_module_path = fs.join(introspection.pts_root_dir,
                                                    "pts/" + subproject + "/config/" + configuration_name + ".py")

                # Add dependencies
                encountered = set()
                introspection.add_dependencies(self.dependencies, logging_path, encountered)
                introspection.add_dependencies(self.dependencies, configuration_module_path, encountered)
                introspection.add_dependencies(self.dependencies, class_module_path, encountered)

    # -----------------------------------------------------------------

    def check_dependencies(self):

        """
        This function ...
        :return:
        """

        # Get the names and versions of all installed python packages
        packages = introspection.installed_python_packages()

        # Get the names of the packages available through Canopy
        if self.config.canopy:

            url = "https://www.enthought.com/products/canopy/package-index/"

            # Try importing lxml
            try: from lxml import html
            except Exception:
                raise RuntimeError("You need the 'lxml' package to be able to specify the '--canopy' option")

            # Try importing requests
            try: import requests
            except Exception:
                raise RuntimeError("You need the 'requests' pacakge to be able to specify the '--canopy' option")

            page_as_string = requests.get(url).content
            tree = html.fromstring(page_as_string)

            tables = [e for e in tree.iter() if e.tag == 'table']
            table = tables[-1]

            table_rows = [e for e in table.iter() if e.tag == 'tr']
            column_headings = [e.text_content() for e in table_rows[0].iter() if e.tag == 'th']

            canopy_packages = []

            for row in table_rows[1:]:
                name = row[0].text_content().strip()
                canopy_packages.append(name.lower())

        else: canopy_packages = None

        # Loop over the packages and check their presence
        for dependency in sorted(self.dependencies, key=str.lower):

            # Get the list of PTS scripts for this dependency
            script_list = self.dependencies[dependency]

            # print(script_list)

            # Skip packages from the standard library, unless the appropriate flag is enabled
            if introspection.is_std_lib(dependency) and not self.config.standard: continue

            # Check presency and version
            if dependency in packages:
                present = True
                version = packages[dependency] if self.config.version else None
            else:
                present = introspection.is_present_package(dependency)
                version = None

            if present: self.present_dependencies.append(dependency)
            else: self.not_present_dependencies.append(dependency)

            # Check whether the package is available through Canopy
            in_canopy = dependency.lower() in canopy_packages if canopy_packages is not None else None

            # Check whether the package is available through pip
            if self.config.pip:

                from pts.core.tools.pypi import search
                results = list(search(dependency))
                in_pip = False
                for result in results:
                    if result["name"] == dependency:
                        in_pip = True
                        break
            else: in_pip = None

            # Check whether the package is available through conda
            if self.config.conda:

                from ..tools import conda
                results = conda.search(dependency)

                if len(results) > 0: in_conda = True
                else: in_conda = False

            else: in_conda = False

            # Get description
            if self.config.description:

                from pts.core.tools.pypi import search
                results = list(search(dependency))
                description = None
                for result in results:
                    if result["name"] == dependency:
                        description = result["summary"]
                        if description.lower().startswith(dependency.lower() + ": "): description = \
                        description.split(": ")[1]
                        break
            else: description = None

            self.data[dependency] = [script_list, present, version, in_canopy, in_pip, in_conda, description]

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        if self.config.subprojects:

            # dependencies_for_subproject = defaultdict(set)
            # for dependency in data:
            #    script_list = data[dependency][0]
            #    for script_path in script_list:
            #        subproject = script_path.split("/pts/")[1].split("/")[0]
            #        if subproject == "do":
            #            subproject = script_path.split("do/")[1].split("/")[0]
            #        dependencies_for_subproject[subproject].add(dependency)

            print("")

            for subproject in self.dependencies_for_subproject:

                print(fmt.underlined + subproject + fmt.reset + ":")
                print("")

                for dependency in self.dependencies_for_subproject[subproject]:

                    if dependency not in self.data: continue  # dependencies that were not added to the data because they belong to the standard library

                    script_list = self.data[dependency][0]
                    present = self.data[dependency][1]
                    version = self.data[dependency][2]
                    in_canopy = self.data[dependency][3]
                    in_pip = self.data[dependency][4]
                    in_conda = self.data[dependency][5]
                    description = self.data[dependency][6]

                    if present: print("   " + fmt.bold + fmt.green + dependency + fmt.reset)
                    else: print("   " + fmt.bold + fmt.red + dependency + fmt.reset)

                print("")

        else:

            print("")

            if len(self.present_dependencies) > 0:

                print(fmt.bold + fmt.green + "present:" + fmt.reset)
                print("")
                for dependency in self.present_dependencies: print_dependency(dependency, self.data, self.config)
                print("")

            if len(self.not_present_dependencies) > 0:

                print(fmt.bold + fmt.red + "not present:" + fmt.reset)
                print("")

                for dependency in self.not_present_dependencies: print_dependency(dependency, self.data, self.config)

                print("")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        self.write_report()

    # -----------------------------------------------------------------

    def write_report(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

def print_dependency(dependency, data, config):

    """
    This function ...
    :param dependency:
    :param data:
    :param config:
    :return:
    """

    from ..tools import types, stringify

    script_list = data[dependency][0]
    version = data[dependency][2]
    in_canopy = data[dependency][3]
    in_pip = data[dependency][4]
    in_conda = data[dependency][5]
    description = data[dependency][6]

    info = []
    if config.version: info.append(("version", version))
    if config.canopy: info.append(("in canopy", in_canopy))
    if config.pip: info.append(("in pip", in_pip))
    if config.conda: info.append(("in conda", in_conda))
    if config.description: info.append(("description", description))

    had_indent = False

    if len(info) == 0: print("   " + fmt.underlined + dependency + fmt.reset)
    elif len(info) == 1:
        info_entry = info[0]
        if info_entry[1] is not None:
            if types.is_boolean_type(info_entry[1]): print("   " + fmt.underlined + dependency + fmt.reset + " (" + info_entry[0] + ": " + stringify.yes_or_no(info_entry[1]) + ")")
            else: print("   " + fmt.underlined + dependency + fmt.reset + " (" + info_entry[0] + ": " + str(info_entry[1]) + ")")
    else:
        print("   " + fmt.underlined + dependency + fmt.reset)
        print("")
        had_indent = True
        for info_entry in info:
            if info_entry[1] is not None:
                if types.is_boolean_type(info_entry[1]): print("    - " + info_entry[0] + ": " + stringify.yes_or_no(info_entry[1]))
                else: print("    - " + info_entry[0] + ": " + str(info_entry[1]))

    # List the PTS modules that have this dependency
    if config.modules:
        print("")
        had_indent = True
        for script in script_list: print("    * " + script.split("PTS/pts/")[1])

    if had_indent: print("")

# -----------------------------------------------------------------
