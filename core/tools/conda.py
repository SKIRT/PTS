#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.conda Interacting with conda.
#  most functions come from: https://github.com/conda/conda-api,
#  with slight modifications.

# -----------------------------------------------------------------

# Import standard modules
import re
import os
import sys
import json
from subprocess import Popen, PIPE
from os.path import basename, isdir, join
from collections import defaultdict

# Import the relevant PTS classes and modules
from . import terminal
from . import filesystem as fs

# -----------------------------------------------------------------

class CondaError(Exception):
    "General Conda error"
    pass

# -----------------------------------------------------------------

class CondaEnvExistsError(CondaError):
    "Conda environment already exists"
    pass

# -----------------------------------------------------------------

def find_conda():

    """
    This function ...
    :return:
    """

    # Find conda path
    if terminal.is_existing_executable("conda"): conda_executable_path = terminal.executable_path("conda", no_pexpect=True)
    else: conda_executable_path = None

    # Find conda installation in the home directory
    if conda_executable_path is None:

        # Search for conda executables in the home directory
        #for path in fs.files_in_path(fs.home, exact_name="conda", extension=""):
        #    print(path)
        conda_path = fs.join(fs.home, "miniconda", "bin", "conda")
        if fs.is_file(conda_path): conda_executable_path = conda_path

    # Return the path to the conda executable
    return conda_executable_path

# -----------------------------------------------------------------

def is_environment(name, conda_path="conda"):

    """
    This function ...
    :param name:
    :param conda_path:
    :return:
    """

    return name in conda_environments(conda_path=conda_path)

# -----------------------------------------------------------------

def conda_environments(conda_path="conda"):

    """
    This function ...
    :return:
    """

    output = terminal.execute_no_pexpect(conda_path + " env list")
    envs = []
    for line in output:
        if line.startswith("#"): continue
        if line.strip() == "": continue
        envs.append(line.split()[0])
    return envs

# -----------------------------------------------------------------

def conda_active_environment(conda_path="conda"):

    """
    This function ...
    :return:
    """

    output = terminal.execute_no_pexpect(conda_path + " env list")
    env = None
    for line in output:
        if line.startswith("#"): continue
        if line.strip() == "": continue
        splitted = line.split()
        if splitted[1].strip() == "*":
            env = splitted[0].strip()
            break
    return env

# -----------------------------------------------------------------

def conda_environment_for_pts():

    """
    This function ...
    :return:
    """

    # Get environment name
    pts_alias = terminal.resolve_alias("pts")
    pts_python_path = pts_alias.split()[0]

    # If just 'python'
    #if pts_python_path == "python": raise Exception("Cannot determine the conda environment used for pts")
    if pts_python_path == "python": return None

    # Get the environment name
    env_path = fs.directory_of(fs.directory_of(pts_python_path))
    environment_name = fs.name(env_path)
    return environment_name

# -----------------------------------------------------------------

def conda_python_path_for_pts():

    """
    This function ...
    :return:
    """

    # Get environment name
    pts_alias = terminal.resolve_alias("pts")
    pts_python_path = pts_alias.split()[0]
    return pts_python_path

# -----------------------------------------------------------------

def activate_environment(environment_name, conda_path="conda", activate_path="activate"):

    """
    This function ...
    :return:
    """

    previous = conda_active_environment(conda_path)
    terminal.execute_no_pexpect("source " + activate_path + " " + environment_name)
    return previous

# -----------------------------------------------------------------

def deactivate(deactivate_path="deactivate"):

    """
    This function ...
    :param deactivate_path:
    :return:
    """

    terminal.execute_no_pexpect("source " + deactivate_path)

# -----------------------------------------------------------------

def is_present_package(name, environment_name=None, conda_path="conda"):

    """
    This function ...
    :param name:
    :param environment_name:
    :param conda_path:
    :return:
    """

    if environment_name is not None: command = conda_path + " list " + name + " --name " + environment_name
    else: command = conda_path + " list " + name

    output = terminal.execute_no_pexpect(command)

    # Find the package
    for line in output:
        if line.startswith("#"): continue
        if line.strip() == "": continue
        if line.split()[0].lower() == name.lower(): return True
    return False

# -----------------------------------------------------------------

def _call_conda(extra_args, conda_path="conda"):

    """
    This function ...
    :param extra_args:
    :param conda_path:
    :return:
    """

    cmd_list = [conda_path]
    cmd_list.extend(extra_args)

    try: p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
    except OSError: raise Exception("could not invoke %r\n" % extra_args)
    return p.communicate()

# -----------------------------------------------------------------

def _call_and_parse(extra_args, conda_path="conda"):

    """
    This function ...
    :param extra_args:
    :param conda_path:
    :return:
    """

    stdout, stderr = _call_conda(extra_args, conda_path=conda_path)
    if stderr.decode().strip():
        raise Exception('conda %r:\nSTDERR:\n%s\nEND' % (extra_args,
                                                         stderr.decode()))
    return json.loads(stdout.decode())

# -----------------------------------------------------------------

def _setup_install_commands_from_kwargs(kwargs, keys=tuple()):

    cmd_list = []
    if kwargs.get('override_channels', False) and 'channel' not in kwargs:
        raise TypeError('conda search: override_channels requires channel')

    if 'env' in kwargs:
        cmd_list.extend(['--name', kwargs.pop('env')])
    if 'prefix' in kwargs:
        cmd_list.extend(['--prefix', kwargs.pop('prefix')])
    if 'channel' in kwargs:
        channel = kwargs.pop('channel')
        if isinstance(channel, str):
            cmd_list.extend(['--channel', channel])
        else:
            cmd_list.append('--channel')
            cmd_list.extend(channel)

    for key in keys:
        if key in kwargs and kwargs[key]:
            cmd_list.append('--' + key.replace('_', '-'))

    return cmd_list

# -----------------------------------------------------------------

def get_conda_version(conda_path="conda"):

    """
    return the version of conda being used (invoked) as a string
    """

    pat = re.compile(r'conda:?\s+(\d+\.\d\S+|unknown)')
    stdout, stderr = _call_conda(['--version'], conda_path=conda_path)
    # argparse outputs version to stderr in Python < 3.4.
    # http://bugs.python.org/issue18920
    m = pat.match(stderr.decode().strip())
    if m is None:
        m = pat.match(stdout.decode().strip())

    if m is None:
        raise Exception('output did not match: %r' % stderr)
    return m.group(1)

# -----------------------------------------------------------------

def get_envs(conda_path="conda"):

    """
    Return all of the (named) environment (this does not include the root
    environment), as a list of absolute path to their prefixes.
    """

    info = _call_and_parse(['info', '--json'], conda_path=conda_path)
    return info['envs']

# -----------------------------------------------------------------

def get_prefix_envname(name):

    """
    Given the name of an environment return its full prefix path, or None
    if it cannot be found.
    """

    #if name == 'root':
    #    return ROOT_PREFIX
    for prefix in get_envs():
        if basename(prefix) == name:
            return prefix
    return None

# -----------------------------------------------------------------

def linked(prefix):

    """
    Return the (set of canonical names) of linked packages in `prefix`.
    """

    if not isdir(prefix):
        raise Exception('no such directory: %r' % prefix)
    meta_dir = join(prefix, 'conda-meta')
    if not isdir(meta_dir):
        # we might have nothing in linked (and no conda-meta directory)
        return set()
    return set(fn[:-5] for fn in os.listdir(meta_dir) if fn.endswith('.json'))

# -----------------------------------------------------------------

def split_canonical_name(cname):

    """
    Split a canonical package name into (name, version, build) strings.
    """

    return tuple(cname.rsplit('-', 2))

# -----------------------------------------------------------------

def info(conda_path="conda"):

    """
    Return a dictionary with configuration information.
    No guarantee is made about which keys exist.  Therefore this function
    should only be used for testing and debugging.
    """

    return _call_and_parse(['info', '--json'], conda_path=conda_path)

# -----------------------------------------------------------------

def package_info(package, conda_path="conda"):

    """
    Return a dictionary with package information.
    """

    return _call_and_parse(['info', package, '--json'], conda_path=conda_path)

# -----------------------------------------------------------------

def search(regex=None, spec=None, **kwargs):

    """
    Search for packages.
    """

    cmd_list = ['search', '--json']

    if regex and spec:
        raise TypeError('conda search: only one of regex or spec allowed')

    if regex:
        cmd_list.append(regex)

    if spec:
        cmd_list.extend(['--spec', spec])

    if 'platform' in kwargs:
        cmd_list.extend(['--platform', kwargs.pop('platform')])

    cmd_list.extend(
        _setup_install_commands_from_kwargs(
            kwargs,
            ('canonical', 'unknown', 'use_index_cache', 'outdated',
             'override_channels')))

    return _call_and_parse(cmd_list, conda_path=kwargs.pop("conda_path", "conda"))

# -----------------------------------------------------------------

def available_versions(regex):

    """
    This function ...
    :param regex:
    :return:
    """

    versions = defaultdict(set)

    results = search(regex)

    for label in results:

        for package_version in results[label]:

            version = package_version["version"]
            versions[label].add(version)

    return versions

# -----------------------------------------------------------------

def create(name=None, prefix=None, pkgs=None):

    """
    Create an environment either by name or path with a specified set of
    packages
    """

    if not pkgs or not isinstance(pkgs, (list, tuple)):
        raise TypeError('must specify a list of one or more packages to '
                        'install into new environment')

    cmd_list = ['create', '--yes', '--quiet']
    if name:
        ref         = name
        search      = [os.path.join(d, name) for d in info()['envs_dirs']]
        cmd_list    = ['create', '--yes', '--quiet', '--name', name]
    elif prefix:
        ref         = prefix
        search      = [prefix]
        cmd_list    = ['create', '--yes', '--quiet', '--prefix', prefix]
    else:
        raise TypeError('must specify either an environment name or a path '
                        'for new environment')

    if any(os.path.exists(prefix) for prefix in search):
        raise CondaEnvExistsError('Conda environment [%s] already exists' % ref)

    cmd_list.extend(pkgs)
    (out, err) = _call_conda(cmd_list)
    if err.decode().strip():
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), err.decode()))
    return out

# -----------------------------------------------------------------

def install(name=None, prefix=None, pkgs=None):

    """
    Install packages into an environment either by name or path with a
    specified set of packages
    :param name:
    :param prefix:
    :param pkgs:
    """

    if not pkgs or not isinstance(pkgs, (list, tuple)):
        raise TypeError('must specify a list of one or more packages to '
                        'install into existing environment')

    cmd_list = ['install', '--yes', '--quiet']
    if name:
        cmd_list.extend(['--name', name])
    elif prefix:
        cmd_list.extend(['--prefix', prefix])
    else: # just install into the current environment, whatever that is
        pass

    cmd_list.extend(pkgs)
    (out, err) = _call_conda(cmd_list)
    if err.decode().strip():
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), err.decode()))
    return out

# -----------------------------------------------------------------

def update(*pkgs, **kwargs):

    """
    Update package(s) (in an environment) by name.
    """

    cmd_list = ['update', '--json', '--quiet', '--yes']

    if not pkgs and not kwargs.get('all'):
        raise TypeError("Must specify at least one package to update, or all=True.")

    if kwargs.get('name') and kwargs.get('path'):
        raise TypeError('conda remove: At most one of name, path allowed')

    if kwargs.get('name'):
        cmd_list.extend(['--name', kwargs.pop('name')])

    if kwargs.get('path'):
        cmd_list.extend(['--prefix', kwargs.pop('path')])

    cmd_list.extend(
        _setup_install_commands_from_kwargs(
            kwargs,
            ('dry_run', 'no_deps', 'override_channels',
             'no_pin', 'force', 'all', 'use_index_cache', 'use_local',
             'alt_hint')))

    cmd_list.extend(pkgs)

    result = _call_and_parse(cmd_list, conda_path=kwargs.pop("conda_path", "conda"))

    if 'error' in result:
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), result['error']))

    return result

# -----------------------------------------------------------------

def remove(*pkgs, **kwargs):

    """
    Remove a package (from an environment) by name.

    Returns {
        success: bool, (this is always true),
        (other information)
    }
    """

    cmd_list = ['remove', '--json', '--quiet', '--yes']

    if not pkgs and not kwargs.get('all'):
        raise TypeError("Must specify at least one package to remove, or all=True.")

    if kwargs.get('name') and kwargs.get('path'):
        raise TypeError('conda remove: At most one of name, path allowed')

    if kwargs.get('name'):
        cmd_list.extend(['--name', kwargs.pop('name')])

    if kwargs.get('path'):
        cmd_list.extend(['--prefix', kwargs.pop('path')])

    cmd_list.extend(
        _setup_install_commands_from_kwargs(
            kwargs,
            ('dry_run', 'features', 'override_channels',
             'no_pin', 'force', 'all')))

    cmd_list.extend(pkgs)

    result = _call_and_parse(cmd_list, conda_path=kwargs.pop("conda_path", "conda"))

    if 'error' in result:
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), result['error']))

    return result

# -----------------------------------------------------------------

def remove_environment(name=None, path=None, **kwargs):

    """
    Remove an environment entirely.
    See ``remove``.
    """

    return remove(name=name, path=path, all=True, **kwargs)

# -----------------------------------------------------------------

def clone_environment(clone, name=None, path=None, **kwargs):

    """
    Clone the environment ``clone`` into ``name`` or ``path``.
    """

    cmd_list = ['create', '--json', '--quiet']

    if (name and path) or not (name or path):
        raise TypeError("conda clone_environment: exactly one of name or path required")

    if name:
        cmd_list.extend(['--name', name])

    if path:
        cmd_list.extend(['--prefix', path])

    cmd_list.extend(['--clone', clone])

    cmd_list.extend(
        _setup_install_commands_from_kwargs(
            kwargs,
            ('dry_run', 'unknown', 'use_index_cache', 'use_local', 'no_pin',
             'force', 'all', 'channel', 'override_channels', 'no_default_packages')))

    result = _call_and_parse(cmd_list, conda_path=kwargs.pop("conda_path", "conda"))

    if 'error' in result:
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), result['error']))

    return result

# -----------------------------------------------------------------

def process(name=None, prefix=None, cmd=None, args=None,
            stdin=None, stdout=None, stderr=None, timeout=None):

    """
    Create a Popen process for cmd using the specified args but in the conda
    environment specified by name or prefix.

    The returned object will need to be invoked with p.communicate() or similar.
    """

    if bool(name) == bool(prefix):
        raise TypeError('exactly one of name or prefix must be specified')

    if not cmd:
        raise TypeError('cmd to execute must be specified')

    if not args:
        args = []

    if name:
        prefix = get_prefix_envname(name)

    conda_env = dict(os.environ)

    if sys.platform == 'win32':
        conda_env['PATH'] = join(prefix, 'Scripts') + os.pathsep + conda_env['PATH']
    else: # Unix
        conda_env['PATH'] = join(prefix, 'bin') + os.pathsep + conda_env['PATH']

    conda_env['PATH'] = prefix + os.pathsep + conda_env['PATH']

    cmd_list = [cmd]
    cmd_list.extend(args)

    try:
        p = Popen(cmd_list, env=conda_env, stdin=stdin, stdout=stdout, stderr=stderr)
    except OSError:
        raise Exception("could not invoke %r\n" % cmd_list)
    return p

# -----------------------------------------------------------------

def _setup_config_from_kwargs(kwargs):

    """
    This function ...
    :param kwargs:
    :return:
    """

    cmd_list = ['--json', '--force']

    if 'file' in kwargs:
        cmd_list.extend(['--file', kwargs['file']])

    if 'system' in kwargs:
        cmd_list.append('--system')

    return cmd_list

# -----------------------------------------------------------------

def config_path(**kwargs):

    """
    Get the path to the config file.
    """

    cmd_list = ['config', '--get']
    cmd_list.extend(_setup_config_from_kwargs(kwargs))

    result = _call_and_parse(cmd_list, conda_path=kwargs.pop("conda_path", "conda"))

    if 'error' in result:
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), result['error']))
    return result['rc_path']

# -----------------------------------------------------------------

def config_get(*keys, **kwargs):

    """
    Get the values of configuration keys.

    Returns a dictionary of values. Note, the key may not be in the
    dictionary if the key wasn't set in the configuration file.
    """

    cmd_list = ['config', '--get']
    cmd_list.extend(keys)
    cmd_list.extend(_setup_config_from_kwargs(kwargs))

    result = _call_and_parse(cmd_list, conda_path=kwargs.pop("conda_path", "conda"))

    if 'error' in result:
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), result['error']))
    return result['get']

# -----------------------------------------------------------------

def config_set(key, value, **kwargs):

    """
    Set a key to a (bool) value.
    Returns a list of warnings Conda may have emitted.
    """

    cmd_list = ['config', '--set', key, str(value)]
    cmd_list.extend(_setup_config_from_kwargs(kwargs))

    result = _call_and_parse(cmd_list, conda_path=kwargs.pop("conda_path", "conda"))

    if 'error' in result:
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), result['error']))
    return result.get('warnings', [])

# -----------------------------------------------------------------

def config_add(key, value, **kwargs):

    """
    Add a value to a key.
    Returns a list of warnings Conda may have emitted.
    """

    cmd_list = ['config', '--add', key, value]
    cmd_list.extend(_setup_config_from_kwargs(kwargs))

    result = _call_and_parse(cmd_list, conda_path=kwargs.pop("conda_path", "conda"))

    if 'error' in result:
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), result['error']))
    return result.get('warnings', [])

# -----------------------------------------------------------------

def config_remove(key, value, **kwargs):

    """
    Remove a value from a key.
    Returns a list of warnings Conda may have emitted.
    """

    cmd_list = ['config', '--remove', key, value]
    cmd_list.extend(_setup_config_from_kwargs(kwargs))

    result = _call_and_parse(cmd_list, conda_path=kwargs.pop("conda_path", "conda"))

    if 'error' in result:
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), result['error']))
    return result.get('warnings', [])

# -----------------------------------------------------------------

def config_delete(key, **kwargs):

    """
    Remove a key entirely.
    Returns a list of warnings Conda may have emitted.
    """

    cmd_list = ['config', '--remove-key', key]
    cmd_list.extend(_setup_config_from_kwargs(kwargs))

    result = _call_and_parse(cmd_list, conda_path=kwargs.pop("conda_path", "conda"))

    if 'error' in result:
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), result['error']))
    return result.get('warnings', [])

# -----------------------------------------------------------------

def run(command, conda_path="conda"):

    """
    Launch the specified app by name or full package name.
    Returns a dictionary containing the key "fn", whose value is the full
    package (ending in ``.tar.bz2``) of the app.
    """

    cmd_list = ['run', '--json', command]

    result = _call_and_parse(cmd_list, conda_path=conda_path)

    if 'error' in result:
        raise CondaError('conda %s: %s' % (" ".join(cmd_list), result['error']))
    return result

# -----------------------------------------------------------------

def test():

    """
    Self-test function, which prints useful debug information.
    This function returns None on success, and will crash the interpreter
    on failure.
    """

    print('sys.version: %r' % sys.version)
    print('sys.prefix : %r' % sys.prefix)
    conda_version = get_conda_version()
    print('conda version: %r' % conda_version)
    print('conda info:')
    d = info()
    for kv in d.items():
        print('\t%s=%r' % kv)
    assert d['conda_version'] == conda_version
    print('OK')

# -----------------------------------------------------------------
