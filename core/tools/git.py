#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.git Provides functions iteracting with git.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from . import introspection
from . import terminal

# -----------------------------------------------------------------

class AuthenticationError(Exception):

    """
    This class ...
    """

    def __init__(self, message, url=None):

        # Call the base class constructor with the parameters it needs
        super(AuthenticationError, self).__init__(message)

        # The URL for which authentication failed
        self.url = url

# -----------------------------------------------------------------

def clone(url, path, show_output=False):

    """
    This function ...
    :param url:
    :param path:
    :param show_output:
    :return:
    """

    command = "git clone " + url + " " + path
    terminal.execute(command, show_output=show_output)

# -----------------------------------------------------------------

def checkout_new_branch(repo_path, branch_name, origin_branch_name, show_output=False):

    """
    Thisn function ...
    :param repo_path:
    :param branch_name:
    :param origin_branch_name:
    :param show_output:
    :return:
    """

    command = "git -C " + repo_path + " checkout -b " + branch_name + " " + origin_branch_name
    terminal.execute(command, show_output=show_output)

# -----------------------------------------------------------------

def check_repository(repo_path, show_output=False):

    """
    This function ...
    :param repo_path:
    :param show_output:
    :return:
    """

    command = "git fsck"
    terminal.execute(command, cwd=repo_path, show_output=show_output)

# -----------------------------------------------------------------

def get_hash_remote_repository(url):

    """
    This function ...
    :return:
    """

    # Import here instead at module level to accomodate clean python installs
    import pexpect

    # Check whether the repo is up-to-date
    command = "git ls-remote " + url + " HEAD"
    # child = pexpect.spawn(command, timeout=30, logfile=sys.stdout)
    child = pexpect.spawn(command, timeout=30)
    index = child.expect([pexpect.EOF, "':"])

    # A prompt appears where username is asked
    if index == 1:

        # Username for 'https://github.ugent.be
        host_url = child.before.split(" for '")[1]
        host = host_url.split("https://")[1]

        # Host info is present
        if introspection.has_account(host):

            username, password = introspection.get_account(host)

            child.sendline(username)
            child.expect("':")
            child.sendline(password)
            child.expect(pexpect.EOF)

            #print(child.before)

            output = child.before.split("\r\n")

        # Error
        else: raise ValueError("Account info is not present for '" + host + "'")

    # No username and password required
    elif index == 0:
        #child.expect(pexpect.EOF)
        #print(child.before)
        output = child.before.split("\r\n")

    # This shouldn't happen
    else: raise RuntimeError("Something went wrong")

    #print(output)

    # Get the latest git hash from the output
    #latest_git_hash = output[0].split("\t")[0]

    latest_git_hash = None
    for line in output:
        if "HEAD" in line:
            latest_git_hash = line.split("\t")[0]

    # Return
    return latest_git_hash

# -----------------------------------------------------------------

def get_url_repository(remote, repo_path, repo_name="origin"):

    """
    This function ...
    :param remote:
    :param repo_path:
    :param repo_name:
    :return:
    """

    # Get the url of the repo from which cloned
    args = ["git", "remote", "show", repo_name]
    show_command = " ".join(args)

    # Change cwd
    original_cwd = remote.change_cwd(repo_path)

    # Send the 'git remote show origin' command and look at what comes out
    #remote.ssh.logfile = sys.stdout
    remote.ssh.sendline(show_command)
    #index = remote.ssh.expect([remote.ssh.PROMPT, "':"])
    index = remote.ssh.expect([remote.ssh.PROMPT, "Username for"])

    # A prompt appears where username is asked
    if index == 1:

        remote.ssh.expect("':")

        #print(remote.ssh.before)

        # Username for 'https://github.ugent.be
        #host_url = remote.ssh.before.split(" for '")[1]
        host_url = remote.ssh.before.split("'")[1].strip()

        #host_url = remote.ssh.before.split("Username for '")[1].split("'")[0]

        #remote.ssh.expect("':")

        #print(remote.ssh.before)
        #print(host_url)

        host = host_url.split("https://")[1]

        if introspection.has_account(host):

            username, password = introspection.get_account(host)

            remote.ssh.sendline(username)
            remote.ssh.expect("':")
            remote.ssh.sendline(password)
            remote.ssh.prompt()

            output = remote.ssh.before.split("\r\n")[1:-1]

        else: raise ValueError("Account info is not present for '" + host + "'")

    # No username and password required
    elif index == 0: output = remote.ssh.before.split("\r\n")[1:]

    # This shouldn't happen
    else: raise RuntimeError("Something went wrong")

    # Reset logfile
    remote.ssh.logfile = None

    # Reset cwd
    remote.change_cwd(original_cwd)

    # Check for errors
    for line in output:
        if "Authentication failed" in line:
            # Try to determine repo URL
            url = line.split("Authentication failed for '")[1][:-1]
            raise AuthenticationError("Connection problem for '" + url + "'", url=url)

    url = None
    for line in output:
        if "Fetch URL" in line: url = line.split(": ")[1]

    # Check
    if url is None: raise RuntimeError("Something went wrong: could not determine the repository url")

    # Return the url
    return url

# -----------------------------------------------------------------

def get_git_hash(remote, repo_path):

    """
    This function ...
    :param remote:
    :param repo_path:
    :return:
    """

    # Check hash
    get_hash_command = ["git", "rev-parse", "HEAD"]
    get_hash_command = " ".join(get_hash_command)
    git_hash = remote.execute(get_hash_command, cwd=repo_path)[0]

    return git_hash

# -----------------------------------------------------------------

def get_short_git_version(repo_path, remote=None):

    """
    This function ...
    :param repo_path:
    :param remote:
    :return:
    """

    if remote is not None:

        # Get the git version
        #first_part_command = "git rev-list --count HEAD" doesn't work for versions git 1.7.x !!
        second_part_command = "git describe --dirty --always"
        #first_part = remote.execute(first_part_command, cwd=repo_path)[0].strip()
        second_part = remote.execute(second_part_command, cwd=repo_path)[0].strip()

        first_part_command = "git rev-list HEAD"
        output = remote.execute(first_part_command, cwd=repo_path)

        first_part = str(len(output))

    else:

        # Get the git version
        #first_part_command = "git rev-list --count HEAD"
        #second_part_command = "git describe --dirty --always"
        #first_part = subprocess.check_output(first_part_command, cwd=repo_path).split("\n")[0].strip()
        #second_part = subprocess.check_output(second_part_command, cwd=repo_path).split("\n")[0].strip()

        first_part_command = "git rev-list --count HEAD"
        second_part_command = "git describe --dirty --always"
        first_part = terminal.execute_no_pexpect(first_part_command, cwd=repo_path)[0].strip()
        second_part = terminal.execute_no_pexpect(second_part_command, cwd=repo_path)[0].strip()

    git_version = first_part + "-" + second_part

    return git_version

# -----------------------------------------------------------------

def decompose_repo_url(url, return_type=False):

    """
    This function ...
    :param url:
    :param return_type:
    :return:
    """

    # Already https link
    if url.startswith("https://"):
        if return_type:
            host, user_or_organization, repo_name, username, password = decompose_https(url)
            return host, user_or_organization, repo_name, username, password, "https"
        else: return decompose_https(url)
    else:
        if return_type:
            host, user_or_organization, repo_name, username, password = decompose_ssh(url)
            return host, user_or_organization, repo_name, username, password, "ssh"
        else: return decompose_ssh(url)

# -----------------------------------------------------------------

def decompose_https(url):

    """
    This function ...
    :param url:
    :return:
    """

    if "@" in url:

        username = url.split("//")[1].split(":")[0]
        password = url.split(username + ":")[1].split("@")[0]
        host = url.split("@")[1].split("/")[0]
        user_or_organization = url.split(host + "/")[1].split("/")[0]
        repo_name = url.split(user_or_organization + "/")[1].split(".git")[0]

    else:

        host = url.split("//")[1].split("/")[0]
        user_or_organization = url.split(host + "/")[1].split("/")[0]
        repo_name = url.split(user_or_organization + "/")[1].split(".git")[0]
        username = None
        password = None

    # Return
    return host, user_or_organization, repo_name, username, password

# -----------------------------------------------------------------

def decompose_ssh(url):

    """
    This function ...
    :param url:
    :return:
    """

    # CONVERT TO HTTPS LINK
    host = url.split("@")[1].split(":")[0]
    user_or_organization = url.split(":")[1].split("/")[0]
    repo_name = url.split("/")[-1].split(".git")[0]

    username = None
    password = None

    # Return
    return host, user_or_organization, repo_name, username, password

# -----------------------------------------------------------------

def compose_repo_url(url_type, host, user_or_organization, repo_name, username=None, password=None):

    """
    This function ...
    :param url_type:
    :param host:
    :param user_or_organization:
    :param repo_name:
    :param username:
    :param password:
    :return:
    """

    # SSH
    if url_type == "ssh":

        # Check input
        if username is not None: raise ValueError("Username cannot be specified for SSH urls")
        if password is not None: raise ValueError("Password cannot be specified for SSH urls")

        # Compose
        return compose_ssh(host, user_or_organization, repo_name)

    # HTTPS
    elif url_type == "https": return compose_https(host, user_or_organization, repo_name, username=username, password=password)

    # Invalid
    else: raise ValueError("Url type '" + url_type + "' not recognized")

# -----------------------------------------------------------------

def compose_https(host, user_or_organization, repo_name, username=None, password=None):

    """
    This function ...
    :param host:
    :param user_or_organization:
    :param repo_name:
    :param username:
    :param password:
    :return:
    """

    # Username is provided
    if username is not None:

        # Add password if specified
        if password is not None: url = "https://" + username + ":" + password + "@" + host + "/" + user_or_organization + "/" + repo_name + ".git"
        else: url = "https://" + username + "@" + host + "/" + user_or_organization + "/" + repo_name + ".git"

    # No username
    else: url = "https://" + host + "/" + user_or_organization + "/" + repo_name + ".git"

    # Return the url
    return url

# -----------------------------------------------------------------

def compose_ssh(host, user_or_organization, repo_name):

    """
    This function ...
    :param host:
    :param user_or_organization:
    :param repo_name:
    :return:
    """

    # COmpose the URL
    url = "git@" + host + ":" + user_or_organization + "/" + repo_name + ".git"
    return url

# -----------------------------------------------------------------

def transform_to_simple_https(url):

    """
    This function ...
    :param url:
    :return:
    """

    host, user_or_organization, repo_name, username, password = decompose_repo_url(url)
    return compose_https(host, user_or_organization, repo_name)

# -----------------------------------------------------------------

def replace_remote_url(url, repo_path, repo_name="origin", remote=None, show_output=False):

    """
    This function ...
    :param url:
    :param repo_path:
    :param repo_name:
    :param remote:
    :param show_output:
    :return:
    """

    # Determine the command
    command = "git remote set-url " + repo_name + " '" + url + "'"

    # Remote
    if remote is not None:
        output = remote.execute(command, cwd=repo_path, show_output=show_output)
        for line in output:
            if "fatal:" in line: raise RuntimeError("Error: " + line)

    # Local
    else:
        output = terminal.execute_no_pexpect(command, cwd=repo_path, show_output=show_output)
        for line in output:
            if "fatal:" in line: raise RuntimeError("Error: " + line)

# -----------------------------------------------------------------
