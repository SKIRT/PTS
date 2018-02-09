#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.handle Contains the ExecutionHandle class.

# -----------------------------------------------------------------

local = "local"
postponed = "postponed"
tty = "tty"
screen = "screen"
job = "job"
group_job = "group-job"
sql = "sql"

# -----------------------------------------------------------------

class ExecutionHandle(object):

    """
    This class ...
    """

    def __init__(self, type, value=None, host_id=None, remote_screen_output_path=None, remote_screen_script_path=None):

        """
        The constructor ...
        :param type:
        :param value:
        :param host_id:
        :param remote_screen_output_path:
        :param remote_screen_script_path:
        """

        # Type (job, screen or tty) and value (job ID, screen name or session number)
        self.type = type
        self.value = value

        # Remote host
        self.host_id = host_id

        # Extra information
        self.remote_screen_output_path = remote_screen_output_path
        self.remote_screen_script_path = remote_screen_script_path

    # -----------------------------------------------------------------

    @classmethod
    def local(cls):

        """
        This function ...
        :return:
        """

        return cls(local)

    # -----------------------------------------------------------------

    @classmethod
    def postponed(cls, host_id):

        """
        This function ...
        :return: 
        """

        return cls(postponed, host_id=host_id)

    # -----------------------------------------------------------------

    @classmethod
    def tty(cls, session_number, host_id):

        """
        This function ...
        :param session_number:
        :param host_id:
        :return:
        """

        return cls(tty, session_number, host_id)

    # -----------------------------------------------------------------

    @classmethod
    def screen(cls, screen_name, host_id, remote_screen_output_path=None, remote_screen_script_path=None):

        """
        This function ...
        :param screen_name:
        :param host_id:
        :param remote_screen_output_path:
        :param remote_screen_script_path:
        :return:
        """

        return cls(screen, screen_name, host_id, remote_screen_output_path=remote_screen_output_path, remote_screen_script_path=remote_screen_script_path)

    # -----------------------------------------------------------------

    @classmethod
    def job(cls, job_id, host_id):

        """
        This function ...
        :param job_id:
        :param host_id:
        :return:
        """

        return cls(job, job_id, host_id)

    # -----------------------------------------------------------------

    @classmethod
    def group_job(cls, job_id, host_id):

        """
        This function ...
        :param job_id:
        :param host_id:
        :return:
        """

        return cls(group_job, job_id, host_id)

    # -----------------------------------------------------------------

    @classmethod
    def sql(cls, name, host_id):

        """
        This function ...
        :param name:
        :param host_id:
        :return:
        """

        return cls(sql, name, host_id)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        string = self.type
        if self.value is not None: string += ": " + str(self.value)
        if self.host_id is not None: string += " [" + self.host_id + "]"
        if self.remote_screen_output_path is not None: string += " (screen output in '" + self.remote_screen_output_path + "')"
        return string

    # -----------------------------------------------------------------

    def to_dict(self):

        """
        This function ...
        :return:
        """

        items = dict()
        items["type"] = self.type

        # Value
        if self.type is not None:

            if self.type == tty: items["session"] = str(self.value)
            elif self.type == screen: items["screen name"] = self.value
            elif self.type == job: items["job ID"] = str(self.value)
            elif self.type == group_job: items["job ID"] = str(self.value)
            elif self.type == sql: items["name"] = str(self.value)
            else: raise ValueError("Type not recognized: '" + self.type + "'")

        # Host iD
        if self.host_id is not None: items["host ID"] = self.host_id

        # Screen output path
        if self.remote_screen_output_path is not None: items["remote screen output path"] = self.remote_screen_output_path

        # Return the dictionary
        return items

    # -----------------------------------------------------------------

    def to_lines(self, line_prefix="", fancy=False):

        """
        This function ...
        :param line_prefix:
        :param fancy:
        :return:
        """

        from ..tools import formatting as fmt

        items = self.to_dict()

        lines = []

        for key in items:

            if fancy: line = line_prefix + fmt.bold + key + fmt.reset + ": " + items[key]
            else: line = line_prefix + key + ": " + items[key]
            lines.append(line)

        return lines

# -----------------------------------------------------------------
