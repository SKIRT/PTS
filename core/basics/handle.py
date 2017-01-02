#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.handle Contains the ExecutionHandle class.

# -----------------------------------------------------------------

class ExecutionHandle(object):

    """
    This class ...
    """

    def __init__(self, type, value, host_id, remote_screen_output_path=None):

        """
        The constructor ...
        :param type:
        :param value:
        :param host_id:
        :param remote_screen_output_path:
        """

        # Type (job, screen or tty) and value (job ID, screen name or session number)
        self.type = type
        self.value = value

        # Remote host
        self.host_id = host_id

        # Extra information
        self.remote_screen_output_path = remote_screen_output_path

    # -----------------------------------------------------------------

    @classmethod
    def tty(cls, session_number, host_id):

        """
        This function ...
        :param session_number:
        :param host_id:
        :return:
        """

        return cls("tty", session_number, host_id)

    # -----------------------------------------------------------------

    @classmethod
    def screen(cls, screen_name, host_id, remote_screen_output_path=None):

        """
        This function ...
        :param screen_name:
        :param host_id:
        :param remote_screen_output_path:
        :return:
        """

        return cls("screen", screen_name, host_id, remote_screen_output_path)

    # -----------------------------------------------------------------

    @classmethod
    def job(cls, job_id, host_id):

        """
        This function ...
        :param job_id:
        :param host_id:
        :return:
        """

        return cls("job", job_id, host_id)

    # -----------------------------------------------------------------

    @classmethod
    def group_job(cls, job_id, host_id):

        """
        This function ...
        :param job_id:
        :param host_id:
        :return:
        """

        return cls("group-job", job_id, host_id)

    # -----------------------------------------------------------------

    @classmethod
    def sql(cls, name, host_id):

        """
        This function ...
        :param name:
        :param host_id:
        :return:
        """

        return cls("sql", name, host_id)

# -----------------------------------------------------------------
