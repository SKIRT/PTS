#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import image modules
from image import Image

# Import standard modules
from cmd import Cmd
from inspect import isfunction

# *****************************************************************

class MyPrompt(Cmd):

    members = Image.__dict__

    for member_name in members:

        pass

    def do_hello(self, args):
        """Says hello. If you provide a name, it will greet you with it."""
        if len(args) == 0:
            name = 'stranger'
        else:
            name = args
        print "Hello, %s" % name

    def do_quit(self, args):
        """Quits the program."""
        print "Quitting."
        raise SystemExit

if __name__ == '__main__':
    
    prompt = MyPrompt()
    prompt.prompt = '> '
    prompt.cmdloop('Starting prompt...')