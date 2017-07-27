#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.utils Contains utilities.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy

# -----------------------------------------------------------------

class lazyproperty(property):

    """
    Works similarly to property(), but computes the value only once.

    This essentially memorizes the value of the property by storing the result
    of its computation in the ``__dict__`` of the object instance.  This is
    useful for computing the value of some property that should otherwise be
    invariant.  For example::

        >>> class LazyTest(object):
        ...     @lazyproperty
        ...     def complicated_property(self):
        ...         print('Computing the value for complicated_property...')
        ...         return 42
        ...
        >>> lt = LazyTest()
        >>> lt.complicated_property
        Computing the value for complicated_property...
        42
        >>> lt.complicated_property
        42

    As the example shows, the second time ``complicated_property`` is accessed,
    the ``print`` statement is not executed.  Only the return value from the
    first access off ``complicated_property`` is returned.

    By default, a setter and deleter are used which simply overwrite and
    delete, respectively, the value stored in ``__dict__``. Any user-specified
    setter or deleter is executed before executing these default actions.
    The one exception is that the default setter is not run if the user setter
    already sets the new value in ``__dict__`` and returns that value and the
    returned value is not ``None``.

    Adapted from the recipe at
    http://code.activestate.com/recipes/363602-lazy-property-evaluation
    """

    def __init__(self, fget, fset=None, fdel=None, doc=None):
        super(lazyproperty, self).__init__(fget, fset, fdel, doc)
        self._key = self.fget.__name__

    def __get__(self, obj, owner=None):
        try:
            #print("TYPE", type(obj))
            #print("DICT", obj.__dict__)
            #print("KEY", self._key)
            return obj.__dict__[self._key]
        except KeyError:
            val = self.fget(obj)
            obj.__dict__[self._key] = val
            return val
        except AttributeError:
            if obj is None:
                return self
            raise

    def __set__(self, obj, val):
        obj_dict = obj.__dict__
        if self.fset:
            ret = self.fset(obj, val)
            if ret is not None and obj_dict.get(self._key) is ret:
                # By returning the value set the setter signals that it took
                # over setting the value in obj.__dict__; this mechanism allows
                # it to override the input value
                return
        obj_dict[self._key] = val

    def __delete__(self, obj):
        if self.fdel:
            self.fdel(obj)
        if self._key in obj.__dict__:
            del obj.__dict__[self._key]

# -----------------------------------------------------------------

def decorate_all_methods(decorator):

    """
    This class decorator decorates all methods

    # USAGE:
    @decorate_all_methods(decorator)
    class C(object):
        def m1(self): pass
        def m2(self, x): pass
        ...

    :param decorator:
    :return:
    """

    def decorate(cls):
        for attr in cls.__dict__: # there's probably a better way to do this
            if callable(getattr(cls, attr)):
                setattr(cls, attr, decorator(getattr(cls, attr)))
        return cls
    return decorate

# -----------------------------------------------------------------

def decorate_all_properties(decorator):

    """
    This class decorator decorates all properties

    # USAGE:
    @decorate_all_properties(decorator)
    class C(object):
        @property
        def m1(self): pass
        @property
        def m2(self): pass
        ...

    :param decorator:
    :return:
    """

    def decorate(cls):
        for attr in cls.__dict__:
            if isinstance(getattr(cls, attr), property):
                setattr(cls, attr, decorator(getattr(cls, attr)))
        return cls
    return decorate

# -----------------------------------------------------------------

def lazify():

    """
    This class decorator decorates all properties with the lazyproperty decorator

    # USAGE:
    @lazify
    class C(object):
        @property
        def m1(self): pass
        @property
        def m2(self): pass
    :return:
    """

    return decorate_all_properties(lazyproperty)

# -----------------------------------------------------------------

def is_hashable(thing):

    """
    Thisf unction ...
    :param thing:
    :return:
    """

    try: hash(thing)
    except TypeError: return False
    else: return True

# -----------------------------------------------------------------

def is_mutable(thing):

    """
    This function ...
    :param thing:
    :return:
    """

    # ACTUALLY, THIS IS NOT ENTIRELY TRUE!
    return not is_hashable(thing)

# -----------------------------------------------------------------

def copy_class(cls, name):

    """
    Thisf unction ...
    :param cls:
    :param name:
    :return:
    """

    return type(name, cls.__bases__, dict(cls.__dict__))

# -----------------------------------------------------------------

def create_lazified_class(cls, name=None):

    """
    This function ...
    :param cls:
    :param name:
    :return:
    """

    bases = cls.__bases__

    #print(bases)

    clsdict = {}

    # Loop over all attributes (methods, properties, ...) of the original class
    for attrname in cls.__dict__:

        # WEIRD, BUT __DICT__ IS ALSO A ATTRNAME IN __DICT__ ????????
        if attrname == "__dict__": continue
        #if attrname.startswith("__") and attrname.endswith("__") and attrname != "__init__": continue # NOT NECESSARY SO STRICT

        # Get the value for the original class
        value = cls.__dict__[attrname]

        # lazify, copy or just plug in
        if isinstance(value, property): clsdict[attrname] = lazyproperty(value.fget)
        elif is_mutable(value): clsdict[attrname] = copy.copy(value)  # FOR EXAMPLE, A CLASS COULD HAVE A CLASS PROPERTY THAT IS A LIST
        else: clsdict[attrname] = value

    # Determine name
    if name is None: name = "Lazy" + cls.__name__

    #print(clsdict)

    # Create the new class
    return type(name, bases, clsdict)

# -----------------------------------------------------------------

class lazyproperties(type):

    """
    This metaclass makes all properties lazy

    # USAGE:
    class myClass(object):
        __metaclass__ = lazyproperties
        def baz(self):
            print self.baz.foo
    """

    def __new__(cls, name, bases, local):

        for attr in local:

            value = local[attr]

            #if callable(value):
            #    local[attr] = myDecorator(value)
            if isinstance(value, property): local[attr] = lazyproperty(value)

        # Create the class
        return type.__new__(cls, name, bases, local)

# -----------------------------------------------------------------

class abstractclassmethod(classmethod):

    """
    This class is a decorator for abstract class methods
    """

    __isabstractmethod__ = True

    def __init__(self, callable):

        callable.__isabstractmethod__ = True
        super(abstractclassmethod, self).__init__(callable)

# -----------------------------------------------------------------

class UserIntervention(Exception):

    """
    This exception should be called when user intervention / inspection is required and therefore
    the execution should be halted or aborted.
    """

    def __init__(self, message, cls, function_name):

        """
        The constructor ...
        :param message:
        :param cls:
        :param function_name:
        """

        # Set the message
        self.message = message

        # Set the class and function name
        self.cls = cls
        self.function_name = function_name

        # Call the constructor of the base class
        super(UserIntervention, self).__init__(message)

# -----------------------------------------------------------------
