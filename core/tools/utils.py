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
import os
import copy
from collections import OrderedDict, Set
from functools import wraps, partial

# -----------------------------------------------------------------

class LazyDictionary(OrderedDict):

    """
    This class ...
    """

    def __init__(self, evaluator, **kwargs):

        """
        Thisf unction ...
        :param evaluator:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(LazyDictionary, self).__init__()

        # The evaluator function
        self.evaluator = evaluator

        # The general evaluator kwargs
        self.general_kwargs = kwargs

        # Keyword arguments for different elements
        self.kwargs = dict()

        # Names of elements that are already evaluated
        self.evaluated = dict()

    # -----------------------------------------------------------------

    def get_kwargs(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Specific kwargs
        if name in self.kwargs:

            # Merge general with unique kwargs
            kwargs = self.general_kwargs.copy()
            for key in self.kwargs[name]: kwargs[key] = self.kwargs[name][key]
            return kwargs

        # Only the general kwargs
        else: return self.general_kwargs

    # -----------------------------------------------------------------

    def set_kwargs(self, name, **kwargs):

        """
        This function ...
        :param name:
        :param kwargs:
        :return:
        """

        if name not in self.kwargs: self.kwargs[name] = dict()
        self.kwargs[name].update(**kwargs)

    # -----------------------------------------------------------------

    def get_raw(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return super(LazyDictionary, self).__getitem__(name)

    # -----------------------------------------------------------------

    def __getitem__(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Evaluated?
        if self.evaluated[name]: return self.get_raw(name)

        # Not yet evaluated
        else:

            # Get value
            value = self.get_raw(name)
            kwargs = self.get_kwargs(name)
            value = self.evaluator(value, **kwargs)

            # Set evaluated value
            self[name] = value

            # Set evaluated flag
            self.evaluated[name] = True

            # Return the value
            return value

    # -----------------------------------------------------------------

    def set(self, name, value, **kwargs):

        """
        Thisj function ...
        :param name:
        :param value:
        :param kwargs:
        :return:
        """

        # Set
        self[name] = value

        # Set kwargs
        if len(kwargs) > 0: self.kwargs[name] = kwargs

    # -----------------------------------------------------------------

    def set_value(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        # Base class implementation
        super(LazyDictionary, self).__setitem__(name, value)

        # Set evaluated flag to TRUE
        self.evaluated[name] = True

    # -----------------------------------------------------------------

    def __setitem__(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        # Base class implementation
        super(LazyDictionary, self).__setitem__(name, value)

        # Set evaluated flag
        self.evaluated[name] = False

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

    def __init__(self, fget, fset=None, fdel=None, doc=None, key=None):

        """
        The constructor ...
        :param fget:
        :param fset:
        :param fdel:
        :param doc:
        :param key:
        """

        # Call the constructor of the base class
        super(lazyproperty, self).__init__(fget, fset, fdel, doc)

        if key is not None: self._key = key
        else: self._key = self.fget.__name__

        #print("KEY:", self._key)

    # -----------------------------------------------------------------

    def __get__(self, obj, owner=None):

        # print("TYPE", type(obj))
        # print("DICT", obj.__dict__)
        # print("KEY", self._key)
        # obj is the instance of the class of which the method is decorated
        # self is the decorator class instance

        try: return obj.__dict__[self._key]
        except KeyError:

            val = self.fget(obj)
            obj.__dict__[self._key] = val
            return val

        except AttributeError:
            if obj is None: return self
            raise

    # -----------------------------------------------------------------

    def __set__(self, obj, val):

        # obj is the instance of the class of which the method is decorated
        # self is the decorator class instance

        obj_dict = obj.__dict__
        if self.fset:
            ret = self.fset(obj, val)
            if ret is not None and obj_dict.get(self._key) is ret:
                # By returning the value set the setter signals that it took
                # over setting the value in obj.__dict__; this mechanism allows
                # it to override the input value
                return
        obj_dict[self._key] = val

    # -----------------------------------------------------------------

    def __delete__(self, obj):

        # obj is the instance of the class of which the method is decorated
        # self is the decorator class instance

        if self.fdel: self.fdel(obj)
        if self._key in obj.__dict__:
            del obj.__dict__[self._key]

    # -----------------------------------------------------------------

    # SOMETHING LIKE THIS CAN NEVER WORK: CALING A METHOD ON THE PROPERTY IS IMPOSSIBLE BECAUSE
    # OBJ.PROPERTY_NAME GIVES THE RETURN VALUE, SO .RESET ON THIS VALUE DOESN'T EXIST
    # USE THE KEYWORD 'DEL' (CALLING __DELETE__) ON THE PROPERTY TO RESET
    # def reset(self, *args, **kwargs):
    #     print(args)
    #     print(kwargs)

# -----------------------------------------------------------------

class lazyfileproperty(object):

    """
    This class ...
    Example:
      >>>  class A(object):
      ...     pass
      ...  filepath = "test.dat"
      >>>  @lazyfileproperty(A, filepath)
      ...  def my_calc(arg=3):
      ...     print("Executing my_calc(%s)" % arg)
      ...     return 1
    """

    def __init__(self, cls, path, is_attribute=False, write=True, fsave=None, fload=None):

        """
        This function ...
        :param cls:
        :param path:
        :param is_attribute:
        :param write:
        :param fsave:
        :param fload:
        """

        # Set things
        self.cls = cls
        self.path = path
        self.is_attribute = is_attribute # is the path actually the name of an attribute of the object?
        self.write = write

        # Set save and load functions
        self.fsave = fsave
        self.fload = fload

    # -----------------------------------------------------------------

    @property
    def has_fsave(self):
        return self.fsave is not None

    # -----------------------------------------------------------------

    @property
    def has_fload(self):
        return self.fload is not None

    # -----------------------------------------------------------------

    def __call__(self, func):

        """
        This function ...
        :param func:
        :return:
        """

        #print("Wrapping function '%s'" % func.__name__)

        # Define function to be decorated
        def wrapper(*args, **kwargs):

            #print("Determining path ...")

            # Determine path
            if self.is_attribute: filepath = getattr(args[0], self.path)
            else: filepath = self.path

            #print("Checking path '" + filepath + "' ...")

            # Check whether file exist
            # If so, load it
            if os.path.isfile(filepath):

                if self.has_fload: out = self.fload(filepath)
                else: out = self.cls.from_file(filepath)

            # File does not exist yet
            else:

                # Calculate
                out = func(*args, **kwargs)

                # Save to path
                if self.write:
                    print("Writing to path '" + filepath + "' ...")
                    if self.has_fsave: self.fsave(filepath, out)
                    else: out.saveto(filepath)

            # Return
            return out

        # Return the actual decorator function
        return lazyproperty(wrapper, doc=func.__doc__, key=func.__name__)

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

def create_static_class(cls, name=None):

    """
    This function ...
    :param cls:
    :param name:
    :return:
    """

    # TODO

    # Lazified + cache function arguments

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

class ForbiddenOperation(Exception):

    """
    This exception ...
    """

    def __init__(self, cls, operation):

        """
        The constructor ...
        :param cls:
        :param operation:
        """

        # Set the class and operation name
        self.cls = cls
        self.operation = operation

        # Set the message
        message = "The " + operation + " operation of class '" + str(type(cls)) + "' is forbidden"
        self.message = message

        # Call the constructor of the base class
        super(ForbiddenOperation, self).__init__(message)

# -----------------------------------------------------------------

class DefaultScope(object):

    """
    This class ...
    """

    #def __init__(self):

    def __enter__(self):
        return None

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is not None: return False # not succesful
        else: pass # succesful
        return None

# -----------------------------------------------------------------

def memoize(function):

    memo = {}

    @wraps(function)
    def wrapper(*args):
        if args in memo:
            return memo[args]
        else:
            rv = function(*args)
            memo[args] = rv
            return rv

    return wrapper

# -----------------------------------------------------------------

# EXAMPLE:
@memoize
def fibonacci(n):
    if n < 2: return n
    return fibonacci(n - 1) + fibonacci(n - 2)

# -----------------------------------------------------------------

class memoize_method(object):

    """
    THIS ONE STORES THE CACHE IN THE CLASS OF WHICH THE METHOD IS DECORATED!

    Cache the return value of a method

    This class is meant to be used as a decorator of methods. The return value
    from a given method invocation will be cached on the instance whose method
    was invoked. All arguments passed to a method decorated with memoize must
    be hashable.

    If a memoized method is invoked directly on its class the result will not
    be cached. Instead the method will be invoked like a static method:
    class Obj(object):

        @memoize
        def add_to(self, arg):
            return self + arg

    Obj.add_to(1) # not enough arguments
    Obj.add_to(1, 2) # returns 3, result is not cached
    """

    def __init__(self, func):

        self.func = func

    def __get__(self, obj, objtype=None):

        if obj is None: return self.func
        return partial(self, obj)

    def __call__(self, *args, **kw):

        obj = args[0]

        try: cache = obj.__cache
        except AttributeError: cache = obj.__cache = {}

        key = (self.func, args[1:], frozenset(kw.items()))

        try: res = cache[key]
        except KeyError: res = cache[key] = self.func(*args, **kw)

        return res

# -----------------------------------------------------------------

class memoize_method_reset(object):

   """
   THIS ONE STORES THE CACHE INSIDE THE DECORATOR CLASS!
   THIS ONE HAS A RESET METHOD

   Decorator that caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned, and
   not re-evaluated.
   """

   def __init__(self, func):

      self.func = func
      self.cache = {}

   def __call__(self, *args):

      try: return self.cache[args]
      except KeyError:
         value = self.func(*args)
         self.cache[args] = value
         return value
      except TypeError:
         # uncachable -- for instance, passing a list as an argument.
         # Better to not cache than to blow up entirely.
         return self.func(*args)

   def __repr__(self):

      """Return the function's docstring."""

      return self.func.__doc__

   def __get__(self, obj, objtype):

      """Support instance methods."""

      fn = partial(self.__call__, obj)
      fn.reset = self._reset
      return fn

   def _reset_for_args(self, *args):
       if args not in self.cache: return
       del self.cache[args]

   def _reset(self):
      self.cache = {}

# -----------------------------------------------------------------

# # example usage
# class Test(object):
#     v = 0
#
#     @memoize
#     def inc_add(self, arg):
#         self.v += 1
#         return self.v + arg
#
#
# t = Test()
# assert t.inc_add(2) == t.inc_add(2)
# assert Test.inc_add(t, 2) != Test.inc_add(t, 2)

# -----------------------------------------------------------------

# Source: https://programmingideaswithjake.wordpress.com/2015/05/23/python-decorator-for-simplifying-delegate-pattern/
class DelegatedAttribute:

    def __init__(self, delegate_name, attr_name):
        self.attr_name = attr_name
        self.delegate_name = delegate_name

    def __get__(self, instance, owner):
        if instance is None:
            return self
        else:
            # return instance.delegate.attr
            return getattr(self.delegate(instance),  self.attr_name)

    def __set__(self, instance, value):
        # instance.delegate.attr = value
        setattr(self.delegate(instance), self.attr_name, value)

    def __delete__(self, instance):
        delattr(self.delegate(instance), self.attr_name)

    def delegate(self, instance):
        return getattr(instance, self.delegate_name)

    def __str__(self):
        return ""

# -----------------------------------------------------------------

# EXAMPLE:

# class Foo:
#     def __init__(self):
#         self.bar = 'bar in foo'
#     def baz(self):
#         return 'baz in foo'
#
# class Baz:
#     def __init__(self):
#         self.foo = Foo()
#     bar = DelegatedAttribute('foo', 'bar')
#     baz = DelegatedAttribute('foo', 'baz')
#
# x = Baz()
# print(x.bar)  # prints 'bar in foo'
# print(x.baz())  # prints 'baz in foo'

# -----------------------------------------------------------------

# You may have noticed the creation of something called a SimpleProperty, which won’t be shown.
# It is simply a descriptor that holds and returns a value given to it, the most basic of properties.
class SimpleProperty:

    def __init__(self, initial_value=None):
        self.value = initial_value

    def __get__(self, obj, objtype):
        return self.value

    def __set__(self, obj, value):
        self.value = value

    def __delete__(self, obj):
        self.value = None

# -----------------------------------------------------------------

def delegate_as(delegate_cls, to='delegate', include=frozenset(), ignore=frozenset()):

    """
    This decorator ...
    :param delegate_cls:
    :param to:
    :param include:
    :param ignore:
    :return:
    """

    # turn include and ignore into sets, if they aren't already
    if not isinstance(include, Set):
        include = set(include)
    if not isinstance(ignore, Set):
        ignore = set(ignore)
    delegate_attrs = set(delegate_cls.__dict__.keys())
    attributes = include | delegate_attrs - ignore

    def inner(cls):

        # create property for storing the delegate
        setattr(cls, to, SimpleProperty())
        # don't bother adding attributes that the class already has
        attrs = attributes - set(cls.__dict__.keys())
        # set all the attributes
        for attr in attrs:
            setattr(cls, attr, DelegatedAttribute(to, attr))
        return cls

    return inner

# -----------------------------------------------------------------

def defaultlen(value, default):

    """
    This function ...
    :param value:
    :param default:
    :return:
    """

    try: return len(value)
    except TypeError: return default

# -----------------------------------------------------------------
