#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.containers Contains container classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
from collections import OrderedDict, Callable

# Import the relevant PTS classes and modules
from ..tools import types
from ..filter.filter import parse_filter, Filter
from ..tools import filesystem as fs

# -----------------------------------------------------------------

# Source: http://stackoverflow.com/a/6190500/562769
class DefaultOrderedDict(OrderedDict):

    """
    This class ...
    """

    def __init__(self, default_factory=None, *a, **kw):

        """
        The cosntructor ...
        :param default_factory: 
        :param a: 
        :param kw: 
        """

        if (default_factory is not None and
           not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    # -----------------------------------------------------------------

    def __getitem__(self, key):

        """
        This function ...
        :param key: 
        :return: 
        """

        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    # -----------------------------------------------------------------

    def __missing__(self, key):

        """
        This fucntion ...
        :param key: 
        :return: 
        """

        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    # -----------------------------------------------------------------

    def __reduce__(self):

        """
        This function ...
        :return: 
        """

        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return: 
        """

        return self.__copy__()

    # -----------------------------------------------------------------

    def __copy__(self):

        """
        This function ...
        :return: 
        """

        return type(self)(self.default_factory, self)

    # -----------------------------------------------------------------

    def __deepcopy__(self, memo):

        """
        This function ...
        :param memo: 
        :return: 
        """

        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return: 
        """

        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory, OrderedDict.__repr__(self))

# -----------------------------------------------------------------

class KeyList(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        self.contents = OrderedDict()

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return: 
        """

        return len(self.contents)

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        THis function ...
        :return: 
        """

        for name in self.contents: yield self.contents[name]

    # -----------------------------------------------------------------

    def to_dictionary(self):

        """
        This function ...
        :return: 
        """

        return copy.copy(self.contents)

    # -----------------------------------------------------------------

    def index(self, value):

        """
        This function ...
        :param value: 
        :return: 
        """

        return self.values.index(value)

    # -----------------------------------------------------------------

    def key_for_value(self, value):

        """
        This function ...
        :param value: 
        :return: 
        """

        index = self.index(value)
        return self.keys[index]

    # -----------------------------------------------------------------

    def remove(self, value):

        """
        This function ...
        :param value: 
        :return: 
        """

        # Find key
        key = self.key_for_value(value)

        # Delete
        del self.contents[key]

    # -----------------------------------------------------------------

    def remove_all(self):

        """
        This function ...
        :return: 
        """

        self.contents = OrderedDict()

    # -----------------------------------------------------------------

    def pop(self, index_or_key):

        """
        This function ...
        :param index_or_key: 
        :return: 
        """

        # Integer type 'list[i]'
        if types.is_integer_type(index_or_key): key = self.keys[index_or_key]
        # Assume it is a proper key
        else: key = index_or_key

        # Pop with key
        return self.contents.pop(key)

    # -----------------------------------------------------------------

    @property
    def keys(self):

        """
        This function ...
        :return: 
        """

        return self.contents.keys()

    # -----------------------------------------------------------------

    @property
    def values(self):

        """
        This function ...
        :return: 
        """

        return self.contents.values()

    # -----------------------------------------------------------------

    @property
    def items(self):

        """
        This function ...
        :return: 
        """

        return self.contents.items()

    # -----------------------------------------------------------------

    def get_key(self, index_or_key):

        """
        This function ...
        :param index_or_key:
        :return:
        """

        # Integer type 'list[i]'
        if types.is_integer_type(index_or_key):

            # Get the key
            key = self.keys[index_or_key]
            return key

        # Assume it is a proper key
        else:

            # Check the key
            if index_or_key not in self.keys: raise ValueError("The element with key '" + str(index_or_key) + "' does not exist")

            # Return the key
            return index_or_key

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_key):

        """
        This function ...
        :param index_or_key:
        :return: 
        """

        # Slice 'list[a:b]'
        if isinstance(index_or_key, slice):

            # Create new class instance
            new = self.__class__()
            indices = index_or_key.indices(len(self))
            for index in range(indices[0], indices[1], indices[2]):
                # Get the key and value
                key = self.keys[index]
                value = self[index]
                new.contents[key] = value
            # Return the new list
            return new

        # Tuple 'list[a,b,c]'
        elif isinstance(index_or_key, tuple) or isinstance(index_or_key, list):

            #print("HERE")

            new = self.__class__()
            for i_or_k in index_or_key:
                # Get the actual key and value
                if types.is_integer_type(i_or_k): key = self.keys[i_or_k]
                else: key = i_or_k
                value = self[key]
                new.contents[key] = value
            # Return the new list
            return new

        # Integer type 'list[i]'
        elif types.is_integer_type(index_or_key):

            # Get the key
            key = self.keys[index_or_key]
            return self.contents[key]

        # Assume it is a proper key
        else: return self.contents[index_or_key]

    # -----------------------------------------------------------------

    def __setitem__(self, index_or_key, value):

        """
        Set an element
        :param index_or_key:
        :param value:
        """

        # Get the name
        if types.is_integer_type(index_or_key): key = self.keys[index_or_key]
        else: key = index_or_key

        # Replace
        self.replace(key, value)

    # -----------------------------------------------------------------

    def append(self, key, element):

        """
        This function ...
        :param key: 
        :param element: 
        :return: 
        """

        # Check the key
        if key in self.keys: raise ValueError("Already an element with the key '" + str(key) + "'")
        self.contents[key] = element

    # -----------------------------------------------------------------

    def replace(self, index_or_key, element, new_key=None):

        """
        This function ...
        :param index_or_key:
        :param element:
        :param new_key:
        :return: 
        """

        # Get the key
        key = self.get_key(index_or_key)

        # New key is given
        if new_key is not None:

            # Get the old value
            old = self.contents[key]

            # Replace the internal ordered dictionary
            self.contents = OrderedDict((new_key, element) if k == key else (k, value) for k, value in self.contents.items())

            # Return the old value
            return old

        # No new key, just replace the value
        else:

            # Get the old value
            old = self.contents[key]

            # Replace
            self.contents[key] = element

            # Return the old value
            return old

    # -----------------------------------------------------------------

    def __contains__(self, key):

        """
        This function ...
        :param key: 
        :return: 
        """

        return key in self.keys

    # -----------------------------------------------------------------

    def sort(self, cmp=None, key=None, reverse=False):

        """
        This function ...
        :param cmp: 
        :param key: 
        :param reverse: 
        :return: 
        """

        # Create key function
        if key is None: key_function = lambda x: x[1]
        else: key_function = lambda x: key(x[1])

        # Create cmp function
        #if cmp is None: cmp_function = None
        #else: cmp_function = lambda x,y: cmp_function(x[1], y[1]) # cmp_function(a,b)
        cmp_function = cmp

        # Sort the contents
        self.contents = OrderedDict(sorted(self.contents.iteritems(), cmp=cmp_function, key=key_function, reverse=reverse))

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return: 
        """

        # Create new
        new = self.__class__()

        # Add elements from this into new
        for key, element in self.items: new.append(key, element)

        # Return the new copy
        return new

# -----------------------------------------------------------------

class FilterBasedList(KeyList):

    """
    This function ...
    """

    @property
    def filters(self):

        """
        THis function ...
        :return: 
        """

        return self.keys

    # -----------------------------------------------------------------

    @property
    def filter_names(self):

        """
        This function ...
        :return: 
        """

        return [str(fltr) for fltr in self.filters]

    # -----------------------------------------------------------------

    def __contains__(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        return fltr in self.keys

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_filter):

        """
        This function ...
        :param index_or_filter:
        :return: 
        """

        # Slice 'list[a:b]'
        if isinstance(index_or_filter, slice):

            # Create new class instance
            new = self.__class__()
            indices = index_or_filter.indices(len(self))
            for index in range(indices[0], indices[1], indices[2]):
                # Get the key and value
                fltr = self.filters[index]
                value = self[fltr]
                new.contents[fltr] = value
            # Return the new list
            return new

        # Tuple 'list[a,b,c]'
        elif isinstance(index_or_filter, tuple) or isinstance(index_or_filter, list):

            # Initialize new
            new = self.__class__()
            for i_or_f in index_or_filter:
                if types.is_integer_type(i_or_f):
                    # Get the key and vlaue
                    fltr = self.filters[i_or_f]
                    value = self[fltr]
                    new.contents[fltr] = value
                elif types.is_string_type(i_or_f):
                    # Convert into filter
                    fltr = parse_filter(i_or_f)
                    value = self[fltr]
                    new.contents[fltr] = value
                elif isinstance(i_or_f, Filter):
                    value = self[i_or_f]
                    new.contents[i_or_f] = value
                else: raise ValueError("Invalid value in tuple: must be index, filter string or filter")
            # Return the new list
            return new

        # Integer type 'list[i]'
        elif types.is_integer_type(index_or_filter):

            # Get the key
            fltr = self.filters[index_or_filter]
            return self.contents[fltr]

        # Get the filter
        elif types.is_string_type(index_or_filter):

            # Convert into filter
            fltr = parse_filter(index_or_filter)
            return self.contents[fltr]

        # Filter
        elif isinstance(index_or_filter, Filter): return self.contents[index_or_filter]

        # Invalid index or filter
        else: raise ValueError("Invalid input: must be slice of indices, tuple of indices and/or filter (strings), single index, filter string or filter")

    # -----------------------------------------------------------------

    def __setitem__(self, index_or_filter, value):

        """
        This functio n...
        :param index_or_filter: 
        :param value: 
        :return: 
        """

        # Get the filters
        if types.is_string_type(index_or_filter): index_or_filter = parse_filter(index_or_filter)

        # Call the function of the base class
        self.replace(index_or_filter, value)

# -----------------------------------------------------------------

class FileList(KeyList):

    """
    This class ...
    :return:
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(FileList, self).__init__()

        # Add
        for key in kwargs: self.append(key, kwargs[key])

    # -----------------------------------------------------------------

    def append(self, key, path):

        """
        This function ...
        :param key:
        :param path:
        :return:
        """

        # Check whether file path
        if not fs.is_file(path): raise ValueError("Not an existing file")

        if key in self.keys: raise ValueError("Already a path with the key '" + str(key) + "'")
        self.contents[key] = path

# -----------------------------------------------------------------

class DirectoryList(KeyList):

    """
    This class ...
    :return:
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DirectoryList, self).__init__()

        # add the paths
        for key in kwargs: self.append(key, kwargs[key])

    # -----------------------------------------------------------------

    def append(self, key, path):

        """
        This function ...
        :param key:
        :param path:
        :return:
        """

        # Check whether file path
        if not fs.is_directory(path): raise ValueError("Not an existing directory")

        if key in self.keys: raise ValueError("Already a path with the key '" + str(key) + "'")
        self.contents[key] = path

# -----------------------------------------------------------------

class NamedList(KeyList):

    """
    This class ...
    """

    @property
    def names(self):

        """
        This function ...
        :return: 
        """

        return self.contents.keys()

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_name):

        """
        THis function ...
        :param index_or_name: 
        :return: 
        """

        # Slice 'list[a:b]'
        if isinstance(index_or_name, slice):

            # Create new class instance
            new = self.__class__()
            indices = index_or_name.indices(len(self))
            for index in range(indices[0], indices[1], indices[2]):
                # Get the key and value
                name = self.names[index]
                value = self[name]
                new.contents[name] = value
            # Return the new list
            return new

        # Tuple 'list[a,b,c]'
        elif isinstance(index_or_name, tuple) or isinstance(index_or_name, list):

            new = self.__class__()
            for i_or_n in index_or_name:
                if types.is_integer_type(i_or_n): name = self.names[i_or_n]
                elif types.is_string_type(i_or_n): name = i_or_n
                else: raise ValueError("Invalid item in tuple: must be index or name")
                value = self.contents[name]
                new[name] = value
            # Return the new list
            return new

        # Integer type 'list[i]'
        elif types.is_integer_type(index_or_name):

            # Get the name
            name = self.names[index_or_name]
            return self.contents[name]

        # String: name
        elif types.is_string_type(index_or_name):
            name = index_or_name
            return self.contents[name]

        # Invalid
        else: raise ValueError("Key must be slice, tuple, index or name")

    # -----------------------------------------------------------------

    def __setitem__(self, index_or_name, value):

        """
        Set an element
        :param index_or_name:
        :param value:
        """

        # Get the name
        if types.is_string_type(index_or_name): name = index_or_name
        elif types.is_integer_type(index_or_name): name = self.names[index_or_name]
        else: raise ValueError("Invalid index or name: must be integer or string")

        # Replace
        self.replace(name, value)

    # -----------------------------------------------------------------

    def append(self, name, element):

        """
        This function ...
        :param name: 
        :param element: 
        :return: 
        """

        #print(name)
        #print(self.names)
        if name in self.names: raise ValueError("Already an element with the name '" + name + "'")
        self.contents[name] = element

    # -----------------------------------------------------------------

    def replace(self, index_or_name, element, new_name=None):

        """
        This fucntion ...
        :param index_or_name:
        :param element:
        :param new_name:
        :return: 
        """

        # Get the name
        name = self.get_name(index_or_name)

        # Check the name
        if name not in self.names: raise ValueError("The element with name '" + name + "' does not exist")

        #return super(NamedList, self).replace(index_or_name, element, new_key=new_name)

        # If the element should get a new name (key)
        if new_name is not None:

            # Get the old value
            old = self.contents[name]

            # Replace the internal ordered dictionary
            self.contents = OrderedDict((new_name, element) if key == name else (key, value) for key, value in self.contents.items())

            # Return the old value
            return old

        # Don't give it a new name, just replace the value
        else:

            # Get the old value
            old = self.contents[name]

            # Replace
            self.contents[name] = element

            # Return the old value
            return old

    # -----------------------------------------------------------------

    def get_name(self, index_or_name):

        """
        This function ...
        :param index_or_name:
        :return:
        """

        # Get name
        if types.is_integer_type(index_or_name): name = self.names[index_or_name]
        elif types.is_string_type(index_or_name): name = index_or_name
        else: raise ValueError("Invalid index or name: " + str(index_or_name))

        return name

    # -----------------------------------------------------------------

    def pop(self, index_or_name):

        """
        This function ...
        :param index_or_name: 
        :return: 
        """

        # Get the name
        name = self.get_name(index_or_name)

        # Pop with name
        return self.contents.pop(name)

# -----------------------------------------------------------------

class NamedFileList(NamedList):

    """
    This class ... 
    :return: 
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(NamedFileList, self).__init__()

        # Add the paths
        for name in kwargs: self.append(name, kwargs[name])

    # -----------------------------------------------------------------

    def append(self, name, path):

        """
        This function ...
        :param name: 
        :param path: 
        :return: 
        """

        # Check whether file path
        if not fs.is_file(path): raise ValueError("Not an existing file")

        if name in self.names: raise ValueError("Already a path with the name '" + name + "'")
        self.contents[name] = path

# -----------------------------------------------------------------

class NamedDirectoryList(NamedList):

    """
    This class ... 
    :return: 
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(NamedDirectoryList, self).__init__()

        # Add the paths
        for name in kwargs: self.append(name, kwargs[name])

    # -----------------------------------------------------------------

    def append(self, name, path):

        """
        This function ...
        :param name: 
        :param path: 
        :return: 
        """

        # Check whether file path
        if not fs.is_directory(path): raise ValueError("Not an existing directory")

        if name in self.names: raise ValueError("Already a path with the path '" + name + "'")
        self.contents[name] = path

# -----------------------------------------------------------------

def equal_dictionaries(dict_a, dict_b):

    """
    This function ...
    :param dict_a:
    :param dict_b:
    :return:
    """

    if len(dict_a) != len(dict_b): return False

    for key in dict_a:

        if key not in dict_b: return False

        value_a = dict_a[key]
        value_b = dict_b[key]

        if value_a != value_b: return False

    return True

# -----------------------------------------------------------------