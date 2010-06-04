# -*- coding: utf-8 -*-
"""
Some utility classes that are useful in `Bio.GO` but are unlikely
to be used in other parts of BioPython.

Maybe this module should be merged with `Bio.utils` eventually.
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"

__all__ = ["pushback_iterator"]

# pylint:disable-msg=C0103
# C0103: invalid name
class pushback_iterator(object):
    """Constructs a pushback iterator from the given iterable.

    A pushback iterator is just like a regular iterable, but it
    has a `push_back` method that can be used to put some of
    the items back into the iterator temporarily.
    """

    def __init__(self, iterable):
        """Constructs a pushback iterator from the given iterable."""
        self.iterator = iter(iterable)
        self.buffer = []

    def __iter__(self):
        return self

    def next(self):
        """Returns the next item from the iterator, or the last
        one that was pushed back if there are items in the internal
        buffer."""
        try:
            return self.buffer.pop()
        except IndexError:
            return self.iterator.next()

    def push_back(self, item):
        """Pushes an item into the internal buffer so it will be
        returned when the user calls `next()` for the next time
        instead of the next element from the original iterable."""
        self.buffer.append(item)
