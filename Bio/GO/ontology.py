#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Classes for the Gene Ontology.

"""

__author__ = 'Chris Lasher'
__email__ = 'chris DOT lasher <AT> gmail DOT com'


# The number of digits a GO ID should be, as determined by the GO
# Consortium
NUM_GOID_DIGITS = 7


def _validate_goid(goid):
    """
    Validates a given GO ID.

    Raises a ValueError if `goid` is not a string of seven digits,
    optionally preceded by the prefix "GO:".

    """

    try:
        if goid.startswith('GO:'):
            digits = goid[3:]
        else:
            digits = goid

        if not digits.isdigit():
            raise ValueError("GO ID should contain only digits "
                    "or optionally digits prefixed with \"GO:\".")
        elif len(digits) != NUM_GOID_DIGITS:
            raise ValueError("GO ID should have precisely %d "
                    "digits." % (NUM_GOID_DIGITS))
    # If the goid doesn't support indexing or .isdigit, the user
    # needs to be told to give a string instead.
    except AttributeError, TypeError:
        raise ValueError("GO ID should be a string.")


class GOTerm(object):
    """
    A class to represent a term in the Gene Ontology.

    """

    def __init__(self):
        """

        :Parameters:
        - `goid`: the GO ID for the term. Should be a string of the form
          "GO:<digits>" or simply "<digits>", where digits is a
          zero-padded seven digit identifier
        - `name`: the name of the GO term
        - `ontology`: the ontology that the term belongs to [should be
          an `Ontology` instance]

        """

        self._validate_goid(goid)
        self.goid = goid
        self.name = name
        self.ontology = ontology


