#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Classes for the Gene Ontology.

"""

__author__ = 'Chris Lasher'
__email__ = 'chris DOT lasher <AT> gmail DOT com'


# The number of digits a GO ID should be, as determined by the GO
# Consortium
NUM_GO_ID_DIGITS = 7


class Term(object):
    """A generic term class."""

    def __init__(self, identifier, name=None, ontology=None):
        """

        :Parameters:
        - `identifier`: a unique identifier for the term
        - `name`: an optional name for the term
        - `ontology`: the ontology that the term belongs to [should be
          an `Ontology` instance]

        """

        self.identifier = identifier
        self.name = name
        self.ontology = ontology


    def __repr__(self):

        outstr = "<%s: %s>" % (self.__class__.__name__, self.identifier)
        return outstr


    def __cmp__(self, other):

        return cmp(self.identifier, other.identifier)


class GOTerm(Term):
    """
    A class to represent a term in the Gene Ontology.

    """

    def __init__(self, identifier, name=None, ontology=None):
        """

        :Parameters:
        - `identifier`: the GO ID for the term. Should be a string of
          the form "GO:<digits>" or simply "<digits>", where digits is a
          zero-padded seven digit identifier
        - `name`: the name of the GO term
        - `ontology`: the ontology that the term belongs to [should be
          an `Ontology` instance]

        """

        identifier = _validate_and_normalize_go_id(identifier)
        super(GOTerm, self).__init__(identifier, name, ontology)


class GeneOntologyNX(object):
    """
    This class represents a gene ontology using NetworkX as the
    underlying graph framework.

    """

    def __init__(self, name, authority=None, identifier=None):
        """
        :Parameters:
        - `name`: name for the ontology
        - `authority`: the name of the authority for this ontology
        - `identifier`: an identifier for the ontology

        """
        self.name = name
        self.authority = authority
        self.identifier = identifier
        # The NetworkX directed graph will serve as the backbone for
        # operations.
        import networkx
        self._internal_dag = networkx.DiGraph()
        # We'll use this so we can retrieve terms by their identifier
        # strings, too.
        self._identifier_dict = {}


    def __repr__(self):

        outstr = "<%s: %s>" % (self.__class__.__name__, self.name)
        return outstr


    def add_term(self, term):
        """Add a term to the ontology.

        :Parameters:
        - `term`: a `GOTerm` instance

        """
        if term.identifier in self._identifier_dict:
            raise ValueError("Term %s already exists in ontology." %
                    term.identifier)
        self._identifier_dict[term.identifier] = term
        self._internal_dag.add_node(term)


    def get_term_by_id(self, term_id):
        """Retrieve a term from the ontology by its identifier.

        :Parameters:
        - `term_id`: a GO identifier (e.g., "GO:1234567")
        """
        return self._identifier_dict[term_id]


    def remove_term(self, term):
        """Add a term to the ontology.

        :Parameters:
        - `term`: a `GOTerm` instance

        """
        del self._identifier_dict[term.identifier]
        self._internal_dag.remove_node(term)


    def add_relationship(self, term1, term2, relationship_type):
        """
        Add a relationship between two terms to the ontology.

        Ontologies are composed of triples in the following form:

            `<SUBJECT> <PREDICATE> <OBJECT>`

        e.g., "mitochondrion is_a organelle"

        We represent this as `term1 relationship term2`.

        :Parameters:
        - `term1`: the subject term
        - `term2`: the object term
        - `relationship`: the predicate term (relationship type)

        """

        pass


    def remove_relationship(self, term1, term2, relationship_type):
        """
        Remove a relationship between two terms from the ontology.

        See `add_relationship()` for an explanation of the relationship
        structure.

        :Parameters:
        - `term1`: the subject term
        - `term2`: the object term
        - `type`

        """

        pass


    def orphaned_terms(self):
        """
        Returns an iterable of terms that have no relationship to any
        other terms in the ontology.

        """

        pass


def _validate_and_normalize_go_id(go_id):
    """
    Validates the format of a given GO identifier.

    Raises a ValueError if `go_id` is not a string of seven digits,
    optionally preceded by the prefix "GO:".

    Returns an identifier guaranteed to be prefixed with "GO:".

    """

    try:
        if go_id.startswith('GO:'):
            digits = go_id[3:]
            normalized_id = go_id
        else:
            digits = go_id
            normalized_id = 'GO:%s' % go_id

        if not digits.isdigit():
            raise ValueError("GO ID %s should contain only digits "
                    "or optionally digits prefixed with \"GO:\"." % (
                    go_id))
        elif len(digits) != NUM_GO_ID_DIGITS:
            raise ValueError("GO ID %s should have precisely %d "
                    "digits." % (go_id, NUM_GO_ID_DIGITS))
    # If the go_id doesn't support indexing or .isdigit, the user
    # needs to be told to give a string instead.
    except AttributeError, TypeError:
        raise ValueError("GO ID %s should be a string." % go_id)

    return normalized_id
