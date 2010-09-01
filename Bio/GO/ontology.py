#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""Classes for the Gene Ontology."""

from Bio.Enum import Enum
import sys

__author__ = 'Chris Lasher'
__email__ = 'chris DOT lasher <AT> gmail DOT com'


# The number of digits a GO ID should be, as determined by the GO
# Consortium
NUM_GO_ID_DIGITS = 7


class InternalStorageInconsistentError(Exception):
    """An exception to be raised when an ontology's internal storage
    data structures show inconsistent (conflicting) states of one or
    more terms being stored.

    """
    pass


class NoSuchTermError(Exception):
    """An exception to be raised when one tries to access or delete
    a term from an ontology when the term does not belong to the
    ontology.

    """

    def __init__(self, term):
        super(NoSuchTermError, self).__init__("no such term: %s" % term.id)


class NoSuchRelationshipError(Exception):
    """An exception to be raised when one tries to access or delete
    a relationship in an ontology when there is no such relationship
    between the two given terms, or when one tries to refer to a
    relationship by name in `GORelationship.from_name` and no such
    relationship name is known.
    """

    def __init__(self, subject_term=None, relation=None, object_term=None):
        if subject_term is None and object_term is None:
            if relation is None:
                message = "no such relationship"
            else:
                message = "no such relationship: %r" % relation
        else:
            if hasattr(subject_term, "id"):
                subject_term = str(subject_term.id)
            else:
                subject_term = repr(subject_term)
            if hasattr(object_term, "id"):
                object_term = str(object_term.id)
            else:
                object_term = repr(object_term)
            message = "no such relationship: %s --> %s" % \
                    (subject_term.id, object_term.id)
        super(NoSuchRelationshipError, self).__init__(message)


# pylint:disable-msg=W0232,R0903
# W0232: class has no __init__ method
# R0903: too few public methods
class Aspect(Enum):
    """Possible aspects of the Gene Ontology"""

    P = "Biological process"
    F = "Molecular function"
    C = "Cellular component"


class GOTerm(object):
    """
    A class to represent a term in the Gene Ontology.

    """

    __slots__ = ("id", "name", "aliases", "tags", "ontology")

    # pylint: disable-msg=C0103,R0913
    # C0103: invalid name
    # R0913: too many arguments
    def __init__(self, identifier, name=None, aliases=None, tags=None, \
                 ontology=None):
        """
        Creates a new Gene Ontology term.

        :Parameters:
        - `identifier`: the GO ID for the term. Should be a string of the
          form `'GO:<digits>'` or simply `'<digits>'`, where `<digits>`
          is a zero-padded seven digit identifier (e.g., `'GO:0006955'`)
        - `name`: the name of the GO term (e.g., `'immune response'`)
        - `aliases`: a list of alternative GO IDs for the term. Each item
          in the list should be a valid GO ID.
        - `tags`: a dict containing the tags associated with this Gene
          Ontology terms. These tags usually come from the stanza that
          defines the term in the ontology file.
        - `ontology`: the ontology which contains this GO term. In general,
          you should not pass anything other than ``None`` here, as the
          preferred way is to construct the term first and then add it
          later to an ontology using `Ontology.add_term()`. The latter
          should take care of setting the `ontology` field of the term
          properly.
        """
        self.id = _validate_and_normalize_go_id(identifier)

        if name:
            self.name = name
        else:
            self.name = ""

        if aliases is None:
            self.aliases = []
        else:
            self.aliases = [_validate_and_normalize_go_id(identifier) \
                            for identifier in aliases]

        if tags:
            self.tags = dict(tags)
        else:
            self.tags = {}

        self.ontology = ontology


    def __repr__(self):
        """String representation of a GO term"""
        return "%s(%r, %r, %r, %r, %r)" % (self.__class__.__name__,
                self.id, self.name, self.aliases, self.tags, self.ontology)

    def __str__(self):
        """Returns just the ID of the GO term"""
        return self.id


class GORelationship(object):
    """A generic class to represent a GO relationship between two terms
    in the ontology.

    """

    __slots__ = ("subject_term", "object_term")
    names = ["any"]

    def __init__(self, subject_term, object_term):
        """Constructs a GO relationship between the given subject
        and object terms.

        :Parameters:
        - `subject`: the subject term. Must be an instance of `GOTerm`.
        - `object`: the object term. Must be an instance of `GOTerm`.
        """
        self.subject_term = subject_term
        self.object_term = object_term

    def __eq__(self, other):
        return self.subject_term == other.subject_term and \
               self.object_term == other.object_term and \
               self.__class__ == other.__class__

    def __repr__(self):
        return "%s(%r, %r)" % (self.__class__.__name__,
                self.subject_term, self.object_term)

    @classmethod
    def from_name(cls, name):
        """Returns a subclass of `GORelationship` based on its
        human-readable name.

        Names are case insensitive. Currently the following names
        are understood: ``is_a``, ``part_of``, ``regulates``,
        ``negatively_regulates``, ``positively_regulates``. Spaces
        and underscores are equivalent.

        There are two additional special relationship names:
        ``any`` and ``inheritable``. ``any`` returns `GORelationship`,
        ``inheritable`` returns `InheritableGORelationship`. Since,
        for instance, `IsARelationship` is a subclass of `GORelationship`,
        it makes sense to perform tests such as::

            >>> if issubclass(my_rel, GORelationship.from_name("inheritable")):
            ...     print "The relationship is inheritable"
        """
        try:
            return cls._registry[name.lower()]
        except KeyError:
            raise NoSuchRelationshipError(relation=name)

    def implies(self, other):
        """Returns ``True`` if this relationship implies `other`.

        For instance, A ``negatively_regulates`` B implies that
        A ``regulates`` B, but not the other way round.

        The formal definition of this method is that the two terms
        involved in `self` and `other` should be equal, while the
        class of `other` must be a superclass of the class of `self`.
        """
        return self.subject_term == other.subject_term and \
               self.object_term == other.object_term and \
               isinstance(self, other.__class__)

    @classmethod
    def _register_relationships(cls, module):
        """Scans the given module for `GORelationship` subclasses
        and registers their names in the class-wide registry
        used by `from_name`.
        
        This method should be called once at the end of the
        initialization of this module.
        """
        if not hasattr(cls, "_registry"):
            cls._registry = {}
        for value in module.__dict__.itervalues():
            if not isinstance(value, type):
                continue
            if issubclass(value, GORelationship) and hasattr(value, "names"):
                for name in value.names:
                    name = name.lower()
                    cls._registry[name] = value
                    cls._registry[name.replace("_", " ")] = value


class InheritableGORelationship(GORelationship):
    """A generic class to represent an inheritable GO relationship.

    An inheritable relationship between two terms implies the subject
    term inherits properties from the object term. Currently, two
    relationships exist for this: the 'is_a' and 'part_of'
    relationships.
    """

    __slots__ = ()
    names = ["inheritable"]


class IsARelationship(InheritableGORelationship):
    """A class representing the 'is_a' GO relationship.

    This relationship is an inheritable relationship. If Term A "is_a"
    Term B, then everything that applies to Term B also applies to Term
    A.

    """

    __slots__ = ()
    names = ["is_a"]


class PartOfRelationship(InheritableGORelationship):
    """A class representing the 'part_of' GO relationship.

    This relationship is an inheritable relationship. If Term A is
    "part_of" Term B, then everything that applies to Term B also
    applies to Term A.

    """

    __slots__ = ()
    names = ["part_of"]


class RegulatesRelationship(GORelationship):
    """A class representing the 'regulates' GO relationship.

    This relationship is not an inheritable relationship.

    """

    __slots__ = ()
    names = ["regulates"]


class PositivelyRegulatesRelationship(RegulatesRelationship):
    """A class representing the 'positively_regulates' GO relationship.

    This relationship is an inheritable relationship. If Term A "is_a"
    Term B, then everything that applies to Term B also applies to Term
    A.

    """

    __slots__ = ()
    names = ["positively_regulates"]


class NegativelyRegulatesRelationship(RegulatesRelationship):
    """A class representing the 'negatively_regulates' GO relationship.

    This relationship is an inheritable relationship. If Term A "is_a"
    Term B, then everything that applies to Term B also applies to Term
    A.

    """

    __slots__ = ()
    names = ["negatively_regulates"]


class Ontology(object):
    """This abstract class specifies the interface that all the
    ontology classes should satisfy.
    """

    def add_term(self, term):
        """Adds the given `term` to this ontology.
        
        :Parameters:
        - `term`: the term to be added; an instance of `GOTerm`.
        """
        raise NotImplementedError

    def ensure_term(self, term_or_id):
        """Given a `GOTerm` or a GO term ID, returns the term itself.

        This method can be used in methods that expect a `GOTerm` to enable
        them to be able to work with GO term IDs as well."""
        if isinstance(term_or_id, GOTerm):
            return term_or_id
        return self.get_term_by_id(term_or_id)

    def get_term_by_id(self, term_id):
        """Retrieves the given term from this ontology using its unique ID.

        Terms in an ontology have a primary ID and may have several alternative
        IDs. This method can be used to look up terms based on both their primary
        or their alternative IDs.

        :Parameters:
        - `term_id`: the primary or an alternative ID of a term we are
          looking for.

        :Returns:
        the `GOTerm` corresponding to the given `term_id`.
        """
        raise NotImplementedError

    def has_term(self, term):
        """Checks whether the given term is in this ontology or not.
        
        :Parameters:
        - `term`: the term to look for; an instance of `GOTerm`.
        """
        raise NotImplementedError

    def remove_term(self, term):
        """Removes the given term from this ontology.

        Raises `NoSuchTermError` if the given term is not in this
        ontology.

        :Parameters:
        - `term`: the term to remove; an instance of `GOTerm`.
        """
        raise NotImplementedError

    def add_relationship(self, term1, term2, relationship):
        """Add a relationship between two terms to the ontology.

        Ontologies are composed of triples in the following form:

            `<SUBJECT> <PREDICATE> <OBJECT>`

        e.g., ``"mitochondrion is_a organelle"``

        We represent this as `term1 relationship term2`.

        :Parameters:
        - `subject_term`: the subject term; an instance of `GOTerm`.
        - `object_term`: the object term; an instance of `GOTerm`.
        - `relationship`: the predicate term (relationship type)
        """
        raise NotImplementedError

    def get_relationships(self, subject_term=None, object_term=None):
        """Returns all the relationships that were added to the two
        given terms.

        :Parameters:
        - `subject_term`: the subject term; an instance of `GOTerm`.
        - `object_term`: the object term; an instance of `GOTerm`.

        :Returns:
        a (possibly empty) list of relationships where `subject_term`
        stands as the subject and `object_term` stands as the object.

        If `subject_term` is ``None`` and `object_term` is not ``None``,
        all relationships will be retrieved where `object_term` is the
        object.

        If `subject_term` is not ``None`` and `object_term` is ``None``,
        all relationships will be retrieved where `subject_term` is
        the subject.

        If both `subject_term` and `object_term` are ``None``, all
        relationships will be retrieved.
        """
        raise NotImplementedError

    def has_relationship(self, subject_term, object_term, \
                         relationship=GORelationship):
        """Checks whether there exists a relationship between the
        two given terms.

        :Parameters:
        - `subject_term`: the subject term; an instance of `GOTerm`.
        - `object_term`: the object term; an instance of `GOTerm`.
        - `relationship`: the type of relationship we are interested in.
        """
        return any(isinstance(rel, relationship) \
                for rel in self.get_relationships(subject_term, object_term))

    def remove_relationship(self, subject_term, object_term, relationship):
        """
        Remove a relationship between two terms from the ontology.

        See `add_relationship()` for an explanation of the relationship
        structure.

        :Parameters:
        - `subject_term`: the subject term; an instance of `GOTerm`.
        - `object_term`: the object term; an instance of `GOTerm`.
        - `relationship`: the type of the relationship to be removed
          from those two terms.
        """
        raise NotImplementedError

    def __contains__(self, term):
        """Checks whether the given term is in this ontology or not.

        This method enables us to use an ontology object in Python
        expressions like ``if term in ontology:``.

        :Parameters:
        - `term`: the term to look for; an instance of `GOTerm`.
        """
        return self.has_term(term)


class GeneOntologyNX(Ontology):
    """This class represents a gene ontology using NetworkX as the
    underlying graph framework.

    """

    def __init__(self, name=None, authority=None, identifier=None):
        """
        :Parameters:
        - `name`: name for the ontology
        - `authority`: the name of the authority for this ontology
        - `identifier`: an identifier for the ontology

        """
        super(GeneOntologyNX, self).__init__()
        self.name = name
        self.authority = authority
        self.identifier = identifier

        # Store a reference to the NetworkX module here so we don't
        # have to import it all the time. We cannot import NetworkX
        # at the module level as the user might not have it.
        try:
            import networkx
            self._nx = networkx
        except ImportError:
            raise ImportError("networkx is required to use %s" % \
                    (self.__class__.__name__, ))

        # The NetworkX directed graph will serve as the backbone for
        # operations.
        self._internal_dag = self._nx.DiGraph()
        # We'll use this so we can retrieve terms by their GO ID
        # strings, too.
        self._goid_dict = {}


    def __repr__(self):
        outstr = "<%s: %s>" % (self.__class__.__name__, self.name)
        return outstr

    # pylint: disable-msg=C0103
    # C0103: invalid name
    def _test_existence_in_internal_storage(self, term):
        """Check on the state of storage of a given term within all the
        internal storage structures.

        Returns a tuple of storage states, where `True` represents that
        a term is stored by a storage structure, and `False` represents
        it is not stored.

        :Parameters:
        - `term`: a `GOTerm` instance

        """
        storage_states = (
                term in self._internal_dag,
                term.id in self._goid_dict
            )
        return storage_states


    def has_term(self, term):
        """Check to see if a term is present in the ontology.

        Raises `InternalStorageInconsistentError` in the event that
        internal storage shows inconsistent states of storage for the
        given term.

        :Parameters:
        - `term`: a `GOTerm` instance

        """
        storage_states = self._test_existence_in_internal_storage(term)
        # if all storage structures report existence, we're in a sane
        # state; return True
        if all(storage_states):
            return True
        # if all storage structures report no existence, we're in a sane
        # state; return False
        elif not any(storage_states):
            return False
        # if neither of those are true, something went horribly awry;
        # raise an error
        else:
            raise InternalStorageInconsistentError("Term %s has"
                    " inconsistent states of storage." % term)


    def add_term(self, term):
        """Add a term to the ontology.

        :Parameters:
        - `term`: a `GOTerm` instance

        """
        if term.id in self._goid_dict:
            raise ValueError("Term %s already exists in ontology." %
                    term.id)
        if term.ontology is not None:
            raise ValueError("Term %s is already added to another ontology." %
                    term.id)

        # Add the term to this ontology
        term.ontology = self
        self._goid_dict[term.id] = term
        self._internal_dag.add_node(term)

        # Register all the alternative IDs of this term in the
        # internal dict
        for alt_id in term.aliases:
            self._goid_dict[alt_id] = term


    def get_term_by_id(self, term_id):
        """Retrieve a term from the ontology by its GO ID.

        This method also supports alternative IDs.

        :Parameters:
        - `term_id`: a GO identifier (e.g., "GO:1234567")
        """
        return self._goid_dict[term_id]


    def remove_term(self, term):
        """Removes the given term from this ontology.

        Raises `NoSuchTermError` if the given term is not in this
        ontology.

        :Parameters:
        - `term`: the term to remove; an instance of `GOTerm`.
        """
        try:
            del self._goid_dict[term.id]
            self._internal_dag.remove_node(term)
            term.ontology = None
        except KeyError:
            raise NoSuchTermError(term)


    def add_relationship(self, subject_term, object_term, relationship):
        """Add a relationship between two terms to the ontology.

        Ontologies are composed of triples in the following form:

            `<SUBJECT> <PREDICATE> <OBJECT>`

        e.g., "mitochondrion is_a organelle"

        We represent this as `term1 relationship term2`.

        :Parameters:
        - `subject_term`: the subject term; an instance of `GOTerm`.
        - `object_term`: the object term; an instance of `GOTerm`.
        - `relationship`: the predicate term (relationship type)

        """
        # add the terms to the internal storage if they're not already
        # there
        for term in (subject_term, object_term):
            if term not in self:
                self.add_term(term)
        self._internal_dag.add_edge(subject_term, object_term, \
                                    relationship=relationship)

    def get_relationships(self, subject_term=None, object_term=None):
        """Returns all the relationships that were added to the two
        given terms.

        :Parameters:
        - `subject_term`: the subject term; an instance of `GOTerm`.
          ``None`` means all terms.
        - `object_term`: the object term; an instance of `GOTerm`.
          ``None`` means all terms.

        :Returns:
        a (possibly empty) list of relationships where `subject_term`
        stands as the subject and `object_term` stands as the object.

        If `subject_term` is ``None`` and `object_term` is not ``None``,
        all relationships will be retrieved where `object_term` is the
        object.

        If `subject_term` is not ``None`` and `object_term` is ``None``,
        all relationships will be retrieved where `subject_term` is
        the subject.

        If both `subject_term` and `object_term` are ``None``, all
        relationships will be retrieved.
        """
        if subject_term is None and object_term is None:
            result = self._internal_dag.edges(data=True)
        elif object_term is None:
            result = self._internal_dag.out_edges(subject_term, data=True)
        elif subject_term is None:
            result = self._internal_dag.in_edges(object_term, data=True)
        else:
            result = self._internal_dag.out_edges(subject_term, data=True)
            result = [edge for edge in result if edge[1] == object_term]

        # Construct the Relationship objects here on-the-fly
        return [data["relationship"](src, dst) for src, dst, data in result]

    def remove_relationship(self, subject_term, object_term, relationship):
        """
        Remove a relationship between two terms from the ontology.

        See `add_relationship()` for an explanation of the relationship
        structure.

        :Parameters:
        - `subject_term`: the subject term; an instance of `GOTerm`.
        - `object_term`: the object term; an instance of `GOTerm`.
        - `relationship`: the type of the relationship to be removed
          from those two terms.
        """
        try:
            self._internal_dag.remove_edge(subject_term, object_term)
        except self._nx.exception.NetworkXError:
            raise NoSuchRelationshipError(subject_term, object_term)


    def orphaned_terms(self):
        """
        Returns an iterable of terms that have no relationship to any
        other terms in the ontology.

        """
        #TODO
        pass


def _validate_and_normalize_go_id(go_id):
    """
    Validates the format of a given GO identifier.

    Raises a ValueError if `go_id` is not a string of seven digits,
    optionally preceded by the prefix "GO:".

    Returns a GO ID guaranteed to be prefixed with "GO:".

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
    except (AttributeError, TypeError):
        raise ValueError("GO ID %s should be a string." % go_id)

    return normalized_id


# Register the GO relationships in this module to the class registry
# of GORelationship
GORelationship._register_relationships(sys.modules[__name__])
