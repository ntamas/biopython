#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""Classes for the Gene Ontology."""

from Bio.Enum import Enum
from Bio.GO.utils import db_cursor_iterator, rewrite_db_query_markers

from collections import defaultdict
from inspect import getmodule
from itertools import izip
import sys

# Based on the work of Chris Lasher (chris DOT lasher <AT> gmail DOT com)
__author__ = 'Tamas Nepusz'
__email__ = 'ntamas <AT> gmail DOT com'

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

    def __init__(self, term_id):
        super(NoSuchTermError, self).__init__("no such term: %s" % term_id)


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

    def __eq__(self, other):
        return self.id == other.id

    def __ne__(self, other):
        return self.id != other.id

    __hash__ = None

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

    Relationships are considered **immutable** objects; equality of two
    relationships is checked by the equality of the corresponding
    `GOTerm` instances and the equality of the class itself. Relationships
    are also hashable, and the hash code depends only on the identity of
    the terms involved and the relationship itself. This is necessary so
    we can detect specific situations in `Bio.GO.inference.InferenceEngine`;
    for instance, when we inferred the same relationship twice in two
    different ways, we do not want to store two copies in the knowledge
    base.

    The bottom line is: **don't** modify `GORelationship` objects after
    they have been created, although theoretically you are able to do
    so.

    """

    # TODO: once BioPython drops support for Python 2.5, we should make
    # this class derived from namedtuple, that would solve all our problems
    # with immutability.

    __slots__ = ("subject_term", "object_term")
    names = ["any"]

    def __init__(self, subject_term, object_term):
        """Constructs a GO relationship between the given subject
        and object terms.

        :Parameters:
        - `subject_term`: the subject term. Must be an instance of `GOTerm`.
        - `object_term`: the object term. Must be an instance of `GOTerm`.
        """
        self.subject_term = subject_term
        self.object_term = object_term

    def __eq__(self, other):
        return self.subject_term == other.subject_term and \
               self.object_term == other.object_term and \
               self.__class__ == other.__class__

    def __hash__(self):
        return id(self.subject_term) ^ id(self.object_term) ^ id(self.__class__)

    def __repr__(self):
        return "%s(%r, %r)" % (self.__class__.__name__,
                self.subject_term, self.object_term)

    def __str__(self):
        return "%r %s %r" % (self.subject_term.name,
                             self.names[0],
                             self.object_term.name)

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


class GORelationshipFactory(object):
    """Factory class that keeps track of a mapping from GO relationship names
    to their corresponding classes.

    This class is usually not instantiated - all the methods are class methods,
    so you can call them without an instance.
    """

    _registry = {}

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

            >>> if issubclass(my_rel, GORelationshipFactory.from_name("inheritable")):
            ...     print "The relationship is inheritable"
        """
        try:
            return cls._registry[name.lower()]
        except KeyError:
            raise NoSuchRelationshipError(relation=name)

    @classmethod
    def get_known_names(cls):
        return cls._registry.keys()

    @classmethod
    def _register_relationships(cls, module):
        """Scans the given module for `GORelationship` subclasses
        and registers their names in the class-wide registry
        used by `from_name`.
        
        This method should be called once at the end of the
        initialization of this module.
        """
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

    def get_number_of_relationships(self):
        """Returns the number of relationships in this ontology."""
        raise NotImplementedError

    def get_number_of_terms(self):
        """Returns the number of terms in this ontology."""
        raise NotImplementedError

    def get_term_by_id(self, term_id):
        """Retrieves the given term from this ontology using its unique ID.

        Terms in an ontology have a primary ID and may have several alternative
        IDs. This method can be used to look up terms based on both their primary
        or their alternative IDs.

        Raises `NoSuchTermError` if the given term is not in this
        ontology.

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

    def terms(self):
        """Returns an iterable of all the terms in the ontology."""
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
            # Check whether networkx is new enough to support dicts
            # for nx.Graph.degree()
            if isinstance(networkx.Graph().degree(), list):
                raise ImportError
            self._nx = networkx
        except ImportError:
            raise ImportError("networkx >= 1.4 is required to use %s" % \
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
                term.id in self._internal_dag,
                term.id in self._goid_dict
            )
        return storage_states


    def has_term(self, term):
        """Checks whether the given term is in this ontology or not.

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
        self._internal_dag.add_node(term.id)

        # Register all the alternative IDs of this term in the
        # internal dict
        for alt_id in term.aliases:
            self._goid_dict[alt_id] = term


    def get_number_of_terms(self):
        """Returns the number of terms in this ontology."""
        return self._internal_dag.number_of_nodes()


    def get_number_of_relationships(self):
        """Returns the number of relationships in this ontology."""
        return self._internal_dag.number_of_edges()


    def get_term_by_id(self, term_id):
        """Retrieve a term from the ontology by its GO ID.

        This method also supports alternative IDs.

        Raises `NoSuchTermError` if the given term is not in this
        ontology.

        :Parameters:
        - `term_id`: a GO identifier (e.g., "GO:1234567")
        """
        try:
            return self._goid_dict[term_id]
        except KeyError:
            raise NoSuchTermError(term_id)


    def remove_term(self, term):
        """Removes the given term from this ontology.

        Raises `NoSuchTermError` if the given term is not in this
        ontology.

        :Parameters:
        - `term`: the term to remove; an instance of `GOTerm`.
        """
        try:
            del self._goid_dict[term.id]
            self._internal_dag.remove_node(term.id)
            term.ontology = None
        except KeyError:
            raise NoSuchTermError(term.id)


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
        self._internal_dag.add_edge(subject_term.id, object_term.id, \
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

        See `Ontology.get_relationships` for more details.
        """
        if subject_term is None and object_term is None:
            result = self._internal_dag.edges(data=True)
        elif object_term is None:
            result = self._internal_dag.out_edges(subject_term.id, data=True)
        elif subject_term is None:
            result = self._internal_dag.in_edges(object_term.id, data=True)
        else:
            result = self._internal_dag.out_edges(subject_term.id, data=True)
            result = [edge for edge in result if edge[1] == object_term.id]

        # Construct the Relationship objects here on-the-fly
        return [data["relationship"](
                self._goid_dict[src], self._goid_dict[dst]
            ) for src, dst, data in result]

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
            self._internal_dag.remove_edge(subject_term.id, object_term.id)
        except self._nx.exception.NetworkXError:
            raise NoSuchRelationshipError(subject_term, object_term)


    def orphaned_terms(self):
        """
        Returns an iterable of terms that have no relationship to any
        other terms in the ontology.
        """
        for term_id, degree in self._internal_dag.degree().iteritems():
            if degree == 0:
                yield self._goid_dict[term_id]


    def summary(self):
        """Returns a summary of the ontology as a string.

        The summary is similar to the contents of the Extended Info panel in
        OBO-Edit 2.1. It can be used to check whether a large ontology file
        has been parsed correctly; if the numbers you see here are equal to
        what you see when loading the ontology in OBO-Edit, chances are that
        the parsing was successful.
        """
        graph = self._internal_dag

        # Calculate the number of terms
        total_terms = len(graph)
        lines = ["Total terms = %d" % total_terms]

        # Find roots and out-degree distribution
        roots, deg_dist = [], {}
        max_degree = 5
        for term_id, degree in graph.out_degree().iteritems():
            degree = min(degree, max_degree+1)
            deg_dist[degree] = deg_dist.get(degree, 0) + 1
            if degree == 0:
                roots.append(term_id)

        # List roots
        for root_id in sorted(roots):
            root = self._goid_dict[root_id]
            rev_graph = graph.reverse()
            num_descendants = len(list(self._nx.dfs_preorder_nodes(rev_graph, root.id))) - 1
            lines.append("  %s (%s) has %d descendants" %
                         (root.name, root.id, num_descendants))

        # Degree distribution
        for degree in range(7):
            if degree > max_degree:
                degree_label = "> %d" % max_degree
            else:
                degree_label = degree
            plural = ["", "s"][degree != 1]
            count = deg_dist.get(degree, 0)
            percentage = "%d" % (count * 100 / total_terms)
            if percentage == "0" and count > 0:
                percentage = "< 1"
            lines.append("  terms with %s parent%s: %d (%s%%)" %
                    (degree_label, plural, count, percentage))

        # How many terms have definitions?
        def_terms = sum(1 for term_id in graph
                        if "def" in self._goid_dict[term_id].tags)
        lines.append("  %d%% of terms have definitions (%d of %d)" %
                ((def_terms * 100 / total_terms), def_terms, total_terms))

        return "\n".join(lines)


    def terms(self):
        """Returns an iterable of all the terms in the ontology."""
        return (self._goid_dict[term_id] for term_id in self._internal_dag)


class GeneOntologySQL(Ontology):
    """This class represents a gene ontology using a read-only database
    as a backend.

    Since the class uses the same cursor for all queries, it is not safe
    to use the same instance of the class from multiple threads.
    """

    def __init__(self, connection):
        """Creates a gene ontology that uses the given connection to a
        database.

        :Parameters:
        - `connection`: a connection to a database. The connection object
          must follow the semantics of the Python Database API
          Specification v2.0.
        """
        self.connection = connection
        self.cursor = connection.cursor()
        self._goid_dict = {}
        self._goid_to_dbid_dict = {}
        self._relid_to_class_dict = {}
        self._class_to_relid_dict = {}
        self._initialize_queries(connection)
        self.clear_cache()

    def _get_dbid_from_goid(self, goid):
        """Returns the database ID (primary key) corresponding to the GO term
        with the given GO ID"""
        try:
            return self._goid_to_dbid_dict[goid]
        except KeyError:
            self.cursor.execute(self._queries["dbid_from_goid"],
                                (goid, ))
            row = self.cursor.fetchone()
            if row is None:
                raise KeyError(goid)
            self._goid_to_dbid_dict[goid] = row[0]
            return row[0]

    def _initialize_queries(self, connection):
        """Initializes the SQL queries that will be used in this class.
        `connection` is the database connection that will be used."""
        module = getmodule(connection.__class__)
        if hasattr(module, "paramstyle"):
            paramstyle = module.paramstyle
        else:
            # Well, let's hope that the standard format strings work.
            # This should be OK for MySQL and PostgreSQL.
            paramstyle = "format"
        self._queries = dict(
            dbid_from_goid="SELECT id FROM term WHERE acc = %s LIMIT 1",
            rel_type_t1="SELECT term.acc, term2term.relationship_type_id "
                        "FROM term2term JOIN term "
                        "ON (term2term.term2_id = term.id) "
                        "WHERE term1_id = %s AND term.is_root != 1",
            rel_type_t1t2="SELECT relationship_type_id FROM term2term "
                          "WHERE term1_id = %s AND term2_id = %s",
            rel_type_t2="SELECT term.acc, term2term.relationship_type_id "
                        "FROM term2term JOIN term "
                        "ON (term2term.term1_id = term.id) "
                        "WHERE term2_id = %s AND term.is_root != 1",
            term_count="SELECT COUNT(*) FROM term WHERE acc = %s",
            term_by_id="SELECT id, acc, name FROM term WHERE acc = %s LIMIT 1",
            term_by_synonym="SELECT term.id, term.acc, term.name "
                            "FROM term JOIN term_synonym "
                            "ON term.id=term_synonym.term_id "
                            "WHERE term_synonym.acc_synonym = %s LIMIT 1"
        )
        for k, v in self._queries.iteritems():
            self._queries[k] = rewrite_db_query_markers(v, paramstyle)

    def _initialize_relid_dicts(self):
        """Initializes the cache mapping relationship type IDs to relationships
        using the database."""
        known_names = set(GORelationshipFactory.get_known_names())
        query = "SELECT id, name FROM term WHERE term_type = \"relationship\""\
                " OR term_type = \"gene_ontology\""
        self.cursor.execute(query)

        self._relid_to_class_dict = {}
        self._class_to_relid_dict = {}
        for row in db_cursor_iterator(self.cursor):
            if str(row[1]) in known_names:
                cls = GORelationshipFactory.from_name(row[1])
                self._class_to_relid_dict[cls] = int(row[0])
                self._relid_to_class_dict[int(row[0])] = cls

    def clear_cache(self):
        """Clears the internal caches of the ontology.

        In ordinary use-cases with an ontology that does not change, you
        should not have to call this method.
        """
        self._goid_dict = {}
        self._goid_to_dbid_dict = {}
        self._initialize_relid_dicts()

    def has_term(self, term):
        """Checks whether the given term is in this ontology or not.

        :Parameters:
        - `term`: a `GOTerm` instance or a GO term ID
        """
        if hasattr(term, "id"):
            term = term.id

        if term in self._goid_dict:
            return True

        self.cursor.execute(self._queries["term_count"], (term, ))
        return self.cursor.fetchone()[0] > 0

    def get_number_of_relationships(self):
        """Returns the number of relationships in this ontology."""
        query = "SELECT COUNT(*) FROM term2term JOIN term AS term1 "\
                "     ON (term2term.term1_id = term1.id) "\
                "     JOIN term AS term2 "\
                "     ON (term2term.term2_id = term2.id) "\
                "WHERE term1.acc LIKE 'GO:%' AND term2.acc LIKE 'GO:%'"
        self.cursor.execute(query)
        return int(self.cursor.fetchone()[0])

    def get_number_of_terms(self):
        """Returns the number of terms in this ontology."""
        query = "SELECT COUNT(*) FROM term WHERE acc LIKE 'GO:%' AND "\
                "is_obsolete=0"
        self.cursor.execute(query)
        return int(self.cursor.fetchone()[0])

    def get_term_by_id(self, term_id):
        """Retrieve a term from the ontology by its GO ID.

        This method does not support alternative IDs as the database schema
        itself does not support it. (In other words, synonyms have no GO IDs,
        therefore you cannot use these IDs when querying the database).

        Raises `NoSuchTermError` if the given term is not in this
        ontology.

        :Parameters:
        - `term_id`: a GO identifier (e.g., "GO:1234567")
        """
        if term_id in self._goid_dict:
            return self._goid_dict[term_id]

        self.cursor.execute(self._queries["term_by_id"], (term_id, ))
        row = self.cursor.fetchone()
        if row is None:
            # Hmmm, not found by primary ID. Maybe using a synonym?
            self.cursor.execute(self._queries["term_by_synonym"],
                    (term_id, ))
            row = self.cursor.fetchone()
            if row is None:
                raise NoSuchTermError(term_id)

        result = GOTerm(row[1], row[2], ontology=self)
        self._goid_dict[term_id] = result
        self._goid_to_dbid_dict[term_id] = row[0]
        return result

    def get_terms_by_ids(self, term_ids):
        """Retrieve multiple terms from the ontology by their GO IDs.

        :Parameters:
        - `term_id`: an iterable of GO identifiers
        
        :Returns:
        - an iterable yielding exactly the same number of elements as there
          are IDs in `term_ids`. The iterator will yield `None` for unknown
          GO identifiers.
        """
        result = []
        unknown_ids = defaultdict(set)
        for idx, term_id in enumerate(term_ids):
            if term_id in self._goid_dict:
                result.append(self._goid_dict[term_id])
            else:
                result.append(None)
                unknown_ids[term_id].add(idx)

        if not unknown_ids:
            return result

        condition = "acc IN (%s)" % (", ".join(repr(str(k))
            for k in unknown_ids.keys()))
        query = "SELECT id, acc, name FROM term WHERE %s" % condition
        self.cursor.execute(query)
        for row in db_cursor_iterator(self.cursor):
            term = GOTerm(row[1], row[2], ontology=self)
            self._goid_dict[term.id] = term
            self._goid_to_dbid_dict[term.id] = row[0]
            for idx in unknown_ids[term.id]:
                result[idx] = term

        for idx, term in enumerate(result):
            if not isinstance(term, GOTerm):
                result[idx] = None

        return result

    def get_relationships(self, subject_term=None, object_term=None):
        """Returns all the relationships that were added to the two
        given terms.

        :Parameters:
        - `subject_term`: the subject term; an instance of `GOTerm`.
        - `object_term`: the object term; an instance of `GOTerm`.

        :Returns:
        a (possibly empty) list of relationships where `subject_term`
        stands as the subject and `object_term` stands as the object.

        See `Ontology.get_relationships` for more details.
        """
        if subject_term is not None:
            if object_term is not None:
                # Querying for relationships between two given terms
                query = self._queries["rel_type_t1t2"]
                subject_dbid = self._get_dbid_from_goid(subject_term.id)
                object_dbid = self._get_dbid_from_goid(object_term.id)
                self.cursor.execute(query, (object_dbid, subject_dbid))
                return [self._relid_to_class_dict[row[0]](subject_term, object_term)
                        for row in db_cursor_iterator(self.cursor)]
            else:
                # Querying for all relationships where the subject is given
                query = self._queries["rel_type_t2"]
                subject_dbid = self._get_dbid_from_goid(subject_term.id)
                self.cursor.execute(query, (subject_dbid, ))

                rows = self.cursor.fetchall()
                object_terms = self.get_terms_by_ids(row[0] for row in rows)
                return [
                    self._relid_to_class_dict[row[1]](subject_term, object_term)
                    for row, object_term in izip(rows, object_terms)
                    if object_term is not None
                ]
        elif object_term is not None:
            # Querying for all relationships where the object is given
            query = self._queries["rel_type_t1"]
            object_dbid = self._get_dbid_from_goid(object_term.id)
            self.cursor.execute(query, (object_dbid, ))

            rows = self.cursor.fetchall()
            subject_terms = self.get_terms_by_ids(row[0] for row in rows)
            return [
                self._relid_to_class_dict[row[1]](subject_term, object_term)
                for row, subject_term in izip(rows, subject_terms)
                if subject_term is not None
            ]
        else:
            # Returning all known relationships.
            query = "SELECT term1.acc, term2.acc, term2term.relationship_type_id "\
                    "FROM term2term JOIN term AS term1 "\
                    "     ON (term2term.term1_id = term1.id) "\
                    "     JOIN term AS term2 "\
                    "     ON (term2term.term2_id = term2.id) "\
                    "WHERE term1.acc LIKE 'GO:%' AND term2.acc LIKE 'GO:%'"
            self.cursor.execute(query)

            rows = self.cursor.fetchall()
            subject_terms = self.get_terms_by_ids(row[1] for row in rows)
            object_terms = self.get_terms_by_ids(row[0] for row in rows)
            rels = [row[2] for row in rows]
            del rows

            return [
                self._relid_to_class_dict[rel](subject_term, object_term)
                for subject_term, object_term, rel in izip(subject_terms, object_terms, rels)
                if subject_term is not None and object_term is not None
            ]
        raise NotImplementedError

    def orphaned_terms(self):
        """Returns an iterable of terms that have no relationship to any
        other terms in the ontology.
        """
        query = "SELECT term.id, COUNT(term2term.term1_id) AS c "\
                "FROM term LEFT OUTER JOIN term2term "\
                "     ON (term.id=term2term.term2_id) "\
                "WHERE term.acc LIKE 'GO:%' "\
                "GROUP BY term.id HAVING c = 0"
        self.cursor.execute(query)
        ids = [row[0] for row in db_cursor_iterator(self.cursor)]
        if ids:
            ids = ", ".join(str(i) for i in ids)
            query = "SELECT term.id, term.acc, COUNT(term2term.term2_id) AS c "\
                    "FROM term LEFT OUTER JOIN term2term "\
                    "     ON (term.id=term2term.term1_id) "\
                    "WHERE term.id IN (%s) "\
                    "GROUP BY term.id HAVING c = 0" % ids
            self.cursor.execute(query)
            orphaned_ids = [row[1] for row in db_cursor_iterator(self.cursor)]
            return self.get_terms_by_ids(orphaned_ids)
        else:
            return []

    def terms(self):
        """Returns an iterable of all the terms in the ontology."""
        query = "SELECT id, acc, name FROM term "\
                "WHERE acc LIKE 'GO:%' AND is_obsolete=0"
        self.cursor.execute(query)
        for row in db_cursor_iterator(self.cursor):
            if row[1] in self._goid_dict:
                yield self._goid_dict[term.id]
            else:
                yield GOTerm(row[1], row[2], ontology=self)


def _validate_and_normalize_go_id(go_id):
    """
    Validates the format of a given GO identifier.

    Raises a ValueError if `go_id` is not a string of seven digits,
    optionally preceded by the prefix "GO:".

    Returns a GO ID guaranteed to be prefixed with "GO:".

    """

    try:
        if go_id.startswith('GO:') or True:
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
# of GORelationshipFactory
GORelationshipFactory._register_relationships(sys.modules[__name__])
