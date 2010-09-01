#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""Inference engine for drawing inferences on the Gene Ontology."""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__version__ = "0.1"


__all__ = ["InferenceEngine", "Rules", "Query", "UNBOUND"]

from Bio.GO.ontology import GORelationship, Ontology
from collections import deque

UNBOUND = None

class Query(object):
    """A query on the Gene Ontology.
    
    Queries consist of three components: a subject term, a relation
    (e.g., ``is_a`` or ``part_of``) and an object term. Either one
    of the subject or the object terms can be unbound (represented
    by the `UNBOUND`_ constant, which maps to the ``None`` constant
    in Python).
    
    If both the subject and the object term is bound to a specific
    GO term, the result of the query will be either ``True`` or
    ``False``, depending on whether the relation can be inferred
    from the Gene Ontology using the standard inference rules or
    not.
    
    If the subject term is unbound, the result of the query will be
    a list of `GOTerm`_ instances that stand in the given relation
    with the given object term.
    
    If the object term is unbound, the result of the query will be
    a list of `GOTerm`_ instances that stand in the given relation
    with the given subject term.

    This class has three instance variables that you may modify:
    `subject_term`, `relation` and `object_term`. They refer to the
    subject term, the relation and the object term, respectively.
    """

    def __init__(self, ontology, subject_term=UNBOUND,
                 relation="is_a", object_term=UNBOUND):
        """Constructs a query on the Gene Ontology with the given
        subject and object terms and the given relation.
        
        :Parameters:
        - `ontology`: the ontology on which the query is formulated.
        - `subject_term`: the subject term as a `GOTerm` instance or its
          GO ID, or ``UNBOUND`` if the subject term is unbound in the
          query.
        - `relation`: the relation (a subclass of `GORelationship`).
        - `object_term`: the object term as a `GOTerm` instance or its
          GO ID, or ``UNBOUND`` if the object term is unbound in the
          query.
        """
        self._subject_term = None
        self._relation = None
        self._object_term = None

        self.ontology = ontology
        self.relation = relation
        self.subject_term = subject_term
        self.object_term = object_term

    def _get_relation(self):
        return self._relation
    def _set_relation(self, relation):
        if relation is UNBOUND:
            raise ValueError("relations may not be unbound")
        if not isinstance(relation, type) or not issubclass(relation, GORelationship):
            relation = GORelationship.from_name(relation)
        self._relation = relation
    relation = property(_get_relation, _set_relation,
                        doc="The relation part of the query")

    def _get_subject_term(self):
        return self._subject_term
    def _set_subject_term(self, subject_term):
        if subject_term is not UNBOUND:
            subject_term = self.ontology.ensure_term(subject_term)
        self._subject_term = subject_term
    subject_term = property(_get_subject_term, _set_subject_term,
                       doc="The subject term of the query")

    def _get_object_term(self):
        return self._object_term
    def _set_object_term(self, object_term):
        if object_term is not UNBOUND:
            object_term = self.ontology.ensure_term(object_term)
        self._object_term = object_term
    object_term = property(_get_object_term, _set_object_term,
                       doc="The object term of the query")


class InvalidQueryError(Exception):
    """Exception thrown when a query is invalid for some reason."""

    def __init__(self, query, reason):
        super.__init__(self, reason)
        self.query = query


class Rules(object):
    """Class that stores information about the rules of inference in the
    Gene Ontology.

    The summary of the `Gene Ontology inference rules`_ defines *rules*
    in the following format: "if we know that A is in a given relationship
    R1 with B, and B is in a given relationship R2 with C, then we can
    infer that A is in a given relationship R3 with C". Therefore, our
    rules will be stored in the form of triplets of R1, R2, and R3. For
    instance, the triplet ``("is_a", "is_a", "is_a")`` means the following
    rule: "if A ``is_a`` B and B ``is_a`` C, then A ``is_a`` C". Similarly,
    the triplet ``("regulates", "part_of", "regulates")`` means that
    "if A ``regulates`` B and B ``part_of`` C, then A ``regulates`` C".

    Besides rules, we also have *rule templates*. Rule templates differ
    from rules in the third element of the triplet, which may be either
    of the integer constants 1 or 2. 1 means that R3 is identical to R1,
    whatever R1 is, while 2 means that R3 is identical to R2. This enables
    us to encode all the five relationships where R1 is ``is_a`` as follows:
    ``("is_a", "any", 2)``, meaning that if A ``is_a`` B and B stands in
    an arbitrary relationship R2 with C, then A stands in the very same
    relationship with C.

    The rules and rule templates are stored in the `rules` member variable.
    There is also a class variable called `default_rules`, which contains
    the default inference rules of the Gene Ontology. If you don't supply
    a ruleset at construction time, the default rules will be used.

    Rules must be listen in the order they should be checked for; in other
    words, more specific rules (such as ``("positively_regulates", "is_a",
    "positively_regulates)``) should come before more general rules (such as
    ``("regulates", "inheritable", "regulates")``).
    """

    # TODO: check whether it's worth indexing the rules based on R1

    default_rules = [
        ("regulates", "is_a", 1),
        ("regulates", "part_of", "regulates"),
        ("part_of", "inheritable", "part_of"),
        ("is_a", "any", 2),
    ]

    def __init__(self, rules=None):
        """Constructs a new ruleset with the given rules.

        If `rules` is ``None``, the default Gene Ontology inference rules
        will be used.
        """
        self.rules = []
        self._initialize(rules)

    def _initialize(self, rules=None):
        """Initializes the `rules` variable by transforming relationship
        names to `GORelationship` instances.

        If `rules` is empty or ``None``, the default rules will be used.
        """
        if not rules:
            rules = self.__class__.default_rules

        new_rules = []
        for rule in rules:
            target_rel = rule[2]
            if target_rel != 1 and target_rel != 2:
                target_rel = GORelationship.from_name(target_rel)
            new_rule = (GORelationship.from_name(rule[0]),
                        GORelationship.from_name(rule[1]),
                        target_rel)
            new_rules.append(new_rule)
        self.rules = new_rules

    def apply(self, rel1, rel2):
        """Applies the inference rules to the two given relationships
        `rel1` and `rel2` and returns the inferred relationship or
        ``None`` if no inference can be drawn.
        """
        if rel1.object_term != rel2.subject_term:
            return None

        for rule in self.rules:
            if isinstance(rel1, rule[0]) and isinstance(rel2, rule[1]):
                if rule[2] == 1:
                    rel3 = rel1.__class__
                elif rule[2] == 2:
                    rel3 = rel2.__class__
                else:
                    rel3 = rule[2]
                return rel3(rel1.subject_term, rel2.object_term)
        return None


class InferenceEngine(object):
    """Inference engine for the Gene Ontology."""

    def solve(query):
        """Solves the given query.
        
        :Parameters:
        - `query`: the query to be solved; an instance of `Query`.
        """
        if query.subject_term is UNBOUND and query.object_term is UNBOUND:
            raise InvalidQueryError(query,
                    "subject and object are both unbound")

        if query.object_term is UNBOUND:
            self.solve_unbound_object(query)

        raise NotImplementedError("queries of this type "
                "are not implemented yet")

    def solve_unbound_object(query):
        """Solves queries where the object is unbound and the subject is
        bound.

        You may call this method directly instead of `solve()`_ if you
        are sure that the query satisfies the preconditions of this method.
        """
        result = []
        # TODO
        return result

