#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Tests for the GO package.

"""

import os
import unittest

from Bio import GO
from Bio.GO.annotation import EvidenceCode
from Bio.GO.inference import InferenceEngine, Rules, Query, UNBOUND
from Bio.GO.ontology import Aspect, GORelationship, \
                            GORelationshipFactory, GOTerm
from Bio.GO.utils import rewrite_db_query_markers
from Bio.GO.Parsers import oboparser as obo
from Bio.GO.Parsers import annotation
from cStringIO import StringIO


def create_sqlite_db(*args):
    """Creates an in-memory SQLite database and loads the fixtures
    given in the arguments. Returns a connection to the database."""
    import sqlite3

    db = sqlite3.connect(":memory:")
    db.executescript(open(
        os.path.join("GO", "sql", "schema.sql")
    ).read())
    db.commit()

    for fixture in args:
        import_fixtures(db, fixture)

    return db

def construct_relationship(term1, rel, term2, ontology=None):
    """Helper function, constructs a GO relationship from
    term IDs and relationship names."""
    if ontology:
        term1 = ontology.ensure_term(term1)
        term2 = ontology.ensure_term(term2)
    return GORelationshipFactory.from_name(rel)(term1, term2)

def first(iterable):
    """Helper function, returns the first element of an iterable or
    raises `IndexError` if the iterable is empty"""
    for elem in iterable:
        return elem
    raise IndexError

def import_fixture_table(db, fixture_name, table_name):
    """Imports the contents of the given table `table_name`
    from the given fixture `fixture_name` to the given
    SQLite database `db`."""
    fp = os.path.join("GO", "sql", "fixtures", fixture_name,
            "%s.txt" % table_name)
    fp = open(fp)

    cursor = db.cursor()
    for line in fp:
        parts = line.strip().split("\t")
        for idx, part in enumerate(parts):
            if part == "\\N":
                parts[idx] = None
        query = "INSERT INTO %s VALUES (%s)" % (
            table_name, ",".join(["?"] * len(parts))
        )
        cursor.execute(query, parts)
    db.commit()

def import_fixtures(db, fixture_name):
    """Imports all the tables from the given fixture `fixture_name`
    to the given SQLite database `db`."""
    fixture_dir = os.path.join("GO", "sql", "fixtures", fixture_name)
    for fname in os.listdir(fixture_dir):
        if fname[-4:].lower() == ".txt":
            import_fixture_table(db, fixture_name, fname[:-4])

def get_mini_ontology():
    """Loads the mini ontology from the test fixtures.
    
    You must *NOT* mutate the returned ontology in any way, doing so would
    result in cross-play between the test cases.
    """
    def load():
        parser = obo.Parser(file("GO/miniontology.obo"))
        return parser.parse()
    if not hasattr(load, "cache"):
        load.cache = load()
    return load.cache

def get_mini_ontology_sql():
    """Creates an in-memory SQLite database storing a mini ontology."""
    db = create_sqlite_db("init", "miniontology")
    return GO.ontology.GeneOntologySQL(db)

###########################################################################
## Tests for Bio.GO.ontology                                             ##
###########################################################################

class OntologyFunctionsTests(unittest.TestCase):

    def test_validate_and_normalize_go_id_raises_error(self):

        invalid_ids = [
                "Australopithecus",
                "GO:1",
                "GO:12345678",
                "1",
                "12345678",
                1234567
        ]
        for invalid_id in invalid_ids:
            self.assertRaises(
                    ValueError,
                    GO.ontology._validate_and_normalize_go_id,
                    invalid_id
            )


    def test_validate_and_normalize_go_id(self):

        cases = (('GO:1234567', '1234567'))
        expected = 'GO:1234567'
        for case in cases:
            self.assertEqual(
                    GO.ontology._validate_and_normalize_go_id(case),
                    expected
            )


class GORelationshipFactoryTests(unittest.TestCase):
    """Test cases for `GO.ontology.GORelationshipFactory`."""

    def test_from_name(self):
        names = ["is a", "is_a", "Is A", "part_of", "part Of",
                 "negatively_rEgUlAtEs", "positively_regulates",
                 "regulates", "negatively regulates",
                 "positively regulates"]
        for name in names:
            rel = GORelationshipFactory.from_name(name) 
            self.failUnless(issubclass(rel, GORelationship))
            self.assertNotEqual(rel, GORelationship)
            self.failUnless(name.lower().replace(" ", "_") in
                    [name.lower() for name in rel.names])

        self.assertRaises(GO.ontology.NoSuchRelationshipError,
                GORelationshipFactory.from_name, "no_such_rel")


class GORelationshipTests(unittest.TestCase):
    """Test cases for `GO.ontology.GORelationship`."""

    def test_implies(self):
        rels = [
            GO.ontology.PositivelyRegulatesRelationship("A", "B"),
            GO.ontology.RegulatesRelationship("A", "B"),
            GO.ontology.IsARelationship("A", "B"),
            GO.ontology.RegulatesRelationship("A", "C"),
            GO.ontology.InheritableGORelationship("A", "B")
        ]
        expected = [[ True,  True, False, False, False],
                    [False,  True, False, False, False],
                    [False, False,  True, False,  True],
                    [False, False, False,  True, False],
                    [False, False, False, False,  True]]
        for i, rel1 in enumerate(rels):
            for j, rel2 in enumerate(rels):
                result = rel1.implies(rel2)
                self.assertEqual(result, expected[i][j],
                        "rels[%d].implies(rels[%d]), "
                        "expected: %s, got: %s" %
                        (i, j, expected[i][j], result))


class GOTermTests(unittest.TestCase):
    """Test cases for GO.ontology.GOTerm."""

    def test_repr(self):

        self.assertEqual(
                GOTerm('GO:1234567', 'somethingase').__repr__(),
                "GOTerm('GO:1234567', 'somethingase', [], {}, None)"
        )


class GenericOntologyTestMixin(object):
    """Test cases for ontologies in general.

    This is a class that can be used as a mixin to test classes for
    descendants of `GO.ontology.Ontology`. The class contains tests for
    the basic interface of `GO.ontology.Ontology` assuming that the
    ontology is read-only, but it does not tell how to prepare an
    ontology instance to be tested. Override `prepare_ontology()` in
    subclasses to specify how to load the terms and relationships into
    the ontology instance that is to be tested. Also override
    `create_empty_ontology()` to specify how to create an empty
    ``self.ontology``.

    The `terms` and `relationships` class variables specify what
    terms and relationships we will need in the ontology.
    This ontology must be built in `prepare_ontology()`.
    """

    #: The terms we need in the ontology
    terms = (
        ('GO:0008150', 'biological_process'),
        ('GO:0051704', 'multi-organism process'),
        ('GO:0043901', 'negative regulation of multi-organism process'),
        ('GO:0003674', 'molecular function', ["GO:0005554"])
    )

    #: The relationships we need in the ontology
    relationships = (
        ("GO:0051704", "GO:0008150", "is_a"),
        ("GO:0043901", "GO:0051704", "negatively_regulates")
    )

    def create_empty_ontology(self):
        """Creates an empty ontology and assigns it to `self.ontology`.
        Should be overridden in child classes."""
        raise NotImplementedError

    def prepare_ontology(self):
        """Prepares an instance of the ontology to be tested. Should be
        overridden in child classes."""
        raise NotImplementedError

    def test_contains(self):
        self.prepare_ontology()

        term = self.ontology.get_term_by_id("GO:0008150")
        self.failUnless(term in self.ontology)

        term = GOTerm('GO:1234567', 'nonexistent term')
        self.failIf(term in self.ontology)

    def test_get_number_of_terms_and_relationships(self):
        self.prepare_ontology()
        self.assertEqual(self.ontology.get_number_of_terms(), 4)
        self.assertEqual(self.ontology.get_number_of_relationships(), 2)

    def test_get_term_by_id(self):
        self.prepare_ontology()
        for term in self.terms:
            result = self.ontology.get_term_by_id(term[0])
            self.failUnless(result.id == term[0])
            # Cannot test "result is term" here because that does not
            # hold for the SQL test cases where terms are constructed
            # on-the-fly as they are retrieved from the DB

    def test_get_term_by_id_aliases(self):
        self.prepare_ontology()
        for term in self.terms:
            if len(term) < 3:
                continue
            for alias in term[2]:
                result = self.ontology.get_term_by_id(alias)
                self.failUnless(result.id == term[0])

    def test_get_relationships_single(self):
        self.prepare_ontology()

        term1 = self.ontology.get_term_by_id("GO:0008150")
        term2 = self.ontology.get_term_by_id("GO:0051704")

        result = self.ontology.get_relationships(term2, term1)
        self.assertEquals(result, [GO.ontology.IsARelationship(term2, term1)])

        result = self.ontology.get_relationships(term1, term2)
        self.assertEquals(result, [])

    def test_get_relationships_given_subject(self):
        self.prepare_ontology()

        term1 = self.ontology.get_term_by_id("GO:0008150")
        term2 = self.ontology.get_term_by_id("GO:0051704")

        result = self.ontology.get_relationships(subject_term=term2)
        self.assertEquals(result, [GO.ontology.IsARelationship(term2, term1)])

        result = self.ontology.get_relationships(subject_term=term1)
        self.assertEquals(result, [])

    def test_get_relationships_given_object(self):
        self.prepare_ontology()

        term1 = self.ontology.get_term_by_id("GO:0008150")
        term2 = self.ontology.get_term_by_id("GO:0051704")
        term3 = self.ontology.get_term_by_id("GO:0043901")

        result = self.ontology.get_relationships(object_term=term1)
        self.assertEquals(result, [GO.ontology.IsARelationship(term2, term1)])

        result = self.ontology.get_relationships(object_term=term2)
        self.assertEquals(result,
                [GO.ontology.NegativelyRegulatesRelationship(term3, term2)])

    def test_get_relationships_all(self):
        self.prepare_ontology()
        
        returned_rels = self.ontology.get_relationships()
        missing_rels = list(self.relationships)
        self.assertEquals(len(returned_rels), len(missing_rels))

        for rel in returned_rels:
            rel = (rel.subject_term.id, rel.object_term.id, \
                   rel.__class__.names[0])
            missing_rels.remove(rel)
        self.assertEquals(missing_rels, [])

    def test_has_term(self):
        self.prepare_ontology()

        term = self.ontology.get_term_by_id("GO:0008150")
        self.failUnless(self.ontology.has_term(term))

        term = GOTerm('GO:1234567', 'nonexistent term')
        self.failIf(self.ontology.has_term(term))


class GeneOntologyNXTests(unittest.TestCase, GenericOntologyTestMixin):
    """Test cases for `GO.ontology.GeneOntologyNX`."""

    def setUp(self):
        self.create_empty_ontology()

    def create_empty_ontology(self):
        self.ontology = GO.ontology.GeneOntologyNX('biological_process')

    def prepare_ontology(self):
        """Prepares ``self.ontology`` by adding all the terms and
        relationships"""
        terms = [GOTerm(*term) for term in self.terms]
        for term in terms:
            self.ontology.add_term(term)
        for subject_id, object_id, rel in self.relationships:
            subject_term = first(term for term in terms \
                                 if term.id == subject_id)
            object_term = first(term for term in terms \
                                if term.id == object_id)
            rel = GORelationshipFactory.from_name(rel)
            self.ontology.add_relationship(subject_term, object_term, rel)

    def test_add_term(self):
        term = GOTerm('GO:0008150', 'biological_process')
        self.failIf(term in self.ontology)

        self.ontology.add_term(term)
        self.failUnless(term in self.ontology)

    def test_term_stored_consistently(self):
        term = GOTerm(*self.terms[0])
        result = self.ontology._test_existence_in_internal_storage(term)
        self.assertEqual(
                result,
                (False, False)
        )

        self.ontology.add_term(term)
        result = self.ontology._test_existence_in_internal_storage(term)
        self.assertEqual(
                result,
                (True, True)
        )

        del self.ontology._goid_dict[term.id]
        result = self.ontology._test_existence_in_internal_storage(term)
        self.assertEqual(
                result,
                (True, False)
        )


    def test_contains_raises_InternalStorageInconsistentError(self):
        term = GOTerm(*self.terms[0])
        self.ontology.add_term(term)
        del self.ontology._goid_dict[term.id]
        self.assertRaises(
                GO.ontology.InternalStorageInconsistentError,
                self.ontology.__contains__,
                term
        )

    def test_remove_term(self):
        term = GOTerm(*self.terms[0])
        self.assertRaises(
                GO.ontology.NoSuchTermError,
                self.ontology.remove_term, term
        )

        self.ontology.add_term(term)
        self.ontology.remove_term(term)
        self.failIf(term in self.ontology)
        self.failUnless(term.ontology is None)


class GeneOntologySQLTests(unittest.TestCase, GenericOntologyTestMixin):
    """Test cases for `GO.ontology.GeneOntologySQL`."""

    def setUp(self):
        self.create_empty_ontology()

    def create_empty_ontology(self):
        self.db = create_sqlite_db("init")
        self.ontology = GO.ontology.GeneOntologySQL(self.db)

    def prepare_ontology(self):
        """Prepares ``self.ontology`` by adding all the terms and
        relationships"""
        import_fixtures(self.db, "small")


#class OboparserFundamentalParsingTests(unittest.TestCase):

#    def test_value_comment(self):
#        case = 'alcohol dehydrogenase ! important for this weekend'
#        expected = 'alcohol dehydrogenase'
#        result = oboparser.value.parseString(case)[0]
#        self.assertEqual(result, expected)


#    def test_value_multi_line(self):
#        case = '''alcohol \\
#                dehydrogenase ! important for this weekend'''
#        expected = 'alcohol dehydrogenase'
#        result = oboparser.value.parseString(case)[0]
#        self.assertEqual(result, expected)


#    def test_tag_value_comment(self):
#        case = ('some_name: alcohol dehydrogenase ! important for this '
#                'weekend ')
#        expected = {
#                'tag': 'some_name',
#                'value': 'alcohol dehydrogenase',
#                'comment': 'important for this weekend'
#        }
#        result = oboparser.tag_value_pair.parseString(case).asDict()
#        self.assertEqual(result, expected)


#    def test_tag_value_comment_multi_line(self):
#        case = '''some_name: alcohol \
#                dehydrogenase ! important for this weekend '''
#        expected = {
#                'tag': 'some_name',
#                'value': 'alcohol dehydrogenase',
#                'comment': 'important for this weekend'
#        }
#        result = oboparser.tag_value_pair.parseString(case).asDict()
#        self.assertEqual(result, expected)

###########################################################################
## Tests for Bio.GO.Parsers.oboparser                                    ##
###########################################################################

class OboParserRealOntologyFileTests(unittest.TestCase):

    def test_mini_ontology(self):
        # Smoke testing
        parser = obo.Parser(file("GO/miniontology.obo"))
        ontology = parser.parse()
        self.assertEquals(list(ontology.orphaned_terms()), [])

        # Get a root term
        term = ontology.get_term_by_id("GO:0003674")
        self.failUnless(term.name == "molecular_function")
        self.failUnless(term.tags["namespace"] == [\
            obo.Value("molecular_function", ())\
        ])

        # Get the same term via one of its aliases
        term2 = ontology.get_term_by_id("GO:0005554")
        self.failUnless(term is term2)

        # Get another term
        term = ontology.get_term_by_id("GO:0003676")
        self.failUnless(term.name == "nucleic acid binding")
        self.failUnless(term.tags["def"] == [\
            obo.Value("Interacting selectively and non-covalently "+\
            "with any nucleic acid.", ("[GOC:jl]", ))\
        ])

        # Testing an is_a relationship
        term2 = ontology.get_term_by_id("GO:0005488")
        self.failUnless(ontology.has_relationship(term, term2))
        rel = GO.ontology.IsARelationship(term, term2)
        self.failUnless(rel in ontology.get_relationships(object_term=term2))
        self.failUnless(rel in ontology.get_relationships(subject_term=term))

        # Testing a part_of relationship
        term = ontology.get_term_by_id("GO:0006350")
        term2 = ontology.get_term_by_id("GO:0008152")
        self.failUnless(ontology.has_relationship(term, term2))
        rel = GO.ontology.PartOfRelationship(term, term2)
        self.failUnless(rel in ontology.get_relationships(object_term=term2))
        self.failUnless(rel in ontology.get_relationships(subject_term=term))


class AnnotationParserRealAnnotationFileTests(unittest.TestCase):

    def setUp(self):
        self.parser = annotation.Parser("GO/gene_association.PAMGO_Ddadantii.gz")
        self.annotations = list(self.parser)

    def test_headers(self):
        self.assertEqual(self.parser.headers["gaf-version"], ["2.0"])
        self.assertEqual(self.parser.headers["CVS Version"], ["Revision: 1.3 $"])
        self.assertEqual(len(self.parser.headers), 4)

    def test_db_and_taxon_fields(self):
        self.failIf(any(ann.db != "ASAP" for ann in self.annotations))
        self.failIf(any(ann.taxons != ("taxon:198628", )
                        for ann in self.annotations))

    def test_simple_annotation(self):
        ann = self.annotations[0]
        self.assertEqual(ann.db_object_id, "ABF-0014586")
        self.assertEqual(ann.db_object_symbol, "pelI")
        self.assertEqual(ann.qualifiers, ())
        self.assertEqual(ann.go_id, "GO:0030570")
        self.assertEqual(ann.db_references, ("PMID:9393696", ))
        self.assertEqual(ann.evidence_code, EvidenceCode.IDA)
        self.assertEqual(ann.froms, ())
        self.assertEqual(ann.aspect, Aspect.F)
        self.assertEqual(ann.db_object_name, "endo-pectate lyase")
        self.assertEqual(ann.db_object_synonyms, ())
        self.assertEqual(ann.db_object_type, "gene")
        self.assertEqual(ann.date, "20080227")
        self.assertEqual(ann.assigned_by, "ASAP")
        self.assertEqual(ann.annotation_extensions, ())
        self.failUnless(ann.gene_product_form_id is None)

    def test_annotation_with_qualifier(self):
        ann = self.annotations[7]
        self.assertEqual(ann.db_object_id, "ABF-0014783")
        self.assertEqual(ann.db_object_symbol, "pelX")
        self.assertEqual(ann.qualifiers, ("NOT", ))
        self.assertEqual(ann.go_id, "GO:0009405")
        self.assertEqual(ann.db_references, ("PMID:10049400", ))
        self.assertEqual(ann.evidence_code, EvidenceCode.IMP)
        self.assertEqual(ann.froms, ())
        self.assertEqual(ann.aspect, Aspect.P)
        self.assertEqual(ann.db_object_name, "exopolygalacturonate lyase PelX")
        self.assertEqual(ann.db_object_synonyms, ())
        self.assertEqual(ann.db_object_type, "gene")
        self.assertEqual(ann.date, "20080229")
        self.assertEqual(ann.assigned_by, "ASAP")
        self.assertEqual(ann.annotation_extensions, ())
        self.failUnless(ann.gene_product_form_id is None)

    def test_annotation_with_tuples(self):
        ann = self.annotations[15]
        self.assertEqual(ann.db_object_id, "ABF-0014958")
        self.assertEqual(ann.db_object_symbol, "pehX")
        self.assertEqual(ann.qualifiers, ())
        self.assertEqual(ann.go_id, "GO:0052179")
        self.assertEqual(ann.db_references, ("PMID:10564505", ))
        self.assertEqual(ann.evidence_code, EvidenceCode.IGI)
        self.assertEqual(ann.froms, ("ASAP:ABF-0014959", "ASAP:ABF-0014963"))
        self.assertEqual(ann.aspect, Aspect.P)
        self.assertEqual(ann.db_object_name, "Exo-poly-alpha-D-galacturonosidase precursor")
        self.assertEqual(ann.db_object_synonyms, ())
        self.assertEqual(ann.db_object_type, "gene")
        self.assertEqual(ann.date, "20080121")
        self.assertEqual(ann.assigned_by, "ASAP")
        self.assertEqual(ann.annotation_extensions, ())
        self.failUnless(ann.gene_product_form_id is None)


###########################################################################
## Tests for Bio.GO.inference                                            ##
###########################################################################

class QueryTests(unittest.TestCase):
    ontology = get_mini_ontology()

    def test_unbound_relation(self):
        self.assertRaises(ValueError, Query, self.ontology,
                          relation=UNBOUND)

    def test_term_id_resolution(self):
        query = Query(self.ontology, "GO:0005488", "is_a", "GO:0003676")
        self.assertEqual(query.relation, GO.ontology.IsARelationship)
        self.assertEqual(query.subject_term,
                         self.ontology.get_term_by_id("GO:0005488"))
        self.assertEqual(query.object_term,
                         self.ontology.get_term_by_id("GO:0003676"))
        self.assertRaises(GO.ontology.NoSuchTermError,
                Query, self.ontology,
                subject_term="no_such_term")


class InferenceRulesTests(unittest.TestCase):
    def get_mini_ontology(self):
        parser = obo.Parser(file("GO/miniontology.obo"))
        ontology = parser.parse()
        ontology.add_term(GOTerm("GO:0040008", "regulation of growth"))
        ontology.add_term(GOTerm("GO:0045926",
            "negative regulation of growth"))
        return ontology

    def test_default_ruleset(self):
        rules = Rules()
        ontology = self.get_mini_ontology()

        vacuole = ontology.get_term_by_id("GO:0005773")
        organelle = ontology.get_term_by_id("GO:0043226")
        cellular_component = ontology.get_term_by_id("GO:0005575")
        cytoplasm = ontology.get_term_by_id("GO:0005737")
        biological_process = ontology.get_term_by_id("GO:0008150")
        reg_biol_proc = ontology.get_term_by_id("GO:0050789")
        neg_reg_biol_proc = ontology.get_term_by_id("GO:0048519")
        growth = ontology.get_term_by_id("GO:0040007")
        reg_growth = ontology.get_term_by_id("GO:0040008")
        neg_reg_growth = ontology.get_term_by_id("GO:0045926")
        intracellular = ontology.get_term_by_id("GO:0005622")

        is_a = "is_a"
        part_of = "part_of"
        regulates = "regulates"
        neg_regulates = "negatively_regulates"

        # is_a o is_a = is_a
        rel1 = construct_relationship(vacuole, is_a, organelle)
        rel2 = construct_relationship(organelle, is_a, cellular_component)
        conclusion = construct_relationship(vacuole, is_a, cellular_component)
        self.assertEquals(rules.apply(rel1, rel2), conclusion)

        # is_a o regulates = regulates
        rel1 = construct_relationship(reg_growth, is_a, reg_biol_proc)
        rel2 = construct_relationship(reg_biol_proc, regulates, biological_process)
        conclusion = construct_relationship(reg_growth, regulates, biological_process)
        self.assertEquals(rules.apply(rel1, rel2), conclusion)

        # is_a o negatively_regulates = negatively_regulates
        rel1 = construct_relationship(neg_reg_growth, is_a, neg_reg_biol_proc)
        rel2 = construct_relationship(neg_reg_biol_proc, neg_regulates, biological_process)
        conclusion = construct_relationship(neg_reg_growth, neg_regulates,
                biological_process)
        self.assertEquals(rules.apply(rel1, rel2), conclusion)

        # part_of o is_a = part_of
        rel1 = construct_relationship(vacuole, part_of, cytoplasm)
        rel2 = construct_relationship(cytoplasm, is_a, cellular_component)
        conclusion = construct_relationship(vacuole, part_of, cellular_component)
        self.assertEquals(rules.apply(rel1, rel2), conclusion)

        # part_of o part_of = part_of
        rel1 = construct_relationship(vacuole, part_of, cytoplasm)
        rel2 = construct_relationship(cytoplasm, part_of, intracellular)
        conclusion = construct_relationship(vacuole, part_of, intracellular)
        self.assertEquals(rules.apply(rel1, rel2), conclusion)

        # negatively_regulates o is_a = negatively_regulates
        rel1 = construct_relationship(neg_reg_growth, neg_regulates, growth)
        rel2 = construct_relationship(growth, is_a, biological_process)
        conclusion = construct_relationship(neg_reg_growth, neg_regulates,
                biological_process)
        self.assertEquals(rules.apply(rel1, rel2), conclusion)

        # negatively_regulates o part_of = regulates
        # This is a fake example -- but the inference should be OK
        rel1 = construct_relationship(neg_reg_growth, neg_regulates, growth)
        rel2 = construct_relationship(growth, part_of, biological_process)
        conclusion = construct_relationship(neg_reg_growth, regulates,
                biological_process)
        self.assertEquals(rules.apply(rel1, rel2), conclusion)

        # is_a o is_a with totally unrelated relationships
        rel1 = construct_relationship(vacuole, is_a, organelle)
        rel2 = construct_relationship(reg_growth, is_a, reg_biol_proc)
        self.assertEquals(rules.apply(rel1, rel2), None)

        # part_of o regulates = ???
        # Again, the relations are fake and do not make sense
        rel1 = construct_relationship(vacuole, part_of, reg_growth)
        rel2 = construct_relationship(reg_growth, regulates, growth)
        self.assertEquals(rules.apply(rel1, rel2), None)

    def test_restrict_to_rules_relevant_for(self):
        rules = Rules()

        expected = {
            "is_a": ([("is_a", "any", 2)], set(["is_a"])),
            "part_of": ([
                ("part_of", "part_of", "part_of"),
                ("part_of", "is_a", "part_of"),
                ("is_a", "any", 2)
            ], set(["is_a", "part_of"])),
            "regulates": (rules.default_rules,
                set(["is_a", "part_of", "regulates"])),
            "negatively_regulates": ([
                ("regulates", "is_a", 1),
                ("is_a", "any", 2)
            ], set(["is_a", "negatively_regulates"])),
            "positively_regulates": ([
                ("regulates", "is_a", 1),
                ("is_a", "any", 2)
            ], set(["is_a", "positively_regulates"]))
        }

        for rel, (exp_rules, exp_rels) in expected.iteritems():
            restricted_rules, restricted_rels = \
                rules.restrict_to_rules_relevant_for(
                    GORelationshipFactory.from_name(rel)
                )
            message = "restrict_to_rules_relevant_for(%r) failed" % rel
            for rule, expected_rule in zip(restricted_rules.rules, exp_rules):
                self.failUnless(expected_rule[0] in rule[0].names, message)
                self.failUnless(expected_rule[1] in rule[1].names, message)
                if rule[2] == 1 or rule[2] == 2:
                    self.assertEqual(expected_rule[2], rule[2], message)
                else:
                    self.failUnless(expected_rule[2] in rule[2].names, message)
            exp_rels = set(GORelationshipFactory.from_name(rel) for rel in exp_rels)
            self.assertEquals(exp_rels, restricted_rels, message)


class InferenceEngineTestMixin(object):
    """Mixin class that should be mixed in to test classes that test
    `InferenceEngine` with various ontology classes.

    Theoretically, a test class that tests `InferenceEngine` has to
    be derived from `unittest.TestCase` as well as `InferenceEngineTests`,
    and must override the `setUp()` method to create ``self.ontology``
    and ``self.engine`` accordingly.
    """

    def check_inference(self, term1_id, relationship, term2_id, negate=False):
        """Checks whether we can infer from the ontology that
        the term with ID `term1_id` stands in the given `relationship`
        with the term having ID `term2_id`.

        The method will check it all the three possible ways: using an
        unbound object query, an unbound subject query and a fully bound
        query as well.
        """
        term1 = self.ontology.get_term_by_id(term1_id)
        term2 = self.ontology.get_term_by_id(term2_id)
        rel = "%r %s %r" % (term1.name, relationship, term2.name)

        # Try unbound object query
        query = Query(self.ontology, term1, relationship, UNBOUND)
        result = any(relationship in result.names and
                     result.subject_term == term1 and
                     result.object_term == term2
                     for result in self.engine.solve(query))

        if negate:
            self.failIf(result, "unbound object query inferred %s (which is false)" % rel)
        else:
            self.failUnless(result, "unbound object query didn't infer that %s" % rel)

        # Try unbound subject query
        query = Query(self.ontology, UNBOUND, relationship, term2)
        result = any(relationship in result.names and
                     result.subject_term == term1 and
                     result.object_term == term2
                     for result in self.engine.solve(query))

        if negate:
            self.failIf(result, "unbound subject query inferred %s (which is false)" % rel)
        else:
            self.failUnless(result, "unbound subject query didn't infer that %s" % rel)

        # Try fully bound query
        query = Query(self.ontology, term1, relationship, term2)
        result = self.engine.solve(query)

        if negate:
            self.failIf(result, "unbound query inferred %s (which is false)" % rel)
        else:
            self.failUnless(result, "unbound query didn't infer that %s" % rel)

    def check_not_provable(self, term1_id, rel, term2_id):
        self.check_inference(term1_id, rel, term2_id, negate=True)

    def test_infer_subject_is_a_simple(self):
        # Try to infer that cytoplasm is_a cellular_component
        self.check_inference("GO:0005737", "is_a", "GO:0005575")

    def test_infer_subject_is_a_twice(self):
        # Try to infer that transcription is_a metabolic process
        # since it is_a biosynthetic process and biosynthetic process
        # is_a metabolic process
        self.check_inference("GO:0006350", "is_a", "GO:0008152")

    def test_infer_subject_is_a_thrice(self):
        # Try to infer that transcription is_a biological process
        # since it is_a biosynthetic process, biosynthetic process
        # is_a metabolic process and metabolic process is_a
        # biological process
        # Bonus: we use an alias of biological_process here
        self.check_inference("GO:0006350", "is_a", "GO:0000004")

    def test_infer_subject_part_of(self):
        # Try to infer that transcription part_of metabolic process
        self.check_inference("GO:0006350", "part_of", "GO:0008152")

    def test_infer_subject_part_of_using_is_a(self):
        # Try to infer that transcription part_of biological process
        # since it part_of metabolic process and metabolic process
        # is_a biological process
        self.check_inference("GO:0006350", "part_of", "GO:0000004")

    def test_infer_subject_part_of_using_is_a_2(self):
        # Try to infer that lysosome part_of intracellular because
        # lysosome is_a vacuole, vacuole part_of cytoplasm and
        # cytoplasm part_of intracellular
        self.check_inference("GO:0005764", "part_of", "GO:0005622")

    def test_infer_subject_regulates(self):
        # Try to infer that regulation of biological process
        # regulates biological process
        self.check_inference("GO:0050789", "regulates", "GO:0008150")

    def test_cannot_infer_lysosome_part_of_vacuole(self):
        # We should not be able to infer that lysosome part_of vacuole,
        # nothing supports it
        self.check_not_provable("GO:0005764", "part_of", "GO:0005773")


class InferenceEngineNXTests(unittest.TestCase, InferenceEngineTestMixin):
    """Test cases for the default inference engine backed by a
    `GeneOntologyNX` instance."""

    def setUp(self):
        self.ontology = get_mini_ontology()
        self.engine = InferenceEngine()


class InferenceEngineSQLTests(unittest.TestCase, InferenceEngineTestMixin):
    """Test cases for the default inference engine backed by a
    `GeneOntologyNX` instance."""

    def setUp(self):
        self.ontology = get_mini_ontology_sql()
        self.engine = InferenceEngine()


###########################################################################
## Tests for Bio.GO.utils                                                ##
###########################################################################

class RewriteDbQueryMarkersTest(unittest.TestCase):
    def setUp(self):
        self.examples = [
            "",
            "%s",
            "SELECT * FROM term WHERE id = %s",
            "SELECT * FROM term WHERE id = %s AND acc = %s",
            "SELECT * FROM term WHERE id = %s AND acc = %s AND 1 = 1"
        ]

    def _test_queries(self, expected, style):
        for idx, (input, exp) in enumerate(zip(self.examples, expected)):
            obs = rewrite_db_query_markers(input, style)
            self.assertEquals(exp, obs,
                "query rewriting failed for style `%s' and query: `%s', "
                "expected: `%s', got: `%s'" % (style, input, exp, obs))

    def test_qmark(self):
        expected = ["", "?", "SELECT * FROM term WHERE id = ?",
                "SELECT * FROM term WHERE id = ? AND acc = ?",
                "SELECT * FROM term WHERE id = ? AND acc = ? AND 1 = 1"]
        self._test_queries(expected, "qmark")

    def test_numeric(self):
        expected = ["", ":1", "SELECT * FROM term WHERE id = :1",
                "SELECT * FROM term WHERE id = :1 AND acc = :2",
                "SELECT * FROM term WHERE id = :1 AND acc = :2 AND 1 = 1"]
        self._test_queries(expected, "numeric")

    def test_format(self):
        self._test_queries(self.examples, "format")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
