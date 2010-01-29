#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Tests for the GO package.

"""

import unittest
from Bio import GO
from Bio.GO.Parsers import oboparser


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


class GOTermTests(unittest.TestCase):
    """Test cases for GO.ontology.GOTerm."""

    def test_repr(self):

        self.assertEqual(
                GO.ontology.GOTerm('GO:1234567', 'somethingase').__repr__(),
                '<GOTerm: GO:1234567-somethingase>'
        )


class GeneOntologyNXTests(unittest.TestCase):
    """Test cases for GO.ontology.GeneOntologyNX."""


    def setUp(self):
        self.ontology = GO.ontology.GeneOntologyNX('biological_process')
        GOTerm = GO.ontology.GOTerm
        self.terms = (
                GOTerm('GO:0008150', 'biological_process',
                    self.ontology),
                GOTerm('GO:0051704', 'multi-organism process',
                    self.ontology),
                GOTerm('GO:0043901', 'negative regulation of '
                    'multi-organism process', self.ontology)
            )


    def test_term_stored_conistently(self):

        term = self.terms[0]
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

        del self.ontology._goid_dict[term.goid]
        result = self.ontology._test_existence_in_internal_storage(term)
        self.assertEqual(
                result,
                (True, False)
        )


    def test_contains(self):

        term = self.terms[0]
        result = self.ontology.__contains__(term)
        self.assertEqual(
                result,
                False
        )

        self.ontology.add_term(term)
        result = self.ontology.__contains__(term)
        self.assertEqual(
                result,
                True
        )



    def test_contains_raises_InternalStorageInconsistentError(self):

        term = self.terms[0]
        self.ontology.add_term(term)
        del self.ontology._goid_dict[term.goid]
        self.assertRaises(
                GO.ontology.InternalStorageInconsistentError,
                self.ontology.__contains__,
                term
        )


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


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
