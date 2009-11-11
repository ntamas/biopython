#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Tests for the GO package.

"""

import unittest
from Bio import GO


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


class TermTests(unittest.TestCase):

    def test_repr(self):

        self.assertEqual(
                GO.ontology.Term('GO:1234567').__repr__(),
                '<Term: GO:1234567>'
        )


    def test_cmp(self):

        cases = (
            ('1', '2', -1),
            ('1', '1', 0),
            ('2', '1', 1)
        )
        for id1, id2, expected in cases:
            term1 = GO.ontology.Term(id1)
            term2 = GO.ontology.Term(id2)
            self.assertEqual(
                    term1.__cmp__(term2),
                    expected
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
