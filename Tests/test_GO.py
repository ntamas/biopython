#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Tests for the GO package.

"""

import unittest
from Bio import GO

class GOTermTests(unittest.TestCase):

    def test_validate_goid_raises_error(self):
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
                    GO.ontology._validate_goid,
                    invalid_id
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
