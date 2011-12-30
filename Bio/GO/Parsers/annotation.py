#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
Parser for Gene Ontology annotation files

Usage example::

    import Bio.GO.Parsers.oboparser as obo
    import Bio.GO.Parsers.annotation as annotation

    parser = obo.Parser(open("gene_ontology.1_2.obo"))
    ontology = parser.parse()

    parser = annotation.Parser(open("gene_association.sgd"))
    annotations = parser.parse(ontology)
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__version__ = "0.1"


__all__ = ["ParseError", "Parser"]

from Bio.GO.annotation import Annotation, AnnotationCollection
from Bio.GO.utils import pushback_iterator

class ParseError(Exception):
    """Exception thrown when a parsing error occurred"""

    def __init__(self, msg, lineno = 1):
        Exception.__init__("%s near line %d" % (msg, lineno))
        self.lineno = lineno

class Parser(object):
    """The main attraction, the annotation file parser.
    
    If you want to create a parser that reads an annotation file, do this:

      >>> import Bio.GO.Parsers.annotation as annotation
      >>> parser = annotation.Parser(open("gene_associations.sgd"))

    Only the headers are read when creating the parser. You can
    access these right after construction as follows:

      >>> parser.headers["gaf-version"]
      ['1.0']

    To read the annotations in the file, you must iterate over the
    parser as if it were a list. The iterator yields `Annotation`_
    objects. If you are not interested in the individual `Annotation`s
    and you only need an `Ontology` object, call the `parse()`
    method which will construct it for you.
    """

    # pylint:disable-msg=C0103
    def __init__(self, fp):
        """Creates an annotation parser that reads the given file-like object.
        """
        self.headers = {}

        if isinstance(fp, (str, unicode)):
            if fp.endswith(".gz"):
                # This is a gzipped file
                from gzip import GzipFile
                fp = GzipFile(fp)
            else:
                fp = open(fp)
        self._line_iterator = pushback_iterator(self._lines(fp))
        self._read_headers()

    @staticmethod
    def _lines(fp):
        """Iterates over the lines of a given file, removing
        comments and empty lines"""
        for line in fp:
            line = line.strip()
            if not line:
                continue
            yield line

    def _read_headers(self):
        """Reads the headers from the annotation file"""
        for line in self._line_iterator:
            if line[0] != '!':
                # We have reached the end of headers
                self._line_iterator.push_back(line)
                return

            parts = [part.rstrip() for part in line[1:].split(":", 1)]
            if len(parts) < 2:
                continue

            key, value = parts
            if key[0] == ' ':
                continue

            value = value.strip()
            try:
                self.headers[key].append(value)
            except KeyError:
                self.headers[key] = [value]

    # pylint:disable-msg=W0142
    # W0142: used * or ** magic
    def annotations(self):
        """Iterates over the annotations in this annotation file,
        yielding an `Annotation`_ object for each annotation."""
        for line in self._line_iterator:
            if line[0] == '!':
                continue
            parts = line.strip().split("\t")
            yield Annotation(*parts)

    def __iter__(self):
        return self.annotations()

    def parse(self, ontology):
        """Parses the file handle given during construction time and
        returns an appropriately constructed `AnnotationCollection`_
        instance.

        Annotations with qualifiers (such as ``NOT``) will be ignored.
        
        :Parameters:
            - `ontology`: the ontology being used to map term IDs to
              term names
        """
        return AnnotationCollection(annotation for annotation in self
                if not annotation.qualifiers)
